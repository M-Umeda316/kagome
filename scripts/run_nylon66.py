"""Nylon-6,6 step-growth polycondensation via TDBB.

System:  n_diamines hexamethylenediamine + n_diacids adipic acid
Backend: OrbMol-v2 (default) or MACE-MP-0 via --backend mace
Ensemble: NPT (Langevin + MC barostat) at 300 K, 1 atm

Paper anchor: arXiv:2511.22874, PDF p.22, Table S2, Fig. S2, Fig. 4.

By default the system is built to paper density (0.5 g/mL, SI S-3) in one run:
direct placement, or dilute placement + classical-FF compression in-process
(--compress-backend classical), mirroring run_vinyl_aibn. Pass --box-size to
place at an explicit (dilute) edge without compression instead.

Usage:
    python scripts/run_nylon66.py --seed 7 --output-dir runs/nylon66
    python scripts/run_nylon66.py --seed 7 --backend mace --output-dir runs/nylon66_mace

    # well-mixed measurement mode (NOT paper-faithful): classical OpenMM/OpenFF
    # mixing after every cycle, refreshing the end-group neighbourhood that the
    # no-mixing paper-scale run exhausted at p~12% (decisions.md 2026-08-04).
    python scripts/run_nylon66.py --seed 7 --mix --mix-ps 25 \
        --mix-platform CUDA --output-dir runs/nylon66_mix
"""
from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path

os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')

import numpy as np

from scripts._mixing_cli import (
    add_mixing_arguments,
    mix_config_from_args,
    mixing_setup_from_args,
    mixing_setup_mismatch,
    resolve_mixing_args,
)
from scripts._systems import (
    _DIACID_SMILES,
    _DIAMINE_SMILES,
    box_from_density,
    build_nylon66_system,
    layout_bonds,
)
from kagome.backends.base import Calculator
from kagome.boost.tdbb import TDBBParams
from kagome.integrators.init_velocities import maxwell_boltzmann_velocities
from kagome.integrators.langevin import LangevinIntegrator, LangevinParams
from kagome.integrators.mc_barostat import MCBarostat, MCBarostatParams
from kagome.reactive.bonds import BondTracker
from kagome.workflows.polymerization import (
    PolymerizationConfig,
    PolymerizationWorkflow,
    SimulationState,
    load_checkpoint,
    masses_from_species,
)

logging.basicConfig(level=logging.INFO, format='%(name)s | %(message)s')
logger = logging.getLogger(__name__)


def _create_backend(backend: str, device: str, model: str,
                    compile_model: bool = False,
                    empty_cache: bool = True) -> Calculator:
    if backend == 'orb':
        from kagome.backends.orb_backend import create_orb_calculator
        return create_orb_calculator(device=device, compile=compile_model,
                                     empty_cache=empty_cache)
    else:
        if compile_model or not empty_cache:
            logger.warning('--compile/--no-empty-cache only apply to the orb '
                           'backend; ignored for %r.', backend)
        from kagome.backends.mace_backend import create_mace_calculator
        return create_mace_calculator(model=model, device=device)


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Nylon-6,6 step-growth polycondensation (TDBB)'
    )
    parser.add_argument('--seed', type=int, default=7)
    parser.add_argument('--output-dir', type=Path, default=Path('runs/nylon66'))
    parser.add_argument('--n-diamines', type=int, default=10)
    parser.add_argument('--n-diacids', type=int, default=10)
    parser.add_argument('--n-cycles', type=int, default=3)
    parser.add_argument('--biased-steps', type=int, default=500)
    parser.add_argument('--unbiased-steps', type=int, default=500)
    parser.add_argument('--f2', type=float, default=2.0,
                        help='TDBB Gaussian width f2 (Å⁻²). Default 2.0 for OrbMol-v2: the '
                             'paper value 10.0 leaves the amine_N–carboxyl_C formation bias '
                             'in a capture-shell dead-zone (force ≈0 across the [3,6] Å '
                             'candidate window) so no amide bond forms; f2=2 bridges it, '
                             'mirroring the vinyl/MA recipe. For a paper-faithful run pass '
                             '--f2 10.0 explicitly. See decisions.md 2026-07-08 / 2026-07-30.')
    # Pre-TDBB relaxation of the classical-compressed dense structure. The
    # compressed 0.5 g/mL box carries close intermolecular contacts (fmax
    # ~60 kcal/mol/Å, unconverged) whose extreme forces segfault the MLIP on the
    # first biased evaluation. FIRE minimize + unbiased equilibration relax them
    # first, exactly as run_vinyl_aibn does (paper anchor PDF p.20: equilibration
    # precedes production reactive MD). Nylon is step-growth condensation with no
    # activation phase, so the order is simply: build → compress → minimize →
    # equilibrate → TDBB. Both steps are driven by PolymerizationConfig and run
    # inside wf.run() (workflows/polymerization.py run() lines 633-636), matching
    # vinyl's non-activation path.
    parser.add_argument('--minimize', dest='minimize', action='store_true', default=True,
                        help='FIRE energy minimization before TDBB (default: on). '
                             'Relaxes close contacts in the compressed structure '
                             '(paper anchor PDF p.20).')
    parser.add_argument('--no-minimize', dest='minimize', action='store_false',
                        help='Skip pre-TDBB energy minimization.')
    parser.add_argument('--minimize-fmax', type=float, default=1.0,
                        help='FIRE convergence threshold (kcal/mol/Å). Default 1.0.')
    parser.add_argument('--equil-steps', type=int, default=2000,
                        help='Unbiased NPT equilibration steps before TDBB '
                             '(paper anchor PDF p.20; length not specified, default 2000 '
                             '= 500 fs matching a TDBB block). 0 disables.')
    parser.add_argument('--box-size', type=float, default=None,
                        help='Box edge (Å). If omitted, computed from --density and '
                             'reached by direct placement or classical compression '
                             '(single-run path to paper density, mirroring run_vinyl_aibn).')
    parser.add_argument('--density', type=float, default=0.5,
                        help='Initial density (g/mL). Paper SI S-3 uses 0.5. Used only '
                             'when --box-size is omitted.')
    parser.add_argument('--temperature', type=float, default=300.0)
    parser.add_argument('--friction-per-fs', type=float, default=0.01,
                        help='Langevin friction (1/fs). Default 0.01 as the cooling '
                             'lever for the OrbMol f2=2 recipe: the TDBB bias work + '
                             'condensation exotherm accumulates under the paper-faithful '
                             '0.001, the same overheating seen in the epoxy 549 K smoke '
                             '(decisions.md 2026-07-07 friction / 2026-07-30). Pass 0.001 '
                             'to restore the paper-faithful value.')
    parser.add_argument('--timestep-fs', type=float, default=0.25,
                        help='MD timestep (fs). Default 0.25 fs (conservative, validated '
                             'for FIRE densification + ML NVT; current behaviour). 1.0 fs '
                             'is standard for organic ML MD and gives 4x speed for the '
                             'same physical time.')
    parser.add_argument('--pressure', type=float, default=1.0)
    parser.add_argument('--no-barostat', action='store_true')
    parser.add_argument('--backend', type=str, default='orb',
                        choices=['orb', 'mace'])
    parser.add_argument('--device', type=str, default='cpu')
    parser.add_argument('--compile', action='store_true', default=False,
                        help='torch.compile the OrbMol model (orb backend only, '
                             'Linux/WSL). Cuts the kernel-launch CPU overhead that '
                             'py-spy measured as the paper-scale bottleneck '
                             '(decisions.md 2026-07-14). First evaluation compiles '
                             'for minutes; graph-size changes may recompile.')
    parser.add_argument('--no-empty-cache', dest='empty_cache',
                        action='store_false', default=True,
                        help='Skip per-step torch.cuda.empty_cache() (orb backend, '
                             'cuda only; ~9%% of CPU time). Safe on >=32 GB GPUs; '
                             'keep the default on 16 GB, where allocator '
                             'fragmentation exhausts VRAM (decisions.md 2026-06-15).')
    parser.add_argument('--model', type=str, default='small',
                        help='MACE model size (only used with --backend mace)')
    parser.add_argument('--compress-backend', type=str, default='classical',
                        choices=['classical', 'ml'],
                        help='Calculator for box compression to paper density. '
                             '"classical" (default) uses OpenMM/OpenFF Sage so densification '
                             'does not consume MLIP GPU time (decision 2026-06-20); "ml" uses '
                             'the production MLIP. Only used when --box-size is omitted.')
    parser.add_argument('--compress-platform', type=str, default='CPU',
                        choices=['CPU', 'CUDA', 'OpenCL', 'Reference'],
                        help='OpenMM platform for --compress-backend classical '
                             '(default CPU, keeps the GPU free for the MLIP MD).')
    # Cycle-boundary checkpointing for crash recovery on long (paper-scale) runs,
    # mirroring run_vinyl_aibn. By default a checkpoint is written every cycle so
    # a killed run can --resume from the last completed cycle. Nylon is step-growth
    # with no activation/spin, so wf.run() restores positions/groups/tracker/rng
    # from the checkpoint directly — no extra state to carry.
    parser.add_argument('--resume', action='store_true', default=False,
                        help='Resume from <output-dir>/checkpoint.pkl if present: '
                             'skip build/minimize/equilibration and continue from the '
                             'next cycle after the last checkpoint (trajectory/bonds/'
                             'topology are appended). No-op if no checkpoint exists.')
    parser.add_argument('--no-checkpoint', action='store_true', default=False,
                        help='Disable writing <output-dir>/checkpoint.pkl each cycle '
                             '(checkpointing is on by default for resumable long runs).')
    # WM-P3 mixing stage (--mix + 6 knobs), shared with run_vinyl_copolymer /
    # run_epoxy_amine via scripts/_mixing_cli.py (decisions.md 2026-08-04).
    add_mixing_arguments(parser)
    args = parser.parse_args()

    resolve_mixing_args(parser, args)

    rng = np.random.default_rng(args.seed)

    counts = {_DIAMINE_SMILES: args.n_diamines, _DIACID_SMILES: args.n_diacids}

    # Backend is created before the build so the 'ml' compress option can reuse it.
    calc = _create_backend(args.backend, args.device, args.model,
                           compile_model=args.compile,
                           empty_cache=args.empty_cache)
    logger.info('Backend: %s', calc.name)

    def _build(edge: float, gen):
        return build_nylon66_system(
            n_diamines=args.n_diamines,
            n_diacids=args.n_diacids,
            box_size=edge,
            rng=gen,
            rdkit_seed=args.seed,  # RF23: conformer geometry follows --seed
        )

    if args.box_size is not None:
        # Explicit box edge (legacy behaviour): direct placement, no compression.
        logger.info(
            'Building nylon-6,6: %d diamines + %d diacids in %.1f Å box (explicit)...',
            args.n_diamines, args.n_diacids, args.box_size,
        )
        positions, species, template, groups = _build(args.box_size, rng)
        cell = np.diag([args.box_size, args.box_size, args.box_size])
    else:
        # Single-run path to paper density (mirrors run_vinyl_aibn): direct
        # placement at the target density, else dilute placement + compression.
        target_edge = box_from_density(counts, args.density)
        logger.info(
            'Building nylon-6,6: %d diamines + %d diacids, target box %.1f Å '
            '(%.2f g/mL, paper SI S-3)...',
            args.n_diamines, args.n_diacids, target_edge, args.density,
        )
        try:
            positions, species, template, groups = _build(target_edge, rng)
            cell = np.diag([target_edge, target_edge, target_edge])
        except RuntimeError:
            logger.warning(
                'Direct placement at %.2f Å failed — placing dilute then compressing.',
                target_edge,
            )
            place_edge = None
            for place_density in (0.25, 0.20, 0.15, 0.10):
                edge = box_from_density(counts, place_density)
                if edge <= target_edge:
                    continue
                try:
                    positions, species, template, groups = _build(
                        edge, np.random.default_rng(args.seed))
                    place_edge = edge
                    logger.info(
                        'Placed at dilute density %.2f g/mL (box %.2f Å); compressing to %.2f Å.',
                        place_density, edge, target_edge,
                    )
                    break
                except RuntimeError:
                    continue
            if place_edge is None:
                raise RuntimeError(
                    'Could not place nylon even at dilute density 0.10 g/mL.'
                )
            from kagome.backends.classical_backend import make_compress_calculator
            from kagome.integrators.minimize import compress_box
            from kagome.prep.openmm_equilibrate import MoleculeSpec
            # Placement order in build_nylon66_system: diamines first (seed),
            # then diacids (seed+1). MoleculeSpec order/seeds must match (RF23).
            specs = [
                MoleculeSpec(_DIAMINE_SMILES, args.n_diamines, rdkit_seed=args.seed),
                MoleculeSpec(_DIACID_SMILES, args.n_diacids, rdkit_seed=args.seed + 1),
            ]
            compress_calc = make_compress_calculator(
                args.compress_backend, specs, calc, platform=args.compress_platform,
                target_edge_A=target_edge,
            )
            place_cell = np.diag([place_edge, place_edge, place_edge])
            result = compress_box(positions, place_cell, target_edge, species, compress_calc)
            positions, cell = result.positions, result.cell

    logger.info(
        'System: %d atoms total  (%d amine_N, %d carboxyl_C, %d amine_H, %d carboxyl_OH)',
        len(species),
        len(groups['amine_N'].atom_indices),
        len(groups['carboxyl_C'].atom_indices),
        len(groups['amine_H'].atom_indices),
        len(groups['carboxyl_OH'].atom_indices),
    )

    # Initial box edge for provenance (state.cell evolves under NPT during the run).
    initial_box_edge_A = float(cell[0, 0])

    langevin_params = LangevinParams(
        temperature_K=args.temperature, friction_per_fs=args.friction_per_fs)
    config = PolymerizationConfig(
        timestep_fs=args.timestep_fs,
        biased_steps=args.biased_steps,
        unbiased_steps=args.unbiased_steps,
        n_cycles=args.n_cycles,
        tdbb=TDBBParams(
            f2=args.f2,
            gamma=1.0,
            f1_max_formation=250.0,
            f1_max_dissociation=125.0,
            lambda_vdw=0.60,
        ),
        seed=args.seed,
        save_interval=50,
        minimize=args.minimize,
        minimize_fmax=args.minimize_fmax,
        equil_steps=args.equil_steps,
        mixing=mix_config_from_args(args, args.temperature),
    )

    integrator = LangevinIntegrator(langevin_params)
    tracker = BondTracker()
    barostat = None if args.no_barostat else MCBarostat(
        MCBarostatParams(pressure_atm=args.pressure, frequency=25)
    )
    if barostat:
        logger.info('NPT barostat enabled: P=%.2f atm', args.pressure)
    else:
        logger.info('Barostat disabled — running NVT.')

    masses = masses_from_species(species)
    velocities = maxwell_boltzmann_velocities(masses, args.temperature, rng)

    state = SimulationState(
        positions=positions,
        velocities=velocities,
        species=species,
        cell=cell,
        masses=masses,
    )

    # Pre-TDBB relaxation runs inside wf.run() from the config above (FIRE
    # minimize then unbiased equilibration), relaxing the classical-compressed
    # dense box before the first biased MLIP step. Nylon has no activation phase,
    # so this is the only pre-production stage (mirrors vinyl's non-activation
    # path; see workflows/polymerization.py run() lines 633-636).
    logger.info(
        'Pre-TDBB: minimize=%s (fmax=%.2f), equilibration=%d steps',
        config.minimize, config.minimize_fmax, config.equil_steps,
    )

    logger.info(
        'Starting TDBB: %d cycles × (%d biased + %d unbiased steps), T=%.0f K',
        config.n_cycles, config.biased_steps, config.unbiased_steps, args.temperature,
    )

    # Initial bond topology for trajectory + Carothers Fig. 4c (measured DPn vs
    # conversion) output. Best-effort — never fail the expensive run over
    # topology extraction. Seeds MUST match build_nylon66_system's placement:
    # diamines first (rdkit_seed=args.seed), then diacids (rdkit_seed=args.seed+1),
    # so the classical topology aligns 1:1 with coordinates/groups (RF23).
    init_bonds = None
    try:
        init_bonds = layout_bonds([
            (_DIAMINE_SMILES, args.n_diamines, args.seed),
            (_DIACID_SMILES, args.n_diacids, args.seed + 1),
        ])
    except Exception as exc:  # noqa: BLE001 — topology output is non-critical
        logger.warning('Bond-topology extraction failed (%s); trajectory will '
                       'carry no explicit bonds and Fig. 4c is unavailable.', exc)
    if args.mix and init_bonds is None:
        # Topology is optional for a plain run but MANDATORY for mixing, which
        # translates the live bond graph into a classical system. Fail here
        # rather than deep inside wf.run (decisions.md 2026-08-04 nylon 固有ガード).
        parser.error('--mix requires the initial bond topology, but its '
                     'extraction failed (see the warning above). Fix the '
                     'topology extraction or drop --mix.')

    wf = PolymerizationWorkflow(
        config, calc, template, groups,
        integrator=integrator,
        bond_tracker=tracker,
        barostat=barostat,
        initial_bonds=init_bonds,
    )

    # Reactive-site count for the alpha(t) denominator. Capture BEFORE the run
    # (the post-cycle updater mutates the groups) and pass it as n_monomers so the
    # trajectory header matches the figure command — no drift (RF22). Note: nylon
    # step-growth conversion is canonically the extent of reaction p (Carothers,
    # src/analysis/carothers.py); alpha(t) here is a secondary view.
    n_reactive_sites = (
        len(groups['amine_N'].atom_indices)
        + len(groups['carboxyl_C'].atom_indices)
    )
    # Checkpointing: on by default (write each cycle); --resume continues from the
    # last completed cycle. --no-checkpoint disables writing unless resuming.
    ckpt_file = args.output_dir / 'checkpoint.pkl'
    resuming = bool(args.resume and ckpt_file.exists())
    run_checkpoint_path = None if (args.no_checkpoint and not args.resume) else ckpt_file
    if args.resume and not ckpt_file.exists():
        logger.warning('--resume given but %s not found; starting a fresh run.', ckpt_file)

    _now_mix = mixing_setup_from_args(args)
    if resuming:
        # Guard the measurement mode across resume: silently switching mixing
        # on/off (or changing its duration) mid-run would corrupt the well-mixed
        # measurement without any recorded reason. The checkpoint records the
        # mixing setup; a mismatch with the current CLI args is a hard error.
        # (Older checkpoints predate this key: absent => the run had mixing off.)
        _ckpt_mix = (load_checkpoint(ckpt_file).get('extra', {}) or {}).get('mixing')
        if mixing_setup_mismatch(_ckpt_mix, _now_mix):
            parser.error(
                f'--mix settings differ from the checkpoint being resumed '
                f'(checkpoint: {_ckpt_mix}, now: {_now_mix}). Resume with the '
                f'same mixing configuration, or start a fresh run.')

    logs = wf.run(
        state,
        output_dir=args.output_dir,
        config_path='configs/boost/paper_faithful.yaml',
        n_monomers=n_reactive_sites,
        checkpoint_path=run_checkpoint_path,
        resume=resuming,
        # Record the mixing setup so resume can detect a mode switch (the guard
        # above compares this against the resume-time CLI args). Same
        # single-source-of-truth builder, so persisted and compared dicts can
        # never drift apart. Nylon has no spin state, so 'mixing' is the only key.
        checkpoint_extra={'mixing': _now_mix},
    )

    # A5: count one condensation per amide bond (amine_N-carboxyl_C). The paired
    # water-forming k-l event carries counts_as_reaction=False and is excluded,
    # so it is not double-counted toward Carothers p / alpha(t).
    all_formations = tracker.confirmed_formations()
    counted_formations = [e for e in all_formations if e.counts_as_reaction]
    n_form = len(counted_formations)
    n_form_all = len(all_formations)
    n_dissoc = len(tracker.confirmed_dissociations())
    logger.info(
        'Confirmed formations: %d counted (%d total incl. water-forming), '
        'dissociations: %d', n_form, n_form_all, n_dissoc,
    )
    summary = {
        'total_steps': state.step,
        'n_diamines': args.n_diamines,
        'n_diacids': args.n_diacids,
        'n_atoms': len(species),
        'box_size_A': initial_box_edge_A,
        'cell_periodic': True,
        'backend': calc.name,
        'compile': args.compile,
        'empty_cache': args.empty_cache,
        'temperature_K': args.temperature,
        'biased_steps': args.biased_steps,
        'unbiased_steps': args.unbiased_steps,
        'n_cycles': args.n_cycles,
        'minimize': args.minimize,
        'minimize_fmax': args.minimize_fmax,
        'equil_steps': args.equil_steps,
        'confirmed_formations': n_form,
        'confirmed_dissociations': n_dissoc,
        'n_reactive_sites': n_reactive_sites,
        'logs': [
            {
                'cycle': log.cycle,
                'phase': log.phase,
                'steps': log.steps,
                'n_candidates': log.n_candidates,
                'n_selected': log.n_selected,
                'bias_energy': log.bias_energy,
            }
            for log in logs
        ],
    }

    # Primary nylon metric is the Carothers extent of reaction p (decisions.md
    # 2026-07-06 / 2026-06-19 RF2:852): one amide bond per counted condensation.
    from kagome.analysis.carothers import dpn_from_bonds
    dpn = dpn_from_bonds(n_bonds=n_form, n_functional_groups=n_reactive_sites)
    p_extent = n_form / (n_reactive_sites / 2) if n_reactive_sites > 0 else 0.0
    summary['carothers_p'] = p_extent
    summary['carothers_dpn'] = dpn
    logger.info('Carothers: p=%.4f, DPn=%.3f (n_bonds=%d, n_groups=%d)',
                p_extent, dpn, n_form, n_reactive_sites)

    out_path = args.output_dir / 'summary.json'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(summary, indent=2), encoding='utf-8')
    logger.info('Done. Results in %s', args.output_dir)

    print('\nPrimary nylon metric is Carothers p (extent of reaction): '
          f'p={p_extent:.4f}, DPn={dpn:.3f}. See summary.json.')
    print(
        '\nTo generate figures (alpha(t)/Eq.11 is a secondary view; its '
        'denominator --n-reactive-sites is the count of reactive end groups, '
        f'amine_N + carboxyl_C = {n_reactive_sites}):'
    )
    print(
        f'  python scripts/reproduce_figures.py '
        f'--trajectory {args.output_dir}/trajectory.jsonl '
        f'--bonds {args.output_dir}/bonds.jsonl '
        f'--n-reactive-sites {n_reactive_sites} '
        f'--target-temperature {args.temperature} '
        f'--timestep-fs {args.timestep_fs} '
        f'--output-dir {args.output_dir}/figures'
    )


if __name__ == '__main__':
    main()
