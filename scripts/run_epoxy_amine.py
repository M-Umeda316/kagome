"""Bulk epoxy-amine curing (ring-opening polyaddition) via TDBB.

System:  n_epoxies DGEBA (2 epoxide rings) + n_amines DETA (2 primary + 1
         secondary N, 5 N-H total) — the paper's resin without the CuO slab
         (organic-only, decisions.md 2026-06-20 / 2026-07-08 E0 design).
Backend: OrbMol-v2 (default) or MACE-MP-0 via --backend mace
Ensemble: NPT (Langevin + MC barostat) at 333 K, 1 atm by default;
         --no-barostat for NVT (paper runs epoxy production in NVT, but on the
         excluded CuO slab — bulk follows the repo NPT densification convention).

Reaction: amine N attacks the epoxide terminal C -> N-C forms + ring C-O breaks
+ one N-H breaks; the H moves to the ring O (beta-hydroxyl). No leaving group
exits the molecule (contrast: nylon condensation releases water). Conjunctive
reaction events (F1', decisions.md 2026-07-08) drive ij AND ik AND jl together;
1° amines react twice (EpoxyAmineAdditionUpdater retires an N only when its
last H is consumed).

Paper anchor: arXiv:2511.22874 SI epoxy template + Table S3 (DGEBA/DETA),
p.24 (epoxy production at 333 K). E0 gate: decisions.md 2026-07-09 (PROCEED;
barrier 36.6 kcal/mol, product held, dE=-29.5 kcal/mol, f2=2 adopted).

Usage:
    python scripts/run_epoxy_amine.py --seed 7 --output-dir runs/epoxy_smoke
    python scripts/run_epoxy_amine.py --seed 7 --n-epoxies 100 --n-amines 50 \
        --device cuda --output-dir runs/epoxy_paper

    # well-mixed measurement mode (NOT paper-faithful): classical OpenMM/OpenFF
    # mixing after every cycle to refresh the reactive neighbourhood
    # (decisions.md 2026-08-04).
    python scripts/run_epoxy_amine.py --seed 7 --mix --mix-ps 25 \
        --mix-platform CUDA --output-dir runs/epoxy_mix
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
    _DETA_SMILES,
    _DGEBA_SMILES,
    box_from_density,
    build_epoxy_amine_system,
    layout_bonds,
)
from kagome.backends.base import Calculator
from kagome.boost.tdbb import TDBBParams
from kagome.integrators.init_velocities import maxwell_boltzmann_velocities
from kagome.integrators.langevin import LangevinIntegrator, LangevinParams
from kagome.integrators.mc_barostat import MCBarostat, MCBarostatParams
from kagome.reactive.bonds import BondTracker
from kagome.workflows.polymerization import (
    EpoxyAmineAdditionUpdater,
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
        description='Bulk epoxy-amine ring-opening curing (TDBB)'
    )
    parser.add_argument('--seed', type=int, default=7)
    parser.add_argument('--output-dir', type=Path, default=Path('runs/epoxy_amine'))
    parser.add_argument('--n-epoxies', type=int, default=10,
                        help='DGEBA molecules (2 epoxides each). Paper Table S3 '
                             'ratio is 100 DGEBA : 50 DETA.')
    parser.add_argument('--n-amines', type=int, default=5,
                        help='DETA molecules (5 N-H each).')
    parser.add_argument('--epoxy-smiles', type=str, default=_DGEBA_SMILES,
                        help='Epoxide monomer SMILES (default: DGEBA).')
    parser.add_argument('--amine-smiles', type=str, default=_DETA_SMILES,
                        help='Amine monomer SMILES (default: DETA).')
    parser.add_argument('--n-cycles', type=int, default=3)
    parser.add_argument('--biased-steps', type=int, default=2000)
    parser.add_argument('--unbiased-steps', type=int, default=1500)
    parser.add_argument('--f2', type=float, default=2.0,
                        help='TDBB Gaussian width f2 (Å⁻²). Paper default is 10.0, '
                             'but the E0 scan measured a capture-shell dead-zone at '
                             'f2=10 on OrbMol-v2 (max bias force ~0 in the [3,6] Å '
                             'candidate window) — f2=2 bridges it and is the adopted '
                             'default for epoxy-amine (decisions.md 2026-07-09, '
                             'user-approved deviation; same as the MA/nylon recipe).')
    # Pre-TDBB relaxation of the compressed dense structure — same rationale as
    # run_nylon66 (segfault fix, decisions.md 2026-07-08; paper anchor PDF p.20:
    # equilibration precedes production reactive MD). No activation phase, so:
    # build → compress → minimize → equilibrate → TDBB.
    parser.add_argument('--minimize', dest='minimize', action='store_true', default=True,
                        help='FIRE energy minimization before TDBB (default: on).')
    parser.add_argument('--no-minimize', dest='minimize', action='store_false',
                        help='Skip pre-TDBB energy minimization.')
    parser.add_argument('--minimize-fmax', type=float, default=1.0,
                        help='FIRE convergence threshold (kcal/mol/Å). Default 1.0.')
    parser.add_argument('--equil-steps', type=int, default=2000,
                        help='Unbiased equilibration steps before TDBB '
                             '(paper anchor PDF p.20; default 2000 as in nylon/vinyl).')
    parser.add_argument('--box-size', type=float, default=None,
                        help='Box edge (Å). If omitted, computed from --density and '
                             'reached by direct placement or classical compression.')
    parser.add_argument('--density', type=float, default=0.5,
                        help='Initial density (g/mL). SI S-3 convention 0.5, as in '
                             'nylon/vinyl. Used only when --box-size is omitted.')
    parser.add_argument('--temperature', type=float, default=333.0,
                        help='Production temperature (K). Paper p.24: epoxy curing '
                             'production at 333 K.')
    parser.add_argument('--friction-per-fs', type=float, default=0.01,
                        help='Langevin friction (1/fs). Paper SI value is 0.001 '
                             '(1.0 ps⁻¹), but the OrbMol f2=2 recipe needs 0.01 as '
                             'the cooling lever: with 0.001 the bias work + ring-'
                             'opening exotherm (-29.5 kcal/mol) accumulates and the '
                             '10+5 smoke overheated to mean 549 K at a 333 K target '
                             '(decisions.md 2026-07-07 friction / 2026-07-09 E1). '
                             'Pass 0.001 to restore the paper-faithful value.')
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
                        help='Calculator for box compression to target density '
                             '(classical = OpenMM/OpenFF Sage, keeps GPU free).')
    parser.add_argument('--compress-platform', type=str, default='CPU',
                        choices=['CPU', 'CUDA', 'OpenCL', 'Reference'])
    parser.add_argument('--resume', action='store_true', default=False,
                        help='Resume from <output-dir>/checkpoint.pkl if present.')
    parser.add_argument('--no-checkpoint', action='store_true', default=False,
                        help='Disable per-cycle checkpoint writing.')
    # WM-P3 mixing stage (--mix + 6 knobs), shared with run_vinyl_copolymer /
    # run_nylon66 via scripts/_mixing_cli.py (decisions.md 2026-08-04).
    add_mixing_arguments(parser)
    args = parser.parse_args()

    resolve_mixing_args(parser, args)

    rng = np.random.default_rng(args.seed)

    counts = {args.epoxy_smiles: args.n_epoxies, args.amine_smiles: args.n_amines}

    calc = _create_backend(args.backend, args.device, args.model,
                           compile_model=args.compile,
                           empty_cache=args.empty_cache)
    logger.info('Backend: %s', calc.name)

    def _build(edge: float, gen):
        return build_epoxy_amine_system(
            n_epoxies=args.n_epoxies,
            n_amines=args.n_amines,
            box_size=edge,
            rng=gen,
            epoxy_smiles=args.epoxy_smiles,
            amine_smiles=args.amine_smiles,
            rdkit_seed=args.seed,  # RF23: conformer geometry follows --seed
        )

    if args.box_size is not None:
        logger.info(
            'Building epoxy-amine: %d epoxies + %d amines in %.1f Å box (explicit)...',
            args.n_epoxies, args.n_amines, args.box_size,
        )
        positions, species, template, groups, amine_h_map = _build(args.box_size, rng)
        cell = np.diag([args.box_size, args.box_size, args.box_size])
    else:
        target_edge = box_from_density(counts, args.density)
        logger.info(
            'Building epoxy-amine: %d epoxies + %d amines, target box %.1f Å '
            '(%.2f g/mL)...',
            args.n_epoxies, args.n_amines, target_edge, args.density,
        )
        try:
            positions, species, template, groups, amine_h_map = _build(
                target_edge, rng)
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
                    positions, species, template, groups, amine_h_map = _build(
                        edge, np.random.default_rng(args.seed))
                    place_edge = edge
                    logger.info(
                        'Placed at dilute density %.2f g/mL (box %.2f Å); '
                        'compressing to %.2f Å.',
                        place_density, edge, target_edge,
                    )
                    break
                except RuntimeError:
                    continue
            if place_edge is None:
                raise RuntimeError(
                    'Could not place epoxy-amine even at dilute density 0.10 g/mL.'
                )
            from kagome.backends.classical_backend import make_compress_calculator
            from kagome.integrators.minimize import compress_box
            from kagome.prep.openmm_equilibrate import MoleculeSpec
            # Placement order in build_epoxy_amine_system: epoxies first (seed),
            # then amines (seed+1). MoleculeSpec order/seeds must match (RF23).
            specs = [
                MoleculeSpec(args.epoxy_smiles, args.n_epoxies, rdkit_seed=args.seed),
                MoleculeSpec(args.amine_smiles, args.n_amines, rdkit_seed=args.seed + 1),
            ]
            compress_calc = make_compress_calculator(
                args.compress_backend, specs, calc, platform=args.compress_platform,
                target_edge_A=target_edge,
            )
            place_cell = np.diag([place_edge, place_edge, place_edge])
            result = compress_box(positions, place_cell, target_edge, species, compress_calc)
            positions, cell = result.positions, result.cell

    n_epoxide_sites = len(groups['epoxy_C'].atom_indices)
    n_amine_h_sites = len(groups['amine_H'].atom_indices)
    logger.info(
        'System: %d atoms total  (%d amine_N, %d epoxy_C, %d amine_H, %d ring_O)',
        len(species),
        len(groups['amine_N'].atom_indices),
        n_epoxide_sites,
        n_amine_h_sites,
        len(groups['ring_O'].atom_indices),
    )

    initial_box_edge_A = float(cell[0, 0])

    langevin_params = LangevinParams(
        temperature_K=args.temperature, friction_per_fs=args.friction_per_fs)
    config = PolymerizationConfig(
        timestep_fs=0.25,
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

    logger.info(
        'Pre-TDBB: minimize=%s (fmax=%.2f), equilibration=%d steps',
        config.minimize, config.minimize_fmax, config.equil_steps,
    )
    logger.info(
        'Starting TDBB: %d cycles × (%d biased + %d unbiased steps), T=%.0f K',
        config.n_cycles, config.biased_steps, config.unbiased_steps, args.temperature,
    )

    # Initial bond topology for trajectory / network analysis output.
    # Best-effort — never fail the expensive run over topology extraction.
    init_bonds = None
    try:
        init_bonds = layout_bonds([
            (args.epoxy_smiles, args.n_epoxies, args.seed),
            (args.amine_smiles, args.n_amines, args.seed + 1),
        ])
    except Exception as exc:  # noqa: BLE001 — topology output is non-critical
        logger.warning('Bond-topology extraction failed (%s); trajectory will '
                       'carry no explicit bonds.', exc)
    if args.mix and init_bonds is None:
        # Topology is optional for a plain run but MANDATORY for mixing, which
        # translates the live bond graph into a classical system. Fail here
        # rather than deep inside wf.run (decisions.md 2026-08-04 nylon 固有ガード;
        # run_epoxy_amine shares the same best-effort extraction).
        parser.error('--mix requires the initial bond topology, but its '
                     'extraction failed (see the warning above). Fix the '
                     'topology extraction or drop --mix.')

    # 1° -> 2° -> 3° amine reassignment: an N stays selectable until its last
    # registered H is consumed (EpoxyAmineAdditionUpdater docstring).
    updater = EpoxyAmineAdditionUpdater(amine_h_map)

    wf = PolymerizationWorkflow(
        config, calc, template, groups,
        integrator=integrator,
        bond_tracker=tracker,
        barostat=barostat,
        initial_bonds=init_bonds,
        updater=updater,
    )

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
        n_monomers=n_epoxide_sites,
        checkpoint_path=run_checkpoint_path,
        resume=resuming,
        # Record the mixing setup so resume can detect a mode switch (the guard
        # above compares this against the resume-time CLI args). Same
        # single-source-of-truth builder, so persisted and compared dicts can
        # never drift apart. Epoxy has no spin state, so 'mixing' is the only key.
        checkpoint_extra={'mixing': _now_mix},
    )

    # One counted N-C formation per opened epoxide; the hydroxyl O-H event
    # carries counts_as_reaction=False and is excluded (no double counting).
    all_formations = tracker.confirmed_formations()
    counted_formations = [e for e in all_formations if e.counts_as_reaction]
    n_form = len(counted_formations)
    n_form_all = len(all_formations)
    n_dissoc = len(tracker.confirmed_dissociations())
    logger.info(
        'Confirmed formations: %d counted (%d total incl. hydroxyl O-H), '
        'dissociations: %d', n_form, n_form_all, n_dissoc,
    )

    # Primary curing metric: epoxide conversion (fraction of rings opened).
    epoxide_conversion = n_form / n_epoxide_sites if n_epoxide_sites > 0 else 0.0
    amine_h_conversion = n_form / n_amine_h_sites if n_amine_h_sites > 0 else 0.0
    logger.info('Epoxide conversion: %.4f (%d/%d rings), amine-H conversion: %.4f',
                epoxide_conversion, n_form, n_epoxide_sites, amine_h_conversion)

    summary = {
        'total_steps': state.step,
        'n_epoxies': args.n_epoxies,
        'n_amines': args.n_amines,
        'epoxy_smiles': args.epoxy_smiles,
        'amine_smiles': args.amine_smiles,
        'n_atoms': len(species),
        'box_size_A': initial_box_edge_A,
        'cell_periodic': True,
        'backend': calc.name,
        'compile': args.compile,
        'empty_cache': args.empty_cache,
        'temperature_K': args.temperature,
        'friction_per_fs': args.friction_per_fs,
        'biased_steps': args.biased_steps,
        'unbiased_steps': args.unbiased_steps,
        'n_cycles': args.n_cycles,
        'f2': args.f2,
        'minimize': args.minimize,
        'minimize_fmax': args.minimize_fmax,
        'equil_steps': args.equil_steps,
        'confirmed_formations': n_form,
        'confirmed_dissociations': n_dissoc,
        'n_epoxide_sites': n_epoxide_sites,
        'n_amine_h_sites': n_amine_h_sites,
        'epoxide_conversion': epoxide_conversion,
        'amine_h_conversion': amine_h_conversion,
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

    out_path = args.output_dir / 'summary.json'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(summary, indent=2), encoding='utf-8')
    logger.info('Done. Results in %s', args.output_dir)

    print('\nPrimary curing metric is epoxide conversion: '
          f'{epoxide_conversion:.4f} ({n_form}/{n_epoxide_sites} rings opened). '
          'See summary.json.')
    print(
        '\nTo generate figures (alpha(t) denominator --n-reactive-sites is the '
        f'epoxide site count = {n_epoxide_sites}):'
    )
    print(
        f'  python scripts/reproduce_figures.py '
        f'--trajectory {args.output_dir}/trajectory.jsonl '
        f'--bonds {args.output_dir}/bonds.jsonl '
        f'--n-reactive-sites {n_epoxide_sites} '
        f'--target-temperature {args.temperature} '
        f'--timestep-fs {config.timestep_fs} '
        f'--output-dir {args.output_dir}/figures'
    )


if __name__ == '__main__':
    main()
