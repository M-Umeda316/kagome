"""Radical vinyl polymerization with methyl acrylate + AIBN initiator.

System:  n_monomers methyl acrylate (C=CC(=O)OC) +
         n_initiators isobutyronitrile radical model (CC(C)C#N)
Backend: OrbMol-v2 (default) or MACE-MP-0 via --backend mace
Ensemble: NPT (Langevin + MC barostat) or NVT with --no-barostat

Design:  specs/decisions.md — "T-G1: vinyl radical polymerization system"
Paper:   arXiv:2511.22874, Section 3, Table S1

Usage:
    python scripts/run_vinyl_aibn.py --seed 7 --output-dir runs/vinyl_aibn
    python scripts/run_vinyl_aibn.py --seed 7 --backend mace --output-dir runs/vinyl_aibn_mace
    python scripts/run_vinyl_aibn.py --seed 7 --no-barostat --output-dir runs/vinyl_aibn_nvt
"""
from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path

os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')
# Long paper-scale runs creep up in VRAM because the OrbMol-v2 neighbour graph
# size varies per step (atoms move), fragmenting the CUDA caching allocator until
# it exhausts memory and the run hangs. Expandable segments defragment this at no
# per-step cost. Must be set before torch is imported (orb import happens lazily
# in _create_backend). See specs/decisions.md 2026-06-15 VRAM record.
os.environ.setdefault('PYTORCH_CUDA_ALLOC_CONF', 'expandable_segments:True')

import numpy as np

from scripts._systems import build_vinyl_aibn_system, build_full_aibn_system, build_activation_template
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


def _create_backend(backend: str, device: str, model: str, spin: int = 1) -> Calculator:
    if backend == 'orb':
        from kagome.backends.orb_backend import create_orb_calculator
        return create_orb_calculator(device=device, spin=spin)
    elif backend == 'aimnet':
        from kagome.backends.aimnet_backend import create_aimnet_calculator
        # model carries the AIMNet2 model name when backend=aimnet (default below)
        aimnet_model = model if model and model.startswith('aimnet') else 'aimnet2-nse'
        return create_aimnet_calculator(model=aimnet_model, device=device, spin=spin)
    else:
        from kagome.backends.mace_backend import create_mace_calculator
        return create_mace_calculator(model=model, device=device)


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Radical vinyl polymerization: methyl acrylate + AIBN'
    )
    parser.add_argument('--seed', type=int, default=7)
    parser.add_argument('--output-dir', type=Path, default=Path('runs/vinyl_aibn'))
    parser.add_argument('--n-monomers', type=int, default=8)
    parser.add_argument('--n-initiators', type=int, default=2)
    parser.add_argument('--n-cycles', type=int, default=3)
    # Defaults are exploratory (fast). Paper-faithful is 2000/2000 (PDF p.7);
    # the S6 production scripts pass 2000 biased / 500 unbiased (decisions.md
    # S2-S4, validated for OrbMol-v2). Effective values are recorded per run in
    # manifest.extra, so the default deviation is auditable (RF22).
    parser.add_argument('--biased-steps', type=int, default=500,
                        help='exploratory default 500; paper-faithful 2000 (PDF p.7)')
    parser.add_argument('--unbiased-steps', type=int, default=500,
                        help='exploratory default 500; paper 2000, S6 production 500')
    parser.add_argument('--box-size', type=float, default=None,
                        help='Box edge (Å). If omitted, computed from --density.')
    parser.add_argument('--density', type=float, default=0.5,
                        help='Initial density (g/mL). Paper SI S-3 uses 0.5 for vinyl. '
                             'Used only when --box-size is omitted.')
    parser.add_argument('--temperature', type=float, default=333.0)
    parser.add_argument('--friction-per-fs', type=float, default=0.001,
                        help='Langevin friction coefficient (1/fs). Default 0.001 '
                             '(weak, τ~1 ps). Increase (e.g. 0.01) to dissipate TDBB '
                             'bias work + exothermic reaction heat faster and keep the '
                             'reactive melt near the target temperature across cycles.')
    parser.add_argument('--pressure', type=float, default=1.0,
                        help='Target pressure (atm). Default 1.0 (assumed, not stated in paper).')
    parser.add_argument('--no-barostat', action='store_true',
                        help='Disable NPT barostat and run NVT instead.')
    parser.add_argument('--backend', type=str, default='orb',
                        choices=['orb', 'mace', 'aimnet'],
                        help='MLIP backend (default: orb = OrbMol-v2; aimnet = AIMNet2-NSE)')
    parser.add_argument('--device', type=str, default='cpu')
    parser.add_argument('--model', type=str, default='small',
                        help='MACE model size (only used with --backend mace)')
    parser.add_argument('--compress-backend', type=str, default='classical',
                        choices=['classical', 'ml'],
                        help='Calculator for box compression to paper density. '
                             '"classical" (default) uses OpenMM/OpenFF Sage so densification '
                             'does not consume MLIP GPU time (decision 2026-06-20); "ml" uses '
                             'the production MLIP. MD always runs on the MLIP.')
    parser.add_argument('--compress-platform', type=str, default='CPU',
                        choices=['CPU', 'CUDA', 'OpenCL', 'Reference'],
                        help='OpenMM platform for --compress-backend classical '
                             '(default CPU, keeps the GPU free for the MLIP MD).')
    parser.add_argument('--initiator-smiles', type=str, default=None,
                        help='Override the initiator SMILES (e.g. "C[C](C)C#N" for the '
                             'real open-shell 2-cyanoprop-2-yl radical). Default: closed-shell model.')
    parser.add_argument('--spin', type=int, default=1,
                        help='Total spin multiplicity (2S+1) passed to OrbMol-v2. '
                             'Use 2 (doublet) for a single radical. Default 1 (singlet).')
    parser.add_argument('--production-spin-cap', type=int, default=None,
                        help='Cap the post-activation production multiplicity at this '
                             'value instead of n_radicals+1 (high-spin sum). Diagnostic '
                             'for isolating whether high multiplicity destabilises the '
                             'dynamics. Default: no cap (high-spin sum).')
    parser.add_argument('--minimize', dest='minimize', action='store_true', default=True,
                        help='FIRE energy minimization before TDBB (default: on). '
                             'Relaxes initial close contacts (paper anchor PDF p.20).')
    parser.add_argument('--no-minimize', dest='minimize', action='store_false',
                        help='Skip pre-TDBB energy minimization.')
    parser.add_argument('--minimize-fmax', type=float, default=1.0,
                        help='FIRE convergence threshold (kcal/mol/Å). Default 1.0.')
    parser.add_argument('--equil-steps', type=int, default=2000,
                        help='Unbiased NPT equilibration steps before TDBB '
                             '(paper anchor PDF p.20; length not specified, default 2000 '
                             '= 500 fs matching a TDBB block). 0 disables.')
    parser.add_argument('--timestep-fs', type=float, default=0.25,
                        help='MD timestep (fs). Default 0.25 fs (conservative, validated '
                             'for FIRE densification + ML NVT). 1.0 fs is standard for '
                             'organic ML MD and gives 4x speed for the same physical time.')
    parser.add_argument('--f2', type=float, default=10.0,
                        help='TDBB Gaussian width parameter f2 (Å⁻²). Paper default 10.0, '
                             'stated robust range 5-20. Lower values widen the bias well '
                             '(capture radius ~1/√f2). Use 5.0 for OrbMol-v2 PES '
                             '(see decisions.md 2026-06-17).')
    parser.add_argument('--select-rmin', type=float, default=None,
                        help='Override candidate selection r_min (Å). Paper Table S1: 3.0. '
                             'For OrbMol-v2 PES-tuned window use 1.5 (see decisions.md 2026-06-17).')
    parser.add_argument('--select-rmax', type=float, default=None,
                        help='Override candidate selection r_max (Å). Paper Table S1: 6.0. '
                             'For OrbMol-v2 PES-tuned window use 3.0 (see decisions.md 2026-06-17).')
    parser.add_argument('--activation', action='store_true', default=False,
                        help='Use full AIBN molecules and run V^d activation phase to '
                             'decompose C-N azo bonds before propagation. Starts with '
                             'spin=1, switches to spin=N_radicals+1 after activation. '
                             'Paper anchor: Table S1 Activation row.')
    parser.add_argument('--activation-steps', type=int, default=3000,
                        help='Max steps for the activation biased phase (default 3000).')
    parser.add_argument('--activation-f2', type=float, default=0.3,
                        help='V^d Gaussian width for activation dissociation (Å⁻²). '
                             'f2=0.3 puts force peak at ~1.29 Å, near C-N bond distance. '
                             'Default 0.3.')
    parser.add_argument('--activation-f1-max', type=float, default=250.0,
                        help='Peak V^d amplitude for activation (kcal/mol). Must exceed '
                             '~200 with f2=0.3 for OrbMol-v2 C-N barrier (~39 kcal/mol). '
                             'Default 250.')
    parser.add_argument('--resume', action='store_true', default=False,
                        help='Resume from <output-dir>/checkpoint.pkl if present: skip '
                             'build-time activation and start at the saved cycle. Bit-exact '
                             'continuation (rng state restored). Use after a long run was '
                             'killed (e.g. crash/swap). No-op if no checkpoint exists.')
    parser.add_argument('--no-checkpoint', action='store_true', default=False,
                        help='Disable writing <output-dir>/checkpoint.pkl each cycle. '
                             'By default a checkpoint is saved at every cycle boundary so '
                             'a killed run can --resume.')
    parser.add_argument('--load-structure', type=Path, default=None,
                        help='Optional: load a classically pre-equilibrated structure (JSON '
                             'from scripts/prep_structure.py) and skip build/place/compress. '
                             'Not required for paper density — by default the run compresses '
                             'in-process (--compress-backend classical). Use this only to '
                             'reuse one prepped structure across many seeds. positions+cell '
                             'come from the file; the short ML re-equil (--minimize/'
                             '--equil-steps) still runs. See decision D-4 (+2026-06-20 WSL).')
    args = parser.parse_args()

    # Cycle-boundary checkpointing for crash recovery (long runs). By default a
    # checkpoint is written every cycle; --resume continues from it (skipping the
    # one-time build-time activation/minimize/equil). See decisions.md 2026-06-26.
    ckpt_file = args.output_dir / 'checkpoint.pkl'
    resuming = bool(args.resume and ckpt_file.exists())
    run_checkpoint_path = None if (args.no_checkpoint and not args.resume) else ckpt_file
    if args.resume and not ckpt_file.exists():
        logger.warning('--resume given but %s not found; starting a fresh run.', ckpt_file)

    rng = np.random.default_rng(args.seed)

    from scripts._systems import (
        _AIBN_SMILES,
        _INITIATOR_SMILES,
        _MONOMER_SMILES,
        box_from_density,
    )

    if args.activation:
        init_smiles_for_density = _AIBN_SMILES
        initial_spin = 1
    else:
        init_smiles_for_density = args.initiator_smiles or _INITIATOR_SMILES
        initial_spin = args.spin

    counts = {_MONOMER_SMILES: args.n_monomers, init_smiles_for_density: args.n_initiators}

    if args.box_size is not None:
        target_edge = args.box_size
    else:
        target_edge = box_from_density(counts, args.density)
        logger.info(
            'Box edge from density %.2f g/mL: %.2f Å (paper SI S-3)',
            args.density, target_edge,
        )

    calc = _create_backend(args.backend, args.device, args.model, spin=initial_spin)
    logger.info('Backend: %s (spin=%d)', calc.name, initial_spin)

    _init_smiles = args.initiator_smiles or _INITIATOR_SMILES

    aibn_azo_bonds = None

    if args.activation:
        def _build_aibn(edge: float, gen: np.random.Generator):
            return build_full_aibn_system(
                n_monomers=args.n_monomers,
                n_aibn=args.n_initiators,
                box_size=edge,
                rng=gen,
                rdkit_seed=args.seed,  # RF23: conformer geometry follows --seed
            )

        logger.info(
            'Building full AIBN system: %d monomers + %d AIBN, target box %.1f Å...',
            args.n_monomers, args.n_initiators, target_edge,
        )
        try:
            positions, species, aibn_azo_bonds, template, groups, propagation_map, chain_c_map = _build_aibn(target_edge, rng)
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
                    positions, species, aibn_azo_bonds, template, groups, propagation_map, chain_c_map = _build_aibn(edge, rng)
                    place_edge = edge
                    break
                except RuntimeError:
                    continue
            if place_edge is None:
                raise RuntimeError('Could not place the system even at dilute density 0.10 g/mL.')
            from kagome.integrators.minimize import compress_box
            from kagome.backends.classical_backend import make_compress_calculator
            from kagome.prep.openmm_equilibrate import MoleculeSpec
            # Placement order in build_full_aibn_system: AIBN first (seed), then
            # monomers (seed+1). MoleculeSpec order/seeds must match (RF23).
            specs = [
                MoleculeSpec(_AIBN_SMILES, args.n_initiators, rdkit_seed=args.seed),
                MoleculeSpec(_MONOMER_SMILES, args.n_monomers, rdkit_seed=args.seed + 1),
            ]
            compress_calc = make_compress_calculator(
                args.compress_backend, specs, calc, platform=args.compress_platform,
                target_edge_A=target_edge,
            )
            place_cell = np.diag([place_edge, place_edge, place_edge])
            result = compress_box(positions, place_cell, target_edge, species, compress_calc)
            positions, cell = result.positions, result.cell
    else:
        def _build(edge: float, gen: np.random.Generator):
            return build_vinyl_aibn_system(
                n_monomers=args.n_monomers,
                n_initiators=args.n_initiators,
                box_size=edge,
                rng=gen,
                initiator_smiles=_init_smiles,
                rdkit_seed=args.seed,  # RF23: conformer geometry follows --seed
            )

        if args.load_structure is not None:
            from kagome.prep.structure_io import PreparedStructure

            prepared = PreparedStructure.load(args.load_structure)
            meta_edge = box_from_density(counts, 0.10)
            _, ref_species, template, groups, propagation_map, chain_c_map = _build(meta_edge, rng)
            if list(ref_species) != list(prepared.species):
                raise ValueError(
                    'Loaded structure species do not match the builder '
                    f'(loaded N={len(prepared.species)}, builder N={len(ref_species)}). '
                    'Ensure --n-monomers/--n-initiators match the prepped structure.'
                )
            positions = prepared.positions
            species = prepared.species
            cell = (prepared.cell if prepared.cell is not None
                    else np.diag([target_edge, target_edge, target_edge]))
            logger.info(
                'Loaded pre-equilibrated structure from %s (%d atoms, box %.2f Å).',
                args.load_structure, len(species), float(cell[0, 0]),
            )
        else:
            logger.info(
                'Building vinyl/AIBN system: %d monomers + %d initiators, target box %.1f Å...',
                args.n_monomers, args.n_initiators, target_edge,
            )
            try:
                positions, species, template, groups, propagation_map, chain_c_map = _build(target_edge, rng)
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
                        positions, species, template, groups, propagation_map, chain_c_map = _build(edge, rng)
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
                        'Could not place the system even at dilute density 0.10 g/mL.'
                    )
                from kagome.integrators.minimize import compress_box
                from kagome.backends.classical_backend import make_compress_calculator
                from kagome.prep.openmm_equilibrate import MoleculeSpec
                # Placement order in build_vinyl_aibn_system: initiators first
                # (seed), then monomers (seed+1). MoleculeSpec order/seeds must
                # match so the classical topology aligns with coords (RF23).
                specs = [
                    MoleculeSpec(_init_smiles, args.n_initiators, rdkit_seed=args.seed),
                    MoleculeSpec(_MONOMER_SMILES, args.n_monomers, rdkit_seed=args.seed + 1),
                ]
                compress_calc = make_compress_calculator(
                    args.compress_backend, specs, calc, platform=args.compress_platform,
                    target_edge_A=target_edge,
                )
                place_cell = np.diag([place_edge, place_edge, place_edge])
                result = compress_box(positions, place_cell, target_edge, species, compress_calc)
                positions, cell = result.positions, result.cell

    if args.select_rmin is not None or args.select_rmax is not None:
        for ps in template.pairs:
            if args.select_rmin is not None:
                ps.r_min = args.select_rmin
            if args.select_rmax is not None:
                ps.r_max = args.select_rmax
        logger.info(
            'Candidate window overridden: [%.1f, %.1f] Å (paper Table S1: [3.0, 6.0])',
            template.pairs[0].r_min, template.pairs[0].r_max,
        )

    logger.info(
        'System: %d atoms total  (%d radical_C, %d vinyl_alpha_C sites), box %.2f Å',
        len(species),
        len(groups['radical_C'].atom_indices),
        len(groups['vinyl_alpha_C'].atom_indices),
        float(cell[0, 0]),
    )
    logger.info('Propagation map: %d entries', len(propagation_map))

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

    wf = PolymerizationWorkflow(
        config, calc, template, groups,
        integrator=integrator,
        bond_tracker=tracker,
        barostat=barostat,
        propagation_map=propagation_map,
        propagation_target_group='radical_C',
        chain_c_map=chain_c_map,
    )

    n_activation_dissoc = 0
    if resuming:
        # The expensive one-time prep (build/compress/minimize/activation/equil)
        # already happened before the checkpoint was written. Restore the
        # post-activation production spin and skip straight to the cycle loop;
        # run() restores positions/groups/tracker/rng from the checkpoint.
        _ck = load_checkpoint(ckpt_file)
        _extra = _ck.get('extra', {}) or {}
        n_activation_dissoc = int(_extra.get('activation_dissociations', 0))
        _spin = _extra.get('spin')
        if _spin is not None and getattr(calc, 'supports_spin', False):
            calc.set_spin(int(_spin))
            logger.info('Resume: restored production spin %d from checkpoint', int(_spin))
        config = PolymerizationConfig(
            timestep_fs=config.timestep_fs,
            biased_steps=config.biased_steps,
            unbiased_steps=config.unbiased_steps,
            n_cycles=config.n_cycles,
            tdbb=config.tdbb,
            seed=config.seed,
            save_interval=config.save_interval,
            minimize=False,
            minimize_fmax=config.minimize_fmax,
            equil_steps=0,
        )
        wf.config = config
        logger.info('Resume mode: skipping build-time activation/minimize/equilibration.')
    elif args.activation and aibn_azo_bonds:
        logger.info(
            'Activation phase: %d C-N azo bonds, %d steps, f2=%.2f, f1_max=%.0f, spin=1',
            len(aibn_azo_bonds), args.activation_steps, args.activation_f2,
            args.activation_f1_max,
        )
        act_template, act_groups = build_activation_template(aibn_azo_bonds)

        # Order matters (see specs/decisions.md 2026-06-24 "活性化と平衡化の順序バグ").
        # minimize (FIRE, 0 K) relaxes placement clashes WITHOUT thermally
        # decomposing the AIBN; activation then dissociates the azo C-N bonds on
        # the still-intact structure. The 333 K equilibration is DEFERRED to
        # after activation — running it before (the old order) thermally
        # decomposed/scattered the azo bonds (C-N stretched past the [0,3] Å
        # activation window), so activation found 0 C-N candidates, no real
        # radicals formed, and the radical_C/chain_C topology was broken,
        # starving downstream candidate selection (qualified candidates 0 vs 16).
        if config.minimize:
            wf._minimize(state)

        dissociated = wf.run_activation(
            state, act_template, act_groups,
            activation_steps=args.activation_steps,
            activation_f2=args.activation_f2,
            activation_f1_max=args.activation_f1_max,
            rng=np.random.default_rng(args.seed + 1),
        )
        n_activation_dissoc = len(dissociated)
        logger.info('Activation result: %d C-N bonds dissociated', n_activation_dissoc)

        if dissociated:
            n_radicals = len(groups['radical_C'].atom_indices)
            production_spin = n_radicals + 1
            if args.production_spin_cap is not None:
                production_spin = min(production_spin, args.production_spin_cap)
                logger.info('Production spin capped at %d (from n_radicals+1=%d)',
                            production_spin, n_radicals + 1)
            calc.set_spin(production_spin)
            if getattr(calc, 'supports_spin', False):
                logger.info('Spin switched: 1 → %d (N_radicals=%d)', production_spin, n_radicals)
            else:
                logger.warning(
                    'Backend %s ignores spin; production spin %d NOT applied '
                    '(N_radicals=%d) — results assume the backend is spin-agnostic.',
                    calc.name, production_spin, n_radicals,
                )

        # Thermalize AFTER activation: radicals are already formed, so there is
        # no intact azo left for the 333 K dynamics to decompose.
        if config.equil_steps > 0:
            act_rng = np.random.default_rng(args.seed)
            wf._run_equilibration_phase(state, act_rng, writer=None)

        config = PolymerizationConfig(
            timestep_fs=config.timestep_fs,
            biased_steps=config.biased_steps,
            unbiased_steps=config.unbiased_steps,
            n_cycles=config.n_cycles,
            tdbb=config.tdbb,
            seed=config.seed,
            save_interval=config.save_interval,
            minimize=False,
            minimize_fmax=config.minimize_fmax,
            equil_steps=0,
        )
        wf.config = config

    logger.info(
        'Starting TDBB: %d cycles × (%d biased + %d unbiased steps), T=%.0f K',
        config.n_cycles, config.biased_steps, config.unbiased_steps, args.temperature,
    )

    logs = wf.run(
        state,
        output_dir=args.output_dir,
        config_path='configs/boost/paper_faithful.yaml',
        n_monomers=args.n_monomers,
        checkpoint_path=run_checkpoint_path,
        resume=resuming,
        checkpoint_extra={
            'spin': getattr(calc, '_spin', None),
            'activation_dissociations': n_activation_dissoc,
        },
    )

    n_form = len(tracker.confirmed_formations())
    n_dissoc = len(tracker.confirmed_dissociations())
    logger.info('Confirmed formations: %d, dissociations: %d', n_form, n_dissoc)

    summary = {
        'total_steps': state.step,
        'n_monomers': args.n_monomers,
        'n_initiators': args.n_initiators,
        'n_atoms': len(species),
        'box_size_A': float(cell[0, 0]),
        'cell_periodic': True,
        'backend': calc.name,
        'temperature_K': langevin_params.temperature_K,
        'biased_steps': args.biased_steps,
        'unbiased_steps': args.unbiased_steps,
        'n_cycles': args.n_cycles,
        'minimize': args.minimize,
        'minimize_fmax': args.minimize_fmax,
        'equil_steps': args.equil_steps,
        'activation': args.activation,
        'activation_steps': args.activation_steps if args.activation else 0,
        'activation_f2': args.activation_f2 if args.activation else None,
        'activation_f1_max': args.activation_f1_max if args.activation else None,
        'activation_dissociations': n_activation_dissoc,
        'confirmed_formations': n_form,
        'confirmed_dissociations': n_dissoc,
        'propagation_events': n_form,
        'logs': [
            {
                'cycle': log.cycle,
                'phase': log.phase,
                'steps': log.steps,
                'n_candidates': log.n_candidates,
                'n_selected': log.n_selected,
                'bias_energy': log.bias_energy,
                'min_pair_distance': (None if log.min_pair_distance == float('inf')
                                      else log.min_pair_distance),
            }
            for log in logs
        ],
    }

    out_path = args.output_dir / 'summary.json'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(summary, indent=2), encoding='utf-8')
    logger.info('Done. Results in %s', args.output_dir)

    print('\nTo generate figures:')
    print(
        f'  python scripts/reproduce_figures.py '
        f'--trajectory {args.output_dir}/trajectory.jsonl '
        f'--bonds {args.output_dir}/bonds.jsonl '
        f'--n-reactive-sites {args.n_monomers} '
        f'--target-temperature {args.temperature} '
        f'--output-dir {args.output_dir}/figures'
    )


if __name__ == '__main__':
    main()
