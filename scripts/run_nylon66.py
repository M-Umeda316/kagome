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
"""
from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path

os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')

import numpy as np

from scripts._systems import (
    _DIACID_SMILES,
    _DIAMINE_SMILES,
    box_from_density,
    build_nylon66_system,
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
    masses_from_species,
)

logging.basicConfig(level=logging.INFO, format='%(name)s | %(message)s')
logger = logging.getLogger(__name__)


def _create_backend(backend: str, device: str, model: str) -> Calculator:
    if backend == 'orb':
        from kagome.backends.orb_backend import create_orb_calculator
        return create_orb_calculator(device=device)
    else:
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
    parser.add_argument('--box-size', type=float, default=None,
                        help='Box edge (Å). If omitted, computed from --density and '
                             'reached by direct placement or classical compression '
                             '(single-run path to paper density, mirroring run_vinyl_aibn).')
    parser.add_argument('--density', type=float, default=0.5,
                        help='Initial density (g/mL). Paper SI S-3 uses 0.5. Used only '
                             'when --box-size is omitted.')
    parser.add_argument('--temperature', type=float, default=300.0)
    parser.add_argument('--pressure', type=float, default=1.0)
    parser.add_argument('--no-barostat', action='store_true')
    parser.add_argument('--backend', type=str, default='orb',
                        choices=['orb', 'mace'])
    parser.add_argument('--device', type=str, default='cpu')
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
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    counts = {_DIAMINE_SMILES: args.n_diamines, _DIACID_SMILES: args.n_diacids}

    # Backend is created before the build so the 'ml' compress option can reuse it.
    calc = _create_backend(args.backend, args.device, args.model)
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

    langevin_params = LangevinParams(temperature_K=args.temperature)
    config = PolymerizationConfig(
        timestep_fs=0.25,
        biased_steps=args.biased_steps,
        unbiased_steps=args.unbiased_steps,
        n_cycles=args.n_cycles,
        tdbb=TDBBParams(
            f2=10.0,
            gamma=1.0,
            f1_max_formation=250.0,
            f1_max_dissociation=125.0,
            lambda_vdw=0.60,
        ),
        seed=args.seed,
        save_interval=50,
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
        'Starting TDBB: %d cycles × (%d biased + %d unbiased steps), T=%.0f K',
        config.n_cycles, config.biased_steps, config.unbiased_steps, args.temperature,
    )

    wf = PolymerizationWorkflow(
        config, calc, template, groups,
        integrator=integrator,
        bond_tracker=tracker,
        barostat=barostat,
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
    logs = wf.run(
        state,
        output_dir=args.output_dir,
        config_path='configs/boost/paper_faithful.yaml',
        n_monomers=n_reactive_sites,
    )

    n_form = len(tracker.confirmed_formations())
    n_dissoc = len(tracker.confirmed_dissociations())
    logger.info('Confirmed formations: %d, dissociations: %d', n_form, n_dissoc)
    summary = {
        'total_steps': state.step,
        'n_diamines': args.n_diamines,
        'n_diacids': args.n_diacids,
        'n_atoms': len(species),
        'box_size_A': initial_box_edge_A,
        'cell_periodic': True,
        'backend': calc.name,
        'temperature_K': args.temperature,
        'biased_steps': args.biased_steps,
        'unbiased_steps': args.unbiased_steps,
        'n_cycles': args.n_cycles,
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

    out_path = args.output_dir / 'summary.json'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(summary, indent=2), encoding='utf-8')
    logger.info('Done. Results in %s', args.output_dir)

    print('\nTo generate figures:')
    print(
        f'  python scripts/reproduce_figures.py '
        f'--trajectory {args.output_dir}/trajectory.jsonl '
        f'--bonds {args.output_dir}/bonds.jsonl '
        f'--n-reactive-sites {n_reactive_sites} '
        f'--target-temperature {args.temperature} '
        f'--output-dir {args.output_dir}/figures'
    )


if __name__ == '__main__':
    main()
