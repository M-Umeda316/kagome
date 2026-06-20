"""Nylon-6,6 step-growth polycondensation via TDBB.

System:  n_diamines hexamethylenediamine + n_diacids adipic acid
Backend: OrbMol-v2 (default) or MACE-MP-0 via --backend mace
Ensemble: NPT (Langevin + MC barostat) at 300 K, 1 atm

Paper anchor: arXiv:2511.22874, PDF p.22, Table S2, Fig. S2, Fig. 4.

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

from scripts._systems import build_nylon66_system
from src.backends.base import Calculator
from src.boost.tdbb import TDBBParams
from src.integrators.init_velocities import maxwell_boltzmann_velocities
from src.integrators.langevin import LangevinIntegrator, LangevinParams
from src.integrators.mc_barostat import MCBarostat, MCBarostatParams
from src.reactive.bonds import BondTracker
from src.workflows.polymerization import (
    PolymerizationConfig,
    PolymerizationWorkflow,
    SimulationState,
    masses_from_species,
)

logging.basicConfig(level=logging.INFO, format='%(name)s | %(message)s')
logger = logging.getLogger(__name__)


def _create_backend(backend: str, device: str, model: str) -> Calculator:
    if backend == 'orb':
        from src.backends.orb_backend import create_orb_calculator
        return create_orb_calculator(device=device)
    else:
        from src.backends.mace_backend import create_mace_calculator
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
    parser.add_argument('--box-size', type=float, default=25.0)
    parser.add_argument('--temperature', type=float, default=300.0)
    parser.add_argument('--pressure', type=float, default=1.0)
    parser.add_argument('--no-barostat', action='store_true')
    parser.add_argument('--backend', type=str, default='orb',
                        choices=['orb', 'mace'])
    parser.add_argument('--device', type=str, default='cpu')
    parser.add_argument('--model', type=str, default='small',
                        help='MACE model size (only used with --backend mace)')
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    logger.info(
        'Building nylon-6,6 system: %d diamines + %d diacids in %.1f Å box...',
        args.n_diamines, args.n_diacids, args.box_size,
    )
    positions, species, template, groups = build_nylon66_system(
        n_diamines=args.n_diamines,
        n_diacids=args.n_diacids,
        box_size=args.box_size,
        rng=rng,
    )
    logger.info(
        'System: %d atoms total  (%d amine_N, %d carboxyl_C, %d amine_H, %d carboxyl_OH)',
        len(species),
        len(groups['amine_N'].atom_indices),
        len(groups['carboxyl_C'].atom_indices),
        len(groups['amine_H'].atom_indices),
        len(groups['carboxyl_OH'].atom_indices),
    )

    cell = np.diag([args.box_size, args.box_size, args.box_size])

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

    calc = _create_backend(args.backend, args.device, args.model)
    logger.info('Backend: %s', calc.name)

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
        'box_size_A': args.box_size,
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
