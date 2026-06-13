"""TDBB polymerization with MACE-MP-0 backend and periodic boundary conditions.

T7.1: Verify that MACE-MP-0 + PBC works end-to-end.
T8.2: Paper-faithful scale with biased/unbiased 2000 steps.

Usage:
    python scripts/run_mace_pbc.py --seed 7 --output-dir runs/mace_pbc
    python scripts/run_mace_pbc.py --seed 7 --output-dir runs/mace_pbc_paper --biased-steps 2000 --unbiased-steps 2000
"""
from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path

os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')

import numpy as np

from scripts._systems import build_ethylene_box, build_template_and_groups
from src.backends.mace_backend import create_mace_calculator
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


def main() -> None:
    parser = argparse.ArgumentParser(description='MACE-MP-0 + PBC TDBB polymerization')
    parser.add_argument('--seed', type=int, default=7)
    parser.add_argument('--output-dir', type=Path, default=Path('runs/mace_pbc'))
    parser.add_argument('--n-molecules', type=int, default=4)
    parser.add_argument('--n-cycles', type=int, default=3)
    parser.add_argument('--biased-steps', type=int, default=500)
    parser.add_argument('--unbiased-steps', type=int, default=500)
    parser.add_argument('--box-size', type=float, default=10.0)
    parser.add_argument('--temperature', type=float, default=500.0)
    parser.add_argument('--pressure', type=float, default=1.0,
                        help='Target pressure (atm). Default 1.0 atm (assumed, not stated in paper).')
    parser.add_argument('--no-barostat', action='store_true',
                        help='Disable NPT barostat and run NVT instead.')
    parser.add_argument('--device', type=str, default='cpu')
    parser.add_argument('--model', type=str, default='small')
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    logger.info('Building %d-ethylene system in %.1f A periodic box...', args.n_molecules, args.box_size)
    positions, species = build_ethylene_box(args.n_molecules, args.box_size, rng)

    # Orthorhombic periodic cell
    cell = np.diag([args.box_size, args.box_size, args.box_size])
    logger.info('System: %d atoms, PBC cell = %.1f x %.1f x %.1f A',
                len(species), args.box_size, args.box_size, args.box_size)

    template, groups = build_template_and_groups(args.n_molecules)

    langevin_params = LangevinParams(
        temperature_K=args.temperature,
        friction_per_fs=0.01,
    )
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

    logger.info('Loading MACE-MP-0 (%s, %s)...', args.model, args.device)
    calc = create_mace_calculator(model=args.model, device=args.device)

    integrator = LangevinIntegrator(langevin_params)
    tracker = BondTracker()
    barostat = None if args.no_barostat else MCBarostat(
        MCBarostatParams(pressure_atm=args.pressure, frequency=25)
    )
    if barostat:
        logger.info('NPT barostat enabled: P=%.2f atm (frequency=25 steps)', args.pressure)
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
        'Starting MACE+PBC TDBB: %d cycles x (%d biased + %d unbiased), T=%.0f K',
        config.n_cycles, config.biased_steps, config.unbiased_steps, args.temperature,
    )

    wf = PolymerizationWorkflow(
        config, calc, template, groups,
        integrator=integrator,
        bond_tracker=tracker,
        barostat=barostat,
    )
    logs = wf.run(
        state,
        output_dir=args.output_dir,
        config_path='configs/boost/paper_faithful.yaml',
    )

    n_form = len(tracker.confirmed_formations())
    n_dissoc = len(tracker.confirmed_dissociations())
    logger.info('Confirmed formations: %d, dissociations: %d', n_form, n_dissoc)

    summary = {
        'total_steps': state.step,
        'n_molecules': args.n_molecules,
        'n_atoms': len(species),
        'box_size_A': args.box_size,
        'cell_periodic': True,
        'backend': calc.name,
        'temperature_K': langevin_params.temperature_K,
        'biased_steps': args.biased_steps,
        'unbiased_steps': args.unbiased_steps,
        'n_cycles': args.n_cycles,
        'confirmed_formations': n_form,
        'confirmed_dissociations': n_dissoc,
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
    logger.info('confirmed_formations=%d', n_form)

    print('\nTo generate figures:')
    print(f'  python scripts/reproduce_figures.py '
          f'--trajectory {args.output_dir}/trajectory.jsonl '
          f'--bonds {args.output_dir}/bonds.jsonl '
          f'--n-reactive-sites {args.n_molecules * 2} '
          f'--target-temperature {args.temperature} '
          f'--output-dir {args.output_dir}/figures')


if __name__ == '__main__':
    main()
