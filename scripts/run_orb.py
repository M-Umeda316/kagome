"""TDBB polymerization with OrbMol-v2 backend.

Runs a small vinyl-like system (ethylene molecules) to demonstrate
TDBB bond formation with OrbMol-v2 (trained on OMol25 + OPoly26).

Usage: python scripts/run_orb.py [--output-dir runs/orb] [--n-molecules 4]
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
from src.backends.orb_backend import create_orb_calculator
from src.boost.tdbb import TDBBParams
from src.integrators.init_velocities import maxwell_boltzmann_velocities
from src.integrators.langevin import LangevinIntegrator, LangevinParams
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
    parser = argparse.ArgumentParser(description='OrbMol-v2 TDBB polymerization')
    parser.add_argument('--output-dir', type=Path, default=Path('runs/orb'))
    parser.add_argument('--seed', type=int, default=7)
    parser.add_argument('--n-molecules', type=int, default=4)
    parser.add_argument('--n-cycles', type=int, default=5)
    parser.add_argument('--biased-steps', type=int, default=100)
    parser.add_argument('--unbiased-steps', type=int, default=100)
    parser.add_argument('--box-size', type=float, default=None)
    parser.add_argument('--device', type=str, default='cpu')
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    logger.info('Building %d-ethylene system...', args.n_molecules)
    box_size = args.box_size or (5.0 + args.n_molecules * 1.5)
    positions, species = build_ethylene_box(args.n_molecules, box_size, rng)

    logger.info('System: %d atoms (%s), box=%.1f A',
                len(species), ', '.join(sorted(set(species))), box_size)

    template, groups = build_template_and_groups(args.n_molecules)

    langevin_params = LangevinParams(temperature_K=333.0)
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
        save_interval=10,
    )

    logger.info('Loading OrbMol-v2 (%s)...', args.device)
    calc = create_orb_calculator(device=args.device)

    langevin = LangevinIntegrator(langevin_params)
    tracker = BondTracker()

    masses = masses_from_species(species)
    velocities = maxwell_boltzmann_velocities(masses, 333.0, rng)

    state = SimulationState(
        positions=positions,
        velocities=velocities,
        species=species,
        cell=None,
        masses=masses,
    )

    logger.info('Starting TDBB polymerization: %d cycles x (%d biased + %d unbiased)',
                config.n_cycles, config.biased_steps, config.unbiased_steps)

    wf = PolymerizationWorkflow(
        config, calc, template, groups,
        integrator=langevin,
        bond_tracker=tracker,
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
        'box_size_A': box_size,
        'backend': calc.name,
        'temperature_K': langevin_params.temperature_K,
        'confirmed_formations': n_form,
        'confirmed_dissociations': n_dissoc,
        'cycles': len(logs) // 2,
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

    logger.info('Done. Outputs in %s', args.output_dir)
    print(f'\nTo generate figures:')
    print(f'  python scripts/reproduce_figures.py '
          f'--trajectory {args.output_dir / "trajectory.jsonl"} '
          f'--bonds {args.output_dir / "bonds.jsonl"} '
          f'--output-dir {args.output_dir / "figures"}')


if __name__ == '__main__':
    main()
