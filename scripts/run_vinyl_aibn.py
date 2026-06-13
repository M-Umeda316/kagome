"""Radical vinyl polymerization with methyl acrylate + AIBN initiator.

System:  n_monomers methyl acrylate (C=CC(=O)OC) +
         n_initiators isobutyronitrile radical model (CC(C)C#N)
Backend: MACE-MP-0 (MIT licence, commercial-safe)
Ensemble: NPT (Langevin + MC barostat) or NVT with --no-barostat

Design:  specs/decisions.md — "T-G1: vinyl radical polymerization system"
Paper:   arXiv:2511.22874, Fig. 1 / Section 2

Usage:
    python scripts/run_vinyl_aibn.py --seed 7 --output-dir runs/vinyl_aibn
    python scripts/run_vinyl_aibn.py --seed 7 --output-dir runs/vinyl_aibn_nvt --no-barostat
"""
from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path

os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')

import numpy as np

from scripts._systems import build_vinyl_aibn_system
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
    parser = argparse.ArgumentParser(
        description='Radical vinyl polymerization: methyl acrylate + AIBN (MACE-MP-0)'
    )
    parser.add_argument('--seed', type=int, default=7)
    parser.add_argument('--output-dir', type=Path, default=Path('runs/vinyl_aibn'))
    parser.add_argument('--n-monomers', type=int, default=8)
    parser.add_argument('--n-initiators', type=int, default=2)
    parser.add_argument('--n-cycles', type=int, default=3)
    parser.add_argument('--biased-steps', type=int, default=500)
    parser.add_argument('--unbiased-steps', type=int, default=500)
    parser.add_argument('--box-size', type=float, default=14.0)
    parser.add_argument('--temperature', type=float, default=500.0)
    parser.add_argument('--pressure', type=float, default=1.0,
                        help='Target pressure (atm). Default 1.0 (assumed, not stated in paper).')
    parser.add_argument('--no-barostat', action='store_true',
                        help='Disable NPT barostat and run NVT instead.')
    parser.add_argument('--device', type=str, default='cpu')
    parser.add_argument('--model', type=str, default='small')
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    logger.info(
        'Building vinyl/AIBN system: %d monomers + %d initiators in %.1f Å box...',
        args.n_monomers, args.n_initiators, args.box_size,
    )
    positions, species, template, groups, propagation_map = build_vinyl_aibn_system(
        n_monomers=args.n_monomers,
        n_initiators=args.n_initiators,
        box_size=args.box_size,
        rng=rng,
    )
    logger.info(
        'System: %d atoms total  (%d radical_C, %d vinyl_alpha_C sites)',
        len(species),
        len(groups['radical_C'].atom_indices),
        len(groups['vinyl_alpha_C'].atom_indices),
    )
    logger.info('Propagation map: %d entries', len(propagation_map))

    cell = np.diag([args.box_size, args.box_size, args.box_size])

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
        propagation_map=propagation_map,
        propagation_target_group='radical_C',
    )
    logs = wf.run(
        state,
        output_dir=args.output_dir,
        config_path='configs/boost/paper_faithful.yaml',
    )

    n_form = len(tracker.confirmed_formations())
    n_dissoc = len(tracker.confirmed_dissociations())
    logger.info('Confirmed formations: %d, dissociations: %d', n_form, n_dissoc)

    n_reactive_sites = len(groups['radical_C'].atom_indices) + len(groups['vinyl_alpha_C'].atom_indices)
    summary = {
        'total_steps': state.step,
        'n_monomers': args.n_monomers,
        'n_initiators': args.n_initiators,
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
        'propagation_events': n_form,  # each formation triggers one propagation
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

    n_reactive = args.n_monomers * 2 + args.n_initiators
    print('\nTo generate figures:')
    print(
        f'  python scripts/reproduce_figures.py '
        f'--trajectory {args.output_dir}/trajectory.jsonl '
        f'--bonds {args.output_dir}/bonds.jsonl '
        f'--n-reactive-sites {n_reactive} '
        f'--target-temperature {args.temperature} '
        f'--output-dir {args.output_dir}/figures'
    )


if __name__ == '__main__':
    main()
