"""T8.1: Demonstrate TDBB bond formation with a toy LJ potential.

Purpose: verify the full TDBB machinery (bias -> attempt -> confirm) produces
confirmed_formations >= 1 in a system where the underlying potential has a well
at the TDBB target distance.

System:
  - 2 "C" atoms starting at ~3.5 Angstrom separation.
  - ToyCalculator with sigma = r0_CC = lambda_vdw * (vdwC + vdwC) = 0.6*3.4 = 2.04 A.
    The LJ minimum sits exactly at the TDBB target distance.
  - epsilon = 10 kcal/mol >> kT(500K)=0.99 kcal/mol, so atoms stay trapped.

Acceptance criterion (T8.1): confirmed_formations >= 1.

Usage:
    python scripts/run_toy_bond_demo.py --seed 7 --output-dir runs/toy_bond_demo
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np

from src.backends.toy import ToyCalculator
from src.boost.tdbb import TDBBParams
from src.integrators.verlet import VelocityVerletIntegrator
from src.reactive.bonds import BondTracker
from src.reactive.groups import PairSpec, ReactiveGroup, ReactionTemplate
from src.workflows.polymerization import (
    PolymerizationConfig,
    PolymerizationWorkflow,
    SimulationState,
)

logging.basicConfig(level=logging.INFO, format='%(name)s | %(message)s')
logger = logging.getLogger(__name__)

# C-C TDBB target: r0 = lambda_vdw * (vdwC + vdwC) = 0.6 * (1.7+1.7) = 2.04 A
_R0_CC = 2.04
_EPSILON = 10.0  # kcal/mol — deep well, >> kT at 500K (0.99 kcal/mol)


def build_two_atom_system(
    initial_distance: float = 3.5,
) -> tuple[np.ndarray, list[str], ReactionTemplate, dict[str, ReactiveGroup]]:
    """Two C atoms separated by initial_distance A."""
    positions = np.array([
        [0.0, 0.0, 0.0],
        [initial_distance, 0.0, 0.0],
    ], dtype=np.float64)
    species = ['C', 'C']

    template = ReactionTemplate(
        name='cc_bond_formation',
        groups=['C_donor', 'C_acceptor'],
        pairs=[
            PairSpec(
                group_a='C_donor',
                group_b='C_acceptor',
                is_formation=True,
                r_min=1.6,
                r_max=5.0,
            ),
        ],
    )
    groups = {
        'C_donor': ReactiveGroup('C_donor', [0]),
        'C_acceptor': ReactiveGroup('C_acceptor', [1]),
    }
    return positions, species, template, groups


def main() -> None:
    parser = argparse.ArgumentParser(description='T8.1 Toy LJ bond formation demo')
    parser.add_argument('--seed', type=int, default=7)
    parser.add_argument('--output-dir', type=Path, default=Path('runs/toy_bond_demo'))
    parser.add_argument('--biased-steps', type=int, default=300)
    parser.add_argument('--unbiased-steps', type=int, default=300)
    parser.add_argument('--n-cycles', type=int, default=2)
    args = parser.parse_args()

    positions, species, template, groups = build_two_atom_system(initial_distance=3.5)

    # ToyCalculator: LJ minimum at sigma = r0_CC = 2.04 A, deep well (epsilon >> kT)
    calc = ToyCalculator(epsilon=_EPSILON, sigma=_R0_CC)

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
        save_interval=20,
    )

    tracker = BondTracker()
    integrator = VelocityVerletIntegrator()

    state = SimulationState(
        positions=positions,
        velocities=np.zeros_like(positions),
        species=species,
        cell=None,
        masses=np.array([12.011, 12.011]),
    )

    logger.info(
        'System: 2 C atoms, initial separation=3.5 A, LJ sigma=%.2f A, epsilon=%.1f kcal/mol',
        _R0_CC, _EPSILON,
    )
    logger.info(
        'TDBB r0=%.2f A, confirm threshold=%.2f A (threshold_fraction=1.0, paper-faithful)',
        _R0_CC, 1.0 * _R0_CC,
    )
    logger.info(
        'Schedule: %d biased + %d unbiased x %d cycles',
        config.biased_steps, config.unbiased_steps, config.n_cycles,
    )

    wf = PolymerizationWorkflow(
        config, calc, template, groups,
        integrator=integrator,
        bond_tracker=tracker,
    )
    logs = wf.run(state, output_dir=args.output_dir, config_path='toy_bond_demo')

    n_form = len(tracker.confirmed_formations())
    logger.info('confirmed_formations=%d (target: >= 1)', n_form)

    for ev in tracker.confirmed_formations():
        logger.info('  Bond confirmed at step=%d, cycle=%d, r=%.3f A (threshold=%.3f A)',
                    ev.step, ev.cycle, ev.distance, 1.0 * ev.r0)

    summary = {
        'backend': calc.name,
        'system': '2-C-atom LJ well',
        'lj_sigma_A': _R0_CC,
        'lj_epsilon_kcal_mol': _EPSILON,
        'confirmed_formations': n_form,
        'confirmed_dissociations': len(tracker.confirmed_dissociations()),
        'total_steps': state.step,
        'logs': [
            {
                'cycle': log.cycle, 'phase': log.phase,
                'n_candidates': log.n_candidates, 'n_selected': log.n_selected,
                'bias_energy': log.bias_energy,
            }
            for log in logs
        ],
    }

    out = args.output_dir / 'summary.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2), encoding='utf-8')

    if n_form >= 1:
        logger.info('T8.1 PASSED: confirmed_formations=%d >= 1', n_form)
    else:
        logger.error('T8.1 FAILED: confirmed_formations=%d < 1', n_form)


if __name__ == '__main__':
    main()
