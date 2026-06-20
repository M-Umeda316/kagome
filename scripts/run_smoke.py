"""Smoke test: run a toy TDBB polymerization end-to-end.

Usage: python scripts/run_smoke.py [--output-dir runs/smoke]
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np

from kagome.backends.toy import ToyCalculator
from kagome.boost.tdbb import TDBBParams
from kagome.reactive.bonds import BondTracker
from kagome.reactive.groups import PairSpec, ReactiveGroup, ReactionTemplate
from kagome.workflows.polymerization import (
    PolymerizationConfig,
    PolymerizationWorkflow,
    SimulationState,
    masses_from_species,
)

logging.basicConfig(level=logging.INFO, format='%(name)s | %(message)s')


def main() -> None:
    parser = argparse.ArgumentParser(description='TDBB smoke test')
    parser.add_argument('--output-dir', type=Path, default=Path('runs/smoke'))
    parser.add_argument('--seed', type=int, default=7)
    args = parser.parse_args()

    config = PolymerizationConfig(
        timestep_fs=0.25,
        biased_steps=200,
        unbiased_steps=200,
        n_cycles=3,
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

    template = ReactionTemplate(
        name='toy_dimerization',
        groups=['A', 'B'],
        pairs=[PairSpec('A', 'B', is_formation=True, r_min=0.5, r_max=4.0)],
    )

    rng = np.random.default_rng(args.seed)
    n_atoms = 8
    positions = rng.uniform(0, 6, size=(n_atoms, 3))
    species = ['C'] * n_atoms

    groups = {
        'A': ReactiveGroup('A', list(range(0, n_atoms, 2))),
        'B': ReactiveGroup('B', list(range(1, n_atoms, 2))),
    }

    state = SimulationState(
        positions=positions,
        velocities=np.zeros((n_atoms, 3)),
        species=species,
        masses=masses_from_species(species),
    )

    calc = ToyCalculator()
    tracker = BondTracker()
    wf = PolymerizationWorkflow(config, calc, template, groups, bond_tracker=tracker)
    logs = wf.run(state, output_dir=args.output_dir, config_path='configs/boost/paper_faithful.yaml')

    summary = {
        'total_steps': state.step,
        'cycles': len(logs) // 2,
        'confirmed_formations': len(tracker.confirmed_formations()),
        'confirmed_dissociations': len(tracker.confirmed_dissociations()),
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

    print(f'Done. Summary written to {out_path}')
    print(f'Trajectory: {args.output_dir / "trajectory.jsonl"}')
    print(f'Bond events: {args.output_dir / "bonds.jsonl"}')
    print(f'\nTo generate figures:')
    print(f'  python scripts/reproduce_figures.py --trajectory {args.output_dir / "trajectory.jsonl"} --output-dir {args.output_dir / "figures"}')


if __name__ == '__main__':
    main()
