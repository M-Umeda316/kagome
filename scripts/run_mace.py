"""Paper-faithful TDBB polymerization with MACE-MP-0 backend.

Runs a small vinyl-like system (ethylene molecules) to demonstrate
TDBB bond formation with a real uMLIP.

Usage: python scripts/run_mace.py [--output-dir runs/mace] [--n-molecules 4]
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np

from src.backends.mace_backend import create_mace_calculator
from src.boost.tdbb import TDBBParams
from src.integrators.langevin import LangevinIntegrator, LangevinParams
from src.reactive.bonds import BondTracker
from src.reactive.groups import PairSpec, ReactiveGroup, ReactionTemplate
from src.workflows.polymerization import (
    PolymerizationConfig,
    PolymerizationWorkflow,
    SimulationState,
    masses_from_species,
)

logging.basicConfig(level=logging.INFO, format='%(name)s | %(message)s')
logger = logging.getLogger(__name__)


def build_ethylene_box(
    n_molecules: int,
    box_size: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, list[str]]:
    """Place n ethylene (C2H4) molecules randomly in a periodic box.

    Each ethylene: C=C bond ~1.34 Å, 4 H atoms at ~1.08 Å from C.
    """
    from ase.build import molecule
    from ase import Atoms

    ethylene = molecule('C2H4')
    template_pos = ethylene.get_positions()
    template_symbols = ethylene.get_chemical_symbols()
    atoms_per_mol = len(template_symbols)

    all_positions = []
    all_species = []

    for i in range(n_molecules):
        offset = rng.uniform(2.0, box_size - 2.0, size=3)
        angle = rng.uniform(0, 2 * np.pi, size=3)

        # Simple rotation around z-axis
        c, s = np.cos(angle[2]), np.sin(angle[2])
        rot = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

        centered = template_pos - template_pos.mean(axis=0)
        rotated = centered @ rot.T
        placed = rotated + offset

        all_positions.append(placed)
        all_species.extend(template_symbols)

    positions = np.vstack(all_positions)
    return positions, all_species


def build_template_and_groups(
    n_molecules: int,
) -> tuple[ReactionTemplate, dict[str, ReactiveGroup]]:
    """Define reactive groups for ethylene C=C bond activation.

    Group A: first C of each C=C pair (carbon index 0 per molecule)
    Group B: second C of each C=C pair (carbon index 1 per molecule)

    The TDBB formation bias will try to bring C atoms from different
    molecules closer, simulating chain growth initiation.
    """
    atoms_per_mol = 6  # C2H4: 2C + 4H

    # C atoms that can form inter-molecular bonds
    group_a_indices = [i * atoms_per_mol + 0 for i in range(n_molecules)]
    group_b_indices = [i * atoms_per_mol + 1 for i in range(n_molecules)]

    template = ReactionTemplate(
        name='vinyl_polymerization',
        groups=['C_donor', 'C_acceptor'],
        pairs=[
            PairSpec(
                group_a='C_donor',
                group_b='C_acceptor',
                is_formation=True,
                r_min=1.6,
                r_max=4.5,
            ),
        ],
    )

    groups = {
        'C_donor': ReactiveGroup('C_donor', group_a_indices),
        'C_acceptor': ReactiveGroup('C_acceptor', group_b_indices),
    }

    return template, groups


def main() -> None:
    parser = argparse.ArgumentParser(description='MACE TDBB polymerization')
    parser.add_argument('--output-dir', type=Path, default=Path('runs/mace'))
    parser.add_argument('--seed', type=int, default=7)
    parser.add_argument('--n-molecules', type=int, default=4)
    parser.add_argument('--n-cycles', type=int, default=5)
    parser.add_argument('--biased-steps', type=int, default=100)
    parser.add_argument('--unbiased-steps', type=int, default=100)
    parser.add_argument('--device', type=str, default='cpu')
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    logger.info('Building %d-ethylene system...', args.n_molecules)
    box_size = 8.0 + args.n_molecules * 1.5
    positions, species = build_ethylene_box(args.n_molecules, box_size, rng)
    cell = np.diag([box_size, box_size, box_size])

    logger.info('System: %d atoms (%s), box=%.1f Å',
                len(species), ', '.join(sorted(set(species))), box_size)

    template, groups = build_template_and_groups(args.n_molecules)

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

    logger.info('Loading MACE-MP-0 (%s)...', args.device)
    calc = create_mace_calculator(model='small', device=args.device)

    langevin = LangevinIntegrator(LangevinParams(
        temperature_K=500.0,
        friction_per_fs=0.01,
    ))

    tracker = BondTracker(threshold_fraction=1.3)

    state = SimulationState(
        positions=positions,
        velocities=np.zeros_like(positions),
        species=species,
        cell=cell,
        masses=masses_from_species(species),
    )

    logger.info('Starting TDBB polymerization: %d cycles × (%d biased + %d unbiased)',
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
        'temperature_K': 500.0,
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
    print(f'  python scripts/reproduce_figures.py --trajectory {args.output_dir / "trajectory.jsonl"} --output-dir {args.output_dir / "figures"}')


if __name__ == '__main__':
    main()
