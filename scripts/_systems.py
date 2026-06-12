"""Shared system builders for polymerization scripts."""
from __future__ import annotations

import numpy as np

from src.reactive.groups import PairSpec, ReactiveGroup, ReactionTemplate


def build_ethylene_box(
    n_molecules: int,
    box_size: float,
    rng: np.random.Generator,
    min_sep: float = 3.0,
) -> tuple[np.ndarray, list[str]]:
    """Place n ethylene (C2H4) molecules with atom-level overlap prevention.

    Uses a for-else pattern: raises RuntimeError if valid placement is not
    found in 200 attempts per molecule (instead of silently accepting overlap).
    Overlap check is atom-vs-atom (not center-vs-atom) so min_sep is exact.
    """
    from ase.build import molecule

    ethylene = molecule('C2H4')
    template_pos = ethylene.get_positions()
    template_symbols = list(ethylene.get_chemical_symbols())

    placed_mols: list[np.ndarray] = []
    all_species: list[str] = []
    centered = template_pos - template_pos.mean(axis=0)

    for mol_idx in range(n_molecules):
        for _attempt in range(200):
            offset = rng.uniform(2.0, box_size - 2.0, size=3)
            angle = rng.uniform(0.0, 2.0 * np.pi)
            c, s = np.cos(angle), np.sin(angle)
            rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
            candidate = centered @ rot.T + offset

            ok = True
            for prev in placed_mols:
                diffs = candidate[:, np.newaxis, :] - prev[np.newaxis, :, :]
                if np.min(np.linalg.norm(diffs, axis=2)) < min_sep:
                    ok = False
                    break
            if ok:
                placed_mols.append(candidate)
                all_species.extend(template_symbols)
                break
        else:
            raise RuntimeError(
                f'Could not place molecule {mol_idx + 1}/{n_molecules} without overlap '
                f'in box_size={box_size:.1f} A after 200 attempts. '
                f'Try increasing --box-size.'
            )

    return np.vstack(placed_mols), all_species


def build_template_and_groups(
    n_molecules: int,
) -> tuple[ReactionTemplate, dict[str, ReactiveGroup]]:
    """Define reactive groups for ethylene C=C bond activation."""
    atoms_per_mol = 6  # C2H4: 2C + 4H

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
