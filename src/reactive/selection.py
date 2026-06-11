"""Candidate reaction group selection with non-overlap constraints.

Paper: arXiv:2511.22874, Mori et al.
Equation 7, scoring d_ijkl, greedy non-overlapping selection.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from src.reactive.groups import ReactiveGroup, ReactionTemplate, PairSpec


@dataclass
class Candidate:
    """A candidate reaction group: one atom per template group."""
    atom_indices: tuple[int, ...]
    pair_distances: dict[tuple[str, str], float]
    score: float = 0.0


def _distance(
    positions: NDArray[np.floating],
    i: int,
    j: int,
) -> float:
    return float(np.linalg.norm(positions[j] - positions[i]))


def find_candidates(
    template: ReactionTemplate,
    groups: dict[str, ReactiveGroup],
    positions: NDArray[np.floating],
) -> list[Candidate]:
    """Enumerate candidate tuples satisfying distance bounds.  Eq. 7.

    For each combination of one atom per group, check all pair distance bounds.
    Returns candidates that pass all constraints.
    """
    group_atoms = [groups[label].atom_indices for label in template.groups]
    label_list = template.groups

    pair_specs: dict[tuple[int, int], PairSpec] = {}
    for ps in template.pairs:
        idx_a = label_list.index(ps.group_a)
        idx_b = label_list.index(ps.group_b)
        pair_specs[(idx_a, idx_b)] = ps

    candidates: list[Candidate] = []
    _enumerate_recursive(
        group_atoms, label_list, pair_specs, positions, candidates,
        depth=0, chosen=[], chosen_distances={},
    )
    return candidates


def _enumerate_recursive(
    group_atoms: list[list[int]],
    label_list: list[str],
    pair_specs: dict[tuple[int, int], PairSpec],
    positions: NDArray[np.floating],
    out: list[Candidate],
    depth: int,
    chosen: list[int],
    chosen_distances: dict[tuple[str, str], float],
) -> None:
    if depth == len(group_atoms):
        out.append(Candidate(
            atom_indices=tuple(chosen),
            pair_distances=dict(chosen_distances),
        ))
        return

    for atom_idx in group_atoms[depth]:
        ok = True
        new_distances: dict[tuple[str, str], float] = {}

        for prev_depth in range(depth):
            key = (prev_depth, depth)
            if key not in pair_specs:
                continue
            ps = pair_specs[key]
            d = _distance(positions, chosen[prev_depth], atom_idx)
            if d < ps.r_min or d > ps.r_max:
                ok = False
                break
            new_distances[(ps.group_a, ps.group_b)] = d

        if not ok:
            continue

        chosen.append(atom_idx)
        merged = {**chosen_distances, **new_distances}
        _enumerate_recursive(
            group_atoms, label_list, pair_specs, positions, out,
            depth + 1, chosen, merged,
        )
        chosen.pop()


def score_candidates(
    candidates: list[Candidate],
    template: ReactionTemplate,
    positions: NDArray[np.floating],
) -> list[Candidate]:
    """Score each candidate: d = sum of all pair distances.  Sort ascending."""
    label_list = template.groups
    for c in candidates:
        total = 0.0
        for ps in template.pairs:
            idx_a = label_list.index(ps.group_a)
            idx_b = label_list.index(ps.group_b)
            d = _distance(positions, c.atom_indices[idx_a], c.atom_indices[idx_b])
            total += d
        c.score = total

    candidates.sort(key=lambda c: c.score)
    return candidates


def select_non_overlapping(candidates: list[Candidate]) -> list[Candidate]:
    """Greedy non-overlapping selection: skip candidates with already-used atoms."""
    used_atoms: set[int] = set()
    selected: list[Candidate] = []

    for c in candidates:
        atoms = set(c.atom_indices)
        if atoms & used_atoms:
            continue
        selected.append(c)
        used_atoms |= atoms

    return selected
