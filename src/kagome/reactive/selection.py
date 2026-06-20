"""Candidate reaction group selection with non-overlap constraints.

Paper: arXiv:2511.22874, Mori et al.
Equation 7, scoring d_ijkl, greedy non-overlapping selection.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from kagome.geometry import minimum_image
from kagome.reactive.groups import ReactiveGroup, ReactionTemplate, PairSpec


@dataclass
class Candidate:
    """A candidate reaction group: one atom per template group."""
    atom_indices: tuple[int, ...]
    score: float = 0.0


def _distance(
    positions: NDArray[np.floating],
    i: int,
    j: int,
    cell: NDArray[np.floating] | None = None,
) -> float:
    r_vec = minimum_image(positions[j] - positions[i], cell)
    return float(np.linalg.norm(r_vec))


def find_candidates(
    template: ReactionTemplate,
    groups: dict[str, ReactiveGroup],
    positions: NDArray[np.floating],
    cell: NDArray[np.floating] | None = None,
) -> list[Candidate]:
    """Enumerate candidate tuples satisfying distance bounds.  Eq. 7.

    Only pairs with score_pair=True participate in candidate identification
    (distance window filtering).  Pairs with score_pair=False (e.g. nylon
    k-l water formation) are bias-only and do not constrain candidate
    selection.
    """
    group_atoms = [groups[label].atom_indices for label in template.groups]
    label_list = template.groups

    pair_specs: dict[tuple[int, int], PairSpec] = {}
    for ps in template.pairs:
        if not ps.score_pair:
            continue
        idx_a = label_list.index(ps.group_a)
        idx_b = label_list.index(ps.group_b)
        pair_specs[(min(idx_a, idx_b), max(idx_a, idx_b))] = ps

    candidates: list[Candidate] = []
    _enumerate_recursive(
        group_atoms, label_list, pair_specs, positions, candidates,
        depth=0, chosen=[], running_score=0.0, cell=cell,
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
    running_score: float,
    cell: NDArray[np.floating] | None = None,
) -> None:
    if depth == len(group_atoms):
        out.append(Candidate(
            atom_indices=tuple(chosen),
            score=running_score,
        ))
        return

    for atom_idx in group_atoms[depth]:
        ok = True
        added_score = 0.0

        for prev_depth in range(depth):
            key = (prev_depth, depth)
            if key not in pair_specs:
                continue
            ps = pair_specs[key]
            d = _distance(positions, chosen[prev_depth], atom_idx, cell)
            if d < ps.r_min or d > ps.r_max:
                ok = False
                break
            added_score += d

        if not ok:
            continue

        chosen.append(atom_idx)
        _enumerate_recursive(
            group_atoms, label_list, pair_specs, positions, out,
            depth + 1, chosen, running_score + added_score, cell=cell,
        )
        chosen.pop()


def score_candidates(
    candidates: list[Candidate],
) -> list[Candidate]:
    """Sort candidates ascending by pre-computed score (d_ijkl, Eq. 7).

    Scores are computed during find_candidates (enumeration phase) as the
    sum of score_pair=True pair distances: d_ijkl = r_ij + r_ik + r_jl.
    """
    candidates.sort(key=lambda c: c.score)
    return candidates


@dataclass
class SelectionDecision:
    """Audit record for one candidate in the greedy non-overlap pass (RF18).

    reason is 'selected' or 'overlap:[...]' listing the already-used atoms that
    caused the candidate to be dropped.
    """
    atom_indices: tuple[int, ...]
    score: float
    selected: bool
    reason: str


def audited_selection(
    candidates: list[Candidate],
) -> tuple[list[Candidate], list[SelectionDecision]]:
    """Greedy non-overlapping selection that also returns a per-candidate audit.

    Single source of the selection logic; ``select_non_overlapping`` is a thin
    wrapper that discards the audit. ``candidates`` is assumed already sorted by
    score (see ``score_candidates``), so the audit reflects the true ranking.
    """
    used_atoms: set[int] = set()
    selected: list[Candidate] = []
    decisions: list[SelectionDecision] = []

    for c in candidates:
        atoms = set(c.atom_indices)
        clash = atoms & used_atoms
        if clash:
            decisions.append(SelectionDecision(
                c.atom_indices, c.score, False, f'overlap:{sorted(clash)}',
            ))
            continue
        selected.append(c)
        used_atoms |= atoms
        decisions.append(SelectionDecision(
            c.atom_indices, c.score, True, 'selected',
        ))

    return selected, decisions


def select_non_overlapping(candidates: list[Candidate]) -> list[Candidate]:
    """Greedy non-overlapping selection: skip candidates with already-used atoms."""
    selected, _ = audited_selection(candidates)
    return selected
