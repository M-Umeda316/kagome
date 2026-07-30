"""Candidate reaction group selection with non-overlap constraints.

Paper: arXiv:2511.22874, Mori et al.
Equation 7, scoring d_ijkl, greedy non-overlapping selection.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from kagome.geometry import minimum_image_fast, validated_box
from kagome.reactive.groups import ReactiveGroup, ReactionTemplate, PairSpec
from kagome.reactive.topology import BondTopology


@dataclass
class Candidate:
    """A candidate reaction group: one atom per template group."""
    atom_indices: tuple[int, ...]
    score: float = 0.0


def _distance(
    positions: NDArray[np.floating],
    i: int,
    j: int,
    box: NDArray[np.floating] | None = None,
) -> float:
    """Minimum-image distance between atoms *i* and *j*.

    ``box`` is a prevalidated (3,) diagonal from ``validated_box`` (the
    orthorhombicity check is hoisted to ``find_candidates`` and done once),
    so this uses ``minimum_image_fast``.  Numerically identical to the previous
    per-pair ``minimum_image`` — same ``r - box*round(r/box)`` rounding.
    """
    r_vec = minimum_image_fast(positions[j] - positions[i], box)
    return float(np.linalg.norm(r_vec))


def find_candidates(
    template: ReactionTemplate,
    groups: dict[str, ReactiveGroup],
    positions: NDArray[np.floating],
    cell: NDArray[np.floating] | None = None,
    topology: BondTopology | None = None,
) -> list[Candidate]:
    """Enumerate candidate tuples satisfying distance bounds.  Eq. 7.

    Only pairs with score_pair=True participate in candidate identification
    (distance window filtering).  Pairs with score_pair=False (e.g. nylon
    k-l water formation) are bias-only and do not constrain candidate
    selection.

    When *topology* is provided, formation pairs whose atoms are already
    bonded are excluded (L6 guard).
    """
    group_atoms = [groups[label].atom_indices for label in template.groups]
    label_list = template.groups

    pair_specs: dict[tuple[int, int], PairSpec] = {}
    formation_group_pairs: set[tuple[int, int]] = set()
    dissociation_group_pairs: set[tuple[int, int]] = set()
    for ps in template.pairs:
        if not ps.score_pair:
            continue
        idx_a = label_list.index(ps.group_a)
        idx_b = label_list.index(ps.group_b)
        key = (min(idx_a, idx_b), max(idx_a, idx_b))
        pair_specs[key] = ps
        if not ps.constraint_only:
            if ps.is_formation:
                formation_group_pairs.add(key)
            else:
                dissociation_group_pairs.add(key)

    # Validate the orthorhombic cell once here; every _distance call reuses the
    # resulting (3,) box via minimum_image_fast (numerically identical).
    box = validated_box(cell)

    candidates: list[Candidate] = []
    _enumerate_recursive(
        group_atoms, label_list, pair_specs, positions, candidates,
        depth=0, chosen=[], running_score=0.0, box=box,
        topology=topology, formation_group_pairs=formation_group_pairs,
        dissociation_group_pairs=dissociation_group_pairs,
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
    box: NDArray[np.floating] | None = None,
    topology: BondTopology | None = None,
    formation_group_pairs: set[tuple[int, int]] | None = None,
    dissociation_group_pairs: set[tuple[int, int]] | None = None,
) -> None:
    if depth == len(group_atoms):
        out.append(Candidate(
            atom_indices=tuple(chosen),
            score=running_score,
        ))
        return

    for atom_idx in group_atoms[depth]:
        if atom_idx in chosen:
            continue
        ok = True
        added_score = 0.0

        for prev_depth in range(depth):
            key = (prev_depth, depth)
            if key not in pair_specs:
                continue
            ps = pair_specs[key]
            d = _distance(positions, chosen[prev_depth], atom_idx, box)
            if d < ps.r_min or d > ps.r_max:
                ok = False
                break
            group_key = (min(prev_depth, depth), max(prev_depth, depth))
            # L6: a formation pair whose atoms are already bonded cannot form.
            if (topology is not None and formation_group_pairs
                    and group_key in formation_group_pairs
                    and topology.has_bond(chosen[prev_depth], atom_idx)):
                ok = False
                break
            # B2 (dual of L6): a dissociation pair can only break a bond that
            # actually exists. Reject candidates that pair a dissociation atom
            # (e.g. nylon amine_N) with a partner it is NOT bonded to (a
            # different molecule's H), which would otherwise spawn spurious
            # inter-molecular N-H / C-OH dissociation candidates.
            if (topology is not None and dissociation_group_pairs
                    and group_key in dissociation_group_pairs
                    and not topology.has_bond(chosen[prev_depth], atom_idx)):
                ok = False
                break
            added_score += d

        if not ok:
            continue

        chosen.append(atom_idx)
        _enumerate_recursive(
            group_atoms, label_list, pair_specs, positions, out,
            depth + 1, chosen, running_score + added_score, box=box,
            topology=topology, formation_group_pairs=formation_group_pairs,
            dissociation_group_pairs=dissociation_group_pairs,
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
    """Audit record for one candidate in the non-overlap selection pass (RF18).

    reason is 'selected' or 'overlap:[...]' listing the already-used atoms that
    caused the candidate to be dropped. ``pool_size`` is set only for
    policy='softmax' picks (WM-P5a): the number of candidates still viable
    (non-overlapping with what was already selected) that competed for that
    particular pick. It stays None for policy='deterministic', so existing
    consumers of this dataclass are unaffected.
    """
    atom_indices: tuple[int, ...]
    score: float
    selected: bool
    reason: str
    pool_size: int | None = None


# Below this softmax_temperature, weights would collapse to (near) a one-hot
# on the best score anyway; treat it as exactly deterministic to avoid 0/0 and
# float-overflow edge cases as T -> 0. This is what makes the T -> 0 limit of
# policy='softmax' recover policy='deterministic' exactly (WM-P5a).
_MIN_SOFTMAX_TEMPERATURE = 1e-9


def _select_deterministic(
    candidates: list[Candidate],
) -> tuple[list[Candidate], list[SelectionDecision]]:
    """Greedy non-overlapping selection that also returns a per-candidate audit.

    ``candidates`` is assumed already sorted by score (see ``score_candidates``),
    so the audit reflects the true ranking. This is the paper-faithful default
    (best-score-first greedy, Eq. 7) and its behaviour is frozen: WM-P5a adds
    policy='softmax' as an opt-in alternative alongside it, never changing this
    path.
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


def _softmax_probs(
    scores: NDArray[np.floating], temperature: float,
) -> NDArray[np.floating]:
    """exp(-(score - min_score)/T), normalised (WM-P5a). Lower score is better,
    so the minimum-score candidate always gets the largest weight; subtracting
    the min keeps the exponent <= 0 for numerical stability."""
    shifted = scores - scores.min()
    weights = np.exp(-shifted / temperature)
    return weights / weights.sum()


def _select_softmax(
    candidates: list[Candidate],
    temperature: float,
    rng: np.random.Generator,
) -> tuple[list[Candidate], list[SelectionDecision]]:
    """Non-overlapping selection where each pick is a softmax draw over the
    candidates currently viable (non-overlapping with what's already been
    selected), instead of always taking the best score first (WM-P5a).

    ``candidates`` is assumed sorted by score, purely so ties and the T -> 0
    limit agree with ``_select_deterministic``; the draw itself does not
    depend on input order.
    """
    if temperature < 0.0:
        raise ValueError(f'softmax_temperature must be >= 0, got {temperature}')

    used_atoms: set[int] = set()
    selected: list[Candidate] = []
    decisions: list[SelectionDecision] = []
    pool = list(candidates)

    while pool:
        viable: list[Candidate] = []
        for c in pool:
            clash = set(c.atom_indices) & used_atoms
            if clash:
                decisions.append(SelectionDecision(
                    c.atom_indices, c.score, False, f'overlap:{sorted(clash)}',
                ))
            else:
                viable.append(c)
        if not viable:
            break

        if len(viable) == 1 or temperature <= _MIN_SOFTMAX_TEMPERATURE:
            pick_idx = int(np.argmin([c.score for c in viable]))
        else:
            scores = np.array([c.score for c in viable], dtype=np.float64)
            probs = _softmax_probs(scores, temperature)
            pick_idx = int(rng.choice(len(viable), p=probs))

        picked = viable[pick_idx]
        selected.append(picked)
        used_atoms |= set(picked.atom_indices)
        decisions.append(SelectionDecision(
            picked.atom_indices, picked.score, True, 'selected',
            pool_size=len(viable),
        ))
        pool = [c for c in viable if c is not picked]

    return selected, decisions


def audited_selection(
    candidates: list[Candidate],
    *,
    policy: str = 'deterministic',
    softmax_temperature: float | None = None,
    rng: np.random.Generator | None = None,
) -> tuple[list[Candidate], list[SelectionDecision]]:
    """Non-overlapping candidate selection that also returns a per-candidate
    audit trail. Single source of the selection logic; ``select_non_overlapping``
    is a thin wrapper that discards the audit.

    policy:
      - 'deterministic' (default): best-score-first greedy (Eq. 7), the
        paper-faithful behaviour. Bit-identical to pre-WM-P5a callers; extra
        kwargs are accepted but ignored on this path.
      - 'softmax': each pick is a draw with probability proportional to
        exp(-(score - min_score)/softmax_temperature) among the candidates
        still viable at that point. Requires ``softmax_temperature`` and
        ``rng``. ``rng`` must be the caller's single shared
        ``np.random.Generator`` so runs stay reproducible per seed.
        softmax_temperature -> 0 recovers the deterministic order exactly
        (see ``_MIN_SOFTMAX_TEMPERATURE``).

    WM-P5a: optional stochastic candidate-selection policy for well-mixed
    measurement mode; see specs/decisions.md "2026-07-17: well-mixed 測定モード"
    item (iv) if present on your branch.
    """
    if policy == 'deterministic':
        return _select_deterministic(candidates)
    if policy == 'softmax':
        if softmax_temperature is None:
            raise ValueError("policy='softmax' requires softmax_temperature")
        if rng is None:
            raise ValueError("policy='softmax' requires rng")
        return _select_softmax(candidates, softmax_temperature, rng)
    raise ValueError(f'Unknown selection policy: {policy!r}')


def select_non_overlapping(
    candidates: list[Candidate],
    *,
    policy: str = 'deterministic',
    softmax_temperature: float | None = None,
    rng: np.random.Generator | None = None,
) -> list[Candidate]:
    """Non-overlapping selection: skip candidates with already-used atoms.
    See ``audited_selection`` for the ``policy`` options."""
    selected, _ = audited_selection(
        candidates, policy=policy, softmax_temperature=softmax_temperature, rng=rng,
    )
    return selected
