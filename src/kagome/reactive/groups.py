"""Reactive group definitions and reaction templates.

Paper: arXiv:2511.22874, Mori et al.
Equations 6-7.
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class ReactiveGroup:
    """A set of chemically equivalent atoms.  Eq. 6: G_X = {a | a ∈ X}."""
    label: str
    atom_indices: list[int] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.atom_indices)

    def remove_atom(self, idx: int) -> None:
        try:
            self.atom_indices.remove(idx)
        except ValueError:
            pass


@dataclass
class PairSpec:
    """Specifies a pair of groups within a reaction template.

    group_a/group_b are labels referencing ReactiveGroup.label.
    is_formation: True for bond-forming, False for bond-breaking.
    r_min, r_max: distance bounds (Å) for candidate selection (Eq. 7).
    constraint_only: if True, participates in candidate identification and
        scoring but no bias force is applied (Table S1: i-k, j-l constraints).
    score_pair: if False, excluded from candidate identification (distance
        window filter) and d_ijkl scoring.  Bias is still applied based on
        is_formation/constraint_only.  Use for nylon k-l (H-OH water
        formation) which is bias-only per Table S2 (Eq. 7: d_ijkl = r_ij +
        r_ik + r_jl, 3 terms fixed).
    count_as_reaction: if False, confirmed formations of this pair are NOT
        counted toward alpha(t) / Carothers p.  The bias and topology effect
        are unchanged.  Use for nylon k-l (amine_H-carboxyl_OH, water O-H
        formation): one condensation = one amide bond (amine_N-carboxyl_C),
        so the paired water-forming event must not be double-counted
        (specs/decisions.md 2026-07-06, A5).
    """
    group_a: str
    group_b: str
    is_formation: bool
    r_min: float = 0.0
    r_max: float = 5.0
    constraint_only: bool = False
    score_pair: bool = True
    count_as_reaction: bool = True


@dataclass
class ReactionTemplate:
    """Defines a reaction type by its participating groups and pair interactions.

    groups: ordered list of group labels [I, J, K, L, ...].
    pairs: which pairs of groups interact and how.
    """
    name: str
    groups: list[str]
    pairs: list[PairSpec]

    def __post_init__(self) -> None:
        if len(self.groups) != len(set(self.groups)):
            dupes = [g for g in self.groups if self.groups.count(g) > 1]
            raise ValueError(
                f'ReactionTemplate {self.name!r}: duplicate group labels {set(dupes)}')
        for ps in self.pairs:
            if ps.group_a not in self.groups:
                raise ValueError(
                    f'ReactionTemplate {self.name!r}: pair group_a {ps.group_a!r} '
                    f'not in groups {self.groups}')
            if ps.group_b not in self.groups:
                raise ValueError(
                    f'ReactionTemplate {self.name!r}: pair group_b {ps.group_b!r} '
                    f'not in groups {self.groups}')
            if ps.group_a == ps.group_b:
                # Symmetric pairs map to a (i,i) key that the enumeration
                # loop (prev_depth < depth) can never match — the distance
                # window would be silently disabled (L4).
                raise ValueError(
                    f'ReactionTemplate {self.name!r}: pair ({ps.group_a!r}, '
                    f'{ps.group_b!r}) references the same group; symmetric '
                    'reactions are not supported by the candidate enumerator')

    def group_labels(self) -> set[str]:
        return set(self.groups)

    def pair_labels(self) -> list[tuple[str, str]]:
        return [(p.group_a, p.group_b) for p in self.pairs]
