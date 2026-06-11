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
    """
    group_a: str
    group_b: str
    is_formation: bool
    r_min: float = 0.0
    r_max: float = 5.0


@dataclass
class ReactionTemplate:
    """Defines a reaction type by its participating groups and pair interactions.

    groups: ordered list of group labels [I, J, K, L, ...].
    pairs: which pairs of groups interact and how.
    """
    name: str
    groups: list[str]
    pairs: list[PairSpec]

    def group_labels(self) -> set[str]:
        return set(self.groups)

    def pair_labels(self) -> list[tuple[str, str]]:
        return [(p.group_a, p.group_b) for p in self.pairs]
