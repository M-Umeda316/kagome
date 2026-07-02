"""Explicit bond topology tracking for reactive MD.

Why this exists
---------------
The coordinate trajectory (``trajectory.jsonl`` / exported XYZ) carries no
connectivity, so viewers (Winmostar / OVITO / VMD) infer bonds from
inter-atomic distance.  Under strong TDBB pair bias (f1 up to 250 kcal/mol)
whole molecules are dragged into near-contact, so distance inference draws
spurious bonds at non-reactive sites; and because the vinyl C=C double bond is
never opened, reacted carbons appear over-valent.  Both are *visualization*
artifacts of missing topology, not recorded chemistry (``bonds.jsonl`` only ever
records the radical_C-vinyl_alpha_C i-j pair).  See specs/decisions.md
2026-07-02.

This module maintains an explicit, chemically-correct bond list: the initial
intramolecular topology from the molecule builder, plus reaction edits applied
on each confirmed formation.  Emitting it lets viewers use real connectivity
instead of distance guessing.

Paper anchor: Table S1 radical addition (radical + C=C -> C-C, unpaired electron
migrates to the beta carbon).  The edit is valence-conserving.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import logging

logger = logging.getLogger(__name__)

BondKey = tuple[int, int]

# Nominal maximum coordination (number of sigma bonds) per element.  Used to
# detect when a closed-shell radical centre is already saturated and must shed a
# placeholder H on reaction (see apply_vinyl_addition).
MAX_COORDINATION: dict[str, int] = {
    'H': 1, 'C': 4, 'N': 3, 'O': 2, 'F': 1, 'S': 2, 'Cl': 1,
}


@dataclass
class BondTopology:
    """Mutable set of bonds with per-bond order, keyed by sorted atom pair."""

    _orders: dict[BondKey, float] = field(default_factory=dict)

    @staticmethod
    def _key(a: int, b: int) -> BondKey:
        return (a, b) if a < b else (b, a)

    @classmethod
    def from_bonds(cls, bonds) -> BondTopology:
        """Build from an iterable of ``(i, j)`` or ``(i, j, order)`` tuples."""
        topo = cls()
        for bond in bonds:
            if len(bond) == 3:
                i, j, order = bond
            else:
                i, j = bond
                order = 1.0
            topo.add_bond(int(i), int(j), float(order))
        return topo

    def add_bond(self, a: int, b: int, order: float = 1.0) -> None:
        if a == b:
            raise ValueError(f'self-bond not allowed: atom {a}')
        self._orders[self._key(a, b)] = float(order)

    def remove_bond(self, a: int, b: int) -> None:
        self._orders.pop(self._key(a, b), None)

    def set_order(self, a: int, b: int, order: float) -> None:
        key = self._key(a, b)
        if key not in self._orders:
            raise KeyError(f'no bond between {a} and {b} to set order on')
        self._orders[key] = float(order)

    def has_bond(self, a: int, b: int) -> bool:
        return self._key(a, b) in self._orders

    def order(self, a: int, b: int) -> float:
        return self._orders.get(self._key(a, b), 0.0)

    def neighbors(self, atom: int) -> list[int]:
        out = []
        for (i, j) in self._orders:
            if i == atom:
                out.append(j)
            elif j == atom:
                out.append(i)
        return out

    def coordination_number(self, atom: int) -> int:
        """Number of bonded neighbours (ignores bond order)."""
        return len(self.neighbors(atom))

    def bonds(self) -> list[tuple[int, int, float]]:
        """All bonds as sorted ``(i, j, order)`` tuples, deterministically ordered."""
        return [(i, j, self._orders[(i, j)]) for (i, j) in sorted(self._orders)]

    def copy(self) -> BondTopology:
        return BondTopology(dict(self._orders))

    def __len__(self) -> int:
        return len(self._orders)


def _spare_hydrogen(
    topology: BondTopology, atom: int, species: list[str] | None,
) -> int | None:
    """Return an H neighbour of *atom* (a placeholder to shed), or None."""
    if species is None:
        return None
    for nbr in topology.neighbors(atom):
        if nbr < len(species) and species[nbr] == 'H':
            return nbr
    return None


def apply_vinyl_addition(
    topology: BondTopology,
    radical_c: int,
    alpha_c: int,
    propagation_map: dict[int, int],
    species: list[str] | None = None,
) -> None:
    """Apply the radical-addition topology edit for one confirmed formation.

    Paper (Table S1): the chain-end radical adds across the monomer vinyl group.
    Valence-conserving edits:

    1. If the radical centre is already 4-coordinate, it is the closed-shell
       initiator model (isobutyronitrile: the radical C carries a placeholder H
       standing in for the unpaired electron; see _systems ``_INITIATOR_SMILES``).
       Shed that H bond so the new C-C bond does not over-coordinate the carbon.
       A propagating radical (beta carbon, 3-coordinate) needs no shedding.
    2. Add the new radical_C--vinyl_alpha_C sigma bond (order 1).
    3. Open the monomer's alpha=beta double bond to a single bond (the unpaired
       electron migrates to beta, which becomes the next chain-end radical).

    ``species`` is required for step 1; without it the shed is skipped and an
    over-valent centre is logged.
    """
    max_c = MAX_COORDINATION.get(
        species[radical_c] if species and radical_c < len(species) else 'C', 4)
    if topology.coordination_number(radical_c) >= max_c:
        spare = _spare_hydrogen(topology, radical_c, species)
        if spare is not None:
            topology.remove_bond(radical_c, spare)
        else:
            logger.warning(
                'apply_vinyl_addition: radical C %d already %d-coordinate and no '
                'spare H to shed; emitted topology will be over-valent at %d',
                radical_c, topology.coordination_number(radical_c), radical_c,
            )

    topology.add_bond(radical_c, alpha_c, 1.0)

    beta = propagation_map.get(alpha_c)
    if beta is not None and topology.has_bond(alpha_c, beta):
        topology.set_order(alpha_c, beta, 1.0)  # C=C -> C-C


def over_coordinated_atoms(
    topology: BondTopology, species: list[str] | None = None,
) -> list[int]:
    """Return atoms whose coordination number exceeds their element maximum."""
    over: list[int] = []
    seen: set[int] = set()
    for (i, j) in topology._orders:
        seen.add(i)
        seen.add(j)
    for atom in seen:
        el = species[atom] if species and atom < len(species) else 'C'
        if topology.coordination_number(atom) > MAX_COORDINATION.get(el, 4):
            over.append(atom)
    return sorted(over)


def vinyl_addition_over_coordinates(
    topology: BondTopology,
    radical_c: int,
    alpha_c: int,
    propagation_map: dict[int, int],
    species: list[str] | None = None,
) -> list[int]:
    """Dry-run the radical-addition edit on a *copy* and report which atoms it
    would leave over-coordinated (empty list => the formation is valence-safe).

    Used by the Layer-2 occupancy guard to refuse biasing / confirming a
    formation that cannot happen without violating valence (e.g. an alpha carbon
    that already reacted, or a radical centre with no spare valence and no H to
    shed).  See specs/decisions.md 2026-07-02.
    """
    trial = topology.copy()
    apply_vinyl_addition(trial, radical_c, alpha_c, propagation_map, species)
    # Only the two reacting centres (and their partners) can change coordination.
    beta = propagation_map.get(alpha_c)
    affected = {radical_c, alpha_c}
    if beta is not None:
        affected.add(beta)
    bad = []
    for atom in sorted(affected):
        el = species[atom] if species and atom < len(species) else 'C'
        if trial.coordination_number(atom) > MAX_COORDINATION.get(el, 4):
            bad.append(atom)
    return bad
