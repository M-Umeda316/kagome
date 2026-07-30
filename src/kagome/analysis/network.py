"""Network analysis for crosslinking polymerization (bulk epoxy-amine curing).

Species are judged purely from the bond graph + element symbols — no MD, no
RDKit, no reactive-group bookkeeping:

- epoxide ring : 3-membered ring of exactly {C, C, O} in the bond graph
- amine rank   : number of H neighbours of an N atom (2 H = primary,
                 1 H = secondary, 0 H = tertiary)
- hydroxyl     : O with exactly one C neighbour and exactly one H neighbour

System context: DGEBA (2 epoxide rings) + DETA (2 primary + 1 secondary
amine N) bulk curing (specs/decisions.md 2026-07-09 Track 2 / E2 design).

Paper anchor: arXiv:2511.22874 Fig. 5 — species concentration traces
c(t)/c(0) during epoxy curing.  The Flory-Stockmayer gel point is a textbook
addition (Flory 1941) used as a comparison baseline for ref#2
(Provenzano 2025); the paper itself does not report a theoretical gel point.
"""
from __future__ import annotations

import json
import logging
import math
from collections.abc import Iterable, Sequence
from itertools import combinations
from pathlib import Path

from kagome.analysis.carothers import monomer_sets_from_bonds

logger = logging.getLogger(__name__)


def _neighbor_map(bonds: Iterable[Sequence[float]]) -> dict[int, set[int]]:
    """Adjacency sets of the bond graph.

    ``bonds`` items are ``(i, j)`` or ``(i, j, order)`` (order ignored),
    matching the ``topology.jsonl`` record format and
    :func:`kagome.analysis.carothers.monomer_sets_from_bonds`.
    """
    adjacency: dict[int, set[int]] = {}
    for bond in bonds:
        i, j = int(bond[0]), int(bond[1])
        adjacency.setdefault(i, set()).add(j)
        adjacency.setdefault(j, set()).add(i)
    return adjacency


def find_epoxide_rings(
    bonds: Iterable[Sequence[float]],
    species: Sequence[str],
) -> list[tuple[int, int, int]]:
    """Closed epoxide rings: 3-membered rings of exactly {C, C, O}.

    Detects triangles in the bond graph whose atoms are two carbons and one
    oxygen.  Returns ``(c1, c2, o)`` tuples with ``c1 < c2``, sorted.  A ring
    that has been opened (any of its three edges removed) no longer forms a
    triangle and is not returned — the count of closed rings per snapshot is
    the epoxide population for Fig. 5-style traces.

    Pure graph pattern: C-C-C triangles (species mismatch) and open C-C-O
    chains or C-O-C ethers (no triangle) are not matched.
    """
    adjacency = _neighbor_map(bonds)
    rings: list[tuple[int, int, int]] = []
    for o, nbrs in adjacency.items():
        if species[o] != 'O':
            continue
        c_neighbors = sorted(n for n in nbrs if species[n] == 'C')
        for c1, c2 in combinations(c_neighbors, 2):
            if c2 in adjacency.get(c1, ()):  # noqa: SIM118 — set membership
                rings.append((c1, c2, o))
    rings.sort()
    return rings


def amine_ranks(
    bonds: Iterable[Sequence[float]],
    species: Sequence[str],
) -> dict[int, int]:
    """H-neighbour count for every N atom in the system.

    Returns ``{n_idx: n_h_neighbours}`` for ALL atoms with species ``'N'``
    (including N with no bonds at all, which reports 0).  Callers interpret
    the rank: 2 H = primary, 1 H = secondary, 0 H = tertiary amine.

    Caveat (documented data property, not a bug): during zwitterion-like
    confirmations (specs/decisions.md 2026-07-09 E1 実走結果) the transferring
    H can be transiently bond-less — the N-H bond is already removed while the
    hydroxyl O-H is not yet formed.  A snapshot taken exactly at such a cycle
    therefore under-counts N-H (an N may appear one rank higher than its
    settled chemical state).
    """
    adjacency = _neighbor_map(bonds)
    ranks: dict[int, int] = {}
    for idx, symbol in enumerate(species):
        if symbol != 'N':
            continue
        ranks[idx] = sum(
            1 for n in adjacency.get(idx, ()) if species[n] == 'H'
        )
    return ranks


def count_hydroxyls(
    bonds: Iterable[Sequence[float]],
    species: Sequence[str],
) -> int:
    """Number of hydroxyl O atoms: exactly one C neighbour and one H neighbour.

    Ring-opening addition turns the epoxide ring O into a beta-hydroxyl
    (C-O-H).  Excluded by the one-C-one-H rule: ether O (2 C, 0 H — DGEBA
    glycidyl/aryl ethers), closed epoxide ring O (2 C), and water O (0 C, 2 H;
    no water exists in this addition chemistry, but the rule is robust to it).
    """
    adjacency = _neighbor_map(bonds)
    count = 0
    for idx, symbol in enumerate(species):
        if symbol != 'O':
            continue
        nbrs = adjacency.get(idx, ())
        n_c = sum(1 for n in nbrs if species[n] == 'C')
        n_h = sum(1 for n in nbrs if species[n] == 'H')
        if n_c == 1 and n_h == 1:
            count += 1
    return count


def species_series(
    snapshots: Sequence[tuple[int, int, list]],
    species: Sequence[str],
) -> list[dict]:
    """Per-snapshot species counts for Fig. 5-style concentration traces.

    ``snapshots`` is ``[(step, cycle, bonds), ...]`` (e.g. from
    :func:`load_topology_snapshots`).  Returns one dict per snapshot:
    ``{'step', 'cycle', 'n_epoxide', 'n_amine_primary', 'n_amine_secondary',
    'n_amine_tertiary', 'n_hydroxyl'}``.

    Amine scope: ALL N atoms in the system are classified by their CURRENT
    H-neighbour count (>= 2 H primary, 1 H secondary, 0 H tertiary).  In the
    DGEBA + DETA system the only N atoms are the DETA amines, so this equals
    the amine population; the >= 2 guard folds hypothetical transient
    over-protonated N into 'primary' rather than dropping it.  See
    :func:`amine_ranks` for the zwitterion-snapshot under-count caveat.
    """
    series: list[dict] = []
    for step, cycle, bonds in snapshots:
        ranks = amine_ranks(bonds, species)
        series.append({
            'step': int(step),
            'cycle': int(cycle),
            'n_epoxide': len(find_epoxide_rings(bonds, species)),
            'n_amine_primary': sum(1 for h in ranks.values() if h >= 2),
            'n_amine_secondary': sum(1 for h in ranks.values() if h == 1),
            'n_amine_tertiary': sum(1 for h in ranks.values() if h == 0),
            'n_hydroxyl': count_hydroxyls(bonds, species),
        })
    return series


def largest_component_fraction(
    bonds: Iterable[Sequence[float]],
    monomer_sets: Sequence[Iterable[int]],
) -> float:
    """Fraction of monomers in the largest connected monomer-level component.

    Gelation indicator: two monomers are connected when any bond joins atoms
    of different monomer sets (``monomer_sets`` as produced by
    :func:`kagome.analysis.carothers.monomer_sets_from_bonds` on the initial
    topology).  Approaches 1.0 as the network percolates; plotted against
    epoxide conversion it locates the simulated gel point next to the
    Flory-Stockmayer prediction (:func:`gel_point_flory_stockmayer`).

    Returns ``(size of largest monomer component) / (number of monomers)``;
    0.0 when there are no monomers.
    """
    n_monomers = len(monomer_sets)
    if n_monomers == 0:
        return 0.0

    atom2mono: dict[int, int] = {}
    for m, atoms in enumerate(monomer_sets):
        for a in atoms:
            atom2mono[int(a)] = m

    # Monomer-level graph: one self-edge per monomer keeps isolated monomers
    # in the component list, inter-monomer bonds are the real edges.  Reuses
    # the carothers component finder instead of a second union-find.
    edges: list[tuple[int, int]] = [(m, m) for m in range(n_monomers)]
    for bond in bonds:
        mi = atom2mono.get(int(bond[0]))
        mj = atom2mono.get(int(bond[1]))
        if mi is None or mj is None or mi == mj:
            continue
        edges.append((mi, mj))

    components = monomer_sets_from_bonds(edges)
    return max(len(c) for c in components) / n_monomers


def max_inter_monomer_degree(
    bonds: Iterable[Sequence[float]],
    monomer_sets: Sequence[Iterable[int]],
) -> int:
    """Largest number of DISTINCT neighbouring monomers any one monomer bonds to.

    Builds the monomer-level graph (nodes = ``monomer_sets`` as produced by
    :func:`kagome.analysis.carothers.monomer_sets_from_bonds` on the initial
    topology; an edge joins two monomers when any bond connects an atom of one
    to an atom of the other) and returns the maximum monomer degree.

    Distinct *neighbours* are counted, not raw bonds, so the nylon amide plus its
    paired water-forming bond — which join the SAME diamine/diacid pair — count
    once. An ideal LINEAR step-growth chain therefore caps this at 2 (interior
    repeat units bond to two neighbours, chain ends to one); a value > 2 means a
    branch point exists, i.e. a monomer of functionality f > 2 (e.g. epoxy-amine
    curing), for which the linear Carothers relation DPn = 1/(1-p) does not hold.

    Returns 0 when there are no monomers or no inter-monomer bonds.
    """
    atom2mono: dict[int, int] = {}
    for m, atoms in enumerate(monomer_sets):
        for a in atoms:
            atom2mono[int(a)] = m

    neighbours: dict[int, set[int]] = {}
    for bond in bonds:
        mi = atom2mono.get(int(bond[0]))
        mj = atom2mono.get(int(bond[1]))
        if mi is None or mj is None or mi == mj:
            continue
        neighbours.setdefault(mi, set()).add(mj)
        neighbours.setdefault(mj, set()).add(mi)

    if not neighbours:
        return 0
    return max(len(n) for n in neighbours.values())


def is_branching_topology(
    bonds: Iterable[Sequence[float]],
    monomer_sets: Sequence[Iterable[int]],
    max_linear_degree: int = 2,
) -> bool:
    """True if any monomer bonds to more than ``max_linear_degree`` neighbours.

    Thin predicate over :func:`max_inter_monomer_degree`; the default threshold
    of 2 flags any network with a branch point (functionality f > 2).
    """
    return max_inter_monomer_degree(bonds, monomer_sets) > max_linear_degree


def gel_point_flory_stockmayer(f: int, g: int, r: float = 1.0) -> float:
    """Flory-Stockmayer gel-point conversion for an f-functional + g-functional
    step-growth system: ``alpha_gel = 1 / sqrt(r * (f - 1) * (g - 1))``.

    Textbook theory (Flory 1941; Flory, "Principles of Polymer Chemistry").
    Not from the reproduced paper — added as a comparison baseline for the
    ref#2 epoxy system (specs/decisions.md 2026-07-09 Track 2 / E2 design).
    For DGEBA (f = 2 epoxide groups) + DETA (g = 5 N-H hydrogens) at
    stoichiometric balance r = 1: alpha_gel = 1/sqrt(1*1*4) = 0.5.

    ``r`` is the stoichiometric imbalance ratio (<= 1 by convention); ``f``
    and ``g`` must both be >= 2 (a monofunctional component never gels).
    """
    if f < 2 or g < 2:
        raise ValueError(f'f and g must be >= 2 for gelation, got f={f}, g={g}')
    if r <= 0.0:
        raise ValueError(f'stoichiometric ratio r must be > 0, got r={r}')
    return 1.0 / math.sqrt(r * (f - 1) * (g - 1))


def crosslink_counts(
    bonds: Iterable[Sequence[float]],
    species: Sequence[str],
    epoxide_rings_initial: Sequence[tuple[int, int, int]],
    monomer_sets: Sequence[Iterable[int]],
) -> dict:
    """Counts of the two crosslink-junction types, kept separate.

    - ``n_tertiary_amine``: N atoms with 0 H neighbours (fully substituted
      amine junction).
    - ``n_fully_reacted_epoxy_monomers``: monomers whose ALL initial epoxide
      rings are now open (every ring lost at least one of its three edges).
      Rings are assigned to monomers by mapping ``epoxide_rings_initial``
      atoms into ``monomer_sets``; monomers with no initial ring (amines) are
      never counted.
    - ``n_monomers``: total monomer count, for normalization by the caller.

    The two junction types are deliberately NOT combined into one density —
    E2 will align the weighting/normalization with ref#2 (Provenzano 2025)
    before reporting a single crosslink density (specs/decisions.md
    2026-07-09 Track 2 / E2 design).
    """
    bond_set: set[tuple[int, int]] = set()
    for bond in bonds:
        i, j = int(bond[0]), int(bond[1])
        bond_set.add((min(i, j), max(i, j)))

    def _has(a: int, b: int) -> bool:
        return (min(a, b), max(a, b)) in bond_set

    ranks = amine_ranks(bonds, species)
    n_tertiary = sum(1 for h in ranks.values() if h == 0)

    atom2mono: dict[int, int] = {}
    for m, atoms in enumerate(monomer_sets):
        for a in atoms:
            atom2mono[int(a)] = m

    rings_open_by_mono: dict[int, list[bool]] = {}
    for c1, c2, o in epoxide_rings_initial:
        mono = atom2mono.get(int(o))
        if mono is None:
            logger.warning(
                'crosslink_counts: epoxide ring O atom %d is not in any '
                'monomer set; ring skipped.', int(o),
            )
            continue
        closed = _has(c1, c2) and _has(c1, o) and _has(c2, o)
        rings_open_by_mono.setdefault(mono, []).append(not closed)

    n_fully_reacted = sum(
        1 for open_flags in rings_open_by_mono.values() if all(open_flags)
    )

    return {
        'n_tertiary_amine': n_tertiary,
        'n_fully_reacted_epoxy_monomers': n_fully_reacted,
        'n_monomers': len(monomer_sets),
    }


def load_topology_snapshots(
    path: Path | str,
) -> list[tuple[int, int, list[tuple[int, int, float]]]]:
    """Thin JSONL loader for ``topology.jsonl`` — the only file-touching
    function in this module.

    Each line is ``{"step": int, "cycle": int, "n_bonds": int,
    "bonds": [[i, j, order], ...]}``; the first record (cycle -1) is the
    initial snapshot.  Returns ``[(step, cycle, bonds), ...]`` sorted by step,
    with ``bonds`` as ``[(i, j, order), ...]``.  Empty list if the file does
    not exist.
    """
    path = Path(path)
    if not path.exists():
        return []
    snapshots: list[tuple[int, int, list[tuple[int, int, float]]]] = []
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            bonds = [(int(i), int(j), float(o)) for i, j, o in rec['bonds']]
            snapshots.append((int(rec['step']), int(rec['cycle']), bonds))
    snapshots.sort(key=lambda s: s[0])
    return snapshots
