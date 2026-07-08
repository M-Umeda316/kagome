"""Carothers equation for step-growth polymerization.

DPn = 1 / (1 - p)

where p is the extent of reaction (conversion) and DPn is the
number-average degree of polymerization.

Paper anchor: arXiv:2511.22874, Fig. 4c — Carothers theoretical curve
compared to TDBB-simulated DPn vs conversion.
"""
from __future__ import annotations

import logging
from collections.abc import Iterable, Sequence

import numpy as np
from numpy.typing import NDArray

logger = logging.getLogger(__name__)


def dpn_carothers(p: NDArray[np.floating] | float) -> NDArray[np.floating]:
    """Carothers equation: DPn = 1 / (1 - p).

    ``p`` (extent of reaction) is clamped to 0.9999 so that p >= 1.0 yields a
    finite DPn (10000) instead of ``inf``/negative, matching ``dpn_from_bonds``.
    A clamp fires a warning because p >= 1.0 signals full conversion or a
    bond-count/denominator miscount (F9, specs/fix-plan-2026-07-06).
    """
    p_arr = np.asarray(p, dtype=np.float64)
    if np.any(p_arr > 0.9999):
        n_clamped = int(np.count_nonzero(p_arr > 0.9999))
        p_max = float(np.max(p_arr))
        logger.warning(
            'dpn_carothers: %d extent-of-reaction value(s) exceed the 0.9999 '
            'clamp (max p=%.6f); DPn capped at 10000. p>=1 would indicate full '
            'conversion or a bond-count/denominator error.',
            n_clamped, p_max,
        )
    p_clamped = np.minimum(p_arr, 0.9999)
    return 1.0 / (1.0 - p_clamped)


def dpn_from_bonds(
    n_bonds: int,
    n_functional_groups: int,
) -> float:
    """Compute DPn from bond count and initial functional group count.

    Equimolar A-A + B-B step-growth only (e.g., nylon-6,6, N_A = N_B):
      p = n_bonds / (n_functional_groups / 2)
      DPn = 1 / (1 - p)

    Stoichiometric imbalance (r != 1) and monofunctional end-capping are NOT
    modeled here. ``n_functional_groups`` is the total number of reactive end
    groups (A + B); for equimolar systems each type contributes
    ``n_functional_groups / 2``.

    Counting convention (A5): ``n_bonds`` is the number of amide bonds, i.e.
    ``confirmed_formation`` events with ``counts_as_reaction=True``. Bias-only
    water-forming k-l events must be excluded by the caller so one condensation
    counts once (specs/decisions.md 2026-07-06).
    """
    if n_functional_groups <= 0:
        return 1.0
    raw_p = n_bonds / (n_functional_groups / 2)
    if raw_p > 0.9999:
        # Surface the clamp instead of silently masking near-complete conversion
        # or a p>1 bond/denominator miscount (RF19b).
        logger.warning(
            'dpn_from_bonds: extent of reaction p=%.6f exceeds the 0.9999 clamp '
            '(n_bonds=%d, n_functional_groups=%d); DPn capped at 10000. p>1 would '
            'indicate a bond-count or denominator error.',
            raw_p, n_bonds, n_functional_groups,
        )
    p = min(raw_p, 0.9999)
    return 1.0 / (1.0 - p)


# ── Measured DPn from the reacted bond graph (Fig. 4c) ───────────────────────
# The functions above are THEORY (DPn = 1/(1-p)).  The two below compute the
# MEASURED number-average DPn directly from the simulated connectivity, so the
# Fig. 4c overlay is a genuine measured-vs-theory comparison rather than
# theory-vs-theory.  Kept pure (bonds + membership in, float out) so they are
# unit-testable without running MD.  Paper anchor: arXiv:2511.22874, Fig. 4c.


def _union_find_components(
    n_nodes: int,
    edges: Iterable[tuple[int, int]],
) -> int:
    """Number of connected components of a graph on ``n_nodes`` nodes.

    Weighted-quotient-free union-find with path halving; edges reference node
    ids in ``range(n_nodes)``.  Isolated nodes each count as one component.
    """
    parent = list(range(n_nodes))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    return len({find(i) for i in range(n_nodes)})


def monomer_sets_from_bonds(
    bonds: Iterable[Sequence[float]],
) -> list[list[int]]:
    """Atom-level connected components of a bond graph, as sorted index lists.

    Applied to the initial (cycle -1) snapshot of ``topology.jsonl`` this
    recovers the atom membership of each initial monomer molecule (every atom of
    a nylon-6,6 monomer is intramolecularly bonded), i.e. the ``monomer_atom_sets``
    argument for :func:`dpn_measured_from_topology`.  Deriving membership from the
    initial topology keeps the Fig. 4c pipeline self-contained (no SMILES/count
    replay needed).

    ``bonds`` items are ``(i, j)`` or ``(i, j, order)`` (order ignored).
    """
    adjacency: dict[int, list[int]] = {}
    for bond in bonds:
        i, j = int(bond[0]), int(bond[1])
        adjacency.setdefault(i, []).append(j)
        adjacency.setdefault(j, []).append(i)

    seen: set[int] = set()
    components: list[list[int]] = []
    for start in adjacency:
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        comp: list[int] = []
        while stack:
            node = stack.pop()
            comp.append(node)
            for nbr in adjacency[node]:
                if nbr not in seen:
                    seen.add(nbr)
                    stack.append(nbr)
        components.append(sorted(comp))

    components.sort(key=lambda c: c[0])
    return components


def dpn_measured_from_topology(
    bonds: Iterable[Sequence[float]],
    monomer_atom_sets: Sequence[Iterable[int]],
) -> float:
    """Measured number-average degree of polymerization from the bond graph.

    ``DPn = (number of monomer units) / (number of molecules)``, where molecules
    are the connected components of the *monomer* graph: each initial monomer is
    one node, and any bond joining an atom of one monomer to an atom of a
    different monomer is an edge.  Intramonomer bonds are ignored (they never
    merge distinct monomers), so passing the full ``topology.jsonl`` bond list is
    correct.  For nylon-6,6 the amide (amine_N-carboxyl_C) and its paired
    water-forming (amine_H-carboxyl_OH) bond join the SAME diamine/diacid pair,
    so the component count is unaffected by whether one or both are present.

    This is the MEASURED counterpart to the theoretical :func:`dpn_carothers`;
    the two agree for an ideal linear step-growth ensemble (Fig. 4c).

    Parameters
    ----------
    bonds:
        Final bond list; each item is ``(i, j)`` or ``(i, j, order)``.  Typically
        the latest snapshot from ``topology.jsonl``.
    monomer_atom_sets:
        ``monomer_atom_sets[m]`` is the atom indices of initial monomer ``m``
        (e.g. from :func:`monomer_sets_from_bonds` on the initial topology).

    Returns
    -------
    float
        DPn >= 1.0.  Returns 1.0 when there are no monomers.
    """
    n_monomers = len(monomer_atom_sets)
    if n_monomers == 0:
        return 1.0

    atom2mono: dict[int, int] = {}
    for m, atoms in enumerate(monomer_atom_sets):
        for a in atoms:
            atom2mono[int(a)] = m

    edges: list[tuple[int, int]] = []
    for bond in bonds:
        mi = atom2mono.get(int(bond[0]))
        mj = atom2mono.get(int(bond[1]))
        if mi is None or mj is None or mi == mj:
            continue
        edges.append((mi, mj))

    n_components = _union_find_components(n_monomers, edges)
    return n_monomers / n_components
