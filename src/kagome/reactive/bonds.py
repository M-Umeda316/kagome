"""Bond event tracking for reaction monitoring.

Records formation/dissociation attempts during biased phases,
detects tentative reactions in-bias (to end the biased segment),
and confirms outcomes after unbiased relaxation.
"""
from __future__ import annotations

import json
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from kagome.boost.tdbb import PairBias
from kagome.geometry import minimum_image


def is_formed(r: float, r0: float, threshold_fraction: float = 1.0) -> bool:
    """True when distance *r* indicates bond formation (r <= threshold)."""
    return r <= threshold_fraction * r0


def is_dissociated(r: float, r0: float, threshold_fraction: float = 1.0) -> bool:
    """True when distance *r* indicates bond dissociation (r > threshold)."""
    return r > threshold_fraction * r0


@dataclass
class BondEvent:
    """結合イベントの記録 (attempt/confirm, formation/dissociation)。"""

    step: int
    cycle: int
    atom_a: int
    atom_b: int
    event_type: str
    distance: float
    r0: float = 0.0
    candidate_id: int = -1
    # Whether a confirmed formation counts toward alpha/Carothers p. Default
    # True keeps old bonds.jsonl (field absent) backward compatible; nylon
    # water-forming k-l events carry False (specs/decisions.md 2026-07-06 A5).
    counts_as_reaction: bool = True


class BondTracker:
    """Tracks bond formation and dissociation events across cycles."""

    def __init__(self, threshold_fraction: float = 1.0) -> None:
        self._threshold_fraction = threshold_fraction
        self._events: list[BondEvent] = []
        self._pending: list[tuple[PairBias, int]] = []
        # Pairs confirmed after unbiased relaxation — keyed by
        # (min_idx, max_idx, is_formation) so the same pair can undergo
        # formation and later dissociation in different cycles.
        self._reacted: set[tuple[int, int, bool]] = set()
        # Pairs tentatively detected during the current biased phase.
        # Cleared at the start of each record_attempts call.
        self._tentative: set[tuple[int, int]] = set()

    @staticmethod
    def _key(a: int, b: int) -> tuple[int, int]:
        return (min(a, b), max(a, b))

    def _pair_satisfied(self, pair: PairBias, r: float) -> bool:
        """True when *pair* meets its bonding condition at distance *r*.

        Paper §2.2 step 3-4: a formation pair is satisfied when the atoms are
        bonded (r <= r0); a dissociation pair when they are NOT bonded (r > r0),
        where r0 = 60% of the vdW-radii sum (encoded in ``pair.r0``).
        """
        if pair.is_formation:
            return is_formed(r, pair.r0, self._threshold_fraction)
        return is_dissociated(r, pair.r0, self._threshold_fraction)

    def _pair_distance(
        self,
        pair: PairBias,
        positions: NDArray[np.floating],
        cell: NDArray[np.floating] | None,
    ) -> float:
        r_vec = minimum_image(
            positions[pair.idx_b] - positions[pair.idx_a], cell,
        )
        return float(np.linalg.norm(r_vec))

    def check_reactions_during_bias(
        self,
        pairs: list[PairBias],
        positions: NDArray[np.floating],
        step: int,
        cycle: int,
        cell: NDArray[np.floating] | None = None,
        pair_dists: list[float] | None = None,
    ) -> list[int]:
        """Detect whether any candidate's full reaction event fires this step.

        Two things happen here (paper §2.2 step 3-4):

        1. *Audit trail* — every pair that individually crosses its bonding
           threshold is recorded once as a ``tentative_formation`` /
           ``tentative_dissociation`` event (deduplicated per pair via
           ``_tentative``).  These are audit records only; confirmation happens
           in ``check_outcomes`` after unbiased relaxation.

        2. *Firing decision* — pairs are grouped by ``candidate_id`` and a
           candidate FIRES iff ALL of its TRIGGER pairs (``is_trigger=True``)
           simultaneously satisfy their bonding condition at the current
           distances (formation bonded AND dissociation broken).  This is the
           paper's conjunction over the identification set P: the biased phase
           must not end until the leaving groups (N-H, C-OH) have dissociated
           together with the amide formation, otherwise the amide reverts.  The
           bias-only water pair (``is_trigger=False``) never gates firing.

        A candidate with ``candidate_id < 0`` (activation / legacy / pre-
        candidate_id paths) falls back to the OLD behaviour: each such pair is
        its own singleton candidate that fires when it alone (being a trigger)
        satisfies.

        ``pair_dists`` (optional) is the per-pair distance list already
        computed by ``total_bias_fast`` for the SAME ``positions`` — reused
        bit-identically when supplied (specs/decisions.md 2026-07-06 S2/W3);
        otherwise the minimum image is recomputed.

        Returns the list of fired ``candidate_id`` values (empty => no candidate
        completed this step); the caller only needs to test it for emptiness.
        """
        dists: list[float] = []
        for i, pair in enumerate(pairs):
            if pair_dists is not None:
                dists.append(pair_dists[i])
            else:
                dists.append(self._pair_distance(pair, positions, cell))

        # 1. Audit trail — record each newly-crossed pair once.
        for i, pair in enumerate(pairs):
            pair_key = self._key(pair.idx_a, pair.idx_b)
            reacted_key = (*pair_key, pair.is_formation)
            if reacted_key in self._reacted:
                continue
            if pair_key in self._tentative:
                continue
            if not self._pair_satisfied(pair, dists[i]):
                continue
            etype = ('tentative_formation' if pair.is_formation
                     else 'tentative_dissociation')
            self._events.append(BondEvent(
                step=step, cycle=cycle,
                atom_a=pair.idx_a, atom_b=pair.idx_b,
                event_type=etype, distance=dists[i], r0=pair.r0,
                candidate_id=pair.candidate_id,
                counts_as_reaction=pair.counts_as_reaction,
            ))
            self._tentative.add(pair_key)

        # 2. Firing decision — conjunction over each candidate's trigger pairs.
        grouped: dict[int, list[tuple[PairBias, float]]] = defaultdict(list)
        singletons: list[tuple[PairBias, float]] = []
        for i, pair in enumerate(pairs):
            if pair.candidate_id < 0:
                singletons.append((pair, dists[i]))
            else:
                grouped[pair.candidate_id].append((pair, dists[i]))

        fired: list[int] = []
        for cid, members in grouped.items():
            triggers = [(p, r) for (p, r) in members if p.is_trigger]
            if not triggers:
                continue
            if all(self._pair_satisfied(p, r) for (p, r) in triggers):
                fired.append(cid)
        for pair, r in singletons:
            if pair.is_trigger and self._pair_satisfied(pair, r):
                fired.append(pair.candidate_id)
        return fired

    def record_attempts(
        self,
        pairs: list[PairBias],
        positions: NDArray[np.floating],
        step: int,
        cycle: int,
        cell: NDArray[np.floating] | None = None,
    ) -> None:
        self._pending.clear()
        self._tentative.clear()
        for pair in pairs:
            r_vec = minimum_image(
                positions[pair.idx_b] - positions[pair.idx_a], cell,
            )
            r = float(np.linalg.norm(r_vec))
            etype = 'attempted_formation' if pair.is_formation else 'attempted_dissociation'
            self._events.append(BondEvent(
                step=step, cycle=cycle,
                atom_a=pair.idx_a, atom_b=pair.idx_b,
                event_type=etype, distance=r, r0=pair.r0,
                candidate_id=pair.candidate_id,
                counts_as_reaction=pair.counts_as_reaction,
            ))
            self._pending.append((pair, cycle))

    def _confirm_pair(
        self,
        pair: PairBias,
        r: float,
        step: int,
        cycle: int,
        confirmed: list[BondEvent],
    ) -> None:
        """Emit a confirmed BondEvent for *pair* if it individually satisfies.

        Formation pairs (incl. the bias-only water pair) emit
        ``confirmed_formation`` only when bonded at the relaxed geometry;
        dissociation pairs emit ``confirmed_dissociation``.  De-duplicated via
        ``_reacted`` on (min, max, is_formation).
        """
        pair_key = self._key(pair.idx_a, pair.idx_b)
        reacted_key = (*pair_key, pair.is_formation)
        if reacted_key in self._reacted:
            return
        if pair.is_formation:
            if not is_formed(r, pair.r0, self._threshold_fraction):
                return
            etype = 'confirmed_formation'
        else:
            etype = 'confirmed_dissociation'
        ev = BondEvent(
            step=step, cycle=cycle,
            atom_a=pair.idx_a, atom_b=pair.idx_b,
            event_type=etype, distance=r, r0=pair.r0,
            candidate_id=pair.candidate_id,
            counts_as_reaction=pair.counts_as_reaction,
        )
        self._events.append(ev)
        self._reacted.add(reacted_key)
        confirmed.append(ev)

    def record_confirmed_dissociation(
        self,
        *,
        step: int,
        cycle: int,
        atom_a: int,
        atom_b: int,
        distance: float,
        r0: float = 0.0,
        candidate_id: int = -1,
        counts_as_reaction: bool = True,
        register_reacted: bool = False,
    ) -> BondEvent:
        """Append a ``confirmed_dissociation`` event through the public API.

        Single source of truth for adding a confirmed dissociation to
        ``_events`` and (optionally) the ``_reacted`` de-duplication set, so
        callers do not reach into the private list by hand (which left
        ``_reacted`` inconsistent with ``_events``).

        The activation path (V^d on azo C-N bonds) historically appended to
        ``_events`` WITHOUT touching ``_reacted``; ``register_reacted`` defaults
        to ``False`` to preserve that behaviour byte-for-byte. Set it ``True``
        for callers that want the ``(min, max, is_formation=False)`` key
        de-duplicated the same way :meth:`_confirm_pair` does.
        """
        ev = BondEvent(
            step=step, cycle=cycle,
            atom_a=atom_a, atom_b=atom_b,
            event_type='confirmed_dissociation',
            distance=distance, r0=r0,
            candidate_id=candidate_id,
            counts_as_reaction=counts_as_reaction,
        )
        self._events.append(ev)
        if register_reacted:
            self._reacted.add((*self._key(atom_a, atom_b), False))
        return ev

    def check_outcomes(
        self,
        positions: NDArray[np.floating],
        step: int,
        cell: NDArray[np.floating] | None = None,
    ) -> list[BondEvent]:
        """Confirm reactions after unbiased relaxation (paper §2.2 step 4).

        The confirmation re-checks the SAME per-candidate trigger conjunction
        used to fire the biased phase: a candidate confirms iff ALL of its
        trigger pairs (``is_trigger=True``) still satisfy their bonding
        condition at the relaxed geometry.  On confirmation the WHOLE candidate
        is committed — each formation pair (incl. the bias-only water pair, if
        it closed) emits ``confirmed_formation`` and each dissociation pair
        emits ``confirmed_dissociation`` — so ``_apply_topology_edits`` can add
        C-N, remove N-H / C-OH and add the water O-H consistently.  A candidate
        whose conjunction does NOT hold confirms NOTHING (no spurious lone
        dissociations).

        Candidates with ``candidate_id < 0`` (activation / legacy) keep the OLD
        per-pair independent confirmation, so those paths are unchanged.
        """
        confirmed: list[BondEvent] = []

        grouped: dict[int, list[tuple[PairBias, int, float]]] = defaultdict(list)
        singletons: list[tuple[PairBias, int, float]] = []
        for pair, cycle in self._pending:
            r = self._pair_distance(pair, positions, cell)
            if pair.candidate_id < 0:
                singletons.append((pair, cycle, r))
            else:
                grouped[pair.candidate_id].append((pair, cycle, r))

        # Per-candidate conjunctive confirmation over trigger pairs.
        for cid, members in grouped.items():
            triggers = [(p, r) for (p, _, r) in members if p.is_trigger]
            if not triggers:
                continue
            if not all(self._pair_satisfied(p, r) for (p, r) in triggers):
                continue
            for pair, cycle, r in members:
                self._confirm_pair(pair, r, step, cycle, confirmed)

        # Legacy / activation: independent per-pair confirmation.
        for pair, cycle, r in singletons:
            self._confirm_pair(pair, r, step, cycle, confirmed)

        self._pending.clear()
        return confirmed

    @property
    def events(self) -> list[BondEvent]:
        return list(self._events)

    def confirmed_formations(self) -> list[BondEvent]:
        return [e for e in self._events if e.event_type == 'confirmed_formation']

    def confirmed_dissociations(self) -> list[BondEvent]:
        return [e for e in self._events if e.event_type == 'confirmed_dissociation']

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, 'w', encoding='utf-8') as f:
            for ev in self._events:
                f.write(json.dumps(asdict(ev)) + '\n')
