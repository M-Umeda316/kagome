"""Bond event tracking for reaction monitoring.

Records formation/dissociation attempts during biased phases,
detects tentative reactions in-bias (to end the biased segment),
and confirms outcomes after unbiased relaxation.
"""
from __future__ import annotations

import json
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

    def check_reactions_during_bias(
        self,
        pairs: list[PairBias],
        positions: NDArray[np.floating],
        step: int,
        cycle: int,
        cell: NDArray[np.floating] | None = None,
    ) -> list[BondEvent]:
        """Detect tentative reaction events DURING the biased phase.

        A formation pair is tentatively detected when its separation falls
        below the vdW bonding threshold; a dissociation pair when it rises
        above it.  Tentative events are recorded for auditing but do NOT
        count as confirmed — confirmation happens only in ``check_outcomes``
        after unbiased relaxation (specs/decisions.md 2026-07-03 D1).

        Returns tentative events so the caller can end the biased segment.
        """
        newly: list[BondEvent] = []
        for pair in pairs:
            pair_key = self._key(pair.idx_a, pair.idx_b)
            reacted_key = (*pair_key, pair.is_formation)
            if reacted_key in self._reacted:
                continue
            if pair_key in self._tentative:
                continue
            r_vec = minimum_image(
                positions[pair.idx_b] - positions[pair.idx_a], cell,
            )
            r = float(np.linalg.norm(r_vec))
            if pair.is_formation:
                reacted = is_formed(r, pair.r0, self._threshold_fraction)
            else:
                reacted = is_dissociated(r, pair.r0, self._threshold_fraction)
            if not reacted:
                continue
            etype = ('tentative_formation' if pair.is_formation
                     else 'tentative_dissociation')
            ev = BondEvent(
                step=step, cycle=cycle,
                atom_a=pair.idx_a, atom_b=pair.idx_b,
                event_type=etype, distance=r, r0=pair.r0,
                candidate_id=pair.candidate_id,
                counts_as_reaction=pair.counts_as_reaction,
            )
            self._events.append(ev)
            self._tentative.add(pair_key)
            newly.append(ev)
        return newly

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

    def check_outcomes(
        self,
        positions: NDArray[np.floating],
        step: int,
        cell: NDArray[np.floating] | None = None,
    ) -> list[BondEvent]:
        confirmed: list[BondEvent] = []
        for pair, cycle in self._pending:
            pair_key = self._key(pair.idx_a, pair.idx_b)
            reacted_key = (*pair_key, pair.is_formation)
            if reacted_key in self._reacted:
                continue
            r_vec = minimum_image(
                positions[pair.idx_b] - positions[pair.idx_a], cell,
            )
            r = float(np.linalg.norm(r_vec))
            if pair.is_formation:
                if is_formed(r, pair.r0, self._threshold_fraction):
                    ev = BondEvent(
                        step=step, cycle=cycle,
                        atom_a=pair.idx_a, atom_b=pair.idx_b,
                        event_type='confirmed_formation',
                        distance=r, r0=pair.r0,
                        candidate_id=pair.candidate_id,
                        counts_as_reaction=pair.counts_as_reaction,
                    )
                    self._events.append(ev)
                    self._reacted.add(reacted_key)
                    confirmed.append(ev)
            else:
                if is_dissociated(r, pair.r0, self._threshold_fraction):
                    ev = BondEvent(
                        step=step, cycle=cycle,
                        atom_a=pair.idx_a, atom_b=pair.idx_b,
                        event_type='confirmed_dissociation',
                        distance=r, r0=pair.r0,
                        candidate_id=pair.candidate_id,
                        counts_as_reaction=pair.counts_as_reaction,
                    )
                    self._events.append(ev)
                    self._reacted.add(reacted_key)
                    confirmed.append(ev)
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
