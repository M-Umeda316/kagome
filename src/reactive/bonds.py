"""Bond event tracking for reaction monitoring.

Records formation/dissociation attempts during biased phases and
confirms outcomes after unbiased relaxation.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from src.boost.tdbb import PairBias
from src.geometry import minimum_image


def is_formed(r: float, r0: float, threshold_fraction: float = 1.0) -> bool:
    """True when distance *r* indicates bond formation (r <= threshold)."""
    return r <= threshold_fraction * r0


def is_dissociated(r: float, r0: float, threshold_fraction: float = 1.0) -> bool:
    """True when distance *r* indicates bond dissociation (r > threshold)."""
    return r > threshold_fraction * r0


@dataclass
class BondEvent:
    step: int
    cycle: int
    atom_a: int
    atom_b: int
    event_type: str  # attempted_formation, attempted_dissociation, confirmed_formation, confirmed_dissociation
    distance: float
    r0: float = 0.0


class BondTracker:
    """Tracks bond formation and dissociation events across cycles."""

    def __init__(self, threshold_fraction: float = 1.0) -> None:
        self._threshold_fraction = threshold_fraction
        self._events: list[BondEvent] = []
        self._pending: list[tuple[PairBias, int]] = []
        # Pairs already confirmed (min,max index) so in-phase detection and the
        # end-of-unbiased check never double-count the same reaction.
        self._reacted: set[tuple[int, int]] = set()

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
        """Detect reaction events DURING the biased phase (paper §2.2 step 3).

        A formation pair reacts when its separation falls below the vdW bonding
        threshold (r ≤ threshold_fraction·r0 = 0.6·Σr_vdw); a dissociation pair
        reacts when it rises above it. Newly reacted pairs are recorded once and
        returned so the caller can end the biased segment (run-until-reaction).
        """
        newly: list[BondEvent] = []
        for pair in pairs:
            key = self._key(pair.idx_a, pair.idx_b)
            if key in self._reacted:
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
            etype = ('confirmed_formation' if pair.is_formation
                     else 'confirmed_dissociation')
            ev = BondEvent(
                step=step, cycle=cycle,
                atom_a=pair.idx_a, atom_b=pair.idx_b,
                event_type=etype, distance=r, r0=pair.r0,
            )
            self._events.append(ev)
            self._reacted.add(key)
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
            if self._key(pair.idx_a, pair.idx_b) in self._reacted:
                continue  # already confirmed during the biased phase
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
                    )
                    self._events.append(ev)
                    self._reacted.add(self._key(pair.idx_a, pair.idx_b))
                    confirmed.append(ev)
            else:
                if is_dissociated(r, pair.r0, self._threshold_fraction):
                    ev = BondEvent(
                        step=step, cycle=cycle,
                        atom_a=pair.idx_a, atom_b=pair.idx_b,
                        event_type='confirmed_dissociation',
                        distance=r, r0=pair.r0,
                    )
                    self._events.append(ev)
                    self._reacted.add(self._key(pair.idx_a, pair.idx_b))
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
