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

    def __init__(self, threshold_fraction: float = 1.2) -> None:
        self._threshold_fraction = threshold_fraction
        self._events: list[BondEvent] = []
        self._pending: list[tuple[PairBias, int]] = []

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
            r_vec = minimum_image(
                positions[pair.idx_b] - positions[pair.idx_a], cell,
            )
            r = float(np.linalg.norm(r_vec))
            if pair.is_formation:
                threshold = self._threshold_fraction * pair.r0
                if r <= threshold:
                    ev = BondEvent(
                        step=step, cycle=cycle,
                        atom_a=pair.idx_a, atom_b=pair.idx_b,
                        event_type='confirmed_formation',
                        distance=r, r0=pair.r0,
                    )
                    self._events.append(ev)
                    confirmed.append(ev)
            else:
                threshold = self._threshold_fraction * pair.r0
                if r > threshold:
                    ev = BondEvent(
                        step=step, cycle=cycle,
                        atom_a=pair.idx_a, atom_b=pair.idx_b,
                        event_type='confirmed_dissociation',
                        distance=r, r0=pair.r0,
                    )
                    self._events.append(ev)
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
