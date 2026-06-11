"""Tests for bond event tracking."""
import numpy as np
import pytest

from src.boost.tdbb import PairBias
from src.reactive.bonds import BondTracker


class TestBondTracker:

    def _make_formation_pair(self, idx_a=0, idx_b=1, r0=1.9):
        return PairBias(idx_a=idx_a, idx_b=idx_b, is_formation=True, r0=r0)

    def _make_dissociation_pair(self, idx_a=0, idx_b=1, r0=1.9):
        return PairBias(idx_a=idx_a, idx_b=idx_b, is_formation=False, r0=r0)

    def test_record_attempts(self):
        tracker = BondTracker()
        positions = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        pair = self._make_formation_pair()

        tracker.record_attempts([pair], positions, step=0, cycle=0)

        events = tracker.events
        assert len(events) == 1
        assert events[0].event_type == 'attempted_formation'
        assert events[0].distance == pytest.approx(3.0)

    def test_confirmed_formation(self):
        tracker = BondTracker(threshold_fraction=1.2)
        pair = self._make_formation_pair(r0=2.0)

        positions_start = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        tracker.record_attempts([pair], positions_start, step=0, cycle=0)

        # atoms moved close: r=2.1 < 1.2 * 2.0 = 2.4
        positions_end = np.array([[0.0, 0.0, 0.0], [2.1, 0.0, 0.0]])
        confirmed = tracker.check_outcomes(positions_end, step=200)

        assert len(confirmed) == 1
        assert confirmed[0].event_type == 'confirmed_formation'
        assert len(tracker.confirmed_formations()) == 1

    def test_unconfirmed_formation(self):
        tracker = BondTracker(threshold_fraction=1.2)
        pair = self._make_formation_pair(r0=2.0)

        positions_start = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        tracker.record_attempts([pair], positions_start, step=0, cycle=0)

        # atoms still far: r=3.5 > 1.2 * 2.0 = 2.4
        positions_end = np.array([[0.0, 0.0, 0.0], [3.5, 0.0, 0.0]])
        confirmed = tracker.check_outcomes(positions_end, step=200)

        assert len(confirmed) == 0
        assert len(tracker.confirmed_formations()) == 0

    def test_confirmed_dissociation(self):
        tracker = BondTracker(threshold_fraction=1.2)
        pair = self._make_dissociation_pair(r0=1.5)

        positions_start = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
        tracker.record_attempts([pair], positions_start, step=0, cycle=0)

        # atoms moved apart: r=2.5 > 1.2 * 1.5 = 1.8
        positions_end = np.array([[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]])
        confirmed = tracker.check_outcomes(positions_end, step=200)

        assert len(confirmed) == 1
        assert confirmed[0].event_type == 'confirmed_dissociation'

    def test_save(self, tmp_path):
        tracker = BondTracker()
        positions = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        pair = self._make_formation_pair()
        tracker.record_attempts([pair], positions, step=0, cycle=0)

        path = tmp_path / 'bonds.jsonl'
        tracker.save(path)

        lines = path.read_text(encoding='utf-8').strip().split('\n')
        assert len(lines) == 1

    def test_multiple_cycles(self):
        tracker = BondTracker(threshold_fraction=1.2)
        pair = self._make_formation_pair(r0=2.0)

        # cycle 0: attempt
        pos = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        tracker.record_attempts([pair], pos, step=0, cycle=0)
        # cycle 0: not confirmed
        tracker.check_outcomes(pos, step=200)

        # cycle 1: attempt again
        tracker.record_attempts([pair], pos, step=400, cycle=1)
        # cycle 1: confirmed
        pos_close = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        confirmed = tracker.check_outcomes(pos_close, step=600)

        assert len(confirmed) == 1
        assert confirmed[0].cycle == 1
        total_formations = tracker.confirmed_formations()
        assert len(total_formations) == 1
