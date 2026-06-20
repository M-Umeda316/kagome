"""Tests for bond event tracking."""
import numpy as np
import pytest

from kagome.boost.tdbb import PairBias
from kagome.reactive.bonds import BondTracker, is_dissociated, is_formed


class TestIsDissociated:

    def test_beyond_threshold(self):
        assert is_dissociated(2.5, r0=1.95, threshold_fraction=1.0) is True

    def test_at_threshold(self):
        assert is_dissociated(1.95, r0=1.95, threshold_fraction=1.0) is False

    def test_below_threshold(self):
        assert is_dissociated(1.5, r0=1.95, threshold_fraction=1.0) is False

    def test_custom_fraction(self):
        assert is_dissociated(2.5, r0=2.0, threshold_fraction=1.2) is True
        assert is_dissociated(2.3, r0=2.0, threshold_fraction=1.2) is False

    def test_default_fraction_is_one(self):
        assert is_dissociated(2.0, r0=1.9) is True
        assert is_dissociated(1.8, r0=1.9) is False

    def test_activation_absolute_threshold(self):
        """Activation uses is_dissociated(r, dissoc_threshold) with absolute 2.5 Å."""
        assert is_dissociated(2.6, r0=2.5) is True
        assert is_dissociated(2.4, r0=2.5) is False


class TestIsFormed:

    def test_below_threshold(self):
        assert is_formed(1.5, r0=1.95, threshold_fraction=1.0) is True

    def test_at_threshold(self):
        assert is_formed(1.95, r0=1.95, threshold_fraction=1.0) is True

    def test_beyond_threshold(self):
        assert is_formed(2.5, r0=1.95, threshold_fraction=1.0) is False

    def test_custom_fraction(self):
        assert is_formed(2.3, r0=2.0, threshold_fraction=1.2) is True
        assert is_formed(2.5, r0=2.0, threshold_fraction=1.2) is False


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

    def test_default_threshold_is_paper_faithful(self):
        # Paper: "60% of the sum of their van der Waals radii" = r0 = lambda * sum(vdW)
        # threshold_fraction=1.0 means threshold = 1.0 * r0 = r0 = 0.6 * sum(vdW)
        tracker = BondTracker()
        assert tracker._threshold_fraction == pytest.approx(1.0), (
            'Default threshold_fraction must be 1.0 (paper: r < 0.6*sum_vdW = r0)'
        )

    # ── check_reactions_during_bias (paper §2.2 step 3) ──────────────

    def test_in_bias_formation_detected(self):
        """Formation pair reacts in-bias when r ≤ threshold·r0."""
        tracker = BondTracker(threshold_fraction=1.0)
        pair = self._make_formation_pair(r0=2.0)
        positions = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
        events = tracker.check_reactions_during_bias(
            [pair], positions, step=10, cycle=0,
        )
        assert len(events) == 1
        assert events[0].event_type == 'confirmed_formation'
        assert events[0].distance == pytest.approx(1.5)
        assert events[0].step == 10
        assert events[0].cycle == 0

    def test_in_bias_dissociation_detected(self):
        """Dissociation pair reacts in-bias when r > threshold·r0."""
        tracker = BondTracker(threshold_fraction=1.0)
        pair = self._make_dissociation_pair(r0=1.5)
        positions = np.array([[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]])
        events = tracker.check_reactions_during_bias(
            [pair], positions, step=10, cycle=0,
        )
        assert len(events) == 1
        assert events[0].event_type == 'confirmed_dissociation'

    def test_in_bias_no_duplicate_detection(self):
        """Once confirmed in-bias, a second call does not re-detect the same pair."""
        tracker = BondTracker(threshold_fraction=1.0)
        pair = self._make_formation_pair(r0=2.0)
        positions = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])

        first = tracker.check_reactions_during_bias(
            [pair], positions, step=10, cycle=0,
        )
        assert len(first) == 1

        second = tracker.check_reactions_during_bias(
            [pair], positions, step=20, cycle=0,
        )
        assert len(second) == 0

    def test_check_outcomes_skips_in_bias_confirmed(self):
        """check_outcomes does not double-count pairs already confirmed during bias."""
        tracker = BondTracker(threshold_fraction=1.0)
        pair = self._make_formation_pair(r0=2.0)
        positions_close = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])

        tracker.record_attempts([pair], positions_close, step=0, cycle=0)
        in_bias = tracker.check_reactions_during_bias(
            [pair], positions_close, step=5, cycle=0,
        )
        assert len(in_bias) == 1

        outcomes = tracker.check_outcomes(positions_close, step=200)
        assert len(outcomes) == 0
        assert len(tracker.confirmed_formations()) == 1

    def test_in_bias_formation_not_detected_when_far(self):
        """Formation pair does NOT react when r > threshold·r0."""
        tracker = BondTracker(threshold_fraction=1.0)
        pair = self._make_formation_pair(r0=2.0)
        positions = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        events = tracker.check_reactions_during_bias(
            [pair], positions, step=10, cycle=0,
        )
        assert len(events) == 0

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
