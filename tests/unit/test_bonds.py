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
    # Return value is now the list of fired candidate_ids; the per-pair
    # tentative crossings are recorded to tracker.events as the audit trail.

    def test_in_bias_formation_tentative(self):
        """Formation pair records a tentative (not confirmed) audit event and,
        as a singleton candidate, fires."""
        tracker = BondTracker(threshold_fraction=1.0)
        pair = self._make_formation_pair(r0=2.0)
        positions = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
        fired = tracker.check_reactions_during_bias(
            [pair], positions, step=10, cycle=0,
        )
        assert fired == [-1]  # singleton (candidate_id<0) trigger crossing fires
        events = tracker.events
        assert len(events) == 1
        assert events[0].event_type == 'tentative_formation'
        assert events[0].distance == pytest.approx(1.5)
        assert events[0].step == 10
        assert events[0].cycle == 0
        assert len(tracker.confirmed_formations()) == 0

    def test_in_bias_dissociation_tentative(self):
        """Dissociation pair records a tentative audit event and fires (singleton)."""
        tracker = BondTracker(threshold_fraction=1.0)
        pair = self._make_dissociation_pair(r0=1.5)
        positions = np.array([[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]])
        fired = tracker.check_reactions_during_bias(
            [pair], positions, step=10, cycle=0,
        )
        assert fired == [-1]
        events = tracker.events
        assert len(events) == 1
        assert events[0].event_type == 'tentative_dissociation'
        assert len(tracker.confirmed_dissociations()) == 0

    def test_in_bias_no_duplicate_audit_record(self):
        """The per-pair audit record is written once; a second in-bias call does
        not append a duplicate tentative event."""
        tracker = BondTracker(threshold_fraction=1.0)
        pair = self._make_formation_pair(r0=2.0)
        positions = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])

        tracker.check_reactions_during_bias([pair], positions, step=10, cycle=0)
        assert len(tracker.events) == 1

        tracker.check_reactions_during_bias([pair], positions, step=20, cycle=0)
        assert len(tracker.events) == 1  # no new audit record

    def test_tentative_confirmed_after_unbiased_relaxation(self):
        """Tentative in-bias detection → confirmed by check_outcomes when still close."""
        tracker = BondTracker(threshold_fraction=1.0)
        pair = self._make_formation_pair(r0=2.0)
        positions_close = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])

        tracker.record_attempts([pair], positions_close, step=0, cycle=0)
        fired = tracker.check_reactions_during_bias(
            [pair], positions_close, step=5, cycle=0,
        )
        assert fired == [-1]
        assert tracker.events[-1].event_type == 'tentative_formation'

        outcomes = tracker.check_outcomes(positions_close, step=200)
        assert len(outcomes) == 1
        assert outcomes[0].event_type == 'confirmed_formation'
        assert len(tracker.confirmed_formations()) == 1

    def test_tentative_not_confirmed_when_drifted_apart(self):
        """Tentative in-bias detection → NOT confirmed if pair drifts apart."""
        tracker = BondTracker(threshold_fraction=1.0)
        pair = self._make_formation_pair(r0=2.0)
        positions_close = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
        positions_far = np.array([[0.0, 0.0, 0.0], [3.5, 0.0, 0.0]])

        tracker.record_attempts([pair], positions_close, step=0, cycle=0)
        fired = tracker.check_reactions_during_bias(
            [pair], positions_close, step=5, cycle=0,
        )
        assert fired == [-1]

        outcomes = tracker.check_outcomes(positions_far, step=200)
        assert len(outcomes) == 0
        assert len(tracker.confirmed_formations()) == 0

    def test_in_bias_formation_not_detected_when_far(self):
        """Formation pair neither fires nor records when r > threshold·r0."""
        tracker = BondTracker(threshold_fraction=1.0)
        pair = self._make_formation_pair(r0=2.0)
        positions = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        fired = tracker.check_reactions_during_bias(
            [pair], positions, step=10, cycle=0,
        )
        assert fired == []
        assert len(tracker.events) == 0

    # ── S2/W3: pair_dists reuse is bit-identical (specs/decisions.md 2026-07-06) ──

    def test_pair_dists_reuse_matches_recompute(self):
        """Passing total_bias_fast pair_dists yields identical in-bias events.

        The biased-phase hot path reuses the distances total_bias_fast already
        computed for the same coordinates instead of recomputing the minimum
        image inside check_reactions_during_bias.  The two paths must produce
        bit-identical distances, identical fired candidate_ids and identical
        audit events (S2/W3).
        """
        from kagome.boost.tdbb import BoostState, TDBBParams, total_bias_fast
        from kagome.geometry import validated_box

        pairs = [
            self._make_formation_pair(idx_a=0, idx_b=1, r0=2.0),
            self._make_formation_pair(idx_a=2, idx_b=3, r0=2.0),
        ]
        # idx 2/3 straddle the periodic boundary so the minimum image matters.
        positions = np.array([
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [0.3, 0.0, 0.0],
            [9.7, 0.0, 0.0],
        ])
        cell = np.diag([10.0, 10.0, 10.0]).astype(float)
        box = validated_box(cell)

        boost = BoostState()
        boost.advance(1.0, 250.0, 250.0)
        _, _, pair_dists = total_bias_fast(
            pairs, positions, boost, TDBBParams(), box,
        )

        tracker_recompute = BondTracker(threshold_fraction=1.0)
        fired_recompute = tracker_recompute.check_reactions_during_bias(
            pairs, positions, step=10, cycle=0, cell=cell,
        )
        tracker_reuse = BondTracker(threshold_fraction=1.0)
        fired_reuse = tracker_reuse.check_reactions_during_bias(
            pairs, positions, step=10, cycle=0, cell=cell,
            pair_dists=pair_dists,
        )

        assert fired_reuse == fired_recompute
        events_recompute = tracker_recompute.events
        events_reuse = tracker_reuse.events
        assert len(events_reuse) == len(events_recompute)
        for ev_r, ev_c in zip(events_reuse, events_recompute):
            assert ev_r.atom_a == ev_c.atom_a
            assert ev_r.atom_b == ev_c.atom_b
            assert ev_r.event_type == ev_c.event_type
            # bit-identical distance (same coords, same MIC formula)
            assert ev_r.distance == ev_c.distance
        # both pairs are within r0 under the minimum image → 2 tentative events
        assert len(events_reuse) == 2

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


class TestConjunctiveReactionEvent:
    """F1' (paper §2.2 step 3-4): a candidate's reaction event fires / confirms
    only when ALL its trigger pairs satisfy simultaneously (formation bonded AND
    every dissociation broken). The bias-only water pair (is_trigger=False) never
    gates the event; it may or may not close."""

    @staticmethod
    def _nylon_pairs(cid: int = 0) -> list[PairBias]:
        """Nylon-like candidate: atoms 0=N, 1=C, 2=H(on N), 3=OH(on C).

        Trigger pairs: amide N-C formation (counted), N-H dissociation, C-OH
        dissociation. Bias-only pair: H-OH water formation (is_trigger=False,
        counts_as_reaction=False).
        """
        return [
            PairBias(idx_a=0, idx_b=1, is_formation=True, r0=2.0,
                     candidate_id=cid, is_trigger=True, counts_as_reaction=True),
            PairBias(idx_a=0, idx_b=2, is_formation=False, r0=1.5,
                     candidate_id=cid, is_trigger=True),
            PairBias(idx_a=1, idx_b=3, is_formation=False, r0=1.5,
                     candidate_id=cid, is_trigger=True),
            PairBias(idx_a=2, idx_b=3, is_formation=True, r0=1.5,
                     candidate_id=cid, is_trigger=False, counts_as_reaction=False),
        ]

    # ── firing (check_reactions_during_bias) ──────────────────────────

    def test_does_not_fire_when_only_formation_satisfied(self):
        """Amide formed but leaving groups still bonded → no fire, even though
        the bias-only water pair is closed."""
        tracker = BondTracker(threshold_fraction=1.0)
        positions = np.array([
            [0.0, 0.0, 0.0],   # N
            [1.8, 0.0, 0.0],   # C  (N-C = 1.8 <= 2.0 : formed)
            [0.0, 1.0, 0.0],   # H  (N-H = 1.0 <= 1.5 : still bonded)
            [0.5, 1.0, 0.0],   # OH (H-OH = 0.5 : water closed)
        ])
        fired = tracker.check_reactions_during_bias(
            self._nylon_pairs(), positions, step=1, cycle=0,
        )
        assert fired == []

    def test_fires_only_when_formation_and_both_dissociations_satisfied(self):
        """Amide formed AND both leaving groups dissociated → fires, even though
        the bias-only water pair is still open."""
        tracker = BondTracker(threshold_fraction=1.0)
        positions = np.array([
            [0.0, 0.0, 0.0],   # N
            [1.8, 0.0, 0.0],   # C  (N-C = 1.8 : formed)
            [0.0, 3.0, 0.0],   # H  (N-H = 3.0 > 1.5 : dissociated)
            [1.8, 3.0, 0.0],   # OH (C-OH = 3.0 dissociated; H-OH = 1.8 : open)
        ])
        fired = tracker.check_reactions_during_bias(
            self._nylon_pairs(), positions, step=1, cycle=0,
        )
        assert fired == [0]

    # ── confirmation (check_outcomes) ─────────────────────────────────

    @staticmethod
    def _start_positions() -> np.ndarray:
        return np.array([
            [0.0, 0.0, 0.0], [1.8, 0.0, 0.0],
            [0.0, 1.0, 0.0], [1.8, 1.0, 0.0],
        ])

    def test_confirms_whole_candidate_when_conjunction_holds(self):
        """All triggers satisfied post-relaxation and water closed: the whole
        candidate commits — counted amide formation, both dissociations, and the
        non-counted water formation."""
        tracker = BondTracker(threshold_fraction=1.0)
        pairs = self._nylon_pairs()
        tracker.record_attempts(pairs, self._start_positions(), step=0, cycle=0)
        relaxed = np.array([
            [0.0, 0.0, 0.0],   # N
            [1.8, 0.0, 0.0],   # C  (amide formed)
            [0.0, 3.0, 0.0],   # H  (N-H = 3.0 dissociated)
            [1.0, 3.0, 0.0],   # OH (C-OH = 3.1 dissociated; H-OH = 1.0 closed)
        ])
        confirmed = tracker.check_outcomes(relaxed, step=200)
        sig = {(e.event_type, e.atom_a, e.atom_b, e.counts_as_reaction)
               for e in confirmed}
        assert ('confirmed_formation', 0, 1, True) in sig       # amide (counted)
        assert ('confirmed_dissociation', 0, 2, True) in sig    # N-H
        assert ('confirmed_dissociation', 1, 3, True) in sig    # C-OH
        assert ('confirmed_formation', 2, 3, False) in sig      # water (uncounted)
        assert len(confirmed) == 4
        assert len(tracker.confirmed_formations()) == 2
        assert len(tracker.confirmed_dissociations()) == 2

    def test_water_pair_confirms_only_if_it_closed(self):
        """Triggers satisfied but water still open: amide + dissociations commit,
        water formation is NOT emitted."""
        tracker = BondTracker(threshold_fraction=1.0)
        pairs = self._nylon_pairs()
        tracker.record_attempts(pairs, self._start_positions(), step=0, cycle=0)
        relaxed = np.array([
            [0.0, 0.0, 0.0],
            [1.8, 0.0, 0.0],
            [0.0, 3.0, 0.0],   # N-H = 3.0 dissociated
            [1.8, 3.0, 0.0],   # C-OH = 3.0 dissociated; H-OH = 1.8 still open
        ])
        confirmed = tracker.check_outcomes(relaxed, step=200)
        sig = {(e.event_type, e.atom_a, e.atom_b) for e in confirmed}
        assert ('confirmed_formation', 0, 1) in sig
        assert ('confirmed_formation', 2, 3) not in sig  # water did not close
        assert len(confirmed) == 3

    def test_confirms_nothing_when_only_some_triggers_hold(self):
        """One leaving group dissociated but the other still bonded → conjunction
        fails → NOTHING confirms (no spurious lone dissociation)."""
        tracker = BondTracker(threshold_fraction=1.0)
        pairs = self._nylon_pairs()
        tracker.record_attempts(pairs, self._start_positions(), step=0, cycle=0)
        relaxed = np.array([
            [0.0, 0.0, 0.0],
            [1.8, 0.0, 0.0],
            [0.0, 3.0, 0.0],   # N-H = 3.0 dissociated
            [1.8, 1.0, 0.0],   # C-OH = 1.0 STILL bonded
        ])
        confirmed = tracker.check_outcomes(relaxed, step=200)
        assert confirmed == []
        assert tracker.confirmed_formations() == []
        assert tracker.confirmed_dissociations() == []

    # ── vinyl (single-formation candidate) unchanged ──────────────────

    def test_vinyl_single_formation_candidate_fires_and_confirms(self):
        tracker = BondTracker(threshold_fraction=1.0)
        pair = PairBias(idx_a=0, idx_b=1, is_formation=True, r0=2.0,
                        candidate_id=0, is_trigger=True)
        close = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
        assert tracker.check_reactions_during_bias(
            [pair], close, step=1, cycle=0) == [0]
        tracker.record_attempts([pair], close, step=0, cycle=0)
        confirmed = tracker.check_outcomes(close, step=200)
        assert len(confirmed) == 1
        assert confirmed[0].event_type == 'confirmed_formation'

    # ── activation / legacy (candidate_id < 0) unchanged ──────────────

    def test_legacy_negative_candidate_id_fires_and_confirms_per_pair(self):
        """A candidate_id<0 dissociation pair (activation style) fires and
        confirms independently, on its own crossing."""
        tracker = BondTracker(threshold_fraction=1.0)
        diss = PairBias(idx_a=0, idx_b=1, is_formation=False, r0=1.5,
                        candidate_id=-1)
        bonded = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])  # not dissociated
        broken = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])  # dissociated
        assert tracker.check_reactions_during_bias(
            [diss], bonded, step=1, cycle=0) == []
        assert tracker.check_reactions_during_bias(
            [diss], broken, step=2, cycle=0) == [-1]
        tracker.record_attempts([diss], broken, step=0, cycle=0)
        confirmed = tracker.check_outcomes(broken, step=200)
        assert len(confirmed) == 1
        assert confirmed[0].event_type == 'confirmed_dissociation'
