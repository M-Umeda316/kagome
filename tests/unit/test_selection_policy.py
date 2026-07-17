"""Tests for the WM-P5a optional stochastic (softmax) selection policy.

Design: specs/decisions.md "2026-07-17: well-mixed 測定モード" item (iv), if
present on your branch. The deterministic policy is the paper-faithful
default (Eq. 7 greedy) and must stay unchanged; 'softmax' is opt-in.
"""
from __future__ import annotations

import numpy as np
import pytest

from kagome.reactive.selection import (
    Candidate,
    audited_selection,
    select_non_overlapping,
)


def _overlapping_candidates() -> list[Candidate]:
    """Four candidates over atoms {0,1,2,3}; every pair overlaps on at least
    one atom with at least one other candidate, so a stochastic run can only
    ever keep one at a time from an overlapping cluster."""
    return [
        Candidate(atom_indices=(0, 1), score=1.0),
        Candidate(atom_indices=(1, 2), score=1.5),
        Candidate(atom_indices=(2, 3), score=2.0),
        Candidate(atom_indices=(0, 3), score=2.5),
    ]


def _disjoint_candidates() -> list[Candidate]:
    """Four mutually non-overlapping candidates, pre-sorted by score ascending
    (matching how score_candidates() would hand them to selection)."""
    return [
        Candidate(atom_indices=(2, 3), score=1.0),
        Candidate(atom_indices=(6, 7), score=2.0),
        Candidate(atom_indices=(4, 5), score=3.0),
        Candidate(atom_indices=(0, 1), score=4.0),
    ]


class TestDefaultPolicyUnchanged:
    """(a) Omitting policy must match today's select_non_overlapping exactly."""

    def test_default_matches_deterministic_helper(self):
        cands = _overlapping_candidates()
        # baseline: no policy kwarg at all (pre-WM-P5a call shape)
        baseline = select_non_overlapping(list(cands))
        explicit = select_non_overlapping(list(cands), policy='deterministic')
        assert [c.atom_indices for c in baseline] == [c.atom_indices for c in explicit]
        assert [c.atom_indices for c in baseline] == [(0, 1), (2, 3)]

    def test_default_audited_selection_unchanged(self):
        cands = _overlapping_candidates()
        selected, decisions = audited_selection(list(cands))
        assert [c.atom_indices for c in selected] == [(0, 1), (2, 3)]
        # deterministic decisions never carry a pool_size
        assert all(d.pool_size is None for d in decisions)

    def test_default_positional_call_still_works(self):
        """Existing positional callers (candidates-only) must keep working."""
        cands = _disjoint_candidates()
        selected = select_non_overlapping(cands)
        assert len(selected) == 4


class TestSoftmaxTinyTemperatureRecoversDeterministic:
    """(b) T -> 0 must reproduce the deterministic ranking."""

    def test_tiny_temperature_matches_deterministic_order(self):
        cands = _overlapping_candidates()
        rng = np.random.default_rng(0)
        det_selected = select_non_overlapping(list(cands), policy='deterministic')
        soft_selected = select_non_overlapping(
            list(cands), policy='softmax', softmax_temperature=1e-12, rng=rng,
        )
        assert [c.atom_indices for c in soft_selected] == \
            [c.atom_indices for c in det_selected]

    def test_zero_temperature_matches_deterministic_order(self):
        cands = _disjoint_candidates()
        rng = np.random.default_rng(1)
        det_selected = select_non_overlapping(list(cands), policy='deterministic')
        soft_selected = select_non_overlapping(
            list(cands), policy='softmax', softmax_temperature=0.0, rng=rng,
        )
        assert [c.atom_indices for c in soft_selected] == \
            [c.atom_indices for c in det_selected]


class TestSoftmaxLargeTemperatureSpread:
    """(c) Large T: over many seeded draws, every viable candidate in an
    overlapping cluster is selected at least once — deterministic via a fixed
    rng sequence, no flakiness."""

    def test_large_temperature_spreads_selection(self):
        winners: set[tuple[int, ...]] = set()
        rng = np.random.default_rng(42)
        for _ in range(200):
            cands = _overlapping_candidates()
            selected = select_non_overlapping(
                cands, policy='softmax', softmax_temperature=1000.0, rng=rng,
            )
            # first pick of the round (the cluster's initial winner)
            winners.add(selected[0].atom_indices)
        expected = {c.atom_indices for c in _overlapping_candidates()}
        assert winners == expected, (
            f'expected every candidate to win at least once over 200 draws, got {winners}'
        )


class TestSoftmaxNonOverlapInvariant:
    """(d) Non-overlap invariant must hold under softmax for a crafted
    overlapping candidate set, across many seeded trials.

    ``_overlapping_candidates`` forms a 4-cycle over atoms {0,1,2,3}
    (0-1, 1-2, 2-3, 0-3): its maximum matching size is 2, so a correct
    non-overlapping selection always keeps exactly 2 of the 4 candidates,
    forming one of the two perfect matchings {(0,1),(2,3)} or {(1,2),(0,3)}.
    """

    _MATCHINGS = (
        {(0, 1), (2, 3)},
        {(1, 2), (0, 3)},
    )

    def test_no_overlap_across_many_trials(self):
        rng = np.random.default_rng(7)
        for _ in range(100):
            cands = _overlapping_candidates()
            selected = select_non_overlapping(
                cands, policy='softmax', softmax_temperature=5.0, rng=rng,
            )
            seen: set[int] = set()
            for c in selected:
                atoms = set(c.atom_indices)
                assert not (atoms & seen), f'overlap detected: {selected}'
                seen |= atoms
            picked = {c.atom_indices for c in selected}
            assert picked in self._MATCHINGS, (
                f'expected one of {self._MATCHINGS}, got {picked}'
            )


class TestSoftmaxReproducibility:
    """(e) Same seed -> identical selections across two independent runs."""

    def test_same_seed_same_selection(self):
        cands = _overlapping_candidates()

        rng_a = np.random.default_rng(123)
        selected_a = select_non_overlapping(
            list(cands), policy='softmax', softmax_temperature=2.0, rng=rng_a,
        )

        rng_b = np.random.default_rng(123)
        selected_b = select_non_overlapping(
            list(cands), policy='softmax', softmax_temperature=2.0, rng=rng_b,
        )

        assert [c.atom_indices for c in selected_a] == \
            [c.atom_indices for c in selected_b]

    def test_same_seed_same_selection_disjoint(self):
        cands = _disjoint_candidates()

        rng_a = np.random.default_rng(9)
        selected_a = select_non_overlapping(
            list(cands), policy='softmax', softmax_temperature=3.0, rng=rng_a,
        )

        rng_b = np.random.default_rng(9)
        selected_b = select_non_overlapping(
            list(cands), policy='softmax', softmax_temperature=3.0, rng=rng_b,
        )

        assert [c.atom_indices for c in selected_a] == \
            [c.atom_indices for c in selected_b]


class TestAuditRecordPolicyFields:
    """(f) Audit decisions carry the pool_size field only for softmax picks."""

    def test_deterministic_decisions_have_no_pool_size(self):
        cands = _overlapping_candidates()
        _, decisions = audited_selection(list(cands))
        assert all(d.pool_size is None for d in decisions)

    def test_softmax_selected_decisions_have_pool_size(self):
        cands = _overlapping_candidates()
        rng = np.random.default_rng(5)
        _, decisions = audited_selection(
            list(cands), policy='softmax', softmax_temperature=2.0, rng=rng,
        )
        selected_decisions = [d for d in decisions if d.selected]
        assert len(selected_decisions) == 2
        for d in selected_decisions:
            assert d.pool_size is not None
            assert d.pool_size >= 1

    def test_softmax_pool_size_matches_viable_count(self):
        """The ring yields exactly 2 picks: the first competes against all 4
        initial candidates, the second is forced (only 1 candidate survives
        once the first pick's atoms are used)."""
        cands = _overlapping_candidates()
        rng = np.random.default_rng(5)
        _, decisions = audited_selection(
            list(cands), policy='softmax', softmax_temperature=2.0, rng=rng,
        )
        selected_decisions = [d for d in decisions if d.selected]
        assert len(selected_decisions) == 2
        assert selected_decisions[0].pool_size == 4
        assert selected_decisions[1].pool_size == 1


class TestSoftmaxValidation:
    """Guard-rail validation for the new kwargs."""

    def test_softmax_requires_temperature(self):
        cands = _disjoint_candidates()
        rng = np.random.default_rng(0)
        with pytest.raises(ValueError):
            select_non_overlapping(cands, policy='softmax', rng=rng)

    def test_softmax_requires_rng(self):
        cands = _disjoint_candidates()
        with pytest.raises(ValueError):
            select_non_overlapping(cands, policy='softmax', softmax_temperature=1.0)

    def test_negative_temperature_rejected(self):
        cands = _disjoint_candidates()
        rng = np.random.default_rng(0)
        with pytest.raises(ValueError):
            select_non_overlapping(
                cands, policy='softmax', softmax_temperature=-1.0, rng=rng,
            )

    def test_unknown_policy_rejected(self):
        cands = _disjoint_candidates()
        with pytest.raises(ValueError):
            select_non_overlapping(cands, policy='bogus')

    def test_empty_candidates_softmax(self):
        rng = np.random.default_rng(0)
        assert select_non_overlapping(
            [], policy='softmax', softmax_temperature=1.0, rng=rng,
        ) == []
