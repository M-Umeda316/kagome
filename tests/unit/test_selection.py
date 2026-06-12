"""Tests for candidate selection (Eq. 7) and non-overlapping filter."""
import numpy as np
import pytest

from src.reactive.groups import PairSpec, ReactiveGroup, ReactionTemplate
from src.reactive.selection import (
    Candidate,
    find_candidates,
    score_candidates,
    select_non_overlapping,
)


def _make_template() -> ReactionTemplate:
    return ReactionTemplate(
        name='test_reaction',
        groups=['I', 'J'],
        pairs=[PairSpec(group_a='I', group_b='J', is_formation=True, r_min=0.5, r_max=3.0)],
    )


def _make_groups() -> dict[str, ReactiveGroup]:
    return {
        'I': ReactiveGroup(label='I', atom_indices=[0, 1]),
        'J': ReactiveGroup(label='J', atom_indices=[2, 3]),
    }


class TestFindCandidates:

    def test_within_bounds(self):
        positions = np.array([
            [0.0, 0.0, 0.0],  # atom 0 (I)
            [5.0, 0.0, 0.0],  # atom 1 (I)
            [1.0, 0.0, 0.0],  # atom 2 (J) — within range of atom 0
            [6.0, 0.0, 0.0],  # atom 3 (J) — within range of atom 1
        ])
        template = _make_template()
        groups = _make_groups()

        cands = find_candidates(template, groups, positions)

        indices_set = {c.atom_indices for c in cands}
        assert (0, 2) in indices_set
        assert (1, 3) in indices_set

    def test_out_of_range_excluded(self):
        positions = np.array([
            [0.0, 0.0, 0.0],
            [100.0, 0.0, 0.0],
            [50.0, 0.0, 0.0],
            [50.0, 50.0, 0.0],
        ])
        template = _make_template()
        groups = _make_groups()

        cands = find_candidates(template, groups, positions)

        assert len(cands) == 0

    def test_empty_groups(self):
        positions = np.array([[0.0, 0.0, 0.0]])
        template = _make_template()
        groups = {
            'I': ReactiveGroup(label='I', atom_indices=[]),
            'J': ReactiveGroup(label='J', atom_indices=[]),
        }
        cands = find_candidates(template, groups, positions)
        assert len(cands) == 0


class TestScoreCandidates:

    def test_sorted_ascending(self):
        positions = np.array([
            [0.0, 0.0, 0.0],
            [5.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [4.5, 0.0, 0.0],
        ])
        template = _make_template()
        groups = _make_groups()

        cands = find_candidates(template, groups, positions)
        scored = score_candidates(cands, template, positions)

        scores = [c.score for c in scored]
        assert scores == sorted(scores)

    def test_score_equals_distance(self):
        positions = np.array([
            [0.0, 0.0, 0.0],
            [5.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
        ])
        template = _make_template()
        groups = _make_groups()

        cands = find_candidates(template, groups, positions)
        scored = score_candidates(cands, template, positions)

        for c in scored:
            expected = np.linalg.norm(
                positions[c.atom_indices[1]] - positions[c.atom_indices[0]]
            )
            assert c.score == pytest.approx(expected)


class TestSelectNonOverlapping:

    def test_no_overlap(self):
        c1 = Candidate(atom_indices=(0, 2), pair_distances={}, score=1.0)
        c2 = Candidate(atom_indices=(1, 3), pair_distances={}, score=2.0)

        selected = select_non_overlapping([c1, c2])
        assert len(selected) == 2

    def test_overlap_rejected(self):
        c1 = Candidate(atom_indices=(0, 2), pair_distances={}, score=1.0)
        c2 = Candidate(atom_indices=(0, 3), pair_distances={}, score=2.0)

        selected = select_non_overlapping([c1, c2])
        assert len(selected) == 1
        assert selected[0] is c1

    def test_greedy_order_matters(self):
        c1 = Candidate(atom_indices=(0, 1), pair_distances={}, score=1.0)
        c2 = Candidate(atom_indices=(1, 2), pair_distances={}, score=2.0)
        c3 = Candidate(atom_indices=(2, 3), pair_distances={}, score=3.0)

        selected = select_non_overlapping([c1, c2, c3])
        assert len(selected) == 2
        assert selected[0] is c1
        assert selected[1] is c3

    def test_empty_input(self):
        assert select_non_overlapping([]) == []


class TestFourGroupTemplate:
    """Test with a 4-group template matching the paper's (I,J,K,L) pattern."""

    def test_four_group_candidates(self):
        template = ReactionTemplate(
            name='vinyl',
            groups=['I', 'J', 'K', 'L'],
            pairs=[
                PairSpec('I', 'J', is_formation=True, r_min=0.5, r_max=3.0),
                PairSpec('I', 'K', is_formation=False, r_min=0.5, r_max=3.0),
                PairSpec('J', 'L', is_formation=False, r_min=0.5, r_max=3.0),
            ],
        )
        groups = {
            'I': ReactiveGroup('I', [0]),
            'J': ReactiveGroup('J', [1]),
            'K': ReactiveGroup('K', [2]),
            'L': ReactiveGroup('L', [3]),
        }
        positions = np.array([
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 1.0, 0.0],
        ])

        cands = find_candidates(template, groups, positions)
        assert len(cands) == 1
        assert cands[0].atom_indices == (0, 1, 2, 3)

    def test_reversed_group_order_in_pair(self):
        """PairSpec with group_b listed before group_a in template.groups
        must still enforce distance constraints."""
        template = ReactionTemplate(
            name='reversed',
            groups=['X', 'Y'],
            pairs=[PairSpec(group_a='Y', group_b='X', is_formation=True,
                            r_min=0.5, r_max=3.0)],
        )
        groups = {
            'X': ReactiveGroup('X', [0]),
            'Y': ReactiveGroup('Y', [1]),
        }
        positions = np.array([
            [0.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],  # far apart — should be excluded by r_max=3.0
        ])
        cands = find_candidates(template, groups, positions)
        assert len(cands) == 0

    def test_reversed_group_order_accepts_valid(self):
        """Reversed group order should accept pairs within range."""
        template = ReactionTemplate(
            name='reversed',
            groups=['X', 'Y'],
            pairs=[PairSpec(group_a='Y', group_b='X', is_formation=True,
                            r_min=0.5, r_max=3.0)],
        )
        groups = {
            'X': ReactiveGroup('X', [0]),
            'Y': ReactiveGroup('Y', [1]),
        }
        positions = np.array([
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],  # within range
        ])
        cands = find_candidates(template, groups, positions)
        assert len(cands) == 1
