"""Tests for candidate selection (Eq. 7) and non-overlapping filter."""
import numpy as np
import pytest

from kagome.reactive.groups import PairSpec, ReactiveGroup, ReactionTemplate
from kagome.reactive.selection import (
    Candidate,
    audited_selection,
    find_candidates,
    score_candidates,
    select_non_overlapping,
)


class TestAuditedSelection:
    """RF18: greedy selection must also report why each candidate was dropped."""

    def test_audit_matches_select_non_overlapping(self):
        cands = [
            Candidate(atom_indices=(0, 2), score=1.0),
            Candidate(atom_indices=(0, 3), score=2.0),  # clashes on atom 0
            Candidate(atom_indices=(1, 4), score=3.0),
        ]
        selected, decisions = audited_selection(cands)
        assert [c.atom_indices for c in selected] == [(0, 2), (1, 4)]
        # wrapper returns the same selected set
        assert [c.atom_indices for c in select_non_overlapping(cands)] == \
            [(0, 2), (1, 4)]
        assert len(decisions) == 3

    def test_audit_records_reasons(self):
        cands = [
            Candidate(atom_indices=(0, 2), score=1.0),
            Candidate(atom_indices=(0, 3), score=2.0),
        ]
        _, decisions = audited_selection(cands)
        assert decisions[0].selected and decisions[0].reason == 'selected'
        assert not decisions[1].selected
        assert decisions[1].reason == 'overlap:[0]'


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


class TestFindCandidatesPBC:
    """#3: the distance fast-path (validated_box once + minimum_image_fast) must
    be numerically identical to the previous per-pair minimum_image call, and
    must still reject triclinic cells."""

    def test_fast_path_matches_minimum_image(self):
        from kagome.geometry import minimum_image

        # atom 2 (J) sits just across the periodic boundary from atom 0 (I), so
        # the minimum image (not the raw separation) decides membership/score.
        positions = np.array([
            [0.3, 0.0, 0.0],   # 0 (I)
            [5.0, 0.0, 0.0],   # 1 (I)
            [9.8, 0.0, 0.0],   # 2 (J) — MIC distance to atom 0 is 0.5, raw is 9.5
            [6.0, 0.0, 0.0],   # 3 (J)
        ])
        cell = np.diag([10.0, 10.0, 10.0]).astype(float)
        template = _make_template()
        groups = _make_groups()

        cands = find_candidates(template, groups, positions, cell=cell)

        # Reference scores via the OLD full-validation minimum_image path.
        def ref_dist(i, j):
            return float(np.linalg.norm(
                minimum_image(positions[j] - positions[i], cell)))

        by_atoms = {c.atom_indices: c.score for c in cands}
        # (0, 2) is only a candidate because the minimum image brings them to 0.5.
        assert (0, 2) in by_atoms
        for atoms, score in by_atoms.items():
            # single-pair template => score is just the I-J distance
            assert score == ref_dist(atoms[0], atoms[1])  # bit-identical

    def test_triclinic_cell_raises(self):
        positions = np.array([
            [0.0, 0.0, 0.0], [5.0, 0.0, 0.0],
            [1.0, 0.0, 0.0], [6.0, 0.0, 0.0],
        ])
        cell = np.array([[10.0, 1.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]])
        with pytest.raises(ValueError, match='orthorhombic'):
            find_candidates(_make_template(), _make_groups(), positions, cell=cell)


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
        scored = score_candidates(cands)

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
        scored = score_candidates(cands)

        for c in scored:
            expected = np.linalg.norm(
                positions[c.atom_indices[1]] - positions[c.atom_indices[0]]
            )
            assert c.score == pytest.approx(expected)


class TestSelectNonOverlapping:

    def test_no_overlap(self):
        c1 = Candidate(atom_indices=(0, 2), score=1.0)
        c2 = Candidate(atom_indices=(1, 3), score=2.0)

        selected = select_non_overlapping([c1, c2])
        assert len(selected) == 2

    def test_overlap_rejected(self):
        c1 = Candidate(atom_indices=(0, 2), score=1.0)
        c2 = Candidate(atom_indices=(0, 3), score=2.0)

        selected = select_non_overlapping([c1, c2])
        assert len(selected) == 1
        assert selected[0] is c1

    def test_greedy_order_matters(self):
        c1 = Candidate(atom_indices=(0, 1), score=1.0)
        c2 = Candidate(atom_indices=(1, 2), score=2.0)
        c3 = Candidate(atom_indices=(2, 3), score=3.0)

        selected = select_non_overlapping([c1, c2, c3])
        assert len(selected) == 2
        assert selected[0] is c1
        assert selected[1] is c3

    def test_empty_input(self):
        assert select_non_overlapping([]) == []


class TestScorePairFlag:
    """RF5: score_pair=False excludes pairs from identification and scoring."""

    def _make_nylon_template(self):
        """Nylon-like template with k-l as bias-only (score_pair=False)."""
        return ReactionTemplate(
            name='nylon_test',
            groups=['amine_N', 'carboxyl_C', 'amine_H', 'carboxyl_OH'],
            pairs=[
                PairSpec('amine_N', 'carboxyl_C', is_formation=True,
                         r_min=3.0, r_max=6.0),
                PairSpec('amine_N', 'amine_H', is_formation=False,
                         r_min=0.0, r_max=3.0),
                PairSpec('carboxyl_C', 'carboxyl_OH', is_formation=False,
                         r_min=0.0, r_max=3.0),
                PairSpec('amine_H', 'carboxyl_OH', is_formation=True,
                         r_min=0.0, r_max=100.0, score_pair=False),
            ],
        )

    def _make_nylon_groups_and_positions(self):
        groups = {
            'amine_N':     ReactiveGroup('amine_N', [0]),
            'carboxyl_C':  ReactiveGroup('carboxyl_C', [1]),
            'amine_H':     ReactiveGroup('amine_H', [2]),
            'carboxyl_OH': ReactiveGroup('carboxyl_OH', [3]),
        }
        positions = np.array([
            [0.0, 0.0, 0.0],   # amine_N (i)
            [4.0, 0.0, 0.0],   # carboxyl_C (j) — 4.0 Å from N, in [3,6]
            [0.0, 1.0, 0.0],   # amine_H (k) — 1.0 Å from N, in [0,3]
            [4.0, 1.0, 0.0],   # carboxyl_OH (l) — 1.0 Å from C, in [0,3]
        ])
        return groups, positions

    def test_nylon_score_is_three_terms(self):
        """Nylon score = r_ij + r_ik + r_jl (3 terms), not 4."""
        template = self._make_nylon_template()
        groups, positions = self._make_nylon_groups_and_positions()

        cands = find_candidates(template, groups, positions)
        assert len(cands) == 1

        scored = score_candidates(cands)
        r_ij = np.linalg.norm(positions[1] - positions[0])  # N-C
        r_ik = np.linalg.norm(positions[2] - positions[0])  # N-H
        r_jl = np.linalg.norm(positions[3] - positions[1])  # C-OH
        expected_3term = r_ij + r_ik + r_jl

        assert scored[0].score == pytest.approx(expected_3term)

    def test_score_pair_false_not_in_candidate_filter(self):
        """score_pair=False pair's distance window does not gate candidates."""
        template = ReactionTemplate(
            name='test_bias_only',
            groups=['A', 'B', 'C'],
            pairs=[
                PairSpec('A', 'B', is_formation=True, r_min=0.5, r_max=5.0),
                PairSpec('B', 'C', is_formation=True, r_min=0.0, r_max=1.0,
                         score_pair=False),
            ],
        )
        groups = {
            'A': ReactiveGroup('A', [0]),
            'B': ReactiveGroup('B', [1]),
            'C': ReactiveGroup('C', [2]),
        }
        positions = np.array([
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [20.0, 0.0, 0.0],  # B-C distance=18 Å, exceeds r_max=1.0
        ])
        cands = find_candidates(template, groups, positions)
        assert len(cands) == 1

    def test_vinyl_score_unchanged(self):
        """Vinyl (all score_pair=True) score = sum of all 3 pair distances."""
        template = ReactionTemplate(
            name='vinyl',
            groups=['I', 'J', 'K', 'L'],
            pairs=[
                PairSpec('I', 'J', is_formation=True, r_min=0.5, r_max=5.0),
                PairSpec('I', 'K', is_formation=False, r_min=0.0, r_max=3.0,
                         constraint_only=True),
                PairSpec('J', 'L', is_formation=False, r_min=0.0, r_max=3.0,
                         constraint_only=True),
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
            [0.0, 1.0, 0.0],
            [1.5, 1.0, 0.0],
        ])
        cands = find_candidates(template, groups, positions)
        scored = score_candidates(cands)

        r_ij = np.linalg.norm(positions[1] - positions[0])
        r_ik = np.linalg.norm(positions[2] - positions[0])
        r_jl = np.linalg.norm(positions[3] - positions[1])
        assert scored[0].score == pytest.approx(r_ij + r_ik + r_jl)


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


class TestTopologyGuardWithConstraints:
    """L6 topology guard must NOT reject candidates based on constraint_only pairs.

    Regression: constraint_only pairs (e.g. radical_C-chain_C) whose atoms are
    bonded in the initial topology were incorrectly included in
    formation_group_pairs, causing ALL candidates to be rejected at cycle 0.
    """

    def test_constraint_bonded_pair_does_not_block_candidates(self):
        """4-group vinyl template: radical_C-chain_C constraint pair is bonded
        in the initial topology, but candidates must still be found."""
        from kagome.reactive.topology import BondTopology

        template = ReactionTemplate(
            name='vinyl_test',
            groups=['radical_C', 'vinyl_alpha_C', 'chain_C', 'vinyl_beta_C'],
            pairs=[
                PairSpec('radical_C', 'vinyl_alpha_C', is_formation=True,
                         r_min=1.0, r_max=6.0),
                PairSpec('radical_C', 'chain_C', is_formation=True,
                         r_min=0.0, r_max=3.0, constraint_only=True),
                PairSpec('vinyl_alpha_C', 'vinyl_beta_C', is_formation=True,
                         r_min=0.0, r_max=3.0, constraint_only=True),
            ],
        )
        groups = {
            'radical_C':      ReactiveGroup('radical_C', [0]),
            'vinyl_alpha_C':  ReactiveGroup('vinyl_alpha_C', [3]),
            'chain_C':        ReactiveGroup('chain_C', [1]),
            'vinyl_beta_C':   ReactiveGroup('vinyl_beta_C', [4]),
        }
        positions = np.array([
            [0.0, 0.0, 0.0],   # 0: radical_C
            [1.5, 0.0, 0.0],   # 1: chain_C (bonded to radical_C)
            [5.0, 5.0, 5.0],   # 2: filler
            [4.0, 0.0, 0.0],   # 3: vinyl_alpha_C (4 Å from radical_C)
            [5.3, 0.0, 0.0],   # 4: vinyl_beta_C (1.3 Å from alpha)
        ])
        topo = BondTopology.from_bonds([(0, 1, 1.0), (3, 4, 2.0)])

        cands_no_topo = find_candidates(template, groups, positions)
        assert len(cands_no_topo) >= 1

        cands_with_topo = find_candidates(template, groups, positions, topology=topo)
        assert len(cands_with_topo) >= 1, (
            'Topology guard on constraint_only pair blocked all candidates'
        )

    def test_actual_formation_pair_bonded_is_excluded(self):
        """The main formation pair (radical_C-alpha_C) should be excluded when
        those atoms are already bonded (L6 original intent)."""
        from kagome.reactive.topology import BondTopology

        template = ReactionTemplate(
            name='simple_formation',
            groups=['A', 'B'],
            pairs=[PairSpec('A', 'B', is_formation=True, r_min=0.5, r_max=5.0)],
        )
        groups = {
            'A': ReactiveGroup('A', [0]),
            'B': ReactiveGroup('B', [1]),
        }
        positions = np.array([
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
        ])
        topo = BondTopology.from_bonds([(0, 1, 1.0)])

        cands = find_candidates(template, groups, positions, topology=topo)
        assert len(cands) == 0, (
            'Already-bonded main formation pair should be excluded'
        )


class TestDissociationTopologyGuard:
    """B2 (dual of L6): a dissociation pair may only break a bond that actually
    exists. Candidates pairing a dissociation atom (nylon amine_N / carboxyl_C)
    with a non-bonded partner from a different molecule must be rejected when a
    topology is supplied. Without topology (activation path) the guard is a
    no-op."""

    @staticmethod
    def _nylon_template():
        return ReactionTemplate(
            name='nylon_b2',
            groups=['amine_N', 'carboxyl_C', 'amine_H', 'carboxyl_OH'],
            pairs=[
                PairSpec('amine_N', 'carboxyl_C', is_formation=True,
                         r_min=3.0, r_max=6.0),
                PairSpec('amine_N', 'amine_H', is_formation=False,
                         r_min=0.0, r_max=3.0),
                PairSpec('carboxyl_C', 'carboxyl_OH', is_formation=False,
                         r_min=0.0, r_max=3.0),
                PairSpec('amine_H', 'carboxyl_OH', is_formation=True,
                         r_min=0.0, r_max=100.0, score_pair=False),
            ],
        )

    def test_bonded_dissociation_allowed_nonbonded_rejected(self):
        from kagome.reactive.topology import BondTopology

        template = self._nylon_template()
        # amine_H has two candidates: atom 2 (bonded to N) and atom 4 (a
        # different molecule's H, geometrically close to N but NOT bonded).
        groups = {
            'amine_N':     ReactiveGroup('amine_N', [0]),
            'carboxyl_C':  ReactiveGroup('carboxyl_C', [1]),
            'amine_H':     ReactiveGroup('amine_H', [2, 4]),
            'carboxyl_OH': ReactiveGroup('carboxyl_OH', [3]),
        }
        positions = np.array([
            [0.0, 0.0, 0.0],   # 0 amine_N
            [4.0, 0.0, 0.0],   # 1 carboxyl_C   (N-C = 4.0 in [3,6])
            [0.0, 1.0, 0.0],   # 2 amine_H bonded to N (N-H = 1.0 in [0,3])
            [4.0, 1.0, 0.0],   # 3 carboxyl_OH bonded to C (C-OH = 1.0 in [0,3])
            [0.0, 0.5, 0.0],   # 4 amine_H of another molecule (N-H4 = 0.5, unbonded)
        ])

        # Without topology BOTH amine_H atoms yield a candidate.
        cands_no_topo = find_candidates(template, groups, positions)
        assert {c.atom_indices[2] for c in cands_no_topo} == {2, 4}

        # With the real intramolecular bonds only the bonded amine_H (2) survives.
        topo = BondTopology.from_bonds([(0, 2, 1.0), (1, 3, 1.0)])
        cands = find_candidates(template, groups, positions, topology=topo)
        assert len(cands) == 1
        assert cands[0].atom_indices == (0, 1, 2, 3)

    def test_no_topology_is_noop(self):
        """Activation path (topology=None): B2 does not filter anything."""
        template = self._nylon_template()
        groups = {
            'amine_N':     ReactiveGroup('amine_N', [0]),
            'carboxyl_C':  ReactiveGroup('carboxyl_C', [1]),
            'amine_H':     ReactiveGroup('amine_H', [2, 4]),
            'carboxyl_OH': ReactiveGroup('carboxyl_OH', [3]),
        }
        positions = np.array([
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [4.0, 1.0, 0.0],
            [0.0, 0.5, 0.0],
        ])
        cands = find_candidates(template, groups, positions)
        assert len(cands) == 2
