"""Tests for reactive group definitions and ReactionTemplate validation (Eq. 6).

ReactionTemplate.__post_init__ enforces three invariants (L4 re-check):
duplicate group labels, pair references to unknown groups, and symmetric
pairs (group_a == group_b) — see src/kagome/reactive/groups.py.
"""
import pytest

from kagome.reactive.groups import PairSpec, ReactiveGroup, ReactionTemplate


class TestReactiveGroup:

    def test_len_reflects_atom_indices(self):
        group = ReactiveGroup(label='I', atom_indices=[0, 1, 2])
        assert len(group) == 3

    def test_remove_atom_present(self):
        group = ReactiveGroup(label='I', atom_indices=[0, 1, 2])
        group.remove_atom(1)
        assert group.atom_indices == [0, 2]

    def test_remove_atom_absent_is_noop(self):
        group = ReactiveGroup(label='I', atom_indices=[0, 1])
        group.remove_atom(99)
        assert group.atom_indices == [0, 1]


class TestReactionTemplateValid:

    def test_valid_template_constructs(self):
        template = ReactionTemplate(
            name='test_reaction',
            groups=['I', 'J'],
            pairs=[PairSpec(group_a='I', group_b='J', is_formation=True,
                             r_min=0.5, r_max=3.0)],
        )
        assert template.group_labels() == {'I', 'J'}
        assert template.pair_labels() == [('I', 'J')]


class TestReactionTemplateDuplicateGroupLabels:

    def test_duplicate_group_labels_raise(self):
        with pytest.raises(ValueError, match='duplicate group labels'):
            ReactionTemplate(
                name='dup_reaction',
                groups=['I', 'J', 'I'],
                pairs=[PairSpec(group_a='I', group_b='J', is_formation=True)],
            )


class TestReactionTemplateUnknownGroupReference:

    def test_unknown_group_a_raises(self):
        with pytest.raises(ValueError, match="group_a 'K' not in groups"):
            ReactionTemplate(
                name='unknown_a',
                groups=['I', 'J'],
                pairs=[PairSpec(group_a='K', group_b='J', is_formation=True)],
            )

    def test_unknown_group_b_raises(self):
        with pytest.raises(ValueError, match="group_b 'K' not in groups"):
            ReactionTemplate(
                name='unknown_b',
                groups=['I', 'J'],
                pairs=[PairSpec(group_a='I', group_b='K', is_formation=True)],
            )


class TestReactionTemplateSymmetricPairRejected:

    def test_group_a_equals_group_b_raises(self):
        # L4 re-check: a pair (I, I) maps to a (i, i) enumeration key that the
        # candidate enumerator's prev_depth < depth loop can never match, so
        # the distance window would be silently disabled instead of erroring.
        with pytest.raises(ValueError, match='symmetric reactions are not supported'):
            ReactionTemplate(
                name='symmetric_reaction',
                groups=['I', 'J'],
                pairs=[PairSpec(group_a='I', group_b='I', is_formation=True)],
            )
