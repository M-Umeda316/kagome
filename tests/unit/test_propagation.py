"""Unit tests for chain-propagation bookkeeping (handoff-plan-v4 S1, T-S1.1).

After a confirmed formation the workflow must: remove the consumed radical (atom_a)
and the reacted vinyl alpha-C (atom_b) from all reactive groups, and add the
monomer's beta-C (propagation_map[alpha]) to the radical group — so the active
radical migrates to the new chain end. Along a single chain there is always
exactly one radical end (the doublet/spin invariant for S1).
"""
from __future__ import annotations

import numpy as np

from src.backends.toy import ToyCalculator
from src.reactive.bonds import BondEvent, BondTracker
from src.reactive.groups import PairSpec, ReactionTemplate, ReactiveGroup
from src.workflows.polymerization import (
    PolymerizationConfig,
    PolymerizationWorkflow,
    SimulationState,
)


def _make_wf():
    template = ReactionTemplate(
        name='radical_vinyl',
        groups=['radical_C', 'vinyl_alpha_C'],
        pairs=[PairSpec('radical_C', 'vinyl_alpha_C', is_formation=True,
                        r_min=1.5, r_max=6.0)],
    )
    groups = {
        'radical_C': ReactiveGroup('radical_C', [0]),
        'vinyl_alpha_C': ReactiveGroup('vinyl_alpha_C', [10, 20]),
    }
    propagation_map = {10: 11, 20: 21}  # each monomer's alpha-C -> beta-C
    tracker = BondTracker()
    wf = PolymerizationWorkflow(
        PolymerizationConfig(), ToyCalculator(), template, groups,
        bond_tracker=tracker, propagation_map=propagation_map,
        propagation_target_group='radical_C',
    )
    return wf, groups, tracker


def _state():
    n = 30
    return SimulationState(
        positions=np.zeros((n, 3)), velocities=np.zeros((n, 3)),
        species=['C'] * n,
    )


def _confirm(tracker, atom_a, atom_b, step, cycle):
    tracker._events.append(BondEvent(
        step=step, cycle=cycle, atom_a=atom_a, atom_b=atom_b,
        event_type='confirmed_formation', distance=1.6, r0=2.04,
    ))


def test_radical_migrates_across_two_additions():
    wf, groups, tracker = _make_wf()
    state = _state()

    # Addition 1: radical(0) + monomer-1 alpha(10) -> new radical at beta(11)
    _confirm(tracker, 0, 10, step=1, cycle=0)
    wf._update_groups_after_cycle(state)
    assert 0 not in groups['radical_C'].atom_indices
    assert groups['radical_C'].atom_indices == [11]
    assert 10 not in groups['vinyl_alpha_C'].atom_indices
    assert 20 in groups['vinyl_alpha_C'].atom_indices

    # Addition 2: new radical(11) + monomer-2 alpha(20) -> radical migrates to beta(21)
    _confirm(tracker, 11, 20, step=2, cycle=1)
    wf._update_groups_after_cycle(state)
    assert groups['radical_C'].atom_indices == [21]
    assert groups['vinyl_alpha_C'].atom_indices == []


def test_single_radical_invariant_holds_each_step():
    """Exactly one radical end at all times along a single chain (S1 doublet)."""
    wf, groups, tracker = _make_wf()
    state = _state()
    assert len(groups['radical_C']) == 1

    _confirm(tracker, 0, 10, step=1, cycle=0)
    wf._update_groups_after_cycle(state)
    assert len(groups['radical_C']) == 1  # radical migrated, not multiplied

    _confirm(tracker, 11, 20, step=2, cycle=1)
    wf._update_groups_after_cycle(state)
    assert len(groups['radical_C']) == 1


def test_processed_formations_advances_no_double_apply():
    wf, groups, tracker = _make_wf()
    state = _state()
    _confirm(tracker, 0, 10, step=1, cycle=0)
    wf._update_groups_after_cycle(state)
    # calling again with no new events must not change groups
    before = list(groups['radical_C'].atom_indices)
    wf._update_groups_after_cycle(state)
    assert groups['radical_C'].atom_indices == before
