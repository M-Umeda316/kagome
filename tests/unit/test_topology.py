"""Unit tests for explicit bond-topology tracking (specs/decisions.md 2026-07-02).

Layer 1 of the chemistry-consistency fix: reacted carbons must not appear
over-valent, and no bond may appear at a non-reactive site.
"""
from __future__ import annotations

import numpy as np
import pytest

from kagome.reactive.bonds import BondEvent
from kagome.reactive.topology import (
    BondTopology,
    apply_vinyl_addition,
    over_coordinated_atoms,
    vinyl_addition_over_coordinates,
)

_ELEM_MAXCOORD = {'H': 1, 'C': 4, 'N': 3, 'O': 2, 'F': 1, 'S': 2, 'Cl': 1}


def test_from_bonds_and_accessors():
    topo = BondTopology.from_bonds([(0, 1), (1, 2, 2.0)])
    assert topo.has_bond(1, 0)  # order-independent key
    assert topo.order(1, 2) == 2.0
    assert topo.order(0, 5) == 0.0
    assert sorted(topo.neighbors(1)) == [0, 2]
    assert topo.coordination_number(1) == 2
    assert topo.bonds() == [(0, 1, 1.0), (1, 2, 2.0)]


def test_set_order_requires_existing_bond():
    topo = BondTopology.from_bonds([(0, 1)])
    topo.set_order(0, 1, 2.0)
    assert topo.order(0, 1) == 2.0
    try:
        topo.set_order(0, 9, 1.0)
    except KeyError:
        pass
    else:
        raise AssertionError('set_order on missing bond should raise')


def test_vinyl_addition_opens_double_bond_propagating_radical():
    """Propagating (3-coordinate) radical: no H shed, C=C opens, new sigma added."""
    # beta radical (atom 3, 3-coordinate) adds to a monomer whose alpha=beta
    # double bond is 10=11.
    topo = BondTopology.from_bonds([
        (3, 4), (3, 5), (3, 6),      # radical carbon: 3 bonds (propagating)
        (10, 11, 2.0),               # monomer vinyl C=C (alpha=beta)
        (10, 12), (10, 13),          # alpha extra H/substituent
    ])
    apply_vinyl_addition(topo, radical_c=3, alpha_c=10,
                         propagation_map={10: 11}, species=['C'] * 20)
    assert topo.has_bond(3, 10)          # new sigma bond
    assert topo.order(3, 10) == 1.0
    assert topo.order(10, 11) == 1.0     # double bond opened
    assert topo.coordination_number(3) == 4   # 3 + new bond, no over-valence
    assert topo.coordination_number(10) == 4


def test_vinyl_addition_sheds_spare_h_on_closed_shell_initiator():
    """Closed-shell initiator radical C (4-coordinate incl. placeholder H) sheds
    one H so the new C-C bond keeps it 4-coordinate, not 5."""
    species = ['C', 'C', 'C', 'C', 'H'] + ['C'] * 15  # atom 4 is the spare H
    # radical C = atom 0, bonded to 3 C (1,2,3) + 1 H (4): closed-shell, 4-coord.
    topo = BondTopology.from_bonds([
        (0, 1), (0, 2), (0, 3), (0, 4),
        (10, 11, 2.0),
    ])
    assert topo.coordination_number(0) == 4
    apply_vinyl_addition(topo, radical_c=0, alpha_c=10,
                         propagation_map={10: 11}, species=species)
    assert not topo.has_bond(0, 4)       # spare H shed
    assert topo.has_bond(0, 10)          # new sigma bond
    assert topo.coordination_number(0) == 4   # still 4, not 5


def test_no_bond_ever_created_at_non_reactive_pair():
    """Only radical_C-alpha_C (and the pre-existing intramolecular bonds) exist;
    the edit never touches unrelated atoms."""
    topo = BondTopology.from_bonds([(0, 1), (0, 2), (0, 3), (10, 11, 2.0),
                                    (10, 12), (10, 13)])
    before = set(k for k in (b[:2] for b in topo.bonds()))
    apply_vinyl_addition(topo, radical_c=0, alpha_c=10,
                         propagation_map={10: 11}, species=['C'] * 20)
    after = set(k for k in (b[:2] for b in topo.bonds()))
    # exactly one new pair (0,10); nothing else added
    assert after - before == {(0, 10)}


def test_two_step_chain_propagation_stays_valence_correct():
    """Radical migrates to beta then adds again; both carbons stay <=4-coordinate."""
    species = ['C'] * 40
    # initiator radical 0 (3C + 1H placeholder -> use 3 bonds to keep it simple:
    # here model as 3-coordinate propagating-style to check migration valence)
    topo = BondTopology.from_bonds([
        (0, 1), (0, 2), (0, 3),
        (10, 11, 2.0), (10, 14), (10, 15),   # monomer 1
        (20, 21, 2.0), (20, 24), (20, 25),   # monomer 2
    ])
    pmap = {10: 11, 20: 21}
    # addition 1: radical 0 + alpha 10 -> beta 11 becomes radical
    apply_vinyl_addition(topo, 0, 10, pmap, species)
    # addition 2: radical 11 + alpha 20 -> beta 21 becomes radical
    apply_vinyl_addition(topo, 11, 20, pmap, species)
    assert topo.order(10, 11) == 1.0
    assert topo.order(20, 21) == 1.0
    assert topo.has_bond(0, 10) and topo.has_bond(11, 20)
    for atom in (0, 10, 11, 20, 21):
        assert topo.coordination_number(atom) <= 4


class TestVinylInitialBondsConsistency:
    """vinyl_initial_bonds must line up 1:1 with build_vinyl_aibn_system's
    groups/positions, and reacting every alpha keeps the whole system
    valence-correct (specs/decisions.md 2026-07-02)."""

    @pytest.fixture(autouse=True)
    def _skip_no_rdkit(self):
        pytest.importorskip('rdkit')

    def _build(self, n_monomers=3, n_initiators=2):
        from scripts._systems import build_vinyl_aibn_system, vinyl_initial_bonds
        rng = np.random.default_rng(0)
        pos, sp, tmpl, groups, pmap, ccmap = build_vinyl_aibn_system(
            n_monomers=n_monomers, n_initiators=n_initiators,
            box_size=30.0, rng=rng,
        )
        bonds = vinyl_initial_bonds(n_monomers=n_monomers, n_initiators=n_initiators)
        return sp, groups, pmap, BondTopology.from_bonds(bonds)

    def test_radical_c_is_closed_shell_four_coordinate(self):
        sp, groups, pmap, topo = self._build()
        for rc in groups['radical_C'].atom_indices:
            nbr_species = sorted(sp[n] for n in topo.neighbors(rc))
            assert nbr_species == ['C', 'C', 'C', 'H']  # 3 C + placeholder H

    def test_every_alpha_beta_is_a_double_bond(self):
        sp, groups, pmap, topo = self._build()
        for a in groups['vinyl_alpha_C'].atom_indices:
            assert topo.order(a, pmap[a]) == 2.0

    def test_initial_topology_has_no_over_valence(self):
        sp, groups, pmap, topo = self._build()
        over = [i for i in range(len(sp))
                if topo.coordination_number(i) > _ELEM_MAXCOORD.get(sp[i], 4)]
        assert over == []

    def test_reacting_all_monomers_stays_valence_correct(self):
        sp, groups, pmap, topo = self._build()
        rc = groups['radical_C'].atom_indices[0]
        for a in list(groups['vinyl_alpha_C'].atom_indices):
            apply_vinyl_addition(topo, rc, a, pmap, sp)
            rc = pmap[a]  # radical migrates to beta
        over = [i for i in range(len(sp))
                if topo.coordination_number(i) > _ELEM_MAXCOORD.get(sp[i], 4)]
        assert over == []


class TestValenceGuard:
    """Layer-2 occupancy guard (specs/decisions.md 2026-07-02)."""

    def test_valid_addition_is_not_flagged(self):
        topo = BondTopology.from_bonds([
            (0, 1), (0, 2), (0, 3),          # 3-coord radical
            (10, 11, 2.0), (10, 12), (10, 13),  # alpha =CH2 (3-coord)
        ])
        assert vinyl_addition_over_coordinates(
            topo, 0, 10, {10: 11}, ['C'] * 20) == []

    def test_saturated_alpha_is_flagged(self):
        # alpha 10 already has 4 single bonds (already reacted / saturated).
        topo = BondTopology.from_bonds([
            (0, 1), (0, 2), (0, 3),
            (10, 11), (10, 12), (10, 13), (10, 14),
        ])
        assert vinyl_addition_over_coordinates(
            topo, 0, 10, {10: 11}, ['C'] * 20) == [10]

    def test_closed_shell_radical_not_flagged_due_to_h_shed(self):
        # radical 0 is 4-coord incl. a spare H (atom 4) -> shed keeps it valid.
        species = ['C', 'C', 'C', 'C', 'H'] + ['C'] * 15
        topo = BondTopology.from_bonds([
            (0, 1), (0, 2), (0, 3), (0, 4),
            (10, 11, 2.0), (10, 12), (10, 13),
        ])
        assert vinyl_addition_over_coordinates(
            topo, 0, 10, {10: 11}, species) == []

    def test_undissociated_aibn_center_flagged_then_cleared_by_azo_removal(self):
        """Minimal reproduction of the incomplete-activation bug
        (specs/decisions.md 2026-07-06).

        An AIBN radical centre whose azo C-N bond did NOT dissociate is
        4-coordinate (3 C + 1 azo-N, no H).  A vinyl addition would push it to
        5-coordinate with no spare H to shed, so the valence guard must flag it.
        After the azo C-N bond is removed (genuine dissociation), the same centre
        is 3-coordinate and the addition is valence-safe.
        """
        species = ['C', 'C', 'C', 'C', 'N'] + ['C'] * 15  # atom 0 = radical C, 4 = azo N
        # radical C (0): 3 C (1,2,3) + 1 azo N (4) -> 4-coordinate, no H.
        intact = BondTopology.from_bonds([
            (0, 1), (0, 2), (0, 3), (0, 4),      # intact AIBN centre (with azo N)
            (10, 11, 2.0), (10, 12), (10, 13),   # monomer alpha =CH2
        ])
        # Un-dissociated: adding a monomer over-coordinates the carbon (-> [0]).
        assert vinyl_addition_over_coordinates(
            intact, 0, 10, {10: 11}, species) == [0]

        # After genuine azo C-N dissociation the centre is a real 3-coord radical.
        dissociated = intact.copy()
        dissociated.remove_bond(0, 4)
        assert vinyl_addition_over_coordinates(
            dissociated, 0, 10, {10: 11}, species) == []

    def test_over_coordinated_atoms_scan(self):
        topo = BondTopology.from_bonds([
            (0, 1), (0, 2), (0, 3), (0, 4), (0, 5),  # atom 0 = 5-coord C
        ])
        assert over_coordinated_atoms(topo, ['C'] * 10) == [0]

    def test_workflow_valence_filter_drops_saturated_candidate(self):
        from kagome.backends.toy import ToyCalculator
        from kagome.reactive.bonds import BondTracker
        from kagome.reactive.groups import (
            PairSpec, ReactionTemplate, ReactiveGroup,
        )
        from kagome.reactive.selection import Candidate
        from kagome.workflows.polymerization import (
            PolymerizationConfig, PolymerizationWorkflow, SimulationState,
        )

        template = ReactionTemplate(
            'radical_vinyl', ['radical_C', 'vinyl_alpha_C'],
            [PairSpec('radical_C', 'vinyl_alpha_C', is_formation=True,
                      r_min=1.5, r_max=6.0)],
        )
        groups = {'radical_C': ReactiveGroup('radical_C', [0]),
                  'vinyl_alpha_C': ReactiveGroup('vinyl_alpha_C', [10, 20])}
        species = ['C'] * 30
        # alpha 10 is saturated (4 single bonds); alpha 20 is a normal =CH2.
        initial_bonds = [(0, 1), (0, 2), (0, 3),
                         (10, 11), (10, 12), (10, 13), (10, 14),
                         (20, 21, 2.0), (20, 22), (20, 23)]
        wf = PolymerizationWorkflow(
            PolymerizationConfig(), ToyCalculator(), template, groups,
            bond_tracker=BondTracker(), propagation_map={10: 11, 20: 21},
            initial_bonds=initial_bonds,
        )
        state = SimulationState(positions=np.zeros((30, 3)),
                                velocities=np.zeros((30, 3)), species=species)
        cands = [Candidate(atom_indices=(0, 10)), Candidate(atom_indices=(0, 20))]
        kept = wf._valence_filter(cands, state, cycle=0)
        kept_pairs = [c.atom_indices for c in kept]
        assert (0, 20) in kept_pairs      # valid addition kept
        assert (0, 10) not in kept_pairs  # saturated alpha dropped


class TestWorkflowTopologyOutput:
    """End-to-end: the workflow emits header bonds + topology.jsonl, and a
    confirmed formation produces a valence-correct connectivity snapshot
    (specs/decisions.md 2026-07-02, Layer 1 minimal reproduction)."""

    def _make(self, tmp_path):
        import json

        from kagome.backends.toy import ToyCalculator
        from kagome.reactive.bonds import BondEvent, BondTracker
        from kagome.reactive.groups import (
            PairSpec, ReactionTemplate, ReactiveGroup,
        )
        from kagome.workflows.polymerization import (
            PolymerizationConfig, PolymerizationWorkflow, SimulationState,
            masses_from_species,
        )

        template = ReactionTemplate(
            'radical_vinyl', ['radical_C', 'vinyl_alpha_C'],
            [PairSpec('radical_C', 'vinyl_alpha_C', is_formation=True,
                      r_min=1.5, r_max=6.0)],
        )
        groups = {
            'radical_C': ReactiveGroup('radical_C', [0]),
            'vinyl_alpha_C': ReactiveGroup('vinyl_alpha_C', [10, 20]),
        }
        propagation_map = {10: 11, 20: 21}
        # atom 0 = closed-shell radical C (3 C + 1 H placeholder at atom 4)
        species = ['C', 'C', 'C', 'C', 'H'] + ['C'] * 26
        initial_bonds = [(0, 1, 1.0), (0, 2, 1.0), (0, 3, 1.0), (0, 4, 1.0),
                         (10, 11, 2.0), (20, 21, 2.0)]
        tracker = BondTracker()
        cfg = PolymerizationConfig(
            n_cycles=1, biased_steps=1, unbiased_steps=1, save_interval=1,
        )
        wf = PolymerizationWorkflow(
            cfg, ToyCalculator(), template, groups, bond_tracker=tracker,
            propagation_map=propagation_map, initial_bonds=initial_bonds,
        )
        n = len(species)
        rng = np.random.default_rng(0)
        state = SimulationState(
            positions=rng.uniform(0, 20, size=(n, 3)),
            velocities=np.zeros((n, 3)), species=species,
            masses=masses_from_species(species),
        )
        return wf, tracker, state, json

    def test_header_carries_initial_bonds(self, tmp_path):
        wf, tracker, state, json = self._make(tmp_path)
        wf.run(state, output_dir=tmp_path)
        header = json.loads((tmp_path / 'trajectory.jsonl').read_text(
            encoding='utf-8').splitlines()[0])
        assert 'bonds' in header
        assert [0, 1, 1.0] in header['bonds']

    def test_topology_jsonl_initial_snapshot_written(self, tmp_path):
        wf, tracker, state, json = self._make(tmp_path)
        wf.run(state, output_dir=tmp_path)
        lines = (tmp_path / 'topology.jsonl').read_text(
            encoding='utf-8').splitlines()
        first = json.loads(lines[0])
        assert first['cycle'] == -1
        assert first['n_bonds'] == 6

    def test_confirmed_formation_emits_valence_correct_snapshot(self, tmp_path):
        wf, tracker, state, json = self._make(tmp_path)
        wf.run(state, output_dir=tmp_path)  # sets up topology_log + initial snapshot
        # Inject a confirmed radical(0) + alpha(10) formation, then run the
        # post-cycle update as the loop would.
        tracker._events.append(BondEvent(
            step=state.step, cycle=0, atom_a=0, atom_b=10,
            event_type='confirmed_formation', distance=1.6, r0=2.04,
        ))
        wf._update_groups_after_cycle(state, cycle=0)
        last = json.loads((tmp_path / 'topology.jsonl').read_text(
            encoding='utf-8').splitlines()[-1])
        pairs = {(i, j) for i, j, _ in last['bonds']}
        assert (0, 10) in pairs           # new sigma bond recorded
        assert (0, 4) not in pairs        # closed-shell placeholder H shed
        # C=C opened
        order_10_11 = next(o for i, j, o in last['bonds'] if (i, j) == (10, 11))
        assert order_10_11 == 1.0
        # no carbon over-valent
        topo = BondTopology.from_bonds([(b[0], b[1], b[2]) for b in last['bonds']])
        over = [a for a in range(len(state.species))
                if topo.coordination_number(a)
                > _ELEM_MAXCOORD.get(state.species[a], 4)]
        assert over == []
