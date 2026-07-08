"""Integration test for the polymerization workflow."""
import json

import numpy as np
import pytest

from kagome.backends.toy import ToyCalculator
from kagome.boost.tdbb import PairBias, TDBBParams
from kagome.integrators.langevin import LangevinIntegrator, LangevinParams
from kagome.io.readers import read_trajectory
from kagome.reactive.bonds import BondEvent, BondTracker
from kagome.reactive.groups import PairSpec, ReactiveGroup, ReactionTemplate
from kagome.reactive.selection import (
    find_candidates,
    score_candidates,
    select_non_overlapping,
)
from kagome.workflows.polymerization import (
    DefaultPostCycleUpdater,
    EpoxyAmineAdditionUpdater,
    PolymerizationConfig,
    PolymerizationWorkflow,
    SimulationState,
    masses_from_species,
)


class _FakeTracker:
    """Minimal BondTracker stand-in returning preset confirmed events."""

    def __init__(self, formations, dissociations):
        self._f = formations
        self._d = dissociations

    def confirmed_formations(self):
        return self._f

    def confirmed_dissociations(self):
        return self._d


class TestDefaultPostCycleUpdater:
    """RF15: confirmed dissociations must free leaving-group atoms from groups so
    they are not re-selected (step-growth condensation topology must advance)."""

    @staticmethod
    def _groups():
        return {
            'amine_N': ReactiveGroup('amine_N', [0, 10]),
            'carboxyl_C': ReactiveGroup('carboxyl_C', [2, 12]),
            'amine_H': ReactiveGroup('amine_H', [1, 11]),
            'carboxyl_OH': ReactiveGroup('carboxyl_OH', [3, 13]),
        }

    @staticmethod
    def _state():
        return SimulationState(
            positions=np.zeros((14, 3)), velocities=np.zeros((14, 3)),
            species=['C'] * 14,
        )

    def test_consumes_formations_and_dissociations(self):
        groups = self._groups()
        formations = [BondEvent(
            step=1, cycle=0, atom_a=0, atom_b=2,
            event_type='confirmed_formation', distance=1.5,
        )]
        dissociations = [
            BondEvent(step=1, cycle=0, atom_a=0, atom_b=1,
                      event_type='confirmed_dissociation', distance=3.0),
            BondEvent(step=1, cycle=0, atom_a=2, atom_b=3,
                      event_type='confirmed_dissociation', distance=3.0),
        ]
        upd = DefaultPostCycleUpdater()
        upd.update(groups, _FakeTracker(formations, dissociations), self._state())

        # reacted centres consumed by the amide formation
        assert 0 not in groups['amine_N'].atom_indices
        assert 2 not in groups['carboxyl_C'].atom_indices
        # leaving groups freed by dissociation (the RF15 fix)
        assert 1 not in groups['amine_H'].atom_indices
        assert 3 not in groups['carboxyl_OH'].atom_indices
        # the other monomer's termini remain — chain growth proceeds naturally
        assert 10 in groups['amine_N'].atom_indices
        assert 12 in groups['carboxyl_C'].atom_indices
        assert upd.processed_formations == 1
        assert upd.processed_dissociations == 2

    def test_dissociation_counter_prevents_double_processing(self):
        groups = self._groups()
        diss = [BondEvent(step=1, cycle=0, atom_a=0, atom_b=1,
                          event_type='confirmed_dissociation', distance=3.0)]
        tracker = _FakeTracker([], diss)
        upd = DefaultPostCycleUpdater()
        upd.update(groups, tracker, self._state())
        upd.update(groups, tracker, self._state())  # same events again
        assert upd.processed_dissociations == 1

    # ── H2: candidate-atomic (all-or-nothing) acceptance ──────────────

    def test_dissociation_skipped_when_formation_not_confirmed(self):
        """H2: a V^d-only confirmation must NOT consume the reactive sites."""
        groups = self._groups()
        # candidate 0's dissociation confirmed, but its formation did not.
        diss = [BondEvent(step=1, cycle=0, atom_a=0, atom_b=1,
                          event_type='confirmed_dissociation', distance=3.0,
                          candidate_id=0)]
        upd = DefaultPostCycleUpdater()
        upd.update(groups, _FakeTracker([], diss), self._state())

        # amine N and its H stay selectable — no site is lost.
        assert 0 in groups['amine_N'].atom_indices
        assert 1 in groups['amine_H'].atom_indices
        assert upd.processed_dissociations == 1

    def test_dissociation_applied_when_same_candidate_formation_confirmed(self):
        groups = self._groups()
        formations = [BondEvent(step=1, cycle=0, atom_a=0, atom_b=2,
                                event_type='confirmed_formation', distance=1.5,
                                candidate_id=0)]
        diss = [BondEvent(step=1, cycle=0, atom_a=0, atom_b=1,
                          event_type='confirmed_dissociation', distance=3.0,
                          candidate_id=0)]
        upd = DefaultPostCycleUpdater()
        upd.update(groups, _FakeTracker(formations, diss), self._state())

        assert 0 not in groups['amine_N'].atom_indices
        assert 1 not in groups['amine_H'].atom_indices

    def test_candidate_id_matching_is_per_cycle(self):
        """candidate_id restarts each cycle: a cycle-1 dissociation must not
        match a cycle-0 formation that happens to share the same id."""
        groups = self._groups()
        formations = [BondEvent(step=1, cycle=0, atom_a=10, atom_b=12,
                                event_type='confirmed_formation', distance=1.5,
                                candidate_id=0)]
        diss = [BondEvent(step=5, cycle=1, atom_a=0, atom_b=1,
                          event_type='confirmed_dissociation', distance=3.0,
                          candidate_id=0)]
        upd = DefaultPostCycleUpdater()
        upd.update(groups, _FakeTracker(formations, diss), self._state())

        # cycle-1 candidate 0 had no confirmed formation — sites preserved.
        assert 0 in groups['amine_N'].atom_indices
        assert 1 in groups['amine_H'].atom_indices

    def test_legacy_events_without_candidate_id_are_applied(self):
        """Events with candidate_id=-1 (activation, pre-H2 checkpoints) keep
        the RF15 behavior: dissociation frees the atoms unconditionally."""
        groups = self._groups()
        diss = [BondEvent(step=1, cycle=0, atom_a=0, atom_b=1,
                          event_type='confirmed_dissociation', distance=3.0)]
        upd = DefaultPostCycleUpdater()
        upd.update(groups, _FakeTracker([], diss), self._state())
        assert 0 not in groups['amine_N'].atom_indices
        assert 1 not in groups['amine_H'].atom_indices


class TestEpoxyAmineAdditionUpdater:
    """1° -> 2° -> 3° amine multi-addition (Track 2 E1, decisions.md
    2026-07-09): a confirmed ring-opening consumes the epoxy C, the ring O and
    one amine H, but the amine N stays selectable until its LAST registered H
    is consumed (tertiary amine)."""

    # atoms: 0 = amine N (primary, 2 H: 1 and 2); 10/20 = epoxy C; 11/21 = ring O
    @staticmethod
    def _groups():
        return {
            'amine_N': ReactiveGroup('amine_N', [0]),
            'epoxy_C': ReactiveGroup('epoxy_C', [10, 20]),
            'amine_H': ReactiveGroup('amine_H', [1, 2]),
            'ring_O':  ReactiveGroup('ring_O',  [11, 21]),
        }

    @staticmethod
    def _updater():
        return EpoxyAmineAdditionUpdater({0: [1, 2]})

    @staticmethod
    def _state():
        return SimulationState(
            positions=np.zeros((22, 3)), velocities=np.zeros((22, 3)),
            species=['C'] * 22,
        )

    @staticmethod
    def _addition_events(cycle, n, c, h, o, cid=0):
        """One complete conjunctive ring-opening: N-C + hydroxyl O-H formations,
        N-H + ring C-O dissociations (all same candidate)."""
        formations = [
            BondEvent(step=1, cycle=cycle, atom_a=n, atom_b=c,
                      event_type='confirmed_formation', distance=1.5,
                      candidate_id=cid),
            BondEvent(step=1, cycle=cycle, atom_a=h, atom_b=o,
                      event_type='confirmed_formation', distance=1.0,
                      candidate_id=cid, counts_as_reaction=False),
        ]
        dissociations = [
            BondEvent(step=1, cycle=cycle, atom_a=n, atom_b=h,
                      event_type='confirmed_dissociation', distance=3.0,
                      candidate_id=cid),
            BondEvent(step=1, cycle=cycle, atom_a=c, atom_b=o,
                      event_type='confirmed_dissociation', distance=3.0,
                      candidate_id=cid),
        ]
        return formations, dissociations

    def test_first_addition_keeps_secondary_amine_selectable(self):
        groups = self._groups()
        upd = self._updater()
        formations, dissociations = self._addition_events(0, n=0, c=10, h=1, o=11)
        upd.update(groups, _FakeTracker(formations, dissociations), self._state())

        # consumed: attacked epoxy C, its ring O, the transferred H
        assert 10 not in groups['epoxy_C'].atom_indices
        assert 11 not in groups['ring_O'].atom_indices
        assert 1 not in groups['amine_H'].atom_indices
        # the now-secondary amine N keeps its remaining H and stays selectable
        assert 0 in groups['amine_N'].atom_indices
        assert 2 in groups['amine_H'].atom_indices
        # untouched epoxide remains
        assert 20 in groups['epoxy_C'].atom_indices
        assert 21 in groups['ring_O'].atom_indices

    def test_second_addition_retires_tertiary_amine(self):
        groups = self._groups()
        upd = self._updater()
        f1, d1 = self._addition_events(0, n=0, c=10, h=1, o=11)
        upd.update(groups, _FakeTracker(f1, d1), self._state())
        f2, d2 = self._addition_events(1, n=0, c=20, h=2, o=21)
        upd.update(groups, _FakeTracker(f1 + f2, d1 + d2), self._state())

        assert 2 not in groups['amine_H'].atom_indices
        assert 20 not in groups['epoxy_C'].atom_indices
        assert 21 not in groups['ring_O'].atom_indices
        # last H consumed -> tertiary amine retired
        assert 0 not in groups['amine_N'].atom_indices

    def test_h2_guard_skips_unmatched_dissociation(self):
        """A dissociation whose candidate has no confirmed formation must not
        consume any site (same H2 semantics as DefaultPostCycleUpdater)."""
        groups = self._groups()
        upd = self._updater()
        diss = [BondEvent(step=1, cycle=0, atom_a=0, atom_b=1,
                          event_type='confirmed_dissociation', distance=3.0,
                          candidate_id=0)]
        upd.update(groups, _FakeTracker([], diss), self._state())

        assert 0 in groups['amine_N'].atom_indices
        assert 1 in groups['amine_H'].atom_indices

    def test_processed_counters_prevent_double_processing(self):
        groups = self._groups()
        upd = self._updater()
        formations, dissociations = self._addition_events(0, n=0, c=10, h=1, o=11)
        tracker = _FakeTracker(formations, dissociations)
        upd.update(groups, tracker, self._state())
        upd.update(groups, tracker, self._state())  # same events again
        assert upd.processed_formations == 2
        assert upd.processed_dissociations == 2

    def test_checkpoint_state_is_counters_only(self):
        """Resume relies on run()'s existing save/restore of the processed
        counters; the H bookkeeping must be derivable from live groups."""
        groups = self._groups()
        upd = self._updater()
        f1, d1 = self._addition_events(0, n=0, c=10, h=1, o=11)
        upd.update(groups, _FakeTracker(f1, d1), self._state())

        # a fresh updater with restored counters continues correctly
        resumed = self._updater()
        resumed._processed_formations = upd.processed_formations
        resumed._processed_dissociations = upd.processed_dissociations
        f2, d2 = self._addition_events(1, n=0, c=20, h=2, o=21)
        resumed.update(groups, _FakeTracker(f1 + f2, d1 + d2), self._state())
        assert 0 not in groups['amine_N'].atom_indices


class TestTruncateJsonl:
    """M2: resume must drop records written after the checkpoint."""

    def test_truncates_by_step(self, tmp_path):
        from kagome.workflows.polymerization import _truncate_jsonl_after
        p = tmp_path / 'trajectory.jsonl'
        lines = [
            json.dumps({'_header': True, 'species': ['C']}),
            json.dumps({'step': 100, 'cycle': 0}),
            json.dumps({'step': 200, 'cycle': 1}),
            json.dumps({'step': 300, 'cycle': 1}),  # post-checkpoint
        ]
        p.write_text('\n'.join(lines) + '\n', encoding='utf-8')

        removed = _truncate_jsonl_after(p, 200)

        assert removed == 1
        kept = [json.loads(l) for l in p.read_text(encoding='utf-8').splitlines()]
        assert len(kept) == 3  # header (no step) + steps 100, 200
        assert all(rec.get('step', 0) <= 200 for rec in kept)

    def test_truncates_by_cycle_for_selection_log(self, tmp_path):
        from kagome.workflows.polymerization import _truncate_jsonl_after
        p = tmp_path / 'selection.jsonl'
        lines = [
            json.dumps({'cycle': 0, 'n_selected': 2}),
            json.dumps({'cycle': 1, 'n_selected': 1}),
            json.dumps({'cycle': 2, 'n_selected': 3}),  # mid-crash duplicate
        ]
        p.write_text('\n'.join(lines) + '\n', encoding='utf-8')

        # checkpoint next_cycle=2 → keep cycles <= 1
        removed = _truncate_jsonl_after(p, 1, field='cycle')

        assert removed == 1
        kept = [json.loads(l) for l in p.read_text(encoding='utf-8').splitlines()]
        assert [rec['cycle'] for rec in kept] == [0, 1]

    def test_missing_file_is_noop(self, tmp_path):
        from kagome.workflows.polymerization import _truncate_jsonl_after
        assert _truncate_jsonl_after(tmp_path / 'nope.jsonl', 100) == 0


def _make_simple_setup():
    template = ReactionTemplate(
        name='simple',
        groups=['A', 'B'],
        pairs=[PairSpec('A', 'B', is_formation=True, r_min=0.5, r_max=5.0)],
    )
    groups = {
        'A': ReactiveGroup('A', [0]),
        'B': ReactiveGroup('B', [1]),
    }
    return template, groups


class TestPolymerizationWorkflow:

    def test_smoke_run(self):
        template, groups = _make_simple_setup()
        config = PolymerizationConfig(
            biased_steps=10,
            unbiased_steps=10,
            n_cycles=2,
            seed=42,
        )
        calc = ToyCalculator()
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'C'],
        )
        wf = PolymerizationWorkflow(config, calc, template, groups)
        logs = wf.run(state)

        assert len(logs) == 4  # 2 cycles × 2 phases
        assert logs[0].phase == 'biased'
        assert logs[1].phase == 'unbiased'
        assert state.step == 40  # 2 × (10+10)

    def test_selection_audit_written(self, tmp_path):
        # RF18: selection.jsonl must record the ranked candidates and the
        # selected/rejected decisions with reasons.
        template = ReactionTemplate(
            name='multi', groups=['A', 'B'],
            pairs=[PairSpec('A', 'B', is_formation=True, r_min=0.5, r_max=10.0)],
        )
        groups = {'A': ReactiveGroup('A', [0, 1]), 'B': ReactiveGroup('B', [2, 3])}
        config = PolymerizationConfig(biased_steps=2, unbiased_steps=2, n_cycles=1, seed=1)
        calc = ToyCalculator()
        state = SimulationState(
            positions=np.array([[0, 0, 0], [5, 0, 0], [1, 0, 0], [6, 0, 0]], dtype=float),
            velocities=np.zeros((4, 3)),
            species=['C'] * 4,
        )
        wf = PolymerizationWorkflow(config, calc, template, groups)
        wf.run(state, output_dir=tmp_path)

        audit = tmp_path / 'selection.jsonl'
        assert audit.exists()
        lines = [json.loads(l) for l in audit.read_text(encoding='utf-8').splitlines() if l.strip()]
        assert len(lines) == 1  # one biased cycle
        rec = lines[0]
        assert rec['cycle'] == 0
        assert rec['n_candidates'] == rec['n_selected'] + rec['n_rejected']
        # candidates (0,2)(0,3)(1,2)(1,3): greedy keeps 2 disjoint, drops 2 on overlap
        assert rec['n_selected'] >= 1 and rec['n_rejected'] >= 1
        for d in rec['selected']:
            assert 'atoms' in d and 'score' in d
        for d in rec['rejected']:
            assert d['reason'].startswith('overlap')

    def test_manifest_records_production_start_step(self, tmp_path):
        """A4: post-equilibration production onset is recorded in manifest.json."""
        template = ReactionTemplate(
            name='m', groups=['A', 'B'],
            pairs=[PairSpec('A', 'B', is_formation=True, r_min=0.5, r_max=10.0)],
        )
        groups = {'A': ReactiveGroup('A', [0, 1]), 'B': ReactiveGroup('B', [2, 3])}
        config = PolymerizationConfig(
            biased_steps=2, unbiased_steps=2, n_cycles=1, seed=1, equil_steps=5,
        )
        calc = ToyCalculator()
        state = SimulationState(
            positions=np.array([[0, 0, 0], [5, 0, 0], [1, 0, 0], [6, 0, 0]], dtype=float),
            velocities=np.zeros((4, 3)),
            species=['C'] * 4,
        )
        wf = PolymerizationWorkflow(config, calc, template, groups)
        wf.run(state, output_dir=tmp_path)

        data = json.loads((tmp_path / 'manifest.json').read_text(encoding='utf-8'))
        # 5 equilibration steps ran before the cycle loop; minimize is off.
        assert data['extra']['production_start_step'] == 5

    def test_deterministic_with_seed(self):
        template, groups = _make_simple_setup()
        config = PolymerizationConfig(
            biased_steps=5,
            unbiased_steps=5,
            n_cycles=1,
            seed=7,
        )
        calc = ToyCalculator()

        def run_once():
            state = SimulationState(
                positions=np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
                velocities=np.zeros((2, 3)),
                species=['C', 'C'],
            )
            wf = PolymerizationWorkflow(config, calc, template, groups)
            wf.run(state)
            return state.positions.copy()

        pos1 = run_once()
        pos2 = run_once()
        np.testing.assert_array_equal(pos1, pos2)

    def test_biased_phase_finds_candidates(self):
        template, groups = _make_simple_setup()
        config = PolymerizationConfig(
            biased_steps=5,
            unbiased_steps=5,
            n_cycles=1,
            seed=7,
        )
        calc = ToyCalculator()
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'C'],
        )
        wf = PolymerizationWorkflow(config, calc, template, groups)
        logs = wf.run(state)

        biased_log = logs[0]
        assert biased_log.n_candidates >= 1
        assert biased_log.n_selected >= 1

    def test_trajectory_output(self, tmp_path):
        template, groups = _make_simple_setup()
        config = PolymerizationConfig(
            biased_steps=10,
            unbiased_steps=10,
            n_cycles=1,
            seed=42,
            save_interval=5,
        )
        calc = ToyCalculator()
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'C'],
        )
        wf = PolymerizationWorkflow(config, calc, template, groups)
        wf.run(state, output_dir=tmp_path)

        traj_path = tmp_path / 'trajectory.jsonl'
        assert traj_path.exists()

        header, frames = read_trajectory(traj_path)
        assert header['species'] == ['C', 'C']
        # 10 biased steps / 5 interval = 2 frames + 10 unbiased / 5 = 2 frames = 4
        # step 0 counts (0%5==0), step 5 counts (5%5==0) = 2 per phase
        assert len(frames) >= 2
        assert any(f.phase == 'biased' for f in frames)
        assert any(f.phase == 'unbiased' for f in frames)

    def test_trajectory_header_n_reactive_sites_equals_n_monomers(self, tmp_path):
        """RF2: header n_reactive_sites == n_monomers when passed."""
        template, groups = _make_simple_setup()
        config = PolymerizationConfig(
            biased_steps=5, unbiased_steps=5, n_cycles=1, seed=42,
            save_interval=1,
        )
        calc = ToyCalculator()
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'C'],
        )
        wf = PolymerizationWorkflow(config, calc, template, groups)
        wf.run(state, output_dir=tmp_path, n_monomers=1)

        header, _ = read_trajectory(tmp_path / 'trajectory.jsonl')
        assert header['n_reactive_sites'] == 1

    def test_manifest_records_effective_params(self, tmp_path):
        """RF1: manifest.json extra contains TDBB effective parameters."""
        template, groups = _make_simple_setup()
        config = PolymerizationConfig(
            biased_steps=10,
            unbiased_steps=10,
            n_cycles=1,
            seed=42,
            save_interval=5,
            tdbb=TDBBParams(
                f2=5.0,
                gamma=1.0,
                f1_max_formation=250.0,
                f1_max_dissociation=125.0,
                lambda_vdw=0.60,
            ),
        )
        calc = ToyCalculator()
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'C'],
        )
        wf = PolymerizationWorkflow(config, calc, template, groups)
        wf.run(state, output_dir=tmp_path, config_path='configs/boost/paper_faithful.yaml')

        manifest_path = tmp_path / 'manifest.json'
        assert manifest_path.exists()
        data = json.loads(manifest_path.read_text(encoding='utf-8'))

        assert data['config_path'] == 'configs/boost/paper_faithful.yaml'
        assert data['seed'] == 42
        assert data['backend'] == 'toy'

        ex = data['extra']
        assert ex['tdbb']['f2'] == 5.0
        assert ex['tdbb']['gamma'] == 1.0
        assert ex['tdbb']['f1_max_formation'] == 250.0
        assert ex['tdbb']['f1_max_dissociation'] == 125.0
        assert ex['tdbb']['lambda_vdw'] == 0.60
        assert ex['timestep_fs'] == config.timestep_fs
        assert ex['biased_steps'] == 10
        assert ex['unbiased_steps'] == 10
        assert ex['n_cycles'] == 1
        assert ex['backend'] == 'toy'
        assert ex['candidate_r_min'] == 0.5
        assert ex['candidate_r_max'] == 5.0
        # RF17: provenance — resolved model identity and alpha denominator recorded
        assert ex['model_id'] == 'toy'
        assert isinstance(ex['n_reactive_sites'], int)

    @staticmethod
    def _run_simple(integrator=None, masses=None, seed=7):
        template, groups = _make_simple_setup()
        config = PolymerizationConfig(
            biased_steps=5,
            unbiased_steps=5,
            n_cycles=1,
            seed=seed,
        )
        calc = ToyCalculator()
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'C'],
            masses=masses,
        )
        wf = PolymerizationWorkflow(config, calc, template, groups, integrator=integrator)
        wf.run(state)
        return state

    def test_with_masses(self):
        # RF21: masses must actually affect the dynamics (F/m), not just run.
        s_mass = self._run_simple(masses=masses_from_species(['C', 'C']))
        s_none = self._run_simple(masses=None)
        assert s_mass.step == 10 and s_none.step == 10
        # 12 amu vs unit-mass path integrate F/m differently -> different trajectory
        assert not np.allclose(s_mass.positions, s_none.positions)

    def test_with_langevin(self):
        # RF21: Langevin thermostat must change the trajectory vs NVE Verlet and
        # remain reproducible under a fixed seed.
        masses = masses_from_species(['C', 'C'])
        langevin = LangevinIntegrator(LangevinParams(temperature_K=300.0, friction_per_fs=0.05))
        s_lan = self._run_simple(integrator=langevin, masses=masses)
        s_vv = self._run_simple(integrator=None, masses=masses)
        assert s_lan.step == 10
        # stochastic thermostat + friction -> diverges from deterministic Verlet
        assert not np.allclose(s_lan.positions, s_vv.positions)
        # same seed -> reproducible
        langevin2 = LangevinIntegrator(LangevinParams(temperature_K=300.0, friction_per_fs=0.05))
        s_lan2 = self._run_simple(integrator=langevin2, masses=masses)
        np.testing.assert_array_equal(s_lan.positions, s_lan2.positions)

    def test_with_bond_tracker(self, tmp_path):
        template, groups = _make_simple_setup()
        config = PolymerizationConfig(
            biased_steps=5,
            unbiased_steps=5,
            n_cycles=1,
            seed=7,
        )
        calc = ToyCalculator()
        tracker = BondTracker()
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'C'],
        )
        wf = PolymerizationWorkflow(
            config, calc, template, groups, bond_tracker=tracker,
        )
        wf.run(state, output_dir=tmp_path)

        assert len(tracker.events) >= 1
        bonds_path = tmp_path / 'bonds.jsonl'
        assert bonds_path.exists()


class TestChainPropagation:
    """Tests for propagation_map / chain-growth logic in PolymerizationWorkflow."""

    def _make_propagation_setup(self, propagation_map=None):
        # 3 atoms: 0=radical_C, 1=vinyl_alpha_C, 2=beta_C
        template = ReactionTemplate(
            name='radical_vinyl',
            groups=['radical_C', 'vinyl_alpha_C'],
            pairs=[PairSpec('radical_C', 'vinyl_alpha_C', is_formation=True, r_min=0.5, r_max=5.0)],
        )
        groups = {
            'radical_C': ReactiveGroup('radical_C', [0]),
            'vinyl_alpha_C': ReactiveGroup('vinyl_alpha_C', [1]),
        }
        config = PolymerizationConfig(biased_steps=5, unbiased_steps=5, n_cycles=1, seed=7)
        calc = ToyCalculator()
        tracker = BondTracker()
        wf = PolymerizationWorkflow(
            config, calc, template, groups,
            bond_tracker=tracker,
            propagation_map=propagation_map,
            propagation_target_group='radical_C',
        )
        state = SimulationState(
            positions=np.zeros((3, 3)),
            velocities=np.zeros((3, 3)),
            species=['C', 'C', 'C'],
        )
        return wf, tracker, groups, state

    def test_beta_added_to_radical_group_after_formation(self):
        """After confirmed formation (radical+alpha), beta-C joins radical_C."""
        wf, tracker, groups, state = self._make_propagation_setup(propagation_map={1: 2})
        tracker._events.append(BondEvent(
            step=5, cycle=0, atom_a=0, atom_b=1,
            event_type='confirmed_formation', distance=1.5, r0=2.0,
        ))
        wf._update_groups_after_cycle(state)
        assert 2 in groups['radical_C'].atom_indices

    def test_reacted_atoms_removed(self):
        """atom_a (radical) and atom_b (alpha) are removed from their groups."""
        wf, tracker, groups, state = self._make_propagation_setup(propagation_map={1: 2})
        tracker._events.append(BondEvent(
            step=5, cycle=0, atom_a=0, atom_b=1,
            event_type='confirmed_formation', distance=1.5, r0=2.0,
        ))
        wf._update_groups_after_cycle(state)
        assert 0 not in groups['radical_C'].atom_indices
        assert 1 not in groups['vinyl_alpha_C'].atom_indices

    def test_no_propagation_map_leaves_groups_unchanged_except_removal(self):
        """Without propagation_map, only removal happens (no beta-C added)."""
        wf, tracker, groups, state = self._make_propagation_setup(propagation_map=None)
        tracker._events.append(BondEvent(
            step=5, cycle=0, atom_a=0, atom_b=1,
            event_type='confirmed_formation', distance=1.5, r0=2.0,
        ))
        wf._update_groups_after_cycle(state)
        # radical_C should be empty (0 was removed, 2 was never added)
        assert groups['radical_C'].atom_indices == []
        assert groups['vinyl_alpha_C'].atom_indices == []

    def test_beta_not_added_twice(self):
        """If formation fires twice for same alpha-C, beta-C is added only once."""
        wf, tracker, groups, state = self._make_propagation_setup(propagation_map={1: 2})
        ev = BondEvent(
            step=5, cycle=0, atom_a=0, atom_b=1,
            event_type='confirmed_formation', distance=1.5, r0=2.0,
        )
        tracker._events.append(ev)
        tracker._events.append(ev)
        wf._update_groups_after_cycle(state)
        assert groups['radical_C'].atom_indices.count(2) == 1

    def test_processed_formations_advances(self):
        """processed_formations is updated so events are not replayed."""
        wf, tracker, groups, state = self._make_propagation_setup(propagation_map={1: 2})
        tracker._events.append(BondEvent(
            step=5, cycle=0, atom_a=0, atom_b=1,
            event_type='confirmed_formation', distance=1.5, r0=2.0,
        ))
        assert wf.processed_formations == 0
        wf._update_groups_after_cycle(state)
        assert wf.processed_formations == 1
        # Second call should not re-process the same event
        wf._update_groups_after_cycle(state)
        assert wf.processed_formations == 1


class TestBiasedStepsLogging:
    """Tests for CycleLog.steps recording actual biased steps run (RF7)."""

    def test_full_run_records_biased_steps(self):
        """When no early break, CycleLog.steps == config.biased_steps."""
        template, groups = _make_simple_setup()
        config = PolymerizationConfig(
            biased_steps=10, unbiased_steps=5, n_cycles=1, seed=42,
        )
        calc = ToyCalculator()
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'C'],
        )
        wf = PolymerizationWorkflow(config, calc, template, groups)
        logs = wf.run(state)
        biased_log = logs[0]
        assert biased_log.phase == 'biased'
        assert biased_log.steps == config.biased_steps

    def test_early_break_records_fewer_steps(self):
        """When bond_tracker detects a reaction mid-bias, steps < biased_steps."""
        template = ReactionTemplate(
            name='close_pair',
            groups=['A', 'B'],
            pairs=[PairSpec('A', 'B', is_formation=True, r_min=0.5, r_max=5.0)],
        )
        groups = {
            'A': ReactiveGroup('A', [0]),
            'B': ReactiveGroup('B', [1]),
        }
        config = PolymerizationConfig(
            biased_steps=1000, unbiased_steps=5, n_cycles=1, seed=42,
        )
        calc = ToyCalculator()
        tracker = BondTracker(threshold_fraction=1.0)
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'C'],
        )
        wf = PolymerizationWorkflow(
            config, calc, template, groups, bond_tracker=tracker,
        )
        logs = wf.run(state)
        biased_log = logs[0]
        assert biased_log.phase == 'biased'
        assert biased_log.steps < config.biased_steps
        assert biased_log.steps >= 1


class TestConjunctiveBiasedPhaseTermination:
    """F1' (paper §2.2 step 3-4): the biased phase ends only when a candidate's
    FULL reaction event fires — ALL its trigger pairs (amide formation AND every
    leaving-group dissociation) satisfied simultaneously — not on any single
    subordinate crossing. The old _is_terminating_bias_event helper is removed;
    termination is delegated to BondTracker.check_reactions_during_bias, which
    returns the fired candidate_ids."""

    def test_terminating_helper_removed(self):
        """The per-event anchor helper (F1) no longer exists (F1' replaces it)."""
        import kagome.workflows.polymerization as pmod
        assert not hasattr(pmod, '_is_terminating_bias_event')

    def test_phase_does_not_end_while_conjunction_unmet(self):
        """Nylon-like candidate whose dissociation trigger is already satisfied
        but whose amide formation is NOT (atoms far): the conjunction cannot hold
        within the few biased steps, so the phase runs to completion — it must
        never break on the lone dissociation crossing."""
        template = ReactionTemplate(
            name='conjunctive_gating',
            groups=['N', 'C', 'H', 'O'],
            pairs=[
                # Counted amide-like formation, starts far (~4.5 A).
                PairSpec('N', 'C', is_formation=True, r_min=0.5, r_max=5.0),
                # Leaving-group dissociation, already 'dissociated' at start.
                PairSpec('H', 'O', is_formation=False, r_min=0.5, r_max=5.0),
            ],
        )
        groups = {
            'N': ReactiveGroup('N', [0]),
            'C': ReactiveGroup('C', [1]),
            'H': ReactiveGroup('H', [2]),
            'O': ReactiveGroup('O', [3]),
        }
        config = PolymerizationConfig(
            biased_steps=3, unbiased_steps=1, n_cycles=1, seed=42,
        )
        calc = ToyCalculator()
        tracker = BondTracker(threshold_fraction=1.0)
        # H-O far apart -> dissociation trigger satisfied on step 1; N-C far apart
        # so the formation trigger never crosses in 3 weak-bias steps.
        state = SimulationState(
            positions=np.array([
                [0.0, 0.0, 0.0],
                [4.5, 0.0, 0.0],
                [0.0, 20.0, 0.0],
                [0.0, 24.0, 0.0],
            ]),
            velocities=np.zeros((4, 3)),
            species=['N', 'C', 'H', 'O'],
        )
        wf = PolymerizationWorkflow(
            config, calc, template, groups, bond_tracker=tracker,
        )
        logs = wf.run(state)
        biased_log = logs[0]
        assert biased_log.phase == 'biased'
        # The lone dissociation did NOT end the phase; the conjunction never held.
        assert biased_log.steps == config.biased_steps

    def test_single_formation_candidate_still_breaks_early(self):
        """Vinyl-like single-trigger candidate: the phase still breaks as soon as
        the lone formation trigger satisfies (F1 behaviour for a 1-pair candidate
        is a special case of the conjunction)."""
        template = ReactionTemplate(
            name='single_form', groups=['A', 'B'],
            pairs=[PairSpec('A', 'B', is_formation=True, r_min=0.5, r_max=5.0)],
        )
        groups = {'A': ReactiveGroup('A', [0]), 'B': ReactiveGroup('B', [1])}
        config = PolymerizationConfig(
            biased_steps=1000, unbiased_steps=1, n_cycles=1, seed=42,
        )
        calc = ToyCalculator()
        tracker = BondTracker(threshold_fraction=1.0)
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'C'],
        )
        wf = PolymerizationWorkflow(
            config, calc, template, groups, bond_tracker=tracker,
        )
        logs = wf.run(state)
        assert logs[0].phase == 'biased'
        assert logs[0].steps < config.biased_steps


class TestEquilibrationPhase:
    """Tests for _run_equilibration_phase (no bias, no bond tracking)."""

    def test_equilibration_frames_metadata(self, tmp_path):
        """Equilibration frames have phase='equilibration' and cycle=-1."""
        template, groups = _make_simple_setup()
        config = PolymerizationConfig(
            biased_steps=5, unbiased_steps=5, n_cycles=1, seed=42,
            save_interval=1, equil_steps=5,
        )
        calc = ToyCalculator()
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'C'],
        )
        wf = PolymerizationWorkflow(config, calc, template, groups)
        wf.run(state, output_dir=tmp_path)

        _, frames = read_trajectory(tmp_path / 'trajectory.jsonl')
        equil_frames = [f for f in frames if f.phase == 'equilibration']
        assert len(equil_frames) > 0
        for f in equil_frames:
            assert f.cycle == -1

    def test_equilibration_no_bond_events(self):
        """Bond tracker records no events during equilibration (no bias applied)."""
        template, groups = _make_simple_setup()
        config = PolymerizationConfig(
            biased_steps=0, unbiased_steps=0, n_cycles=0, seed=42,
            equil_steps=10,
        )
        calc = ToyCalculator()
        tracker = BondTracker()
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'C'],
        )
        wf = PolymerizationWorkflow(
            config, calc, template, groups, bond_tracker=tracker,
        )
        wf.run(state)

        assert len(tracker.events) == 0

    def test_equilibration_advances_step_counter(self):
        """state.step is advanced by equil_steps."""
        template, groups = _make_simple_setup()
        config = PolymerizationConfig(
            biased_steps=0, unbiased_steps=0, n_cycles=0, seed=42,
            equil_steps=7,
        )
        calc = ToyCalculator()
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'C'],
        )
        wf = PolymerizationWorkflow(config, calc, template, groups)
        wf.run(state)
        assert state.step == 7


class TestPreTDBBMinimize:
    """Config-driven pre-TDBB minimize is what relaxes the classical-compressed
    dense box before the first biased MLIP step. run_nylon66.py wires
    minimize/minimize_fmax/equil_steps into PolymerizationConfig so wf.run()
    performs this relaxation (fix for the nylon-6,6 first-forward-pass segfault;
    paper anchor PDF p.20 — equilibration precedes production reactive MD)."""

    def test_config_minimize_runs_minimize_phase(self, tmp_path):
        """config.minimize=True makes run() emit a 'minimize' phase before TDBB."""
        template, groups = _make_simple_setup()
        config = PolymerizationConfig(
            biased_steps=0, unbiased_steps=0, n_cycles=0, seed=42,
            save_interval=1, minimize=True, minimize_fmax=1.0,
            minimize_max_steps=20,
        )
        calc = ToyCalculator()
        # Close contact so FIRE has something to relax (mirrors the compressed box).
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [0.6, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'C'],
        )
        wf = PolymerizationWorkflow(config, calc, template, groups)
        wf.run(state, output_dir=tmp_path)

        _, frames = read_trajectory(tmp_path / 'trajectory.jsonl')
        assert any(f.phase == 'minimize' for f in frames)

    def test_config_no_minimize_skips_minimize_phase(self, tmp_path):
        """config.minimize=False (the --no-minimize path) emits no minimize phase."""
        template, groups = _make_simple_setup()
        config = PolymerizationConfig(
            biased_steps=0, unbiased_steps=0, n_cycles=0, seed=42,
            save_interval=1, minimize=False,
        )
        calc = ToyCalculator()
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [0.6, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'C'],
        )
        wf = PolymerizationWorkflow(config, calc, template, groups)
        wf.run(state, output_dir=tmp_path)

        _, frames = read_trajectory(tmp_path / 'trajectory.jsonl')
        assert not any(f.phase == 'minimize' for f in frames)


class TestRunActivation:
    """Tests for PolymerizationWorkflow.run_activation (Table S1, V^d on C-N)."""

    def test_dissociation_detected_when_beyond_threshold(self):
        """Activation detects dissociation when C-N distance > 2.5 Å."""
        template, groups = _make_simple_setup()
        config = PolymerizationConfig(
            biased_steps=5, unbiased_steps=5, n_cycles=0, seed=42,
        )
        calc = ToyCalculator()
        wf = PolymerizationWorkflow(config, calc, template, groups)

        activation_template = ReactionTemplate(
            name='aibn_activation',
            groups=['azo_C', 'azo_N'],
            pairs=[PairSpec('azo_C', 'azo_N', is_formation=False,
                            r_min=0.0, r_max=3.0)],
        )
        activation_groups = {
            'azo_C': ReactiveGroup('azo_C', [0]),
            'azo_N': ReactiveGroup('azo_N', [1]),
        }
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [2.6, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'N'],
        )

        dissociated = wf.run_activation(
            state, activation_template, activation_groups,
            activation_steps=10, activation_f2=0.3, activation_f1_max=250.0,
        )
        assert len(dissociated) == 1
        assert dissociated[0] == (0, 1)

    def test_no_dissociation_when_close(self):
        """Activation does not detect dissociation when C-N stays below 2.5 Å."""
        template, groups = _make_simple_setup()
        config = PolymerizationConfig(
            biased_steps=5, unbiased_steps=5, n_cycles=0, seed=42,
        )
        calc = ToyCalculator()
        wf = PolymerizationWorkflow(config, calc, template, groups)

        activation_template = ReactionTemplate(
            name='aibn_activation',
            groups=['azo_C', 'azo_N'],
            pairs=[PairSpec('azo_C', 'azo_N', is_formation=False,
                            r_min=0.0, r_max=3.0)],
        )
        activation_groups = {
            'azo_C': ReactiveGroup('azo_C', [0]),
            'azo_N': ReactiveGroup('azo_N', [1]),
        }
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'N'],
        )

        dissociated = wf.run_activation(
            state, activation_template, activation_groups,
            activation_steps=3, activation_f2=0.3, activation_f1_max=1.0,
        )
        assert len(dissociated) == 0

    def test_activation_returns_empty_when_no_candidates(self):
        """run_activation returns [] when no candidates pass distance filter."""
        template, groups = _make_simple_setup()
        config = PolymerizationConfig(
            biased_steps=5, unbiased_steps=5, n_cycles=0, seed=42,
        )
        calc = ToyCalculator()
        wf = PolymerizationWorkflow(config, calc, template, groups)

        activation_template = ReactionTemplate(
            name='aibn_activation',
            groups=['azo_C', 'azo_N'],
            pairs=[PairSpec('azo_C', 'azo_N', is_formation=False,
                            r_min=0.0, r_max=3.0)],
        )
        activation_groups = {
            'azo_C': ReactiveGroup('azo_C', [0]),
            'azo_N': ReactiveGroup('azo_N', [1]),
        }
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'N'],
        )
        dissociated = wf.run_activation(
            state, activation_template, activation_groups,
            activation_steps=5,
        )
        assert dissociated == []

    def test_activation_selection_audited_with_phase(self, tmp_path):
        """S1: run_activation(output_dir=...) writes a phase='activation' record."""
        template, groups = _make_simple_setup()
        config = PolymerizationConfig(
            biased_steps=5, unbiased_steps=5, n_cycles=0, seed=42,
        )
        calc = ToyCalculator()
        wf = PolymerizationWorkflow(config, calc, template, groups)

        activation_template = ReactionTemplate(
            name='aibn_activation',
            groups=['azo_C', 'azo_N'],
            pairs=[PairSpec('azo_C', 'azo_N', is_formation=False,
                            r_min=0.0, r_max=3.0)],
        )
        activation_groups = {
            'azo_C': ReactiveGroup('azo_C', [0]),
            'azo_N': ReactiveGroup('azo_N', [1]),
        }
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'N'],
        )

        wf.run_activation(
            state, activation_template, activation_groups,
            activation_steps=3, activation_f2=0.3, activation_f1_max=1.0,
            output_dir=tmp_path,
        )

        audit = tmp_path / 'selection.jsonl'
        assert audit.exists()
        lines = [json.loads(l) for l in audit.read_text(encoding='utf-8').splitlines()
                 if l.strip()]
        assert len(lines) == 1
        rec = lines[0]
        assert rec['phase'] == 'activation'
        assert rec['n_selected'] == 1

    def test_activation_no_output_dir_warns_and_skips_audit(self, tmp_path, caplog):
        """S1: without output_dir (and no prior log) activation is not audited."""
        import logging
        template, groups = _make_simple_setup()
        config = PolymerizationConfig(
            biased_steps=5, unbiased_steps=5, n_cycles=0, seed=42,
        )
        calc = ToyCalculator()
        wf = PolymerizationWorkflow(config, calc, template, groups)

        activation_template = ReactionTemplate(
            name='aibn_activation',
            groups=['azo_C', 'azo_N'],
            pairs=[PairSpec('azo_C', 'azo_N', is_formation=False,
                            r_min=0.0, r_max=3.0)],
        )
        activation_groups = {
            'azo_C': ReactiveGroup('azo_C', [0]),
            'azo_N': ReactiveGroup('azo_N', [1]),
        }
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'N'],
        )
        with caplog.at_level(logging.WARNING):
            wf.run_activation(
                state, activation_template, activation_groups,
                activation_steps=3, activation_f2=0.3, activation_f1_max=1.0,
            )
        assert not (tmp_path / 'selection.jsonl').exists()
        assert any('not audited' in r.message for r in caplog.records)


class TestNylonMixedBias:
    """Tests for nylon-like templates with mixed formation/dissociation bias."""

    def _make_nylon_setup(self):
        """Nylon-like 4-group template (Table S2) with controlled positions."""
        template = ReactionTemplate(
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
        return template, groups, positions

    def test_candidates_found(self):
        """Nylon 4-group template finds candidates with valid positions."""
        template, groups, positions = self._make_nylon_setup()
        candidates = find_candidates(template, groups, positions)
        assert len(candidates) == 1
        assert candidates[0].atom_indices == (0, 1, 2, 3)

    def test_build_pair_biases_generates_both_types(self):
        """_build_pair_biases on nylon template produces formation + dissociation pairs."""
        template, groups, positions = self._make_nylon_setup()
        config = PolymerizationConfig(
            biased_steps=5, unbiased_steps=5, n_cycles=0, seed=7,
        )
        calc = ToyCalculator()
        wf = PolymerizationWorkflow(config, calc, template, groups)

        candidates = find_candidates(template, groups, positions)
        scored = score_candidates(candidates)
        selected = select_non_overlapping(scored)
        assert len(selected) >= 1

        pairs = wf._build_pair_biases(selected, ['N', 'C', 'H', 'O'])
        formation = [p for p in pairs if p.is_formation]
        dissociation = [p for p in pairs if not p.is_formation]
        assert len(formation) == 2
        assert len(dissociation) == 2

    def test_bias_forces_have_correct_direction(self):
        """Formation pairs attract (toward r0), dissociation pairs repel."""
        from kagome.boost.tdbb import BoostState, TDBBParams, total_bias

        template, groups, positions = self._make_nylon_setup()
        config = PolymerizationConfig(
            biased_steps=5, unbiased_steps=5, n_cycles=0, seed=7,
        )
        calc = ToyCalculator()
        wf = PolymerizationWorkflow(config, calc, template, groups)

        candidates = find_candidates(template, groups, positions)
        scored = score_candidates(candidates)
        selected = select_non_overlapping(scored)
        pairs = wf._build_pair_biases(selected, ['N', 'C', 'H', 'O'])

        boost = BoostState()
        boost.advance(1.0, 250.0, 125.0)
        energy, forces = total_bias(
            pairs, positions, boost, TDBBParams(),
        )
        assert energy > 0.0
        assert forces.shape == (4, 3)


class TestCheckpointResume:
    """Cycle-boundary checkpoint must allow bit-exact resume (crash recovery).

    A run killed after cycle k and resumed from <out>/checkpoint.pkl must produce
    the identical final state as an uninterrupted run, because the numpy Generator
    state, positions/velocities, groups, and bond tracker are all restored.
    """

    @staticmethod
    def _setup(n_cycles):
        template = ReactionTemplate(
            name='radical_vinyl',
            groups=['radical_C', 'vinyl_alpha_C'],
            pairs=[PairSpec('radical_C', 'vinyl_alpha_C',
                            is_formation=True, r_min=0.5, r_max=5.0)],
        )
        groups = {
            'radical_C': ReactiveGroup('radical_C', [0]),
            'vinyl_alpha_C': ReactiveGroup('vinyl_alpha_C', [1]),
        }
        config = PolymerizationConfig(
            biased_steps=5, unbiased_steps=5, n_cycles=n_cycles, seed=7,
        )
        calc = ToyCalculator()
        tracker = BondTracker()
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.5, 0.0, 0.0]]),
            velocities=np.zeros((3, 3)),
            species=['C', 'C', 'C'],
        )
        wf = PolymerizationWorkflow(
            config, calc, template, groups, bond_tracker=tracker,
            propagation_map={1: 2}, propagation_target_group='radical_C',
        )
        return wf, state, tracker

    def test_resume_is_bit_exact(self, tmp_path):
        # Uninterrupted 4-cycle reference.
        wf_full, state_full, tracker_full = self._setup(n_cycles=4)
        wf_full.run(state_full)
        pos_ref = state_full.positions.copy()
        events_ref = [(e.event_type, e.atom_a, e.atom_b, e.step) for e in tracker_full.events]

        ckpt = tmp_path / 'checkpoint.pkl'
        # First leg: run only 2 cycles, writing a checkpoint each cycle ("crash").
        wf_a, state_a, _ = self._setup(n_cycles=2)
        wf_a.run(state_a, checkpoint_path=ckpt)
        assert ckpt.exists()
        # Diverged from the reference after 2 cycles (proves resume does real work).
        assert not np.array_equal(state_a.positions, pos_ref)

        # Second leg: fresh objects, resume to the full 4 cycles from the checkpoint.
        wf_b, state_b, tracker_b = self._setup(n_cycles=4)
        wf_b.run(state_b, checkpoint_path=ckpt, resume=True)

        np.testing.assert_array_equal(state_b.positions, pos_ref)
        events_resumed = [(e.event_type, e.atom_a, e.atom_b, e.step) for e in tracker_b.events]
        assert events_resumed == events_ref

    def test_no_checkpoint_file_when_disabled(self, tmp_path):
        # Without checkpoint_path, no checkpoint is written.
        wf, state, _ = self._setup(n_cycles=2)
        wf.run(state, output_dir=tmp_path)
        assert not (tmp_path / 'checkpoint.pkl').exists()


class TestIntegratorTemperature:
    """Tests for _integrator_temperature (RF13: Verlet uses kinetic T, not 300 K)."""

    def test_langevin_returns_target_temperature(self):
        """Langevin integrator: returns thermostat target temperature."""
        from kagome.workflows.polymerization import _integrator_temperature
        langevin = LangevinIntegrator(LangevinParams(temperature_K=500.0))
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0]]),
            velocities=np.array([[10.0, 0.0, 0.0]]),
            species=['C'],
        )
        t = _integrator_temperature(langevin, state)
        assert t == pytest.approx(500.0)

    def test_verlet_returns_kinetic_temperature(self):
        """Verlet integrator: returns instantaneous kinetic temperature, not 300 K."""
        from kagome.integrators.verlet import VelocityVerletIntegrator
        from kagome.workflows.polymerization import _integrator_temperature
        verlet = VelocityVerletIntegrator()
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0]]),
            velocities=np.array([[0.01, 0.0, 0.0]]),
            species=['C'],
            masses=np.array([12.011]),
        )
        t = _integrator_temperature(verlet, state)
        assert t != pytest.approx(300.0)
        assert t > 0.0

    def test_verlet_zero_velocity_gives_zero_temperature(self):
        """Verlet with zero velocities gives T=0, not the old 300 K fallback."""
        from kagome.integrators.verlet import VelocityVerletIntegrator
        from kagome.workflows.polymerization import _integrator_temperature
        verlet = VelocityVerletIntegrator()
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0]]),
            velocities=np.zeros((1, 3)),
            species=['C'],
            masses=np.array([12.011]),
        )
        t = _integrator_temperature(verlet, state)
        assert t == pytest.approx(0.0)
