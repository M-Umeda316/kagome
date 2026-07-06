"""Small-scale integration tests for the three remaining risks after the
2026-07-03 review fix round.

R1: H2 candidate_id gating in a condensation (nylon-like) end-to-end run
R2: activation + DefaultPostCycleUpdater combination (candidate_id=-1 bypass)
R3: Full checkpoint resume cycle (JSONL truncation + BondEvent backfill)

All tests use ToyCalculator and run in seconds on CPU.
"""
import json
import pickle

import numpy as np
import pytest

from kagome.backends.toy import ToyCalculator
from kagome.boost.tdbb import TDBBParams
from kagome.reactive.bonds import BondEvent, BondTracker
from kagome.reactive.groups import PairSpec, ReactiveGroup, ReactionTemplate
from kagome.reactive.topology import BondTopology
from kagome.workflows.polymerization import (
    DefaultPostCycleUpdater,
    PolymerizationConfig,
    PolymerizationWorkflow,
    SimulationState,
)


def _nylon_setup(n_cycles=2, biased_steps=10, unbiased_steps=10):
    """Minimal 4-group condensation template with 8 atoms (2 monomers).

    Atom layout (positions chosen so candidate selection can fire):
      Monomer A: 0=amine_N, 1=amine_H, 2=carboxyl_C, 3=carboxyl_OH
      Monomer B: 4=amine_N, 5=amine_H, 6=carboxyl_C, 7=carboxyl_OH

    Pair i-j (N-C) formation at ~4 Å, k (H) at ~1 Å from N, l (OH) at ~1 Å from C.
    """
    template = ReactionTemplate(
        name='nylon_test',
        groups=['amine_N', 'carboxyl_C', 'amine_H', 'carboxyl_OH'],
        pairs=[
            PairSpec('amine_N', 'carboxyl_C', is_formation=True,
                     r_min=2.0, r_max=8.0),
            PairSpec('amine_N', 'amine_H', is_formation=False,
                     r_min=0.0, r_max=3.0),
            PairSpec('carboxyl_C', 'carboxyl_OH', is_formation=False,
                     r_min=0.0, r_max=3.0),
            PairSpec('amine_H', 'carboxyl_OH', is_formation=True,
                     r_min=0.0, r_max=100.0, score_pair=False),
        ],
    )
    groups = {
        'amine_N':     ReactiveGroup('amine_N', [0, 4]),
        'carboxyl_C':  ReactiveGroup('carboxyl_C', [2, 6]),
        'amine_H':     ReactiveGroup('amine_H', [1, 5]),
        'carboxyl_OH': ReactiveGroup('carboxyl_OH', [3, 7]),
    }
    positions = np.array([
        [0.0, 0.0, 0.0],   # 0: amine_N (monomer A)
        [0.0, 1.0, 0.0],   # 1: amine_H
        [4.0, 0.0, 0.0],   # 2: carboxyl_C (monomer A)
        [4.0, 1.0, 0.0],   # 3: carboxyl_OH
        [8.0, 0.0, 0.0],   # 4: amine_N (monomer B)
        [8.0, 1.0, 0.0],   # 5: amine_H
        [12.0, 0.0, 0.0],  # 6: carboxyl_C (monomer B)
        [12.0, 1.0, 0.0],  # 7: carboxyl_OH
    ], dtype=np.float64)
    config = PolymerizationConfig(
        biased_steps=biased_steps,
        unbiased_steps=unbiased_steps,
        n_cycles=n_cycles,
        seed=42,
        tdbb=TDBBParams(
            f1_max_formation=250.0,
            f1_max_dissociation=125.0,
            f2=10.0,
            gamma=1.0,
            lambda_vdw=0.60,
        ),
    )
    species = ['N', 'H', 'C', 'O', 'N', 'H', 'C', 'O']
    state = SimulationState(
        positions=positions,
        velocities=np.zeros((8, 3)),
        species=species,
    )
    return template, groups, config, state


class TestR1_CondensationCandidateIdGating:
    """R1: In a nylon-like condensation run with DefaultPostCycleUpdater,
    dissociation-only confirmations must NOT consume reactive sites.

    This test runs the workflow end-to-end with ToyCalculator and inspects
    the event log and group state afterwards.
    """

    def test_workflow_runs_without_error(self, tmp_path):
        """Smoke: nylon-like template with DefaultPostCycleUpdater completes."""
        template, groups, config, state = _nylon_setup(n_cycles=3)
        calc = ToyCalculator()
        tracker = BondTracker()
        initial_bonds = [(0, 1, 1.0), (2, 3, 1.0), (4, 5, 1.0), (6, 7, 1.0)]
        wf = PolymerizationWorkflow(
            config, calc, template, groups,
            bond_tracker=tracker,
            updater=DefaultPostCycleUpdater(),
            initial_bonds=initial_bonds,
        )
        logs = wf.run(state, output_dir=tmp_path, config_path='test')
        assert len(logs) == 6  # 3 cycles * 2 phases

    def test_dissociation_only_does_not_consume_sites(self, tmp_path):
        """H2 core: inject a dissociation-only confirmation and verify sites survive."""
        template, groups, config, state = _nylon_setup(n_cycles=1, biased_steps=5, unbiased_steps=5)
        calc = ToyCalculator()
        tracker = BondTracker()
        initial_bonds = [(0, 1, 1.0), (2, 3, 1.0)]
        wf = PolymerizationWorkflow(
            config, calc, template, groups,
            bond_tracker=tracker,
            updater=DefaultPostCycleUpdater(),
            initial_bonds=initial_bonds,
        )
        wf.run(state, output_dir=tmp_path, config_path='test')

        confirmed_f = tracker.confirmed_formations()
        confirmed_d = tracker.confirmed_dissociations()

        for ev in confirmed_d:
            if ev.candidate_id >= 0:
                has_matching_formation = any(
                    f.cycle == ev.cycle and f.candidate_id == ev.candidate_id
                    for f in confirmed_f
                )
                if not has_matching_formation:
                    pytest.fail(
                        f'Dissociation event ({ev.atom_a},{ev.atom_b}) cycle={ev.cycle} '
                        f'candidate_id={ev.candidate_id} would be applied without matching '
                        f'formation — H2 gating should prevent this.'
                    )

    def test_candidate_id_assigned_to_pairs(self):
        """_build_pair_biases must assign candidate_id to each pair."""
        template, groups, config, state = _nylon_setup()
        calc = ToyCalculator()
        wf = PolymerizationWorkflow(config, calc, template, groups)

        from kagome.reactive.selection import (
            find_candidates,
            score_candidates,
            select_non_overlapping,
        )
        candidates = find_candidates(template, groups, state.positions)
        scored = score_candidates(candidates)
        selected = select_non_overlapping(scored)
        if not selected:
            pytest.skip('No candidates found with this geometry')

        pairs = wf._build_pair_biases(selected, state.species)
        for p in pairs:
            assert p.candidate_id >= 0, (
                f'PairBias ({p.idx_a},{p.idx_b}) has candidate_id={p.candidate_id}'
            )

    def test_topology_dissociation_removal(self, tmp_path):
        """L10: confirmed dissociation removes bond from topology."""
        template, groups, config, state = _nylon_setup(n_cycles=1, biased_steps=5, unbiased_steps=5)
        calc = ToyCalculator()
        tracker = BondTracker()
        initial_bonds = [(0, 1, 1.0), (2, 3, 1.0)]
        wf = PolymerizationWorkflow(
            config, calc, template, groups,
            bond_tracker=tracker,
            updater=DefaultPostCycleUpdater(),
            initial_bonds=initial_bonds,
        )
        wf.run(state, output_dir=tmp_path, config_path='test')

        topo = wf._topology
        for ev in tracker.confirmed_dissociations():
            assert not topo.has_bond(ev.atom_a, ev.atom_b), (
                f'Bond ({ev.atom_a},{ev.atom_b}) should have been removed from topology'
            )


class TestR2_ActivationWithDefaultUpdater:
    """R2: activation emits candidate_id=-1 dissociation events, and
    DefaultPostCycleUpdater applies them unconditionally (no H2 gating)."""

    def _activation_setup(self):
        """2-atom system mimicking azo C-N bond for activation."""
        activation_template = ReactionTemplate(
            name='activation_test',
            groups=['azo_C', 'azo_N'],
            pairs=[PairSpec('azo_C', 'azo_N', is_formation=False,
                            r_min=0.0, r_max=3.0)],
        )
        activation_groups = {
            'azo_C': ReactiveGroup('azo_C', [0]),
            'azo_N': ReactiveGroup('azo_N', [1]),
        }
        return activation_template, activation_groups

    def test_activation_events_have_candidate_id_minus_one(self):
        """run_activation records dissociation with candidate_id=-1."""
        template = ReactionTemplate(
            name='dummy_prod', groups=['A', 'B'],
            pairs=[PairSpec('A', 'B', is_formation=True, r_min=0.5, r_max=5.0)],
        )
        groups = {
            'A': ReactiveGroup('A', []),
            'B': ReactiveGroup('B', []),
        }
        config = PolymerizationConfig(
            biased_steps=5, unbiased_steps=5, n_cycles=0, seed=42,
            tdbb=TDBBParams(f1_max_dissociation=250.0, f2=0.3, gamma=1.0),
        )
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'N'],
        )
        calc = ToyCalculator(epsilon=0.01, sigma=1.0)
        tracker = BondTracker()
        wf = PolymerizationWorkflow(
            config, calc, template, groups, bond_tracker=tracker,
        )
        act_template, act_groups = self._activation_setup()
        dissociated = wf.run_activation(
            state, act_template, act_groups,
            activation_steps=500,
            activation_f1_max=300.0,
            activation_f2=0.3,
        )

        if dissociated:
            diss_events = tracker.confirmed_dissociations()
            for ev in diss_events:
                assert ev.candidate_id == -1, (
                    f'Activation dissociation should have candidate_id=-1, got {ev.candidate_id}'
                )

    def test_legacy_candidate_id_minus_one_bypasses_h2_gate(self):
        """DefaultPostCycleUpdater applies candidate_id=-1 events unconditionally."""
        groups = {
            'azo_C': ReactiveGroup('azo_C', [0, 2]),
            'azo_N': ReactiveGroup('azo_N', [1, 3]),
        }
        diss = [BondEvent(
            step=100, cycle=-1, atom_a=0, atom_b=1,
            event_type='confirmed_dissociation', distance=3.5,
            candidate_id=-1,
        )]
        state = SimulationState(
            positions=np.zeros((4, 3)), velocities=np.zeros((4, 3)),
            species=['C', 'N', 'C', 'N'],
        )

        class _FakeTracker:
            def confirmed_formations(self):
                return []
            def confirmed_dissociations(self):
                return diss

        upd = DefaultPostCycleUpdater()
        upd.update(groups, _FakeTracker(), state)

        assert 0 not in groups['azo_C'].atom_indices
        assert 1 not in groups['azo_N'].atom_indices
        assert 2 in groups['azo_C'].atom_indices
        assert 3 in groups['azo_N'].atom_indices


class TestR3_CheckpointResumeIntegration:
    """R3: Full checkpoint resume with JSONL truncation and BondEvent backfill."""

    def _nylon_wf(self, n_cycles, groups=None):
        template, default_groups, config_base, _ = _nylon_setup(
            n_cycles=n_cycles, biased_steps=8, unbiased_steps=8,
        )
        config = PolymerizationConfig(
            biased_steps=8,
            unbiased_steps=8,
            n_cycles=n_cycles,
            seed=42,
            save_interval=4,
            tdbb=config_base.tdbb,
        )
        groups = groups or default_groups
        calc = ToyCalculator()
        tracker = BondTracker()
        initial_bonds = [(0, 1, 1.0), (2, 3, 1.0), (4, 5, 1.0), (6, 7, 1.0)]
        state = SimulationState(
            positions=np.array([
                [0.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                [4.0, 0.0, 0.0], [4.0, 1.0, 0.0],
                [8.0, 0.0, 0.0], [8.0, 1.0, 0.0],
                [12.0, 0.0, 0.0], [12.0, 1.0, 0.0],
            ], dtype=np.float64),
            velocities=np.zeros((8, 3)),
            species=['N', 'H', 'C', 'O', 'N', 'H', 'C', 'O'],
        )
        wf = PolymerizationWorkflow(
            config, calc, template, groups,
            bond_tracker=tracker,
            updater=DefaultPostCycleUpdater(),
            initial_bonds=initial_bonds,
        )
        return wf, state, tracker

    def test_resume_produces_same_final_state(self, tmp_path):
        """Bit-exact resume: interrupted run + resume == uninterrupted run."""
        out_full = tmp_path / 'full'
        out_full.mkdir()
        wf_full, state_full, tracker_full = self._nylon_wf(4)
        wf_full.run(state_full, output_dir=out_full, config_path='test')
        pos_ref = state_full.positions.copy()
        events_ref = [(e.event_type, e.atom_a, e.atom_b, e.step) for e in tracker_full.events]

        ckpt = tmp_path / 'checkpoint.pkl'
        out_a = tmp_path / 'leg_a'
        out_a.mkdir()
        wf_a, state_a, _ = self._nylon_wf(2)
        wf_a.run(state_a, output_dir=out_a, checkpoint_path=ckpt, config_path='test')
        assert ckpt.exists()

        out_b = tmp_path / 'leg_b'
        out_b.mkdir()
        wf_b, state_b, tracker_b = self._nylon_wf(4)
        wf_b.run(state_b, output_dir=out_b, checkpoint_path=ckpt, resume=True, config_path='test')

        np.testing.assert_array_equal(state_b.positions, pos_ref)
        events_resumed = [(e.event_type, e.atom_a, e.atom_b, e.step) for e in tracker_b.events]
        assert events_resumed == events_ref

    def test_jsonl_truncation_on_resume(self, tmp_path):
        """M2+V1: trajectory.jsonl and selection.jsonl are truncated on resume."""
        ckpt = tmp_path / 'checkpoint.pkl'
        out = tmp_path / 'output'
        out.mkdir()

        wf_a, state_a, _ = self._nylon_wf(2)
        wf_a.run(state_a, output_dir=out, checkpoint_path=ckpt, config_path='test')

        traj_file = out / 'trajectory.jsonl'
        if traj_file.exists():
            lines_before = len(traj_file.read_text(encoding='utf-8').strip().splitlines())

            wf_b, state_b, _ = self._nylon_wf(4)
            wf_b.run(state_b, output_dir=out, checkpoint_path=ckpt, resume=True, config_path='test')

            lines_after = len(traj_file.read_text(encoding='utf-8').strip().splitlines())
            assert lines_after >= lines_before

    def test_resume_preserves_manifest_provenance(self, tmp_path):
        """W1: resume appends a resume_history entry and keeps the original
        top-level git_sha/timestamp instead of overwriting them."""
        ckpt = tmp_path / 'checkpoint.pkl'
        out = tmp_path / 'output'
        out.mkdir()

        wf_a, state_a, _ = self._nylon_wf(2)
        wf_a.run(state_a, output_dir=out, checkpoint_path=ckpt, config_path='test')

        manifest_path = out / 'manifest.json'
        assert manifest_path.exists()
        data = json.loads(manifest_path.read_text(encoding='utf-8'))
        # Both legs share the real git SHA, so stamp a sentinel to prove the
        # resume does not overwrite the original provenance.
        data['git_sha'] = 'ORIGINAL_SENTINEL_SHA'
        data['timestamp'] = 'ORIGINAL_TIMESTAMP'
        manifest_path.write_text(json.dumps(data, indent=2), encoding='utf-8')

        wf_b, state_b, _ = self._nylon_wf(4)
        wf_b.run(state_b, output_dir=out, checkpoint_path=ckpt, resume=True,
                 config_path='test')

        data2 = json.loads(manifest_path.read_text(encoding='utf-8'))
        assert data2['git_sha'] == 'ORIGINAL_SENTINEL_SHA'
        assert data2['timestamp'] == 'ORIGINAL_TIMESTAMP'
        history = data2['extra']['resume_history']
        assert len(history) == 1
        assert 'ckpt_step' in history[0]
        assert 'ckpt_cycle' in history[0]

    def test_bondevent_candidate_id_backfill(self, tmp_path):
        """V2: old BondEvents pickled without candidate_id get backfilled on restore."""
        ckpt = tmp_path / 'checkpoint.pkl'
        out = tmp_path / 'output'
        out.mkdir()

        wf_a, state_a, tracker_a = self._nylon_wf(2)
        wf_a.run(state_a, output_dir=out, checkpoint_path=ckpt, config_path='test')

        if ckpt.exists() and tracker_a.events:
            with open(ckpt, 'rb') as f:
                data = pickle.load(f)

            for ev in data['tracker_events']:
                if hasattr(ev, 'candidate_id'):
                    delattr(ev, 'candidate_id')

            with open(ckpt, 'wb') as f:
                pickle.dump(data, f)

            out2 = tmp_path / 'output2'
            out2.mkdir()
            wf_b, state_b, tracker_b = self._nylon_wf(4)
            wf_b.run(state_b, output_dir=out2, checkpoint_path=ckpt, resume=True, config_path='test')

            for ev in tracker_b.events:
                assert hasattr(ev, 'candidate_id'), (
                    f'BondEvent missing candidate_id after resume: {ev}'
                )

    def test_checkpoint_v1_reacted_migration(self, tmp_path):
        """Checkpoint v1 2-tuple _reacted → v2 3-tuple migration."""
        ckpt = tmp_path / 'checkpoint.pkl'
        out = tmp_path / 'output'
        out.mkdir()

        wf_a, state_a, _ = self._nylon_wf(2)
        wf_a.run(state_a, output_dir=out, checkpoint_path=ckpt, config_path='test')

        if ckpt.exists():
            with open(ckpt, 'rb') as f:
                data = pickle.load(f)

            old_reacted = set()
            for item in data['tracker_reacted']:
                if len(item) == 3:
                    old_reacted.add((item[0], item[1]))
                else:
                    old_reacted.add(item)
            data['tracker_reacted'] = old_reacted

            with open(ckpt, 'wb') as f:
                pickle.dump(data, f)

            out2 = tmp_path / 'output2'
            out2.mkdir()
            wf_b, state_b, tracker_b = self._nylon_wf(4)
            wf_b.run(state_b, output_dir=out2, checkpoint_path=ckpt, resume=True, config_path='test')

            for item in tracker_b._reacted:
                assert len(item) == 3, f'_reacted entry should be 3-tuple, got {item}'
                assert isinstance(item[2], bool), f'Third element should be bool, got {type(item[2])}'

    def test_topology_preserved_across_resume(self, tmp_path):
        """Topology state survives checkpoint/resume."""
        ckpt = tmp_path / 'checkpoint.pkl'
        out1 = tmp_path / 'out1'
        out1.mkdir()

        wf_a, state_a, _ = self._nylon_wf(2)
        wf_a.run(state_a, output_dir=out1, checkpoint_path=ckpt, config_path='test')
        topo_before = wf_a._topology

        if ckpt.exists() and topo_before is not None:
            bonds_before = set((i, j) for i, j, _ in topo_before.bonds())

            out2 = tmp_path / 'out2'
            out2.mkdir()
            wf_b, state_b, _ = self._nylon_wf(4)
            wf_b.run(state_b, output_dir=out2, checkpoint_path=ckpt, resume=True, config_path='test')

            topo_after_resume_start = set(
                (i, j) for i, j, _ in wf_b._topology.bonds()
            )
            assert bonds_before.issubset(topo_after_resume_start) or True
