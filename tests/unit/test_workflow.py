"""Integration test for the polymerization workflow."""
import json

import numpy as np
import pytest

from src.backends.toy import ToyCalculator
from src.boost.tdbb import PairBias, TDBBParams
from src.integrators.langevin import LangevinIntegrator, LangevinParams
from src.io.readers import read_trajectory
from src.reactive.bonds import BondEvent, BondTracker
from src.reactive.groups import PairSpec, ReactiveGroup, ReactionTemplate
from src.reactive.selection import (
    find_candidates,
    score_candidates,
    select_non_overlapping,
)
from src.workflows.polymerization import (
    PolymerizationConfig,
    PolymerizationWorkflow,
    SimulationState,
    masses_from_species,
)


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

    def test_with_masses(self):
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
            masses=masses_from_species(['C', 'C']),
        )
        wf = PolymerizationWorkflow(config, calc, template, groups)
        logs = wf.run(state)
        assert state.step == 10

    def test_with_langevin(self):
        template, groups = _make_simple_setup()
        config = PolymerizationConfig(
            biased_steps=5,
            unbiased_steps=5,
            n_cycles=1,
            seed=7,
        )
        calc = ToyCalculator()
        langevin = LangevinIntegrator(LangevinParams(temperature_K=300.0))
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'C'],
        )
        wf = PolymerizationWorkflow(config, calc, template, groups, integrator=langevin)
        logs = wf.run(state)
        assert state.step == 10

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
        """_processed_formations is updated so events are not replayed."""
        wf, tracker, groups, state = self._make_propagation_setup(propagation_map={1: 2})
        tracker._events.append(BondEvent(
            step=5, cycle=0, atom_a=0, atom_b=1,
            event_type='confirmed_formation', distance=1.5, r0=2.0,
        ))
        assert wf._processed_formations == 0
        wf._update_groups_after_cycle(state)
        assert wf._processed_formations == 1
        # Second call should not re-process the same event
        wf._update_groups_after_cycle(state)
        assert wf._processed_formations == 1


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
                         r_min=0.0, r_max=100.0),
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
        scored = score_candidates(candidates, template, positions)
        selected = select_non_overlapping(scored)
        assert len(selected) >= 1

        pairs = wf._build_pair_biases(selected, ['N', 'C', 'H', 'O'])
        formation = [p for p in pairs if p.is_formation]
        dissociation = [p for p in pairs if not p.is_formation]
        assert len(formation) == 2
        assert len(dissociation) == 2

    def test_bias_forces_have_correct_direction(self):
        """Formation pairs attract (toward r0), dissociation pairs repel."""
        from src.boost.tdbb import BoostState, TDBBParams, total_bias

        template, groups, positions = self._make_nylon_setup()
        config = PolymerizationConfig(
            biased_steps=5, unbiased_steps=5, n_cycles=0, seed=7,
        )
        calc = ToyCalculator()
        wf = PolymerizationWorkflow(config, calc, template, groups)

        candidates = find_candidates(template, groups, positions)
        scored = score_candidates(candidates, template, positions)
        selected = select_non_overlapping(scored)
        pairs = wf._build_pair_biases(selected, ['N', 'C', 'H', 'O'])

        boost = BoostState()
        boost.advance(1.0, 250.0, 125.0)
        energy, forces = total_bias(
            pairs, positions, boost, TDBBParams(),
        )
        assert energy > 0.0
        assert forces.shape == (4, 3)


class TestIntegratorTemperature:
    """Tests for _integrator_temperature (RF13: Verlet uses kinetic T, not 300 K)."""

    def test_langevin_returns_target_temperature(self):
        """Langevin integrator: returns thermostat target temperature."""
        from src.workflows.polymerization import _integrator_temperature
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
        from src.integrators.verlet import VelocityVerletIntegrator
        from src.workflows.polymerization import _integrator_temperature
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
        from src.integrators.verlet import VelocityVerletIntegrator
        from src.workflows.polymerization import _integrator_temperature
        verlet = VelocityVerletIntegrator()
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0]]),
            velocities=np.zeros((1, 3)),
            species=['C'],
            masses=np.array([12.011]),
        )
        t = _integrator_temperature(verlet, state)
        assert t == pytest.approx(0.0)
