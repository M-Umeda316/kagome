"""Integration test for the polymerization workflow."""
import json

import numpy as np
import pytest

from src.backends.toy import ToyCalculator
from src.integrators.langevin import LangevinIntegrator, LangevinParams
from src.io.readers import read_trajectory
from src.reactive.bonds import BondEvent, BondTracker
from src.reactive.groups import PairSpec, ReactiveGroup, ReactionTemplate
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
