"""Integration test for the polymerization workflow."""
import json

import numpy as np
import pytest

from src.backends.toy import ToyCalculator
from src.integrators.langevin import LangevinIntegrator, LangevinParams
from src.io.readers import read_trajectory
from src.reactive.bonds import BondTracker
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
