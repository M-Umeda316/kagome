"""Polymerization workflow: biased/unbiased alternation loop.

Paper: arXiv:2511.22874, Mori et al.
Fig. 1: biased (2000 steps) → unbiased (2000 steps) → repeat.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from src.backends.base import Calculator
from src.boost.tdbb import BoostState, PairBias, TDBBParams, target_distance, total_bias
from src.integrators.verlet import Integrator, VelocityVerletIntegrator
from src.io.trajectory import TrajectoryFrame, TrajectoryWriter
from src.reactive.bonds import BondTracker
from src.reactive.groups import ReactiveGroup, ReactionTemplate
from src.reactive.selection import (
    Candidate,
    find_candidates,
    score_candidates,
    select_non_overlapping,
)
from src.workflows.manifest import RunManifest

logger = logging.getLogger(__name__)

VDW_RADII: dict[str, float] = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47,
    'S': 1.80, 'Cl': 1.75, 'Cu': 1.40,
}

ATOMIC_MASSES: dict[str, float] = {
    'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
    'F': 18.998, 'S': 32.06, 'Cl': 35.45, 'Cu': 63.546,
}


@dataclass
class SimulationState:
    positions: NDArray[np.floating]
    velocities: NDArray[np.floating]
    species: list[str]
    cell: NDArray[np.floating] | None = None
    masses: NDArray[np.floating] | None = None
    step: int = 0


@dataclass
class CycleLog:
    cycle: int
    phase: str
    steps: int
    n_candidates: int = 0
    n_selected: int = 0
    bias_energy: float = 0.0


@dataclass
class PolymerizationConfig:
    timestep_fs: float = 0.25
    biased_steps: int = 2000
    unbiased_steps: int = 2000
    n_cycles: int = 10
    tdbb: TDBBParams = field(default_factory=TDBBParams)
    seed: int = 7
    save_interval: int = 0


def masses_from_species(species: list[str]) -> NDArray[np.floating]:
    return np.array([ATOMIC_MASSES.get(s, 12.0) for s in species], dtype=np.float64)


class PolymerizationWorkflow:
    """Alternating biased/unbiased MD loop for polymerization."""

    def __init__(
        self,
        config: PolymerizationConfig,
        calculator: Calculator,
        template: ReactionTemplate,
        groups: dict[str, ReactiveGroup],
        integrator: Integrator | None = None,
        bond_tracker: BondTracker | None = None,
    ) -> None:
        self.config = config
        self.calculator = calculator
        self.template = template
        self.groups = groups
        self.integrator = integrator or VelocityVerletIntegrator()
        self.bond_tracker = bond_tracker
        self.logs: list[CycleLog] = []

    def run(
        self,
        state: SimulationState,
        output_dir: Path | None = None,
        config_path: str = '',
    ) -> list[CycleLog]:
        if output_dir:
            manifest = RunManifest(
                config_path=config_path,
                seed=self.config.seed,
                backend=self.calculator.name,
                output_dir=str(output_dir),
            )
            manifest.save(output_dir / 'manifest.json')

        writer: TrajectoryWriter | None = None
        if output_dir and self.config.save_interval > 0:
            writer = TrajectoryWriter(
                output_dir / 'trajectory.jsonl',
                species=state.species,
                save_interval=self.config.save_interval,
                metadata={'config_path': config_path, 'seed': self.config.seed},
            )

        rng = np.random.default_rng(self.config.seed)

        try:
            for cycle in range(self.config.n_cycles):
                log_biased = self._run_biased_phase(state, cycle, rng, writer)
                self.logs.append(log_biased)
                logger.info(
                    'Cycle %d biased: %d candidates, %d selected, bias_E=%.2f',
                    cycle, log_biased.n_candidates, log_biased.n_selected,
                    log_biased.bias_energy,
                )

                log_unbiased = self._run_unbiased_phase(state, cycle, rng, writer)
                self.logs.append(log_unbiased)
                logger.info('Cycle %d unbiased: %d steps', cycle, log_unbiased.steps)

                self._update_groups_after_cycle(state)
        finally:
            if writer:
                writer.close()
            if self.bond_tracker and output_dir:
                self.bond_tracker.save(output_dir / 'bonds.jsonl')

        return self.logs

    def _run_biased_phase(
        self,
        state: SimulationState,
        cycle: int,
        rng: np.random.Generator,
        writer: TrajectoryWriter | None,
    ) -> CycleLog:
        candidates = find_candidates(self.template, self.groups, state.positions)
        scored = score_candidates(candidates, self.template, state.positions)
        selected = select_non_overlapping(scored)

        active_pairs = self._build_pair_biases(selected, state.species)

        if self.bond_tracker:
            self.bond_tracker.record_attempts(
                active_pairs, state.positions, state.step, cycle,
            )

        boost = BoostState()
        total_bias_energy = 0.0

        for step_in_phase in range(self.config.biased_steps):
            boost.advance(
                self.config.tdbb.gamma,
                self.config.tdbb.f1_max_formation,
                self.config.tdbb.f1_max_dissociation,
            )

            base_energy, base_forces = self.calculator.compute(
                state.positions, state.species, state.cell,
            )
            bias_energy, bias_forces = total_bias(
                active_pairs, state.positions, boost, self.config.tdbb,
            )

            total_forces = base_forces + bias_forces
            total_bias_energy = bias_energy

            self.integrator.step(
                state.positions, state.velocities, total_forces,
                state.masses, self.config.timestep_fs, rng,
            )
            state.step += 1

            if writer and writer.should_write(step_in_phase):
                writer.write_frame(TrajectoryFrame(
                    step=state.step,
                    time_fs=state.step * self.config.timestep_fs,
                    phase='biased',
                    cycle=cycle,
                    energy_base=base_energy,
                    energy_bias=bias_energy,
                    energy_total=base_energy + bias_energy,
                    positions=state.positions.tolist(),
                    n_candidates=len(candidates),
                    n_selected=len(selected),
                ))

        return CycleLog(
            cycle=cycle, phase='biased',
            steps=self.config.biased_steps,
            n_candidates=len(candidates),
            n_selected=len(selected),
            bias_energy=total_bias_energy,
        )

    def _run_unbiased_phase(
        self,
        state: SimulationState,
        cycle: int,
        rng: np.random.Generator,
        writer: TrajectoryWriter | None,
    ) -> CycleLog:
        for step_in_phase in range(self.config.unbiased_steps):
            energy, forces = self.calculator.compute(
                state.positions, state.species, state.cell,
            )

            self.integrator.step(
                state.positions, state.velocities, forces,
                state.masses, self.config.timestep_fs, rng,
            )
            state.step += 1

            if writer and writer.should_write(step_in_phase):
                writer.write_frame(TrajectoryFrame(
                    step=state.step,
                    time_fs=state.step * self.config.timestep_fs,
                    phase='unbiased',
                    cycle=cycle,
                    energy_base=energy,
                    energy_bias=0.0,
                    energy_total=energy,
                    positions=state.positions.tolist(),
                ))

        if self.bond_tracker:
            self.bond_tracker.check_outcomes(state.positions, state.step)

        return CycleLog(
            cycle=cycle, phase='unbiased',
            steps=self.config.unbiased_steps,
        )

    def _update_groups_after_cycle(self, state: SimulationState) -> None:
        if not self.bond_tracker:
            return
        for ev in self.bond_tracker.confirmed_formations():
            for group in self.groups.values():
                if ev.atom_a in group.atom_indices:
                    group.remove_atom(ev.atom_a)
                if ev.atom_b in group.atom_indices:
                    group.remove_atom(ev.atom_b)

    def _build_pair_biases(
        self,
        selected: list[Candidate],
        species: list[str],
    ) -> list[PairBias]:
        pairs: list[PairBias] = []
        label_list = self.template.groups

        for cand in selected:
            for ps in self.template.pairs:
                idx_a_pos = label_list.index(ps.group_a)
                idx_b_pos = label_list.index(ps.group_b)
                atom_a = cand.atom_indices[idx_a_pos]
                atom_b = cand.atom_indices[idx_b_pos]

                vdw_a = VDW_RADII.get(species[atom_a], 1.5)
                vdw_b = VDW_RADII.get(species[atom_b], 1.5)
                r0 = target_distance(
                    np.array([vdw_a, vdw_b]),
                    self.config.tdbb.lambda_vdw,
                )

                pairs.append(PairBias(
                    idx_a=atom_a,
                    idx_b=atom_b,
                    is_formation=ps.is_formation,
                    r0=r0,
                ))
        return pairs
