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
from src.integrators.mc_barostat import MCBarostat
from src.integrators.minimize import FireParams, fire_minimize
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
from src.units import FORCE_CONV, KB
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
    # Pre-TDBB relaxation (paper anchor: PDF p.20 — equilibration precedes
    # reactive production).  Disabled by default to preserve legacy behaviour;
    # run scripts opt in.
    minimize: bool = False
    minimize_fmax: float = 1.0
    minimize_max_steps: int = 500
    equil_steps: int = 0


def masses_from_species(species: list[str]) -> NDArray[np.floating]:
    masses = []
    for s in species:
        if s not in ATOMIC_MASSES:
            logger.warning('Unknown element %r — using fallback mass 12.0 amu', s)
        masses.append(ATOMIC_MASSES.get(s, 12.0))
    return np.array(masses, dtype=np.float64)


def _instant_temperature(
    velocities: NDArray[np.floating],
    masses: NDArray[np.floating] | None,
) -> float:
    """Instantaneous kinetic temperature from velocities.

    KE[kcal/mol] = 0.5 * sum(m * v^2) / FORCE_CONV
    T[K] = 2 * KE / (3 * N * KB)
    """
    n = velocities.shape[0]
    if n == 0:
        return 0.0
    if masses is not None:
        ke_amu = 0.5 * float(np.sum(masses[:, np.newaxis] * velocities ** 2))
    else:
        ke_amu = 0.5 * float(np.sum(velocities ** 2))
    ke_kcal = ke_amu / FORCE_CONV
    return 2.0 * ke_kcal / (3.0 * n * KB)


def _integrator_temperature(integrator: object) -> float:
    """Extract target temperature from a Langevin integrator, or return 300 K."""
    from src.integrators.langevin import LangevinIntegrator
    if isinstance(integrator, LangevinIntegrator):
        return integrator.params.temperature_K
    return 300.0


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
        barostat: MCBarostat | None = None,
        propagation_map: dict[int, int] | None = None,
        propagation_target_group: str = 'radical_C',
    ) -> None:
        self.config = config
        self.calculator = calculator
        self.template = template
        self.groups = groups
        self.integrator = integrator or VelocityVerletIntegrator()
        self.bond_tracker = bond_tracker
        self.barostat = barostat
        self.propagation_map = propagation_map or {}
        self.propagation_target_group = propagation_target_group
        self.logs: list[CycleLog] = []
        self._processed_formations: int = 0

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

        n_reactive_sites = sum(len(g.atom_indices) for g in self.groups.values())

        writer: TrajectoryWriter | None = None
        if output_dir and self.config.save_interval > 0:
            writer = TrajectoryWriter(
                output_dir / 'trajectory.jsonl',
                species=state.species,
                save_interval=self.config.save_interval,
                metadata={'config_path': config_path, 'seed': self.config.seed},
                n_reactive_sites=n_reactive_sites,
            )

        rng = np.random.default_rng(self.config.seed)

        try:
            if self.config.minimize:
                self._minimize(state)
            if self.config.equil_steps > 0:
                self._run_equilibration_phase(state, rng, writer)

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

    def _minimize(self, state: SimulationState) -> None:
        """Relax close contacts in the initial structure before dynamics.

        Paper anchor: PDF p.20 — production reactive MD follows equilibration.
        Grid-packed structures carry intermolecular clashes whose large forces
        would spike the temperature; FIRE minimization removes them first.
        """
        logger.info(
            'Pre-TDBB energy minimization (FIRE, fmax=%.2f, max_steps=%d)...',
            self.config.minimize_fmax, self.config.minimize_max_steps,
        )
        result = fire_minimize(
            state.positions, state.species, state.cell, self.calculator,
            FireParams(
                fmax_kcal_mol_A=self.config.minimize_fmax,
                max_steps=self.config.minimize_max_steps,
            ),
        )
        state.positions[:] = result.positions

    def _run_equilibration_phase(
        self,
        state: SimulationState,
        rng: np.random.Generator,
        writer: TrajectoryWriter | None,
    ) -> None:
        """Unbiased NPT/NVT equilibration before the first TDBB cycle.

        Paper anchor: PDF p.20 — "Equilibration simulations were performed in
        the NPT ensemble ... Production simulations using reactive acceleration
        MD were then carried out".  No bias and no bond tracking; frames are
        labelled phase='equilibration', cycle=-1.
        """
        dt = self.config.timestep_fs
        energy, forces = self.calculator.compute(
            state.positions, state.species, state.cell,
        )
        current_forces = forces

        for step_in_phase in range(self.config.equil_steps):
            self.integrator.pre_force(
                state.positions, state.velocities, current_forces,
                state.masses, dt, rng, state.cell,
            )
            energy, forces = self.calculator.compute(
                state.positions, state.species, state.cell,
            )
            current_forces = forces
            self.integrator.post_force(
                state.velocities, current_forces, state.masses, dt,
            )
            state.step += 1

            if (self.barostat is not None
                    and state.cell is not None
                    and self.barostat.should_attempt(step_in_phase)):
                accepted, new_e, new_f = self.barostat.try_step(
                    state.positions, state.species, state.cell,
                    energy, self.calculator, rng,
                    _integrator_temperature(self.integrator),
                )
                if accepted and new_f is not None:
                    energy, current_forces = new_e, new_f

            if writer and writer.should_write(step_in_phase):
                writer.write_frame(TrajectoryFrame(
                    step=state.step,
                    time_fs=state.step * dt,
                    phase='equilibration',
                    cycle=-1,
                    energy_base=energy,
                    energy_bias=0.0,
                    energy_total=energy,
                    positions=state.positions.tolist(),
                    temperature_K=_instant_temperature(state.velocities, state.masses),
                ))

        logger.info('Equilibration: %d steps complete', self.config.equil_steps)

    def _run_biased_phase(
        self,
        state: SimulationState,
        cycle: int,
        rng: np.random.Generator,
        writer: TrajectoryWriter | None,
    ) -> CycleLog:
        candidates = find_candidates(
            self.template, self.groups, state.positions, state.cell,
        )
        scored = score_candidates(candidates, self.template, state.positions, state.cell)
        selected = select_non_overlapping(scored)

        active_pairs = self._build_pair_biases(selected, state.species)

        if self.bond_tracker:
            self.bond_tracker.record_attempts(
                active_pairs, state.positions, state.step, cycle, state.cell,
            )

        boost = BoostState()
        dt = self.config.timestep_fs
        last_bias_energy = 0.0

        base_energy, base_forces = self.calculator.compute(
            state.positions, state.species, state.cell,
        )
        bias_energy, bias_forces = total_bias(
            active_pairs, state.positions, boost, self.config.tdbb, state.cell,
        )
        current_forces = base_forces + bias_forces

        for step_in_phase in range(self.config.biased_steps):
            boost.advance(
                self.config.tdbb.gamma,
                self.config.tdbb.f1_max_formation,
                self.config.tdbb.f1_max_dissociation,
            )

            self.integrator.pre_force(
                state.positions, state.velocities, current_forces,
                state.masses, dt, rng, state.cell,
            )

            base_energy, base_forces = self.calculator.compute(
                state.positions, state.species, state.cell,
            )
            bias_energy, bias_forces = total_bias(
                active_pairs, state.positions, boost, self.config.tdbb, state.cell,
            )
            current_forces = base_forces + bias_forces
            last_bias_energy = bias_energy

            self.integrator.post_force(
                state.velocities, current_forces, state.masses, dt,
            )
            state.step += 1

            if (self.barostat is not None
                    and state.cell is not None
                    and self.barostat.should_attempt(step_in_phase)):
                accepted, new_base_e, new_base_f = self.barostat.try_step(
                    state.positions, state.species, state.cell,
                    base_energy, self.calculator, rng,
                    _integrator_temperature(self.integrator),
                )
                if accepted and new_base_f is not None:
                    base_energy, base_forces = new_base_e, new_base_f
                    bias_energy, bias_forces = total_bias(
                        active_pairs, state.positions, boost, self.config.tdbb, state.cell,
                    )
                    current_forces = base_forces + bias_forces
                    last_bias_energy = bias_energy

            if writer and writer.should_write(step_in_phase):
                writer.write_frame(TrajectoryFrame(
                    step=state.step,
                    time_fs=state.step * dt,
                    phase='biased',
                    cycle=cycle,
                    energy_base=base_energy,
                    energy_bias=bias_energy,
                    energy_total=base_energy + bias_energy,
                    positions=state.positions.tolist(),
                    n_candidates=len(candidates),
                    n_selected=len(selected),
                    temperature_K=_instant_temperature(state.velocities, state.masses),
                ))

        return CycleLog(
            cycle=cycle, phase='biased',
            steps=self.config.biased_steps,
            n_candidates=len(candidates),
            n_selected=len(selected),
            bias_energy=last_bias_energy,
        )

    def _run_unbiased_phase(
        self,
        state: SimulationState,
        cycle: int,
        rng: np.random.Generator,
        writer: TrajectoryWriter | None,
    ) -> CycleLog:
        dt = self.config.timestep_fs

        energy, forces = self.calculator.compute(
            state.positions, state.species, state.cell,
        )
        current_forces = forces

        for step_in_phase in range(self.config.unbiased_steps):
            self.integrator.pre_force(
                state.positions, state.velocities, current_forces,
                state.masses, dt, rng, state.cell,
            )

            energy, forces = self.calculator.compute(
                state.positions, state.species, state.cell,
            )
            current_forces = forces

            self.integrator.post_force(
                state.velocities, current_forces, state.masses, dt,
            )
            state.step += 1

            if (self.barostat is not None
                    and state.cell is not None
                    and self.barostat.should_attempt(step_in_phase)):
                accepted, new_e, new_f = self.barostat.try_step(
                    state.positions, state.species, state.cell,
                    energy, self.calculator, rng,
                    _integrator_temperature(self.integrator),
                )
                if accepted and new_f is not None:
                    energy, current_forces = new_e, new_f

            if writer and writer.should_write(step_in_phase):
                writer.write_frame(TrajectoryFrame(
                    step=state.step,
                    time_fs=state.step * dt,
                    phase='unbiased',
                    cycle=cycle,
                    energy_base=energy,
                    energy_bias=0.0,
                    energy_total=energy,
                    positions=state.positions.tolist(),
                    temperature_K=_instant_temperature(state.velocities, state.masses),
                ))

        if self.bond_tracker:
            self.bond_tracker.check_outcomes(state.positions, state.step, state.cell)

        return CycleLog(
            cycle=cycle, phase='unbiased',
            steps=self.config.unbiased_steps,
        )

    def _update_groups_after_cycle(self, state: SimulationState) -> None:
        if not self.bond_tracker:
            return
        formations = self.bond_tracker.confirmed_formations()
        for ev in formations[self._processed_formations:]:
            for group in self.groups.values():
                if ev.atom_a in group.atom_indices:
                    group.remove_atom(ev.atom_a)
                if ev.atom_b in group.atom_indices:
                    group.remove_atom(ev.atom_b)
            # Chain propagation: beta-C of reacted monomer becomes new radical site.
            # ev.atom_b = vinyl_alpha_C; propagation_map[alpha] = beta.
            if self.propagation_map and ev.atom_b in self.propagation_map:
                beta_idx = self.propagation_map[ev.atom_b]
                target = self.groups.get(self.propagation_target_group)
                if target is not None and beta_idx not in target.atom_indices:
                    target.atom_indices.append(beta_idx)
                    logger.info(
                        'Chain propagation: atom %d (beta-C) → %s',
                        beta_idx, self.propagation_target_group,
                    )
        self._processed_formations = len(formations)

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

                sp_a, sp_b = species[atom_a], species[atom_b]
                if sp_a not in VDW_RADII:
                    logger.warning('Unknown element %r — using fallback vdW radius 1.5 Å', sp_a)
                if sp_b not in VDW_RADII:
                    logger.warning('Unknown element %r — using fallback vdW radius 1.5 Å', sp_b)
                vdw_a = VDW_RADII.get(sp_a, 1.5)
                vdw_b = VDW_RADII.get(sp_b, 1.5)
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
