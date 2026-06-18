"""Polymerization workflow: biased/unbiased alternation loop.

Paper: arXiv:2511.22874, Mori et al.
Fig. 1: biased (2000 steps) → unbiased (2000 steps) → repeat.
"""
from __future__ import annotations

import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Protocol

import numpy as np
from numpy.typing import NDArray

from src.backends.base import Calculator
from src.boost.tdbb import BoostState, PairBias, TDBBParams, target_distance, total_bias
from src.geometry import minimum_image
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
from src.integrators.init_velocities import instant_temperature_K
from src.workflows.manifest import RunManifest, _normalize_value

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
    # Closest any biased formation pair got during the phase (Å). Evidence for
    # whether [3,6]-window pairs reach the bias-capture shell / bonding distance.
    min_pair_distance: float = float('inf')


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


def _integrator_temperature(integrator: object, state: SimulationState) -> float:
    """Return the temperature for MC barostat acceptance.

    Langevin: use the thermostat target temperature (constant).
    Other (e.g. Verlet/NVE): use the instantaneous kinetic temperature.
    """
    from src.integrators.langevin import LangevinIntegrator
    if isinstance(integrator, LangevinIntegrator):
        return integrator.params.temperature_K
    return instant_temperature_K(state.velocities, state.masses)


# ---------------------------------------------------------------------------
# PostCycleUpdater: injectable group-update strategy (RF4)
# ---------------------------------------------------------------------------

class PostCycleUpdater(Protocol):
    """Protocol for post-cycle group updates after confirmed formations."""

    def update(
        self,
        groups: dict[str, ReactiveGroup],
        tracker: BondTracker | None,
        state: SimulationState,
    ) -> None: ...


class DefaultPostCycleUpdater:
    """Remove reacted atoms from their groups after confirmed formation."""

    def __init__(self) -> None:
        self._processed_formations: int = 0

    def update(
        self,
        groups: dict[str, ReactiveGroup],
        tracker: BondTracker | None,
        state: SimulationState,
    ) -> None:
        if not tracker:
            return
        formations = tracker.confirmed_formations()
        for ev in formations[self._processed_formations:]:
            for group in groups.values():
                if ev.atom_a in group.atom_indices:
                    group.remove_atom(ev.atom_a)
                if ev.atom_b in group.atom_indices:
                    group.remove_atom(ev.atom_b)
        self._processed_formations = len(formations)


class VinylChainPropagationUpdater:
    """Vinyl radical chain propagation: beta-C becomes new radical after formation.

    Paper anchor: Table S1 — after radical + vinyl_alpha_C formation,
    the monomer's beta-C becomes the new radical site.
    """

    def __init__(
        self,
        propagation_map: dict[int, int],
        propagation_target_group: str = 'radical_C',
        chain_c_map: dict[int, int] | None = None,
    ) -> None:
        self.propagation_map = propagation_map
        self.propagation_target_group = propagation_target_group
        self.chain_c_map = chain_c_map if chain_c_map is not None else {}
        self._processed_formations: int = 0

    def update(
        self,
        groups: dict[str, ReactiveGroup],
        tracker: BondTracker | None,
        state: SimulationState,
    ) -> None:
        if not tracker:
            return
        formations = tracker.confirmed_formations()
        for ev in formations[self._processed_formations:]:
            for group in groups.values():
                if ev.atom_a in group.atom_indices:
                    group.remove_atom(ev.atom_a)
                if ev.atom_b in group.atom_indices:
                    group.remove_atom(ev.atom_b)

            chain_c_group = groups.get('chain_C')
            if chain_c_group is not None and ev.atom_a in self.chain_c_map:
                old_chain_c = self.chain_c_map.pop(ev.atom_a)
                chain_c_group.remove_atom(old_chain_c)

            beta_c_group = groups.get('vinyl_beta_C')
            if beta_c_group is not None and ev.atom_b in self.propagation_map:
                beta_idx = self.propagation_map[ev.atom_b]
                beta_c_group.remove_atom(beta_idx)

            if self.propagation_map and ev.atom_b in self.propagation_map:
                beta_idx = self.propagation_map[ev.atom_b]
                target = groups.get(self.propagation_target_group)
                if target is not None and beta_idx not in target.atom_indices:
                    target.atom_indices.append(beta_idx)
                    logger.info(
                        'Chain propagation: atom %d (beta-C) → %s',
                        beta_idx, self.propagation_target_group,
                    )

                if chain_c_group is not None:
                    chain_c_group.atom_indices.append(ev.atom_b)
                    self.chain_c_map[beta_idx] = ev.atom_b

        self._processed_formations = len(formations)


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
        chain_c_map: dict[int, int] | None = None,
        updater: PostCycleUpdater | None = None,
    ) -> None:
        self.config = config
        self.calculator = calculator
        self.template = template
        self.groups = groups
        self.integrator = integrator or VelocityVerletIntegrator()
        self.bond_tracker = bond_tracker
        self.barostat = barostat
        self.logs: list[CycleLog] = []

        if updater is not None:
            self._updater = updater
        elif propagation_map:
            self._updater = VinylChainPropagationUpdater(
                propagation_map=propagation_map,
                propagation_target_group=propagation_target_group,
                chain_c_map=chain_c_map,
            )
        else:
            self._updater = DefaultPostCycleUpdater()

    @property
    def _processed_formations(self) -> int:
        return self._updater._processed_formations

    def run(
        self,
        state: SimulationState,
        output_dir: Path | None = None,
        config_path: str = '',
        n_monomers: int | None = None,
    ) -> list[CycleLog]:
        if output_dir:
            effective_params = _normalize_value(asdict(self.config))
            effective_params['backend'] = self.calculator.name
            effective_params['candidate_r_min'] = self.template.pairs[0].r_min if self.template.pairs else None
            effective_params['candidate_r_max'] = self.template.pairs[0].r_max if self.template.pairs else None
            manifest = RunManifest(
                config_path=config_path,
                seed=self.config.seed,
                backend=self.calculator.name,
                output_dir=str(output_dir),
                extra=effective_params,
            )
            manifest.save(output_dir / 'manifest.json')

        n_reactive_sites = n_monomers if n_monomers is not None else sum(len(g.atom_indices) for g in self.groups.values())

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
                    'Cycle %d biased: %d candidates, %d selected, bias_E=%.2f, '
                    'min_pair_dist=%.2f A',
                    cycle, log_biased.n_candidates, log_biased.n_selected,
                    log_biased.bias_energy, log_biased.min_pair_distance,
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

    # ------------------------------------------------------------------
    # Consolidated MD step (RF4): single implementation for all phases
    # ------------------------------------------------------------------

    def _md_step(
        self,
        state: SimulationState,
        current_forces: NDArray[np.floating],
        dt: float,
        rng: np.random.Generator,
        step_in_phase: int,
        *,
        active_pairs: list[PairBias] | None = None,
        boost: BoostState | None = None,
        tdbb: TDBBParams | None = None,
        enable_barostat: bool = True,
    ) -> tuple[float, float, NDArray[np.floating]]:
        """Single MD step: pre_force → compute → [bias] → post_force → [barostat].

        Returns (base_energy, bias_energy, current_forces).
        """
        self.integrator.pre_force(
            state.positions, state.velocities, current_forces,
            state.masses, dt, rng, state.cell,
        )

        base_energy, base_forces = self.calculator.compute(
            state.positions, state.species, state.cell,
        )

        bias_energy = 0.0
        if active_pairs is not None and boost is not None and tdbb is not None:
            bias_energy, bias_forces = total_bias(
                active_pairs, state.positions, boost, tdbb, state.cell,
            )
            current_forces = base_forces + bias_forces
        else:
            current_forces = base_forces

        self.integrator.post_force(
            state.velocities, current_forces, state.masses, dt,
        )
        state.step += 1

        if (enable_barostat
                and self.barostat is not None
                and state.cell is not None
                and self.barostat.should_attempt(step_in_phase)):
            accepted, new_base_e, new_base_f = self.barostat.try_step(
                state.positions, state.species, state.cell,
                base_energy, self.calculator, rng,
                _integrator_temperature(self.integrator, state),
            )
            if accepted and new_base_f is not None:
                base_energy = new_base_e
                if active_pairs is not None and boost is not None and tdbb is not None:
                    bias_energy, bias_forces = total_bias(
                        active_pairs, state.positions, boost, tdbb, state.cell,
                    )
                    current_forces = new_base_f + bias_forces
                else:
                    current_forces = new_base_f

        return base_energy, bias_energy, current_forces

    # ------------------------------------------------------------------
    # Phase implementations — each delegates to _md_step
    # ------------------------------------------------------------------

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
            energy, _, current_forces = self._md_step(
                state, current_forces, dt, rng, step_in_phase,
            )

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
                    temperature_K=instant_temperature_K(state.velocities, state.masses),
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

        form_pairs = [p for p in active_pairs if p.is_formation]
        min_form_dist = float('inf')

        steps_run = 0
        for step_in_phase in range(self.config.biased_steps):
            steps_run = step_in_phase + 1
            boost.advance(
                self.config.tdbb.gamma,
                self.config.tdbb.f1_max_formation,
                self.config.tdbb.f1_max_dissociation,
            )

            base_energy, bias_energy, current_forces = self._md_step(
                state, current_forces, dt, rng, step_in_phase,
                active_pairs=active_pairs, boost=boost, tdbb=self.config.tdbb,
            )
            last_bias_energy = bias_energy

            for _p in form_pairs:
                _r = float(np.linalg.norm(minimum_image(
                    state.positions[_p.idx_b] - state.positions[_p.idx_a], state.cell)))
                if _r < min_form_dist:
                    min_form_dist = _r

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
                    temperature_K=instant_temperature_K(state.velocities, state.masses),
                ))

            if self.bond_tracker is not None:
                events = self.bond_tracker.check_reactions_during_bias(
                    active_pairs, state.positions, state.step, cycle, state.cell,
                )
                if events:
                    logger.info(
                        'Cycle %d biased: reaction event at step %d (%d pair(s)) '
                        '- ending biased phase', cycle, step_in_phase + 1, len(events),
                    )
                    break

        return CycleLog(
            cycle=cycle, phase='biased',
            steps=steps_run,
            n_candidates=len(candidates),
            n_selected=len(selected),
            bias_energy=last_bias_energy,
            min_pair_distance=min_form_dist,
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
            energy, _, current_forces = self._md_step(
                state, current_forces, dt, rng, step_in_phase,
            )

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
                    temperature_K=instant_temperature_K(state.velocities, state.masses),
                ))

        if self.bond_tracker:
            self.bond_tracker.check_outcomes(state.positions, state.step, state.cell)

        return CycleLog(
            cycle=cycle, phase='unbiased',
            steps=self.config.unbiased_steps,
        )

    def _update_groups_after_cycle(self, state: SimulationState) -> None:
        self._updater.update(self.groups, self.bond_tracker, state)

    def run_activation(
        self,
        state: SimulationState,
        activation_template: ReactionTemplate,
        activation_groups: dict[str, ReactiveGroup],
        activation_steps: int = 2000,
        activation_f2: float = 0.3,
        activation_f1_max: float = 250.0,
        rng: np.random.Generator | None = None,
    ) -> list[tuple[int, int]]:
        """Run AIBN activation phase: V^d on azo C-N bonds until dissociation.

        activation_f2: V^d Gaussian width — must be small enough that the
        force peak (r=1/√(2·f2)) overlaps the C-N bond distance (~1.49 Å).
        activation_f1_max: peak V^d amplitude.  Must be large enough that the
        effective potential PES+V^d is monotonically repulsive at the C-N bond.
        With OrbMol-v2 C-N barrier ≈39 kcal/mol, f2=0.3+f1≥200 guarantees this.

        Returns list of (C_idx, N_idx) pairs that dissociated.

        Paper anchor: Table S1 — Activation row, V^d applied to C-N azo bonds.
        """
        if rng is None:
            rng = np.random.default_rng(self.config.seed)

        candidates = find_candidates(
            activation_template, activation_groups, state.positions, state.cell,
        )
        scored = score_candidates(
            candidates, activation_template, state.positions, state.cell,
        )
        selected = select_non_overlapping(scored)

        if not selected:
            logger.warning('Activation: no C-N candidates found.')
            return []

        pairs = self._build_pair_biases(
            selected, state.species, template=activation_template,
        )

        logger.info(
            'Activation: %d candidates, %d selected, %d V^d pairs',
            len(candidates), len(selected), len(pairs),
        )
        for i, p in enumerate(pairs):
            r = float(np.linalg.norm(minimum_image(
                state.positions[p.idx_b] - state.positions[p.idx_a], state.cell)))
            logger.info('  V^d pair %d: atoms (%d, %d) r=%.2f A', i, p.idx_a, p.idx_b, r)

        boost = BoostState()
        dt = self.config.timestep_fs
        tdbb = TDBBParams(
            f1_max_formation=self.config.tdbb.f1_max_formation,
            f1_max_dissociation=activation_f1_max,
            f2=activation_f2,
            gamma=self.config.tdbb.gamma,
            lambda_vdw=self.config.tdbb.lambda_vdw,
        )

        base_energy, base_forces = self.calculator.compute(
            state.positions, state.species, state.cell,
        )
        bias_energy, bias_forces = total_bias(
            pairs, state.positions, boost, tdbb, state.cell,
        )
        current_forces = base_forces + bias_forces

        dissoc_threshold = 2.5
        dissociated: list[tuple[int, int]] = []

        for step_in_phase in range(activation_steps):
            boost.advance(tdbb.gamma, tdbb.f1_max_formation, tdbb.f1_max_dissociation)

            base_energy, bias_energy, current_forces = self._md_step(
                state, current_forces, dt, rng, step_in_phase,
                active_pairs=pairs, boost=boost, tdbb=tdbb,
                enable_barostat=False,
            )

            for p in pairs:
                r = float(np.linalg.norm(minimum_image(
                    state.positions[p.idx_b] - state.positions[p.idx_a], state.cell)))
                if (step_in_phase + 1) % 500 == 0 or step_in_phase == 0:
                    logger.info(
                        'Activation step %d: (%d,%d) r=%.3f A, f1_d=%.1f, bias_E=%.1f',
                        step_in_phase + 1, p.idx_a, p.idx_b, r,
                        boost.f1_dissociation, bias_energy,
                    )
                if r > dissoc_threshold and (p.idx_a, p.idx_b) not in dissociated:
                    dissociated.append((p.idx_a, p.idx_b))
                    logger.info(
                        'Activation: C-N dissociation at step %d, atoms (%d, %d), r=%.2f A',
                        step_in_phase + 1, p.idx_a, p.idx_b, r,
                    )

            if len(dissociated) == len(pairs):
                logger.info(
                    'Activation: all %d C-N bonds dissociated at step %d',
                    len(pairs), step_in_phase + 1,
                )
                break

        logger.info('Activation done: %d/%d dissociated', len(dissociated), len(pairs))
        return dissociated

    def _build_pair_biases(
        self,
        selected: list[Candidate],
        species: list[str],
        *,
        template: ReactionTemplate | None = None,
    ) -> list[PairBias]:
        template = template or self.template
        pairs: list[PairBias] = []
        label_list = template.groups

        for cand in selected:
            for ps in template.pairs:
                if ps.constraint_only:
                    continue
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
