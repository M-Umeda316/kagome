"""Polymerization workflow: biased/unbiased alternation loop.

Paper: arXiv:2511.22874, Mori et al.
Fig. 1: biased (2000 steps) → unbiased (2000 steps) → repeat.
"""
from __future__ import annotations

import json
import logging
import pickle
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Protocol

import numpy as np
from numpy.typing import NDArray

from kagome.backends.base import Calculator
from kagome.boost.tdbb import BoostState, PairBias, TDBBParams, target_distance, total_bias_fast
from kagome.diagnostics import StepWatchdog
from kagome.geometry import minimum_image, validated_box
from kagome.integrators.mc_barostat import MCBarostat
from kagome.integrators.minimize import FireParams, fire_minimize
from kagome.integrators.verlet import Integrator, VelocityVerletIntegrator
from kagome.io.trajectory import TrajectoryFrame, TrajectoryWriter
from kagome.reactive.bonds import BondTracker, is_dissociated
from kagome.reactive.groups import ReactiveGroup, ReactionTemplate
from kagome.reactive.topology import (
    BondTopology, apply_vinyl_addition, vinyl_addition_over_coordinates,
)
from kagome.reactive.selection import (
    Candidate,
    audited_selection,
    find_candidates,
    score_candidates,
    select_non_overlapping,
)
from kagome.integrators.init_velocities import instant_temperature_K
from kagome.workflows.manifest import RunManifest, _normalize_value

logger = logging.getLogger(__name__)

VDW_RADII: dict[str, float] = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47,
    'S': 1.80, 'Cl': 1.75, 'Cu': 1.40,
}

ATOMIC_MASSES: dict[str, float] = {
    'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
    'F': 18.998, 'S': 32.06, 'Cl': 35.45, 'Cu': 63.546,
}

# ---------------------------------------------------------------------------
# Cycle-boundary checkpointing (crash recovery for long TDBB runs)
# ---------------------------------------------------------------------------
# A checkpoint captures everything the cycle loop mutates so a killed run can
# resume at the next cycle instead of from scratch (build + activation + all
# completed cycles). The static parts (species, template, calculator/weights)
# are NOT stored — the resuming process rebuilds them and overrides the dynamic
# state below. Resume is bit-exact because the numpy Generator state is saved.
CHECKPOINT_VERSION = 2


def load_checkpoint(path: Path) -> dict:
    """Load a cycle-boundary checkpoint written by :func:`save_checkpoint`."""
    with open(path, 'rb') as f:
        data = pickle.load(f)
    if data.get('version') != CHECKPOINT_VERSION:
        logger.warning(
            'Checkpoint version %s != current %s; resume may be incompatible.',
            data.get('version'), CHECKPOINT_VERSION,
        )
    return data


def save_checkpoint(
    path: Path,
    *,
    next_cycle: int,
    state: SimulationState,
    groups: dict[str, ReactiveGroup],
    updater: object,
    tracker: BondTracker | None,
    rng: np.random.Generator,
    logs: list[CycleLog],
    extra: dict | None = None,
    topology: BondTopology | None = None,
    topo_processed_formations: int = 0,
) -> None:
    """Atomically write a checkpoint at a cycle boundary.

    Written to ``path.with_suffix('.tmp')`` then renamed, so a crash mid-write
    never corrupts the previous good checkpoint.
    """
    data = {
        'version': CHECKPOINT_VERSION,
        'next_cycle': next_cycle,
        'positions': np.asarray(state.positions),
        'velocities': np.asarray(state.velocities),
        'cell': None if state.cell is None else np.asarray(state.cell),
        'step': state.step,
        'groups': {label: list(g.atom_indices) for label, g in groups.items()},
        'updater_processed_formations': getattr(updater, '_processed_formations', 0),
        'updater_processed_dissociations': getattr(updater, '_processed_dissociations', 0),
        'updater_chain_c_map': dict(getattr(updater, 'chain_c_map', {})),
        'tracker_events': list(tracker._events) if tracker else [],
        'tracker_reacted': set(tracker._reacted) if tracker else set(),
        'rng_state': rng.bit_generator.state,
        'logs': list(logs),
        'extra': extra or {},
        # Explicit bond topology (specs/decisions.md 2026-07-02): resume must
        # continue the same connectivity, not re-derive it from scratch.
        'topology_bonds': (topology.bonds() if topology is not None else None),
        'topo_processed_formations': topo_processed_formations,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    with open(tmp, 'wb') as f:
        pickle.dump(data, f)
    tmp.replace(path)


def _pair_distance(
    positions: NDArray[np.floating],
    idx_a: int,
    idx_b: int,
    cell: NDArray[np.floating] | None,
) -> float:
    """Minimum-image distance (Å) between two atoms. Single source for the
    several diagnostic distance recomputations in this module (RF24)."""
    return float(np.linalg.norm(
        minimum_image(positions[idx_b] - positions[idx_a], cell)
    ))


@dataclass
class SimulationState:
    """MD シミュレーションの瞬間状態。

    ワークフロー全体を通じて更新され、積分器・バックエンド・ワークフローの
    間で受け渡される中心データ構造。
    """

    positions: NDArray[np.floating]
    velocities: NDArray[np.floating]
    species: list[str]
    cell: NDArray[np.floating] | None = None
    masses: NDArray[np.floating] | None = None
    step: int = 0


@dataclass
class CycleLog:
    """1 サイクル中の 1 フェーズ (biased/unbiased) の診断情報。"""

    cycle: int
    phase: str
    steps: int
    n_candidates: int = 0
    n_selected: int = 0
    bias_energy: float = 0.0
    min_pair_distance: float = float('inf')


@dataclass
class PolymerizationConfig:
    """重合ワークフローの全パラメータ。

    TDBB パラメータ (tdbb)、積分ステップ数、前処理 (minimize/equil) を含む。
    run スクリプトから CLI 引数を経由して構築される。
    """

    timestep_fs: float = 0.25
    biased_steps: int = 2000
    unbiased_steps: int = 2000
    n_cycles: int = 10
    tdbb: TDBBParams = field(default_factory=TDBBParams)
    seed: int = 7
    save_interval: int = 0
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
    from kagome.integrators.langevin import LangevinIntegrator
    if isinstance(integrator, LangevinIntegrator):
        return integrator.params.temperature_K
    return instant_temperature_K(state.velocities, state.masses)


# ---------------------------------------------------------------------------
# PostCycleUpdater: injectable group-update strategy (RF4)
# ---------------------------------------------------------------------------

class PostCycleUpdater(Protocol):
    """Protocol for post-cycle group updates after confirmed formations."""

    @property
    def processed_formations(self) -> int: ...

    def update(
        self,
        groups: dict[str, ReactiveGroup],
        tracker: BondTracker | None,
        state: SimulationState,
    ) -> None: ...


class DefaultPostCycleUpdater:
    """Remove reacted atoms from their groups after confirmed formation AND
    dissociation.

    For step-growth condensation (nylon-6,6, Table S2) the leaving groups
    (amine_H, carboxyl_OH) and reacted centres are freed by V^d *dissociation*
    events; without consuming those events the same atoms stay selectable in
    later cycles and the condensation topology never advances (RF15). Chain
    extension itself needs no special handling here: build_nylon66_system
    registers *both* termini of every monomer, so growth continues on the
    remaining ends. See specs/decisions.md 2026-06-20 RF15.
    """

    def __init__(self) -> None:
        self._processed_formations: int = 0
        self._processed_dissociations: int = 0

    @property
    def processed_formations(self) -> int:
        return self._processed_formations

    @property
    def processed_dissociations(self) -> int:
        return self._processed_dissociations

    @staticmethod
    def _remove_pair(groups: dict[str, ReactiveGroup], ev) -> None:
        for group in groups.values():
            if ev.atom_a in group.atom_indices:
                group.remove_atom(ev.atom_a)
            if ev.atom_b in group.atom_indices:
                group.remove_atom(ev.atom_b)

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
            self._remove_pair(groups, ev)
        self._processed_formations = len(formations)

        # Consume dissociations so freed leaving groups / reacted centres are not
        # re-selected next cycle (remove_atom is idempotent if also consumed above).
        dissociations = tracker.confirmed_dissociations()
        for ev in dissociations[self._processed_dissociations:]:
            self._remove_pair(groups, ev)
        self._processed_dissociations = len(dissociations)


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

    @property
    def processed_formations(self) -> int:
        return self._processed_formations

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
        initial_bonds: list[tuple[int, int, float]] | None = None,
    ) -> None:
        self.config = config
        self.calculator = calculator
        self.template = template
        self.groups = groups
        self.integrator = integrator or VelocityVerletIntegrator()
        self.bond_tracker = bond_tracker
        self.barostat = barostat
        self.logs: list[CycleLog] = []
        # Per-cycle candidate-selection audit (RF18); set in run() when output_dir given.
        self._selection_log: Path | None = None

        # Explicit bond topology (specs/decisions.md 2026-07-02): when initial
        # bonds are supplied we track connectivity through the run so the
        # trajectory carries real bonds (no distance-inferred spurious/over-
        # valent bonds in the viewer). None => topology output disabled.
        self._propagation_map = propagation_map or {}
        self._topology = (
            BondTopology.from_bonds(initial_bonds) if initial_bonds else None)
        self._topo_processed_formations = 0
        self._topology_log: Path | None = None

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
    def processed_formations(self) -> int:
        return self._updater.processed_formations

    def run(
        self,
        state: SimulationState,
        output_dir: Path | None = None,
        config_path: str = '',
        n_monomers: int | None = None,
        checkpoint_path: Path | None = None,
        resume: bool = False,
        checkpoint_extra: dict | None = None,
    ) -> list[CycleLog]:
        """biased/unbiased 交互ループを n_cycles 回実行する。

        output_dir を指定すると trajectory.jsonl, bonds.jsonl, manifest.json,
        summary.json を出力する。戻り値は各フェーズの CycleLog リスト。

        checkpoint_path を指定すると各サイクル境界で状態を保存する。resume=True
        かつ checkpoint が存在すれば minimize/equil/初期frame を skip し、保存
        サイクルの次から継続する（trajectory/selection は追記）。checkpoint_extra
        は checkpoint に保存する追加情報（例: production spin）。
        """
        n_reactive_sites = n_monomers if n_monomers is not None else sum(len(g.atom_indices) for g in self.groups.values())

        # Resume from a cycle-boundary checkpoint: restore the mutated dynamic
        # state and continue at the saved cycle. The static parts (species,
        # template, calculator) were rebuilt by the caller; here we override the
        # parts the cycle loop mutates. minimize/equil/initial-frame are skipped.
        _ckpt = None
        start_cycle = 0
        resuming = False
        if resume and checkpoint_path is not None and Path(checkpoint_path).exists():
            _ckpt = load_checkpoint(Path(checkpoint_path))
            resuming = True
            start_cycle = int(_ckpt['next_cycle'])
            state.positions = np.array(_ckpt['positions'])
            state.velocities = np.array(_ckpt['velocities'])
            state.cell = None if _ckpt['cell'] is None else np.array(_ckpt['cell'])
            state.step = int(_ckpt['step'])
            for label, idxs in _ckpt['groups'].items():
                if label in self.groups:
                    self.groups[label].atom_indices[:] = list(idxs)
            self._updater._processed_formations = _ckpt['updater_processed_formations']
            if hasattr(self._updater, '_processed_dissociations'):
                self._updater._processed_dissociations = _ckpt['updater_processed_dissociations']
            if hasattr(self._updater, 'chain_c_map'):
                self._updater.chain_c_map = dict(_ckpt['updater_chain_c_map'])
            if self.bond_tracker is not None:
                self.bond_tracker._events = list(_ckpt['tracker_events'])
                raw_reacted = _ckpt['tracker_reacted']
                # Migrate v1 checkpoints: (int, int) -> (int, int, True)
                migrated: set[tuple[int, int, bool]] = set()
                for item in raw_reacted:
                    if len(item) == 2:
                        migrated.add((*item, True))
                    else:
                        migrated.add(tuple(item))
                self.bond_tracker._reacted = migrated
            if self._topology is not None and _ckpt.get('topology_bonds') is not None:
                self._topology = BondTopology.from_bonds(_ckpt['topology_bonds'])
                self._topo_processed_formations = _ckpt.get(
                    'topo_processed_formations', 0)
            self.logs = list(_ckpt['logs'])
            logger.info(
                'Resuming from checkpoint: starting at cycle %d of %d '
                '(%d confirmed formations so far)',
                start_cycle, self.config.n_cycles,
                len(self.bond_tracker.confirmed_formations()) if self.bond_tracker else 0,
            )

        if output_dir:
            effective_params = _normalize_value(asdict(self.config))
            effective_params['backend'] = self.calculator.name
            # Resolved weights identity (RF17): two runs with the same backend
            # name but different weights are otherwise indistinguishable.
            effective_params['model_id'] = getattr(
                self.calculator, 'model_id', self.calculator.name)
            # alpha(t) denominator (RF17): also lives in the trajectory header, but
            # record it in the manifest so provenance is complete without the JSONL.
            effective_params['n_reactive_sites'] = n_reactive_sites
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

            # Per-cycle candidate-selection audit (RF18): "why was X dropped for Y"
            # must be reconstructable from artifacts, not just n_candidates counts.
            self._selection_log = output_dir / 'selection.jsonl'
            if not resuming:
                self._selection_log.write_text('', encoding='utf-8')  # truncate prior runs

            # Explicit bond-topology history (specs/decisions.md 2026-07-02).
            if self._topology is not None:
                self._topology_log = output_dir / 'topology.jsonl'
                if not resuming:
                    self._topology_log.write_text('', encoding='utf-8')
                    self._write_topology_snapshot(state, cycle=-1)  # initial

        writer: TrajectoryWriter | None = None
        if output_dir and self.config.save_interval > 0:
            writer = TrajectoryWriter(
                output_dir / 'trajectory.jsonl',
                species=state.species,
                save_interval=self.config.save_interval,
                metadata={'config_path': config_path, 'seed': self.config.seed},
                n_reactive_sites=n_reactive_sites,
                append=resuming,  # keep prior frames when resuming
                initial_bonds=(self._topology.bonds()
                               if self._topology is not None else None),
            )

        rng = np.random.default_rng(self.config.seed)
        if _ckpt is not None:
            # Bit-exact continuation: restore the generator to its cycle-boundary state.
            rng.bit_generator.state = _ckpt['rng_state']

        try:
            if not resuming:
                if writer:
                    writer.write_frame(TrajectoryFrame(
                        step=state.step,
                        time_fs=0.0,
                        phase='initial',
                        cycle=-1,
                        energy_base=0.0,
                        energy_bias=0.0,
                        energy_total=0.0,
                        positions=state.positions.tolist(),
                    ))
                if self.config.minimize:
                    self._minimize(state, writer)
                if self.config.equil_steps > 0:
                    self._run_equilibration_phase(state, rng, writer)

            for cycle in range(start_cycle, self.config.n_cycles):
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

                self._update_groups_after_cycle(state, cycle)

                if checkpoint_path is not None:
                    # Save AFTER the group update so resume starts cleanly at the
                    # next cycle. Atomic write — a crash mid-save keeps the prior one.
                    save_checkpoint(
                        Path(checkpoint_path),
                        next_cycle=cycle + 1,
                        state=state,
                        groups=self.groups,
                        updater=self._updater,
                        tracker=self.bond_tracker,
                        rng=rng,
                        logs=self.logs,
                        extra=checkpoint_extra,
                        topology=self._topology,
                        topo_processed_formations=self._topo_processed_formations,
                    )
        finally:
            if writer:
                writer.close()
            if self.bond_tracker and output_dir:
                self.bond_tracker.save(output_dir / 'bonds.jsonl')

        return self.logs

    def _minimize(
        self,
        state: SimulationState,
        writer: TrajectoryWriter | None = None,
    ) -> None:
        """Relax close contacts in the initial structure before dynamics.

        Paper anchor: PDF p.20 — production reactive MD follows equilibration.
        Grid-packed structures carry intermolecular clashes whose large forces
        would spike the temperature; FIRE minimization removes them first.
        """
        logger.info(
            'Pre-TDBB energy minimization (FIRE, fmax=%.2f, max_steps=%d)...',
            self.config.minimize_fmax, self.config.minimize_max_steps,
        )
        save_interval = max(1, self.config.save_interval)

        def _on_step(step: int, pos, energy: float, fmax: float) -> None:
            if writer and step % save_interval == 0:
                writer.write_frame(TrajectoryFrame(
                    step=step,
                    time_fs=0.0,
                    phase='minimize',
                    cycle=-1,
                    energy_base=energy,
                    energy_bias=0.0,
                    energy_total=energy,
                    positions=pos.tolist(),
                ))

        result = fire_minimize(
            state.positions, state.species, state.cell, self.calculator,
            FireParams(
                fmax_kcal_mol_A=self.config.minimize_fmax,
                max_steps=self.config.minimize_max_steps,
            ),
            on_step=_on_step,
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
        box: NDArray[np.floating] | None = None,
    ) -> tuple[float, float, NDArray[np.floating], list[float] | None]:
        """Single MD step: pre_force → compute → [bias] → post_force → [barostat].

        Returns (base_energy, bias_energy, current_forces, pair_distances).
        """
        self.integrator.pre_force(
            state.positions, state.velocities, current_forces,
            state.masses, dt, rng, state.cell,
        )

        base_energy, base_forces = self.calculator.compute(
            state.positions, state.species, state.cell,
        )

        bias_energy = 0.0
        pair_dists: list[float] | None = None
        if active_pairs is not None and boost is not None and tdbb is not None:
            bias_energy, bias_forces, pair_dists = total_bias_fast(
                active_pairs, state.positions, boost, tdbb, box,
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
            bias_fn = None
            if active_pairs is not None and boost is not None and tdbb is not None:
                def bias_fn(pos, bx):
                    return total_bias_fast(active_pairs, pos, boost, tdbb, bx)[0]
            accepted, new_base_e, new_base_f = self.barostat.try_step(
                state.positions, state.species, state.cell,
                base_energy, self.calculator, rng,
                _integrator_temperature(self.integrator, state),
                bias_energy_fn=bias_fn,
            )
            if accepted and new_base_f is not None:
                base_energy = new_base_e
                if box is not None and state.cell is not None:
                    box[0] = state.cell[0, 0]
                    box[1] = state.cell[1, 1]
                    box[2] = state.cell[2, 2]
                if active_pairs is not None and boost is not None and tdbb is not None:
                    bias_energy, bias_forces, pair_dists = total_bias_fast(
                        active_pairs, state.positions, boost, tdbb, box,
                    )
                    current_forces = new_base_f + bias_forces
                else:
                    current_forces = new_base_f

        return base_energy, bias_energy, current_forces, pair_dists

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
            energy, _, current_forces, _ = self._md_step(
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

    def _log_starvation_diag(self, cycle, positions, cell, n_qualified) -> None:
        """Decompose the formation criterion on the LIVE groups (authoritative).

        Counts raw i-j geometric pairs in the formation window, then how many of
        those i (j) also have a valid constraint partner k (l) within its window,
        and how many satisfy both — vs the fully-qualified 4-tuple candidates.
        Isolates whether candidate scarcity is geometric (raw_ij ~ 0) or
        constraint-driven (raw_ij large but both ~ 0).
        """
        score_pairs = [p for p in self.template.pairs if p.score_pair]
        form = next((p for p in score_pairs if p.is_formation
                     and p.group_a in self.groups and p.group_b in self.groups), None)
        if form is None:
            return
        A = self.groups[form.group_a].atom_indices
        B = self.groups[form.group_b].atom_indices

        def _partner(anchor):
            for p in score_pairs:
                if p is form:
                    continue
                if p.group_a == anchor:
                    return p.group_b, p
                if p.group_b == anchor:
                    return p.group_a, p
            return None, None

        kc_label, kc = _partner(form.group_a)   # radical_C -> chain_C
        lc_label, lc = _partner(form.group_b)   # vinyl_alpha_C -> vinyl_beta_C
        Kset = self.groups[kc_label].atom_indices if kc_label else []
        Lset = self.groups[lc_label].atom_indices if lc_label else []

        def _has_partner(i, partners, ps):
            for q in partners:
                d = float(np.linalg.norm(minimum_image(positions[q] - positions[i], cell)))
                if ps.r_min <= d <= ps.r_max:
                    return True
            return False

        raw = vk = vl = both = 0
        for i in A:
            ik = _has_partner(i, Kset, kc) if kc else True
            for j in B:
                d = float(np.linalg.norm(minimum_image(positions[j] - positions[i], cell)))
                if not (form.r_min <= d <= form.r_max):
                    continue
                raw += 1
                jl = _has_partner(j, Lset, lc) if lc else True
                vk += ik
                vl += jl
                both += (ik and jl)
        logger.info(
            'DIAG-STARV cycle %d: raw_ij=%d (window [%.1f,%.1f]) '
            'with_validk=%d with_validl=%d both=%d qualified=%d',
            cycle, raw, form.r_min, form.r_max, vk, vl, both, n_qualified,
        )

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
        scored = score_candidates(candidates)

        # DIAG (candidate-starvation, 2026-06-24): decompose the formation
        # criterion live. Counts raw i-j geometric pairs (formation window only)
        # and how many of those also satisfy each constraint (i-k, j-l), vs the
        # fully-qualified 4-tuple candidates. Isolates geometry vs constraint.
        import os
        if os.environ.get('KAGOME_DIAG_STARVATION'):
            self._log_starvation_diag(cycle, state.positions, state.cell,
                                      len(candidates))
        if self._selection_log is not None:
            selected, decisions = audited_selection(scored)
            self._write_selection_audit(cycle, decisions)
        else:
            selected = select_non_overlapping(scored)

        selected = self._valence_filter(selected, state, cycle)

        active_pairs = self._build_pair_biases(selected, state.species)

        if self.bond_tracker:
            self.bond_tracker.record_attempts(
                active_pairs, state.positions, state.step, cycle, state.cell,
            )

        boost = BoostState()
        dt = self.config.timestep_fs
        last_bias_energy = 0.0

        box = validated_box(state.cell)

        base_energy, base_forces = self.calculator.compute(
            state.positions, state.species, state.cell,
        )
        bias_energy, bias_forces, _ = total_bias_fast(
            active_pairs, state.positions, boost, self.config.tdbb, box,
        )
        current_forces = base_forces + bias_forces

        form_indices = {i for i, p in enumerate(active_pairs) if p.is_formation}
        min_form_dist = float('inf')

        watchdog = StepWatchdog()
        steps_run = 0
        for step_in_phase in range(self.config.biased_steps):
            steps_run = step_in_phase + 1
            boost.advance(
                self.config.tdbb.gamma,
                self.config.tdbb.f1_max_formation,
                self.config.tdbb.f1_max_dissociation,
            )

            watchdog.arm()
            base_energy, bias_energy, current_forces, pair_dists = self._md_step(
                state, current_forces, dt, rng, step_in_phase,
                active_pairs=active_pairs, boost=boost, tdbb=self.config.tdbb,
                box=box,
            )
            watchdog.step_done(phase='biased', cycle=cycle, step=state.step)
            last_bias_energy = bias_energy

            if pair_dists is not None:
                for i in form_indices:
                    if pair_dists[i] < min_form_dist:
                        min_form_dist = pair_dists[i]

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
                        'Cycle %d biased: tentative reaction at step %d (%d pair(s)) '
                        '- ending biased phase for unbiased confirmation',
                        cycle, step_in_phase + 1, len(events),
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

    def _formation_pair_positions(self) -> tuple[int, int] | None:
        """Index positions (in template.groups) of the biased formation pair's
        two groups — the atoms whose bond would form.  None if not resolvable."""
        label_list = self.template.groups
        for ps in self.template.pairs:
            if ps.is_formation and not ps.constraint_only:
                try:
                    return label_list.index(ps.group_a), label_list.index(ps.group_b)
                except ValueError:
                    return None
        return None

    def _valence_filter(
        self, selected: list[Candidate], state: SimulationState, cycle: int,
    ) -> list[Candidate]:
        """Layer-2 occupancy guard: drop selected candidates whose formation
        would over-coordinate an atom in the current topology (specs/decisions.md
        2026-07-02).  Makes valence-safety a guaranteed invariant rather than an
        emergent property of group bookkeeping; a no-op when topology tracking is
        off or for non-vinyl chemistries.
        """
        if self._topology is None or not self._propagation_map:
            return selected
        pos = self._formation_pair_positions()
        if pos is None:
            return selected
        ia, ib = pos
        kept: list[Candidate] = []
        dropped = 0
        for cand in selected:
            radical_c = cand.atom_indices[ia]
            alpha_c = cand.atom_indices[ib]
            bad = vinyl_addition_over_coordinates(
                self._topology, radical_c, alpha_c,
                self._propagation_map, state.species,
            )
            if bad:
                dropped += 1
                logger.warning(
                    'Cycle %d valence guard: dropped candidate (radical %d + '
                    'alpha %d) — would over-coordinate atoms %s',
                    cycle, radical_c, alpha_c, bad,
                )
            else:
                kept.append(cand)
        if dropped:
            logger.info('Cycle %d valence guard: dropped %d/%d candidate(s)',
                        cycle, dropped, len(selected))
        return kept

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

        watchdog = StepWatchdog()
        for step_in_phase in range(self.config.unbiased_steps):
            watchdog.arm()
            energy, _, current_forces, _ = self._md_step(
                state, current_forces, dt, rng, step_in_phase,
            )
            watchdog.step_done(phase='unbiased', cycle=cycle, step=state.step)

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

    def _update_groups_after_cycle(
        self, state: SimulationState, cycle: int | None = None,
    ) -> None:
        self._updater.update(self.groups, self.bond_tracker, state)
        changed = self._apply_topology_updates(state)
        if changed and cycle is not None:
            self._write_topology_snapshot(state, cycle)

    def _apply_topology_updates(self, state: SimulationState) -> bool:
        """Apply reaction edits to the explicit bond topology for newly confirmed
        formations. Returns True if the topology changed this cycle.

        Vinyl radical addition (propagation_map present): new radical_C-alpha_C
        sigma bond + C=C opening + closed-shell placeholder-H shed
        (apply_vinyl_addition).  Other chemistries: add the recorded formation
        bond generically (condensation leaving-group edits are a follow-up).
        """
        if self._topology is None or self.bond_tracker is None:
            return False
        formations = self.bond_tracker.confirmed_formations()
        new = formations[self._topo_processed_formations:]
        for ev in new:
            if self._propagation_map:
                # Defensive: the selection valence guard should already prevent
                # over-coordinating formations, but never emit an over-valent
                # topology if one slips through — skip the edit and flag it.
                bad = vinyl_addition_over_coordinates(
                    self._topology, ev.atom_a, ev.atom_b,
                    self._propagation_map, state.species,
                )
                if bad:
                    logger.error(
                        'Confirmed formation (radical %d + alpha %d) would '
                        'over-coordinate %s; topology edit skipped (recorded '
                        'formation count is unaffected).',
                        ev.atom_a, ev.atom_b, bad,
                    )
                    continue
                apply_vinyl_addition(
                    self._topology, ev.atom_a, ev.atom_b,
                    self._propagation_map, state.species,
                )
            else:
                self._topology.add_bond(ev.atom_a, ev.atom_b, 1.0)
        self._topo_processed_formations = len(formations)
        return bool(new)

    def _write_topology_snapshot(self, state: SimulationState, cycle: int) -> None:
        """Append the current full bond list to topology.jsonl (one JSON line).

        Written only when connectivity changed, so the file is a compact,
        replayable history: line 0 is the initial topology (cycle -1), each later
        line the connectivity after a cycle that formed a bond.  The XYZ exporter
        replays it to attach the correct bonds to each frame.
        """
        if self._topology_log is None or self._topology is None:
            return
        record = {
            'step': state.step,
            'cycle': cycle,
            'n_bonds': len(self._topology),
            'bonds': [[i, j, o] for i, j, o in self._topology.bonds()],
        }
        with open(self._topology_log, 'a', encoding='utf-8') as f:
            f.write(json.dumps(record) + '\n')

    # Cap on rejected candidates written per cycle; large systems can enumerate
    # many. Selected candidates are always written in full.
    _MAX_REJECTED_LOGGED = 500

    def _write_selection_audit(self, cycle: int, decisions: list) -> None:
        """Append one JSON line recording the cycle's candidate ranking and the
        non-overlap selection/rejection decisions (RF18)."""
        if self._selection_log is None:
            return
        selected = [d for d in decisions if d.selected]
        rejected = [d for d in decisions if not d.selected]
        shown = rejected[:self._MAX_REJECTED_LOGGED]
        record = {
            'cycle': cycle,
            'n_candidates': len(decisions),
            'n_selected': len(selected),
            'n_rejected': len(rejected),
            'rejected_truncated': len(rejected) - len(shown),
            'selected': [
                {'atoms': list(d.atom_indices), 'score': d.score} for d in selected
            ],
            'rejected': [
                {'atoms': list(d.atom_indices), 'score': d.score, 'reason': d.reason}
                for d in shown
            ],
        }
        with open(self._selection_log, 'a', encoding='utf-8') as f:
            f.write(json.dumps(record) + '\n')
        if record['rejected_truncated']:
            logger.info(
                'Cycle %d selection audit: logged %d/%d rejected candidates '
                '(truncated %d)',
                cycle, len(shown), len(rejected), record['rejected_truncated'],
            )

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
        scored = score_candidates(candidates)
        selected = select_non_overlapping(scored)

        import os
        if os.environ.get('KAGOME_DIAG_STARVATION'):
            # Measure the ACTUAL bonded azo C-N distances (azo_C[k] is bonded to
            # azo_N[k]). If these are stretched beyond the activation window the
            # structure prep distorted the AIBN molecules -> activation can't fire.
            c_idx = activation_groups['azo_C'].atom_indices
            n_idx = activation_groups['azo_N'].atom_indices
            bonded = [
                float(np.linalg.norm(minimum_image(
                    state.positions[n] - state.positions[c], state.cell)))
                for c, n in zip(c_idx, n_idx)
            ]
            ps = activation_template.pairs[0]
            in_win = sum(ps.r_min <= d <= ps.r_max for d in bonded)
            logger.info(
                'DIAG-ACTIV: %d azo C-N bonds; dist min=%.2f max=%.2f mean=%.2f; '
                'in window [%.1f,%.1f]=%d; find_candidates=%d',
                len(bonded), min(bonded), max(bonded), sum(bonded) / len(bonded),
                ps.r_min, ps.r_max, in_win, len(candidates),
            )

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
            r = _pair_distance(state.positions, p.idx_a, p.idx_b, state.cell)
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

        box = validated_box(state.cell)

        base_energy, base_forces = self.calculator.compute(
            state.positions, state.species, state.cell,
        )
        bias_energy, bias_forces, _ = total_bias_fast(
            pairs, state.positions, boost, tdbb, box,
        )
        current_forces = base_forces + bias_forces

        # Absolute dissociation threshold for azo C-N bonds (decisions.md
        # 2026-06-19 RF8).  C-N equilibrium ~1.49 Å; 2.5 Å is well beyond
        # the barrier and indicates a clearly broken bond.  Expressed via
        # is_dissociated(r, r0=threshold, fraction=1.0) for API consistency
        # with BondTracker's r0-relative convention.
        dissoc_threshold = 2.5
        dissociated: list[tuple[int, int]] = []

        watchdog = StepWatchdog()
        for step_in_phase in range(activation_steps):
            boost.advance(tdbb.gamma, tdbb.f1_max_formation, tdbb.f1_max_dissociation)

            watchdog.arm()
            base_energy, bias_energy, current_forces, pair_dists = self._md_step(
                state, current_forces, dt, rng, step_in_phase,
                active_pairs=pairs, boost=boost, tdbb=tdbb,
                enable_barostat=False, box=box,
            )
            watchdog.step_done(phase='activation', cycle=-1, step=state.step)

            for i, p in enumerate(pairs):
                r = pair_dists[i] if pair_dists else _pair_distance(
                    state.positions, p.idx_a, p.idx_b, state.cell)
                if (step_in_phase + 1) % 500 == 0 or step_in_phase == 0:
                    logger.info(
                        'Activation step %d: (%d,%d) r=%.3f A, f1_d=%.1f, bias_E=%.1f',
                        step_in_phase + 1, p.idx_a, p.idx_b, r,
                        boost.f1_dissociation, bias_energy,
                    )
                if is_dissociated(r, dissoc_threshold) and (p.idx_a, p.idx_b) not in dissociated:
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
