"""Classical NVT mixing-MD kernel for the well-mixed measurement mode (WM-P3).

Part of the "well-mixed 測定モード" (specs/decisions.md 2026-07-17, 追補
2026-07-18 (l)). The workflow's mixing phase builds a :class:`ClassicalMix`
via :func:`kagome.prep.mixing.build_classical_mix` and hands it here;
:func:`run_mix_md` honours the translator's ``requires_minimization=True``
contract (LocalEnergyMinimizer before any dynamics), seeds Maxwell-Boltzmann
velocities OpenMM-side, then runs Langevin NVT in chunks with an early
non-finite-energy bail-out (same pattern as
:mod:`kagome.prep.openmm_equilibrate`).

This module is the numerical kernel only — it never touches workflow state.
The caller draws ``seed`` from the workflow rng (determinism, decisions.md
追補 2026-07-18 (i)) and maps the returned OpenMM-order coordinates back to
the MLIP order via ``ClassicalMix.write_back``.

All OpenMM imports are deferred to call time so this module imports cleanly
in ML-only environments.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from kagome.prep.mixing import ClassicalMix
from kagome.units import ANGSTROM_PER_NM, NM_PER_ANGSTROM

logger = logging.getLogger(__name__)


@dataclass
class MixMDConfig:
    """Parameters for one classical mixing-MD segment.

    ``n_steps`` (from the workflow's ``mix_time_ps``) carries no default —
    the mixing duration must be chosen explicitly until the P4 sweep
    establishes measured defaults (specs/decisions.md 2026-07-17 (v)).
    ``timestep_fs`` / ``friction_per_ps`` / ``platform`` mirror
    :class:`kagome.prep.openmm_equilibrate.ClassicalPrepConfig` precedent.

    Soft-start warm-up (WM-P4 robustness, specs/decisions.md 追補 2026-07-18):
    the classical mixing PES (Sage LJ) has a much harder repulsive wall than the
    MLIP (orb) that produced the handed-over coordinates, so the translated
    system intermittently carries a close contact — an orb-equilibrated atom pair
    (or an injected cap/placeholder H) closer than Sage's LJ core tolerates. At
    production scale (~500+ atoms, 0.5 g/mL) the ``LocalEnergyMinimizer``'s
    RMS-force tolerance is met globally while such a contact stays hot (its force
    is diluted across all atoms), and stepping the configured ``timestep_fs``
    Langevin straight away diverges to NaN. Before the main mixing we therefore
    run a short warm-up MD segment at a much smaller timestep with high friction,
    ramping up to ``timestep_fs``, which bleeds off the residual hot contacts —
    including contacts between Sage-*constrained* X–H atoms, which a geometric
    push-apart cannot fix (the constraint snaps a displaced H back) but bounded,
    overdamped dynamics can (it moves whole rigid groups apart). The warm-up steps
    are counted and logged SEPARATELY from the reported mixing ``n_steps`` (they
    are not part of the mixing measurement). These are numerical-stability knobs,
    so — unlike the mixing *duration* (decisions.md (v)) — they carry documented
    engineering defaults. The bridge is needed EVERY mixing cycle (the wall
    mismatch is reaction-independent, specs/decisions.md 追補 2026-07-18 (WM-P4
    bridge)), not only reacted ones, so it runs whenever ``warmup_steps > 0``.
    Set ``warmup_steps=0`` to disable it.

    ``minimize_tolerance_kj_mol_nm`` is the RMS-force tolerance for the
    pre-dynamics ``LocalEnergyMinimizer`` (the ``requires_minimization``
    contract); ``warmup_minimize_tolerance_kj_mol_nm`` (< the above) drives a
    tighter second minimization pass — run every cycle (when ``warmup_steps > 0``)
    — squeezing the worst contacts before the warm-up dynamics.
    """

    temperature_K: float
    n_steps: int
    timestep_fs: float = 0.5
    friction_per_ps: float = 1.0
    platform: str = 'CPU'                    # 'CUDA'|'OpenCL'|'CPU'|'Reference'
    minimize_tolerance_kj_mol_nm: float = 10.0
    # ── soft-start warm-up (see class docstring) ──
    warmup_steps: int = 2000                  # 0 disables the warm-up
    # Start the ramp very gently (0.05 fs): the close contact can be between ANY
    # handed-over pair (not just a moderate injected one), so the first stage must
    # be a near-steepest-descent crawl to survive a hard orb->Sage contact.
    warmup_timestep_fs: float = 0.05          # starting (smallest) ramp timestep
    warmup_friction_per_ps: float = 50.0      # overdamp the hot contacts
    warmup_stages: int = 6                    # timestep-ramp stages -> timestep_fs
    warmup_minimize_tolerance_kj_mol_nm: float = 1.0  # tighter 2nd-pass min tol

    def __post_init__(self) -> None:
        if self.n_steps <= 0:
            raise ValueError(f'n_steps must be positive; got {self.n_steps}')
        if self.temperature_K <= 0:
            raise ValueError(
                f'temperature_K must be positive; got {self.temperature_K}')
        if self.timestep_fs <= 0:
            raise ValueError(
                f'timestep_fs must be positive; got {self.timestep_fs}')
        if self.warmup_steps < 0:
            raise ValueError(
                f'warmup_steps must be >= 0; got {self.warmup_steps}')
        if self.warmup_steps > 0:
            if self.warmup_timestep_fs <= 0:
                raise ValueError(
                    'warmup_timestep_fs must be positive; got '
                    f'{self.warmup_timestep_fs}')
            if self.warmup_timestep_fs > self.timestep_fs:
                raise ValueError(
                    'warmup_timestep_fs must be <= timestep_fs (the warm-up '
                    f'ramps UP to the mixing timestep); got '
                    f'{self.warmup_timestep_fs} > {self.timestep_fs}')
            if self.warmup_friction_per_ps < 0:
                raise ValueError(
                    'warmup_friction_per_ps must be >= 0; got '
                    f'{self.warmup_friction_per_ps}')
            if self.warmup_stages < 1:
                raise ValueError(
                    f'warmup_stages must be >= 1; got {self.warmup_stages}')


@dataclass
class MixMDResult:
    """Outcome of one classical mixing segment (OpenMM atom order, Å).

    ``n_steps`` / ``energies_kj_mol`` / ``final_energy_kj_mol`` describe the
    reported *mixing* segment only. The soft-start warm-up (WM-P4) is tracked
    separately in ``n_warmup_steps`` / ``warmup_energies_kj_mol`` so it never
    contaminates the mixing measurement.
    """

    positions_A: NDArray[np.floating]        # (M, 3) final coordinates
    minimized_energy_kj_mol: float           # after LocalEnergyMinimizer
    final_energy_kj_mol: float               # after the last MD chunk
    energies_kj_mol: list[float]             # per-chunk potential energies
    n_steps: int
    n_warmup_steps: int = 0                   # soft-start steps (outside n_steps)
    warmup_energies_kj_mol: list[float] = field(default_factory=list)


def run_mix_md(
    mix: ClassicalMix,
    cfg: MixMDConfig,
    seed: int,
) -> MixMDResult:
    """Minimize + Langevin-NVT-mix a translated classical system.

    ``seed`` seeds both the OpenMM Langevin integrator and the
    Maxwell-Boltzmann velocity draw; the caller derives it from the workflow
    rng so the whole run stays reproducible from one seed.

    Raises RuntimeError when the potential energy turns non-finite mid-run
    (diverged classical MD — e.g. timestep too large for the dense box).
    """
    import openmm
    from openmm import unit as ommunit

    integrator = openmm.LangevinMiddleIntegrator(
        cfg.temperature_K * ommunit.kelvin,
        cfg.friction_per_ps / ommunit.picosecond,
        cfg.timestep_fs * ommunit.femtosecond,
    )
    integrator.setRandomNumberSeed(int(seed))
    platform = openmm.Platform.getPlatformByName(cfg.platform)
    context = openmm.Context(mix.system, integrator, platform)

    box_nm = np.asarray(mix.box_vectors_A, dtype=np.float64) * NM_PER_ANGSTROM
    context.setPeriodicBoxVectors(*[openmm.Vec3(*row) * ommunit.nanometer
                                    for row in box_nm])
    context.setPositions(
        (mix.positions_A * NM_PER_ANGSTROM) * ommunit.nanometer)

    # requires_minimization contract (ClassicalMix.metadata): cap H sit at ideal
    # guesses and placeholder H may clash with their former parent — always
    # minimize before dynamics.
    openmm.LocalEnergyMinimizer.minimize(
        context, cfg.minimize_tolerance_kj_mol_nm, 0,
    )
    # Tighter second pass — run EVERY mixing cycle, not only reacted ones. The
    # orb<->Sage repulsive-wall mismatch that produces a Sage-intolerable close
    # contact is reaction-independent (specs/decisions.md 追補 2026-07-18 (WM-P4
    # bridge)): a cycle with NO injected atoms (production cycle 2) can still hand
    # over a hard contact, so the classical<->MLIP bridge (tighter minimize +
    # warm-up) is needed regardless of injection. Cheap (the first pass already
    # did the bulk) and never loosens the tolerance.
    if (cfg.warmup_steps > 0
            and cfg.warmup_minimize_tolerance_kj_mol_nm
            < cfg.minimize_tolerance_kj_mol_nm):
        openmm.LocalEnergyMinimizer.minimize(
            context, cfg.warmup_minimize_tolerance_kj_mol_nm, 0,
        )
    minimized_e = context.getState(getEnergy=True).getPotentialEnergy() \
        .value_in_unit(ommunit.kilojoule_per_mole)
    if not np.isfinite(minimized_e):
        raise RuntimeError(
            'Classical mixing: non-finite energy after minimization — the '
            'translated system is unphysical (check the bond graph handed to '
            'build_classical_mix).'
        )
    context.setVelocitiesToTemperature(
        cfg.temperature_K * ommunit.kelvin, int(seed))

    # ── soft-start warm-up (classical<->MLIP bridge; every cycle) ──
    # Bleeds off any residual hot contact left after minimization, whether from an
    # injected atom or an orb-handed-over pair inside Sage's repulsive wall.
    n_warmup = 0
    warmup_energies: list[float] = []
    if cfg.warmup_steps > 0:
        n_warmup, warmup_energies = _run_warmup(
            context, integrator, cfg, ommunit, np,
        )
        # Restore the configured mixing timestep / friction for the measured
        # segment (the warm-up left them at their ramped / high-friction values).
        integrator.setStepSize(cfg.timestep_fs * ommunit.femtosecond)
        integrator.setFriction(cfg.friction_per_ps / ommunit.picosecond)

    chunk = max(1, cfg.n_steps // 20)
    done = 0
    energies: list[float] = []
    while done < cfg.n_steps:
        take = min(chunk, cfg.n_steps - done)
        integrator.step(take)
        done += take
        e = context.getState(getEnergy=True).getPotentialEnergy() \
            .value_in_unit(ommunit.kilojoule_per_mole)
        energies.append(float(e))
        if not np.isfinite(e):
            raise RuntimeError(
                f'Classical mixing: NVT became non-finite after {done} steps '
                f'(timestep {cfg.timestep_fs} fs may be too large for the '
                'dense box even after the soft-start warm-up).'
            )

    out_pos_nm = np.array(
        context.getState(getPositions=True).getPositions()
        .value_in_unit(ommunit.nanometer),
        dtype=np.float64,
    )
    logger.info(
        'Classical mixing: %d steps done (warm-up %d; E %.1f -> %.1f kJ/mol).',
        cfg.n_steps, n_warmup, minimized_e, energies[-1],
    )
    return MixMDResult(
        positions_A=out_pos_nm * ANGSTROM_PER_NM,
        minimized_energy_kj_mol=float(minimized_e),
        final_energy_kj_mol=energies[-1],
        energies_kj_mol=energies,
        n_steps=cfg.n_steps,
        n_warmup_steps=n_warmup,
        warmup_energies_kj_mol=warmup_energies,
    )


def _run_warmup(context, integrator, cfg, ommunit, np) -> tuple[int, list[float]]:
    """Soft-start MD before the main mixing (WM-P4; see :class:`MixMDConfig`).

    Ramps the integrator timestep geometrically from ``warmup_timestep_fs`` up
    to the configured ``timestep_fs`` over ``warmup_stages`` stages, all at the
    high ``warmup_friction_per_ps`` so the injected hot contacts are overdamped
    rather than kicked. Deterministic (the integrator's seed is already set) and
    stepped in small chunks with the same non-finite bail-out as the main loop,
    so a genuinely unrecoverable start still fails informatively instead of
    silently producing NaN. Returns ``(steps_run, per-chunk energies)``.
    """
    integrator.setFriction(cfg.warmup_friction_per_ps / ommunit.picosecond)
    stages = cfg.warmup_stages
    base = cfg.warmup_steps // stages
    remainder = cfg.warmup_steps - base * stages
    dt0, dt1 = cfg.warmup_timestep_fs, cfg.timestep_fs

    done = 0
    energies: list[float] = []
    for k in range(stages):
        frac = k / (stages - 1) if stages > 1 else 1.0
        dt = dt0 * (dt1 / dt0) ** frac          # geometric ramp dt0 -> dt1
        integrator.setStepSize(dt * ommunit.femtosecond)
        stage_steps = base + (1 if k < remainder else 0)
        # Check finiteness a few times per stage so a divergence is caught early.
        sub = max(1, stage_steps // 4)
        left = stage_steps
        while left > 0:
            take = min(sub, left)
            integrator.step(take)
            left -= take
            done += take
            e = context.getState(getEnergy=True).getPotentialEnergy() \
                .value_in_unit(ommunit.kilojoule_per_mole)
            energies.append(float(e))
            if not np.isfinite(e):
                raise RuntimeError(
                    'Classical mixing: soft-start warm-up became non-finite '
                    f'after {done} warm-up steps (dt {dt:.3f} fs) — the '
                    'translated system is too hot even for the warm-up (a '
                    'handed-over pair or injected H is inside Sage\'s repulsive '
                    'wall and could not be bridged).'
                )
    return done, energies
