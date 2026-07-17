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
from dataclasses import dataclass

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
    The remaining values mirror :class:`kagome.prep.openmm_equilibrate.
    ClassicalPrepConfig` precedent (no new engineering constants).
    """

    temperature_K: float
    n_steps: int
    timestep_fs: float = 0.5
    friction_per_ps: float = 1.0
    platform: str = 'CPU'                    # 'CUDA'|'OpenCL'|'CPU'|'Reference'
    minimize_tolerance_kj_mol_nm: float = 10.0

    def __post_init__(self) -> None:
        if self.n_steps <= 0:
            raise ValueError(f'n_steps must be positive; got {self.n_steps}')
        if self.temperature_K <= 0:
            raise ValueError(
                f'temperature_K must be positive; got {self.temperature_K}')
        if self.timestep_fs <= 0:
            raise ValueError(
                f'timestep_fs must be positive; got {self.timestep_fs}')


@dataclass
class MixMDResult:
    """Outcome of one classical mixing segment (OpenMM atom order, Å)."""

    positions_A: NDArray[np.floating]        # (M, 3) final coordinates
    minimized_energy_kj_mol: float           # after LocalEnergyMinimizer
    final_energy_kj_mol: float               # after the last MD chunk
    energies_kj_mol: list[float]             # per-chunk potential energies
    n_steps: int


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
                'dense box).'
            )

    out_pos_nm = np.array(
        context.getState(getPositions=True).getPositions()
        .value_in_unit(ommunit.nanometer),
        dtype=np.float64,
    )
    logger.info(
        'Classical mixing: %d steps done (E %.1f -> %.1f kJ/mol).',
        cfg.n_steps, minimized_e, energies[-1],
    )
    return MixMDResult(
        positions_A=out_pos_nm * ANGSTROM_PER_NM,
        minimized_energy_kj_mol=float(minimized_e),
        final_energy_kj_mol=energies[-1],
        energies_kj_mol=energies,
        n_steps=cfg.n_steps,
    )
