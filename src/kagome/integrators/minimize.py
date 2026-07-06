"""FIRE energy minimizer (numerical kernel).

Paper anchor: PDF p.20 (SI) — production reactive-acceleration MD is preceded by
equilibration.  Before any thermalization, the freshly placed (grid-packed)
structure carries close intermolecular contacts whose large forces would
otherwise convert to kinetic energy and spike the temperature.  Energy
minimization removes those contacts so the subsequent NPT equilibration can
thermalize to the setpoint.

FIRE (Fast Inertial Relaxation Engine), Bitzek et al., PRL 97, 170201 (2006).
Operates purely through the backend-agnostic ``Calculator.compute`` interface;
forces are in kcal/mol/Å, positions in Å.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Callable

import numpy as np
from numpy.typing import NDArray

from kagome.backends.base import Calculator
from kagome.geometry import wrap_positions

logger = logging.getLogger(__name__)


@dataclass
class FireParams:
    """FIRE minimizer parameters.

    fmax is the convergence threshold on the largest per-atom force magnitude.
    The default 1.0 kcal/mol/Å ≈ 0.043 eV/Å is slightly *stricter* than ASE's
    default of 0.05 eV/Å (≈ 1.15 kcal/mol/Å); combined with a max_steps cap it
    is sufficient for relaxing the close contacts of a freshly placed structure.
    """
    fmax_kcal_mol_A: float = 1.0
    max_steps: int = 500
    dt_start: float = 0.1
    dt_max: float = 1.0
    maxstep_A: float = 0.2          # per-step displacement clamp (Å) — robustness
    n_min: int = 5
    f_inc: float = 1.1
    f_dec: float = 0.5
    alpha_start: float = 0.1
    f_alpha: float = 0.99


@dataclass
class MinimizeResult:
    positions: NDArray[np.floating]
    energy: float
    fmax: float
    n_steps: int
    converged: bool


def fire_minimize(
    positions: NDArray[np.floating],
    species: list[str],
    cell: NDArray[np.floating] | None,
    calculator: Calculator,
    params: FireParams | None = None,
    on_step: Callable[[int, NDArray, float, float], None] | None = None,
) -> MinimizeResult:
    """Relax ``positions`` toward a local energy minimum via FIRE.

    Returns a MinimizeResult with the relaxed positions.  Unit masses are used
    (FIRE is a relaxation, not physical dynamics); a per-step displacement clamp
    (``maxstep_A``) keeps the very first steps stable even when the initial
    structure has severe clashes.

    on_step: optional callback(step, positions, energy, fmax) invoked after
    each FIRE iteration.
    """
    p = params or FireParams()
    pos = np.array(positions, dtype=np.float64, copy=True)
    # ASE sentinel: velocity is None until the first kick so the power/mixing
    # branch is skipped on step 0 (v is zero there).
    vel: NDArray[np.floating] | None = None

    dt = p.dt_start
    alpha = p.alpha_start
    n_positive = 0

    energy, forces = calculator.compute(pos, species, cell)
    fmax = float(np.sqrt((forces ** 2).sum(axis=1).max()))

    if on_step is not None:
        on_step(0, pos, energy, fmax)

    converged = fmax < p.fmax_kcal_mol_A
    step = 0
    while not converged and step < p.max_steps:
        # Canonical FIRE update order (Bitzek et al., PRL 97, 170201 (2006) /
        # ase.optimize.fire.FIRE.step): power test -> velocity mixing on the
        # *pre-kick* v and |v| -> dt/alpha update (check N_min, then count) ->
        # inertial kick.  On the first step v is zero, so — as in ASE — the
        # whole mixing/power branch is skipped and only the kick is applied.
        if vel is None:
            vel = np.zeros_like(pos)
        else:
            power = float(np.sum(forces * vel))
            if power > 0.0:
                fnorm = float(np.linalg.norm(forces))
                if fnorm > 1e-12:
                    vnorm = float(np.linalg.norm(vel))
                    vel = (1.0 - alpha) * vel + alpha * (forces / fnorm) * vnorm
                if n_positive > p.n_min:
                    dt = min(dt * p.f_inc, p.dt_max)
                    alpha *= p.f_alpha
                n_positive += 1
            else:
                n_positive = 0
                dt *= p.f_dec
                alpha = p.alpha_start
                vel[:] = 0.0

        # Inertial kick (semi-implicit Euler, unit mass) — after the mixing.
        vel += dt * forces

        dx = dt * vel
        dxmax = float(np.sqrt((dx ** 2).sum(axis=1).max()))
        if dxmax > p.maxstep_A:
            dx *= p.maxstep_A / dxmax
        pos += dx
        wrap_positions(pos, cell)

        energy, forces = calculator.compute(pos, species, cell)
        fmax = float(np.sqrt((forces ** 2).sum(axis=1).max()))
        converged = fmax < p.fmax_kcal_mol_A
        step += 1

        if on_step is not None:
            on_step(step, pos, energy, fmax)

    logger.info(
        'FIRE minimize: %d steps, fmax=%.3f kcal/mol/Å (target %.3f), '
        'E=%.2f kcal/mol, converged=%s',
        step, fmax, p.fmax_kcal_mol_A, energy, converged,
    )
    return MinimizeResult(
        positions=pos, energy=energy, fmax=fmax, n_steps=step, converged=converged,
    )


@dataclass
class CompressResult:
    positions: NDArray[np.floating]
    cell: NDArray[np.floating]
    n_stages: int


def compress_box(
    positions: NDArray[np.floating],
    cell: NDArray[np.floating],
    target_edge_A: float,
    species: list[str],
    calculator: Calculator,
    n_stages: int = 20,
    fire_params: FireParams | None = None,
) -> CompressResult:
    """Deterministically compress a cubic box to ``target_edge_A``.

    Used to reach paper-relevant densities at scale: the greedy placer can only
    seed molecules at a dilute density (large box), so the structure is then
    compressed in small geometric steps, FIRE-relaxing the close contacts that
    each shrink introduces.  This avoids relying on the slow MC barostat
    (max 1% volume change/move) to recover liquid density during equilibration.

    Paper anchor: density is the paper-specified initial condition (SI S-3); the
    compression path is a non-physical preparation device and does not bias the
    subsequent biased/unbiased dynamics.  Only compression is performed
    (target must be smaller than the current edge); otherwise it is a no-op.

    Assumes an isotropic (diagonal, cubic) cell; positions scale affinely with
    the cell at each stage.
    """
    p = fire_params or FireParams(fmax_kcal_mol_A=2.0, max_steps=200)
    pos = np.array(positions, dtype=np.float64, copy=True)
    cur_cell = np.array(cell, dtype=np.float64, copy=True)
    cur_edge = float(cur_cell[0, 0])

    if target_edge_A >= cur_edge:
        logger.info(
            'compress_box: target edge %.2f Å >= current %.2f Å — no compression.',
            target_edge_A, cur_edge,
        )
        return CompressResult(positions=pos, cell=cur_cell, n_stages=0)

    ratio = (target_edge_A / cur_edge) ** (1.0 / n_stages)  # per-stage linear factor
    logger.info(
        'compress_box: %.2f Å -> %.2f Å in %d stages (%.4f/stage)...',
        cur_edge, target_edge_A, n_stages, ratio,
    )
    for stage in range(1, n_stages + 1):
        pos *= ratio
        cur_cell *= ratio
        result = fire_minimize(pos, species, cur_cell, calculator, p)
        pos = result.positions
        logger.info(
            'compress_box stage %d/%d: edge=%.2f Å, fmax=%.2f kcal/mol/Å',
            stage, n_stages, float(cur_cell[0, 0]), result.fmax,
        )

    return CompressResult(positions=pos, cell=cur_cell, n_stages=n_stages)
