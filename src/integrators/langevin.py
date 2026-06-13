"""Langevin thermostat integrator (BAOAB splitting) for NVT dynamics.

Paper anchor: Section 2 describes NVT ensemble simulations.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from src.geometry import wrap_positions
from src.units import FORCE_CONV, KB


@dataclass
class LangevinParams:
    temperature_K: float = 300.0
    friction_per_fs: float = 0.01


class LangevinIntegrator:
    """BAOAB Langevin integrator.

    Split: pre_force  = B-A-O-A  (half-kick, half-drift, thermostat, half-drift)
           post_force = B         (half-kick with new forces)
    """

    def __init__(self, params: LangevinParams) -> None:
        self.params = params

    def pre_force(
        self,
        positions: NDArray[np.floating],
        velocities: NDArray[np.floating],
        forces: NDArray[np.floating],
        masses: NDArray[np.floating] | None,
        dt: float,
        rng: np.random.Generator,
        cell: NDArray[np.floating] | None = None,
    ) -> None:
        gamma = self.params.friction_per_fs
        kT = KB * self.params.temperature_K

        accel = _accel(forces, masses)
        c1 = np.exp(-gamma * dt)

        if masses is not None:
            m_col = masses[:, np.newaxis]
            c2 = np.sqrt(kT * FORCE_CONV * (1.0 - c1 ** 2) / m_col)
        else:
            c2 = np.sqrt(kT * FORCE_CONV * (1.0 - c1 ** 2))

        # B: half-kick
        velocities += 0.5 * dt * accel
        # A: half-drift
        positions += 0.5 * dt * velocities
        # O: Ornstein-Uhlenbeck thermostat
        noise = rng.standard_normal(velocities.shape)
        velocities[:] = c1 * velocities + c2 * noise
        # A: half-drift
        positions += 0.5 * dt * velocities
        wrap_positions(positions, cell)

    def post_force(
        self,
        velocities: NDArray[np.floating],
        forces: NDArray[np.floating],
        masses: NDArray[np.floating] | None,
        dt: float,
    ) -> None:
        accel = _accel(forces, masses)
        velocities += 0.5 * dt * accel


def _accel(
    forces: NDArray[np.floating],
    masses: NDArray[np.floating] | None,
) -> NDArray[np.floating]:
    if masses is not None:
        inv_m = (1.0 / masses)[:, np.newaxis]
    else:
        inv_m = 1.0
    return FORCE_CONV * forces * inv_m
