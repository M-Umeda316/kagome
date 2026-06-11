"""Langevin thermostat integrator for NVT dynamics.

Paper anchor: Section 2 describes NVT ensemble simulations.
Standard MD methodology — not a paper-specific equation.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

# Boltzmann constant in kcal/(mol·K)
KB_KCAL_MOL_K = 0.001987204


@dataclass
class LangevinParams:
    temperature_K: float = 300.0
    friction_per_fs: float = 0.01


class LangevinIntegrator:
    """BAOAB-style Langevin integrator with velocity randomization.

    c1 = exp(-γ·dt)
    c2 = sqrt(kT·(1 - c1²)/m)
    """

    def __init__(self, params: LangevinParams) -> None:
        self.params = params

    def step(
        self,
        positions: NDArray[np.floating],
        velocities: NDArray[np.floating],
        forces: NDArray[np.floating],
        masses: NDArray[np.floating] | None,
        dt: float,
        rng: np.random.Generator,
    ) -> None:
        gamma = self.params.friction_per_fs
        kT = KB_KCAL_MOL_K * self.params.temperature_K

        if masses is not None:
            inv_m = (1.0 / masses)[:, np.newaxis]
            m_col = masses[:, np.newaxis]
        else:
            inv_m = 1.0
            m_col = 1.0

        c1 = np.exp(-gamma * dt)
        c2 = np.sqrt(kT * (1.0 - c1 ** 2) / m_col) if masses is not None else np.sqrt(kT * (1.0 - c1 ** 2))

        accel = forces * inv_m

        # B: half-kick
        velocities += 0.5 * dt * accel
        # O: Ornstein-Uhlenbeck (thermostat)
        noise = rng.standard_normal(velocities.shape)
        velocities[:] = c1 * velocities + c2 * noise
        # A: drift
        positions += dt * velocities
        # B: half-kick (reusing same forces — single-force-eval variant)
        velocities += 0.5 * dt * accel
