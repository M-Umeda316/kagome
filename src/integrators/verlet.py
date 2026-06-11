"""Velocity Verlet integrator."""
from __future__ import annotations

from typing import Protocol

import numpy as np
from numpy.typing import NDArray


class Integrator(Protocol):
    def step(
        self,
        positions: NDArray[np.floating],
        velocities: NDArray[np.floating],
        forces: NDArray[np.floating],
        masses: NDArray[np.floating] | None,
        dt: float,
        rng: np.random.Generator,
    ) -> None:
        """Integrate one step in-place."""
        ...


class VelocityVerletIntegrator:
    """Standard velocity Verlet. Uses masses when provided, unit masses otherwise."""

    def step(
        self,
        positions: NDArray[np.floating],
        velocities: NDArray[np.floating],
        forces: NDArray[np.floating],
        masses: NDArray[np.floating] | None,
        dt: float,
        rng: np.random.Generator,
    ) -> None:
        if masses is not None:
            inv_m = (1.0 / masses)[:, np.newaxis]
        else:
            inv_m = 1.0

        accel = forces * inv_m
        velocities += 0.5 * dt * accel
        positions += dt * velocities
        velocities += 0.5 * dt * accel
