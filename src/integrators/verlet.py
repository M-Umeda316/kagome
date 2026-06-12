"""Velocity Verlet integrator with correct unit conversion."""
from __future__ import annotations

from typing import Protocol

import numpy as np
from numpy.typing import NDArray

from src.units import FORCE_CONV


class Integrator(Protocol):
    def pre_force(
        self,
        positions: NDArray[np.floating],
        velocities: NDArray[np.floating],
        forces: NDArray[np.floating],
        masses: NDArray[np.floating] | None,
        dt: float,
        rng: np.random.Generator,
    ) -> None:
        """Half-kick + drift.  Call BEFORE force evaluation at new positions."""
        ...

    def post_force(
        self,
        velocities: NDArray[np.floating],
        forces: NDArray[np.floating],
        masses: NDArray[np.floating] | None,
        dt: float,
    ) -> None:
        """Second half-kick with newly computed forces."""
        ...


class VelocityVerletIntegrator:
    """Standard velocity Verlet (leapfrog split).

    pre_force:  v += 0.5·dt·a(old),  x += dt·v
    post_force: v += 0.5·dt·a(new)
    """

    def pre_force(
        self,
        positions: NDArray[np.floating],
        velocities: NDArray[np.floating],
        forces: NDArray[np.floating],
        masses: NDArray[np.floating] | None,
        dt: float,
        rng: np.random.Generator,
    ) -> None:
        accel = self._accel(forces, masses)
        velocities += 0.5 * dt * accel
        positions += dt * velocities

    def post_force(
        self,
        velocities: NDArray[np.floating],
        forces: NDArray[np.floating],
        masses: NDArray[np.floating] | None,
        dt: float,
    ) -> None:
        accel = self._accel(forces, masses)
        velocities += 0.5 * dt * accel

    @staticmethod
    def _accel(
        forces: NDArray[np.floating],
        masses: NDArray[np.floating] | None,
    ) -> NDArray[np.floating]:
        if masses is not None:
            inv_m = (1.0 / masses)[:, np.newaxis]
        else:
            inv_m = 1.0
        return FORCE_CONV * forces * inv_m
