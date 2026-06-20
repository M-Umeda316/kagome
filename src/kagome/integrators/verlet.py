"""Velocity Verlet integrator with correct unit conversion."""
from __future__ import annotations

from typing import Protocol

import numpy as np
from numpy.typing import NDArray

from kagome.geometry import wrap_positions
from kagome.units import force_to_accel


class Integrator(Protocol):
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

    pre_force:  v += 0.5·dt·a(old),  x += dt·v,  wrap(x)
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
        cell: NDArray[np.floating] | None = None,
    ) -> None:
        accel = force_to_accel(forces, masses)
        velocities += 0.5 * dt * accel
        positions += dt * velocities
        wrap_positions(positions, cell)

    def post_force(
        self,
        velocities: NDArray[np.floating],
        forces: NDArray[np.floating],
        masses: NDArray[np.floating] | None,
        dt: float,
    ) -> None:
        accel = force_to_accel(forces, masses)
        velocities += 0.5 * dt * accel
