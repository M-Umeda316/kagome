"""Velocity Verlet integrator with correct unit conversion."""
from __future__ import annotations

from typing import Protocol

import numpy as np
from numpy.typing import NDArray

from kagome.geometry import wrap_positions_fast
from kagome.units import force_to_accel_fast, precompute_inv_masses


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

    def __init__(self) -> None:
        self._inv_masses: NDArray[np.floating] | None = None
        self._masses_id: int | None = None
        self._box: NDArray[np.floating] | None = None

    def _get_inv_masses(self, masses: NDArray[np.floating] | None) -> NDArray[np.floating] | None:
        mid = id(masses) if masses is not None else None
        if mid != self._masses_id:
            self._inv_masses = precompute_inv_masses(masses)
            self._masses_id = mid
        return self._inv_masses

    def _get_box(self, cell: NDArray[np.floating] | None) -> NDArray[np.floating] | None:
        if cell is None:
            return None
        d0, d1, d2 = cell[0, 0], cell[1, 1], cell[2, 2]
        if self._box is None or self._box[0] != d0 or self._box[1] != d1 or self._box[2] != d2:
            self._box = np.array([d0, d1, d2])
        return self._box

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
        accel = force_to_accel_fast(forces, self._get_inv_masses(masses))
        velocities += 0.5 * dt * accel
        positions += dt * velocities
        wrap_positions_fast(positions, self._get_box(cell))

    def post_force(
        self,
        velocities: NDArray[np.floating],
        forces: NDArray[np.floating],
        masses: NDArray[np.floating] | None,
        dt: float,
    ) -> None:
        accel = force_to_accel_fast(forces, self._get_inv_masses(masses))
        velocities += 0.5 * dt * accel
