"""Velocity Verlet integrator with correct unit conversion."""
from __future__ import annotations

from typing import Protocol

import numpy as np
from numpy.typing import NDArray

from kagome.geometry import validated_box_cached, wrap_positions_fast
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
        self._cached_masses: NDArray[np.floating] | None = None
        self._box: NDArray[np.floating] | None = None
        self._cached_cell: NDArray[np.floating] | None = None

    def _get_inv_masses(self, masses: NDArray[np.floating] | None) -> NDArray[np.floating] | None:
        if masses is not self._cached_masses:
            self._inv_masses = precompute_inv_masses(masses)
            self._cached_masses = masses
        return self._inv_masses

    def _get_box(self, cell: NDArray[np.floating] | None) -> NDArray[np.floating] | None:
        # Validate only on a cell content change; unchanged steps reuse the
        # cached box (see geometry.validated_box_cached).
        self._box, self._cached_cell = validated_box_cached(
            cell, self._cached_cell, self._box,
        )
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
