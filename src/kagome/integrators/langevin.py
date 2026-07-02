"""Langevin thermostat integrator (BAOAB splitting) for NVT dynamics.

Paper anchor: Section 2 describes NVT ensemble simulations.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from kagome.geometry import validated_box, wrap_positions_fast
from kagome.units import FORCE_CONV, KB, force_to_accel_fast, precompute_inv_masses


@dataclass
class LangevinParams:
    """Langevin サーモスタットのパラメータ (論文 SI: coupling 1.0 ps^-1)。"""

    temperature_K: float = 300.0
    friction_per_fs: float = 0.001


class LangevinIntegrator:
    """BAOAB Langevin integrator.

    Split: pre_force  = B-A-O-A  (half-kick, half-drift, thermostat, half-drift)
           post_force = B         (half-kick with new forces)
    """

    def __init__(self, params: LangevinParams) -> None:
        self.params = params
        self._cached_dt: float | None = None
        self._cached_masses: NDArray[np.floating] | None = None
        self._cached_temperature_K: float | None = None
        self._cached_friction: float | None = None
        self._c1: float = 0.0
        self._c2: NDArray[np.floating] | float = 0.0
        self._inv_masses: NDArray[np.floating] | None = None
        self._box: NDArray[np.floating] | None = None

    def _update_cache(self, dt: float, masses: NDArray[np.floating] | None) -> None:
        if (dt == self._cached_dt
                and masses is self._cached_masses
                and self.params.temperature_K == self._cached_temperature_K
                and self.params.friction_per_fs == self._cached_friction):
            return
        gamma = self.params.friction_per_fs
        kT = KB * self.params.temperature_K
        self._c1 = float(np.exp(-gamma * dt))
        if masses is not None:
            m_col = masses[:, np.newaxis]
            self._c2 = np.sqrt(kT * FORCE_CONV * (1.0 - self._c1 ** 2) / m_col)
        else:
            self._c2 = np.sqrt(kT * FORCE_CONV * (1.0 - self._c1 ** 2))
        self._inv_masses = precompute_inv_masses(masses)
        self._cached_dt = dt
        self._cached_masses = masses
        self._cached_temperature_K = self.params.temperature_K
        self._cached_friction = self.params.friction_per_fs

    def _get_box(self, cell: NDArray[np.floating] | None) -> NDArray[np.floating] | None:
        if cell is None:
            return None
        box = validated_box(cell)
        if box is None:
            return None
        if self._box is None or self._box[0] != box[0] or self._box[1] != box[1] or self._box[2] != box[2]:
            self._box = box
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
        self._update_cache(dt, masses)
        accel = force_to_accel_fast(forces, self._inv_masses)
        c1 = self._c1
        c2 = self._c2

        # B: half-kick
        velocities += 0.5 * dt * accel
        # A: half-drift
        positions += 0.5 * dt * velocities
        # O: Ornstein-Uhlenbeck thermostat
        noise = rng.standard_normal(velocities.shape)
        velocities[:] = c1 * velocities + c2 * noise
        # A: half-drift
        positions += 0.5 * dt * velocities
        wrap_positions_fast(positions, self._get_box(cell))

    def post_force(
        self,
        velocities: NDArray[np.floating],
        forces: NDArray[np.floating],
        masses: NDArray[np.floating] | None,
        dt: float,
    ) -> None:
        self._update_cache(dt, masses)
        accel = force_to_accel_fast(forces, self._inv_masses)
        velocities += 0.5 * dt * accel
