"""Langevin thermostat integrator (BAOAB splitting) for NVT dynamics.

Paper anchor: Section 2 describes NVT ensemble simulations.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from kagome.geometry import wrap_positions_fast
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
        self._cached_masses_id: int | None = None
        self._c1: float = 0.0
        self._c2: NDArray[np.floating] | float = 0.0
        self._inv_masses: NDArray[np.floating] | None = None
        self._box: NDArray[np.floating] | None = None

    def _update_cache(self, dt: float, masses: NDArray[np.floating] | None) -> None:
        masses_id = id(masses) if masses is not None else None
        if dt == self._cached_dt and masses_id == self._cached_masses_id:
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
        self._cached_masses_id = masses_id

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
