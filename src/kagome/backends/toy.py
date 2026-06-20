"""Toy calculator for testing the TDBB workflow without a real MLIP.

Returns a simple Lennard-Jones-like energy/force for smoke tests.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from kagome.backends.base import Calculator
from kagome.geometry import minimum_image


class ToyCalculator(Calculator):
    """Pairwise LJ-like potential: V = ε((σ/r)^12 - 2(σ/r)^6)."""

    def __init__(self, epsilon: float = 0.1, sigma: float = 1.5) -> None:
        self._epsilon = epsilon
        self._sigma = sigma

    @property
    def name(self) -> str:
        return 'toy'

    def compute(
        self,
        positions: NDArray[np.floating],
        species: list[str],
        cell: NDArray[np.floating] | None = None,
    ) -> tuple[float, NDArray[np.floating]]:
        n = positions.shape[0]
        energy = 0.0
        forces = np.zeros_like(positions, dtype=np.float64)
        eps = self._epsilon
        sig = self._sigma

        for i in range(n):
            for j in range(i + 1, n):
                r_vec = minimum_image(positions[j] - positions[i], cell)
                r = np.linalg.norm(r_vec)
                if r < 1e-12:
                    continue
                sr6 = (sig / r) ** 6
                sr12 = sr6 ** 2
                energy += eps * (sr12 - 2.0 * sr6)

                dv_dr = eps * (-12.0 * sr12 / r + 12.0 * sr6 / r)
                f_vec = dv_dr * (r_vec / r)
                forces[i] += f_vec
                forces[j] -= f_vec

        return energy, forces
