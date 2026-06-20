"""Maxwell-Boltzmann velocity initialization for NVT/NPT MD.

Unit system (LAMMPS 'real'):
    mass: amu, velocity: Å/fs, energy: kcal/mol, temperature: K
    KE[kcal/mol] = 0.5 * m[amu] * v[Å/fs]^2 / FORCE_CONV
    Equipartition: 0.5 * m * v^2 / FORCE_CONV = 0.5 * KB * T (per DOF)
    => sigma = sqrt(KB * T * FORCE_CONV / m)  [Å/fs]
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from kagome.units import FORCE_CONV, KB


def maxwell_boltzmann_velocities(
    masses: NDArray[np.floating],
    temperature_K: float,
    rng: np.random.Generator,
    remove_com: bool = True,
) -> NDArray[np.floating]:
    """Assign Maxwell-Boltzmann velocities at target temperature.

    Args:
        masses:        per-atom masses (amu), shape (N,)
        temperature_K: target temperature (K)
        rng:           numpy random generator (for reproducibility)
        remove_com:    subtract center-of-mass drift (default True)

    Returns:
        velocities: shape (N, 3), units Å/fs
    """
    n = len(masses)
    sigmas = np.sqrt(KB * temperature_K * FORCE_CONV / masses)  # (N,)
    velocities = rng.standard_normal((n, 3)) * sigmas[:, np.newaxis]

    if remove_com and n > 1:
        total_mass = np.sum(masses)
        v_com = np.sum(masses[:, np.newaxis] * velocities, axis=0) / total_mass
        velocities -= v_com

    return velocities


def instant_temperature_K(
    velocities: NDArray[np.floating],
    masses: NDArray[np.floating] | None = None,
) -> float:
    """Kinetic temperature from velocities.

    T = 2*KE / (dof*KB),  KE[kcal/mol] = 0.5*sum(m*v^2)/FORCE_CONV.
    dof = 3N - 3 for N > 1 (center-of-mass translation removed; the MB
    initializer removes COM and NVE Verlet conserves it), else 3N. The COM
    correction is an O(1/N) diagnostic adjustment; see specs/decisions.md
    2026-06-20 RF19b.
    """
    n = velocities.shape[0]
    if n == 0:
        return 0.0
    if masses is not None:
        ke_amu = 0.5 * float(np.sum(masses[:, np.newaxis] * velocities ** 2))
    else:
        ke_amu = 0.5 * float(np.sum(velocities ** 2))
    ke_kcal = ke_amu / FORCE_CONV
    dof = 3 * n - 3 if n > 1 else 3 * n
    return 2.0 * ke_kcal / (dof * KB)
