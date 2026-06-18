"""Physical constants and unit conversions for MD in 'real' units.

Unit system: energy=kcal/mol, distance=Å, time=fs, mass=amu(g/mol).
Matches LAMMPS 'real' style.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

# Boltzmann constant: kcal/(mol·K)
KB = 0.001987204

# F[kcal/(mol·Å)] / m[amu] → a[Å/fs²]
# Derived: 4184 / (N_A × 1e-10 × m_u) × 1e-20 ≈ 4.184e-4
FORCE_CONV = 4.184e-4

# Pressure: 1 atm in kcal/(mol·Å³)
# 1 kcal/(mol·Å³) = (4184/N_A) J / (1e-10)³ m³ ≈ 6.947e9 Pa = 6.947 GPa
# 1 atm = 101325 Pa → 101325 / 6.947e9 ≈ 1.4596e-5 kcal/(mol·Å³)
ATM_TO_KCAL_MOL_A3 = 1.4596e-5

# 1 eV = 23.060548 kcal/mol (NIST CODATA 2018)
EV_TO_KCAL_MOL = 23.060548

# ── Boundary conversions for OpenMM (nm, kJ/mol) ↔ this repo (Å, kcal/mol) ─────
# Used only at the classical structure-prep boundary (src/prep). OpenMM's native
# unit system is nm / kJ/mol / ps; the rest of the repo uses Å / kcal/mol / fs.
ANGSTROM_PER_NM = 10.0
NM_PER_ANGSTROM = 0.1
KJ_PER_KCAL = 4.184
KCAL_PER_KJ = 1.0 / 4.184


def force_to_accel(
    forces: NDArray[np.floating],
    masses: NDArray[np.floating] | None,
) -> NDArray[np.floating]:
    """Convert forces [kcal/(mol·Å)] to accelerations [Å/fs²]."""
    if masses is not None:
        inv_m = (1.0 / masses)[:, np.newaxis]
    else:
        inv_m = 1.0
    return FORCE_CONV * forces * inv_m
