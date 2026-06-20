"""Lock the physical constants in src.units to first-principles derivations.

These guard against transcription errors like the RF19 ATM_TO_KCAL_MOL_A3 fix.
All derivations use CODATA 2018 values; the unit system is LAMMPS 'real'
(energy=kcal/mol, distance=Å, time=fs, mass=amu).
"""
from __future__ import annotations

import numpy as np

from src.units import (
    ATM_TO_KCAL_MOL_A3,
    EV_TO_KCAL_MOL,
    FORCE_CONV,
    KB,
    force_to_accel,
)

# CODATA 2018
_N_A = 6.02214076e23          # mol^-1
_KB_J = 1.380649e-23          # J/K
_EV_J = 1.602176634e-19       # J
_AMU_KG = 1.66053906660e-27   # kg
_CAL = 4184.0                 # J per kcal


def test_boltzmann_kcal_mol_K():
    # kB·N_A / (J per kcal) = kcal/(mol·K)
    expected = _KB_J * _N_A / _CAL
    np.testing.assert_allclose(KB, expected, rtol=1e-5)


def test_force_conv_kcal_mol_A_amu_to_A_fs2():
    # a[Å/fs²] = FORCE_CONV · F[kcal/(mol·Å)] / m[amu]
    # Derive FORCE_CONV from first principles:
    #   F = 1 kcal/(mol·Å) = (_CAL/_N_A) J / 1e-10 m  [N]
    #   a = F / (1 amu)  [m/s²]  → convert to Å/fs² (×1e-20)
    f_newton = (_CAL / _N_A) / 1e-10
    a_m_s2 = f_newton / _AMU_KG
    a_A_fs2 = a_m_s2 * 1e-20
    np.testing.assert_allclose(FORCE_CONV, a_A_fs2, rtol=1e-3)


def test_atm_to_kcal_mol_A3():
    # 1 kcal/(mol·Å³) in Pa, then 1 atm in those units (RF19 transcription guard)
    kcal_mol_A3_in_Pa = (_CAL / _N_A) / 1e-30
    expected = 101325.0 / kcal_mol_A3_in_Pa
    np.testing.assert_allclose(ATM_TO_KCAL_MOL_A3, expected, rtol=1e-4)
    # The old wrong literal (1.4596e-5) must NOT pass this guard.
    assert not np.isclose(1.4596e-5, expected, rtol=1e-4)


def test_ev_to_kcal_mol():
    expected = _EV_J * _N_A / _CAL
    np.testing.assert_allclose(EV_TO_KCAL_MOL, expected, rtol=1e-5)


def test_force_to_accel_uses_force_conv():
    forces = np.array([[2.0, 0.0, 0.0]])
    masses = np.array([4.0])
    accel = force_to_accel(forces, masses)
    np.testing.assert_allclose(accel[0, 0], FORCE_CONV * 2.0 / 4.0, rtol=1e-12)


def test_force_to_accel_unit_mass_when_none():
    forces = np.array([[3.0, 0.0, 0.0]])
    accel = force_to_accel(forces, None)
    np.testing.assert_allclose(accel[0, 0], FORCE_CONV * 3.0, rtol=1e-12)
