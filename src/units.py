"""Physical constants and unit conversions for MD in 'real' units.

Unit system: energy=kcal/mol, distance=Å, time=fs, mass=amu(g/mol).
Matches LAMMPS 'real' style.
"""

# Boltzmann constant: kcal/(mol·K)
KB = 0.001987204

# F[kcal/(mol·Å)] / m[amu] → a[Å/fs²]
# Derived: 4184 / (N_A × 1e-10 × m_u) × 1e-20 ≈ 4.184e-4
FORCE_CONV = 4.184e-4

# Pressure: 1 atm in kcal/(mol·Å³)
# 1 kcal/(mol·Å³) = (4184/N_A) J / (1e-10)³ m³ ≈ 6.947e9 Pa = 6.947 GPa
# 1 atm = 101325 Pa → 101325 / 6.947e9 ≈ 1.4596e-5 kcal/(mol·Å³)
ATM_TO_KCAL_MOL_A3 = 1.4596e-5
