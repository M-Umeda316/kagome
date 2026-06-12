"""Physical constants and unit conversions for MD in 'real' units.

Unit system: energy=kcal/mol, distance=Å, time=fs, mass=amu(g/mol).
Matches LAMMPS 'real' style.
"""

# Boltzmann constant: kcal/(mol·K)
KB = 0.001987204

# F[kcal/(mol·Å)] / m[amu] → a[Å/fs²]
# Derived: 4184 / (N_A × 1e-10 × m_u) × 1e-20 ≈ 4.184e-4
FORCE_CONV = 4.184e-4
