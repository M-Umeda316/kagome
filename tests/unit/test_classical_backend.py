"""Classical (OpenMM/OpenFF) Calculator backend.

Guarded with importorskip so the suite still runs in ML-only environments that
lack the OpenFF/OpenMM classical stack.
"""
from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip('openmm')
pytest.importorskip('openff.toolkit')

from scripts._systems import _INITIATOR_SMILES, _MONOMER_SMILES, build_vinyl_aibn_system
from kagome.backends.classical_backend import create_classical_calculator
from kagome.prep.openmm_equilibrate import MoleculeSpec


def _tiny_system():
    rng = np.random.default_rng(7)
    positions, species, _, _, _, _ = build_vinyl_aibn_system(
        n_monomers=2, n_initiators=1, box_size=25.0, rng=rng, rdkit_seed=42,
    )
    cell = np.diag([25.0, 25.0, 25.0])
    # MoleculeSpec order/seeds MUST match the builder: initiators first
    # (seed=42), then monomers (seed=43); see scripts/prep_structure.py.
    specs = [
        MoleculeSpec(_INITIATOR_SMILES, 1, rdkit_seed=42),
        MoleculeSpec(_MONOMER_SMILES, 2, rdkit_seed=43),
    ]
    return positions, species, cell, specs


def test_compute_returns_finite_energy_and_forces() -> None:
    positions, species, cell, specs = _tiny_system()
    calc = create_classical_calculator(specs, platform='CPU')

    energy, forces = calc.compute(positions, species, cell)

    assert np.isfinite(energy)
    assert forces.shape == (len(species), 3)
    assert np.all(np.isfinite(forces))
    assert calc.name.startswith('classical-')


def test_compress_box_shrinks_with_classical_calculator() -> None:
    from kagome.integrators.minimize import FireParams, compress_box

    positions, species, cell, specs = _tiny_system()
    calc = create_classical_calculator(specs, platform='CPU')

    target_edge = 18.0  # > 2×cutoff(0.8 nm=8 Å) so OpenMM accepts the dense box
    result = compress_box(
        positions, cell, target_edge, species, calc,
        n_stages=2, fire_params=FireParams(fmax_kcal_mol_A=5.0, max_steps=5),
    )

    assert result.cell[0, 0] < cell[0, 0]
    assert result.cell[0, 0] == pytest.approx(target_edge, abs=1e-6)
    assert result.positions.shape == positions.shape
    assert np.all(np.isfinite(result.positions))
