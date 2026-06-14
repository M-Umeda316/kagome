"""Unit tests for the classical structure-prep scaffolding (Phase 1).

Covers the env-agnostic pieces: atom-order sharing between the system builder
and prep (decision D-4), the portable structure handoff file, unit-conversion
constants, and the prep config defaults. The OpenMM/OpenFF body is exercised in
the Phase 2/3 integration tests (it requires those packages to be installed).
"""
from __future__ import annotations

import numpy as np
import pytest

from scripts._systems import _MONOMER_SMILES, _rdkit_3d, _rdkit_mol
from src.prep.openmm_equilibrate import ClassicalPrepConfig, MoleculeSpec
from src.prep.structure_io import PreparedStructure
from src.units import (
    ANGSTROM_PER_NM,
    KCAL_PER_KJ,
    KJ_PER_KCAL,
    NM_PER_ANGSTROM,
)


# ── decision D-4: shared atom ordering ────────────────────────────────────────

def test_rdkit_mol_and_rdkit_3d_share_atom_order() -> None:
    """_rdkit_3d species/positions must match _rdkit_mol's AddHs atom order."""
    mol = _rdkit_mol(_MONOMER_SMILES, seed=43)
    positions, species = _rdkit_3d(_MONOMER_SMILES, seed=43)

    mol_symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    assert species == mol_symbols
    assert positions.shape == (len(species), 3)
    np.testing.assert_allclose(positions, mol.GetConformer().GetPositions())


def test_rdkit_mol_ordering_is_seed_independent() -> None:
    """Conformer seed changes coordinates but not atom indexing/ordering."""
    sym_a = [a.GetSymbol() for a in _rdkit_mol(_MONOMER_SMILES, seed=1).GetAtoms()]
    sym_b = [a.GetSymbol() for a in _rdkit_mol(_MONOMER_SMILES, seed=2).GetAtoms()]
    assert sym_a == sym_b


# ── structure handoff file ────────────────────────────────────────────────────

def test_prepared_structure_roundtrip(tmp_path) -> None:
    rng = np.random.default_rng(0)
    positions = rng.standard_normal((5, 3))
    species = ['C', 'H', 'H', 'O', 'N']
    cell = np.diag([10.0, 10.0, 10.0])
    s = PreparedStructure(positions, species, cell, metadata={'k': 'v'})

    path = tmp_path / 's.json'
    s.save(path)
    loaded = PreparedStructure.load(path)

    np.testing.assert_allclose(loaded.positions, positions)
    assert loaded.species == species
    np.testing.assert_allclose(loaded.cell, cell)
    assert loaded.metadata == {'k': 'v'}
    assert loaded.n_atoms == 5


def test_prepared_structure_nonperiodic_roundtrip(tmp_path) -> None:
    s = PreparedStructure(np.zeros((2, 3)), ['C', 'C'], cell=None)
    path = tmp_path / 's.json'
    s.save(path)
    loaded = PreparedStructure.load(path)
    assert loaded.cell is None


def test_prepared_structure_species_length_mismatch() -> None:
    with pytest.raises(ValueError, match='species length'):
        PreparedStructure(np.zeros((3, 3)), ['C', 'C'])


def test_prepared_structure_bad_cell_shape() -> None:
    with pytest.raises(ValueError, match='cell must be'):
        PreparedStructure(np.zeros((2, 3)), ['C', 'C'], cell=np.zeros((2, 2)))


def test_prepared_structure_bad_positions_shape() -> None:
    with pytest.raises(ValueError, match='positions must be'):
        PreparedStructure(np.zeros((2, 4)), ['C', 'C'])


def test_prepared_structure_rejects_unknown_schema(tmp_path) -> None:
    path = tmp_path / 's.json'
    path.write_text('{"schema_version": 999, "species": [], "positions_A": []}',
                    encoding='utf-8')
    with pytest.raises(ValueError, match='schema_version'):
        PreparedStructure.load(path)


# ── unit conversions ──────────────────────────────────────────────────────────

def test_length_conversions_are_inverses() -> None:
    assert ANGSTROM_PER_NM == pytest.approx(10.0)
    assert NM_PER_ANGSTROM == pytest.approx(0.1)
    assert ANGSTROM_PER_NM * NM_PER_ANGSTROM == pytest.approx(1.0)


def test_energy_conversions_are_inverses() -> None:
    assert KJ_PER_KCAL == pytest.approx(4.184)
    assert KJ_PER_KCAL * KCAL_PER_KJ == pytest.approx(1.0)


# ── config defaults ───────────────────────────────────────────────────────────

def test_prep_config_defaults_match_paper_anchors() -> None:
    cfg = ClassicalPrepConfig()
    assert cfg.target_density_g_per_ml == 0.5   # paper SI S-3 initial density
    assert cfg.temperature_K == 333.0           # production setpoint
    assert cfg.protocol == 'simple'             # decision D-3
    assert cfg.charge_method == 'nagl'          # decision D-2


def test_molecule_spec_fields() -> None:
    spec = MoleculeSpec(smiles=_MONOMER_SMILES, count=8, rdkit_seed=43)
    assert spec.count == 8
    assert spec.rdkit_seed == 43
