"""Tests for the minimal LAMMPS data-file parser (E2 external comparison)."""
import pytest

from kagome.io.lammps_data import (
    LammpsData,
    element_from_mass,
    read_lammps_data,
)

# Synthetic minimal data file in the exact shape of the Provenzano 2025
# structures: title line, header counts/bounds, 'Atoms # full' style comment,
# unordered atom ids, skipped sections (Pair Coeffs, Velocities).
_MINIMAL = """\
Test LAMMPS data file via write_data, units = real

4 atoms
2 bonds
3 atom types

-1.0 9.0 xlo xhi
0.0 10.0 ylo yhi
0.5 10.5 zlo zhi

Masses

1 12.01115
2 1.00797
3 15.9994

Pair Coeffs # lj/class2/coul/long

1 0.068 3.915 12
2 0.02 2.99 12
3 0.08 3.6 12

Atoms # full

2 1 2 0.100 1.0 2.0 3.0 0 0 0
1 1 1 -0.200 0.0 0.0 0.0 0 0 0
4 1 3 -0.300 2.0 2.0 2.0
3 1 1 0.000 1.0 1.0 1.0 0 0 0

Velocities

1 0.1 0.2 0.3
2 0.0 0.0 0.0
3 0.0 0.0 0.0
4 0.0 0.0 0.0

Bonds

1 1 1 2
2 2 3 4
"""


class TestElementFromMass:

    def test_class_ii_masses(self):
        # Class II force-field masses from the real Provenzano files.
        assert element_from_mass(12.01115) == 'C'
        assert element_from_mass(1.00797) == 'H'
        assert element_from_mass(14.0067) == 'N'
        assert element_from_mass(15.9994) == 'O'

    def test_unknown_mass_raises(self):
        # Fe (55.85 amu) matches nothing within the tolerance.
        with pytest.raises(ValueError, match='matches no element'):
            element_from_mass(55.85)

    def test_borderline_outside_tolerance_raises(self):
        # Halfway between C (12.011) and N (14.007) is > 0.5 amu from both.
        with pytest.raises(ValueError):
            element_from_mass(13.0)


class TestReadLammpsData:

    @pytest.fixture()
    def data(self, tmp_path) -> LammpsData:
        path = tmp_path / 'data.test'
        path.write_text(_MINIMAL, encoding='utf-8')
        return read_lammps_data(path)

    def test_species_inferred_and_ordered_by_atom_id(self, data):
        # Atom rows are deliberately unordered (2, 1, 4, 3); species must be
        # indexed by 0-based atom id: id1=C, id2=H, id3=C, id4=O.
        assert data.species == ['C', 'H', 'C', 'O']

    def test_bonds_converted_to_zero_indexed(self, data):
        assert data.bonds == [(0, 1), (2, 3)]

    def test_counts_and_box(self, data):
        assert data.n_atoms == 4
        assert data.box_bounds == ((-1.0, 9.0), (0.0, 10.0), (0.5, 10.5))

    def test_skipped_sections_do_not_leak(self, data):
        # Pair Coeffs and Velocities rows must not be parsed as atoms/bonds.
        assert len(data.bonds) == 2
        assert len(data.species) == 4

    def test_unknown_mass_in_file_raises(self, tmp_path):
        path = tmp_path / 'data.bad'
        path.write_text(
            _MINIMAL.replace('1 12.01115', '1 55.85'), encoding='utf-8')
        with pytest.raises(ValueError, match='matches no element'):
            read_lammps_data(path)

    def test_atom_count_mismatch_raises(self, tmp_path):
        path = tmp_path / 'data.short'
        path.write_text(
            _MINIMAL.replace('3 1 1 0.000 1.0 1.0 1.0 0 0 0\n', ''),
            encoding='utf-8')
        with pytest.raises(ValueError, match='atoms'):
            read_lammps_data(path)

    def test_bond_count_mismatch_raises(self, tmp_path):
        path = tmp_path / 'data.bonds'
        path.write_text(
            _MINIMAL.replace('2 2 3 4\n', ''), encoding='utf-8')
        with pytest.raises(ValueError, match='[Bb]onds'):
            read_lammps_data(path)

    def test_non_full_atom_row_raises(self, tmp_path):
        # 'atomic' style rows (id type x y z) have < 7 columns.
        path = tmp_path / 'data.atomic'
        path.write_text(
            _MINIMAL.replace('2 1 2 0.100 1.0 2.0 3.0 0 0 0',
                             '2 2 1.0 2.0 3.0'),
            encoding='utf-8')
        with pytest.raises(ValueError, match='full'):
            read_lammps_data(path)
