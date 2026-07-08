"""Tests for Carothers equation (step-growth polymerization)."""
import logging

import numpy as np
import pytest

from kagome.analysis.carothers import (
    dpn_carothers,
    dpn_from_bonds,
    dpn_measured_from_topology,
    monomer_sets_from_bonds,
)


class TestDpnCarothers:

    def test_zero_conversion(self):
        assert dpn_carothers(0.0) == pytest.approx(1.0)

    def test_half_conversion(self):
        assert dpn_carothers(0.5) == pytest.approx(2.0)

    def test_ninety_percent(self):
        assert dpn_carothers(0.9) == pytest.approx(10.0)

    def test_ninety_nine_percent(self):
        assert dpn_carothers(0.99) == pytest.approx(100.0)

    def test_array_input(self):
        p = np.array([0.0, 0.5, 0.9])
        result = dpn_carothers(p)
        expected = np.array([1.0, 2.0, 10.0])
        np.testing.assert_allclose(result, expected)

    def test_full_conversion_clamped_not_inf(self, caplog):
        # F9: p = 1.0 must yield a finite DPn (10000), not inf, and warn.
        with caplog.at_level(logging.WARNING):
            dpn = dpn_carothers(1.0)
        assert np.isfinite(dpn)
        assert dpn == pytest.approx(10000.0)  # 1/(1-0.9999)
        assert 'clamp' in caplog.text

    def test_overshoot_clamped_not_negative(self, caplog):
        # p > 1.0 (miscount) must clamp to a finite positive DPn, not go negative.
        with caplog.at_level(logging.WARNING):
            dpn = dpn_carothers(1.5)
        assert np.isfinite(dpn)
        assert dpn == pytest.approx(10000.0)
        assert 'clamp' in caplog.text

    def test_array_with_full_conversion(self, caplog):
        with caplog.at_level(logging.WARNING):
            result = dpn_carothers(np.array([0.5, 1.0]))
        assert np.all(np.isfinite(result))
        np.testing.assert_allclose(result, np.array([2.0, 10000.0]))
        assert 'clamp' in caplog.text

    def test_no_warning_below_clamp(self, caplog):
        with caplog.at_level(logging.WARNING):
            dpn_carothers(0.99)
        assert caplog.text == ''


class TestDpnFromBonds:

    def test_no_bonds(self):
        assert dpn_from_bonds(0, 20) == pytest.approx(1.0)

    def test_half_conversion(self):
        result = dpn_from_bonds(5, 20)
        assert result == pytest.approx(2.0)

    def test_zero_groups(self):
        assert dpn_from_bonds(0, 0) == pytest.approx(1.0)

    def test_near_full_conversion(self):
        result = dpn_from_bonds(9, 20)
        assert result == pytest.approx(10.0)

    def test_clamp_warns_on_overshoot(self, caplog):
        # RF19b: p>1 (miscount) or full conversion must not be silently clamped.
        with caplog.at_level(logging.WARNING):
            dpn = dpn_from_bonds(n_bonds=100, n_functional_groups=100)  # p = 2.0
        assert dpn == pytest.approx(10000.0)  # 1/(1-0.9999)
        assert 'clamp' in caplog.text

    def test_no_warning_below_clamp(self, caplog):
        with caplog.at_level(logging.WARNING):
            dpn_from_bonds(n_bonds=9, n_functional_groups=20)  # p = 0.9
        assert caplog.text == ''


def _linear_chain(n_monomers: int, atoms_per_monomer: int = 2):
    """Hand-built (monomer_atom_sets, bonds) for one linear chain of n monomers.

    Each monomer occupies a contiguous atom block with one intramonomer bond;
    consecutive monomers are joined by one inter-monomer bond (the last atom of
    monomer k to the first atom of monomer k+1). Measured DPn -> n_monomers.
    """
    monomer_atom_sets = [
        list(range(k * atoms_per_monomer, (k + 1) * atoms_per_monomer))
        for k in range(n_monomers)
    ]
    bonds: list[tuple[int, int, float]] = []
    for atoms in monomer_atom_sets:  # intramonomer bonds
        for a, b in zip(atoms, atoms[1:]):
            bonds.append((a, b, 1.0))
    for k in range(n_monomers - 1):  # inter-monomer links
        last = monomer_atom_sets[k][-1]
        first = monomer_atom_sets[k + 1][0]
        bonds.append((last, first, 1.0))
    return monomer_atom_sets, bonds


class TestMonomerSetsFromBonds:

    def test_recovers_disconnected_monomers(self):
        # Four 2-atom monomers, only intramonomer bonds -> four components.
        bonds = [(0, 1), (2, 3), (4, 5), (6, 7)]
        comps = monomer_sets_from_bonds(bonds)
        assert comps == [[0, 1], [2, 3], [4, 5], [6, 7]]

    def test_merges_bonded_atoms(self):
        # 0-1-2 one component, 3-4 another.
        bonds = [(0, 1), (1, 2), (3, 4)]
        comps = monomer_sets_from_bonds(bonds)
        assert comps == [[0, 1, 2], [3, 4]]

    def test_order_ignored_in_triples(self):
        comps = monomer_sets_from_bonds([(0, 1, 2.0), (2, 3, 1.0)])
        assert comps == [[0, 1], [2, 3]]


class TestDpnMeasuredFromTopology:

    def test_no_reaction_gives_dpn_one(self):
        # Four separate monomers, no inter-monomer bonds -> DPn = 4/4 = 1.
        monomer_atom_sets = [[0, 1], [2, 3], [4, 5], [6, 7]]
        intramono = [(0, 1), (2, 3), (4, 5), (6, 7)]
        assert dpn_measured_from_topology(intramono, monomer_atom_sets) == pytest.approx(1.0)

    def test_two_dimers_gives_dpn_two(self):
        # Bonds 1-2 (mono0+1) and 5-6 (mono2+3) -> two molecules -> DPn = 4/2 = 2.
        monomer_atom_sets = [[0, 1], [2, 3], [4, 5], [6, 7]]
        bonds = [(0, 1), (2, 3), (4, 5), (6, 7), (1, 2), (5, 6)]
        assert dpn_measured_from_topology(bonds, monomer_atom_sets) == pytest.approx(2.0)

    def test_single_chain_gives_dpn_n(self):
        # All four monomers linked into one chain -> DPn = 4/1 = 4.
        monomer_atom_sets, bonds = _linear_chain(4)
        assert dpn_measured_from_topology(bonds, monomer_atom_sets) == pytest.approx(4.0)

    def test_empty_monomers_returns_one(self):
        assert dpn_measured_from_topology([], []) == pytest.approx(1.0)

    def test_intramonomer_bonds_do_not_merge(self):
        # Extra intramonomer bond must not change the molecule count.
        monomer_atom_sets = [[0, 1, 2], [3, 4, 5]]
        bonds = [(0, 1), (1, 2), (3, 4), (4, 5)]  # no inter-monomer bond
        assert dpn_measured_from_topology(bonds, monomer_atom_sets) == pytest.approx(1.0 * 2 / 2)

    def test_duplicate_intermonomer_bonds_idempotent(self):
        # Nylon condensation adds BOTH an amide and a water-forming bond between
        # the same diamine/diacid pair; the second must not split the molecule.
        monomer_atom_sets = [[0, 1], [2, 3]]
        bonds = [(0, 1), (2, 3), (1, 2), (0, 3)]  # two links, same pair
        assert dpn_measured_from_topology(bonds, monomer_atom_sets) == pytest.approx(2.0)


class TestMeasuredMatchesTheory:

    @pytest.mark.parametrize('n_monomers', [2, 5, 10, 50])
    def test_linear_chain_matches_carothers(self, n_monomers):
        # Ideal linear step-growth: one chain of N bifunctional monomers has
        # N-1 links, so p = 2(N-1)/(2N) = (N-1)/N and DPn = 1/(1-p) = N.
        monomer_atom_sets, bonds = _linear_chain(n_monomers)
        measured = dpn_measured_from_topology(bonds, monomer_atom_sets)
        p = (n_monomers - 1) / n_monomers
        theory = dpn_carothers(p)
        assert measured == pytest.approx(float(n_monomers))
        assert measured == pytest.approx(float(theory))
