"""Tests for reproduce_figures helpers.

Focus: the branching-detection gate for the Carothers Fig. 4c DPn plot
(:func:`kagome.analysis.network.max_inter_monomer_degree` /
:func:`is_branching_topology`). Carothers DPn = 1/(1-p) is linear step-growth
theory, so a branched network (functionality f>2) must skip the plot while a
linear chain (nylon) still generates it (approved 2026-07-30).

Pure graph tests — no MD, torch, or matplotlib needed. Each monomer is a pair
of atoms {2m, 2m+1} joined by one intramonomer bond; monomer membership is
recovered from that initial topology exactly as the figure pipeline does.
"""
from kagome.analysis.carothers import monomer_sets_from_bonds
from kagome.analysis.network import (
    is_branching_topology,
    max_inter_monomer_degree,
)

# Four 2-atom monomers: M0={0,1}, M1={2,3}, M2={4,5}, M3={6,7}.
_INITIAL_4 = [(0, 1), (2, 3), (4, 5), (6, 7)]
# Two 2-atom monomers: M0={0,1}, M1={2,3}.
_INITIAL_2 = [(0, 1), (2, 3)]


class TestMaxInterMonomerDegree:

    def test_linear_chain_max_degree_two(self):
        # Line M0-M1-M2-M3: interior monomers bond to two neighbours.
        monomer_sets = monomer_sets_from_bonds(_INITIAL_4)
        final = _INITIAL_4 + [(1, 2), (3, 4), (5, 6)]
        assert max_inter_monomer_degree(final, monomer_sets) == 2

    def test_star_branch_max_degree_three(self):
        # Star: M0 bonds to M1, M2, M3 -> central monomer degree 3.
        monomer_sets = monomer_sets_from_bonds(_INITIAL_4)
        final = _INITIAL_4 + [(0, 2), (0, 4), (0, 6)]
        assert max_inter_monomer_degree(final, monomer_sets) == 3

    def test_double_bond_same_pair_counts_once(self):
        # Two inter-monomer bonds between the SAME pair (nylon amide + paired
        # water-forming bond) -> distinct neighbours, degree 1 not 2.
        monomer_sets = monomer_sets_from_bonds(_INITIAL_2)
        final = _INITIAL_2 + [(1, 2), (0, 3)]
        assert max_inter_monomer_degree(final, monomer_sets) == 1

    def test_no_inter_monomer_bonds_is_zero(self):
        monomer_sets = monomer_sets_from_bonds(_INITIAL_4)
        assert max_inter_monomer_degree(_INITIAL_4, monomer_sets) == 0

    def test_no_monomers_is_zero(self):
        assert max_inter_monomer_degree([], []) == 0


class TestIsBranchingTopology:

    def test_linear_chain_not_branching(self):
        monomer_sets = monomer_sets_from_bonds(_INITIAL_4)
        final = _INITIAL_4 + [(1, 2), (3, 4), (5, 6)]
        assert is_branching_topology(final, monomer_sets) is False

    def test_star_is_branching(self):
        monomer_sets = monomer_sets_from_bonds(_INITIAL_4)
        final = _INITIAL_4 + [(0, 2), (0, 4), (0, 6)]
        assert is_branching_topology(final, monomer_sets) is True

    def test_nylon_double_bond_pair_not_branching(self):
        monomer_sets = monomer_sets_from_bonds(_INITIAL_2)
        final = _INITIAL_2 + [(1, 2), (0, 3)]
        assert is_branching_topology(final, monomer_sets) is False
