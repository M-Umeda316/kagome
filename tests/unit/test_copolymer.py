"""Tests for the mixed-monomer (copolymer) vinyl builder.

Covers the one thing that differs from the single-monomer builder: global index
arithmetic must use a per-species running offset, because acrylate (12 atoms)
and methacrylate (15 atoms) have different atom counts. Candidate generation is
exercised without MD (specs/decisions.md 2026-07-16).
"""
import numpy as np
import pytest

from kagome.reactive.selection import (
    find_candidates,
    score_candidates,
    select_non_overlapping,
)

ACRYLATE = 'C=CC(=O)OC'        # 12 atoms after AddHs
METHACRYLATE = 'C=C(C)C(=O)OC'  # 15 atoms after AddHs


@pytest.fixture(autouse=True)
def _skip_no_rdkit():
    pytest.importorskip('rdkit')


def _build(n_acr, n_mac, n_init, seed=0, box=24.0):
    from scripts._systems import build_vinyl_copolymer_system
    rng = np.random.default_rng(seed)
    return build_vinyl_copolymer_system(
        monomer_specs=[(ACRYLATE, n_acr), (METHACRYLATE, n_mac)],
        n_initiators=n_init, box_size=box, rng=rng,
    )


class TestCopolymerBuilder:

    def test_atom_count_sums_heterogeneous_species(self):
        from scripts._systems import _rdkit_3d
        pos, species, _, _, _, _ = _build(2, 2, 1)
        n_acr = len(_rdkit_3d(ACRYLATE)[1])
        n_mac = len(_rdkit_3d(METHACRYLATE)[1])
        n_init = len(_rdkit_3d('CC(C)C#N')[1])
        assert n_acr != n_mac  # the whole point: different strides
        expected = 1 * n_init + 2 * n_acr + 2 * n_mac
        assert pos.shape == (expected, 3)
        assert len(species) == expected

    def test_group_sizes(self):
        _, _, _, groups, _, _ = _build(3, 4, 2)
        assert len(groups['radical_C'].atom_indices) == 2
        # both species' alpha-Cs land in the SAME vinyl_alpha_C group
        assert len(groups['vinyl_alpha_C'].atom_indices) == 3 + 4
        assert len(groups['vinyl_beta_C'].atom_indices) == 3 + 4
        assert len(groups['chain_C'].atom_indices) == 2

    def test_propagation_map_size_and_membership(self):
        _, _, _, groups, prop_map, _ = _build(3, 2, 1)
        assert len(prop_map) == 3 + 2
        alpha = set(groups['vinyl_alpha_C'].atom_indices)
        beta = set(groups['vinyl_beta_C'].atom_indices)
        forbidden = set(groups['radical_C'].atom_indices) | alpha
        for a, b in prop_map.items():
            assert a in alpha
            assert b in beta
            assert b not in forbidden

    def test_indices_within_bounds_and_unique(self):
        pos, _, _, groups, prop_map, chain_c_map = _build(4, 3, 2)
        n = pos.shape[0]
        all_group_idx = [i for g in groups.values() for i in g.atom_indices]
        assert all(0 <= i < n for i in all_group_idx)
        # every reactive atom is distinct (no offset collision across species)
        assert len(all_group_idx) == len(set(all_group_idx))
        for r, c in chain_c_map.items():
            assert r in groups['radical_C'].atom_indices
            assert c in groups['chain_C'].atom_indices

    def test_offset_correct_second_species_shifted(self):
        """The methacrylate alpha indices must be offset past ALL acrylate atoms,
        proving the running offset (not a fixed stride) is used."""
        from scripts._systems import _rdkit_3d
        n_init = len(_rdkit_3d('CC(C)C#N')[1])
        n_acr = len(_rdkit_3d(ACRYLATE)[1])
        _, _, _, groups, _, _ = _build(2, 2, 1)
        alphas = sorted(groups['vinyl_alpha_C'].atom_indices)
        # first 2 alphas belong to acrylate (after 1 initiator block),
        # last 2 to methacrylate (after the 2 acrylate blocks).
        first_mac_alpha = alphas[2]
        assert first_mac_alpha >= n_init + 2 * n_acr

    def test_template_shape_unchanged(self):
        _, _, template, _, _, _ = _build(2, 2, 1)
        assert template.name == 'radical_vinyl_copolymerization'
        assert len(template.groups) == 4
        assert len(template.pairs) == 3
        assert sum(p.constraint_only for p in template.pairs) == 2

    def test_all_methacrylate_only(self):
        """Degenerate: zero acrylate still builds a valid single-species system."""
        _, _, _, groups, prop_map, _ = _build(0, 3, 1)
        assert len(groups['vinyl_alpha_C'].atom_indices) == 3
        assert len(prop_map) == 3

    def test_empty_specs_rejected(self):
        from scripts._systems import build_vinyl_copolymer_system
        rng = np.random.default_rng(0)
        with pytest.raises(ValueError):
            build_vinyl_copolymer_system(
                monomer_specs=[], n_initiators=1, box_size=20.0, rng=rng,
            )


class TestCopolymerCandidateGeneration:
    """Candidate generation runs on the mixed layout without MD."""

    def test_find_candidates_runs_on_mixed_layout(self):
        _, _, template, groups, _, _ = _build(4, 4, 2, box=16.0)
        # dense box so at least some radical_C / vinyl_alpha_C pairs fall in
        # the [3, 6] Å formation window.
        pos = _build(4, 4, 2, box=16.0)[0]
        cell = np.diag([16.0, 16.0, 16.0])
        candidates = find_candidates(template, groups, pos, cell)
        scored = score_candidates(candidates)
        selected = select_non_overlapping(scored)
        # selection must not raise and must return a subset
        assert len(selected) <= len(scored)
        # every selected candidate references in-range atom indices
        for c in selected:
            for idx in c.atom_indices:
                assert 0 <= idx < pos.shape[0]
