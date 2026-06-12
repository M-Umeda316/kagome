"""Tests for shared ethylene system builders."""
import numpy as np
import pytest


class TestBuildEthyleneBox:

    @pytest.fixture(autouse=True)
    def _skip_no_ase(self):
        pytest.importorskip('ase')

    def test_correct_atom_count(self):
        from scripts._systems import build_ethylene_box
        rng = np.random.default_rng(0)
        positions, species = build_ethylene_box(3, 12.0, rng)
        assert positions.shape == (18, 3)  # 3 molecules × 6 atoms (C2H4)
        assert len(species) == 18
        assert species.count('C') == 6
        assert species.count('H') == 12

    def test_no_overlap(self):
        from scripts._systems import build_ethylene_box
        rng = np.random.default_rng(42)
        atoms_per_mol = 6
        n_mol = 4
        positions, _ = build_ethylene_box(n_mol, 14.0, rng, min_sep=2.5)
        # Check inter-molecular distances only (intra-mol bonds are ~1.3 Å)
        for mol_i in range(n_mol):
            for mol_j in range(mol_i + 1, n_mol):
                s_i = mol_i * atoms_per_mol
                s_j = mol_j * atoms_per_mol
                mol_i_pos = positions[s_i:s_i + atoms_per_mol]
                mol_j_pos = positions[s_j:s_j + atoms_per_mol]
                diffs = mol_i_pos[:, np.newaxis, :] - mol_j_pos[np.newaxis, :, :]
                d_min = float(np.min(np.linalg.norm(diffs, axis=2)))
                assert d_min >= 2.5, f'Inter-mol overlap: mol {mol_i} and {mol_j} at {d_min:.3f} Å'

    def test_too_small_box_raises(self):
        from scripts._systems import build_ethylene_box
        rng = np.random.default_rng(0)
        with pytest.raises(RuntimeError, match='without overlap'):
            build_ethylene_box(10, 4.0, rng)

    def test_reproducible_with_seed(self):
        from scripts._systems import build_ethylene_box
        pos1, sp1 = build_ethylene_box(2, 10.0, np.random.default_rng(7))
        pos2, sp2 = build_ethylene_box(2, 10.0, np.random.default_rng(7))
        np.testing.assert_array_equal(pos1, pos2)
        assert sp1 == sp2


class TestBuildTemplateAndGroups:

    def test_group_sizes(self):
        from scripts._systems import build_template_and_groups
        template, groups = build_template_and_groups(4)
        assert len(groups['C_donor'].atom_indices) == 4
        assert len(groups['C_acceptor'].atom_indices) == 4

    def test_correct_indices(self):
        from scripts._systems import build_template_and_groups
        _, groups = build_template_and_groups(3)
        assert groups['C_donor'].atom_indices == [0, 6, 12]
        assert groups['C_acceptor'].atom_indices == [1, 7, 13]
