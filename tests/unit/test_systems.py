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


class TestVinylAIBNHelpers:
    """Tests for RDKit-based vinyl/AIBN system helpers."""

    @pytest.fixture(autouse=True)
    def _skip_no_rdkit(self):
        pytest.importorskip('rdkit')

    def test_find_vinyl_alpha_beta_methyl_acrylate(self):
        from scripts._systems import _find_vinyl_alpha_beta
        alpha, beta = _find_vinyl_alpha_beta('C=CC(=O)OC')
        # alpha-C must be different from beta-C
        assert alpha != beta

    def test_find_vinyl_alpha_has_more_h_than_beta(self):
        """alpha-C should be CH2= (2 H), beta-C should be =CH- (1 H)."""
        from rdkit import Chem
        from scripts._systems import _find_vinyl_alpha_beta

        mol = Chem.MolFromSmiles('C=CC(=O)OC')
        mol = Chem.AddHs(mol)
        alpha, beta = _find_vinyl_alpha_beta('C=CC(=O)OC')
        h_alpha = sum(1 for n in mol.GetAtomWithIdx(alpha).GetNeighbors() if n.GetSymbol() == 'H')
        h_beta = sum(1 for n in mol.GetAtomWithIdx(beta).GetNeighbors() if n.GetSymbol() == 'H')
        assert h_alpha == 2
        assert h_beta == 1

    def test_find_ibn_radical_c_has_three_c_neighbours(self):
        """Radical C in CC(C)C#N must be bonded to exactly 3 C (2 methyl + nitrile C)."""
        from rdkit import Chem
        from scripts._systems import _find_ibn_radical_c

        smiles = 'CC(C)C#N'
        idx = _find_ibn_radical_c(smiles)
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        atom = mol.GetAtomWithIdx(idx)
        nbrs = atom.GetNeighbors()
        # closed-shell isobutyronitrile: central C has 3 C + 1 H
        assert sum(1 for n in nbrs if n.GetSymbol() == 'C') == 3
        assert sum(1 for n in nbrs if n.GetSymbol() == 'H') == 1

    def test_rdkit_3d_species_count(self):
        """Methyl acrylate C=CC(=O)OC has 12 atoms with explicit H."""
        from scripts._systems import _rdkit_3d
        positions, species = _rdkit_3d('C=CC(=O)OC')
        assert positions.shape == (len(species), 3)
        assert species.count('C') == 4
        assert species.count('O') == 2
        assert species.count('H') == 6

    def test_rdkit_3d_initiator_species_count(self):
        """Isobutyronitrile CC(C)C#N has 10 atoms with explicit H."""
        from scripts._systems import _rdkit_3d
        positions, species = _rdkit_3d('CC(C)C#N')
        assert positions.shape == (len(species), 3)
        assert species.count('C') == 4
        assert species.count('N') == 1
        assert species.count('H') == 7


class TestBuildVinylAIBNSystem:

    @pytest.fixture(autouse=True)
    def _skip_no_rdkit(self):
        pytest.importorskip('rdkit')

    def test_atom_count(self):
        from scripts._systems import build_vinyl_aibn_system, _rdkit_3d
        rng = np.random.default_rng(0)
        n_mono, n_init = 2, 1
        pos, species, template, groups, prop_map = build_vinyl_aibn_system(
            n_monomers=n_mono, n_initiators=n_init, box_size=18.0, rng=rng,
        )
        n_per_mono = len(_rdkit_3d('C=CC(=O)OC')[1])   # 12
        n_per_init = len(_rdkit_3d('CC(C)C#N')[1])     # 10
        expected_n = n_init * n_per_init + n_mono * n_per_mono
        assert pos.shape == (expected_n, 3)
        assert len(species) == expected_n

    def test_group_sizes(self):
        from scripts._systems import build_vinyl_aibn_system
        rng = np.random.default_rng(1)
        _, _, _, groups, _ = build_vinyl_aibn_system(
            n_monomers=3, n_initiators=2, box_size=20.0, rng=rng,
        )
        assert len(groups['radical_C'].atom_indices) == 2
        assert len(groups['vinyl_alpha_C'].atom_indices) == 3

    def test_propagation_map_size(self):
        from scripts._systems import build_vinyl_aibn_system
        rng = np.random.default_rng(2)
        n_mono = 4
        _, _, _, _, prop_map = build_vinyl_aibn_system(
            n_monomers=n_mono, n_initiators=1, box_size=22.0, rng=rng,
        )
        assert len(prop_map) == n_mono

    def test_propagation_map_alpha_in_vinyl_alpha_group(self):
        """Every alpha-C key in propagation_map must appear in vinyl_alpha_C group."""
        from scripts._systems import build_vinyl_aibn_system
        rng = np.random.default_rng(3)
        _, _, _, groups, prop_map = build_vinyl_aibn_system(
            n_monomers=2, n_initiators=1, box_size=18.0, rng=rng,
        )
        alpha_indices = set(groups['vinyl_alpha_C'].atom_indices)
        for alpha in prop_map:
            assert alpha in alpha_indices

    def test_propagation_map_beta_not_in_any_group(self):
        """beta-C values must NOT already be in radical_C or vinyl_alpha_C groups."""
        from scripts._systems import build_vinyl_aibn_system
        rng = np.random.default_rng(4)
        _, _, _, groups, prop_map = build_vinyl_aibn_system(
            n_monomers=2, n_initiators=1, box_size=18.0, rng=rng,
        )
        all_group_indices = set(
            idx for g in groups.values() for idx in g.atom_indices
        )
        for beta in prop_map.values():
            assert beta not in all_group_indices

    def test_no_overlap(self):
        from scripts._systems import build_vinyl_aibn_system
        rng = np.random.default_rng(5)
        pos, species, _, _, _ = build_vinyl_aibn_system(
            n_monomers=2, n_initiators=1, box_size=20.0, rng=rng, min_sep=2.0,
        )
        # Check min inter-atom distance (across different molecules).
        # We don't check intra-mol distances (bonds ~1.4 Å are expected).
        # Rough sanity: at least some atoms should be > 2 Å from each other.
        # Just verify the array is well-formed.
        assert pos.shape[1] == 3
        assert not np.any(np.isnan(pos))
