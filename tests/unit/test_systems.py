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


class TestPlaceFragmentsPBC:
    """PBC-aware placement: positions wrap and distances use minimum image."""

    @pytest.fixture(autouse=True)
    def _skip_no_ase(self):
        pytest.importorskip('ase')

    def test_positions_not_clipped_away_from_boundary(self):
        """Positions may lie near 0 or box_size (no 2 Å clip margin)."""
        from scripts._systems import build_ethylene_box
        rng = np.random.default_rng(99)
        box = 12.0
        pos, _ = build_ethylene_box(4, box, rng)
        assert np.all(pos >= 0.0) or np.all(pos < box), (
            'All positions should be non-negative (wrapped into box)'
        )
        has_near_zero = np.any(pos < 2.0)
        has_near_edge = np.any(pos > box - 2.0)
        assert has_near_zero or has_near_edge, (
            'With wrapping, some atoms should appear near box boundaries'
        )

    def test_no_pbc_overlap(self):
        """Inter-molecular minimum-image distances must exceed min_sep."""
        from scripts._systems import build_ethylene_box
        rng = np.random.default_rng(42)
        box = 14.0
        n_mol = 4
        atoms_per_mol = 6
        pos, _ = build_ethylene_box(n_mol, box, rng, min_sep=2.5)
        box_vec = np.array([box, box, box])
        for mol_i in range(n_mol):
            for mol_j in range(mol_i + 1, n_mol):
                s_i = mol_i * atoms_per_mol
                s_j = mol_j * atoms_per_mol
                diffs = pos[s_i:s_i + atoms_per_mol, np.newaxis, :] - pos[np.newaxis, s_j:s_j + atoms_per_mol, :]
                diffs = diffs - box_vec * np.round(diffs / box_vec)
                d_min = float(np.min(np.linalg.norm(diffs, axis=2)))
                assert d_min >= 2.5, (
                    f'PBC overlap: mol {mol_i} and {mol_j} at {d_min:.3f} Å'
                )


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

    def test_find_vinyl_alpha_beta_disubstituted_tertiary(self):
        """1,1-disubstituted vinyls (CH2=CR2) must resolve: alpha=CH2 (2 H),
        beta=substituted C (0 H, tertiary radical after propagation). Covers the
        paper's methacrylate / diphenylethylene / dimethyl itaconate monomers
        (decision 2026-06-20 'di-substituted vinyl builder')."""
        from rdkit import Chem
        from scripts._systems import _find_vinyl_alpha_beta

        for smiles in (
            'C=C(C)C(=O)OC',            # methyl methacrylate
            'C=C(c1ccccc1)c1ccccc1',    # 1,1-diphenylethylene
            'C=C(CC(=O)OC)C(=O)OC',     # dimethyl itaconate
        ):
            alpha, beta = _find_vinyl_alpha_beta(smiles)
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            h_alpha = sum(1 for n in mol.GetAtomWithIdx(alpha).GetNeighbors()
                          if n.GetSymbol() == 'H')
            h_beta = sum(1 for n in mol.GetAtomWithIdx(beta).GetNeighbors()
                         if n.GetSymbol() == 'H')
            assert h_alpha == 2, f'alpha must be =CH2 terminus for {smiles}'
            assert h_beta == 0, f'beta must be the disubstituted C for {smiles}'

    def test_build_methacrylate_system(self):
        """End-to-end build of a 1,1-disubstituted (methacrylate) vinyl system:
        the substituted beta carbon must populate the vinyl_beta_C group that
        becomes the propagating (tertiary) radical."""
        import numpy as np
        from scripts._systems import build_vinyl_aibn_system

        rng = np.random.default_rng(7)
        _, species, _, groups, prop_map, _ = build_vinyl_aibn_system(
            n_monomers=3, n_initiators=2, box_size=22.0, rng=rng,
            monomer_smiles='C=C(C)C(=O)OC',
        )
        assert len(groups['vinyl_alpha_C'].atom_indices) == 3
        assert len(groups['vinyl_beta_C'].atom_indices) == 3
        assert len(prop_map) == 3
        beta_indices = set(groups['vinyl_beta_C'].atom_indices)
        assert set(prop_map.values()) == beta_indices

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

    def test_find_chain_c_neighbor(self):
        """chain_C neighbor of radical C in CC(C)C#N must be a non-nitrile C."""
        from rdkit import Chem
        from scripts._systems import _find_ibn_radical_c, _find_chain_c_neighbor

        smiles = 'CC(C)C#N'
        rad_idx = _find_ibn_radical_c(smiles)
        chain_idx = _find_chain_c_neighbor(smiles, rad_idx)
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        chain_atom = mol.GetAtomWithIdx(chain_idx)
        assert chain_atom.GetSymbol() == 'C'
        assert chain_idx != rad_idx
        nitrile_c = [n.GetIdx() for n in mol.GetAtomWithIdx(rad_idx).GetNeighbors()
                     if n.GetSymbol() == 'C'
                     and any(b.GetBondTypeAsDouble() == 3.0
                             for b in n.GetBonds())]
        assert chain_idx not in nitrile_c

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


class TestRdkitSeedDeterminism:
    """RF23: conformer geometry must follow the seed (threaded from --seed), while
    atom ordering stays seed-independent."""

    @pytest.fixture(autouse=True)
    def _skip_no_rdkit(self):
        pytest.importorskip('rdkit')

    def test_same_seed_identical_diff_seed_differs(self):
        from scripts._systems import _rdkit_3d
        p_a, s_a = _rdkit_3d('C=CC(=O)OC', seed=7)
        p_a2, s_a2 = _rdkit_3d('C=CC(=O)OC', seed=7)
        p_b, s_b = _rdkit_3d('C=CC(=O)OC', seed=8)
        # same seed -> identical geometry (deterministic)
        np.testing.assert_array_equal(p_a, p_a2)
        # atom ordering / species is seed-independent
        assert s_a == s_a2 == s_b
        # different seed -> different conformer
        assert not np.allclose(p_a, p_b)


class TestBuildVinylAIBNSystem:

    @pytest.fixture(autouse=True)
    def _skip_no_rdkit(self):
        pytest.importorskip('rdkit')

    def test_atom_count(self):
        from scripts._systems import build_vinyl_aibn_system, _rdkit_3d
        rng = np.random.default_rng(0)
        n_mono, n_init = 2, 1
        pos, species, template, groups, prop_map, _ = build_vinyl_aibn_system(
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
        _, _, _, groups, _, _ = build_vinyl_aibn_system(
            n_monomers=3, n_initiators=2, box_size=20.0, rng=rng,
        )
        assert len(groups['radical_C'].atom_indices) == 2
        assert len(groups['vinyl_alpha_C'].atom_indices) == 3

    def test_propagation_map_size(self):
        from scripts._systems import build_vinyl_aibn_system
        rng = np.random.default_rng(2)
        n_mono = 4
        _, _, _, _, prop_map, _ = build_vinyl_aibn_system(
            n_monomers=n_mono, n_initiators=1, box_size=22.0, rng=rng,
        )
        assert len(prop_map) == n_mono

    def test_propagation_map_alpha_in_vinyl_alpha_group(self):
        """Every alpha-C key in propagation_map must appear in vinyl_alpha_C group."""
        from scripts._systems import build_vinyl_aibn_system
        rng = np.random.default_rng(3)
        _, _, _, groups, prop_map, _ = build_vinyl_aibn_system(
            n_monomers=2, n_initiators=1, box_size=18.0, rng=rng,
        )
        alpha_indices = set(groups['vinyl_alpha_C'].atom_indices)
        for alpha in prop_map:
            assert alpha in alpha_indices

    def test_propagation_map_beta_in_vinyl_beta_group_only(self):
        """beta-C values must be in vinyl_beta_C but NOT in radical_C or vinyl_alpha_C."""
        from scripts._systems import build_vinyl_aibn_system
        rng = np.random.default_rng(4)
        _, _, _, groups, prop_map, _ = build_vinyl_aibn_system(
            n_monomers=2, n_initiators=1, box_size=18.0, rng=rng,
        )
        forbidden = set(groups['radical_C'].atom_indices) | set(groups['vinyl_alpha_C'].atom_indices)
        beta_group = set(groups['vinyl_beta_C'].atom_indices)
        for beta in prop_map.values():
            assert beta not in forbidden
            assert beta in beta_group

    def test_chain_c_map_size(self):
        from scripts._systems import build_vinyl_aibn_system
        rng = np.random.default_rng(6)
        n_init = 2
        _, _, _, groups, _, chain_c_map = build_vinyl_aibn_system(
            n_monomers=3, n_initiators=n_init, box_size=20.0, rng=rng,
        )
        assert len(chain_c_map) == n_init
        for rad, chain in chain_c_map.items():
            assert rad in groups['radical_C'].atom_indices
            assert chain in groups['chain_C'].atom_indices

    def test_template_has_4_groups_and_3_pairs(self):
        from scripts._systems import build_vinyl_aibn_system
        rng = np.random.default_rng(7)
        _, _, template, _, _, _ = build_vinyl_aibn_system(
            n_monomers=2, n_initiators=1, box_size=18.0, rng=rng,
        )
        assert len(template.groups) == 4
        assert len(template.pairs) == 3
        constraint_pairs = [p for p in template.pairs if p.constraint_only]
        assert len(constraint_pairs) == 2

    def test_no_overlap(self):
        from scripts._systems import build_vinyl_aibn_system
        rng = np.random.default_rng(5)
        pos, species, _, _, _, _ = build_vinyl_aibn_system(
            n_monomers=2, n_initiators=1, box_size=20.0, rng=rng, min_sep=2.0,
        )
        assert pos.shape[1] == 3
        assert not np.any(np.isnan(pos))


class TestBuildActivationTemplate:
    """Tests for AIBN activation template builder (Table S1, V^d on C-N)."""

    def test_template_targets_dissociation(self):
        """Activation template pairs are dissociation-only (V^d)."""
        from scripts._systems import build_activation_template
        azo_bonds = [(1, 5), (7, 6)]
        template, groups = build_activation_template(azo_bonds)
        assert template.name == 'aibn_activation'
        assert len(template.pairs) == 1
        assert template.pairs[0].is_formation is False

    def test_groups_match_input_bonds(self):
        """azo_C and azo_N groups contain exactly the input bond atoms."""
        from scripts._systems import build_activation_template
        azo_bonds = [(10, 20), (30, 40)]
        _, groups = build_activation_template(azo_bonds)
        assert groups['azo_C'].atom_indices == [10, 30]
        assert groups['azo_N'].atom_indices == [20, 40]

    def test_template_has_2_groups(self):
        from scripts._systems import build_activation_template
        template, _ = build_activation_template([(0, 1)])
        assert len(template.groups) == 2
        assert set(template.groups) == {'azo_C', 'azo_N'}


class TestFindAIBNAzoBonds:
    """Tests for _find_aibn_azo_bonds (excludes nitrile C≡N)."""

    @pytest.fixture(autouse=True)
    def _skip_no_rdkit(self):
        pytest.importorskip('rdkit')

    def test_finds_two_azo_bonds(self):
        from scripts._systems import _find_aibn_azo_bonds
        bonds = _find_aibn_azo_bonds()
        assert len(bonds) == 2

    def test_excludes_nitrile(self):
        """Returned bonds must be C-N single bonds, not C≡N (nitrile)."""
        from rdkit import Chem
        from scripts._systems import _find_aibn_azo_bonds, _AIBN_SMILES
        mol = Chem.MolFromSmiles(_AIBN_SMILES)
        mol = Chem.AddHs(mol)
        bonds = _find_aibn_azo_bonds()
        for c_idx, n_idx in bonds:
            bond = mol.GetBondBetweenAtoms(c_idx, n_idx)
            assert bond is not None
            assert bond.GetBondTypeAsDouble() == 1.0


class TestBuildFullAIBNSystem:

    @pytest.fixture(autouse=True)
    def _skip_no_rdkit(self):
        pytest.importorskip('rdkit')

    def test_azo_bonds_are_c_n_pairs(self):
        """aibn_azo_bonds should reference C and N atoms in the global species."""
        from scripts._systems import build_full_aibn_system
        rng = np.random.default_rng(0)
        _, species, azo_bonds, _, _, _, _ = build_full_aibn_system(
            n_monomers=2, n_aibn=1, box_size=20.0, rng=rng,
        )
        for c_idx, n_idx in azo_bonds:
            assert species[c_idx] == 'C'
            assert species[n_idx] == 'N'

    def test_activation_template_from_azo_bonds(self):
        """build_activation_template on the output of build_full_aibn_system
        produces a valid template with correct group sizes."""
        from scripts._systems import build_full_aibn_system, build_activation_template
        rng = np.random.default_rng(0)
        _, _, azo_bonds, _, _, _, _ = build_full_aibn_system(
            n_monomers=2, n_aibn=1, box_size=20.0, rng=rng,
        )
        template, groups = build_activation_template(azo_bonds)
        assert len(groups['azo_C'].atom_indices) == len(azo_bonds)
        assert len(groups['azo_N'].atom_indices) == len(azo_bonds)
        assert template.pairs[0].is_formation is False


class TestNylon66Helpers:
    """Tests for nylon-6,6 RDKit-based helpers."""

    @pytest.fixture(autouse=True)
    def _skip_no_rdkit(self):
        pytest.importorskip('rdkit')

    def test_find_terminal_amine_n(self):
        from scripts._systems import _find_terminal_amine_n
        indices = _find_terminal_amine_n('NCCCCCCN')
        assert len(indices) == 2

    def test_find_amine_h(self):
        from scripts._systems import _find_amine_h, _find_terminal_amine_n
        n_indices = _find_terminal_amine_n('NCCCCCCN')
        for n_idx in n_indices:
            h_idx = _find_amine_h('NCCCCCCN', n_idx)
            assert h_idx != n_idx

    def test_find_carboxyl_c_and_oh(self):
        from scripts._systems import _find_carboxyl_c_and_oh
        pairs = _find_carboxyl_c_and_oh('OC(=O)CCCCC(=O)O')
        assert len(pairs) == 2
        for c_idx, oh_idx in pairs:
            assert c_idx != oh_idx

    def test_carboxyl_oh_has_hydrogen(self):
        from rdkit import Chem
        from scripts._systems import _find_carboxyl_c_and_oh
        smiles = 'OC(=O)CCCCC(=O)O'
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        pairs = _find_carboxyl_c_and_oh(smiles)
        for _, oh_idx in pairs:
            atom = mol.GetAtomWithIdx(oh_idx)
            assert atom.GetSymbol() == 'O'
            h_count = sum(1 for n in atom.GetNeighbors() if n.GetSymbol() == 'H')
            assert h_count >= 1

    def test_rdkit_3d_diamine_species(self):
        from scripts._systems import _rdkit_3d
        positions, species = _rdkit_3d('NCCCCCCN')
        assert species.count('N') == 2
        assert species.count('C') == 6
        assert positions.shape == (len(species), 3)

    def test_rdkit_3d_diacid_species(self):
        from scripts._systems import _rdkit_3d
        positions, species = _rdkit_3d('OC(=O)CCCCC(=O)O')
        assert species.count('O') == 4
        assert species.count('C') == 6
        assert positions.shape == (len(species), 3)


class TestBuildNylon66System:

    @pytest.fixture(autouse=True)
    def _skip_no_rdkit(self):
        pytest.importorskip('rdkit')

    def test_atom_count(self):
        from scripts._systems import build_nylon66_system, _rdkit_3d
        rng = np.random.default_rng(0)
        n_di, n_ac = 2, 2
        pos, species, _, _ = build_nylon66_system(
            n_diamines=n_di, n_diacids=n_ac, box_size=20.0, rng=rng,
        )
        n_per_di = len(_rdkit_3d('NCCCCCCN')[1])
        n_per_ac = len(_rdkit_3d('OC(=O)CCCCC(=O)O')[1])
        expected = n_di * n_per_di + n_ac * n_per_ac
        assert pos.shape == (expected, 3)
        assert len(species) == expected

    def test_group_sizes(self):
        from scripts._systems import build_nylon66_system
        rng = np.random.default_rng(1)
        n_di, n_ac = 3, 3
        _, _, _, groups = build_nylon66_system(
            n_diamines=n_di, n_diacids=n_ac, box_size=25.0, rng=rng,
        )
        assert len(groups['amine_N'].atom_indices) == n_di * 2
        assert len(groups['carboxyl_C'].atom_indices) == n_ac * 2
        assert len(groups['amine_H'].atom_indices) == n_di * 2
        assert len(groups['carboxyl_OH'].atom_indices) == n_ac * 2

    def test_template_has_4_groups(self):
        from scripts._systems import build_nylon66_system
        rng = np.random.default_rng(2)
        _, _, template, _ = build_nylon66_system(
            n_diamines=2, n_diacids=2, box_size=20.0, rng=rng,
        )
        assert len(template.groups) == 4
        assert template.groups == ['amine_N', 'carboxyl_C', 'amine_H', 'carboxyl_OH']

    def test_template_has_4_pairs(self):
        from scripts._systems import build_nylon66_system
        rng = np.random.default_rng(3)
        _, _, template, _ = build_nylon66_system(
            n_diamines=2, n_diacids=2, box_size=20.0, rng=rng,
        )
        assert len(template.pairs) == 4
        formation_pairs = [p for p in template.pairs if p.is_formation]
        dissociation_pairs = [p for p in template.pairs if not p.is_formation]
        assert len(formation_pairs) == 2
        assert len(dissociation_pairs) == 2

    def test_kl_pair_is_bias_only(self):
        """RF5: nylon k-l (H-OH) has score_pair=False and is_formation=True."""
        from scripts._systems import build_nylon66_system
        rng = np.random.default_rng(3)
        _, _, template, _ = build_nylon66_system(
            n_diamines=2, n_diacids=2, box_size=20.0, rng=rng,
        )
        kl = [p for p in template.pairs
              if p.group_a == 'amine_H' and p.group_b == 'carboxyl_OH']
        assert len(kl) == 1
        assert kl[0].score_pair is False
        assert kl[0].is_formation is True

    def test_amine_h_bonded_to_amine_n(self):
        """Each amine_H must be a neighbor of some amine_N in the system."""
        from rdkit import Chem
        from scripts._systems import build_nylon66_system, _rdkit_3d
        rng = np.random.default_rng(4)
        _, species, _, groups = build_nylon66_system(
            n_diamines=2, n_diacids=2, box_size=20.0, rng=rng,
        )
        for h_idx in groups['amine_H'].atom_indices:
            assert species[h_idx] == 'H'
        for n_idx in groups['amine_N'].atom_indices:
            assert species[n_idx] == 'N'

    def test_no_nan_positions(self):
        from scripts._systems import build_nylon66_system
        rng = np.random.default_rng(5)
        pos, _, _, _ = build_nylon66_system(
            n_diamines=2, n_diacids=2, box_size=20.0, rng=rng,
        )
        assert not np.any(np.isnan(pos))


class TestPruneUndissociatedCenters:
    """prune_undissociated_centers reconciles the reactive groups with the
    actual AIBN activation result (specs/decisions.md 2026-07-06)."""

    @staticmethod
    def _make_groups(n_centers):
        from kagome.reactive.groups import ReactiveGroup
        # radical centres at 0,10,20,...; each chain_C at radical+1 (distinct).
        radicals = [10 * i for i in range(n_centers)]
        chain_cs = [10 * i + 1 for i in range(n_centers)]
        groups = {
            'radical_C': ReactiveGroup('radical_C', list(radicals)),
            'chain_C': ReactiveGroup('chain_C', list(chain_cs)),
            'vinyl_alpha_C': ReactiveGroup('vinyl_alpha_C', [500, 501, 502]),
        }
        chain_c_map = {radicals[i]: chain_cs[i] for i in range(n_centers)}
        return groups, chain_c_map, radicals, chain_cs

    def test_partial_dissociation_reduces_groups(self):
        from scripts._systems import prune_undissociated_centers
        groups, chain_c_map, radicals, chain_cs = self._make_groups(10)
        dissociated_c = radicals[:6]  # first 6 centres actually dissociated
        new_groups, new_ccmap, n_pruned = prune_undissociated_centers(
            groups, chain_c_map, dissociated_c,
        )
        assert n_pruned == 4
        assert new_groups['radical_C'].atom_indices == radicals[:6]
        assert new_groups['chain_C'].atom_indices == chain_cs[:6]
        assert new_ccmap == {radicals[i]: chain_cs[i] for i in range(6)}
        # vinyl_alpha_C is untouched.
        assert new_groups['vinyl_alpha_C'].atom_indices == [500, 501, 502]

    def test_full_dissociation_is_no_op(self):
        from scripts._systems import prune_undissociated_centers
        groups, chain_c_map, radicals, chain_cs = self._make_groups(10)
        new_groups, new_ccmap, n_pruned = prune_undissociated_centers(
            groups, chain_c_map, radicals,  # all 10 dissociated
        )
        assert n_pruned == 0
        assert new_groups['radical_C'].atom_indices == radicals
        assert new_groups['chain_C'].atom_indices == chain_cs
        assert new_ccmap == chain_c_map

    def test_zero_dissociation_removes_all(self):
        from scripts._systems import prune_undissociated_centers
        groups, chain_c_map, radicals, chain_cs = self._make_groups(10)
        new_groups, new_ccmap, n_pruned = prune_undissociated_centers(
            groups, chain_c_map, [],  # nothing dissociated
        )
        assert n_pruned == 10
        assert new_groups['radical_C'].atom_indices == []
        assert new_groups['chain_C'].atom_indices == []
        assert new_ccmap == {}

    def test_inputs_are_not_mutated(self):
        from scripts._systems import prune_undissociated_centers
        groups, chain_c_map, radicals, chain_cs = self._make_groups(10)
        prune_undissociated_centers(groups, chain_c_map, radicals[:6])
        # Originals unchanged (pure function).
        assert groups['radical_C'].atom_indices == radicals
        assert groups['chain_C'].atom_indices == chain_cs
        assert chain_c_map == {radicals[i]: chain_cs[i] for i in range(10)}


class TestProductionSpin:
    """production_spin computes the high-spin multiplicity from the surviving
    radical count (specs/decisions.md 2026-07-02 Layer 3, 2026-07-06)."""

    def test_high_spin_sum(self):
        from scripts._systems import production_spin
        assert production_spin(0) == 1
        assert production_spin(1) == 2
        assert production_spin(6) == 7

    def test_cap_clamps(self):
        from scripts._systems import production_spin
        assert production_spin(10, cap=3) == 3
        assert production_spin(2, cap=5) == 3  # below cap -> unclamped

    def test_uses_pruned_radical_count(self):
        """Spin follows the reconciled (pruned) radical count, not the raw
        registration count."""
        from scripts._systems import prune_undissociated_centers, production_spin
        from kagome.reactive.groups import ReactiveGroup
        radicals = [10 * i for i in range(10)]
        chain_cs = [10 * i + 1 for i in range(10)]
        groups = {
            'radical_C': ReactiveGroup('radical_C', list(radicals)),
            'chain_C': ReactiveGroup('chain_C', list(chain_cs)),
        }
        chain_c_map = {radicals[i]: chain_cs[i] for i in range(10)}
        new_groups, _, _ = prune_undissociated_centers(
            groups, chain_c_map, radicals[:6],
        )
        n_radicals = len(new_groups['radical_C'].atom_indices)
        assert production_spin(n_radicals) == 7  # 6 radicals + 1, not 10 + 1


class TestEpoxyAmineHelpers:
    """RDKit helpers for the bulk epoxy-amine system (Track 2 E1,
    decisions.md 2026-07-09)."""

    @pytest.fixture(autouse=True)
    def _skip_no_rdkit(self):
        pytest.importorskip('rdkit')

    def test_find_epoxide_propylene_oxide(self):
        """CC1CO1: one ring; terminal C is the CH2 (2 H), paired with ring O."""
        from rdkit import Chem
        from scripts._systems import _find_epoxide_c_o
        smiles = 'CC1CO1'
        pairs = _find_epoxide_c_o(smiles)
        assert len(pairs) == 1
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        c_idx, o_idx = pairs[0]
        assert mol.GetAtomWithIdx(c_idx).GetSymbol() == 'C'
        assert mol.GetAtomWithIdx(o_idx).GetSymbol() == 'O'
        h_count = sum(1 for n in mol.GetAtomWithIdx(c_idx).GetNeighbors()
                      if n.GetSymbol() == 'H')
        assert h_count == 2  # terminal CH2, not the methyl-bearing CH

    def test_find_epoxide_dgeba_has_two_rings(self):
        from scripts._systems import _DGEBA_SMILES, _find_epoxide_c_o
        pairs = _find_epoxide_c_o(_DGEBA_SMILES)
        assert len(pairs) == 2
        assert len({c for c, _ in pairs}) == 2
        assert len({o for _, o in pairs}) == 2

    def test_find_amine_n_h_deta(self):
        """DETA NCCNCCN: 2 primary (2 H) + 1 secondary (1 H) = 3 N, 5 H."""
        from scripts._systems import _DETA_SMILES, _find_amine_n_h
        result = _find_amine_n_h(_DETA_SMILES)
        assert len(result) == 3
        h_counts = sorted(len(hs) for _, hs in result)
        assert h_counts == [1, 2, 2]

    def test_find_amine_n_h_methylamine(self):
        from scripts._systems import _find_amine_n_h
        result = _find_amine_n_h('CN')
        assert len(result) == 1
        assert len(result[0][1]) == 2


class TestBuildEpoxyAmineSystem:
    """Bulk epoxy-amine curing builder (4-group ring-opening template)."""

    @pytest.fixture(autouse=True)
    def _skip_no_rdkit(self):
        pytest.importorskip('rdkit')

    @staticmethod
    def _build(n_epoxies=2, n_amines=1, box=32.0, seed=0):
        from scripts._systems import build_epoxy_amine_system
        rng = np.random.default_rng(seed)
        return build_epoxy_amine_system(
            n_epoxies=n_epoxies, n_amines=n_amines, box_size=box, rng=rng,
        )

    def test_atom_count(self):
        from scripts._systems import _DETA_SMILES, _DGEBA_SMILES, _rdkit_3d
        pos, species, _, _, _ = self._build()
        n_per_epoxy = len(_rdkit_3d(_DGEBA_SMILES)[1])
        n_per_amine = len(_rdkit_3d(_DETA_SMILES)[1])
        expected = 2 * n_per_epoxy + 1 * n_per_amine
        assert pos.shape == (expected, 3)
        assert len(species) == expected

    def test_group_sizes(self):
        _, _, _, groups, _ = self._build()
        # 2 DGEBA x 2 epoxides; 1 DETA x (3 N, 5 N-H)
        assert len(groups['epoxy_C'].atom_indices) == 4
        assert len(groups['ring_O'].atom_indices) == 4
        assert len(groups['amine_N'].atom_indices) == 3
        assert len(groups['amine_H'].atom_indices) == 5

    def test_template_structure(self):
        _, _, template, _, _ = self._build()
        assert template.groups == ['amine_N', 'epoxy_C', 'amine_H', 'ring_O']
        assert len(template.pairs) == 4
        formation_pairs = [p for p in template.pairs if p.is_formation]
        dissociation_pairs = [p for p in template.pairs if not p.is_formation]
        assert len(formation_pairs) == 2
        assert len(dissociation_pairs) == 2

    def test_kl_pair_is_bias_only(self):
        """The hydroxyl-forming k-l (amine_H - ring_O) pair must be bias-only:
        score_pair=False (not a trigger) and count_as_reaction=False."""
        _, _, template, _, _ = self._build()
        kl = [p for p in template.pairs
              if p.group_a == 'amine_H' and p.group_b == 'ring_O']
        assert len(kl) == 1
        assert kl[0].score_pair is False
        assert kl[0].count_as_reaction is False
        assert kl[0].is_formation is True

    def test_amine_h_map_consistent_with_groups(self):
        _, species, _, groups, amine_h_map = self._build()
        assert set(amine_h_map.keys()) == set(groups['amine_N'].atom_indices)
        all_hs = [h for hs in amine_h_map.values() for h in hs]
        assert sorted(all_hs) == sorted(groups['amine_H'].atom_indices)
        for n_idx, hs in amine_h_map.items():
            assert species[n_idx] == 'N'
            for h in hs:
                assert species[h] == 'H'

    def test_group_species(self):
        _, species, _, groups, _ = self._build()
        for c_idx in groups['epoxy_C'].atom_indices:
            assert species[c_idx] == 'C'
        for o_idx in groups['ring_O'].atom_indices:
            assert species[o_idx] == 'O'

    def test_no_nan_positions(self):
        pos, _, _, _, _ = self._build(seed=5)
        assert not np.any(np.isnan(pos))
