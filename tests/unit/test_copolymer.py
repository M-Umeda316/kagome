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


class TestCopolymerAlphaSpecies:
    """copolymer_alpha_species must reconstruct the SAME alpha indices as the
    builder, so the reactivity analysis labels the right species."""

    def test_matches_builder_alpha_indices(self):
        from scripts._systems import copolymer_alpha_species
        specs = [(ACRYLATE, 3), (METHACRYLATE, 2)]
        _, _, _, groups, _, _ = _build(3, 2, 1)
        mapping = copolymer_alpha_species(specs, n_initiators=1)
        # every alpha index in the builder's group is mapped, and vice versa
        assert set(mapping) == set(groups['vinyl_alpha_C'].atom_indices)

    def test_species_labels_partition_correctly(self):
        from scripts._systems import copolymer_alpha_species
        specs = [(ACRYLATE, 3), (METHACRYLATE, 2)]
        mapping = copolymer_alpha_species(specs, n_initiators=1)
        n_acr = sum(1 for s in mapping.values() if s == ACRYLATE)
        n_mac = sum(1 for s in mapping.values() if s == METHACRYLATE)
        assert n_acr == 3
        assert n_mac == 2

    def test_acrylate_indices_precede_methacrylate(self):
        from scripts._systems import copolymer_alpha_species
        mapping = copolymer_alpha_species([(ACRYLATE, 2), (METHACRYLATE, 2)],
                                          n_initiators=1)
        acr = sorted(i for i, s in mapping.items() if s == ACRYLATE)
        mac = sorted(i for i, s in mapping.items() if s == METHACRYLATE)
        assert max(acr) < min(mac)  # placement order preserved


class TestCopolymerAtomSpecies:
    """copolymer_atom_species must label every atom with its block, on the
    same offsets as the builder, so radical endpoints classify correctly."""

    def test_covers_all_atoms_and_agrees_with_alpha_map(self):
        from scripts._systems import (
            copolymer_alpha_species,
            copolymer_atom_species,
        )
        specs = [(ACRYLATE, 3), (METHACRYLATE, 2)]
        pos, _, _, _, _, _ = _build(3, 2, 1)
        atom_map = copolymer_atom_species(specs, n_initiators=1)
        assert set(atom_map) == set(range(pos.shape[0]))
        # every alpha-C's block label matches the alpha map's species label
        alpha_map = copolymer_alpha_species(specs, n_initiators=1)
        for idx, smi in alpha_map.items():
            assert atom_map[idx] == smi

    def test_initiator_atoms_labeled_initiator(self):
        from scripts._systems import _rdkit_3d, copolymer_atom_species
        atom_map = copolymer_atom_species(
            [(ACRYLATE, 2), (METHACRYLATE, 2)], n_initiators=2)
        n_init_atoms = 2 * len(_rdkit_3d('CC(C)C#N')[1])
        for i in range(n_init_atoms):
            assert atom_map[i] == 'CC(C)C#N'
        assert atom_map[n_init_atoms] == ACRYLATE


class TestCopolymerInitialBonds:
    """copolymer_initial_bonds (WM-P1, specs/decisions.md 2026-07-17 前提工事)
    must line up 1:1 with build_vinyl_copolymer_system's running-offset layout,
    cover every fragment with no cross-fragment bonds, carry the alpha=beta
    C=C double bonds, be valence-clean, and support apply_vinyl_addition for
    BOTH monomer species (including the initiator's placeholder-H shed)."""

    _SPECS = [(ACRYLATE, 3), (METHACRYLATE, 2)]

    def test_bond_count_and_fragment_coverage(self):
        """Returned bonds equal the union of per-fragment RDKit bonds shifted
        by the same offsets as copolymer_atom_species: every bond's endpoints
        map to the SAME fragment, and no cross-fragment bonds exist."""
        from scripts._systems import (
            _rdkit_local_bonds,
            copolymer_atom_species,
            copolymer_initial_bonds,
        )
        atom_map = copolymer_atom_species(self._SPECS, n_initiators=1)
        bonds = copolymer_initial_bonds(self._SPECS, n_initiators=1)
        assert bonds

        n_init_bonds = len(_rdkit_local_bonds('CC(C)C#N'))
        n_acr_bonds = len(_rdkit_local_bonds(ACRYLATE))
        n_mac_bonds = len(_rdkit_local_bonds(METHACRYLATE))
        expected_count = 1 * n_init_bonds + 3 * n_acr_bonds + 2 * n_mac_bonds
        assert len(bonds) == expected_count

        for i, j, order in bonds:
            assert atom_map[i] == atom_map[j]  # no cross-fragment bond
            assert order > 0.0

    def test_alpha_beta_double_bond_present_for_every_monomer(self):
        from kagome.reactive.topology import BondTopology
        from scripts._systems import copolymer_initial_bonds

        _, _, _, _, pmap, _ = _build(3, 2, 1)
        bonds = copolymer_initial_bonds(self._SPECS, n_initiators=1)
        topo = BondTopology.from_bonds(bonds)
        for alpha, beta in pmap.items():
            assert topo.has_bond(alpha, beta)
            assert topo.order(alpha, beta) == 2.0

    def test_no_over_coordinated_atoms_in_initial_topology(self):
        from kagome.reactive.topology import BondTopology, over_coordinated_atoms
        from scripts._systems import copolymer_initial_bonds

        _, species, _, _, _, _ = _build(3, 2, 1)
        bonds = copolymer_initial_bonds(self._SPECS, n_initiators=1)
        topo = BondTopology.from_bonds(bonds)
        assert over_coordinated_atoms(topo, species) == []

    def test_reaction_round_trip_for_both_species(self):
        """apply_vinyl_addition for one initiator radical + one ACRYLATE alpha,
        and separately one METHACRYLATE alpha: new sigma bond, C=C order drops
        to 1.0, no over-coordination — including the initiator placeholder-H
        shed (the radical C is 4-coordinate closed-shell before addition)."""
        from kagome.reactive.topology import (
            BondTopology,
            apply_vinyl_addition,
            over_coordinated_atoms,
        )
        from scripts._systems import copolymer_alpha_species, copolymer_initial_bonds

        _, species, _, groups, pmap, _ = _build(3, 2, 1)
        bonds = copolymer_initial_bonds(self._SPECS, n_initiators=1)
        alpha_species = copolymer_alpha_species(self._SPECS, n_initiators=1)
        radical_c = groups['radical_C'].atom_indices[0]

        acr_alpha = next(i for i, s in alpha_species.items() if s == ACRYLATE)
        mac_alpha = next(i for i, s in alpha_species.items() if s == METHACRYLATE)

        for alpha in (acr_alpha, mac_alpha):
            topo = BondTopology.from_bonds(bonds)
            beta = pmap[alpha]
            assert topo.coordination_number(radical_c) == 4  # closed-shell + placeholder H
            apply_vinyl_addition(topo, radical_c, alpha, pmap, species)
            assert topo.has_bond(radical_c, alpha)            # new sigma bond
            assert topo.order(alpha, beta) == 1.0              # C=C -> C-C
            assert topo.coordination_number(radical_c) == 4    # placeholder H shed, not 5
            assert over_coordinated_atoms(topo, species) == []


class TestCrossPropagationAnalysis:
    """analyze() classifies each confirmed formation as
    (terminal species -> incorporated species) from the radical endpoint."""

    def _fake_run_dir(self, tmp_path, events):
        import json
        (tmp_path / 'summary.json').write_text(json.dumps({
            'n_acrylate': 3, 'n_methacrylate': 2, 'n_initiators': 1,
        }), encoding='utf-8')
        lines = [json.dumps(ev) for ev in events]
        (tmp_path / 'bonds.jsonl').write_text('\n'.join(lines), encoding='utf-8')
        return tmp_path

    def test_initiator_then_chain_terminal(self, tmp_path):
        from scripts._systems import (
            _rdkit_3d,
            copolymer_alpha_species,
            copolymer_atom_species,
        )
        from scripts.analyze_copolymer_reactivity import analyze

        specs = [(ACRYLATE, 3), (METHACRYLATE, 2)]
        alpha = copolymer_alpha_species(specs, n_initiators=1)
        acr_alphas = sorted(i for i, s in alpha.items() if s == ACRYLATE)
        mac_alphas = sorted(i for i, s in alpha.items() if s == METHACRYLATE)
        n_init_atoms = len(_rdkit_3d('CC(C)C#N')[1])

        # event 1: initiator radical (atom 1, inside the initiator block)
        # adds the first acrylate; event 2: the radical now sits on that
        # acrylate's beta-C (same block as its alpha) and adds a methacrylate.
        first_acr_alpha = acr_alphas[0]
        assert first_acr_alpha >= n_init_atoms
        events = [
            {'event_type': 'confirmed_formation', 'cycle': 0,
             'atom_a': 1, 'atom_b': first_acr_alpha},
            {'event_type': 'confirmed_formation', 'cycle': 1,
             'atom_a': mac_alphas[0], 'atom_b': first_acr_alpha + 1},
            {'event_type': 'tentative_formation', 'cycle': 2,
             'atom_a': 1, 'atom_b': acr_alphas[1]},  # must be ignored
        ]
        result = analyze(self._fake_run_dir(tmp_path, events))

        cross = result['cross_propagation']
        assert cross['initiator'] == {'acrylate': 1, 'methacrylate': 0}
        assert cross['acrylate'] == {'acrylate': 0, 'methacrylate': 1}
        assert cross['methacrylate'] == {'acrylate': 0, 'methacrylate': 0}
        assert result['total_incorporations'] == 2
        est = result['reactivity_ratio_estimates']
        assert est['monomer_terminal_events'] == 1
        # acr->mac exists but acr->acr is 0: point estimate degenerates to 0.0
        assert est['r_acrylate'] == 0.0
        assert est['r_methacrylate'] is None  # no mac->acr denominator event

    def test_reactivity_ratios_from_full_2x2(self, tmp_path):
        from scripts._systems import copolymer_alpha_species
        from scripts.analyze_copolymer_reactivity import analyze

        specs = [(ACRYLATE, 3), (METHACRYLATE, 2)]
        alpha = copolymer_alpha_species(specs, n_initiators=1)
        acr = sorted(i for i, s in alpha.items() if s == ACRYLATE)
        mac = sorted(i for i, s in alpha.items() if s == METHACRYLATE)

        # radical endpoints chosen INSIDE monomer blocks (alpha+1 = same block):
        # acr*->acr x2, acr*->mac x1, mac*->mac x2, mac*->acr x1
        events = [
            {'event_type': 'confirmed_formation', 'cycle': 0,
             'atom_a': acr[0] + 1, 'atom_b': acr[1]},
            {'event_type': 'confirmed_formation', 'cycle': 1,
             'atom_a': acr[1] + 1, 'atom_b': acr[2]},
            {'event_type': 'confirmed_formation', 'cycle': 2,
             'atom_a': acr[2] + 1, 'atom_b': mac[0]},
            {'event_type': 'confirmed_formation', 'cycle': 3,
             'atom_a': mac[0] + 1, 'atom_b': mac[1]},
            {'event_type': 'confirmed_formation', 'cycle': 4,
             'atom_a': mac[1] + 1, 'atom_b': mac[1]},
            {'event_type': 'confirmed_formation', 'cycle': 5,
             'atom_a': mac[1] + 2, 'atom_b': acr[0]},
        ]
        # NOTE: event 5 reuses mac[1] as alpha; fine — analyze() does not
        # dedupe, it just classifies endpoints.
        result = analyze(self._fake_run_dir(tmp_path, events))
        est = result['reactivity_ratio_estimates']
        # r_acr = (2/1) * (n_mac/n_acr = 2/3), r_mac = (2/1) * (3/2)
        assert est['r_acrylate'] == pytest.approx(2 * 2 / 3)
        assert est['r_methacrylate'] == pytest.approx(2 * 3 / 2)
        assert est['monomer_terminal_events'] == 6


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
