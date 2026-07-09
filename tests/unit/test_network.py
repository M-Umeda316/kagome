"""Tests for network analysis (epoxy-amine curing species / gel point)."""
import json
from pathlib import Path

import pytest

from kagome.analysis.carothers import monomer_sets_from_bonds
from kagome.analysis.network import (
    amine_ranks,
    count_hydroxyls,
    crosslink_counts,
    find_epoxide_rings,
    gel_point_flory_stockmayer,
    largest_component_fraction,
    load_topology_snapshots,
    species_series,
)


class TestFindEpoxideRings:

    def test_triangle_with_decoys(self):
        # 0-1-2: closed C-C-O triangle (the only epoxide).
        # 3-4-5: OPEN C-C-O chain (no 3-5 bond) -> not a ring.
        # 6-7-8: C-O-C ether (no C-C bond) -> not a ring.
        species = ['C', 'C', 'O', 'C', 'C', 'O', 'C', 'O', 'C']
        bonds = [
            (0, 1), (0, 2), (1, 2),   # closed epoxide
            (3, 4), (4, 5),           # open chain
            (6, 7), (7, 8),           # ether
        ]
        assert find_epoxide_rings(bonds, species) == [(0, 1, 2)]

    def test_ccc_triangle_species_mismatch(self):
        # Cyclopropane-like C-C-C triangle must not match {C, C, O}.
        species = ['C', 'C', 'C']
        bonds = [(0, 1), (0, 2), (1, 2)]
        assert find_epoxide_rings(bonds, species) == []

    def test_c1_less_than_c2_and_sorted(self):
        # Two rings; carbons returned as (min, max) and list sorted.
        species = ['O', 'C', 'C', 'C', 'C', 'O']
        bonds = [(4, 3), (3, 5), (4, 5), (2, 1), (1, 0), (2, 0)]
        assert find_epoxide_rings(bonds, species) == [(1, 2, 0), (3, 4, 5)]

    def test_bond_order_triples_accepted(self):
        species = ['C', 'C', 'O']
        bonds = [(0, 1, 1.0), (0, 2, 1.0), (1, 2, 1.0)]
        assert find_epoxide_rings(bonds, species) == [(0, 1, 2)]


class TestAmineRanks:

    def test_primary_secondary_tertiary(self):
        # N0: 2 H (primary), N4: 1 H (secondary), N8: 0 H (tertiary).
        species = ['N', 'H', 'H', 'C',
                   'N', 'H', 'C', 'C',
                   'N', 'C', 'C', 'C']
        bonds = [
            (0, 1), (0, 2), (0, 3),
            (4, 5), (4, 6), (4, 7),
            (8, 9), (8, 10), (8, 11),
        ]
        assert amine_ranks(bonds, species) == {0: 2, 4: 1, 8: 0}

    def test_bondless_n_reports_zero(self):
        # Every N in the species list is reported, even with no bonds at all.
        species = ['N', 'C']
        assert amine_ranks([], species) == {0: 0}

    def test_no_nitrogen(self):
        assert amine_ranks([(0, 1)], ['C', 'H']) == {}


class TestCountHydroxyls:

    def test_hydroxyl_on_carbon_counted(self):
        # C-O-H: one C neighbour, one H neighbour -> counted.
        species = ['C', 'O', 'H']
        bonds = [(0, 1), (1, 2)]
        assert count_hydroxyls(bonds, species) == 1

    def test_ether_not_counted(self):
        # C-O-C: two C, zero H.
        species = ['C', 'O', 'C']
        bonds = [(0, 1), (1, 2)]
        assert count_hydroxyls(bonds, species) == 0

    def test_closed_epoxide_ring_o_not_counted(self):
        # Ring O has two C neighbours.
        species = ['C', 'C', 'O']
        bonds = [(0, 1), (0, 2), (1, 2)]
        assert count_hydroxyls(bonds, species) == 0

    def test_water_not_counted(self):
        # Water O has 0 C and 2 H -> fails the one-C-one-H rule. (No water
        # exists in the epoxy-amine addition, but the rule must be robust.)
        species = ['H', 'O', 'H']
        bonds = [(0, 1), (1, 2)]
        assert count_hydroxyls(bonds, species) == 0


class TestSpeciesSeries:

    def test_ring_opening_species_balance(self):
        # Minimal ring-opening addition: epoxide (C0, C1, O2) + primary amine
        # (N3 with H4, H5, backbone C6).  The event removes the terminal ring
        # C-O bond (0-2), forms N-C (3-0), transfers H4 from N to O
        # (remove 3-4, add 2-4).  Balance: epoxide -1, primary -1,
        # secondary +1, hydroxyl +1.
        species = ['C', 'C', 'O', 'N', 'H', 'H', 'C']
        bonds_before = [(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (3, 6)]
        bonds_after = [(0, 1), (1, 2), (0, 3), (3, 5), (3, 6), (2, 4)]
        snapshots = [(0, -1, bonds_before), (100, 0, bonds_after)]

        series = species_series(snapshots, species)

        assert series[0] == {
            'step': 0, 'cycle': -1, 'n_epoxide': 1, 'n_amine_primary': 1,
            'n_amine_secondary': 0, 'n_amine_tertiary': 0, 'n_hydroxyl': 0,
        }
        assert series[1] == {
            'step': 100, 'cycle': 0, 'n_epoxide': 0, 'n_amine_primary': 0,
            'n_amine_secondary': 1, 'n_amine_tertiary': 0, 'n_hydroxyl': 1,
        }


class TestGelPointFloryStockmayer:

    def test_dgeba_deta(self):
        # DGEBA f=2, DETA g=5, r=1 -> 1/sqrt(1*1*4) = 0.5.
        assert gel_point_flory_stockmayer(2, 5, 1.0) == pytest.approx(0.5)

    def test_f2_g4(self):
        assert gel_point_flory_stockmayer(2, 4, 1.0) == pytest.approx(
            1.0 / 3.0 ** 0.5)

    def test_imbalance_ratio(self):
        assert gel_point_flory_stockmayer(2, 5, 0.5) == pytest.approx(
            1.0 / 2.0 ** 0.5)

    def test_invalid_functionality_raises(self):
        with pytest.raises(ValueError):
            gel_point_flory_stockmayer(1, 5)
        with pytest.raises(ValueError):
            gel_point_flory_stockmayer(2, 1)

    def test_invalid_ratio_raises(self):
        with pytest.raises(ValueError):
            gel_point_flory_stockmayer(2, 5, 0.0)
        with pytest.raises(ValueError):
            gel_point_flory_stockmayer(2, 5, -1.0)


class TestLargestComponentFraction:

    def test_one_bridge_gives_half(self):
        # 4 monomers of 2 atoms; one inter-monomer bond (1-2) bridges two of
        # them -> largest component 2/4 = 0.5.
        monomer_sets = [[0, 1], [2, 3], [4, 5], [6, 7]]
        bonds = [(0, 1), (2, 3), (4, 5), (6, 7), (1, 2)]
        assert largest_component_fraction(bonds, monomer_sets) == pytest.approx(0.5)

    def test_all_disconnected(self):
        monomer_sets = [[0, 1], [2, 3], [4, 5], [6, 7]]
        bonds = [(0, 1), (2, 3), (4, 5), (6, 7)]
        assert largest_component_fraction(bonds, monomer_sets) == pytest.approx(0.25)

    def test_fully_connected(self):
        monomer_sets = [[0, 1], [2, 3]]
        bonds = [(0, 1), (2, 3), (1, 2)]
        assert largest_component_fraction(bonds, monomer_sets) == pytest.approx(1.0)

    def test_no_monomers(self):
        assert largest_component_fraction([], []) == pytest.approx(0.0)


class TestCrosslinkCounts:

    def test_synthetic_two_monomer_system(self):
        # Monomer 0: one epoxide ring (0, 1, 2), fully opened in the current
        # bonds.  Monomer 1: amine N3 with both H transferred -> tertiary.
        species = ['C', 'C', 'O', 'N', 'H', 'H', 'C']
        monomer_sets = [[0, 1, 2], [3, 4, 5, 6]]
        rings_initial = [(0, 1, 2)]
        bonds = [(1, 2), (0, 3), (3, 6), (2, 4)]  # ring open, N has 0 H

        counts = crosslink_counts(bonds, species, rings_initial, monomer_sets)

        assert counts == {
            'n_tertiary_amine': 1,
            'n_fully_reacted_epoxy_monomers': 1,
            'n_monomers': 2,
        }

    def test_partially_reacted_monomer_not_counted(self):
        # Monomer 0 has TWO initial rings; only one is open -> not fully
        # reacted.  Ring A (0, 1, 2) open, ring B (3, 4, 5) still closed.
        species = ['C', 'C', 'O', 'C', 'C', 'O']
        monomer_sets = [[0, 1, 2, 3, 4, 5]]
        rings_initial = [(0, 1, 2), (3, 4, 5)]
        bonds = [(1, 2), (3, 4), (3, 5), (4, 5)]  # A open, B closed

        counts = crosslink_counts(bonds, species, rings_initial, monomer_sets)

        assert counts['n_fully_reacted_epoxy_monomers'] == 0
        assert counts['n_tertiary_amine'] == 0
        assert counts['n_monomers'] == 1


class TestLoadTopologySnapshots:

    def test_load_and_sort(self, tmp_path):
        path = tmp_path / 'topology.jsonl'
        rec_late = {'step': 100, 'cycle': 0, 'n_bonds': 1, 'bonds': [[0, 2, 1.0]]}
        rec_init = {'step': 0, 'cycle': -1, 'n_bonds': 1, 'bonds': [[0, 1, 1.0]]}
        path.write_text(
            json.dumps(rec_late) + '\n' + json.dumps(rec_init) + '\n',
            encoding='utf-8',
        )

        snapshots = load_topology_snapshots(path)

        assert snapshots == [
            (0, -1, [(0, 1, 1.0)]),
            (100, 0, [(0, 2, 1.0)]),
        ]

    def test_missing_file_returns_empty(self, tmp_path):
        assert load_topology_snapshots(tmp_path / 'nope.jsonl') == []


# ── Real-data regression (runs/epoxy_amine_smoke10) ─────────────────────────

_SMOKE10_DIR = Path('runs/epoxy_amine_smoke10')


@pytest.fixture(scope='module')
def smoke10_species():
    pytest.importorskip('rdkit')
    from scripts._systems import layout_species

    summary = json.loads(
        (_SMOKE10_DIR / 'summary.json').read_text(encoding='utf-8'))
    manifest = json.loads(
        (_SMOKE10_DIR / 'manifest.json').read_text(encoding='utf-8'))
    return layout_species([
        (summary['epoxy_smiles'], summary['n_epoxies'], manifest['seed']),
        (summary['amine_smiles'], summary['n_amines'], manifest['seed'] + 1),
    ])


@pytest.fixture(scope='module')
def smoke10_snapshots():
    return load_topology_snapshots(_SMOKE10_DIR / 'topology.jsonl')


@pytest.mark.skipif(
    not (_SMOKE10_DIR / 'topology.jsonl').exists(),
    reason='smoke10 run data not available',
)
class TestSmoke10Regression:
    """Measured smoke10 values, asserted as a regression guard.

    Species are rebuilt from summary.json SMILES/counts + manifest.json seed
    via scripts._systems.layout_species (epoxies first then amines, seeds
    (seed, seed + 1)) — the same layout as build_epoxy_amine_system.  The
    rebuild was verified identical to the species recorded in the
    trajectory.jsonl header.
    """

    @pytest.fixture()
    def species(self, smoke10_species):
        return smoke10_species

    @pytest.fixture()
    def snapshots(self, smoke10_snapshots):
        return smoke10_snapshots

    def test_snapshot_structure(self, snapshots):
        assert [(s, c) for s, c, _ in snapshots] == [(0, -1), (4468, 0), (8577, 2)]

    def test_initial_species_counts(self, snapshots, species):
        # 10 DGEBA x 2 epoxide rings, 5 DETA x (2 primary + 1 secondary N);
        # DGEBA glycidyl/aryl ethers are excluded by the one-C-one-H hydroxyl
        # rule, so the initial hydroxyl count is 0.
        row = species_series(snapshots[:1], species)[0]
        assert row['n_epoxide'] == 20
        assert row['n_amine_primary'] == 10
        assert row['n_amine_secondary'] == 5
        assert row['n_amine_tertiary'] == 0
        assert row['n_hydroxyl'] == 0
        # 5 DETA x (2+2+1) N-H hydrogens = 25 (matches summary n_amine_h_sites).
        ranks = amine_ranks(snapshots[0][2], species)
        assert len(ranks) == 15
        assert sum(ranks.values()) == 25

    def test_final_species_counts(self, snapshots, species):
        # 2 ring openings confirmed (summary: confirmed_formations=2,
        # epoxide_conversion=0.1) -> 18 closed epoxides, 2 amines promoted
        # 1° -> 2°.  Hydroxyl is exactly 1, not 2: the cycle-0 opening was
        # zwitterion-like (decisions.md 2026-07-09 E1 実走結果) — its
        # transferring H is bond-less in the recorded snapshot, so only the
        # cycle-2 opening carries the full C-O-H set.  Measured empirically;
        # documented data property, not a bug (see amine_ranks docstring).
        row = species_series(snapshots[-1:], species)[0]
        assert row['n_epoxide'] == 18
        assert row['n_amine_primary'] == 8
        assert row['n_amine_secondary'] == 7
        assert row['n_amine_tertiary'] == 0
        assert row['n_hydroxyl'] == 1

    def test_gel_indicator_and_crosslinks(self, snapshots, species):
        monomer_sets = monomer_sets_from_bonds(snapshots[0][2])
        assert len(monomer_sets) == 15  # 10 DGEBA + 5 DETA

        # Initial: all monomers separate -> 1/15.  Final: one epoxy-amine
        # dimer -> 2/15 (far below the FS gel point at conversion 0.1).
        assert largest_component_fraction(
            snapshots[0][2], monomer_sets) == pytest.approx(1 / 15)
        assert largest_component_fraction(
            snapshots[-1][2], monomer_sets) == pytest.approx(2 / 15)

        rings0 = find_epoxide_rings(snapshots[0][2], species)
        assert len(rings0) == 20
        counts = crosslink_counts(
            snapshots[-1][2], species, rings0, monomer_sets)
        assert counts == {
            'n_tertiary_amine': 0,
            'n_fully_reacted_epoxy_monomers': 0,
            'n_monomers': 15,
        }
