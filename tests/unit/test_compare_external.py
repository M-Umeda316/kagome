"""Tests for the E2 external comparison helpers (Provenzano 2025, ref#2)."""
from pathlib import Path

import pytest

from scripts.compare_epoxy_external import (
    external_metrics,
    network_metrics,
    parse_xl_trend,
    run_series,
)

_REPO_ROOT = Path(__file__).resolve().parents[2]
_XL_TREND_FIXTURE = _REPO_ROOT / 'tests' / 'data' / 'provenzano2025' / 'xl_trend.txt'
_EXTERNAL_DIR = _REPO_ROOT / 'data' / 'external' / 'provenzano2025' / 'xlinker'


class TestParseXlTrend:

    def test_committed_fixture(self):
        # Redistributed CC-BY-4.0 fixture (tests/data/provenzano2025/README.md):
        # 50 rows of (radius A, iteration, crosslink %), 45% example run.
        rows = parse_xl_trend(_XL_TREND_FIXTURE)
        assert len(rows) == 50
        assert rows[0] == (2.0, 1, 0.7)
        assert rows[-1] == (3.0, 10, 45.0)
        # Schedule: radii 2.0 -> 3.0 A in 0.25 steps, 10 iterations each.
        assert sorted({r for r, _i, _p in rows}) == [2.0, 2.25, 2.5, 2.75, 3.0]
        assert all(1 <= i <= 10 for _r, i, _p in rows)
        # Crosslink degree is monotonically non-decreasing (no smoothing).
        percents = [p for _r, _i, p in rows]
        assert all(a <= b for a, b in zip(percents, percents[1:]))

    def test_malformed_row_raises(self, tmp_path):
        path = tmp_path / 'xl_trend.txt'
        path.write_text('2.0 1 0.7 extra\n', encoding='utf-8')
        with pytest.raises(ValueError, match='3 columns'):
            parse_xl_trend(path)

    def test_header_and_comments_skipped(self, tmp_path):
        path = tmp_path / 'xl_trend.txt'
        path.write_text(
            '#Trend of crosslink \n'
            'Radi    Iter    %\n'
            '----    ----    ----\n'
            '2.0 1 0.7\n',
            encoding='utf-8')
        assert parse_xl_trend(path) == [(2.0, 1, 0.7)]


# Minimal ring-opening addition (same fixture style as test_network.py):
# epoxide (C0, C1, O2) + primary amine (N3 with H4, H5, backbone C6).
# The event opens the ring (remove 0-2), forms N-C (0-3), and transfers H4
# from N to O (remove 3-4, add 2-4).
_SPECIES = ['C', 'C', 'O', 'N', 'H', 'H', 'C']
_BONDS_BEFORE = [(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (3, 6)]
_BONDS_AFTER = [(0, 1), (1, 2), (0, 3), (3, 5), (3, 6), (2, 4)]


class TestNetworkMetrics:

    def test_unreacted_state_is_all_zero(self):
        metrics = network_metrics(_BONDS_BEFORE, _BONDS_BEFORE, _SPECIES)
        assert metrics['epoxide_conversion'] == pytest.approx(0.0)
        assert metrics['amine_h_conversion'] == pytest.approx(0.0)
        assert metrics['n_epoxide_initial'] == 1
        assert metrics['amine_h_initial'] == 2
        assert metrics['n_monomers'] == 2
        assert metrics['n_epoxy_monomers'] == 1
        assert metrics['n_amine_n'] == 1
        # Two disconnected monomers -> largest component 1/2.
        assert metrics['largest_component_fraction'] == pytest.approx(0.5)
        assert metrics['tertiary_amine_fraction'] == pytest.approx(0.0)
        assert metrics['fully_reacted_epoxy_fraction'] == pytest.approx(0.0)

    def test_ring_opening_addition(self):
        metrics = network_metrics(_BONDS_BEFORE, _BONDS_AFTER, _SPECIES)
        # The single ring opened: epoxide conversion 1.0 (ring basis).
        assert metrics['epoxide_conversion'] == pytest.approx(1.0)
        # One of two N-H consumed: amine-H conversion 0.5 (ref#2 basis).
        assert metrics['amine_h_conversion'] == pytest.approx(0.5)
        # N-C bond bridges the two monomers -> largest component 2/2.
        assert metrics['largest_component_fraction'] == pytest.approx(1.0)
        # N3 still has H5 -> secondary, not tertiary.
        assert metrics['n_tertiary_amine'] == 0
        assert metrics['tertiary_amine_fraction'] == pytest.approx(0.0)
        # The only epoxy monomer has its only ring open -> fraction 1.0.
        assert metrics['n_fully_reacted_epoxy_monomers'] == 1
        assert metrics['fully_reacted_epoxy_fraction'] == pytest.approx(1.0)

    def test_run_series_rows(self):
        snapshots = [(0, -1, _BONDS_BEFORE), (100, 0, _BONDS_AFTER)]
        series = run_series(snapshots, _SPECIES)
        assert [row['cycle'] for row in series] == [-1, 0]
        assert series[0]['epoxide_conversion'] == pytest.approx(0.0)
        assert series[1]['epoxide_conversion'] == pytest.approx(1.0)
        assert series[1]['amine_h_conversion'] == pytest.approx(0.5)
        assert series[1]['largest_component_fraction'] == pytest.approx(1.0)

    def test_run_series_empty(self):
        assert run_series([], _SPECIES) == []


# ── Real-data regression (data/external/provenzano2025, not committed) ──────

@pytest.mark.skipif(
    not (_EXTERNAL_DIR / 'data.xlinker').exists(),
    reason='Provenzano 2025 dataset not fetched '
           '(run scripts/fetch_provenzano2025.py)',
)
class TestProvenzanoRegression:
    """Metrics computed from the published structures (45% example).

    data.relaxed00 = initial 800 DGEBA + 320 DETA structure (~45,600 atoms);
    data.xlinker = the published 45%-conversion crosslinked example.  The
    amine-H conversion computed by us must land near the 45% reported in
    xl_trend.txt (their basis: reacted amine H / initial reactive sites).
    """

    @pytest.fixture(scope='class')
    def metrics(self):
        return external_metrics(
            _EXTERNAL_DIR / 'data.relaxed00',
            _EXTERNAL_DIR / 'data.xlinker',
        )

    def test_system_size(self, metrics):
        assert metrics['n_atoms'] == 45600, (
            f'n_atoms = {metrics["n_atoms"]}')
        assert metrics['n_monomers'] == 1120, (
            f'n_monomers = {metrics["n_monomers"]} (expected 800 DGEBA '
            f'+ 320 DETA)')
        assert metrics['n_epoxy_monomers'] == 800, (
            f'n_epoxy_monomers = {metrics["n_epoxy_monomers"]}')
        assert metrics['n_amine_n'] == 960, (
            f'n_amine_n = {metrics["n_amine_n"]} (expected 320 DETA x 3 N)')

    def test_initial_reactive_sites(self, metrics):
        # 800 DGEBA x 2 epoxide rings; 320 DETA x 5 N-H (stoichiometric 1:1).
        assert metrics['n_epoxide_initial'] == 1600, (
            f'n_epoxide_initial = {metrics["n_epoxide_initial"]}')
        assert metrics['amine_h_initial'] == 1600, (
            f'amine_h_initial = {metrics["amine_h_initial"]}')

    def test_amine_h_conversion_near_45_percent(self, metrics):
        assert metrics['amine_h_conversion'] == pytest.approx(0.45, abs=0.05), (
            f'amine_h_conversion = {metrics["amine_h_conversion"]:.4f}, '
            f'epoxide_conversion = {metrics["epoxide_conversion"]:.4f}')

    def test_largest_component_fraction_in_unit_interval(self, metrics):
        lcf = metrics['largest_component_fraction']
        assert 0.0 < lcf <= 1.0, f'largest_component_fraction = {lcf}'

    def test_crosslink_fractions_recorded(self, metrics):
        # Learn/record the actual numbers via assert messages.
        tert = metrics['tertiary_amine_fraction']
        full = metrics['fully_reacted_epoxy_fraction']
        assert 0.0 <= tert <= 1.0, f'tertiary_amine_fraction = {tert}'
        assert 0.0 <= full <= 1.0, f'fully_reacted_epoxy_fraction = {full}'
