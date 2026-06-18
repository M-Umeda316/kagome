"""Tests for conversion tracking (Eq. 11-12)."""
import numpy as np
import pytest

from src.analysis.conversion import conversion, conversion_timeseries, fit_conversion_exponential, monomer_site_count
from src.reactive.bonds import BondEvent
from src.reactive.groups import ReactiveGroup


class TestMonomerSiteCount:
    """RF2: denominator = initial monomer count, not all groups."""

    def test_vinyl_returns_monomer_count(self):
        """vinyl_alpha_C group size = n_monomers."""
        groups = {
            'radical_C': ReactiveGroup('radical_C', [0, 1]),
            'vinyl_alpha_C': ReactiveGroup('vinyl_alpha_C', [2, 3, 4]),
            'chain_C': ReactiveGroup('chain_C', [5]),
            'vinyl_beta_C': ReactiveGroup('vinyl_beta_C', [6, 7, 8]),
        }
        assert monomer_site_count(groups) == 3

    def test_excludes_constraint_groups(self):
        """chain_C and vinyl_beta_C are excluded from denominator."""
        groups = {
            'radical_C': ReactiveGroup('radical_C', [0]),
            'vinyl_alpha_C': ReactiveGroup('vinyl_alpha_C', [1, 2]),
            'chain_C': ReactiveGroup('chain_C', [3, 4, 5, 6]),
            'vinyl_beta_C': ReactiveGroup('vinyl_beta_C', [7, 8, 9, 10]),
        }
        assert monomer_site_count(groups) == 2

    def test_preserves_count_after_propagation(self):
        """Even after atoms are removed from vinyl_alpha_C, the
        caller should capture the initial count before any updates."""
        groups = {
            'vinyl_alpha_C': ReactiveGroup('vinyl_alpha_C', [1, 2, 3]),
        }
        initial = monomer_site_count(groups)
        groups['vinyl_alpha_C'].atom_indices.remove(1)
        assert initial == 3
        assert monomer_site_count(groups) == 2

    def test_missing_group_returns_zero(self):
        groups = {'radical_C': ReactiveGroup('radical_C', [0])}
        assert monomer_site_count(groups) == 0

    def test_custom_monomer_group(self):
        groups = {'amine_N': ReactiveGroup('amine_N', [0, 1, 2])}
        assert monomer_site_count(groups, monomer_group='amine_N') == 3


class TestConversion:

    def test_zero(self):
        assert conversion(0, 10) == 0.0

    def test_full(self):
        assert conversion(10, 10) == 1.0

    def test_partial(self):
        assert conversion(3, 10) == pytest.approx(0.3)

    def test_empty(self):
        assert conversion(0, 0) == 0.0


class TestConversionTimeseries:

    def test_simple(self):
        events = [
            BondEvent(step=100, cycle=0, atom_a=0, atom_b=1,
                      event_type='confirmed_formation', distance=1.8),
            BondEvent(step=300, cycle=1, atom_a=2, atom_b=3,
                      event_type='confirmed_formation', distance=1.9),
        ]
        steps = np.array([0, 100, 200, 300, 400])
        steps_out, alpha = conversion_timeseries(events, n_total_sites=4, step_range=steps)

        assert alpha[0] == 0.0
        assert alpha[1] == pytest.approx(0.25)
        assert alpha[2] == pytest.approx(0.25)
        assert alpha[3] == pytest.approx(0.50)
        assert alpha[4] == pytest.approx(0.50)

    def test_no_events(self):
        steps = np.array([0, 100, 200])
        steps_out, alpha = conversion_timeseries([], n_total_sites=10, step_range=steps)
        np.testing.assert_array_equal(alpha, 0.0)

    def test_filters_non_formation(self):
        events = [
            BondEvent(step=100, cycle=0, atom_a=0, atom_b=1,
                      event_type='attempted_formation', distance=3.0),
        ]
        steps = np.array([0, 100, 200])
        _, alpha = conversion_timeseries(events, n_total_sites=4, step_range=steps)
        np.testing.assert_array_equal(alpha, 0.0)


class TestFitConversionExponential:

    def _synthetic_alpha(self, kp_true: float, steps: np.ndarray, timestep_fs: float = 0.25) -> np.ndarray:
        t = steps.astype(np.float64) * timestep_fs
        return 1.0 - np.exp(-kp_true * t)

    def test_recovers_known_kp(self):
        pytest.importorskip('scipy')
        kp_true = 5e-4  # 1/fs
        steps = np.arange(0, 10000, 100, dtype=np.int64)
        alpha = self._synthetic_alpha(kp_true, steps)
        kp_fit, r2 = fit_conversion_exponential(steps, alpha, timestep_fs=0.25)
        assert kp_fit == pytest.approx(kp_true, rel=0.05), f'kp_fit={kp_fit:.4e} vs kp_true={kp_true:.4e}'
        assert r2 > 0.99

    def test_zero_alpha_returns_zeros(self):
        pytest.importorskip('scipy')
        steps = np.arange(0, 1000, 10, dtype=np.int64)
        alpha = np.zeros(len(steps))
        kp, r2 = fit_conversion_exponential(steps, alpha)
        assert kp == pytest.approx(0.0)
        assert r2 == pytest.approx(0.0)

    def test_too_few_points(self):
        pytest.importorskip('scipy')
        steps = np.array([0, 1], dtype=np.int64)
        alpha = np.array([0.0, 0.5])
        kp, r2 = fit_conversion_exponential(steps, alpha)
        assert kp == pytest.approx(0.0)

    def test_missing_scipy_raises(self, monkeypatch):
        import sys
        # Simulate scipy not installed
        monkeypatch.setitem(sys.modules, 'scipy', None)
        monkeypatch.setitem(sys.modules, 'scipy.optimize', None)
        steps = np.arange(0, 1000, 10, dtype=np.int64)
        alpha = self._synthetic_alpha(1e-4, steps)
        with pytest.raises((ImportError, TypeError)):
            fit_conversion_exponential(steps, alpha)
