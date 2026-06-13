"""Tests for conversion tracking (Eq. 11-12)."""
import numpy as np
import pytest

from src.analysis.conversion import conversion, conversion_timeseries, fit_conversion_exponential
from src.reactive.bonds import BondEvent


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
