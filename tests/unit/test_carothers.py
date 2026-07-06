"""Tests for Carothers equation (step-growth polymerization)."""
import logging

import numpy as np
import pytest

from kagome.analysis.carothers import dpn_carothers, dpn_from_bonds


class TestDpnCarothers:

    def test_zero_conversion(self):
        assert dpn_carothers(0.0) == pytest.approx(1.0)

    def test_half_conversion(self):
        assert dpn_carothers(0.5) == pytest.approx(2.0)

    def test_ninety_percent(self):
        assert dpn_carothers(0.9) == pytest.approx(10.0)

    def test_ninety_nine_percent(self):
        assert dpn_carothers(0.99) == pytest.approx(100.0)

    def test_array_input(self):
        p = np.array([0.0, 0.5, 0.9])
        result = dpn_carothers(p)
        expected = np.array([1.0, 2.0, 10.0])
        np.testing.assert_allclose(result, expected)

    def test_full_conversion_clamped_not_inf(self, caplog):
        # F9: p = 1.0 must yield a finite DPn (10000), not inf, and warn.
        with caplog.at_level(logging.WARNING):
            dpn = dpn_carothers(1.0)
        assert np.isfinite(dpn)
        assert dpn == pytest.approx(10000.0)  # 1/(1-0.9999)
        assert 'clamp' in caplog.text

    def test_overshoot_clamped_not_negative(self, caplog):
        # p > 1.0 (miscount) must clamp to a finite positive DPn, not go negative.
        with caplog.at_level(logging.WARNING):
            dpn = dpn_carothers(1.5)
        assert np.isfinite(dpn)
        assert dpn == pytest.approx(10000.0)
        assert 'clamp' in caplog.text

    def test_array_with_full_conversion(self, caplog):
        with caplog.at_level(logging.WARNING):
            result = dpn_carothers(np.array([0.5, 1.0]))
        assert np.all(np.isfinite(result))
        np.testing.assert_allclose(result, np.array([2.0, 10000.0]))
        assert 'clamp' in caplog.text

    def test_no_warning_below_clamp(self, caplog):
        with caplog.at_level(logging.WARNING):
            dpn_carothers(0.99)
        assert caplog.text == ''


class TestDpnFromBonds:

    def test_no_bonds(self):
        assert dpn_from_bonds(0, 20) == pytest.approx(1.0)

    def test_half_conversion(self):
        result = dpn_from_bonds(5, 20)
        assert result == pytest.approx(2.0)

    def test_zero_groups(self):
        assert dpn_from_bonds(0, 0) == pytest.approx(1.0)

    def test_near_full_conversion(self):
        result = dpn_from_bonds(9, 20)
        assert result == pytest.approx(10.0)

    def test_clamp_warns_on_overshoot(self, caplog):
        # RF19b: p>1 (miscount) or full conversion must not be silently clamped.
        with caplog.at_level(logging.WARNING):
            dpn = dpn_from_bonds(n_bonds=100, n_functional_groups=100)  # p = 2.0
        assert dpn == pytest.approx(10000.0)  # 1/(1-0.9999)
        assert 'clamp' in caplog.text

    def test_no_warning_below_clamp(self, caplog):
        with caplog.at_level(logging.WARNING):
            dpn_from_bonds(n_bonds=9, n_functional_groups=20)  # p = 0.9
        assert caplog.text == ''
