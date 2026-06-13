"""Tests for Carothers equation (step-growth polymerization)."""
import numpy as np
import pytest

from src.analysis.carothers import dpn_carothers, dpn_from_bonds


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
