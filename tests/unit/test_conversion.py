"""Tests for conversion tracking (Eq. 11-12)."""
import numpy as np
import pytest

from src.analysis.conversion import conversion, conversion_timeseries
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
