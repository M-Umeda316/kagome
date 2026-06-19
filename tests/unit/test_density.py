"""Tests for depth-resolved reaction density (Eq. 13)."""
import numpy as np
import pytest

from src.analysis.density import reaction_density_profile
from src.reactive.bonds import BondEvent


class TestReactionDensityProfile:

    def test_single_event(self):
        events = [
            BondEvent(step=100, cycle=0, atom_a=0, atom_b=1,
                      event_type='confirmed_formation', distance=1.8),
        ]
        positions = {
            100: np.array([[0.0, 0.0, 2.5], [0.0, 0.0, 3.5]]),  # midpoint z=3.0
        }
        z_bins = np.array([0.0, 2.0, 4.0, 6.0])
        area_xy = 10.0

        density = reaction_density_profile(events, positions, z_bins, area_xy)

        assert density.shape == (3,)
        # event at z=3.0 falls in bin [2.0, 4.0) → index 1
        assert density[0] == 0.0
        assert density[1] > 0.0
        assert density[2] == 0.0

    def test_normalization(self):
        events = [
            BondEvent(step=100, cycle=0, atom_a=0, atom_b=1,
                      event_type='confirmed_formation', distance=1.8),
        ]
        positions = {100: np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]])}
        z_bins = np.array([0.0, 2.0])
        area_xy = 5.0
        n_frames = 2

        density = reaction_density_profile(events, positions, z_bins, area_xy, n_frames)

        # 1 event / (5.0 * 2.0 * 2) = 0.05
        assert density[0] == pytest.approx(0.05)

    def test_no_events(self):
        z_bins = np.array([0.0, 5.0, 10.0])
        density = reaction_density_profile([], {}, z_bins, area_xy=10.0)
        np.testing.assert_array_equal(density, 0.0)

    def test_missing_positions_skipped(self):
        events = [
            BondEvent(step=999, cycle=0, atom_a=0, atom_b=1,
                      event_type='confirmed_formation', distance=1.8),
        ]
        z_bins = np.array([0.0, 5.0])
        density = reaction_density_profile(events, {}, z_bins, area_xy=10.0)
        assert density[0] == 0.0

    def test_pbc_midpoint_correction(self):
        """Atoms straddling the PBC boundary should have midpoint near z=0."""
        cell = np.diag([30.0, 30.0, 30.0])
        events = [
            BondEvent(step=100, cycle=0, atom_a=0, atom_b=1,
                      event_type='confirmed_formation', distance=1.8),
        ]
        positions = {
            100: np.array([[5.0, 5.0, 0.5], [5.0, 5.0, 29.5]]),
        }
        z_bins = np.array([0.0, 5.0, 15.0, 25.0, 30.0])
        area_xy = 10.0

        density_pbc = reaction_density_profile(
            events, positions, z_bins, area_xy, cell=cell,
        )
        density_no_pbc = reaction_density_profile(
            events, positions, z_bins, area_xy, cell=None,
        )
        assert density_pbc[0] > 0.0, 'PBC midpoint should be near z=0'
        assert density_no_pbc[2] > 0.0, 'naive midpoint at z=15 without PBC'
        assert density_pbc[2] == 0.0, 'PBC should not place event at z=15'
