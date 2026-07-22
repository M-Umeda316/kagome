"""Tests for depth-resolved reaction density (Eq. 13)."""
import numpy as np
import pytest

from kagome.analysis.density import reaction_density_profile
from kagome.io.trajectory import TrajectoryFrame
from kagome.reactive.bonds import BondEvent


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

    def test_per_event_cell_used(self):
        """A1: cells_at_event supplies each event's own box (NPT varies)."""
        events = [
            BondEvent(step=10, cycle=0, atom_a=0, atom_b=1,
                      event_type='confirmed_formation', distance=1.8),
        ]
        positions = {10: np.array([[5.0, 5.0, 0.5], [5.0, 5.0, 29.5]])}
        # Wide box (30): minimum image wraps the pair -> midpoint near z=0.
        cells_at_event = {10: np.diag([30.0, 30.0, 30.0]).astype(float)}
        z_bins = np.array([0.0, 5.0, 15.0, 25.0, 30.0])

        density = reaction_density_profile(
            events, positions, z_bins, area_xy=10.0,
            cells_at_event=cells_at_event,
        )
        assert density[0] > 0.0, 'per-event cell should wrap pair to z~0'
        assert density[2] == 0.0

    def test_n_frames_halves_density(self):
        """A2: doubling N_frames halves the normalized density."""
        events = [
            BondEvent(step=1, cycle=0, atom_a=0, atom_b=1,
                      event_type='confirmed_formation', distance=1.8),
        ]
        positions = {1: np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]])}
        z_bins = np.array([0.0, 2.0])

        d1 = reaction_density_profile(events, positions, z_bins, 5.0, n_frames=10)
        d2 = reaction_density_profile(events, positions, z_bins, 5.0, n_frames=20)
        assert d2[0] == pytest.approx(0.5 * d1[0])


class TestDensityProfilePlotting:
    """A3: the script must not silently invent an area; it needs a cell or
    an explicit --cell-xy-area."""

    def _write_jsonl(self, path, records):
        import json
        with open(path, 'w', encoding='utf-8') as f:
            for r in records:
                f.write(json.dumps(r) + '\n')

    def test_missing_cell_and_area_raises(self, tmp_path):
        pytest.importorskip('matplotlib')
        from scripts.reproduce_figures import plot_density_profile

        traj = tmp_path / 'trajectory.jsonl'
        bonds = tmp_path / 'bonds.jsonl'
        # Frame with no 'cell' field (old schema) and no area supplied.
        self._write_jsonl(traj, [
            {'_header': True, 'schema_version': 1, 'species': ['C', 'C'],
             'n_atoms': 2, 'save_interval': 1},
            {'step': 100, 'time_fs': 0.0, 'phase': 'unbiased', 'cycle': 0,
             'energy_base': 0.0, 'energy_bias': 0.0, 'energy_total': 0.0,
             'positions': [[0.0, 0.0, 1.0], [0.0, 0.0, 2.0]]},
        ])
        self._write_jsonl(bonds, [
            {'step': 100, 'cycle': 0, 'atom_a': 0, 'atom_b': 1,
             'event_type': 'confirmed_formation', 'distance': 1.5, 'r0': 1.5},
        ])

        with pytest.raises(ValueError, match='cross-sectional area'):
            plot_density_profile(bonds, traj, tmp_path / 'figs')

    def test_density_excludes_non_counting_formations(self, tmp_path, monkeypatch):
        """A5: rho_rxn excludes water-forming (counts_as_reaction=False) events.

        A nylon-type run mixes the primary amide bond (counts) with the water
        O-H formation (does not count). A vinyl-only run has no such events, so
        adding one must not change which formations reach the density kernel.
        """
        pytest.importorskip('matplotlib')
        import scripts.reproduce_figures as rf

        captured = {}

        def _fake_density(formations, positions, z_bins, area, n_frames,
                          cells_at_event=None):
            captured['pairs'] = {(e.atom_a, e.atom_b) for e in formations}
            return np.zeros(len(z_bins) - 1)

        monkeypatch.setattr(rf, 'reaction_density_profile', _fake_density)

        traj = tmp_path / 'trajectory.jsonl'
        bonds = tmp_path / 'bonds.jsonl'
        self._write_jsonl(traj, [
            {'_header': True, 'schema_version': 1, 'species': ['C'] * 4,
             'n_atoms': 4, 'save_interval': 1},
            {'step': 100, 'time_fs': 0.0, 'phase': 'unbiased', 'cycle': 0,
             'energy_base': 0.0, 'energy_bias': 0.0, 'energy_total': 0.0,
             'positions': [[0.0, 0.0, 1.0]] * 4,
             'cell': [[10.0, 0, 0], [0, 10.0, 0], [0, 0, 10.0]]},
        ])
        self._write_jsonl(bonds, [
            # Primary amide bond (counts) — the vinyl-equivalent real reaction.
            {'step': 100, 'cycle': 0, 'atom_a': 0, 'atom_b': 1,
             'event_type': 'confirmed_formation', 'distance': 1.5, 'r0': 1.5,
             'counts_as_reaction': True},
            # Water O-H formation (nylon-only) — must be excluded.
            {'step': 100, 'cycle': 0, 'atom_a': 2, 'atom_b': 3,
             'event_type': 'confirmed_formation', 'distance': 1.5, 'r0': 1.5,
             'counts_as_reaction': False},
        ])

        rf.plot_density_profile(bonds, traj, tmp_path / 'figs')
        # Only the counting formation reaches the density kernel; the water
        # event (2,3) is filtered out — identical to a vinyl-only result.
        assert captured['pairs'] == {(0, 1)}

    def test_density_excludes_mix_settle_frames_from_n_frames(self, tmp_path, monkeypatch):
        """WM-P2 (decisions.md 追補 2026-07-18 (k) extension): mix_settle
        frames are real MLIP frames written by _run_plain_md during a
        classical-mixing run's settle stage. No reaction can fire during
        them, so they must not inflate N_frames in
        rho_rxn(z) = N_rxn(z)/(A*dz*N_frames) -- mirroring the energy-plot's
        unbiased_mask exclusion (reproduce_figures.py ~line 68-69).

        Asserts plot_density_profile passes the filtered frame count
        (excluding mix_settle) to reaction_density_profile, and that using
        the correct (smaller) N_frames yields a larger density than counting
        every trajectory frame would (the pre-fix behavior).
        """
        pytest.importorskip('matplotlib')
        import scripts.reproduce_figures as rf

        captured = {}
        real_density = rf.reaction_density_profile

        def _spy_density(formations, positions, z_bins, area, n_frames,
                         cells_at_event=None):
            captured['n_frames'] = n_frames
            return real_density(formations, positions, z_bins, area,
                                 n_frames=n_frames, cells_at_event=cells_at_event)

        monkeypatch.setattr(rf, 'reaction_density_profile', _spy_density)

        traj = tmp_path / 'trajectory.jsonl'
        bonds = tmp_path / 'bonds.jsonl'
        cell = [[10.0, 0, 0], [0, 10.0, 0], [0, 0, 10.0]]

        def _rec(step, phase):
            return {'step': step, 'time_fs': 0.0, 'phase': phase, 'cycle': 0,
                     'energy_base': 0.0, 'energy_bias': 0.0, 'energy_total': 0.0,
                     'positions': [[0.0, 0.0, 3.0], [0.0, 0.0, 3.0]], 'cell': cell}

        self._write_jsonl(traj, [
            {'_header': True, 'schema_version': 1, 'species': ['C', 'C'],
             'n_atoms': 2, 'save_interval': 1},
            # 3 classical mix_settle frames (would inflate N_frames pre-fix).
            _rec(1, 'mix_settle'),
            _rec(2, 'mix_settle'),
            _rec(3, 'mix_settle'),
            # 2 MLIP-analysis-window frames; the event lands on step 100.
            _rec(100, 'unbiased'),
            _rec(101, 'unbiased'),
        ])
        self._write_jsonl(bonds, [
            {'step': 100, 'cycle': 0, 'atom_a': 0, 'atom_b': 1,
             'event_type': 'confirmed_formation', 'distance': 1.5, 'r0': 1.5},
        ])

        rf.plot_density_profile(bonds, traj, tmp_path / 'figs', n_z_bins=5)

        # 5 total frames, 3 of them mix_settle -> only the 2 unbiased frames
        # should reach N_frames. (Pre-fix this would have been 5.)
        assert captured['n_frames'] == 2

        # Directly confirm the practical consequence: normalizing by the
        # correct (filtered) N_frames gives a larger density than the buggy
        # all-frames count would have.
        events = [
            BondEvent(step=100, cycle=0, atom_a=0, atom_b=1,
                      event_type='confirmed_formation', distance=1.5),
        ]
        positions = {100: np.array([[0.0, 0.0, 3.0], [0.0, 0.0, 3.0]])}
        z_bins = np.linspace(0.0, 10.0, 6)
        density_correct = reaction_density_profile(
            events, positions, z_bins, area_xy=100.0, n_frames=2)
        density_buggy_all_frames = reaction_density_profile(
            events, positions, z_bins, area_xy=100.0, n_frames=5)
        assert np.sum(density_correct) > np.sum(density_buggy_all_frames)


class TestAnalysisFramesHelper:
    """WM-P2 (decisions.md 追補 2026-07-18 (k) extension): the pure helper
    behind the density-profile N_frames fix, tested in isolation from
    matplotlib/figure-saving so the phase-filtering logic itself is covered.
    """

    @staticmethod
    def _frame(step: int, phase: str) -> TrajectoryFrame:
        return TrajectoryFrame(
            step=step, time_fs=0.0, phase=phase, cycle=0,
            energy_base=0.0, energy_bias=0.0, energy_total=0.0,
            positions=[[0.0, 0.0, 0.0]],
        )

    def test_excludes_mix_settle_and_mixing(self):
        from scripts.reproduce_figures import _analysis_frames

        frames = [
            self._frame(1, 'mixing'),
            self._frame(2, 'mix_settle'),
            self._frame(3, 'biased'),
            self._frame(4, 'unbiased'),
            self._frame(5, 'mix_settle'),
        ]
        result = _analysis_frames(frames)
        assert [f.step for f in result] == [3, 4]

    def test_no_op_when_no_mixing_phases(self):
        """Runs without mixing carry neither phase -> filter is a no-op."""
        from scripts.reproduce_figures import _analysis_frames

        frames = [self._frame(1, 'biased'), self._frame(2, 'unbiased')]
        result = _analysis_frames(frames)
        assert result == frames
