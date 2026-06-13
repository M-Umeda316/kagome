"""Tests for MC barostat (NPT ensemble, T-B)."""
import numpy as np
import pytest

from src.integrators.mc_barostat import MCBarostat, MCBarostatParams
from src.units import ATM_TO_KCAL_MOL_A3, KB


class _FlatCalculator:
    """Fake calculator: returns zero energy and zero forces regardless of configuration."""
    name = 'flat'

    def compute(self, positions, species, cell):
        n = positions.shape[0]
        return 0.0, np.zeros((n, 3))


class _WellCalculator:
    """Fake calculator: energy proportional to volume deviation from a reference."""
    name = 'well'

    def __init__(self, V_ref: float, k: float = 0.1):
        self.V_ref = V_ref
        self.k = k

    def compute(self, positions, species, cell):
        box = np.array([cell[0, 0], cell[1, 1], cell[2, 2]])
        V = float(np.prod(box))
        energy = 0.5 * self.k * (V - self.V_ref) ** 2
        n = positions.shape[0]
        return energy, np.zeros((n, 3))


class TestMCBarostatParams:

    def test_default_pressure_is_1atm(self):
        p = MCBarostatParams()
        assert p.pressure_atm == pytest.approx(1.0)

    def test_pressure_conversion(self):
        # 1 atm -> kcal/(mol*A^3): known from src/units.py
        assert ATM_TO_KCAL_MOL_A3 == pytest.approx(1.4596e-5, rel=1e-3)


class TestMCBarostatShouldAttempt:

    def setup_method(self):
        self.b = MCBarostat(MCBarostatParams(frequency=25))

    def test_step_0_is_false(self):
        assert not self.b.should_attempt(0)

    def test_multiple_of_frequency_is_true(self):
        assert self.b.should_attempt(25)
        assert self.b.should_attempt(50)
        assert self.b.should_attempt(100)

    def test_non_multiple_is_false(self):
        assert not self.b.should_attempt(24)
        assert not self.b.should_attempt(26)


class TestMCBarostatStep:

    def _cell(self, box: float = 10.0) -> np.ndarray:
        return np.diag([box, box, box])

    def _positions(self, n: int = 4, box: float = 10.0) -> np.ndarray:
        rng = np.random.default_rng(0)
        return rng.uniform(0, box, (n, 3))

    def test_always_accepts_zero_energy_change(self):
        """With a flat potential, delta_E=0, so acceptance depends on P*dV and Jacobian."""
        barostat = MCBarostat(MCBarostatParams(pressure_atm=1.0, max_volume_change_frac=1e-6))
        rng = np.random.default_rng(42)
        cell = self._cell(10.0)
        positions = self._positions(4, 10.0)
        calc = _FlatCalculator()

        # With tiny volume change and flat potential, almost always accepted
        n_accepted = sum(
            barostat.try_step(positions.copy(), ['C'] * 4, cell.copy(), 0.0, calc, rng, 300.0)[0]
            for _ in range(50)
        )
        assert n_accepted > 0, 'At least some steps should be accepted'

    def test_accepted_step_changes_cell(self):
        barostat = MCBarostat(MCBarostatParams(pressure_atm=0.0, max_volume_change_frac=0.05))
        rng = np.random.default_rng(7)
        calc = _FlatCalculator()

        cell = self._cell(10.0)
        positions = self._positions(4, 10.0)
        original_cell = cell.copy()
        original_pos = positions.copy()

        accepted, _, _ = barostat.try_step(positions, ['C'] * 4, cell, 0.0, calc, rng, 300.0)
        if accepted:
            assert not np.allclose(cell, original_cell), 'Cell should change on acceptance'

    def test_rejected_step_preserves_cell(self):
        """Force rejection by making volume expansion very costly."""
        # Use extremely high pressure so volume expansion is always rejected
        barostat = MCBarostat(MCBarostatParams(pressure_atm=1e12, max_volume_change_frac=0.1))
        rng = np.random.default_rng(99)
        calc = _FlatCalculator()

        cell = self._cell(5.0)
        positions = self._positions(4, 5.0)
        original_cell = cell.copy()

        # Most steps will be rejected at ultra-high pressure (volume expansion is costly)
        rejections = 0
        for _ in range(20):
            cell_copy = cell.copy()
            pos_copy = positions.copy()
            accepted, _, _ = barostat.try_step(pos_copy, ['C'] * 4, cell_copy, 0.0, calc, rng, 300.0)
            if not accepted:
                rejections += 1
                np.testing.assert_allclose(cell_copy, original_cell)

        assert rejections > 0

    def test_stats_track_attempts_and_acceptances(self):
        barostat = MCBarostat(MCBarostatParams(pressure_atm=1.0, max_volume_change_frac=0.01))
        rng = np.random.default_rng(42)
        calc = _FlatCalculator()
        cell = self._cell(10.0)
        positions = self._positions(4)

        for _ in range(10):
            barostat.try_step(positions.copy(), ['C'] * 4, cell.copy(), 0.0, calc, rng, 300.0)

        assert barostat.stats.attempts == 10
        assert 0 <= barostat.stats.accepted <= 10
        assert 0.0 <= barostat.stats.acceptance_rate <= 1.0

    def test_positions_scale_with_cell(self):
        """When accepted, positions should scale by the same factor as the cell."""
        # P=0, flat potential, tiny move -> high acceptance
        barostat = MCBarostat(MCBarostatParams(pressure_atm=0.0, max_volume_change_frac=0.02))
        rng = np.random.default_rng(5)
        calc = _FlatCalculator()

        cell = np.diag([10.0, 10.0, 10.0])
        positions = np.array([[5.0, 5.0, 5.0], [2.0, 3.0, 4.0]], dtype=float)
        orig_pos = positions.copy()
        orig_box = cell[0, 0]

        accepted, _, _ = barostat.try_step(positions, ['C', 'C'], cell, 0.0, calc, rng, 300.0)

        if accepted:
            scale = cell[0, 0] / orig_box
            np.testing.assert_allclose(positions, orig_pos * scale, rtol=1e-10)
