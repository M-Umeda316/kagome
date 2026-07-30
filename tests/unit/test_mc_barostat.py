"""Tests for MC barostat (NPT ensemble, T-B)."""
import math

import numpy as np
import pytest

from kagome.integrators.mc_barostat import MCBarostat, MCBarostatParams
from kagome.units import ATM_TO_KCAL_MOL_A3, KB


class _FixedRng:
    """Deterministic stand-in: fixed volume proposal and acceptance draw."""

    def __init__(self, delta_ln_v: float, accept_draw: float):
        self._dlnv = delta_ln_v
        self._r = accept_draw

    def uniform(self, lo, hi):
        return self._dlnv

    def random(self):
        return self._r


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


class _ConstantBiasFn:
    """Deterministic stand-in for bias_energy_fn: returns fixed bias energies
    regardless of the configuration passed in.

    try_step calls bias_energy_fn twice per attempt, in this order (see
    mc_barostat.py ~121-125): first with (new_positions, new_box), then with
    (positions, box). So the first call returns ``new_bias`` and the second
    returns ``old_bias``, giving an exact, deterministic ΔE_bias = new_bias -
    old_bias regardless of the actual scaled positions/cell.
    """

    def __init__(self, new_bias: float, old_bias: float):
        self._values = [new_bias, old_bias]
        self._calls = 0

    def __call__(self, positions, box):
        value = self._values[self._calls]
        self._calls += 1
        return value


class TestMCBarostatParams:

    def test_default_pressure_is_1atm(self):
        p = MCBarostatParams()
        assert p.pressure_atm == pytest.approx(1.0)

    def test_pressure_conversion(self):
        # 1 atm -> kcal/(mol*A^3): see src/units.py (canonical first-principles
        # check lives in tests/unit/test_units.py; RF19 corrected this value).
        assert ATM_TO_KCAL_MOL_A3 == pytest.approx(1.4584e-5, rel=1e-3)


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
        # RF21: must not pass vacuously — require at least one acceptance and
        # assert the cell changed on every accepted move.
        barostat = MCBarostat(MCBarostatParams(pressure_atm=0.0, max_volume_change_frac=0.05))
        rng = np.random.default_rng(7)
        calc = _FlatCalculator()

        original_cell = self._cell(10.0)
        positions = self._positions(4, 10.0)

        n_accepted = 0
        for _ in range(50):
            cell = original_cell.copy()
            accepted, _, _ = barostat.try_step(
                positions.copy(), ['C'] * 4, cell, 0.0, calc, rng, 300.0
            )
            if accepted:
                n_accepted += 1
                assert not np.allclose(cell, original_cell), 'Cell should change on acceptance'

        assert n_accepted > 0, 'expected at least one accepted volume move'

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

    def test_jacobian_uses_n_plus_one(self):
        # RF19b: ln(V)-uniform proposal => acceptance Jacobian is (N+1), not N.
        # Flat potential, P=0: delta_H = -(N+1)*kT*delta_ln_V. For a volume
        # *decrease* (delta_ln_V<0) accept iff draw < exp(-(N+1)|delta_ln_V|).
        n = 4
        dlnv = -0.01
        thresh_n1 = math.exp(-(n + 1) * abs(dlnv))   # ~0.95123
        thresh_n = math.exp(-n * abs(dlnv))          # ~0.96079
        assert thresh_n1 < thresh_n  # the two rules give different thresholds

        cell = self._cell(10.0)
        positions = self._positions(n, 10.0)
        calc = _FlatCalculator()

        def attempt(draw):
            b = MCBarostat(MCBarostatParams(pressure_atm=0.0, max_volume_change_frac=0.01))
            accepted, _, _ = b.try_step(
                positions.copy(), ['C'] * n, cell.copy(), 0.0, calc,
                _FixedRng(dlnv, draw), 300.0,
            )
            return accepted

        # draw between the two thresholds: rejected under (N+1), accepted under N
        assert attempt(0.5 * (thresh_n1 + thresh_n)) is False
        # draw below the (N+1) threshold: accepted
        assert attempt(thresh_n1 - 0.01) is True

    def test_positions_scale_with_cell(self):
        """When accepted, positions should scale by the same factor as the cell.

        RF21: loop until at least one acceptance so the invariant is actually
        exercised (no vacuous pass).
        """
        # P=0, flat potential, tiny move -> high acceptance
        barostat = MCBarostat(MCBarostatParams(pressure_atm=0.0, max_volume_change_frac=0.02))
        rng = np.random.default_rng(5)
        calc = _FlatCalculator()

        orig_pos = np.array([[5.0, 5.0, 5.0], [2.0, 3.0, 4.0]], dtype=float)
        orig_box = 10.0

        n_accepted = 0
        for _ in range(50):
            cell = np.diag([orig_box, orig_box, orig_box]).astype(float)
            positions = orig_pos.copy()
            accepted, _, _ = barostat.try_step(positions, ['C', 'C'], cell, 0.0, calc, rng, 300.0)
            if accepted:
                n_accepted += 1
                scale = cell[0, 0] / orig_box
                np.testing.assert_allclose(positions, orig_pos * scale, rtol=1e-10)

        assert n_accepted > 0, 'expected at least one accepted volume move'


class TestMCBarostatBiasEnergy:
    """bias_energy_fn must fold ΔE_bias into the ΔH acceptance criterion
    (specs/decisions.md 2026-07-03 D2). Covers: (a) a fixed bias shifts the
    accept/reject boundary by exactly ΔE_bias, (b) omitting bias_energy_fn
    (default None) reproduces the pre-bias behaviour exactly, (c) a bias that
    grows with volume makes expansion less likely to be accepted."""

    def _cell(self, box: float = 10.0) -> np.ndarray:
        return np.diag([box, box, box])

    def _positions(self, n: int = 4, box: float = 10.0) -> np.ndarray:
        rng = np.random.default_rng(0)
        return rng.uniform(0, box, (n, 3))

    def test_bias_shifts_delta_h_by_delta_e_bias(self):
        # (a) P=0, flat calculator => dE=0, so without bias
        # delta_H = -(n+1)*kT*delta_ln_V exactly (RF19b Jacobian term).
        n = 4
        dlnv = -0.01
        kT = KB * 300.0
        delta_h_no_bias = -(n + 1) * kT * dlnv

        delta_e_bias = 0.5  # kcal/mol, arbitrary fixed shift added by the bias
        delta_h_with_bias = delta_h_no_bias + delta_e_bias

        thresh_no_bias = math.exp(-delta_h_no_bias / kT)
        thresh_with_bias = math.exp(-delta_h_with_bias / kT)
        assert thresh_with_bias < thresh_no_bias  # adding positive bias tightens acceptance

        # Draw strictly between the two thresholds: accepted without bias,
        # rejected once ΔE_bias is folded in — the boundary moved by exactly
        # delta_e_bias.
        draw = 0.5 * (thresh_with_bias + thresh_no_bias)

        cell = self._cell(10.0)
        positions = self._positions(n, 10.0)
        calc = _FlatCalculator()

        def run(bias_fn):
            b = MCBarostat(MCBarostatParams(pressure_atm=0.0, max_volume_change_frac=0.01))
            accepted, _, _ = b.try_step(
                positions.copy(), ['C'] * n, cell.copy(), 0.0, calc,
                _FixedRng(dlnv, draw), 300.0, bias_energy_fn=bias_fn,
            )
            return accepted

        assert run(None) is True
        assert run(_ConstantBiasFn(new_bias=delta_e_bias, old_bias=0.0)) is False

    def test_bias_energy_fn_none_matches_omitted_argument(self):
        # (b) explicit bias_energy_fn=None must be indistinguishable from
        # leaving the argument out (the pre-bias code path).
        n = 4
        dlnv = -0.01
        cell = self._cell(10.0)
        positions = self._positions(n, 10.0)
        calc = _WellCalculator(V_ref=900.0, k=0.05)  # nonzero dE, exercises dE + bias path

        for draw in (0.1, 0.5, 0.9, 0.99):
            b_omitted = MCBarostat(MCBarostatParams(pressure_atm=1.0, max_volume_change_frac=0.01))
            b_explicit = MCBarostat(MCBarostatParams(pressure_atm=1.0, max_volume_change_frac=0.01))

            accepted_omitted, e_omitted, _ = b_omitted.try_step(
                positions.copy(), ['C'] * n, cell.copy(), 0.0, calc,
                _FixedRng(dlnv, draw), 300.0,
            )
            accepted_explicit, e_explicit, _ = b_explicit.try_step(
                positions.copy(), ['C'] * n, cell.copy(), 0.0, calc,
                _FixedRng(dlnv, draw), 300.0, bias_energy_fn=None,
            )
            assert accepted_omitted == accepted_explicit
            assert e_omitted == e_explicit

    def test_bias_growing_with_volume_penalizes_expansion(self):
        # (c) directional check: a bias energy that increases with volume
        # should make the barostat accept volume-expanding moves less often
        # overall, since ΔE_bias adds a positive penalty on expansion and a
        # negative (favorable) term on contraction.
        n = 4
        cell = self._cell(10.0)
        positions = self._positions(n, 10.0)
        calc = _FlatCalculator()

        def growing_bias(pos, box):
            return 5.0 * float(np.prod(box))

        barostat_no_bias = MCBarostat(MCBarostatParams(pressure_atm=0.0, max_volume_change_frac=0.05))
        barostat_with_bias = MCBarostat(MCBarostatParams(pressure_atm=0.0, max_volume_change_frac=0.05))

        rng_no_bias = np.random.default_rng(3)
        rng_with_bias = np.random.default_rng(3)  # same seed => same proposal sequence

        for _ in range(200):
            barostat_no_bias.try_step(
                positions.copy(), ['C'] * n, cell.copy(), 0.0, calc, rng_no_bias, 300.0,
            )
            barostat_with_bias.try_step(
                positions.copy(), ['C'] * n, cell.copy(), 0.0, calc, rng_with_bias, 300.0,
                bias_energy_fn=growing_bias,
            )

        assert barostat_with_bias.stats.acceptance_rate < barostat_no_bias.stats.acceptance_rate
