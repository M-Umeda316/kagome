"""Tests for TDBB potentials and forces.

Validates against paper equations 2-5, 8.
"""
import numpy as np
import pytest

from kagome.boost.tdbb import (
    BoostState,
    PairBias,
    TDBBParams,
    boost_amplitude,
    dissociation_force_magnitude,
    dissociation_potential,
    formation_force_magnitude,
    formation_potential,
    target_distance,
    total_bias,
    total_bias_fast,
)
from kagome.geometry import validated_box


class TestBoostAmplitude:
    """Eq. 5: f1(t) = min(γt, f1_max)."""

    def test_linear_ramp(self):
        assert boost_amplitude(100, gamma=1.0, f1_max=250.0) == 100.0

    def test_saturation(self):
        assert boost_amplitude(300, gamma=1.0, f1_max=250.0) == 250.0

    def test_exact_boundary(self):
        assert boost_amplitude(250, gamma=1.0, f1_max=250.0) == 250.0

    def test_zero(self):
        assert boost_amplitude(0, gamma=1.0, f1_max=250.0) == 0.0

    def test_custom_gamma(self):
        assert boost_amplitude(50, gamma=2.0, f1_max=250.0) == 100.0


class TestTargetDistance:
    """Eq. 4: r0 = λ Σ r_a^vdw."""

    def test_two_atoms(self):
        vdw = np.array([1.70, 1.52])  # C, O
        r0 = target_distance(vdw, lambda_vdw=0.60)
        assert r0 == pytest.approx(0.60 * (1.70 + 1.52))

    def test_single_atom(self):
        r0 = target_distance(np.array([1.5]), lambda_vdw=0.60)
        assert r0 == pytest.approx(0.90)


class TestFormationPotential:
    """Eq. 2: V^f = f1·(1 - exp(-f2·(r-r0)²))."""

    def test_at_r0_is_zero(self):
        r = np.array([2.0])
        v = formation_potential(r, r0=2.0, f1=250.0, f2=10.0)
        assert v[0] == pytest.approx(0.0)

    def test_far_from_r0_approaches_f1(self):
        r = np.array([100.0])
        v = formation_potential(r, r0=2.0, f1=250.0, f2=10.0)
        assert v[0] == pytest.approx(250.0, rel=1e-6)

    def test_symmetric_around_r0(self):
        v_above = formation_potential(np.array([2.5]), r0=2.0, f1=1.0, f2=10.0)
        v_below = formation_potential(np.array([1.5]), r0=2.0, f1=1.0, f2=10.0)
        assert v_above[0] == pytest.approx(v_below[0])

    def test_vectorized(self):
        r = np.array([1.0, 2.0, 3.0])
        v = formation_potential(r, r0=2.0, f1=100.0, f2=10.0)
        assert v.shape == (3,)
        assert v[1] == pytest.approx(0.0)


def _central_difference(func, r, h=1e-6, **kwargs):
    """Central finite-difference dV/dr of a scalar potential at distance r."""
    up = func(np.array([r + h]), **kwargs)[0]
    dn = func(np.array([r - h]), **kwargs)[0]
    return (up - dn) / (2.0 * h)


class TestFormationForce:
    """dV^f/dr = f1·2f2·(r-r0)·exp(-f2·(r-r0)²)."""

    def test_zero_at_r0(self):
        dv = formation_force_magnitude(np.array([2.0]), r0=2.0, f1=250.0, f2=10.0)
        assert dv[0] == pytest.approx(0.0)

    def test_positive_above_r0(self):
        dv = formation_force_magnitude(np.array([2.1]), r0=2.0, f1=250.0, f2=10.0)
        assert dv[0] > 0.0

    def test_negative_below_r0(self):
        dv = formation_force_magnitude(np.array([1.9]), r0=2.0, f1=250.0, f2=10.0)
        assert dv[0] < 0.0

    @pytest.mark.parametrize('r', [1.5, 1.9, 2.0, 2.1, 2.5, 3.0])
    def test_magnitude_matches_analytic(self, r):
        """RF14: lock the numeric magnitude, not just the sign.

        A regression that drops the analytic factor of 2 or mis-scales f2 keeps
        the sign but changes this value, so it would be caught here.
        """
        r0, f1, f2 = 2.0, 250.0, 10.0
        dr = r - r0
        expected = f1 * 2.0 * f2 * dr * np.exp(-f2 * dr ** 2)
        got = formation_force_magnitude(np.array([r]), r0=r0, f1=f1, f2=f2)[0]
        np.testing.assert_allclose(got, expected, rtol=1e-12)

    @pytest.mark.parametrize('r', [1.5, 1.9, 2.1, 2.5, 3.0])
    def test_matches_finite_difference_of_potential(self, r):
        """dV^f/dr must equal the central difference of formation_potential."""
        r0, f1, f2 = 2.0, 250.0, 10.0
        fd = _central_difference(formation_potential, r, r0=r0, f1=f1, f2=f2)
        got = formation_force_magnitude(np.array([r]), r0=r0, f1=f1, f2=f2)[0]
        np.testing.assert_allclose(got, fd, rtol=1e-5)


class TestDissociationPotential:
    """Eq. 3: V^d = f1·exp(-f2·r²)."""

    def test_at_origin(self):
        v = dissociation_potential(np.array([0.0]), f1=125.0, f2=10.0)
        assert v[0] == pytest.approx(125.0)

    def test_decays_with_distance(self):
        v = dissociation_potential(np.array([0.0, 1.0, 2.0]), f1=125.0, f2=10.0)
        assert v[0] > v[1] > v[2]

    def test_far_approaches_zero(self):
        v = dissociation_potential(np.array([10.0]), f1=125.0, f2=10.0)
        assert v[0] == pytest.approx(0.0, abs=1e-10)


class TestDissociationForce:
    """dV^d/dr = -2f1·f2·r·exp(-f2·r²)."""

    def test_zero_at_origin(self):
        dv = dissociation_force_magnitude(np.array([0.0]), f1=125.0, f2=10.0)
        assert dv[0] == pytest.approx(0.0)

    def test_always_negative(self):
        r = np.array([0.5, 1.0, 2.0])
        dv = dissociation_force_magnitude(r, f1=125.0, f2=10.0)
        assert np.all(dv < 0)

    @pytest.mark.parametrize('r', [0.2, 0.5, 1.0, 1.5, 2.0])
    def test_magnitude_matches_analytic(self, r):
        f1, f2 = 125.0, 10.0
        expected = -2.0 * f1 * f2 * r * np.exp(-f2 * r ** 2)
        got = dissociation_force_magnitude(np.array([r]), f1=f1, f2=f2)[0]
        np.testing.assert_allclose(got, expected, rtol=1e-12)

    @pytest.mark.parametrize('r', [0.2, 0.5, 1.0, 1.5])
    def test_matches_finite_difference_of_potential(self, r):
        f1, f2 = 125.0, 10.0
        fd = _central_difference(dissociation_potential, r, f1=f1, f2=f2)
        got = dissociation_force_magnitude(np.array([r]), f1=f1, f2=f2)[0]
        np.testing.assert_allclose(got, fd, rtol=1e-5)


class TestBoostState:

    def test_advance_ramps(self):
        state = BoostState()
        state.advance(gamma=1.0, f1_max_form=250.0, f1_max_dissoc=125.0)
        assert state.step == 1
        assert state.f1_formation == 1.0
        assert state.f1_dissociation == 1.0

    def test_advance_saturates(self):
        state = BoostState()
        for _ in range(300):
            state.advance(gamma=1.0, f1_max_form=250.0, f1_max_dissoc=125.0)
        assert state.f1_formation == 250.0
        assert state.f1_dissociation == 125.0

    def test_reset(self):
        state = BoostState(step=100, f1_formation=100.0, f1_dissociation=50.0)
        state.reset()
        assert state.step == 0
        assert state.f1_formation == 0.0


class TestTotalBias:
    """Eq. 8: total bias energy and forces."""

    def test_single_formation_pair(self):
        positions = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        pair = PairBias(idx_a=0, idx_b=1, is_formation=True, r0=1.9)
        state = BoostState(step=100, f1_formation=100.0, f1_dissociation=50.0)
        params = TDBBParams(f2=10.0)

        energy, forces = total_bias([pair], positions, state, params)

        assert energy > 0.0
        # force on atom 0 should point toward atom 1 (positive x)
        assert forces[0, 0] > 0.0
        # Newton's third law
        np.testing.assert_allclose(forces[0], -forces[1])

    def test_single_dissociation_pair(self):
        positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        pair = PairBias(idx_a=0, idx_b=1, is_formation=False, r0=0.0)
        state = BoostState(step=100, f1_formation=100.0, f1_dissociation=50.0)
        params = TDBBParams(f2=10.0)

        energy, forces = total_bias([pair], positions, state, params)

        assert energy > 0.0
        # dissociation: force on atom 0 should push away from atom 1 (negative x)
        assert forces[0, 0] < 0.0
        np.testing.assert_allclose(forces[0], -forces[1])

    def test_no_pairs_returns_zero(self):
        positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        state = BoostState(step=100, f1_formation=100.0, f1_dissociation=50.0)
        params = TDBBParams()

        energy, forces = total_bias([], positions, state, params)

        assert energy == 0.0
        np.testing.assert_array_equal(forces, 0.0)

    def test_zero_boost_gives_zero(self):
        positions = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        pair = PairBias(idx_a=0, idx_b=1, is_formation=True, r0=1.5)
        state = BoostState()  # step=0, f1=0
        params = TDBBParams()

        energy, forces = total_bias([pair], positions, state, params)

        assert energy == 0.0
        np.testing.assert_array_equal(forces, 0.0)

    def test_force_magnitude_equals_dv_dr(self):
        """RF14: |force on atom| must equal |dV/dr| and point along the bond.

        With atoms on the x-axis, e_ij = +x, so forces[0] = dV/dr · e_ij.
        """
        r = 3.0
        positions = np.array([[0.0, 0.0, 0.0], [r, 0.0, 0.0]])
        r0, f1, f2 = 1.9, 100.0, 10.0
        pair = PairBias(idx_a=0, idx_b=1, is_formation=True, r0=r0)
        state = BoostState(step=100, f1_formation=f1, f1_dissociation=50.0)
        params = TDBBParams(f2=f2)

        _, forces = total_bias([pair], positions, state, params)

        dv_dr = formation_force_magnitude(np.array([r]), r0=r0, f1=f1, f2=f2)[0]
        # along +x, other components zero
        np.testing.assert_allclose(forces[0], [dv_dr, 0.0, 0.0], atol=1e-12)
        np.testing.assert_allclose(np.linalg.norm(forces[0]), abs(dv_dr), rtol=1e-12)

    def test_momentum_conserved_three_atoms(self):
        """Sum of bias forces is zero (Newton's third law over all pairs)."""
        positions = np.array([
            [0.0, 0.0, 0.0],
            [2.3, 0.1, 0.0],
            [0.2, 2.4, -0.3],
        ])
        pairs = [
            PairBias(idx_a=0, idx_b=1, is_formation=True, r0=1.9),
            PairBias(idx_a=1, idx_b=2, is_formation=False, r0=0.0),
            PairBias(idx_a=0, idx_b=2, is_formation=True, r0=2.0),
        ]
        state = BoostState(step=100, f1_formation=80.0, f1_dissociation=40.0)
        params = TDBBParams(f2=10.0)

        _, forces = total_bias(pairs, positions, state, params)

        np.testing.assert_allclose(forces.sum(axis=0), [0.0, 0.0, 0.0], atol=1e-12)

    def test_uses_minimum_image_across_pbc(self):
        """A pair separated across a periodic boundary uses the wrapped distance.

        atoms at x=0.5 and x=9.5 in a 10 Å box wrap to a 1.0 Å separation.
        With r0=1.0 the formation potential is 0 at the wrapped distance but
        ~f1 at the naive 9.0 Å distance, so the two paths are clearly distinct.
        """
        box = 10.0
        cell = np.diag([box, box, box]).astype(np.float64)
        positions = np.array([[0.5, 0.0, 0.0], [9.5, 0.0, 0.0]])
        pair = PairBias(idx_a=0, idx_b=1, is_formation=True, r0=1.0)
        state = BoostState(step=100, f1_formation=100.0, f1_dissociation=50.0)
        params = TDBBParams(f2=10.0)

        e_pbc, f_pbc = total_bias([pair], positions, state, params, cell=cell)
        e_naive, _ = total_bias([pair], positions, state, params, cell=None)

        # minimum image: r=1.0=r0 → V^f≈0; naive: r=9.0 → V^f≈f1
        assert e_pbc == pytest.approx(0.0, abs=1e-6)
        assert e_naive == pytest.approx(100.0, rel=1e-3)


class TestTotalBiasFiniteDifference:
    """RF14: 3D finite difference validation of TDBB forces.

    For every atom i and dimension d, verify:
        forces[i, d] ≈ -(E(x_i+h) - E(x_i-h)) / (2h)

    This is the gold-standard test for force correctness: it catches sign errors,
    missing factors, and incorrect chain-rule applications that unit tests on
    individual potentials might miss.
    """

    @staticmethod
    def _finite_diff_forces(pairs, positions, state, params, box, h=1e-5):
        """Compute numerical forces via central finite difference."""
        n_atoms = positions.shape[0]
        num_forces = np.zeros_like(positions)
        for i in range(n_atoms):
            for d in range(3):
                pos_plus = positions.copy()
                pos_plus[i, d] += h
                e_plus, _, _ = total_bias_fast(pairs, pos_plus, state, params, box=box)

                pos_minus = positions.copy()
                pos_minus[i, d] -= h
                e_minus, _, _ = total_bias_fast(pairs, pos_minus, state, params, box=box)

                num_forces[i, d] = -(e_plus - e_minus) / (2.0 * h)
        return num_forces

    def _check(self, pairs, positions, state, params, cell=None, rtol=1e-6):
        box = validated_box(cell)
        _, forces, _ = total_bias_fast(pairs, positions, state, params, box=box)
        num_forces = self._finite_diff_forces(pairs, positions, state, params, box)
        np.testing.assert_allclose(forces, num_forces, rtol=rtol, atol=1e-10)

    def test_formation_on_axis(self):
        positions = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        pair = PairBias(idx_a=0, idx_b=1, is_formation=True, r0=1.9)
        state = BoostState(step=100, f1_formation=100.0, f1_dissociation=50.0)
        params = TDBBParams(f2=10.0)
        self._check([pair], positions, state, params)

    def test_formation_off_axis(self):
        positions = np.array([[0.5, -0.3, 0.7], [2.1, 1.4, -0.2]])
        pair = PairBias(idx_a=0, idx_b=1, is_formation=True, r0=1.9)
        state = BoostState(step=150, f1_formation=150.0, f1_dissociation=75.0)
        params = TDBBParams(f2=10.0)
        self._check([pair], positions, state, params)

    def test_dissociation_off_axis(self):
        positions = np.array([[1.0, 0.5, -0.3], [1.8, 1.2, 0.4]])
        pair = PairBias(idx_a=0, idx_b=1, is_formation=False, r0=0.0)
        state = BoostState(step=80, f1_formation=80.0, f1_dissociation=80.0)
        params = TDBBParams(f2=10.0)
        self._check([pair], positions, state, params)

    def test_multiple_pairs_three_atoms(self):
        positions = np.array([
            [0.0, 0.0, 0.0],
            [2.3, 0.1, 0.5],
            [0.2, 2.4, -0.3],
        ])
        pairs = [
            PairBias(idx_a=0, idx_b=1, is_formation=True, r0=1.9),
            PairBias(idx_a=1, idx_b=2, is_formation=False, r0=0.0),
            PairBias(idx_a=0, idx_b=2, is_formation=True, r0=2.0),
        ]
        state = BoostState(step=100, f1_formation=80.0, f1_dissociation=40.0)
        params = TDBBParams(f2=10.0)
        self._check(pairs, positions, state, params)

    def test_with_pbc(self):
        cell = np.diag([10.0, 10.0, 10.0]).astype(np.float64)
        positions = np.array([[0.5, 5.0, 5.0], [9.5, 5.0, 5.0]])
        pair = PairBias(idx_a=0, idx_b=1, is_formation=True, r0=1.0)
        state = BoostState(step=50, f1_formation=200.0, f1_dissociation=100.0)
        params = TDBBParams(f2=10.0)
        self._check([pair], positions, state, params, cell=cell)

    def test_pbc_off_axis(self):
        cell = np.diag([8.0, 8.0, 8.0]).astype(np.float64)
        positions = np.array([
            [0.3, 7.5, 4.0],
            [7.8, 0.6, 4.2],
            [4.0, 4.0, 0.5],
        ])
        pairs = [
            PairBias(idx_a=0, idx_b=1, is_formation=True, r0=2.0),
            PairBias(idx_a=0, idx_b=2, is_formation=False, r0=0.0),
        ]
        state = BoostState(step=200, f1_formation=250.0, f1_dissociation=125.0)
        params = TDBBParams(f2=10.0)
        self._check(pairs, positions, state, params, cell=cell)

    @pytest.mark.parametrize('f2', [5.0, 10.0, 20.0])
    def test_varying_f2(self, f2):
        positions = np.array([[0.0, 0.0, 0.0], [2.5, 1.0, -0.5]])
        pair = PairBias(idx_a=0, idx_b=1, is_formation=True, r0=1.9)
        state = BoostState(step=100, f1_formation=100.0, f1_dissociation=50.0)
        params = TDBBParams(f2=f2)
        self._check([pair], positions, state, params)

    def test_near_r0_formation(self):
        """Near r0 the gradient is small but must still be accurate."""
        r0 = 1.932
        positions = np.array([[0.0, 0.0, 0.0], [r0 + 0.01, 0.005, -0.003]])
        pair = PairBias(idx_a=0, idx_b=1, is_formation=True, r0=r0)
        state = BoostState(step=250, f1_formation=250.0, f1_dissociation=125.0)
        params = TDBBParams(f2=10.0)
        self._check([pair], positions, state, params)


class TestTotalBiasFast:
    """Verify total_bias_fast matches total_bias for all existing scenarios."""

    def _compare(self, pairs, positions, state, params, cell=None):
        e_ref, f_ref = total_bias(pairs, positions, state, params, cell=cell)
        box = validated_box(cell)
        e_fast, f_fast, dists = total_bias_fast(pairs, positions, state, params, box=box)
        np.testing.assert_allclose(e_fast, e_ref, rtol=1e-12)
        np.testing.assert_allclose(f_fast, f_ref, atol=1e-12)
        assert len(dists) == len(pairs)
        return e_fast, f_fast, dists

    def test_single_formation(self):
        positions = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        pair = PairBias(idx_a=0, idx_b=1, is_formation=True, r0=1.9)
        state = BoostState(step=100, f1_formation=100.0, f1_dissociation=50.0)
        params = TDBBParams(f2=10.0)
        e, f, d = self._compare([pair], positions, state, params)
        assert e > 0.0
        assert d[0] == pytest.approx(3.0)

    def test_single_dissociation(self):
        positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        pair = PairBias(idx_a=0, idx_b=1, is_formation=False, r0=0.0)
        state = BoostState(step=100, f1_formation=100.0, f1_dissociation=50.0)
        params = TDBBParams(f2=10.0)
        self._compare([pair], positions, state, params)

    def test_no_pairs(self):
        positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        state = BoostState(step=100, f1_formation=100.0, f1_dissociation=50.0)
        params = TDBBParams()
        e, f, d = self._compare([], positions, state, params)
        assert e == 0.0
        assert d == []

    def test_zero_boost(self):
        positions = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        pair = PairBias(idx_a=0, idx_b=1, is_formation=True, r0=1.5)
        state = BoostState()
        params = TDBBParams()
        self._compare([pair], positions, state, params)

    def test_momentum_conserved_three_atoms(self):
        positions = np.array([
            [0.0, 0.0, 0.0],
            [2.3, 0.1, 0.0],
            [0.2, 2.4, -0.3],
        ])
        pairs = [
            PairBias(idx_a=0, idx_b=1, is_formation=True, r0=1.9),
            PairBias(idx_a=1, idx_b=2, is_formation=False, r0=0.0),
            PairBias(idx_a=0, idx_b=2, is_formation=True, r0=2.0),
        ]
        state = BoostState(step=100, f1_formation=80.0, f1_dissociation=40.0)
        params = TDBBParams(f2=10.0)
        _, forces, _ = self._compare(pairs, positions, state, params)
        np.testing.assert_allclose(forces.sum(axis=0), [0.0, 0.0, 0.0], atol=1e-12)

    def test_pbc_minimum_image(self):
        box_l = 10.0
        cell = np.diag([box_l, box_l, box_l]).astype(np.float64)
        positions = np.array([[0.5, 0.0, 0.0], [9.5, 0.0, 0.0]])
        pair = PairBias(idx_a=0, idx_b=1, is_formation=True, r0=1.0)
        state = BoostState(step=100, f1_formation=100.0, f1_dissociation=50.0)
        params = TDBBParams(f2=10.0)
        e, _, d = self._compare([pair], positions, state, params, cell=cell)
        assert e == pytest.approx(0.0, abs=1e-6)
        assert d[0] == pytest.approx(1.0)
