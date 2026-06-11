"""Tests for TDBB potentials and forces.

Validates against paper equations 2-5, 8.
"""
import numpy as np
import pytest

from src.boost.tdbb import (
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
)


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
