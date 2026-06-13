"""Tests for Maxwell-Boltzmann velocity initialization."""
import numpy as np
import pytest

from src.integrators.init_velocities import instant_temperature_K, maxwell_boltzmann_velocities
from src.units import FORCE_CONV, KB


class TestMaxwellBoltzmannVelocities:

    def _masses(self, n: int = 10, element: str = 'C') -> np.ndarray:
        return np.full(n, 12.011)

    def test_shape(self):
        masses = self._masses(8)
        rng = np.random.default_rng(42)
        v = maxwell_boltzmann_velocities(masses, 300.0, rng)
        assert v.shape == (8, 3)

    def test_temperature_within_tolerance(self):
        # Statistical test: N=100 atoms, 5 independent seeds, all within ±30%
        rng = np.random.default_rng(0)
        masses = np.full(100, 12.011)
        target_T = 500.0
        for seed in range(5):
            rng = np.random.default_rng(seed)
            v = maxwell_boltzmann_velocities(masses, target_T, rng)
            T_inst = instant_temperature_K(v, masses)
            assert abs(T_inst - target_T) / target_T < 0.3, (
                f'Temperature {T_inst:.1f} K deviates more than 30% from target {target_T:.1f} K'
            )

    def test_com_velocity_removed(self):
        masses = self._masses(20)
        rng = np.random.default_rng(7)
        v = maxwell_boltzmann_velocities(masses, 300.0, rng, remove_com=True)
        total_mass = np.sum(masses)
        v_com = np.sum(masses[:, np.newaxis] * v, axis=0) / total_mass
        np.testing.assert_allclose(v_com, np.zeros(3), atol=1e-12)

    def test_com_not_removed_when_false(self):
        masses = self._masses(20)
        rng = np.random.default_rng(99)
        v = maxwell_boltzmann_velocities(masses, 300.0, rng, remove_com=False)
        total_mass = np.sum(masses)
        v_com = np.sum(masses[:, np.newaxis] * v, axis=0) / total_mass
        # With remove_com=False the COM velocity is generally nonzero
        assert np.linalg.norm(v_com) > 1e-10

    def test_single_atom(self):
        masses = np.array([12.011])
        rng = np.random.default_rng(1)
        v = maxwell_boltzmann_velocities(masses, 300.0, rng, remove_com=False)
        assert v.shape == (1, 3)

    def test_reproducible_with_same_seed(self):
        masses = self._masses(10)
        v1 = maxwell_boltzmann_velocities(masses, 300.0, np.random.default_rng(42))
        v2 = maxwell_boltzmann_velocities(masses, 300.0, np.random.default_rng(42))
        np.testing.assert_array_equal(v1, v2)


class TestInstantTemperature:

    def test_zero_velocities(self):
        masses = np.array([12.011, 12.011])
        v = np.zeros((2, 3))
        assert instant_temperature_K(v, masses) == pytest.approx(0.0)

    def test_no_atoms(self):
        masses = np.array([], dtype=float)
        v = np.zeros((0, 3))
        assert instant_temperature_K(v, masses) == pytest.approx(0.0)

    def test_known_temperature(self):
        # Single atom, 1D: KE = 0.5*m*v^2/FORCE_CONV, T = 2*KE/(3*1*KB)
        # Set v so T = 300 K exactly for 1 atom
        # sigma = sqrt(KB * T * FORCE_CONV / m) per DOF
        m = 12.011
        T = 300.0
        sigma = (KB * T * FORCE_CONV / m) ** 0.5
        # Assign equal velocity to all 3 DOF — then average KE = 3 * 0.5 * m * sigma^2 / FC
        v = np.array([[sigma, sigma, sigma]])
        masses = np.array([m])
        T_measured = instant_temperature_K(v, masses)
        assert T_measured == pytest.approx(T, rel=1e-6)
