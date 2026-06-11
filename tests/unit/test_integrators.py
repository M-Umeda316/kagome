"""Tests for Verlet and Langevin integrators."""
import numpy as np
import pytest

from src.integrators.verlet import VelocityVerletIntegrator
from src.integrators.langevin import LangevinIntegrator, LangevinParams, KB_KCAL_MOL_K


class TestVelocityVerlet:

    def test_free_particle(self):
        pos = np.array([[0.0, 0.0, 0.0]])
        vel = np.array([[1.0, 0.0, 0.0]])
        forces = np.array([[0.0, 0.0, 0.0]])
        rng = np.random.default_rng(42)
        integrator = VelocityVerletIntegrator()

        integrator.step(pos, vel, forces, None, dt=1.0, rng=rng)

        np.testing.assert_allclose(pos[0, 0], 1.0)
        np.testing.assert_allclose(vel[0, 0], 1.0)

    def test_constant_force(self):
        pos = np.array([[0.0, 0.0, 0.0]])
        vel = np.array([[0.0, 0.0, 0.0]])
        forces = np.array([[2.0, 0.0, 0.0]])
        rng = np.random.default_rng(42)
        integrator = VelocityVerletIntegrator()

        integrator.step(pos, vel, forces, None, dt=1.0, rng=rng)

        # With unit mass: a=2, v_half=1, x=1, v_final=2
        np.testing.assert_allclose(pos[0, 0], 1.0)
        np.testing.assert_allclose(vel[0, 0], 2.0)

    def test_mass_affects_acceleration(self):
        """Heavier atoms should move less under the same force."""
        rng = np.random.default_rng(42)
        integrator = VelocityVerletIntegrator()

        # light atom (mass=1)
        pos_light = np.array([[0.0, 0.0, 0.0]])
        vel_light = np.array([[0.0, 0.0, 0.0]])
        masses_light = np.array([1.0])
        forces = np.array([[1.0, 0.0, 0.0]])
        integrator.step(pos_light, vel_light, forces.copy(), masses_light, dt=1.0, rng=rng)

        # heavy atom (mass=12)
        pos_heavy = np.array([[0.0, 0.0, 0.0]])
        vel_heavy = np.array([[0.0, 0.0, 0.0]])
        masses_heavy = np.array([12.0])
        integrator.step(pos_heavy, vel_heavy, forces.copy(), masses_heavy, dt=1.0, rng=rng)

        assert pos_light[0, 0] > pos_heavy[0, 0]
        assert vel_light[0, 0] > vel_heavy[0, 0]

    def test_backward_compat_none_masses(self):
        """None masses should behave like unit masses."""
        rng = np.random.default_rng(42)
        integrator = VelocityVerletIntegrator()

        pos_none = np.array([[0.0, 0.0, 0.0]])
        vel_none = np.array([[0.0, 0.0, 0.0]])
        forces = np.array([[1.0, 0.0, 0.0]])
        integrator.step(pos_none, vel_none, forces.copy(), None, dt=0.5, rng=rng)

        pos_unit = np.array([[0.0, 0.0, 0.0]])
        vel_unit = np.array([[0.0, 0.0, 0.0]])
        masses_unit = np.array([1.0])
        integrator.step(pos_unit, vel_unit, forces.copy(), masses_unit, dt=0.5, rng=rng)

        np.testing.assert_allclose(pos_none, pos_unit)
        np.testing.assert_allclose(vel_none, vel_unit)


class TestLangevin:

    def test_deterministic_with_seed(self):
        params = LangevinParams(temperature_K=300.0, friction_per_fs=0.01)
        integrator = LangevinIntegrator(params)

        def run_one():
            rng = np.random.default_rng(42)
            pos = np.array([[0.0, 0.0, 0.0]])
            vel = np.array([[0.0, 0.0, 0.0]])
            forces = np.array([[0.0, 0.0, 0.0]])
            for _ in range(10):
                integrator.step(pos, vel, forces, None, dt=0.25, rng=rng)
            return pos.copy(), vel.copy()

        p1, v1 = run_one()
        p2, v2 = run_one()
        np.testing.assert_array_equal(p1, p2)
        np.testing.assert_array_equal(v1, v2)

    def test_temperature_equilibration(self):
        """Check that kinetic temperature stabilizes near target for many particles."""
        params = LangevinParams(temperature_K=300.0, friction_per_fs=0.05)
        integrator = LangevinIntegrator(params)
        rng = np.random.default_rng(7)

        n_atoms = 100
        pos = rng.uniform(0, 10, size=(n_atoms, 3))
        vel = np.zeros((n_atoms, 3))
        masses = np.full(n_atoms, 12.0)
        forces = np.zeros((n_atoms, 3))

        for _ in range(2000):
            integrator.step(pos, vel, forces, masses, dt=0.25, rng=rng)

        # Kinetic temperature: T = (2/3) * KE / (n * kB)
        ke = 0.5 * np.sum(masses[:, np.newaxis] * vel ** 2)
        T_kinetic = (2.0 / 3.0) * ke / (n_atoms * KB_KCAL_MOL_K)

        # Allow ±50% tolerance for statistical fluctuations
        assert 150 < T_kinetic < 450, f'T_kinetic = {T_kinetic:.1f} K, expected ~300 K'

    def test_mass_dependence(self):
        params = LangevinParams(temperature_K=300.0, friction_per_fs=0.01)
        integrator = LangevinIntegrator(params)

        rng = np.random.default_rng(42)
        pos = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        vel = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        masses = np.array([1.0, 100.0])
        forces = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])

        integrator.step(pos, vel, forces, masses, dt=0.25, rng=rng)

        assert abs(pos[0, 0]) > abs(pos[1, 0])
