"""Tests for Verlet and Langevin integrators."""
import numpy as np
import pytest

from src.integrators.verlet import VelocityVerletIntegrator
from src.integrators.langevin import LangevinIntegrator, LangevinParams
from src.units import FORCE_CONV, KB


class TestWrapPositions:
    """Position wrapping via integrators (T-D)."""

    def test_verlet_wraps_with_cell(self):
        rng = np.random.default_rng(42)
        integrator = VelocityVerletIntegrator()
        cell = np.diag([10.0, 10.0, 10.0])

        # Atom at 9.5 with velocity 1.0 A/fs over dt=1.0 -> new pos = 10.5 -> should wrap to 0.5
        pos = np.array([[9.5, 0.0, 0.0]])
        vel = np.array([[1.0, 0.0, 0.0]])
        forces = np.zeros((1, 3))
        integrator.pre_force(pos, vel, forces, None, dt=1.0, rng=rng, cell=cell)

        assert 0.0 <= pos[0, 0] < 10.0, f'Expected wrapped position, got {pos[0, 0]}'

    def test_verlet_no_wrap_without_cell(self):
        rng = np.random.default_rng(42)
        integrator = VelocityVerletIntegrator()

        pos = np.array([[9.5, 0.0, 0.0]])
        vel = np.array([[1.0, 0.0, 0.0]])
        forces = np.zeros((1, 3))
        integrator.pre_force(pos, vel, forces, None, dt=1.0, rng=rng, cell=None)

        assert pos[0, 0] > 10.0, 'Without cell, position should not be wrapped'

    def test_langevin_wraps_with_cell(self):
        params = LangevinParams(temperature_K=300.0, friction_per_fs=0.001)
        integrator = LangevinIntegrator(params)
        rng = np.random.default_rng(42)
        cell = np.diag([5.0, 5.0, 5.0])

        # Run 100 steps from zero velocities — all positions must stay inside [0, 5)
        pos = np.array([[4.9, 0.0, 0.0], [0.1, 0.0, 0.0]])
        vel = np.zeros((2, 3))
        forces = np.zeros((2, 3))
        masses = np.array([12.0, 12.0])

        for _ in range(100):
            integrator.pre_force(pos, vel, forces, masses, dt=0.25, rng=rng, cell=cell)
            integrator.post_force(vel, forces, masses, dt=0.25)

        assert np.all(pos[:, 0] >= 0.0) and np.all(pos[:, 0] < 5.0)

    def test_all_positions_in_cell_after_many_steps(self):
        rng = np.random.default_rng(7)
        integrator = VelocityVerletIntegrator()
        box = 8.0
        cell = np.diag([box, box, box])

        n = 20
        pos = rng.uniform(0, box, (n, 3))
        vel = rng.standard_normal((n, 3)) * 0.1
        forces = np.zeros((n, 3))

        for _ in range(500):
            integrator.pre_force(pos, vel, forces, None, dt=0.25, rng=rng, cell=cell)
            integrator.post_force(vel, forces, None, dt=0.25)

        assert np.all(pos >= 0.0) and np.all(pos < box)


class TestVelocityVerlet:

    def test_free_particle(self):
        pos = np.array([[0.0, 0.0, 0.0]])
        vel = np.array([[1.0, 0.0, 0.0]])
        forces = np.array([[0.0, 0.0, 0.0]])
        rng = np.random.default_rng(42)
        integrator = VelocityVerletIntegrator()

        integrator.pre_force(pos, vel, forces, None, dt=1.0, rng=rng)
        integrator.post_force(vel, forces, None, dt=1.0)

        np.testing.assert_allclose(pos[0, 0], 1.0)
        np.testing.assert_allclose(vel[0, 0], 1.0)

    def test_constant_force(self):
        pos = np.array([[0.0, 0.0, 0.0]])
        vel = np.array([[0.0, 0.0, 0.0]])
        forces = np.array([[2.0, 0.0, 0.0]])
        rng = np.random.default_rng(42)
        integrator = VelocityVerletIntegrator()

        integrator.pre_force(pos, vel, forces, None, dt=1.0, rng=rng)
        integrator.post_force(vel, forces, None, dt=1.0)

        a = FORCE_CONV * 2.0
        np.testing.assert_allclose(pos[0, 0], 0.5 * a, rtol=1e-10)
        np.testing.assert_allclose(vel[0, 0], a, rtol=1e-10)

    def test_proper_vv_with_different_forces(self):
        """pre_force uses old forces, post_force uses new forces."""
        pos = np.array([[0.0, 0.0, 0.0]])
        vel = np.array([[0.0, 0.0, 0.0]])
        old_forces = np.array([[2.0, 0.0, 0.0]])
        new_forces = np.array([[4.0, 0.0, 0.0]])
        rng = np.random.default_rng(42)
        integrator = VelocityVerletIntegrator()

        integrator.pre_force(pos, vel, old_forces, None, dt=1.0, rng=rng)
        integrator.post_force(vel, new_forces, None, dt=1.0)

        a_old = FORCE_CONV * 2.0
        a_new = FORCE_CONV * 4.0
        np.testing.assert_allclose(vel[0, 0], 0.5 * a_old + 0.5 * a_new)
        np.testing.assert_allclose(pos[0, 0], 0.5 * a_old)

    def test_mass_affects_acceleration(self):
        rng = np.random.default_rng(42)
        integrator = VelocityVerletIntegrator()

        pos_light = np.array([[0.0, 0.0, 0.0]])
        vel_light = np.array([[0.0, 0.0, 0.0]])
        masses_light = np.array([1.0])
        forces = np.array([[1.0, 0.0, 0.0]])
        integrator.pre_force(pos_light, vel_light, forces.copy(), masses_light, dt=1.0, rng=rng)
        integrator.post_force(vel_light, forces.copy(), masses_light, dt=1.0)

        pos_heavy = np.array([[0.0, 0.0, 0.0]])
        vel_heavy = np.array([[0.0, 0.0, 0.0]])
        masses_heavy = np.array([12.0])
        integrator.pre_force(pos_heavy, vel_heavy, forces.copy(), masses_heavy, dt=1.0, rng=rng)
        integrator.post_force(vel_heavy, forces.copy(), masses_heavy, dt=1.0)

        assert pos_light[0, 0] > pos_heavy[0, 0]
        assert vel_light[0, 0] > vel_heavy[0, 0]

    def test_backward_compat_none_masses(self):
        rng = np.random.default_rng(42)
        integrator = VelocityVerletIntegrator()

        pos_none = np.array([[0.0, 0.0, 0.0]])
        vel_none = np.array([[0.0, 0.0, 0.0]])
        forces = np.array([[1.0, 0.0, 0.0]])
        integrator.pre_force(pos_none, vel_none, forces.copy(), None, dt=0.5, rng=rng)
        integrator.post_force(vel_none, forces.copy(), None, dt=0.5)

        pos_unit = np.array([[0.0, 0.0, 0.0]])
        vel_unit = np.array([[0.0, 0.0, 0.0]])
        masses_unit = np.array([1.0])
        integrator.pre_force(pos_unit, vel_unit, forces.copy(), masses_unit, dt=0.5, rng=rng)
        integrator.post_force(vel_unit, forces.copy(), masses_unit, dt=0.5)

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
                integrator.pre_force(pos, vel, forces, None, dt=0.25, rng=rng)
                integrator.post_force(vel, forces, None, dt=0.25)
            return pos.copy(), vel.copy()

        p1, v1 = run_one()
        p2, v2 = run_one()
        np.testing.assert_array_equal(p1, p2)
        np.testing.assert_array_equal(v1, v2)

    def test_temperature_equilibration(self):
        """Kinetic temperature should stabilize near target with correct units."""
        params = LangevinParams(temperature_K=300.0, friction_per_fs=0.05)
        integrator = LangevinIntegrator(params)
        rng = np.random.default_rng(7)

        n_atoms = 100
        pos = rng.uniform(0, 10, size=(n_atoms, 3))
        vel = np.zeros((n_atoms, 3))
        masses = np.full(n_atoms, 12.0)
        forces = np.zeros((n_atoms, 3))

        for _ in range(2000):
            integrator.pre_force(pos, vel, forces, masses, dt=0.25, rng=rng)
            integrator.post_force(vel, forces, masses, dt=0.25)

        ke_amu = 0.5 * np.sum(masses[:, np.newaxis] * vel ** 2)
        ke_kcal = ke_amu / FORCE_CONV
        T_kinetic = (2.0 / 3.0) * ke_kcal / (n_atoms * KB)

        assert 150 < T_kinetic < 450, f'T_kinetic = {T_kinetic:.1f} K, expected ~300 K'

    def test_mass_dependence(self):
        params = LangevinParams(temperature_K=300.0, friction_per_fs=0.01)
        integrator = LangevinIntegrator(params)

        rng = np.random.default_rng(42)
        pos = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        vel = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        masses = np.array([1.0, 100.0])
        forces = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])

        integrator.pre_force(pos, vel, forces, masses, dt=0.25, rng=rng)
        integrator.post_force(vel, forces, masses, dt=0.25)

        assert abs(pos[0, 0]) > abs(pos[1, 0])


class TestVelocityVerletConservation:
    """RF21: symplectic energy conservation and time-reversibility under a
    conservative (harmonic) force. Free-particle / constant-force tests keep the
    force constant, so they cannot catch a regression that uses the new force in
    pre_force or applies a half-kick twice; these do."""

    K = 10.0          # spring constant, kcal/(mol·Å²)
    DT = 0.25         # fs
    MASS = 12.0       # amu

    def _energy(self, pos, vel):
        ke = 0.5 * self.MASS * float(np.sum(vel ** 2)) / FORCE_CONV
        pe = 0.5 * self.K * float(np.sum(pos ** 2))
        return ke + pe

    def test_energy_conserved_harmonic(self):
        integrator = VelocityVerletIntegrator()
        rng = np.random.default_rng(0)
        masses = np.array([self.MASS])
        pos = np.array([[1.0, 0.0, 0.0]])
        vel = np.zeros((1, 3))
        forces = -self.K * pos
        e0 = self._energy(pos, vel)

        max_dev = 0.0
        for _ in range(4000):  # ~12 oscillation periods
            integrator.pre_force(pos, vel, forces, masses, dt=self.DT, rng=rng)
            forces = -self.K * pos
            integrator.post_force(vel, forces, masses, dt=self.DT)
            max_dev = max(max_dev, abs(self._energy(pos, vel) - e0))

        # symplectic: bounded oscillation around e0, no secular drift
        assert max_dev / e0 < 1e-3, f'energy drift {max_dev / e0:.2e}'

    def test_time_reversibility_harmonic(self):
        integrator = VelocityVerletIntegrator()
        rng = np.random.default_rng(0)
        masses = np.array([self.MASS])
        pos = np.array([[1.0, 0.0, 0.0]])
        vel = np.array([[0.2, -0.1, 0.05]])
        x0, v0 = pos.copy(), vel.copy()
        forces = -self.K * pos

        for _ in range(500):
            integrator.pre_force(pos, vel, forces, masses, dt=self.DT, rng=rng)
            forces = -self.K * pos
            integrator.post_force(vel, forces, masses, dt=self.DT)

        vel[:] = -vel            # reverse time
        forces = -self.K * pos
        for _ in range(500):
            integrator.pre_force(pos, vel, forces, masses, dt=self.DT, rng=rng)
            forces = -self.K * pos
            integrator.post_force(vel, forces, masses, dt=self.DT)

        np.testing.assert_allclose(pos, x0, atol=1e-6)
        np.testing.assert_allclose(vel, -v0, atol=1e-6)
