"""Tests for calculator backends."""
import logging

import numpy as np
import pytest

from src.backends.toy import ToyCalculator


class TestToyCalculator:

    def test_energy_and_forces_shape(self):
        calc = ToyCalculator()
        positions = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        energy, forces = calc.compute(positions, ['C', 'C'])
        assert isinstance(energy, float)
        assert forces.shape == (2, 3)

    def test_newtons_third_law(self):
        calc = ToyCalculator()
        positions = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        _, forces = calc.compute(positions, ['C', 'C'])
        np.testing.assert_allclose(forces[0], -forces[1], atol=1e-12)

    def test_name(self):
        assert ToyCalculator().name == 'toy'

    def test_model_id_defaults_to_name(self):
        assert ToyCalculator().model_id == 'toy'


class TestSpinAndPeriodicGuards:
    """RF20: spin capability is explicit; spin-agnostic backends warn; periodic
    OrbMol-v2 fails clearly when nvalchemiops is unavailable."""

    def test_toy_supports_spin_false(self):
        assert ToyCalculator().supports_spin is False

    def test_set_spin_warns_when_unsupported(self, caplog):
        calc = ToyCalculator()
        with caplog.at_level(logging.WARNING):
            calc.set_spin(3)
        assert 'ignores spin' in caplog.text

    def test_orb_supports_spin_true(self):
        pytest.importorskip('ase')
        from src.backends.orb_backend import OrbCalculatorAdapter
        adapter = OrbCalculatorAdapter(model=object(), atoms_adapter=object())
        assert adapter.supports_spin is True

    def test_orb_periodic_guard_raises_without_nvalchemiops(self, monkeypatch):
        pytest.importorskip('ase')
        import importlib.util
        from src.backends.orb_backend import OrbCalculatorAdapter

        adapter = OrbCalculatorAdapter(model=object(), atoms_adapter=object())
        real = importlib.util.find_spec
        monkeypatch.setattr(
            importlib.util, 'find_spec',
            lambda name, *a, **k: None if name == 'nvalchemiops' else real(name, *a, **k),
        )
        cell = np.diag([10.0, 10.0, 10.0]).astype(float)
        with pytest.raises(RuntimeError, match='nvalchemiops'):
            adapter.compute(np.zeros((1, 3)), ['C'], cell=cell)


class TestASEAdapter:

    @pytest.fixture
    def _skip_no_ase(self):
        pytest.importorskip('ase')

    @pytest.mark.usefixtures('_skip_no_ase')
    def test_adapter_with_lj(self):
        from ase.calculators.lj import LennardJones
        from src.backends.ase_adapter import ASECalculatorAdapter

        calc = ASECalculatorAdapter(LennardJones(), name='ase-lj')
        positions = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        energy, forces = calc.compute(positions, ['Ar', 'Ar'])

        assert isinstance(energy, float)
        assert forces.shape == (2, 3)
        assert calc.name == 'ase-lj'


class TestMACEBackend:

    @pytest.fixture
    def _skip_no_mace(self):
        pytest.importorskip('mace')

    @pytest.mark.usefixtures('_skip_no_mace')
    @pytest.mark.slow
    def test_mace_compute(self):
        from src.backends.mace_backend import create_mace_calculator

        calc = create_mace_calculator(model='small', device='cpu')
        positions = np.array([
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
        ])
        energy, forces = calc.compute(positions, ['C', 'C'])

        assert isinstance(energy, float)
        assert forces.shape == (2, 3)
        # C-C at 1.5 Å should have repulsive forces
        np.testing.assert_allclose(forces[0], -forces[1], atol=1e-6)
        assert 'mace' in calc.name


class TestOrbBackend:

    @pytest.fixture
    def _skip_no_orb(self):
        pytest.importorskip('orb_models')

    @pytest.mark.usefixtures('_skip_no_orb')
    @pytest.mark.slow
    def test_orb_compute(self):
        from src.backends.orb_backend import create_orb_calculator

        calc = create_orb_calculator(device='cpu')
        positions = np.array([
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
        ])
        energy, forces = calc.compute(positions, ['C', 'C'])

        assert isinstance(energy, float)
        assert forces.shape == (2, 3)
        np.testing.assert_allclose(forces[0], -forces[1], atol=1e-4)
        assert 'orb' in calc.name
