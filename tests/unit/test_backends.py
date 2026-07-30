"""Tests for calculator backends."""
import logging

import numpy as np
import pytest

from kagome.backends.toy import ToyCalculator


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
        from kagome.backends.orb_backend import OrbCalculatorAdapter
        adapter = OrbCalculatorAdapter(model=object(), atoms_adapter=object())
        assert adapter.supports_spin is True

    def test_orb_periodic_guard_raises_without_nvalchemiops(self, monkeypatch):
        pytest.importorskip('ase')
        import importlib.util
        from kagome.backends.orb_backend import OrbCalculatorAdapter

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
        from kagome.backends.ase_adapter import ASECalculatorAdapter

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
        from kagome.backends.mace_backend import create_mace_calculator

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


class TestLocalModelPaths:
    """B9: _LOCAL_MODELS must resolve to <repo-root>/models, not <repo-root>/src.

    parents[2] pointed at src/ (old layout), so local checkpoints were never
    found and the loaders silently fell back to downloaded pretrained weights.
    We verify path *structure* (not file existence) so the check passes in CI
    environments without the weight files.
    """

    def _assert_under_models_root(self, path):
        from pathlib import Path
        path = Path(path)
        # <root>/models/<file>
        assert path.parent.name == 'models', f'{path} not under a models/ dir'
        root = path.parent.parent
        # The repo root must contain src/ (i.e. we did not resolve into src/).
        assert (root / 'src').is_dir(), (
            f'model path root {root} has no src/ dir — _PROJECT_ROOT likely '
            f'points into src/ (the parents[2] B9 bug)'
        )
        assert root.name != 'src', f'model path root resolved into src/: {root}'

    def test_orb_local_models_under_repo_models(self):
        from kagome.backends import orb_backend
        assert orb_backend._LOCAL_MODELS
        for path in orb_backend._LOCAL_MODELS.values():
            self._assert_under_models_root(path)

    def test_mace_local_models_under_repo_models(self):
        from kagome.backends import mace_backend
        assert mace_backend._LOCAL_MODELS
        for path in mace_backend._LOCAL_MODELS.values():
            self._assert_under_models_root(path)


class TestOrbBackend:

    @pytest.fixture
    def _skip_no_orb(self):
        pytest.importorskip('orb_models')

    @pytest.mark.usefixtures('_skip_no_orb')
    @pytest.mark.slow
    def test_orb_compute(self):
        from kagome.backends.orb_backend import create_orb_calculator

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


class TestOrbTorchEnvConfig:
    """compile=True must not set TORCHDYNAMO_DISABLE=1, which disables
    torch.compile process-wide (decisions.md 2026-07-14). compile=False keeps
    the Windows cl.exe workaround."""

    def test_no_compile_sets_dynamo_disable(self, monkeypatch):
        import os
        from kagome.backends.orb_backend import _configure_torch_env
        monkeypatch.delenv('TORCHDYNAMO_DISABLE', raising=False)
        _configure_torch_env(compile=False)
        assert os.environ.get('TORCHDYNAMO_DISABLE') == '1'

    def test_compile_leaves_dynamo_enabled(self, monkeypatch):
        import os
        from kagome.backends.orb_backend import _configure_torch_env
        monkeypatch.delenv('TORCHDYNAMO_DISABLE', raising=False)
        _configure_torch_env(compile=True)
        assert 'TORCHDYNAMO_DISABLE' not in os.environ

    def test_compile_warns_when_dynamo_already_disabled(self, monkeypatch, caplog):
        from kagome.backends.orb_backend import _configure_torch_env
        monkeypatch.setenv('TORCHDYNAMO_DISABLE', '1')
        with caplog.at_level(logging.WARNING):
            _configure_torch_env(compile=True)
        assert any('no-op' in r.message for r in caplog.records)


class TestOrbEmptyCacheGating:
    """Per-step torch.cuda.empty_cache() costs ~9% CPU at paper scale (py-spy,
    decisions.md 2026-07-14); empty_cache=False must skip it while the default
    keeps the 16 GB fragmentation mitigation (decisions.md 2026-06-15)."""

    def _make_adapter(self, empty_cache: bool):
        pytest.importorskip('torch')
        pytest.importorskip('ase')
        import torch
        from kagome.backends.orb_backend import OrbCalculatorAdapter

        class _FakeModel:
            def parameters(self):
                yield torch.zeros(1)

            def __call__(self, batch):
                return {
                    'energy': torch.tensor(-1.0),
                    'grad_forces': torch.zeros((2, 3)),
                }

        class _FakeBatch:
            positions = None

        class _FakeAtomsAdapter:
            def from_ase_atoms(self, atoms, device=None):
                return _FakeBatch()

        # device='cuda' exercises the cleanup branch; the fake adapter/model
        # never touch a real GPU.
        return OrbCalculatorAdapter(
            _FakeModel(), _FakeAtomsAdapter(),
            device='cuda', empty_cache=empty_cache,
        )

    def _compute(self, calc):
        positions = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
        return calc.compute(positions, ['C', 'C'])

    def test_empty_cache_called_by_default(self, monkeypatch):
        torch = pytest.importorskip('torch')
        calls = []
        monkeypatch.setattr(torch.cuda, 'empty_cache', lambda: calls.append(1))
        calc = self._make_adapter(empty_cache=True)
        self._compute(calc)
        assert calls == [1]

    def test_empty_cache_skipped_when_disabled(self, monkeypatch):
        torch = pytest.importorskip('torch')
        calls = []
        monkeypatch.setattr(torch.cuda, 'empty_cache', lambda: calls.append(1))
        calc = self._make_adapter(empty_cache=False)
        energy, forces = self._compute(calc)
        assert calls == []
        assert forces.shape == (2, 3)


class TestAimnetBackend:
    """AIMNet2 backend (code+weights MIT, verified 2026-06-25; specs/
    dependency-license-matrix.md). Complements OrbMol-v2 for open-shell/
    radical chemistry (specs/decisions.md 2026-06-25)."""

    @pytest.fixture
    def _skip_no_aimnet(self):
        pytest.importorskip('aimnet')

    @pytest.mark.usefixtures('_skip_no_aimnet')
    @pytest.mark.slow
    def test_aimnet_compute(self):
        from kagome.backends.aimnet_backend import create_aimnet_calculator

        calc = create_aimnet_calculator(model='aimnet2-nse', device='cpu', spin=2)
        positions = np.array([
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
        ])
        energy, forces = calc.compute(positions, ['C', 'C'])

        assert isinstance(energy, float)
        assert forces.shape == (2, 3)
        assert 'aimnet' in calc.name

    def test_create_aimnet_calculator_raises_without_aimnet_package(self, monkeypatch):
        """Error path: missing 'aimnet' package must raise a clear ImportError
        (mirrors the mace-torch / orb-models missing-dependency messages)."""
        import sys

        monkeypatch.setitem(sys.modules, 'aimnet', None)
        monkeypatch.setitem(sys.modules, 'aimnet.calculators', None)

        from kagome.backends.aimnet_backend import create_aimnet_calculator

        with pytest.raises(ImportError, match='pip install aimnet'):
            create_aimnet_calculator()

    @pytest.fixture
    def _skip_no_ase(self):
        pytest.importorskip('ase')

    @pytest.mark.usefixtures('_skip_no_ase')
    def test_adapter_compute_with_fake_ase_calculator(self):
        """Adapter generation + compute() path, exercised with a fake
        ASE-calculator-like stand-in so this test does not require the real
        aimnet package (only ase, which is already exercised elsewhere)."""
        from ase.calculators.lj import LennardJones

        from kagome.backends.aimnet_backend import AimnetCalculatorAdapter

        class _FakeAimnetAseCalc(LennardJones):
            def __init__(self):
                super().__init__()
                self.mult_calls: list[int] = []
                self.charge_calls: list[int] = []

            def set_mult(self, mult):
                self.mult_calls.append(mult)

            def set_charge(self, charge):
                self.charge_calls.append(charge)

        fake_calc = _FakeAimnetAseCalc()
        adapter = AimnetCalculatorAdapter(fake_calc, name='aimnet-test', spin=3, charge=0)

        positions = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        energy, forces = adapter.compute(positions, ['Ar', 'Ar'])

        assert isinstance(energy, float)
        assert forces.shape == (2, 3)
        assert fake_calc.mult_calls == [3]
        assert fake_calc.charge_calls == [0]
        assert adapter.supports_spin is True
        assert adapter.model_id == 'aimnet-test'
        assert adapter.magnetic_moments() is None  # LJ has no magnetic moments
