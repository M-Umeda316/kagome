"""Tests for RunManifest save/load and effective parameter recording (RF1)."""
import json

import numpy as np
import pytest

from kagome.workflows.manifest import RunManifest, _normalize_value


class TestNormalizeValue:

    def test_numpy_float(self):
        assert isinstance(_normalize_value(np.float64(1.5)), float)

    def test_numpy_int(self):
        assert isinstance(_normalize_value(np.int64(3)), int)

    def test_nested_dict(self):
        d = {'a': np.float64(1.0), 'b': {'c': np.int32(2)}}
        result = _normalize_value(d)
        assert result == {'a': 1.0, 'b': {'c': 2}}
        assert isinstance(result['a'], float)
        assert isinstance(result['b']['c'], int)

    def test_list(self):
        result = _normalize_value([np.float64(1.0), np.int64(2)])
        assert result == [1.0, 2]

    def test_plain_values_unchanged(self):
        assert _normalize_value(42) == 42
        assert _normalize_value('hello') == 'hello'
        assert _normalize_value(3.14) == 3.14


class TestRunManifest:

    def test_save_and_load(self, tmp_path):
        manifest = RunManifest(
            config_path='configs/boost/paper_faithful.yaml',
            seed=7,
            backend='toy',
            output_dir=str(tmp_path),
        )
        out = tmp_path / 'manifest.json'
        manifest.save(out)
        assert out.exists()
        data = json.loads(out.read_text(encoding='utf-8'))
        assert data['seed'] == 7
        assert data['backend'] == 'toy'
        assert data['config_path'] == 'configs/boost/paper_faithful.yaml'
        assert data['git_sha']
        assert data['timestamp']

    def test_extra_contains_effective_params(self, tmp_path):
        extra = {
            'timestep_fs': 0.25,
            'biased_steps': 2000,
            'unbiased_steps': 2000,
            'n_cycles': 10,
            'seed': 7,
            'tdbb': {
                'f2': 10.0,
                'gamma': 1.0,
                'f1_max_formation': 250.0,
                'f1_max_dissociation': 125.0,
                'lambda_vdw': 0.60,
            },
            'backend': 'toy',
            'candidate_r_min': 3.0,
            'candidate_r_max': 6.0,
        }
        manifest = RunManifest(
            config_path='configs/boost/paper_faithful.yaml',
            seed=7,
            backend='toy',
            output_dir=str(tmp_path),
            extra=extra,
        )
        out = tmp_path / 'manifest.json'
        manifest.save(out)
        data = json.loads(out.read_text(encoding='utf-8'))

        assert 'extra' in data
        ex = data['extra']
        assert ex['tdbb']['f2'] == 10.0
        assert ex['tdbb']['gamma'] == 1.0
        assert ex['tdbb']['f1_max_formation'] == 250.0
        assert ex['tdbb']['f1_max_dissociation'] == 125.0
        assert ex['tdbb']['lambda_vdw'] == 0.60
        assert ex['timestep_fs'] == 0.25
        assert ex['biased_steps'] == 2000
        assert ex['unbiased_steps'] == 2000
        assert ex['n_cycles'] == 10
        assert ex['candidate_r_min'] == 3.0
        assert ex['candidate_r_max'] == 6.0
        assert ex['backend'] == 'toy'


class TestGitProvenance:
    """RF17: dirty-working-tree detection so a recorded SHA is trustworthy."""

    def test_git_dirty_field_is_serialized(self, tmp_path):
        manifest = RunManifest(
            config_path='c', seed=1, backend='toy', output_dir=str(tmp_path),
        )
        out = tmp_path / 'manifest.json'
        manifest.save(out)
        data = json.loads(out.read_text(encoding='utf-8'))
        assert 'git_dirty' in data
        assert isinstance(data['git_dirty'], bool)

    def test_get_git_dirty_true_when_porcelain_nonempty(self, monkeypatch):
        import kagome.workflows.manifest as mod

        class _R:
            returncode = 0
            stdout = ' M src/foo.py\n?? new.py\n'

        monkeypatch.setattr(mod.subprocess, 'run', lambda *a, **k: _R())
        assert mod._get_git_dirty() is True

    def test_get_git_dirty_false_when_clean(self, monkeypatch):
        import kagome.workflows.manifest as mod

        class _R:
            returncode = 0
            stdout = '   \n'

        monkeypatch.setattr(mod.subprocess, 'run', lambda *a, **k: _R())
        assert mod._get_git_dirty() is False

    def test_get_git_dirty_false_on_error(self, monkeypatch):
        import kagome.workflows.manifest as mod

        def _boom(*a, **k):
            raise OSError('no git')

        monkeypatch.setattr(mod.subprocess, 'run', _boom)
        assert mod._get_git_dirty() is False
