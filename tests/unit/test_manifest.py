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


class TestAppendResume:
    """W1: resume must preserve the original run's provenance, not overwrite it."""

    def _fresh(self, tmp_path):
        manifest = RunManifest(
            config_path='c', seed=7, backend='toy', output_dir=str(tmp_path),
            extra={'equil_steps': 500},
        )
        out = tmp_path / 'manifest.json'
        manifest.save(out)
        return out, manifest

    def test_original_git_sha_preserved(self, tmp_path):
        out, manifest = self._fresh(tmp_path)
        original_sha = manifest.git_sha
        original_ts = manifest.timestamp

        ok = RunManifest.append_resume(out, ckpt_step=1234, ckpt_cycle=5)
        assert ok is True
        data = json.loads(out.read_text(encoding='utf-8'))
        # Top-level provenance of the original run is untouched.
        assert data['git_sha'] == original_sha
        assert data['timestamp'] == original_ts
        assert data['seed'] == 7

    def test_resume_history_grows_by_one(self, tmp_path):
        out, _ = self._fresh(tmp_path)
        RunManifest.append_resume(out, ckpt_step=100, ckpt_cycle=1)
        RunManifest.append_resume(out, ckpt_step=200, ckpt_cycle=2)
        data = json.loads(out.read_text(encoding='utf-8'))
        history = data['extra']['resume_history']
        assert len(history) == 2
        assert history[0]['ckpt_step'] == 100
        assert history[0]['ckpt_cycle'] == 1
        assert history[1]['ckpt_step'] == 200
        for entry in history:
            assert 'timestamp' in entry
            assert 'git_sha' in entry
            assert isinstance(entry['git_dirty'], bool)

    def test_original_extra_preserved(self, tmp_path):
        out, _ = self._fresh(tmp_path)
        RunManifest.append_resume(out, ckpt_step=10, ckpt_cycle=0)
        data = json.loads(out.read_text(encoding='utf-8'))
        assert data['extra']['equil_steps'] == 500

    def test_result_is_valid_json(self, tmp_path):
        out, _ = self._fresh(tmp_path)
        RunManifest.append_resume(out, ckpt_step=10, ckpt_cycle=0)
        # Must parse cleanly (no truncation/corruption).
        json.loads(out.read_text(encoding='utf-8'))

    def test_missing_file_returns_false(self, tmp_path):
        missing = tmp_path / 'nope.json'
        assert RunManifest.append_resume(missing, ckpt_step=1, ckpt_cycle=1) is False


class TestProductionStartStep:
    """A4: record the post-equilibration production onset for k_p fitting."""

    def test_records_step_into_extra(self, tmp_path):
        manifest = RunManifest(
            config_path='c', seed=1, backend='toy', output_dir=str(tmp_path),
        )
        out = tmp_path / 'manifest.json'
        manifest.save(out)
        RunManifest.record_production_start_step(out, 2048)
        data = json.loads(out.read_text(encoding='utf-8'))
        assert data['extra']['production_start_step'] == 2048

    def test_no_op_when_file_missing(self, tmp_path):
        missing = tmp_path / 'nope.json'
        # Should not raise.
        RunManifest.record_production_start_step(missing, 100)
        assert not missing.exists()


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
