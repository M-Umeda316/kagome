"""Unit tests for quality gate scripts.

Tests that:
- check_dependency_licenses.py fails on blocked imports, passes otherwise.
- validate_configs.py passes valid configs, fails on missing/bad keys.
- check_seed_defined.py requires --seed.
- check_output_path.py requires --output-dir.
"""
from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from scripts.check_dependency_licenses import check as check_licenses, _find_imports_in_source
from scripts.check_seed_defined import check as check_seed
from scripts.check_output_path import check as check_output
from scripts.validate_configs import validate_file


# ---------- check_dependency_licenses ----------

def _write_yaml(tmp_path: Path, content: str) -> Path:
    p = tmp_path / 'approved.yaml'
    p.write_text(textwrap.dedent(content), encoding='utf-8')
    return p


def test_license_check_passes_when_blocked_not_imported(tmp_path):
    yaml_file = _write_yaml(tmp_path, """
        software:
          - name: pfp
            status: blocked_pending_review
            reason: proprietary
        models: []
    """)
    assert check_licenses(yaml_file, [tmp_path]) == 0


def test_license_check_fails_when_blocked_is_imported(tmp_path):
    yaml_file = _write_yaml(tmp_path, """
        software:
          - name: pfp
            status: blocked_pending_review
            reason: proprietary
        models: []
    """)
    # Create a fake source file that imports pfp
    src_dir = tmp_path / 'src'
    src_dir.mkdir()
    (src_dir / 'bad.py').write_text('import pfp\n', encoding='utf-8')
    assert check_licenses(yaml_file, [src_dir]) == 1


def test_license_check_warns_review_required_not_imported(tmp_path, capsys):
    yaml_file = _write_yaml(tmp_path, """
        software:
          - name: openmm
            status: review_required
            required_evidence: upstream_license_url
        models: []
    """)
    result = check_licenses(yaml_file, [tmp_path])
    assert result == 0
    captured = capsys.readouterr()
    assert 'WARNING' in captured.out


def test_license_check_all_approved(tmp_path):
    yaml_file = _write_yaml(tmp_path, """
        software:
          - name: numpy
            status: approved
            evidence: BSD-3-Clause
        models:
          - name: toy
            status: approved
    """)
    assert check_licenses(yaml_file, [tmp_path]) == 0


# ---------- check_seed_defined ----------

def test_seed_check_passes():
    assert check_seed(['--output-dir', 'runs/test', '--seed', '7']) == 0


def test_seed_check_fails_missing():
    assert check_seed(['--output-dir', 'runs/test']) == 1


def test_seed_check_fails_no_value():
    assert check_seed(['--seed']) == 1


def test_seed_check_fails_negative():
    assert check_seed(['--seed', '-1']) == 1


def test_seed_check_fails_non_integer():
    assert check_seed(['--seed', 'abc']) == 1


# ---------- check_output_path ----------

def test_output_check_passes():
    assert check_output(['--seed', '7', '--output-dir', 'runs/test']) == 0


def test_output_check_fails_missing():
    assert check_output(['--seed', '7']) == 1


def test_output_check_fails_no_value():
    assert check_output(['--output-dir']) == 1


# ---------- validate_configs ----------

def test_validate_boost_config_passes(tmp_path):
    boost_dir = tmp_path / 'boost'
    boost_dir.mkdir()
    config = boost_dir / 'test.yaml'
    config.write_text(
        'boost:\n'
        '  timestep_fs: 0.25\n'
        '  biased_steps: 2000\n'
        '  unbiased_steps: 2000\n'
        '  lambda_vdw: 0.60\n'
        '  f2_invA2: 10.0\n'
        '  gamma: 1.0\n'
        '  f1_max_form_kcal_mol: 250.0\n'
        '  f1_max_break_kcal_mol: 125.0\n',
        encoding='utf-8',
    )
    errors = validate_file(config)
    assert errors == [], errors


def test_validate_boost_config_fails_missing_key(tmp_path):
    boost_dir = tmp_path / 'boost'
    boost_dir.mkdir()
    config = boost_dir / 'test.yaml'
    config.write_text(
        'boost:\n'
        '  biased_steps: 2000\n',  # missing most required keys
        encoding='utf-8',
    )
    errors = validate_file(config)
    assert len(errors) > 0


def test_validate_boost_config_fails_invalid_value(tmp_path):
    boost_dir = tmp_path / 'boost'
    boost_dir.mkdir()
    config = boost_dir / 'test.yaml'
    config.write_text(
        'boost:\n'
        '  timestep_fs: -1.0\n'  # invalid: must be > 0
        '  biased_steps: 2000\n'
        '  unbiased_steps: 2000\n'
        '  lambda_vdw: 0.60\n'
        '  f2_invA2: 10.0\n'
        '  gamma: 1.0\n'
        '  f1_max_form_kcal_mol: 250.0\n'
        '  f1_max_break_kcal_mol: 125.0\n',
        encoding='utf-8',
    )
    errors = validate_file(config)
    assert any('timestep_fs' in e for e in errors)
