"""Unit tests for scripts/_mixing_cli.py — the shared ``--mix`` CLI machinery.

Extracted from run_vinyl_copolymer and wired into run_nylon66 /
run_epoxy_amine (specs/decisions.md 2026-08-04). No MD, no MLIP: the driver
tests replace the backend with the toy calculator and the workflow with a
recorder, so ``main()`` stops at wf.run.

Covers:
* ``mixing_setup_from_args`` — the checkpoint/resume record (all seven fields);
* ``mixing_setup_mismatch`` — identical / differing knob / exactly one None /
  both None, and the OLD-checkpoint case whose stored dict predates
  friction_per_ps + temperature_K (absent keys cannot be verified);
* ``MIX_KNOB_DEFAULTS`` and its resolution by ``resolve_mixing_args``;
* the symmetric stray-knob guard (a knob without --mix is a usage error);
* ``mix_config_from_args`` -> MixConfig field mapping;
* ``collect_mixing_skips`` / ``mixing_summary_fields`` — the summary.json
  provenance block, including the WM-P4 de-clash red-flag counter;
* end-to-end wiring in run_nylon66 / run_epoxy_amine: PolymerizationConfig.mixing,
  checkpoint_extra['mixing'] and the summary.json mixing block, with and
  without --mix;
* the nylon/epoxy guard that --mix with a failed bond-topology extraction is a
  hard error (mixing needs the topology; decisions.md 2026-08-04 nylon 固有ガード).

The tests that need ``--mix`` to survive the availability check are skipped when
OpenMM/OpenFF are absent (same convention as test_workflow_mixing.py).
"""
from __future__ import annotations

import argparse
import json
import logging
import sys

import pytest

from kagome.backends.toy import ToyCalculator
from kagome.workflows.polymerization import MixConfig
from scripts._mixing_cli import (
    MIX_KNOB_DEFAULTS,
    add_mixing_arguments,
    collect_mixing_skips,
    mix_config_from_args,
    mixing_setup_from_args,
    mixing_setup_mismatch,
    mixing_summary_fields,
    resolve_mixing_args,
)

_LOG = logging.getLogger('test_mixing_cli')

# The summary.json mixing block, in order. Pinned because run_vinyl_copolymer
# splices it with ** where it used to write the keys inline: a reorder or rename
# would silently change that script's long-standing summary.json layout.
MIXING_SUMMARY_KEYS = (
    'mixing_enabled', 'mix_ps', 'mix_settle_steps', 'mix_timestep_fs',
    'mix_platform', 'mix_charge_method', 'mixing_skipped_cycles',
    'n_mixing_skipped',
)


def _has_openmm() -> bool:
    try:
        import openmm  # noqa: F401
        import openff.toolkit  # noqa: F401
        return True
    except Exception:
        return False


requires_openmm = pytest.mark.skipif(
    not _has_openmm(), reason='OpenMM/OpenFF not installed',
)


def _args(**overrides) -> argparse.Namespace:
    """A resolved (post-sentinel) args namespace with --mix ON by default."""
    base = dict(
        mix=True,
        mix_ps=1.0,
        mix_settle_steps=100,
        mix_timestep_fs=0.5,
        mix_platform='CPU',
        mix_charge_method='nagl',
        mix_friction_per_ps=1.0,
        temperature=333.0,
    )
    base.update(overrides)
    return argparse.Namespace(**base)


def _parser() -> argparse.ArgumentParser:
    """A bare parser carrying only the mixing block plus --temperature."""
    p = argparse.ArgumentParser()
    p.add_argument('--temperature', type=float, default=333.0)
    add_mixing_arguments(p)
    return p


# ── mixing_setup_from_args: the checkpoint record ─────────────────────────────

class TestMixingSetupFromArgs:

    def test_off_returns_none(self):
        assert mixing_setup_from_args(_args(mix=False)) is None

    def test_includes_friction_and_temperature(self):
        setup = mixing_setup_from_args(_args())
        # the two fields the review found missing from the guarded setup
        assert setup['friction_per_ps'] == 1.0
        assert setup['temperature_K'] == 333.0
        # and the pre-existing five
        assert setup['mix_ps'] == 1.0
        assert setup['mix_settle_steps'] == 100
        assert setup['mix_timestep_fs'] == 0.5
        assert setup['mix_platform'] == 'CPU'
        assert setup['mix_charge_method'] == 'nagl'


# ── mixing_setup_mismatch: the resume mode-switch guard ───────────────────────

class TestResumeGuard:

    def test_matched_config_resumes_cleanly(self):
        ckpt = mixing_setup_from_args(_args())
        now = mixing_setup_from_args(_args())
        assert mixing_setup_mismatch(ckpt, now) is False

    def test_changed_mix_ps_is_mismatch(self):
        ckpt = mixing_setup_from_args(_args(mix_ps=1.0))
        now = mixing_setup_from_args(_args(mix_ps=2.0))
        assert mixing_setup_mismatch(ckpt, now) is True

    def test_changed_friction_is_mismatch(self):
        ckpt = mixing_setup_from_args(_args(mix_friction_per_ps=1.0))
        now = mixing_setup_from_args(_args(mix_friction_per_ps=5.0))
        assert mixing_setup_mismatch(ckpt, now) is True

    def test_changed_temperature_is_mismatch(self):
        ckpt = mixing_setup_from_args(_args(temperature=333.0))
        now = mixing_setup_from_args(_args(temperature=300.0))
        assert mixing_setup_mismatch(ckpt, now) is True

    def test_turning_mixing_on_or_off_is_mismatch(self):
        on = mixing_setup_from_args(_args())
        assert mixing_setup_mismatch(None, on) is True     # was off, now on
        assert mixing_setup_mismatch(on, None) is True     # was on, now off
        assert mixing_setup_mismatch(None, None) is False  # off both runs

    def test_old_checkpoint_missing_new_keys_does_not_spuriously_mismatch(self):
        """An OLD checkpoint predates friction_per_ps/temperature_K; those keys
        are absent from its stored dict and cannot be verified, so they must be
        SKIPPED (not compared against a default). Shared keys still match -> no
        spurious error."""
        old = mixing_setup_from_args(_args())
        del old['friction_per_ps']
        del old['temperature_K']
        now = mixing_setup_from_args(_args())
        assert mixing_setup_mismatch(old, now) is False

    def test_absent_checkpoint_key_cannot_flag_a_change(self):
        """Even if the current friction differs, an OLD checkpoint that never
        stored friction_per_ps cannot flag it (the field predates the
        checkpoint) — only shared keys are compared."""
        old = mixing_setup_from_args(_args(mix_friction_per_ps=1.0))
        del old['friction_per_ps']
        del old['temperature_K']
        now = mixing_setup_from_args(_args(mix_friction_per_ps=99.0))
        assert mixing_setup_mismatch(old, now) is False
        # but a change in a shared key IS still caught
        now_bad = mixing_setup_from_args(
            _args(mix_ps=2.0, mix_friction_per_ps=99.0))
        assert mixing_setup_mismatch(old, now_bad) is True


# ── knob defaults (WM-P5b: 25 ps default) ─────────────────────────────────────

class TestKnobDefaults:

    def test_documented_default_values(self):
        """The WM-P5b sweep result (decisions.md 追補 2026-07-22)."""
        assert MIX_KNOB_DEFAULTS['mix_ps'] == 25.0
        assert MIX_KNOB_DEFAULTS['mix_settle_steps'] == 500
        assert MIX_KNOB_DEFAULTS['mix_timestep_fs'] == 0.5
        assert MIX_KNOB_DEFAULTS['mix_friction_per_ps'] == 1.0
        assert MIX_KNOB_DEFAULTS['mix_platform'] == 'CPU'
        assert MIX_KNOB_DEFAULTS['mix_charge_method'] == 'nagl'

    def test_knobs_parse_to_none_sentinels(self):
        """"not given" must stay distinguishable from "given" — the stray-knob
        guard depends on it."""
        args = _parser().parse_args([])
        assert args.mix is False
        for knob in MIX_KNOB_DEFAULTS:
            assert getattr(args, knob) is None

    @requires_openmm
    def test_unset_knobs_resolve_to_defaults(self):
        p = _parser()
        args = p.parse_args(['--mix'])
        resolve_mixing_args(p, args)
        for knob, default in MIX_KNOB_DEFAULTS.items():
            assert getattr(args, knob) == default

    @requires_openmm
    def test_explicit_values_are_preserved(self):
        p = _parser()
        args = p.parse_args(['--mix', '--mix-ps', '100', '--mix-settle-steps', '750'])
        resolve_mixing_args(p, args)
        assert args.mix_ps == 100.0        # not overwritten by the default
        assert args.mix_settle_steps == 750
        assert args.mix_timestep_fs == MIX_KNOB_DEFAULTS['mix_timestep_fs']


# ── symmetric stray-knob guard ────────────────────────────────────────────────

@pytest.mark.parametrize('flag', [
    ['--mix-ps', '1.0'],
    ['--mix-settle-steps', '100'],
    ['--mix-friction-per-ps', '1.0'],
    ['--mix-timestep-fs', '0.5'],
    ['--mix-platform', 'CPU'],
    ['--mix-charge-method', 'nagl'],
])
def test_stray_knob_without_mix_errors(flag):
    p = _parser()
    args = p.parse_args(flag)
    with pytest.raises(SystemExit):
        resolve_mixing_args(p, args)


def test_no_mixing_flags_is_clean():
    p = _parser()
    args = p.parse_args([])
    resolve_mixing_args(p, args)         # must not raise
    assert mixing_setup_from_args(args) is None


# ── mix_config_from_args ──────────────────────────────────────────────────────

class TestMixConfigFromArgs:

    def test_off_returns_none(self):
        assert mix_config_from_args(_args(mix=False), 333.0) is None

    def test_field_mapping(self):
        cfg = mix_config_from_args(_args(), 300.0)
        assert isinstance(cfg, MixConfig)
        assert cfg.mix_time_ps == 1.0
        assert cfg.settle_steps == 100
        assert cfg.temperature_K == 300.0   # the production setpoint, not a knob
        assert cfg.timestep_fs == 0.5
        assert cfg.friction_per_ps == 1.0
        assert cfg.platform == 'CPU'
        assert cfg.charge_method == 'nagl'


# ── summary.json provenance (decisions.md 2026-08-04 追補) ────────────────────

def _write_mixing_log(tmp_path, records) -> None:
    (tmp_path / 'mixing.jsonl').write_text(
        ''.join(json.dumps(r) + '\n' for r in records), encoding='utf-8')


class TestCollectMixingSkips:

    def test_no_log_file_is_empty(self, tmp_path):
        assert collect_mixing_skips(_args(n_cycles=3), tmp_path, _LOG) == []

    def test_mixing_off_ignores_an_existing_log(self, tmp_path):
        _write_mixing_log(tmp_path, [{'cycle': 0, 'skipped': True}])
        args = _args(mix=False, n_cycles=3)
        assert collect_mixing_skips(args, tmp_path, _LOG) == []

    def test_collects_only_skipped_cycles(self, tmp_path):
        _write_mixing_log(tmp_path, [
            {'cycle': 0, 'skipped': False},
            {'cycle': 1, 'skipped': True, 'skip_reason': 'diverged'},
            {'cycle': 2},                                   # no key => not skipped
            {'cycle': 3, 'skipped': True},
        ])
        args = _args(n_cycles=4)
        assert collect_mixing_skips(args, tmp_path, _LOG) == [1, 3]

    def test_blank_and_malformed_lines_are_tolerated(self, tmp_path):
        (tmp_path / 'mixing.jsonl').write_text(
            '\n{"cycle": 0, "skipped": true}\nnot json\n\n', encoding='utf-8')
        args = _args(n_cycles=2)
        assert collect_mixing_skips(args, tmp_path, _LOG) == [0]

    def test_skips_are_warned(self, tmp_path, caplog):
        _write_mixing_log(tmp_path, [{'cycle': 1, 'skipped': True}])
        with caplog.at_level(logging.WARNING, logger='test_mixing_cli'):
            collect_mixing_skips(_args(n_cycles=2), tmp_path, _LOG)
        assert 'SKIPPED' in caplog.text


class TestMixingSummaryFields:

    def test_key_order_is_pinned(self):
        assert tuple(mixing_summary_fields(_args(), [])) == MIXING_SUMMARY_KEYS

    def test_mixing_on(self):
        fields = mixing_summary_fields(_args(), [2, 5])
        assert fields['mixing_enabled'] is True
        assert fields['mix_ps'] == 1.0
        assert fields['mix_settle_steps'] == 100
        assert fields['mix_timestep_fs'] == 0.5
        assert fields['mix_platform'] == 'CPU'
        assert fields['mix_charge_method'] == 'nagl'
        assert fields['mixing_skipped_cycles'] == [2, 5]
        assert fields['n_mixing_skipped'] == 2

    def test_mixing_off_is_all_none(self):
        """--mix off leaves every knob at its None sentinel, so the block stays
        explicit (recorded as off) rather than absent."""
        off = argparse.Namespace(
            mix=False, temperature=300.0,
            **{knob: None for knob in MIX_KNOB_DEFAULTS})
        fields = mixing_summary_fields(off, [])
        assert fields['mixing_enabled'] is False
        assert all(fields[k] is None for k in MIXING_SUMMARY_KEYS[1:6])
        assert fields['mixing_skipped_cycles'] == []
        assert fields['n_mixing_skipped'] == 0


# ── driver wiring: run_nylon66 / run_epoxy_amine ──────────────────────────────

def _run_driver(monkeypatch, module, argv: list[str]) -> dict:
    """Run ``module.main()`` with the backend and the workflow stubbed out.

    Returns the captured PolymerizationWorkflow construction/run arguments; the
    real system build still runs (cheap at 1+1 molecules) so the wiring is
    exercised end to end up to wf.run.
    """
    captured: dict = {}

    class _RecordingWorkflow:
        def __init__(self, config, *args, **kwargs):
            captured['config'] = config
            captured['initial_bonds'] = kwargs.get('initial_bonds')
            self.config = config

        def run(self, state, **kwargs):
            captured['run_kwargs'] = kwargs
            return []

    monkeypatch.setattr(module, 'PolymerizationWorkflow', _RecordingWorkflow)
    monkeypatch.setattr(module, '_create_backend',
                        lambda *a, **kw: ToyCalculator())
    monkeypatch.setattr(sys, 'argv', argv)
    module.main()
    return captured


def _summary(tmp_path) -> dict:
    return json.loads((tmp_path / 'summary.json').read_text(encoding='utf-8'))


def _nylon_argv(tmp_path, *extra: str) -> list[str]:
    return [
        'run_nylon66.py', '--output-dir', str(tmp_path),
        '--n-diamines', '1', '--n-diacids', '1', '--n-cycles', '1',
        '--box-size', '25', *extra,
    ]


def _epoxy_argv(tmp_path, *extra: str) -> list[str]:
    return [
        'run_epoxy_amine.py', '--output-dir', str(tmp_path),
        '--n-epoxies', '1', '--n-amines', '1', '--n-cycles', '1',
        '--box-size', '30', *extra,
    ]


class TestNylonWiring:

    def test_without_mix_is_off(self, monkeypatch, tmp_path):
        from scripts import run_nylon66
        cap = _run_driver(monkeypatch, run_nylon66, _nylon_argv(tmp_path))
        assert cap['config'].mixing is None
        assert cap['run_kwargs']['checkpoint_extra'] == {'mixing': None}
        summary = _summary(tmp_path)
        assert set(MIXING_SUMMARY_KEYS) <= set(summary)
        assert summary['mixing_enabled'] is False
        assert summary['mix_ps'] is None
        assert summary['mixing_skipped_cycles'] == []
        assert summary['n_mixing_skipped'] == 0

    @requires_openmm
    def test_mix_builds_mixconfig_and_checkpoint_extra(self, monkeypatch, tmp_path):
        from scripts import run_nylon66
        cap = _run_driver(monkeypatch, run_nylon66, _nylon_argv(
            tmp_path, '--mix', '--mix-ps', '2.0', '--mix-settle-steps', '10',
            '--mix-platform', 'Reference', '--mix-charge-method', 'gasteiger',
            '--temperature', '300'))
        mixing = cap['config'].mixing
        assert isinstance(mixing, MixConfig)
        assert mixing.mix_time_ps == 2.0
        assert mixing.settle_steps == 10
        assert mixing.platform == 'Reference'
        assert mixing.charge_method == 'gasteiger'
        assert mixing.temperature_K == 300.0     # follows --temperature
        assert mixing.timestep_fs == MIX_KNOB_DEFAULTS['mix_timestep_fs']
        assert cap['run_kwargs']['checkpoint_extra']['mixing'] == {
            'mix_ps': 2.0,
            'mix_settle_steps': 10,
            'mix_timestep_fs': 0.5,
            'mix_platform': 'Reference',
            'mix_charge_method': 'gasteiger',
            'friction_per_ps': 1.0,
            'temperature_K': 300.0,
        }
        # mixing needs the bond topology, which must have been extracted
        assert cap['initial_bonds'] is not None
        summary = _summary(tmp_path)
        assert summary['mixing_enabled'] is True
        assert summary['mix_ps'] == 2.0
        assert summary['mix_settle_steps'] == 10
        assert summary['mix_timestep_fs'] == 0.5
        assert summary['mix_platform'] == 'Reference'
        assert summary['mix_charge_method'] == 'gasteiger'
        assert summary['mixing_skipped_cycles'] == []
        assert summary['n_mixing_skipped'] == 0

    @requires_openmm
    def test_mix_without_bond_topology_is_a_hard_error(self, monkeypatch, tmp_path):
        """init_bonds extraction is best-effort for a plain run, but mixing
        cannot work without it — fail at the CLI, not inside wf.run."""
        from scripts import run_nylon66

        def _boom(*_a, **_kw):
            raise RuntimeError('synthetic topology failure')

        monkeypatch.setattr(run_nylon66, 'layout_bonds', _boom)
        with pytest.raises(SystemExit):
            _run_driver(monkeypatch, run_nylon66,
                        _nylon_argv(tmp_path, '--mix', '--mix-ps', '2.0'))

    def test_topology_failure_without_mix_still_runs(self, monkeypatch, tmp_path):
        """The best-effort path is unchanged when mixing is off."""
        from scripts import run_nylon66

        def _boom(*_a, **_kw):
            raise RuntimeError('synthetic topology failure')

        monkeypatch.setattr(run_nylon66, 'layout_bonds', _boom)
        cap = _run_driver(monkeypatch, run_nylon66, _nylon_argv(tmp_path))
        assert cap['initial_bonds'] is None
        assert cap['config'].mixing is None

    def test_stray_knob_without_mix_errors(self, monkeypatch, tmp_path):
        from scripts import run_nylon66
        with pytest.raises(SystemExit):
            _run_driver(monkeypatch, run_nylon66,
                        _nylon_argv(tmp_path, '--mix-ps', '2.0'))


class TestEpoxyWiring:

    def test_without_mix_is_off(self, monkeypatch, tmp_path):
        from scripts import run_epoxy_amine
        cap = _run_driver(monkeypatch, run_epoxy_amine, _epoxy_argv(tmp_path))
        assert cap['config'].mixing is None
        assert cap['run_kwargs']['checkpoint_extra'] == {'mixing': None}
        summary = _summary(tmp_path)
        assert set(MIXING_SUMMARY_KEYS) <= set(summary)
        assert summary['mixing_enabled'] is False
        assert summary['mix_ps'] is None
        assert summary['n_mixing_skipped'] == 0

    @requires_openmm
    def test_mix_builds_mixconfig_and_checkpoint_extra(self, monkeypatch, tmp_path):
        from scripts import run_epoxy_amine
        cap = _run_driver(monkeypatch, run_epoxy_amine, _epoxy_argv(
            tmp_path, '--mix', '--mix-ps', '3.0', '--mix-settle-steps', '20',
            '--mix-platform', 'Reference', '--mix-charge-method', 'gasteiger'))
        mixing = cap['config'].mixing
        assert isinstance(mixing, MixConfig)
        assert mixing.mix_time_ps == 3.0
        assert mixing.settle_steps == 20
        assert mixing.temperature_K == 333.0     # epoxy production setpoint
        assert cap['run_kwargs']['checkpoint_extra']['mixing']['mix_ps'] == 3.0
        assert cap['initial_bonds'] is not None
        summary = _summary(tmp_path)
        assert summary['mixing_enabled'] is True
        assert summary['mix_ps'] == 3.0
        assert summary['mix_settle_steps'] == 20
        assert summary['mix_platform'] == 'Reference'
        assert summary['mix_charge_method'] == 'gasteiger'
        assert summary['mixing_skipped_cycles'] == []
        assert summary['n_mixing_skipped'] == 0

    @requires_openmm
    def test_mix_without_bond_topology_is_a_hard_error(self, monkeypatch, tmp_path):
        from scripts import run_epoxy_amine

        def _boom(*_a, **_kw):
            raise RuntimeError('synthetic topology failure')

        monkeypatch.setattr(run_epoxy_amine, 'layout_bonds', _boom)
        with pytest.raises(SystemExit):
            _run_driver(monkeypatch, run_epoxy_amine,
                        _epoxy_argv(tmp_path, '--mix', '--mix-ps', '3.0'))

    def test_stray_knob_without_mix_errors(self, monkeypatch, tmp_path):
        from scripts import run_epoxy_amine
        with pytest.raises(SystemExit):
            _run_driver(monkeypatch, run_epoxy_amine,
                        _epoxy_argv(tmp_path, '--mix-ps', '3.0'))
