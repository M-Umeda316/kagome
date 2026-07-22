"""Unit tests for scripts/run_vinyl_copolymer.py CLI guards.

Pure-Python (no orb/GPU/openmm): all covered logic fires in argparse / the
mixing-setup guard helpers before any system build or backend creation.

Covers (review-driven mixing correctness batch):
* the resume mode-switch guard (``_mixing_setup_from_args`` /
  ``_mixing_setup_mismatch``): matched configs resume cleanly; a changed knob
  (mix_ps / friction_per_ps / temperature) is a hard mismatch; an OLD checkpoint
  whose stored mixing dict predates friction_per_ps/temperature_K does NOT
  spuriously mismatch (a field absent from the checkpoint cannot be verified);
* the symmetric CLI guard: any mixing knob given WITHOUT --mix is a usage error
  (parser.error -> SystemExit).
"""
from __future__ import annotations

import argparse
import sys

import pytest

from scripts.run_vinyl_copolymer import (
    MIX_KNOB_DEFAULTS,
    _apply_mix_defaults,
    _mixing_setup_from_args,
    _mixing_setup_mismatch,
    main,
)


def _args(**overrides):
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


# ── (a) resume mode-switch guard ──────────────────────────────────────────────

class TestMixingSetupFromArgs:

    def test_off_returns_none(self):
        assert _mixing_setup_from_args(_args(mix=False)) is None

    def test_includes_friction_and_temperature(self):
        setup = _mixing_setup_from_args(_args())
        # the two fields the review found missing from the guarded setup
        assert setup['friction_per_ps'] == 1.0
        assert setup['temperature_K'] == 333.0
        # and the pre-existing five
        assert setup['mix_ps'] == 1.0
        assert setup['mix_settle_steps'] == 100
        assert setup['mix_timestep_fs'] == 0.5
        assert setup['mix_platform'] == 'CPU'
        assert setup['mix_charge_method'] == 'nagl'


class TestResumeGuard:

    def test_matched_config_resumes_cleanly(self):
        ckpt = _mixing_setup_from_args(_args())
        now = _mixing_setup_from_args(_args())
        assert _mixing_setup_mismatch(ckpt, now) is False

    def test_changed_mix_ps_is_mismatch(self):
        ckpt = _mixing_setup_from_args(_args(mix_ps=1.0))
        now = _mixing_setup_from_args(_args(mix_ps=2.0))
        assert _mixing_setup_mismatch(ckpt, now) is True

    def test_changed_friction_is_mismatch(self):
        ckpt = _mixing_setup_from_args(_args(mix_friction_per_ps=1.0))
        now = _mixing_setup_from_args(_args(mix_friction_per_ps=5.0))
        assert _mixing_setup_mismatch(ckpt, now) is True

    def test_changed_temperature_is_mismatch(self):
        ckpt = _mixing_setup_from_args(_args(temperature=333.0))
        now = _mixing_setup_from_args(_args(temperature=300.0))
        assert _mixing_setup_mismatch(ckpt, now) is True

    def test_turning_mixing_on_or_off_is_mismatch(self):
        on = _mixing_setup_from_args(_args())
        assert _mixing_setup_mismatch(None, on) is True     # was off, now on
        assert _mixing_setup_mismatch(on, None) is True     # was on, now off
        assert _mixing_setup_mismatch(None, None) is False  # off both runs

    def test_old_checkpoint_missing_new_keys_does_not_spuriously_mismatch(self):
        """An OLD checkpoint predates friction_per_ps/temperature_K; those keys
        are absent from its stored dict and cannot be verified, so they must be
        SKIPPED (not compared against a default). Shared keys still match -> no
        spurious error."""
        old = _mixing_setup_from_args(_args())
        del old['friction_per_ps']
        del old['temperature_K']
        now = _mixing_setup_from_args(_args())
        assert _mixing_setup_mismatch(old, now) is False

    def test_absent_checkpoint_key_cannot_flag_a_change(self):
        """Even if the current friction differs, an OLD checkpoint that never
        stored friction_per_ps cannot flag it (the field predates the
        checkpoint) — only shared keys are compared."""
        old = _mixing_setup_from_args(_args(mix_friction_per_ps=1.0))
        del old['friction_per_ps']
        del old['temperature_K']
        now = _mixing_setup_from_args(_args(mix_friction_per_ps=99.0))
        assert _mixing_setup_mismatch(old, now) is False
        # but a change in a shared key IS still caught
        now_bad = _mixing_setup_from_args(
            _args(mix_ps=2.0, mix_friction_per_ps=99.0))
        assert _mixing_setup_mismatch(old, now_bad) is True


# ── (c) mixing-knob default resolution (WM-P5b: 25 ps default) ────────────────

class TestApplyMixDefaults:
    def test_mix_ps_default_is_25(self):
        """The WM-P5b sweep result (decisions.md 追補 2026-07-22)."""
        assert MIX_KNOB_DEFAULTS['mix_ps'] == 25.0
        assert MIX_KNOB_DEFAULTS['mix_settle_steps'] == 500

    def test_unset_knobs_resolve_to_defaults(self):
        ns = argparse.Namespace(
            mix=True, mix_ps=None, mix_settle_steps=None, mix_timestep_fs=None,
            mix_friction_per_ps=None, mix_platform=None, mix_charge_method=None)
        _apply_mix_defaults(ns)
        assert ns.mix_ps == 25.0
        assert ns.mix_settle_steps == 500
        assert ns.mix_timestep_fs == 0.5
        assert ns.mix_friction_per_ps == 1.0
        assert ns.mix_platform == 'CPU'
        assert ns.mix_charge_method == 'nagl'

    def test_explicit_values_are_preserved(self):
        ns = argparse.Namespace(
            mix=True, mix_ps=100.0, mix_settle_steps=750, mix_timestep_fs=None,
            mix_friction_per_ps=None, mix_platform=None, mix_charge_method=None)
        _apply_mix_defaults(ns)
        assert ns.mix_ps == 100.0        # not overwritten by the default
        assert ns.mix_settle_steps == 750


# ── (b) symmetric CLI guard: mixing knob without --mix errors ─────────────────

@pytest.mark.parametrize('flag', [
    ['--mix-ps', '1.0'],                 # pre-existing behavior, pinned
    ['--mix-settle-steps', '100'],       # pre-existing behavior, pinned
    ['--mix-friction-per-ps', '1.0'],    # newly guarded
    ['--mix-timestep-fs', '0.5'],        # newly guarded
    ['--mix-platform', 'CPU'],           # newly guarded
    ['--mix-charge-method', 'nagl'],     # newly guarded
])
def test_mix_knob_without_mix_errors(monkeypatch, tmp_path, flag):
    """Any mixing knob supplied WITHOUT --mix is a usage error. The guard fires
    right after argparse, before any system build or backend creation."""
    monkeypatch.setattr(sys, 'argv', [
        'run_vinyl_copolymer.py',
        '--output-dir', str(tmp_path),
        '--backend', 'toy',
        *flag,
    ])
    with pytest.raises(SystemExit):
        main()
