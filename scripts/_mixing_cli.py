"""Shared CLI wiring for the optional WM-P3 classical mixing stage.

Extracted verbatim from ``scripts/run_vinyl_copolymer.py`` so the vinyl,
nylon-6,6 and epoxy-amine drivers expose ONE set of ``--mix*`` flags with
identical semantics (specs/decisions.md 2026-08-04 "混合(well-mixed)CLI 配線を
nylon-6,6 / epoxy-amine スクリプトへ拡張"). The mixing stage itself is unchanged:
it stays opt-in and OFF by default, so the paper-faithful path is untouched
(specs/decisions.md 2026-07-17 "well-mixed 測定モード" (i), 追補 2026-07-18).

Consumers wire it in five places:

1. ``add_mixing_arguments(parser)``       — the ``--mix`` flag + 6 knobs
2. ``resolve_mixing_args(parser, args)``  — default resolution + guards
3. ``mix_config_from_args(args, T)``      — the ``MixConfig`` for the workflow
4. ``mixing_setup_from_args(args)``       — checkpoint record + resume guard,
   compared via ``mixing_setup_mismatch``
5. ``collect_mixing_skips`` + ``mixing_summary_fields`` — summary.json
   provenance, including the de-clash red-flag counter
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from kagome.workflows.polymerization import MixConfig

__all__ = [
    'MIX_KNOB_DEFAULTS',
    'add_mixing_arguments',
    'collect_mixing_skips',
    'mix_config_from_args',
    'mixing_setup_from_args',
    'mixing_setup_mismatch',
    'mixing_summary_fields',
    'resolve_mixing_args',
]


# Documented defaults for the optional mixing knobs. Every knob carries an
# argparse sentinel of None so "not given" is distinguishable from "given"
# (the symmetric stray-knob guard depends on this); the sentinels are resolved
# to these values only when --mix is present. mix_ps=25.0 is the WM-P5b sweep
# result (decisions.md 追補 2026-07-22: refresh saturates by 25 ps, yield
# unresolved across 25/50/100, 25 ps cheapest); mix_settle_steps=500 is the
# established WM-P4/P5b campaign value. The rest mirror MixConfig's defaults.
MIX_KNOB_DEFAULTS = {
    'mix_ps': 25.0,
    'mix_settle_steps': 500,
    'mix_timestep_fs': 0.5,
    'mix_friction_per_ps': 1.0,
    'mix_platform': 'CPU',
    'mix_charge_method': 'nagl',
}


def mixing_setup_from_args(args: argparse.Namespace) -> dict | None:
    """The mixing setup persisted in the checkpoint and compared on resume.

    Returns ``None`` when ``--mix`` is off. Every field here is part of the
    measured mixing diffusion, so a change across resume must be rejected. Single
    source of truth so ``checkpoint_extra['mixing']`` and the resume-time
    comparison dict can never drift apart. ``args.mix_*`` sentinels must already
    be resolved to their real defaults (see ``resolve_mixing_args``).
    """
    if not args.mix:
        return None
    return {
        'mix_ps': args.mix_ps,
        'mix_settle_steps': args.mix_settle_steps,
        'mix_timestep_fs': args.mix_timestep_fs,
        'mix_platform': args.mix_platform,
        'mix_charge_method': args.mix_charge_method,
        # friction_per_ps and temperature_K both change the measured mixing
        # diffusion, so they belong to the guarded setup too.
        'friction_per_ps': args.mix_friction_per_ps,
        'temperature_K': args.temperature,
    }


def mixing_setup_mismatch(ckpt_mix: dict | None, now_mix: dict | None) -> bool:
    """True if the resumed run's mixing setup differs from the checkpoint's.

    Turning mixing on/off (exactly one side ``None``) is always a mismatch. When
    both are dicts, compare only keys present in BOTH: an OLD checkpoint predates
    ``friction_per_ps``/``temperature_K``, and a field that predates the
    checkpoint cannot be verified — skip it rather than flag a spurious mismatch.
    New checkpoints carry every key.
    """
    if (ckpt_mix is None) != (now_mix is None):
        return True
    if ckpt_mix is None:                     # both None: mixing off both runs
        return False
    shared = set(ckpt_mix) & set(now_mix)
    return any(ckpt_mix[k] != now_mix[k] for k in shared)


def add_mixing_arguments(parser: argparse.ArgumentParser) -> None:
    """Register ``--mix`` and its six knobs on ``parser``.

    WM-P3 (specs/decisions.md "2026-07-17: well-mixed 測定モード" item (i),
    追補 2026-07-18): optional per-cycle classical mixing stage. Off by
    default — the paper-faithful loop is unchanged. Mixing durations now have
    documented defaults resolved from MIX_KNOB_DEFAULTS (mix_ps=25 ps is the
    WM-P5b sweep result, decisions.md 追補 2026-07-22, lifting the earlier
    "no defaults until measured" rule of decision (v)).
    """
    parser.add_argument('--mix', action='store_true',
                        help='enable the classical (OpenFF/OpenMM) mixing '
                             'stage after each cycle (well-mixed measurement '
                             'mode; NOT paper-faithful).')
    parser.add_argument('--mix-ps', type=float, default=None,
                        help='classical mixing duration per cycle in ps '
                             '(default 25, WM-P5b); requires --mix.')
    parser.add_argument('--mix-settle-steps', type=int, default=None,
                        help='MLIP settle steps after coordinate write-back '
                             '(default 500, WM-P4/P5b); requires --mix.')
    # Like every mixing knob these carry a sentinel default of None so we can
    # tell "not given" from "given" and (a) error if any is passed WITHOUT --mix
    # (the symmetric stray-knob guard in resolve_mixing_args) and (b) resolve
    # None to the documented default (MIX_KNOB_DEFAULTS) only when --mix is
    # present.
    parser.add_argument('--mix-timestep-fs', type=float, default=None,
                        help='classical mixing timestep in fs (default 0.5, '
                             'ClassicalPrepConfig precedent); requires --mix.')
    parser.add_argument('--mix-friction-per-ps', type=float, default=None,
                        help='Langevin friction for the mixing MD in ps^-1 '
                             '(default 1.0); requires --mix.')
    parser.add_argument('--mix-platform',
                        choices=['CUDA', 'OpenCL', 'CPU', 'Reference'],
                        default=None,
                        help='OpenMM platform for the mixing MD (default CPU); '
                             'requires --mix.')
    parser.add_argument('--mix-charge-method', choices=['nagl', 'gasteiger'],
                        default=None,
                        help='partial-charge method for the mixing force field '
                             '(default nagl; NAGL falls back to Gasteiger); '
                             'requires --mix.')


def resolve_mixing_args(parser: argparse.ArgumentParser,
                        args: argparse.Namespace) -> None:
    """Resolve the mixing sentinels in place and enforce the CLI guards.

    With ``--mix``: every unset knob is resolved to its documented default now,
    so each downstream consumer (MixConfig, the resume guard, checkpoint_extra,
    summary.json) sees a real value rather than None; then OpenMM/OpenFF
    availability is checked here rather than deep inside the first mixing cycle.

    Without ``--mix``: symmetric guard — ANY mixing knob given without --mix is
    a usage error. Previously only --mix-ps/--mix-settle-steps errored; the four
    other knobs were silently ignored because they carried real defaults.
    """
    if args.mix:
        for knob, default in MIX_KNOB_DEFAULTS.items():
            if getattr(args, knob) is None:
                setattr(args, knob, default)
        try:
            import openff.toolkit  # noqa: F401
            import openmm  # noqa: F401
        except ImportError as exc:
            parser.error(f'--mix needs OpenMM + OpenFF in this environment '
                         f'(prep extras): {exc}')
    else:
        stray = [
            flag for flag, val in (
                ('--mix-ps', args.mix_ps),
                ('--mix-settle-steps', args.mix_settle_steps),
                ('--mix-timestep-fs', args.mix_timestep_fs),
                ('--mix-friction-per-ps', args.mix_friction_per_ps),
                ('--mix-platform', args.mix_platform),
                ('--mix-charge-method', args.mix_charge_method),
            ) if val is not None
        ]
        if stray:
            parser.error(f'{"/".join(stray)} given without --mix.')


def mix_config_from_args(args: argparse.Namespace,
                         temperature_K: float) -> MixConfig | None:
    """The workflow's ``MixConfig``, or ``None`` when ``--mix`` is off.

    ``temperature_K`` must be the production thermostat setpoint: the mixing
    thermostat and the post-mixing Maxwell-Boltzmann re-draw both use it
    (decisions.md 2026-07-17 (g)/(i)). Call only after ``resolve_mixing_args``.
    """
    if not args.mix:
        return None
    return MixConfig(
        mix_time_ps=args.mix_ps,
        settle_steps=args.mix_settle_steps,
        temperature_K=temperature_K,
        timestep_fs=args.mix_timestep_fs,
        friction_per_ps=args.mix_friction_per_ps,
        platform=args.mix_platform,
        charge_method=args.mix_charge_method,
    )


def collect_mixing_skips(args: argparse.Namespace, output_dir: Path,
                         log: logging.Logger) -> list[int]:
    """Cycle indices whose classical mixing was gracefully skipped.

    Surfaces graceful mixing skips (WM-P4 de-clash safety net, specs/decisions.md
    追補 2026-07-18): a cycle whose classical mixing diverged even after de-clash
    + minimize + warm-up is skipped (pre-mixing MLIP state kept) rather than
    crashing the run. De-clash should make this near-zero; a non-trivial tally
    (say > 2 of n_cycles) is a RED FLAG that the measurement is compromised, so
    it is recorded here for the analysis to check.

    ``log`` is the caller's logger so the warning keeps the driver's name.
    """
    skipped: list[int] = []
    mixing_log = output_dir / 'mixing.jsonl'
    if args.mix and mixing_log.exists():
        for raw in mixing_log.read_text(encoding='utf-8').splitlines():
            if not raw.strip():
                continue
            try:
                rec = json.loads(raw)
            except json.JSONDecodeError:
                continue
            if rec.get('skipped'):
                skipped.append(rec.get('cycle'))
    if skipped:
        log.warning(
            'Mixing was SKIPPED on %d/%d cycle(s): %s. De-clash should make this '
            'near-zero — investigate the packing/geometry if this is non-trivial.',
            len(skipped), args.n_cycles, skipped,
        )
    return skipped


def mixing_summary_fields(args: argparse.Namespace,
                          skipped_cycles: list[int]) -> dict:
    """The mixing block of a driver's summary.json (decisions.md 2026-08-04).

    Every knob is None when ``--mix`` is off (unresolved sentinels), so a
    non-mixing run's provenance stays explicit rather than silently absent.
    Key order is part of the contract: splice this with ``**`` where the vinyl
    driver already writes these keys so its summary.json stays byte-identical.
    """
    return {
        'mixing_enabled': args.mix,
        'mix_ps': args.mix_ps,
        'mix_settle_steps': args.mix_settle_steps,
        'mix_timestep_fs': args.mix_timestep_fs if args.mix else None,
        'mix_platform': args.mix_platform if args.mix else None,
        'mix_charge_method': args.mix_charge_method if args.mix else None,
        # Graceful mixing-skip tally (WM-P4 de-clash safety net): cycle indices
        # whose classical mixing diverged and was skipped. Empty is the healthy
        # case; a non-trivial count flags a compromised measurement.
        'mixing_skipped_cycles': skipped_cycles,
        'n_mixing_skipped': len(skipped_cycles),
    }
