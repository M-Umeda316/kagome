#!/usr/bin/env python3
'''WM-P5b production-campaign analysis: mix_time_ps sweep x multi-seed.

Design ref: specs/decisions.md '追補 2026-07-19 — WM-P5b 実験計画
(sweep×多シード、ユーザー承認スコープ C)'. The campaign
(scripts/run_wm_p5b_campaign.sh) runs the deterministic-selection arm over
mix_time_ps in {0, 25, 50, 100} (0 = no-mix baseline) x seed in
{7, 11, 17, 23, 42} = 20 runs, each written to
`<runs-dir>/s{seed}_mix{ps}`.

This script answers the two questions the WM-P4 preliminary campaign left
open (specs/decisions.md '追補 2026-07-19 — WM-P4 混合実装 フラットレビュー
と対応', item 4/6):

  Q1 - is the WM-P4 confirmed-yield drop under mixing (A1->A2: 5->2 confirmed
       formations) signal or single-seed noise?
  Q2 - does a shorter mix_time_ps keep the candidate-refresh benefit
       (re-selection rate / Jaccard) while protecting yield?

Per-run metrics are NOT recomputed here -- they are imported directly from
`analyze_wm_p4` (`compute_arm_metrics`, which itself derives the complete-
cycle cap, the re-selection rate, longest stuck-pair chain, mean
cycle-to-cycle Jaccard, confirmed yield, and mixing rms/skip diagnostics from
each run's own jsonl logs). This script's job is purely the seed-aggregation
layer on top: mean +/- standard error per mix_ps, figures, and a report.

Pure data analysis: no orb/GPU/MD is invoked. A run whose directory is
missing or whose selection.jsonl is missing/empty is skipped with a logged
note rather than crashing the rest of the analysis.
'''
from __future__ import annotations

import argparse
import math
import statistics
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

# analyze_wm_p4.py lives next to this file in scripts/; import it directly
# rather than duplicating its metric functions (compute_arm_metrics already
# derives the cycle cap and computes re-selection/Jaccard/stuck-chain/yield/
# mixing-rms from a run's jsonl logs -- see that module's docstring for the
# exact definitions).
sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_wm_p4 import ArmMetrics, compute_arm_metrics  # noqa: E402

DEFAULT_SEEDS = [7, 11, 17, 23, 42]
DEFAULT_MIX_PS = [0, 25, 50, 100]


# --------------------------------------------------------------------------
# Discovery
# --------------------------------------------------------------------------


def run_dir_name(seed: int, mix_ps: int) -> str:
    '''The launcher's naming convention: `s{seed}_mix{ps}` (integer ps, no
    decimal point -- matches scripts/run_wm_p5b_campaign.sh's `${RUN}`).'''
    return f's{seed}_mix{mix_ps}'


def _has_usable_data(run_dir: Path) -> bool:
    '''Same "is there anything to analyze" gate as analyze_wm_p4.main(): the
    run directory must exist and selection.jsonl must exist and be non-empty.
    '''
    sel = run_dir / 'selection.jsonl'
    return run_dir.exists() and sel.exists() and sel.stat().st_size > 0


def load_run_metrics(
    runs_dir: Path,
    seed: int,
    mix_ps: int,
    notes: list[str],
) -> Optional[ArmMetrics]:
    '''Load one run's metrics via analyze_wm_p4.compute_arm_metrics, or
    return None (with a note appended to `notes`) if the run is missing or
    has no usable data. Never raises.'''
    name = run_dir_name(seed, mix_ps)
    run_dir = runs_dir / name
    if not _has_usable_data(run_dir):
        notes.append(f'{name}: missing or empty selection.jsonl -- excluded from aggregation.')
        return None
    try:
        return compute_arm_metrics(name, run_dir)
    except Exception as exc:  # defensive: a malformed log must not crash the campaign analysis
        notes.append(f'{name}: failed to compute metrics ({exc!r}) -- excluded from aggregation.')
        return None


# --------------------------------------------------------------------------
# Aggregation
# --------------------------------------------------------------------------


@dataclass
class MetricAgg:
    '''Mean +/- standard error of one metric across the seeds present for a
    given mix_ps. `se` is None when fewer than 2 seeds are present (standard
    error over a single sample is undefined).'''
    mean: float
    se: Optional[float]
    n: int
    values: list[float] = field(default_factory=list)


def _agg(values: list[float]) -> MetricAgg:
    n = len(values)
    if n == 0:
        return MetricAgg(mean=float('nan'), se=None, n=0, values=[])
    mean = statistics.fmean(values)
    if n < 2:
        return MetricAgg(mean=mean, se=None, n=n, values=list(values))
    se = statistics.stdev(values) / math.sqrt(n)
    return MetricAgg(mean=mean, se=se, n=n, values=list(values))


@dataclass
class MixPsAggregate:
    mix_ps: int
    n_seeds_present: int
    seeds_present: list[int]
    seeds_missing: list[int]
    reselection_rate: MetricAgg
    mean_jaccard: MetricAgg
    longest_stuck_chain: MetricAgg
    yield_per_cycle: MetricAgg
    yield_per_slot: MetricAgg
    mixing_rms_A: MetricAgg  # empty (n=0) for mix_ps == 0
    mixing_skip_rate: MetricAgg  # empty (n=0) for mix_ps == 0


def aggregate_mix_ps(
    mix_ps: int,
    seeds: list[int],
    runs_dir: Path,
    notes: list[str],
) -> MixPsAggregate:
    '''Collect every present seed's metrics for one mix_ps value and reduce
    them to mean +/- SE per metric.'''
    present: dict[int, ArmMetrics] = {}
    for seed in seeds:
        m = load_run_metrics(runs_dir, seed, mix_ps, notes)
        if m is not None:
            present[seed] = m

    seeds_present = sorted(present)
    seeds_missing = sorted(set(seeds) - set(present))

    rms_vals = [m.mean_rms_A for m in present.values() if m.mean_rms_A is not None]
    skip_vals = [m.skip_rate for m in present.values() if m.mixing_present]

    return MixPsAggregate(
        mix_ps=mix_ps,
        n_seeds_present=len(seeds_present),
        seeds_present=seeds_present,
        seeds_missing=seeds_missing,
        reselection_rate=_agg([present[s].reselection_rate for s in seeds_present]),
        mean_jaccard=_agg([present[s].mean_jaccard for s in seeds_present]),
        longest_stuck_chain=_agg(
            [float(present[s].longest_stuck_chain) for s in seeds_present]
        ),
        yield_per_cycle=_agg([present[s].yield_per_cycle for s in seeds_present]),
        yield_per_slot=_agg([present[s].yield_per_slot for s in seeds_present]),
        mixing_rms_A=_agg(rms_vals),
        mixing_skip_rate=_agg(skip_vals),
    )


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------


def _errorbar_series(aggs: list[MixPsAggregate], field_name: str) -> tuple[
    list[int], list[float], list[float]
]:
    xs = [a.mix_ps for a in aggs]
    metric_aggs = [getattr(a, field_name) for a in aggs]
    ys = [m.mean for m in metric_aggs]
    yerr = [m.se if m.se is not None else 0.0 for m in metric_aggs]
    return xs, ys, yerr


def fig_yield_vs_mix_ps(aggs: list[MixPsAggregate], out: Path) -> None:
    '''Headline figure for Q1/Q2: confirmed yield (per cycle) vs mix_ps.'''
    xs, ys, yerr = _errorbar_series(aggs, 'yield_per_cycle')
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(xs, ys, yerr=yerr, fmt='o-', color='#55A868', capsize=4,
                markersize=7, linewidth=1.5)
    for x, y, a in zip(xs, ys, aggs):
        ax.annotate(f'n={a.n_seeds_present}', (x, y), textcoords='offset points',
                    xytext=(6, 6), fontsize=8)
    ax.set_xlabel('mix_time_ps (0 = no-mix baseline)')
    ax.set_ylabel('Confirmed formations per cycle')
    ax.set_title('Confirmed yield vs mix_time_ps (mean +/- SE across seeds)')
    ax.set_xticks(xs)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    _save(fig, out, 'yield_vs_mix_ps')


def fig_reselection_vs_mix_ps(aggs: list[MixPsAggregate], out: Path) -> None:
    xs, ys, yerr = _errorbar_series(aggs, 'reselection_rate')
    ys_pct = [y * 100 for y in ys]
    yerr_pct = [e * 100 for e in yerr]
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(xs, ys_pct, yerr=yerr_pct, fmt='o-', color='#4C72B0', capsize=4,
                markersize=7, linewidth=1.5)
    for x, y, a in zip(xs, ys_pct, aggs):
        ax.annotate(f'n={a.n_seeds_present}', (x, y), textcoords='offset points',
                    xytext=(6, 6), fontsize=8)
    ax.set_xlabel('mix_time_ps (0 = no-mix baseline)')
    ax.set_ylabel('Re-selection rate (% of selected slots)')
    ax.set_title('Consecutive re-selection rate vs mix_time_ps (mean +/- SE)')
    ax.set_xticks(xs)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    _save(fig, out, 'reselection_vs_mix_ps')


def fig_jaccard_vs_mix_ps(aggs: list[MixPsAggregate], out: Path) -> None:
    xs, ys, yerr = _errorbar_series(aggs, 'mean_jaccard')
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(xs, ys, yerr=yerr, fmt='o-', color='#DD8452', capsize=4,
                markersize=7, linewidth=1.5)
    for x, y, a in zip(xs, ys, aggs):
        ax.annotate(f'n={a.n_seeds_present}', (x, y), textcoords='offset points',
                    xytext=(6, 6), fontsize=8)
    ax.set_xlabel('mix_time_ps (0 = no-mix baseline)')
    ax.set_ylabel('Mean cycle-to-cycle candidate-set Jaccard')
    ax.set_title('Candidate-neighbourhood refresh vs mix_time_ps (lower = more refresh)')
    ax.set_xticks(xs)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    _save(fig, out, 'jaccard_vs_mix_ps')


def _save(fig: 'plt.Figure', out: Path, name: str) -> None:
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f'{name}.png', dpi=150)
    fig.savefig(out / f'{name}.pdf')
    plt.close(fig)


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------


def _fmt_mean_se(agg: MetricAgg, digits: int = 3, pct: bool = False) -> str:
    if agg.n == 0:
        return 'n/a (0 seeds)'
    scale = 100.0 if pct else 1.0
    suffix = '%' if pct else ''
    mean_s = f'{agg.mean * scale:.{digits}f}{suffix}'
    if agg.se is None:
        return f'{mean_s} (n={agg.n}, SE n/a)'
    return f'{mean_s} +/- {agg.se * scale:.{digits}f}{suffix} (n={agg.n})'


def _intervals_overlap(a: MetricAgg, b: MetricAgg) -> bool:
    '''Whether a's and b's mean +/- 1 SE intervals overlap. Treated as
    "overlapping" (conservative / not resolved) whenever either side lacks an
    SE (fewer than 2 seeds) -- there is no basis to claim resolution.'''
    if a.n == 0 or b.n == 0:
        return True
    if a.se is None or b.se is None:
        return True
    a_lo, a_hi = a.mean - a.se, a.mean + a.se
    b_lo, b_hi = b.mean - b.se, b.mean + b.se
    return a_lo <= b_hi and b_lo <= a_hi


def _yield_resolution_note(aggs: list[MixPsAggregate]) -> str:
    '''Q1: is the mixing yield difference statistically resolved (to 1 SE)
    relative to the mix_ps=0 baseline? Never overclaims -- overlapping
    intervals (or insufficient seeds) are reported as unresolved.'''
    baseline = next((a for a in aggs if a.mix_ps == 0), None)
    if baseline is None or baseline.n_seeds_present == 0:
        return (
            'No mix_ps=0 baseline data available -- Q1 (signal vs. noise) cannot be '
            'assessed at all.\n'
        )
    lines = [
        f'Baseline mix_ps=0 yield/cycle: {_fmt_mean_se(baseline.yield_per_cycle)}.\n'
    ]
    any_resolved = False
    for a in aggs:
        if a.mix_ps == 0:
            continue
        overlap = _intervals_overlap(baseline.yield_per_cycle, a.yield_per_cycle)
        verdict = (
            'overlaps the baseline (NOT statistically resolved to 1 SE)'
            if overlap
            else 'does NOT overlap the baseline (resolved to 1 SE)'
        )
        if not overlap:
            any_resolved = True
        lines.append(
            f'- mix_ps={a.mix_ps}: yield/cycle {_fmt_mean_se(a.yield_per_cycle)} -- {verdict}.\n'
        )
    if not any_resolved:
        lines.append(
            '\n**Q1 verdict: the yield difference across mix_time_ps is NOT '
            'statistically resolved** at the 1-SE level with the seeds currently '
            'present -- consistent with either a real but small/noisy effect or no '
            'effect at all. Do not claim a confirmed yield penalty from mixing on '
            'this evidence alone.\n'
        )
    else:
        lines.append(
            '\n**Q1 note: at least one mix_ps arm shows a yield difference that does '
            'not overlap the mix_ps=0 baseline at 1 SE.** This is a first statistical '
            'signal, not proof at conventional confidence levels (1 SE is a weak '
            'bar); treat as motivating, not confirmatory.\n'
        )
    return ''.join(lines)


def write_report(
    aggs: list[MixPsAggregate],
    out: Path,
    notes: list[str],
    seeds: list[int],
    mix_ps_values: list[int],
) -> Path:
    lines: list[str] = []
    lines.append('# WM-P5b analysis: mix_time_ps sweep x multi-seed\n')
    lines.append(
        'Generated by `scripts/analyze_wm_p5b.py` (figures + numbers reproducible from '
        'the jsonl logs via `analyze_wm_p4.compute_arm_metrics`; no orb/GPU/MD invoked).\n'
    )
    lines.append(
        f'Requested matrix: mix_time_ps in {mix_ps_values} x seed in {seeds} '
        f'({len(mix_ps_values) * len(seeds)} runs).\n'
    )

    lines.append('\n## Limitations (read first)\n')
    lines.append(
        '- **Single arm: deterministic selection only.** This campaign does not include '
        'the softmax selection policy (WM-P4 A3/A4); it isolates the mixing effect on '
        'the deterministic-selection axis only (specs/decisions.md 追補 2026-07-19 '
        'WM-P5b 実験計画). Findings here should not be generalized to stochastic '
        'selection without a separate campaign.\n'
    )
    lines.append(
        '- **Small system: 40 monomers (20 acrylate + 20 methacrylate + 2 initiators, '
        '~564 atoms), 15 cycles.** Per-run confirmed-formation counts are single-digit; '
        'seed-to-seed variance can be large relative to the mean. Paper-scale (hundreds '
        'of monomers) is explicitly out of scope for this campaign.\n'
    )
    missing_runs = [
        run_dir_name(seed, ps) for ps in mix_ps_values for seed in seeds
        if seed not in next((a.seeds_present for a in aggs if a.mix_ps == ps), [])
    ]
    if missing_runs:
        lines.append(
            f'- **Missing/excluded runs ({len(missing_runs)}/{len(mix_ps_values) * len(seeds)}):** '
            + ', '.join(missing_runs) + '.\n'
        )
    else:
        lines.append('- All 20 planned runs have usable data.\n')
    if notes:
        lines.append('- Per-run notes:\n')
        for note in notes:
            lines.append(f'  - {note}\n')

    lines.append('\n## Per-mix_ps aggregate table (mean +/- SE across seeds)\n')
    lines.append(
        '| mix_time_ps | n seeds | Yield/cycle | Yield/slot | Re-selection rate | '
        'Mean Jaccard | Longest stuck chain | Mixing RMS (A) | Mixing skip rate |\n'
        '|---|---|---|---|---|---|---|---|---|\n'
    )
    for a in aggs:
        lines.append(
            f'| {a.mix_ps} | {a.n_seeds_present}/{len(seeds)} | '
            f'{_fmt_mean_se(a.yield_per_cycle)} | {_fmt_mean_se(a.yield_per_slot, pct=True)} | '
            f'{_fmt_mean_se(a.reselection_rate, pct=True)} | {_fmt_mean_se(a.mean_jaccard)} | '
            f'{_fmt_mean_se(a.longest_stuck_chain, digits=2)} | '
            f'{_fmt_mean_se(a.mixing_rms_A, digits=2) if a.mix_ps != 0 else "n/a (no mixing)"} | '
            f'{_fmt_mean_se(a.mixing_skip_rate, pct=True) if a.mix_ps != 0 else "n/a (no mixing)"} |\n'
        )

    lines.append('\n## Q1 -- is the confirmed-yield drop under mixing signal or noise?\n')
    lines.append(_yield_resolution_note(aggs))

    lines.append(
        '\n## Q2 -- does a shorter mix_time_ps trade off refresh (re-selection/Jaccard) '
        'against yield?\n'
    )
    lines.append(
        'See `reselection_vs_mix_ps` and `jaccard_vs_mix_ps` alongside `yield_vs_mix_ps`: '
        'a favourable mix_time_ps would show materially lower re-selection/Jaccard than '
        'mix_ps=0 without a resolved drop in yield/cycle. Read the three figures together '
        'rather than any single one; per-mix_ps sample sizes above indicate how much '
        'weight each point should be given.\n'
    )

    lines.append('\n## Metric definitions\n')
    lines.append(
        'All metrics are as defined in `scripts/analyze_wm_p4.py` (re-selection rate, '
        'longest stuck-pair chain, cycle-to-cycle candidate-set Jaccard, confirmed yield '
        'per cycle/per selected slot, mixing RMS displacement, mixing skip rate) and are '
        'imported, not recomputed, here. Standard error = sample stdev / sqrt(n_seeds); '
        'reported as "SE n/a" when fewer than 2 seeds are present for a mix_ps.\n'
    )

    text = ''.join(lines)
    out.mkdir(parents=True, exist_ok=True)
    report_path = out / 'WM-P5b-analysis.md'
    report_path.write_text(text, encoding='utf-8')
    return report_path


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description='WM-P5b mix_time_ps sweep x multi-seed analysis')
    parser.add_argument(
        '--runs-dir',
        type=Path,
        default=Path('runs/wm_p5b'),
        help='directory containing the s{seed}_mix{ps} run subdirectories',
    )
    parser.add_argument(
        '--seeds',
        nargs='+',
        type=int,
        default=DEFAULT_SEEDS,
        help=f'seeds in the matrix (default: {DEFAULT_SEEDS})',
    )
    parser.add_argument(
        '--mix-ps',
        nargs='+',
        type=int,
        default=DEFAULT_MIX_PS,
        help=f'mix_time_ps values in the matrix (default: {DEFAULT_MIX_PS})',
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=None,
        help='output dir for figures + report (default: <runs-dir>/analysis)',
    )
    args = parser.parse_args()

    out = args.output_dir or (args.runs_dir / 'analysis')
    notes: list[str] = []

    aggs: list[MixPsAggregate] = []
    for ps in args.mix_ps:
        agg = aggregate_mix_ps(ps, args.seeds, args.runs_dir, notes)
        aggs.append(agg)
        print(
            f'[ok] mix_ps={ps}: n_seeds={agg.n_seeds_present}/{len(args.seeds)} '
            f'yield/cycle={_fmt_mean_se(agg.yield_per_cycle)} '
            f'resel={_fmt_mean_se(agg.reselection_rate, pct=True)} '
            f'jaccard={_fmt_mean_se(agg.mean_jaccard)}'
        )

    for note in notes:
        print(f'[note] {note}', file=sys.stderr)

    if all(a.n_seeds_present == 0 for a in aggs):
        print('No runs with usable data anywhere in the matrix. Nothing to do.')
        return

    fig_yield_vs_mix_ps(aggs, out)
    fig_reselection_vs_mix_ps(aggs, out)
    fig_jaccard_vs_mix_ps(aggs, out)
    report_path = write_report(aggs, out, notes, args.seeds, args.mix_ps)

    print(f'\nFigures + report written to: {out}')
    print(f'Report: {report_path}')


if __name__ == '__main__':
    main()
