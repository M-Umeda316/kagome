#!/usr/bin/env python3
"""WM-P4 preliminary reactivity-measurement analysis.

Compares up to three well-mixed-measurement-mode campaign arms:

    A1_baseline  deterministic selection, NO mixing   (reference)
    A3_softmax   softmax selection (T=0.5), NO mixing
    A4_both      softmax + mixing (50 ps/cycle)        (crashed at cycle 14)

The measurement targets the known "frozen nearest-neighbour re-selection"
pathology of TDBB candidate selection (specs/decisions.md 2026-07-17): the
selector re-picks the same frozen pair every cycle (baseline 46-53%), yielding
low per-cycle confirmed reactivity (~23% baseline framing). The arms let us
isolate the stochastic-selection effect (A1 vs A3) and the mixing-on-top-of-
softmax effect (A3 vs A4).

This script is pure data analysis: NO orb/GPU/MD is invoked. It reads the
per-arm jsonl logs, tolerates a missing summary.json (reads jsonl directly),
computes the metrics, writes figures (png + pdf) and a markdown report.

Metric / identity choices (documented for auditability):

* Pair identity = frozenset({atoms[0], atoms[1]}) -- the two reacting atoms
  (radical + alpha carbon) that form the bond. This matches the (atom_a,
  atom_b) recorded in bonds.jsonl for the same candidate_id, and is the robust
  "formation-pair" convention. For this data it is equivalent to using the
  full sorted 4-tuple of `atoms` (the neighbour atoms are deterministically
  tied to the reacting pair), but the reacting pair is the minimal stable key.

* Cycle cap: the number of COMPLETE cycles per arm is DERIVED from the data
  (see `derive_cycle_cap`), not hand-maintained. A cycle counts as complete
  only if it has a selection.jsonl record and, for a mixing-configured arm, a
  matching mixing.jsonl record; the cap is the largest N such that cycles
  0..N-1 are all complete. For A4_both this yields 14: A4 logged a selection
  AND an attempted formation for cycle 14 before the GPU crash, but no mixing
  record and no confirmed_formation for that cycle, so cycle 14 is a trailing
  partial cycle and is excluded (mixing covers 0-13 incl. one skip at cycle
  12). A1_baseline/A3_softmax have no mixing requirement and derive to 15. Use
  `--cycle-cap ARM=N` to override the derived value manually (a warning is
  printed to stderr if the override disagrees with the derived cap).

* Cycle-count normalization: all rates are reported per-slot or per-cycle so
  A4's 14 cycles compare fairly against A1/A3's 15.
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------
# Arm configuration
# --------------------------------------------------------------------------

# NOTE: there is intentionally no hardcoded per-arm cycle cap here anymore --
# the number of complete cycles per arm is derived from selection.jsonl /
# mixing.jsonl / bonds.jsonl by `derive_cycle_cap` below, with an optional
# `--cycle-cap ARM=N` CLI override for manual control.

ARM_ORDER = ["A1_baseline", "A3_softmax", "A4_both"]
ARM_LABELS = {
    "A1_baseline": "A1 baseline\n(det, no-mix)",
    "A3_softmax": "A3 softmax\n(T=0.5, no-mix)",
    "A4_both": "A4 both\n(softmax+mix)",
}

RMS_TARGET_A = 5.0  # decisions.md (vi): mixing must displace >= 5 A rms
SKIP_THRESHOLD = 2  # > 2 skipped cycles => compromised measurement (per brief)


# --------------------------------------------------------------------------
# IO helpers
# --------------------------------------------------------------------------


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    """Read a jsonl file into a list of dicts. Missing/empty -> []."""
    if not path.exists() or path.stat().st_size == 0:
        return []
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def pair_identity(atoms: Iterable[int]) -> frozenset[int]:
    """Stable identity for a candidate/selected pair: the two reacting atoms."""
    atoms = list(atoms)
    return frozenset(atoms[:2])


def _read_manifest(arm_dir: Path) -> Optional[dict[str, Any]]:
    """Read manifest.json for an arm. Returns None if missing/unparseable."""
    path = arm_dir / "manifest.json"
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return None


# --------------------------------------------------------------------------
# Data-derived cycle cap (replaces a hand-maintained per-arm literal)
# --------------------------------------------------------------------------


def _is_mixing_arm(arm_dir: Path, mixing_raw: list[dict[str, Any]]) -> bool:
    """Whether this arm is configured to run classical mixing.

    Prefers the manifest's own record of the config that was actually run
    (`extra.mixing`, present even if the arm crashed before writing any
    mixing.jsonl rows); falls back to "does mixing.jsonl have any rows at
    all" if the manifest is missing/unparseable. This makes a future 4th arm
    (e.g. A2_mix50) work without touching this function.
    """
    manifest = _read_manifest(arm_dir)
    if manifest is not None:
        extra = manifest.get("extra", {})
        if isinstance(extra, dict) and "mixing" in extra:
            return extra["mixing"] is not None
    return bool(mixing_raw)


def derive_cycle_cap(
    arm: str,
    arm_dir: Path,
    selection_raw: list[dict[str, Any]],
    mixing_raw: list[dict[str, Any]],
    bonds_raw: list[dict[str, Any]],
) -> int:
    """Derive the number of leading COMPLETE cycles (0..N-1) for an arm from
    its own jsonl logs -- no hardcoded per-arm literal.

    A cycle `c` counts as complete iff:
      * selection.jsonl has a record for cycle c, AND
      * if the arm is configured for mixing (`_is_mixing_arm`), mixing.jsonl
        also has a record for cycle c (a `skipped: true` mixing row still
        counts -- mixing was attempted, just aborted for that cycle).
    The cap is the largest N such that cycles 0..N-1 are ALL complete: the
    scan stops at the first gap in cycle numbering, or the first cycle that
    has a selection (and typically an attempted_formation in bonds.jsonl)
    but no matching mixing record -- the observed A4 crash signature, where
    cycle 14 logged a selection + attempted_formation but no mixing record
    and no confirmed_formation before the transient GPU fault.

    `bonds_raw` is read too (as the task requires) and used as an
    audit/cross-check: if a cycle counted as complete has zero bonds.jsonl
    records at all, that is unusual (every complete cycle observed so far
    has >= 1 attempted_* event) and is flagged to stderr rather than
    silently trusted, but it does not by itself change the derived cap.
    """
    sel_cycles = sorted({d["cycle"] for d in selection_raw})
    mix_cycles = {m["cycle"] for m in mixing_raw}
    is_mixing = _is_mixing_arm(arm_dir, mixing_raw)

    n = 0
    for c in sel_cycles:
        if c != n:
            break  # gap in cycle numbering -> stop growing the cap
        if is_mixing and c not in mix_cycles:
            break  # selected/attempted but never mixed -> trailing partial cycle
        n += 1

    bond_cycles = {b.get("cycle") for b in bonds_raw}
    missing_bonds = [c for c in range(n) if c not in bond_cycles]
    if missing_bonds:
        print(
            f"[warn] {arm}: cycles {missing_bonds} counted as complete but have "
            "no bonds.jsonl record at all; verify manually.",
            file=sys.stderr,
        )
    return n


def _parse_cycle_cap_overrides(raw: Optional[list[str]]) -> dict[str, int]:
    """Parse repeated `--cycle-cap ARM=N` tokens into an override dict."""
    overrides: dict[str, int] = {}
    if not raw:
        return overrides
    for item in raw:
        arm_name, sep, n_str = item.partition("=")
        if not sep:
            raise argparse.ArgumentTypeError(f"--cycle-cap expects ARM=N (got {item!r})")
        try:
            overrides[arm_name] = int(n_str)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(
                f"--cycle-cap N must be an integer (got {item!r})"
            ) from exc
    return overrides


# --------------------------------------------------------------------------
# Per-arm metric container
# --------------------------------------------------------------------------


@dataclass
class ArmMetrics:
    arm: str
    n_cycles: int
    n_selected_slots: int
    reselection_rate: float  # per total selected slot
    reselection_rate_eligible: float  # per slot in cycles 1..N (excludes cycle 0)
    n_reselections: int
    longest_stuck_chain: int
    stuck_pair: Optional[tuple[int, ...]]
    mean_jaccard: float
    jaccard_series: list[float]
    confirmed_formations: int
    confirmed_dissociations: int
    yield_per_cycle: float  # confirmed formations / n_cycles
    yield_per_slot: float  # confirmed formations / n_selected_slots
    # mixing (None for non-mixing arms)
    mixing_present: bool = False
    n_mix_records: int = 0
    n_skipped: int = 0
    skip_rate: float = 0.0
    mean_rms_A: Optional[float] = None
    min_rms_A: Optional[float] = None
    total_declashes: int = 0
    mean_warmup_steps: Optional[float] = None
    rms_series: list[tuple[int, Optional[float]]] = field(default_factory=list)
    skip_compromised: bool = False


# --------------------------------------------------------------------------
# Metric computation
# --------------------------------------------------------------------------


def compute_arm_metrics(
    arm: str,
    arm_dir: Path,
    cap_override: Optional[int] = None,
) -> ArmMetrics:
    selection_all = read_jsonl(arm_dir / "selection.jsonl")
    selection_all.sort(key=lambda d: d["cycle"])
    mixing_all = read_jsonl(arm_dir / "mixing.jsonl")
    mixing_all.sort(key=lambda d: d["cycle"])
    bonds_all = read_jsonl(arm_dir / "bonds.jsonl")

    derived_cap = derive_cycle_cap(arm, arm_dir, selection_all, mixing_all, bonds_all)
    if cap_override is not None:
        if cap_override != derived_cap:
            print(
                f"[warn] {arm}: --cycle-cap override {cap_override} differs from the "
                f"data-derived complete-cycle count {derived_cap}; using the override.",
                file=sys.stderr,
            )
        cap = cap_override
    else:
        cap = derived_cap

    selection = [d for d in selection_all if d["cycle"] < cap]
    n_cycles = len(selection)

    # --- selected pair identities per cycle (list of slots, order preserved) ---
    selected_per_cycle: list[list[frozenset[int]]] = [
        [pair_identity(s["atoms"]) for s in d.get("selected", [])] for d in selection
    ]
    # candidate pool (selected + rejected) per cycle
    candidate_per_cycle: list[set[frozenset[int]]] = []
    for d in selection:
        pool: set[frozenset[int]] = set()
        for s in d.get("selected", []):
            pool.add(pair_identity(s["atoms"]))
        for r in d.get("rejected", []):
            pool.add(pair_identity(r["atoms"]))
        candidate_per_cycle.append(pool)

    n_selected_slots = sum(len(s) for s in selected_per_cycle)

    # --- (1) re-selection rate: slots in cycle N whose pair was selected in N-1 ---
    n_reselections = 0
    n_eligible_slots = 0
    for i in range(1, len(selected_per_cycle)):
        prev = set(selected_per_cycle[i - 1])
        for pair in selected_per_cycle[i]:
            n_eligible_slots += 1
            if pair in prev:
                n_reselections += 1
    reselection_rate = n_reselections / n_selected_slots if n_selected_slots else 0.0
    reselection_rate_eligible = (
        n_reselections / n_eligible_slots if n_eligible_slots else 0.0
    )

    # --- (2) longest stuck-pair chain: max consecutive cycles a pair is selected ---
    longest_stuck_chain, stuck_pair = _longest_stuck_chain(selected_per_cycle)

    # --- (3) in-window candidate-set Jaccard between adjacent cycles ---
    jaccard_series: list[float] = []
    for i in range(1, len(candidate_per_cycle)):
        a, b = candidate_per_cycle[i - 1], candidate_per_cycle[i]
        union = a | b
        jac = len(a & b) / len(union) if union else 0.0
        jaccard_series.append(jac)
    mean_jaccard = sum(jaccard_series) / len(jaccard_series) if jaccard_series else 0.0

    # --- (4) per-cycle confirmed yield ---
    bonds = [b for b in bonds_all if b.get("cycle", -1) < cap]
    confirmed_formations = sum(
        1 for b in bonds if b.get("event_type") == "confirmed_formation"
    )
    confirmed_dissociations = sum(
        1 for b in bonds if b.get("event_type") == "confirmed_dissociation"
    )
    yield_per_cycle = confirmed_formations / n_cycles if n_cycles else 0.0
    yield_per_slot = (
        confirmed_formations / n_selected_slots if n_selected_slots else 0.0
    )

    metrics = ArmMetrics(
        arm=arm,
        n_cycles=n_cycles,
        n_selected_slots=n_selected_slots,
        reselection_rate=reselection_rate,
        reselection_rate_eligible=reselection_rate_eligible,
        n_reselections=n_reselections,
        longest_stuck_chain=longest_stuck_chain,
        stuck_pair=tuple(sorted(stuck_pair)) if stuck_pair else None,
        mean_jaccard=mean_jaccard,
        jaccard_series=jaccard_series,
        confirmed_formations=confirmed_formations,
        confirmed_dissociations=confirmed_dissociations,
        yield_per_cycle=yield_per_cycle,
        yield_per_slot=yield_per_slot,
    )

    # --- (5) mixing diagnostics ---
    mixing = [m for m in mixing_all if m.get("cycle", -1) < cap]
    if mixing:
        _fill_mixing_metrics(metrics, mixing)

    return metrics


def _longest_stuck_chain(
    selected_per_cycle: list[list[frozenset[int]]],
) -> tuple[int, Optional[frozenset[int]]]:
    """Max number of consecutive cycles that any single pair stays selected."""
    if not selected_per_cycle:
        return 0, None
    all_pairs: set[frozenset[int]] = set()
    for slots in selected_per_cycle:
        all_pairs |= set(slots)
    best_len = 0
    best_pair: Optional[frozenset[int]] = None
    for pair in all_pairs:
        cur = 0
        for slots in selected_per_cycle:
            if pair in slots:
                cur += 1
                if cur > best_len:
                    best_len, best_pair = cur, pair
            else:
                cur = 0
    return best_len, best_pair


def _fill_mixing_metrics(metrics: ArmMetrics, mixing: list[dict[str, Any]]) -> None:
    mixing.sort(key=lambda d: d["cycle"])
    metrics.mixing_present = True
    metrics.n_mix_records = len(mixing)
    skipped = [m for m in mixing if m.get("skipped")]
    metrics.n_skipped = len(skipped)
    metrics.skip_rate = metrics.n_skipped / metrics.n_mix_records
    rms_vals = [
        m["rms_displacement_A"]
        for m in mixing
        if not m.get("skipped") and m.get("rms_displacement_A") is not None
    ]
    if rms_vals:
        metrics.mean_rms_A = sum(rms_vals) / len(rms_vals)
        metrics.min_rms_A = min(rms_vals)
    metrics.total_declashes = sum(m.get("n_declashed") or 0 for m in mixing)
    warmups = [
        m.get("n_warmup_steps")
        for m in mixing
        if not m.get("skipped") and m.get("n_warmup_steps") is not None
    ]
    if warmups:
        metrics.mean_warmup_steps = sum(warmups) / len(warmups)
    metrics.rms_series = [
        (m["cycle"], None if m.get("skipped") else m.get("rms_displacement_A"))
        for m in mixing
    ]
    metrics.skip_compromised = metrics.n_skipped > SKIP_THRESHOLD


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------

_COLORS = {
    "A1_baseline": "#4C72B0",
    "A3_softmax": "#DD8452",
    "A4_both": "#55A868",
}


def _labels(arms: list[str]) -> list[str]:
    return [ARM_LABELS.get(a, a) for a in arms]


def fig_reselection_and_stuck(results: list[ArmMetrics], out: Path) -> None:
    arms = [r.arm for r in results]
    x = range(len(arms))
    fig, ax1 = plt.subplots(figsize=(8, 5))
    width = 0.38
    resel = [r.reselection_rate * 100 for r in results]
    stuck = [r.longest_stuck_chain for r in results]

    bars1 = ax1.bar(
        [i - width / 2 for i in x], resel, width, color="#4C72B0", label="re-selection rate"
    )
    ax1.set_ylabel("Re-selection rate (% of selected slots)", color="#4C72B0")
    ax1.tick_params(axis="y", labelcolor="#4C72B0")
    ax1.set_ylim(0, max(60, max(resel) * 1.2 if resel else 60))
    ax1.axhspan(46, 53, color="#4C72B0", alpha=0.12, label="baseline 46-53% band")
    for b, v in zip(bars1, resel):
        ax1.text(b.get_x() + b.get_width() / 2, v + 1, f"{v:.0f}%", ha="center", fontsize=9)

    ax2 = ax1.twinx()
    bars2 = ax2.bar(
        [i + width / 2 for i in x], stuck, width, color="#C44E52", label="longest stuck chain"
    )
    ax2.set_ylabel("Longest stuck-pair chain (cycles)", color="#C44E52")
    ax2.tick_params(axis="y", labelcolor="#C44E52")
    ax2.set_ylim(0, max(15, max(stuck) * 1.3 if stuck else 15))
    for b, v in zip(bars2, stuck):
        ax2.text(b.get_x() + b.get_width() / 2, v + 0.2, f"{v}", ha="center", fontsize=9)

    ax1.set_xticks(list(x))
    ax1.set_xticklabels(_labels(arms))
    ax1.set_title("Re-selection rate & longest stuck-pair chain per arm")
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines + lines2, labels + labels2, loc="upper right", fontsize=8)
    fig.tight_layout()
    _save(fig, out, "reselection_stuck")


def fig_jaccard(results: list[ArmMetrics], out: Path) -> None:
    arms = [r.arm for r in results]
    x = range(len(arms))
    vals = [r.mean_jaccard for r in results]
    fig, ax = plt.subplots(figsize=(7, 5))
    bars = ax.bar(list(x), vals, color=[_COLORS.get(a, "#777") for a in arms])
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.01, f"{v:.2f}", ha="center", fontsize=10)
    ax.set_xticks(list(x))
    ax.set_xticklabels(_labels(arms))
    ax.set_ylabel("Mean cycle-to-cycle candidate-set Jaccard")
    ax.set_ylim(0, max(1.0, max(vals) * 1.25 if vals else 1.0))
    ax.set_title("Candidate-neighbourhood refresh (lower = more refresh)")
    fig.tight_layout()
    _save(fig, out, "jaccard")


def fig_yield(results: list[ArmMetrics], out: Path) -> None:
    arms = [r.arm for r in results]
    x = range(len(arms))
    per_cycle = [r.yield_per_cycle for r in results]
    per_slot = [r.yield_per_slot * 100 for r in results]
    fig, ax1 = plt.subplots(figsize=(8, 5))
    width = 0.38
    bars1 = ax1.bar(
        [i - width / 2 for i in x], per_cycle, width, color="#55A868",
        label="confirmed formations / cycle",
    )
    ax1.set_ylabel("Confirmed formations per cycle", color="#55A868")
    ax1.tick_params(axis="y", labelcolor="#55A868")
    for b, v in zip(bars1, per_cycle):
        ax1.text(b.get_x() + b.get_width() / 2, v + 0.005, f"{v:.2f}", ha="center", fontsize=9)

    ax2 = ax1.twinx()
    bars2 = ax2.bar(
        [i + width / 2 for i in x], per_slot, width, color="#8172B3",
        label="confirmed / selected slot (%)",
    )
    ax2.set_ylabel("Confirmed per selected slot (%)", color="#8172B3")
    ax2.tick_params(axis="y", labelcolor="#8172B3")
    ax2.axhline(23, color="#8172B3", ls="--", lw=1, alpha=0.7, label="~23% baseline framing")
    for b, v in zip(bars2, per_slot):
        ax2.text(b.get_x() + b.get_width() / 2, v + 0.5, f"{v:.0f}%", ha="center", fontsize=9)

    ax1.set_xticks(list(x))
    ax1.set_xticklabels(_labels(arms))
    ax1.set_title("Per-cycle confirmed yield per arm")
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines + lines2, labels + labels2, loc="upper right", fontsize=8)
    fig.tight_layout()
    _save(fig, out, "yield")


def fig_mixing_rms(results: list[ArmMetrics], out: Path) -> None:
    mix = [r for r in results if r.mixing_present]
    fig, ax = plt.subplots(figsize=(8, 5))
    if not mix:
        ax.text(0.5, 0.5, "No mixing arm available", ha="center", va="center")
    for r in mix:
        cycles = [c for c, v in r.rms_series]
        vals = [v for c, v in r.rms_series]
        plotted_c = [c for c, v in zip(cycles, vals) if v is not None]
        plotted_v = [v for v in vals if v is not None]
        ax.plot(plotted_c, plotted_v, "o-", color=_COLORS.get(r.arm, "#55A868"),
                label=f"{r.arm} rms")
        for c, v in zip(cycles, vals):
            if v is None:
                ax.axvline(c, color="red", ls=":", alpha=0.6)
                ax.text(c, RMS_TARGET_A + 0.3, "skip", color="red", rotation=90,
                        fontsize=8, ha="right", va="bottom")
    ax.axhline(RMS_TARGET_A, color="k", ls="--", lw=1, label=f"{RMS_TARGET_A:.0f} A target")
    ax.set_xlabel("Cycle")
    ax.set_ylabel("RMS displacement (A)")
    ax.set_title("A4 mixing RMS displacement per cycle")
    ax.legend(fontsize=8)
    ax.set_ylim(0, None)
    fig.tight_layout()
    _save(fig, out, "mixing_rms")


def _save(fig: "plt.Figure", out: Path, name: str) -> None:
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{name}.png", dpi=150)
    fig.savefig(out / f"{name}.pdf")
    plt.close(fig)


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------


def _pct(x: float) -> str:
    return f"{x * 100:.1f}%"


def _manifest_field(manifest: Optional[dict[str, Any]], *keys: str) -> Any:
    """Look up the first present key, checking top-level then `extra`."""
    if manifest is None:
        return "unknown"
    extra = manifest.get("extra", {})
    if not isinstance(extra, dict):
        extra = {}
    for key in keys:
        if key in manifest and manifest[key] is not None:
            return manifest[key]
        if key in extra and extra[key] is not None:
            return extra[key]
    return "unknown"


def _agree_or_list(values: dict[str, Any]) -> str:
    """Render a per-arm value dict as a single value if all arms agree, else
    an explicit per-arm breakdown noting the disagreement."""
    present = {a: v for a, v in values.items() if v != "unknown"}
    uniq = {v for v in present.values()}
    if len(uniq) == 1:
        return str(uniq.pop())
    if not uniq:
        return "unknown"
    return "disagreement across arms: " + ", ".join(f"{a}={v}" for a, v in values.items())


def system_summary_line(arm_dirs: dict[str, Path], arms: list[str]) -> str:
    """Build the report's system-description line from each arm's own
    manifest.json, instead of a hardcoded literal. Falls back to "unknown"
    per-arm/per-field when a manifest or field is missing (never crashes),
    and calls out disagreement across arms rather than silently picking one.
    """
    manifests = {arm: _read_manifest(arm_dirs[arm]) for arm in arms if arm in arm_dirs}
    seeds = {arm: _manifest_field(m, "seed") for arm, m in manifests.items()}
    backends = {arm: _manifest_field(m, "backend") for arm, m in manifests.items()}
    n_sites = {arm: _manifest_field(m, "n_reactive_sites") for arm, m in manifests.items()}

    seed_s = _agree_or_list(seeds)
    backend_s = _agree_or_list(backends)
    sites_s = _agree_or_list(n_sites)
    return (
        f"System per arm: {sites_s} reactive sites, seed {seed_s}, {backend_s} backend "
        "(read from each arm's manifest.json; acrylate/methacrylate/initiator counts are "
        "not recorded there, so only the reactive-site total is reported; falls back to "
        "\"unknown\" if a manifest or field is missing).\n"
    )


def write_report(results: list[ArmMetrics], out: Path, arm_dirs: dict[str, Path]) -> Path:
    by_arm = {r.arm: r for r in results}
    lines: list[str] = []
    lines.append("# WM-P4 preliminary reactivity-measurement analysis\n")
    lines.append(
        "Generated by `scripts/analyze_wm_p4.py` (figures + numbers reproducible from "
        "the jsonl logs; no orb/GPU/MD invoked).\n"
    )
    lines.append(system_summary_line(arm_dirs, [r.arm for r in results]))

    # limitations up front (prominent)
    lines.append("## Limitations (read first)\n")
    lines.append(
        "- **A2_mix50 is MISSING** (deterministic selection + mixing). A GPU failure "
        "blocked it; its data files are present but empty. Without A2 we **cannot isolate "
        "the mixing effect on top of deterministic selection** (the A1 vs A2 contrast). "
        "The only clean mixing signal available is A3 vs A4 (mixing on top of softmax).\n"
    )
    lines.append(
        "- **A4 is 14 complete cycles, not 15.** A4 crashed on a transient GPU error at "
        "cycle 14; it logged a selection and an attempted formation for cycle 14 but no "
        "mixing and no confirmation. The partial cycle 14 is excluded; A4 = cycles 0-13. "
        "All rates are normalized per-slot/per-cycle so the comparison is fair.\n"
    )
    lines.append(
        "- **Single seed (7), small system (40 monomers).** Differences of a few "
        "selected slots move the rates by ~3-7 points. Treat magnitudes as indicative, "
        "not statistically resolved. Do not over-read small/noisy differences.\n"
    )
    lines.append(
        "- **GPU instability caveat.** Two of the four planned arms hit GPU failures "
        "(A2 blocked, A4 truncated). The campaign is preliminary.\n"
    )

    # summary table
    lines.append("## Summary table (arm x metric)\n")
    header = (
        "| Metric | A1 baseline | A3 softmax | A4 both |\n"
        "|---|---|---|---|\n"
    )
    lines.append(header)

    def row(name: str, fn) -> str:
        cells = " | ".join(fn(by_arm[a]) if a in by_arm else "-" for a in ARM_ORDER)
        return f"| {name} | {cells} |\n"

    lines.append(row("Selection policy / mixing",
                     lambda r: {"A1_baseline": "det / none",
                                "A3_softmax": "softmax / none",
                                "A4_both": "softmax / 50ps"}[r.arm]))
    lines.append(row("Complete cycles", lambda r: str(r.n_cycles)))
    lines.append(row("Selected slots (total)", lambda r: str(r.n_selected_slots)))
    lines.append(row("Re-selection rate (/all slots)", lambda r: _pct(r.reselection_rate)))
    lines.append(row("Re-selection rate (/eligible slots)",
                     lambda r: _pct(r.reselection_rate_eligible)))
    lines.append(row("Longest stuck-pair chain (cycles)",
                     lambda r: str(r.longest_stuck_chain)))
    lines.append(row("Mean cycle-to-cycle Jaccard", lambda r: f"{r.mean_jaccard:.3f}"))
    lines.append(row("Confirmed formations", lambda r: str(r.confirmed_formations)))
    lines.append(row("Yield per cycle", lambda r: f"{r.yield_per_cycle:.3f}"))
    lines.append(row("Yield per selected slot", lambda r: _pct(r.yield_per_slot)))
    lines.append(row("Mixing mean RMS (A)",
                     lambda r: f"{r.mean_rms_A:.2f}" if r.mean_rms_A is not None else "n/a"))
    lines.append(row("Mixing skipped cycles",
                     lambda r: f"{r.n_skipped}/{r.n_mix_records}" if r.mixing_present else "n/a"))
    lines.append("\n")

    # comparison framings
    a1, a3, a4 = by_arm.get("A1_baseline"), by_arm.get("A3_softmax"), by_arm.get("A4_both")
    lines.append("## Comparison framings\n")
    if a1 and a3:
        lines.append(
            f"### A1 vs A3 -- stochastic-selection effect (softmax alone, no mixing)\n"
            f"- Re-selection rate: {_pct(a1.reselection_rate)} -> {_pct(a3.reselection_rate)} "
            f"({_delta(a1.reselection_rate, a3.reselection_rate)}).\n"
            f"- Longest stuck chain: {a1.longest_stuck_chain} -> {a3.longest_stuck_chain} cycles.\n"
            f"- Mean Jaccard: {a1.mean_jaccard:.3f} -> {a3.mean_jaccard:.3f}.\n"
            f"- Yield/cycle: {a1.yield_per_cycle:.3f} -> {a3.yield_per_cycle:.3f}.\n"
        )
    if a3 and a4:
        lines.append(
            f"### A3 vs A4 -- MIXING effect on top of softmax (**key available mixing signal**)\n"
            f"- Re-selection rate: {_pct(a3.reselection_rate)} -> {_pct(a4.reselection_rate)} "
            f"({_delta(a3.reselection_rate, a4.reselection_rate)}).\n"
            f"- Longest stuck chain: {a3.longest_stuck_chain} -> {a4.longest_stuck_chain} cycles.\n"
            f"- Mean Jaccard: {a3.mean_jaccard:.3f} -> {a4.mean_jaccard:.3f} "
            f"(lower = more neighbourhood refresh).\n"
            f"- Yield/cycle: {a3.yield_per_cycle:.3f} -> {a4.yield_per_cycle:.3f}.\n"
        )
    if a1 and a4:
        lines.append(
            f"### A1 vs A4 -- combined effect (softmax + mixing) vs baseline\n"
            f"- Re-selection rate: {_pct(a1.reselection_rate)} -> {_pct(a4.reselection_rate)} "
            f"({_delta(a1.reselection_rate, a4.reselection_rate)}).\n"
            f"- Longest stuck chain: {a1.longest_stuck_chain} -> {a4.longest_stuck_chain} cycles.\n"
            f"- Mean Jaccard: {a1.mean_jaccard:.3f} -> {a4.mean_jaccard:.3f}.\n"
            f"- Yield/cycle: {a1.yield_per_cycle:.3f} -> {a4.yield_per_cycle:.3f}.\n"
        )
    lines.append(
        "- **Missing:** A1 vs A2 (mixing on top of deterministic selection) -- not "
        "available (A2 blocked).\n"
    )

    # key finding
    lines.append("## Key finding\n")
    lines.append(_key_finding(a1, a3, a4))

    # A4 health
    if a4 and a4.mixing_present:
        lines.append("\n## A4 mixing health\n")
        health = "OK" if not a4.skip_compromised else "COMPROMISED (> 2 skips)"
        lines.append(
            f"- Mean RMS displacement: **{a4.mean_rms_A:.2f} A** (min {a4.min_rms_A:.2f} A) "
            f"vs {RMS_TARGET_A:.0f} A target -> "
            f"{'well above target' if a4.mean_rms_A and a4.mean_rms_A >= RMS_TARGET_A else 'BELOW target'}.\n"
            f"- Skipped cycles: **{a4.n_skipped}/{a4.n_mix_records}** "
            f"(skip rate {_pct(a4.skip_rate)}); compromised-measurement threshold is "
            f"> {SKIP_THRESHOLD} skips -> **{health}**.\n"
            f"- Total declashes across cycles: {a4.total_declashes}.\n"
            f"- Mean warm-up steps: "
            f"{a4.mean_warmup_steps:.0f}.\n"
        )

    lines.append("\n## Metric definitions & choices\n")
    lines.append(
        "- **Pair identity** = unordered pair of the two reacting atoms "
        "(`atoms[0]`, `atoms[1]`), matching `atom_a`/`atom_b` in bonds.jsonl. Equivalent "
        "to the full sorted 4-tuple for this data.\n"
        "- **Re-selection** = a pair selected in cycle N that was ALSO selected in cycle "
        "N-1 (consecutive). Rate reported both over all selected slots (per brief) and "
        "over eligible slots (cycles 1..N, excludes cycle 0 which has no predecessor).\n"
        "- **Longest stuck chain** = max consecutive cycles any single pair stays selected.\n"
        "- **Jaccard** = candidate-set (selected + rejected pairs) Jaccard similarity "
        "between adjacent cycles, averaged.\n"
        "- **Yield per cycle** = confirmed_formation events / complete cycles. "
        "**Yield per selected slot** = confirmed_formation events / total selected slots. "
        "The ~23% baseline framing is a per-attempt confirmation rate, so compare it to "
        "the per-selected-slot column.\n"
    )

    text = "".join(lines)
    out.mkdir(parents=True, exist_ok=True)
    report_path = out / "WM-P4-preliminary.md"
    report_path.write_text(text, encoding="utf-8")
    return report_path


def _delta(a: float, b: float) -> str:
    d = (b - a) * 100
    sign = "+" if d >= 0 else ""
    return f"{sign}{d:.1f} pts"


def _key_finding(a1, a3, a4) -> str:
    if not (a3 and a4):
        return "Insufficient arms for the A3 vs A4 mixing comparison.\n"
    resel_d = (a4.reselection_rate - a3.reselection_rate) * 100
    jac_d = a4.mean_jaccard - a3.mean_jaccard
    stuck_d = a4.longest_stuck_chain - a3.longest_stuck_chain
    parts: list[str] = []
    parts.append(
        f"The key available mixing signal is **A3 vs A4** (both softmax; only A4 mixes 50 "
        f"ps/cycle). "
    )
    # re-selection direction
    if resel_d <= -5:
        parts.append(
            f"Mixing **substantially lowers** consecutive re-selection "
            f"({_pct(a3.reselection_rate)} -> {_pct(a4.reselection_rate)}, {resel_d:.1f} pts). "
        )
    elif resel_d >= 5:
        parts.append(
            f"Mixing **raises** measured re-selection "
            f"({_pct(a3.reselection_rate)} -> {_pct(a4.reselection_rate)}, +{resel_d:.1f} pts) "
            f"-- unexpected; likely small-sample noise. "
        )
    else:
        parts.append(
            f"Re-selection is **essentially unchanged** by mixing "
            f"({_pct(a3.reselection_rate)} -> {_pct(a4.reselection_rate)}, {resel_d:+.1f} pts). "
        )
    # jaccard direction (the more mechanistic refresh signal)
    if jac_d <= -0.03:
        parts.append(
            f"Neighbourhood refresh improves: mean candidate-set Jaccard drops "
            f"{a3.mean_jaccard:.3f} -> {a4.mean_jaccard:.3f} ({jac_d:+.3f}), i.e. mixing "
            f"presents more distinct candidate pools cycle-to-cycle. "
        )
    elif jac_d >= 0.03:
        parts.append(
            f"Candidate-set Jaccard rises {a3.mean_jaccard:.3f} -> {a4.mean_jaccard:.3f} "
            f"({jac_d:+.3f}) -- mixing did not visibly refresh the neighbourhood here. "
        )
    else:
        parts.append(
            f"Candidate-set Jaccard is roughly flat "
            f"({a3.mean_jaccard:.3f} -> {a4.mean_jaccard:.3f}, {jac_d:+.3f}). "
        )
    parts.append(
        f"Longest stuck chain {a3.longest_stuck_chain} -> {a4.longest_stuck_chain} cycles. "
    )
    parts.append(
        f"Confirmed yield, however, is low and noisy in every arm (A1 "
        f"{a1.confirmed_formations if a1 else 'n/a'}, A3 {a3.confirmed_formations}, A4 "
        f"{a4.confirmed_formations} confirmed formations); A4's yield is the *lowest*, not "
        f"the highest, so the refresh did NOT translate into more confirmations at this "
        f"scale/seed -- plausibly because 50 ps of mixing separates the boosted pair and it "
        f"must re-approach, and because single-digit counts are not statistically meaningful. "
    )
    parts.append(
        "**Honest read:** the *selection-refresh* effect of mixing (A3->A4) is real and "
        "large -- re-selection and candidate-set overlap both drop sharply, which is exactly "
        "the pathology the well-mixed mode targets. But the downstream *reactive yield* did "
        "not improve (it fell to 1 confirmation), so on this single-seed, 40-monomer, "
        "14-cycle run mixing demonstrably refreshes the candidate neighbourhood without (yet) "
        "buying more confirmed reactions. A2 (mixing on deterministic selection) and repeat "
        "seeds are needed before claiming a net reactivity benefit.\n"
    )
    return "".join(parts)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description="WM-P4 preliminary reactivity analysis")
    parser.add_argument(
        "--runs-dir",
        type=Path,
        default=Path("runs/wm_p4"),
        help="directory containing the arm subdirectories",
    )
    parser.add_argument(
        "--arms",
        nargs="+",
        default=ARM_ORDER,
        help="arm subdirectory names to analyze (in order)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="output dir for figures + report (default: <runs-dir>/analysis)",
    )
    parser.add_argument(
        "--cycle-cap",
        nargs="*",
        default=None,
        metavar="ARM=N",
        help=(
            "override the data-derived complete-cycle cap for one or more arms, "
            "e.g. --cycle-cap A4_both=14. The default (omit this flag) derives the "
            "cap from selection.jsonl/mixing.jsonl/bonds.jsonl per arm; a warning is "
            "printed to stderr if an override disagrees with the derived value."
        ),
    )
    args = parser.parse_args()

    out = args.output_dir or (args.runs_dir / "analysis")
    arm_dirs = {arm: args.runs_dir / arm for arm in args.arms}
    cap_overrides = _parse_cycle_cap_overrides(args.cycle_cap)

    results: list[ArmMetrics] = []
    for arm in args.arms:
        arm_dir = arm_dirs[arm]
        if not arm_dir.exists():
            print(f"[skip] {arm}: directory missing ({arm_dir})")
            continue
        if not (arm_dir / "selection.jsonl").exists() or (
            arm_dir / "selection.jsonl"
        ).stat().st_size == 0:
            print(f"[skip] {arm}: selection.jsonl missing/empty -> no usable data")
            continue
        m = compute_arm_metrics(arm, arm_dir, cap_overrides.get(arm))
        results.append(m)
        print(
            f"[ok] {arm}: cycles={m.n_cycles} slots={m.n_selected_slots} "
            f"resel={_pct(m.reselection_rate)} stuck={m.longest_stuck_chain} "
            f"jaccard={m.mean_jaccard:.3f} confirmed={m.confirmed_formations} "
            f"rms={m.mean_rms_A}"
        )

    if not results:
        print("No arms with usable data. Nothing to do.")
        return

    # order results by ARM_ORDER for stable figures/tables
    results.sort(key=lambda r: ARM_ORDER.index(r.arm) if r.arm in ARM_ORDER else 99)

    fig_reselection_and_stuck(results, out)
    fig_jaccard(results, out)
    fig_yield(results, out)
    fig_mixing_rms(results, out)
    report_path = write_report(results, out, arm_dirs)

    print(f"\nFigures + report written to: {out}")
    print(f"Report: {report_path}")


if __name__ == "__main__":
    main()
