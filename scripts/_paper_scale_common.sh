#!/usr/bin/env bash
# Shared boilerplate for the paper-scale run launchers (epoxy-amine, nylon66,
# S6 vinyl, vinyl copolymer). Source from a launcher (after `set -euo
# pipefail`, before defining launcher-specific parameters) via:
#   source "$(dirname "$0")/_paper_scale_common.sh"
#
# This file only sets defaults / defines helpers; it does not read CLI args
# or run anything itself. Launcher-specific parameters, rationale comments,
# and the actual python invocation stay in each launcher.

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PYTHONPATH="${REPO_ROOT}:${REPO_ROOT}/src:${PYTHONPATH:-}"
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"
export KMP_DUPLICATE_LIB_OK="${KMP_DUPLICATE_LIB_OK:-TRUE}"

# Interpreter: override with PYTHON=/abs/path/to/python for a specific conda env.
PYTHON="${PYTHON:-python}"

# ─── Perf flags (decisions.md 追補 2026-07-22/23, runs/scaleup_a + scaleup_matrix) ─
# NO_EMPTY_CACHE=1 (default): skip per-step torch.cuda.empty_cache(). On WSL,
#   expandable_segments absorbs fragmentation (reserved flat, 1.29-1.41x faster)
#   — the measured recommended operation for WSL runs. Set NO_EMPTY_CACHE=0 only
#   on Windows-native python, where expandable_segments is a no-op.
# COMPILE=1 (opt-in, default 0): torch.compile — 1.24x alone, 1.65x combined
#   with NO_EMPTY_CACHE=1, and ~8% less VRAM. Forces match eager within the
#   TF32 noise floor. Kept opt-in until the first long soak on a production
#   workload (barostat/bond formation + resume) per decisions.md 2026-07-23.
NO_EMPTY_CACHE="${NO_EMPTY_CACHE:-1}"
COMPILE="${COMPILE:-0}"

# build_perf_flags: sets PERF_FLAGS from NO_EMPTY_CACHE / COMPILE.
build_perf_flags() {
    PERF_FLAGS=""
    if [ "${NO_EMPTY_CACHE}" = "1" ]; then PERF_FLAGS="--no-empty-cache"; fi
    if [ "${COMPILE}" = "1" ]; then PERF_FLAGS="${PERF_FLAGS} --compile"; fi
}

# build_resume_flag: sets RESUME_FLAG from RESUME (1 -> --resume). RESUME=1
# continues from ${OUTPUT_DIR}/checkpoint.pkl after a killed run; checkpoints
# are written every cycle by default.
build_resume_flag() {
    RESUME_FLAG=""
    if [ "${RESUME:-0}" = "1" ]; then RESUME_FLAG="--resume"; fi
}

# build_barostat_flag <default_ensemble>
#   default_ensemble=npt : NPT-by-default launchers. NO_BAROSTAT=1 -> --no-barostat.
#   default_ensemble=nvt : NVT-by-default launchers. BAROSTAT=1 -> --barostat
#                          --pressure "${PRESSURE:-1.0}".
build_barostat_flag() {
    local default_ensemble="$1"
    BAROSTAT_FLAG=""
    if [ "${default_ensemble}" = "npt" ]; then
        if [ "${NO_BAROSTAT:-0}" = "1" ]; then BAROSTAT_FLAG="--no-barostat"; fi
    else
        if [ "${BAROSTAT:-0}" = "1" ]; then BAROSTAT_FLAG="--barostat --pressure ${PRESSURE:-1.0}"; fi
    fi
}

# build_mix_flags: sets MIX_FLAGS from MIX / MIX_PS / MIX_SETTLE_STEPS /
# MIX_PLATFORM (specs/decisions.md 2026-08-04 — the --mix CLI is shared by the
# vinyl copolymer, nylon-6,6 and epoxy-amine drivers).
#
# MIX=1 enables the well-mixed measurement mode: one classical OpenMM/OpenFF
# mixing segment per cycle, which refreshes the reactive neighbourhood that
# step-growth/curing runs deplete. OFF by default — the paper-faithful loop is
# unchanged. Every knob is forwarded ONLY when set, so an unset knob keeps the
# CLI's documented default (mix_ps=25 ps from WM-P5b, mix_settle_steps=500,
# mix_platform=CPU; the WM campaigns ran mix_platform=CUDA).
#
# A knob set without MIX=1 aborts here rather than reaching the driver: under
# `set -euo pipefail` the CLI's own stray-knob error is equally fatal but its
# message points at a flag the user never typed.
build_mix_flags() {
    MIX_FLAGS=""
    if [ "${MIX:-0}" = "1" ]; then
        MIX_FLAGS="--mix"
        if [ -n "${MIX_PS:-}" ]; then MIX_FLAGS="${MIX_FLAGS} --mix-ps ${MIX_PS}"; fi
        if [ -n "${MIX_SETTLE_STEPS:-}" ]; then MIX_FLAGS="${MIX_FLAGS} --mix-settle-steps ${MIX_SETTLE_STEPS}"; fi
        if [ -n "${MIX_PLATFORM:-}" ]; then MIX_FLAGS="${MIX_FLAGS} --mix-platform ${MIX_PLATFORM}"; fi
        return 0
    fi
    local stray=""
    local knob
    for knob in MIX_PS MIX_SETTLE_STEPS MIX_PLATFORM; do
        if [ -n "${!knob:-}" ]; then stray="${stray:+${stray}, }${knob}"; fi
    done
    if [ -n "${stray}" ]; then
        echo "ERROR: ${stray} set without MIX=1 — mixing is off, so the knob(s) would do nothing." >&2
        echo "       Add MIX=1 to enable the classical mixing stage, or unset the knob(s)." >&2
        exit 1
    fi
}

# mix_status: one-line description of the mixing configuration for the launcher
# header ('off', or 'on' plus the flags actually forwarded).
mix_status() {
    if [ -z "${MIX_FLAGS:-}" ]; then
        echo "off"
    else
        echo "on (${MIX_FLAGS}; unset knobs use CLI defaults)"
    fi
}

# check_vram <warn_threshold_mb> <description>
# Prints the detected total GPU VRAM and a warning if it is below the given
# threshold. `description` is echoed alongside the reading and reused in the
# warning line, so it should carry both the "what this system wants" context
# and a shrink suggestion.
#
# Never aborts the script under `set -euo pipefail`: if nvidia-smi is absent
# this is a silent no-op; if nvidia-smi is present but the query fails or
# returns nothing, that is reported as a failed detection and the launcher
# continues.
check_vram() {
    local warn_threshold_mb="$1"
    local description="$2"
    if ! command -v nvidia-smi >/dev/null 2>&1; then
        return 0
    fi
    local vram_mb=""
    vram_mb="$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d ' ')" || true
    if [ -z "${vram_mb}" ]; then
        echo "  GPU VRAM: (detection failed; ${description})"
        echo ""
        return 0
    fi
    echo "  GPU VRAM: ${vram_mb} MB (${description})"
    if [ "${vram_mb}" -lt "${warn_threshold_mb}" ] 2>/dev/null; then
        echo "  WARNING: <${warn_threshold_mb} MB VRAM — ${description}"
    fi
    echo ""
}
