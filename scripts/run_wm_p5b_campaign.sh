#!/usr/bin/env bash
# WM-P5b production campaign — mix_time_ps sweep x multi-seed, deterministic
# selection only.
#
# Design ref: specs/decisions.md "追補 2026-07-19 — WM-P5b 実験計画
# (sweep×多シード、ユーザー承認スコープ C)". Answers two questions left open by
# WM-P4 item 4/6 (specs/decisions.md "追補 2026-07-19 — WM-P4 混合実装 フラット
# レビューと対応"):
#   Q1 - is the WM-P4 confirmed-yield drop under mixing (A1->A2: 5->2) signal
#        or single-seed noise?
#   Q2 - does a shorter mix_time_ps preserve the candidate-refresh benefit
#        (re-selection rate / Jaccard) while protecting yield?
#
# ─── Matrix (FIXED by decisions.md — do not change here) ───────────────────
#   Deterministic selection only (isolates the mixing effect; no softmax —
#   extends the WM-P4 A1/A2 axis). mix_time_ps in {0, 25, 50, 100} (0 = no-mix
#   baseline) x seed in {7, 11, 17, 23, 42} = 20 runs.
#   System/schedule identical to WM-P4 A2_mix50 (scripts/run_wm_p4_campaign.sh):
#   20 acrylate + 20 methacrylate + 2 initiators (564 atoms), 15 cycles,
#   biased 2000 / unbiased 1500 / equil 2000 steps, timestep 0.25 fs, 333 K,
#   orb/CUDA backend. mix_ps > 0 arms add classical mixing: 500 settle steps,
#   CUDA mixing platform, NAGL charges (--mix-ps/--mix-settle-steps/
#   --mix-platform/--mix-charge-method copied VERBATIM from the A2_mix50 arm
#   flags in run_wm_p4_campaign.sh). mix_ps == 0 runs WITHOUT any --mix* flag
#   at all (the paper-faithful no-mix baseline), not `--mix --mix-ps 0`.
#   Output dir per run: runs/wm_p5b/s{seed}_mix{ps}.
#
# ─── Retry mechanism (REUSED from scripts/run_wm_p4_resilient.sh) ──────────
# The per-run attempt/resume/stuck-detection/GPU-drain loop below is a
# faithful line-for-line port of run_wm_p4_resilient.sh's per-arm loop,
# generalized from "arm" to "run" (mix_ps/seed instead of arm name). It is
# duplicated rather than sourced because run_wm_p4_resilient.sh is a
# top-level script that executes its own 4-arm campaign as soon as it is
# read (no functions are factored out to `return` before that happens), so
# `source`-ing it here would immediately kick off the WM-P4 campaign too.
# See run_wm_p4_resilient.sh for the canonical version and its design notes;
# any behavioural change to the retry policy should be made in both places.
#
# Retry policy per run (identical to run_wm_p4_resilient.sh):
#   * Up to MAX_ATTEMPTS (default 8 -- higher than P4's 6: 20 runs over ~3.5
#     days give more opportunities to hit a transient CUDA fault).
#   * Attempt 1 runs fresh if no checkpoint exists, else --resume. Every later
#     attempt passes --resume (gated on checkpoint.pkl existence).
#   * After a crash, compare checkpoint next_cycle before (prev) vs. after
#     (now):
#       - now > prev  => made progress => transient => stuck counter reset, retry.
#       - now <= prev => no progress   => stuck++. Two consecutive no-progress
#         attempts => FAIL-STUCK, stop retrying this run.
#   * Between attempts: reap any hung process for THIS run, wait for the GPU
#     to drain (< 2000 MiB) or ~120 s, then a short settle sleep.
#   * MAX_ATTEMPTS reached with an advancing checkpoint => FAIL-EXHAUSTED.
#   * A failed/exhausted/stuck run does NOT abort the campaign — the next run
#     in the matrix still runs.
#
# ─── Idempotency (campaign-restart-safe, run-level) ────────────────────────
# Before starting a run, if <output-dir>/summary.json already exists, is
# non-empty, and parses as JSON, the run is SKIPPED (logged as SKIP-DONE).
# summary.json is written by run_vinyl_copolymer.py only once, at the very
# end of a fully successful run (see scripts/run_vinyl_copolymer.py, the
# `out_path.write_text(...)` near the bottom of `main()`), so its presence is
# a reliable "this run already completed" signal — a crash never leaves a
# partial summary.json behind. Restarting this launcher after a kill/reboot
# therefore resumes at the run level (already-finished runs are skipped
# outright) as well as the per-run cycle level (an in-progress run resumes
# from its checkpoint via the retry loop above).
#
# ─── RUN IT (detached, on the GPU box / WSL pfpoly-gpu) ─────────────────────
#   nohup bash scripts/run_wm_p5b_campaign.sh > runs/wm_p5b/campaign.out 2>&1 &
#
#   # Only a subset of the matrix / fewer attempts:
#   SEEDS="7 11" MIX_PS="0 50" MAX_ATTEMPTS=4 bash scripts/run_wm_p5b_campaign.sh
#
# ─── Estimated wall-clock (specs/decisions.md, measured unit costs) ────────
#   mix_ps=0 ~2.8h, 25 ~3.4h, 50 ~3.7h, 100 ~4.2h per seed (GPU sequential).
#   1 seed, full mix_ps sweep: ~14.1h. x5 seeds: ~70h nominal; with retry
#   headroom for the intermittent WSL/CUDA instability, budget ~3.5-4 days.
#
# ─── Testing hooks (inert in production; all env vars unset by default) ────
#   WM_P5B_STUB_CMD    : if set, each attempt runs `eval "$WM_P5B_STUB_CMD"`
#                        instead of the real python command. Lets the retry +
#                        matrix + idempotency logic be tested with no
#                        orb/GPU/python. The stub's exit code drives the loop
#                        exactly like the real command's would. The stub can
#                        inspect SEED/MIX_PS/OUT/CKPT/RESUME_FLAG/EXTRA_ARGS
#                        (exported) to assert on the exact flags this launcher
#                        would have passed for a given (seed, mix_ps).
#   WM_P4_STUB_CMD     : accepted as a fallback alias for WM_P5B_STUB_CMD (for
#                        parity with run_wm_p4_resilient.sh); WM_P5B_STUB_CMD
#                        takes precedence if both are set.
#   WM_P5B_GPU_WAIT    : max seconds to wait for the GPU to drain (default 120).
#   WM_P5B_POLL_SEC    : GPU poll interval (default 5).
#   WM_P5B_SETTLE_SEC  : settle sleep after the GPU drains (default 10).
#   WM_P5B_CAMPAIGN_DIR: output/log root (default runs/wm_p5b). Overridden by
#                        the test harness so it never writes under the live run.
#   GPU memory threshold for "drained" is < 2000 MiB.
# nvidia-smi / the checkpoint python one-liner can be shadowed on PATH for tests.
#
#   SEEDS / MIX_PS     : override the matrix (space-separated), e.g. for a
#                        partial re-run. Order in the env var is preserved.
#   MAX_ATTEMPTS       : override the per-run retry cap (default 8).

set -uo pipefail   # NOT -e: a failing run must not abort the whole campaign.

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "${REPO_ROOT}"

export PYTHONPATH="${REPO_ROOT}:${REPO_ROOT}/src:${PYTHONPATH:-}"
export KMP_DUPLICATE_LIB_OK="${KMP_DUPLICATE_LIB_OK:-TRUE}"
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"

# Interpreter: the MLIP conda env. Override with PYTHON=/abs/path/to/python.
PYTHON="${PYTHON:-${HOME}/miniconda3/envs/pfpoly-gpu/bin/python}"

CAMPAIGN_DIR="${WM_P5B_CAMPAIGN_DIR:-runs/wm_p5b}"
mkdir -p "${CAMPAIGN_DIR}"
CAMPAIGN_LOG="${CAMPAIGN_DIR}/campaign.log"

MAX_ATTEMPTS="${MAX_ATTEMPTS:-8}"

# GPU-recovery tunables (overridable so tests run fast).
GPU_WAIT="${WM_P5B_GPU_WAIT:-120}"
POLL_SEC="${WM_P5B_POLL_SEC:-5}"
SETTLE_SEC="${WM_P5B_SETTLE_SEC:-10}"
GPU_FREE_MIB=2000

# Stub hook: WM_P5B_STUB_CMD takes precedence; WM_P4_STUB_CMD accepted as an
# alias so the same env-var name works with either launcher.
STUB_CMD="${WM_P5B_STUB_CMD:-${WM_P4_STUB_CMD:-}}"

# Shared composition / schedule (identical across every run) — copied
# VERBATIM from the A2_mix50 arm in scripts/run_wm_p4_campaign.sh, except
# --seed (which varies per run and is appended below).
COMMON=(
    --n-acrylate 20
    --n-methacrylate 20
    --n-initiators 2
    --n-cycles 15
    --biased-steps 2000
    --unbiased-steps 1500
    --equil-steps 2000
    --timestep-fs 0.25
    --temperature 333
    --backend orb
    --device cuda
)
# Deterministic selection is the CLI default (no --selection-policy flag) --
# WM-P5b fixes the deterministic arm on purpose (decisions.md: isolate the
# mixing effect, no softmax stochasticity).

# Mixing flags for mix_ps > 0 — copied VERBATIM from the A2_mix50 arm in
# scripts/run_wm_p4_campaign.sh, with --mix-ps substituted per run.
mix_flags_for() {
    local ps="$1"
    if [ "${ps}" = "0" ]; then
        return 0   # no --mix* flags at all: the no-mix baseline
    fi
    printf '%s\n' --mix --mix-ps "${ps}" --mix-settle-steps 500 \
        --mix-platform CUDA --mix-charge-method nagl
}

# Read checkpoint next_cycle for a run; prints 0 if absent/corrupt. Never fails.
# (Ported from run_wm_p4_resilient.sh's read_next_cycle.)
read_next_cycle() {
    local ckpt="$1"
    if [ ! -f "${ckpt}" ]; then
        echo 0
        return 0
    fi
    local val
    val="$("${PYTHON}" -c "import pickle
try:
    print(pickle.load(open('${ckpt}','rb')).get('next_cycle', 0))
except Exception:
    print(0)" 2>/dev/null)"
    case "${val}" in
        ''|*[!0-9]*) echo 0 ;;
        *) echo "${val}" ;;
    esac
}

# Whether <out>/summary.json exists, is non-empty, and parses as JSON — the
# run-level idempotency check (a crashed run never leaves a valid
# summary.json behind; see the header comment).
has_valid_summary() {
    local path="$1"
    if [ ! -s "${path}" ]; then
        return 1
    fi
    "${PYTHON}" -c "import json,sys
try:
    json.load(open(sys.argv[1], encoding='utf-8'))
except Exception:
    sys.exit(1)" "${path}" 2>/dev/null
}

# Wait for the GPU to drain the dead process's CUDA context (or GPU_WAIT secs).
# (Ported from run_wm_p4_resilient.sh's wait_for_gpu.)
wait_for_gpu() {
    local waited=0 used
    while [ "${waited}" -lt "${GPU_WAIT}" ]; do
        used="$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null | head -n1 | tr -d '[:space:]')"
        case "${used}" in
            ''|*[!0-9]*) used=999999 ;;   # unreadable -> keep waiting up to cap
        esac
        if [ "${used}" -lt "${GPU_FREE_MIB}" ]; then
            return 0
        fi
        sleep "${POLL_SEC}"
        waited=$((waited + POLL_SEC))
    done
    return 0
}

# Matrix (overridable for partial/debug runs). Order is preserved: mix_ps is
# the outer loop, seed the inner loop, so the default enumeration is
# s7_mix0, s11_mix0, ..., s42_mix0, s7_mix25, ..., s42_mix100 (20 runs).
MIX_PS="${MIX_PS:-0 25 50 100}"
SEEDS="${SEEDS:-7 11 17 23 42}"

echo "=== WM-P5b campaign start $(date '+%Y-%m-%d %H:%M:%S') ===" | tee -a "${CAMPAIGN_LOG}"
echo "  Python:       ${PYTHON}"       | tee -a "${CAMPAIGN_LOG}"
echo "  mix_ps:       ${MIX_PS}"       | tee -a "${CAMPAIGN_LOG}"
echo "  seeds:        ${SEEDS}"        | tee -a "${CAMPAIGN_LOG}"
echo "  MAX_ATTEMPTS: ${MAX_ATTEMPTS}" | tee -a "${CAMPAIGN_LOG}"
[ -n "${STUB_CMD}" ] && echo "  STUB MODE:    ${STUB_CMD}" | tee -a "${CAMPAIGN_LOG}"

for ps in ${MIX_PS}; do
    mapfile -t EXTRA < <(mix_flags_for "${ps}")
    for seed in ${SEEDS}; do
        RUN="s${seed}_mix${ps}"
        OUT="${CAMPAIGN_DIR}/${RUN}"
        RUN_LOG="${CAMPAIGN_DIR}/${RUN}.log"
        CKPT="${OUT}/checkpoint.pkl"
        SUMMARY="${OUT}/summary.json"
        mkdir -p "${OUT}"

        # Run-level idempotency: a campaign restart resumes past runs that
        # already finished successfully.
        if has_valid_summary "${SUMMARY}"; then
            echo "SKIP-DONE ${RUN} (valid summary.json already present) $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "${CAMPAIGN_LOG}"
            continue
        fi

        echo "START ${RUN} $(date '+%Y-%m-%d %H:%M:%S')  -> ${OUT}" | tee -a "${CAMPAIGN_LOG}"

        attempts=0
        stuck=0
        run_start_cycle="$(read_next_cycle "${CKPT}")"
        while [ "${attempts}" -lt "${MAX_ATTEMPTS}" ]; do
            prev="$(read_next_cycle "${CKPT}")"

            # Resume iff a checkpoint already exists (crash recovery or
            # pre-existing).
            RESUME_FLAG=""
            if [ -f "${CKPT}" ]; then RESUME_FLAG="--resume"; fi

            n=$((attempts + 1))
            echo "RUN   ${RUN} attempt=${n}/${MAX_ATTEMPTS} resume=$( [ -n "${RESUME_FLAG}" ] && echo yes || echo no ) from cycle=${prev} $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "${CAMPAIGN_LOG}"

            SECONDS=0
            if [ -n "${STUB_CMD}" ]; then
                # Test hook: run the stub instead of the real command. Inert
                # in production (WM_P5B_STUB_CMD/WM_P4_STUB_CMD unset).
                # SEED/MIX_PS/OUT/CKPT/RESUME_FLAG/EXTRA_ARGS are exported so
                # the stub can inspect and assert on them.
                SEED="${seed}" MIX_PS_VAL="${ps}" OUT="${OUT}" CKPT="${CKPT}" \
                    RESUME_FLAG="${RESUME_FLAG}" EXTRA_ARGS="${EXTRA[*]:-}" \
                    eval "${STUB_CMD}" 2>&1 | tee -a "${RUN_LOG}"
                RC=${PIPESTATUS[0]}
            else
                "${PYTHON}" scripts/run_vinyl_copolymer.py \
                    --output-dir "${OUT}" \
                    "${COMMON[@]}" \
                    --seed "${seed}" \
                    "${EXTRA[@]}" \
                    ${RESUME_FLAG} 2>&1 | tee -a "${RUN_LOG}"
                RC=${PIPESTATUS[0]}
            fi
            DUR=${SECONDS}

            if [ "${RC}" -eq 0 ]; then
                echo "END   ${RUN} rc=0 dur=${DUR}s cycle=${run_start_cycle}->$(read_next_cycle "${CKPT}") $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "${CAMPAIGN_LOG}"
                break
            fi

            now="$(read_next_cycle "${CKPT}")"
            echo "RETRY ${RUN} attempt=${n} rc=${RC} cycle ${prev}->${now} dur=${DUR}s $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "${CAMPAIGN_LOG}"

            # Progress check: did we complete at least one new cycle this attempt?
            if [ "${now}" -gt "${prev}" ]; then
                stuck=0
            else
                stuck=$((stuck + 1))
                if [ "${stuck}" -ge 2 ]; then
                    echo "FAIL-STUCK ${RUN} (no cycle progress across 2 attempts — persistent error, not transient) $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "${CAMPAIGN_LOG}"
                    break
                fi
            fi

            # GPU recovery: reap any hung process for THIS run, then wait for
            # the GPU to release the dead CUDA context before the next attempt.
            pkill -f "output-dir ${OUT}" 2>/dev/null || true
            wait_for_gpu
            sleep "${SETTLE_SEC}"

            attempts=$((attempts + 1))
        done

        if [ "${attempts}" -ge "${MAX_ATTEMPTS}" ]; then
            echo "FAIL-EXHAUSTED ${RUN} after ${MAX_ATTEMPTS} attempts $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "${CAMPAIGN_LOG}"
        fi
    done
done

echo "=== WM-P5b campaign done $(date '+%Y-%m-%d %H:%M:%S') ===" | tee -a "${CAMPAIGN_LOG}"
