#!/usr/bin/env bash
# WM-P4 resilient auto-resume launcher — well-mixed measurement campaign.
#
# Same four arms / same composition / same schedule as scripts/run_wm_p4_campaign.sh
# (arm flags copied VERBATIM from that script), but wrapped in a per-arm retry
# loop that auto-resumes an arm past *transient* CUDA crashes.
#
# Why this exists: the orb/CUDA stack is intermittently unstable in this WSL box.
# Two arms have died on transient errors that are NOT code bugs:
#   - CUBLAS_STATUS_INTERNAL_ERROR at startup
#   - torch.AcceleratorError: CUDA error: unknown error (mid-run)
# A fresh process succeeds. Each arm checkpoints every cycle to
# <output-dir>/checkpoint.pkl and supports bit-exact --resume (rng restored), so
# a crashed arm can pick up from its last completed cycle. This launcher retries
# a crashed arm, resuming from its checkpoint, until it either finishes or proves
# the error is persistent (no cycle progress).
#
# ─── RUN IT (detached, on the GPU box / WSL pfpoly-gpu) ─────────────────────
#   nohup bash scripts/run_wm_p4_resilient.sh > runs/wm_p4/campaign.out 2>&1 &
#
#   # Only a subset / fewer attempts:
#   ARMS="A2_mix50 A4_both" MAX_ATTEMPTS=4 bash scripts/run_wm_p4_resilient.sh
#
# ─── Retry policy per arm ───────────────────────────────────────────────────
#   * Up to MAX_ATTEMPTS (default 6) attempts.
#   * Attempt 1 runs fresh if no checkpoint exists, else --resume. Every later
#     attempt passes --resume (gated on checkpoint.pkl existence).
#   * After a crash we compare checkpoint next_cycle before (prev) vs. after
#     (now):
#       - now > prev  => made progress => transient => stuck counter reset, retry.
#       - now <= prev => no progress   => stuck++. Two consecutive no-progress
#         attempts => FAIL-STUCK, stop retrying this arm (persistent error, not
#         worth burning attempts on).
#   * Between attempts: reap any hung process for THIS arm, wait for the GPU to
#     drain (CUDA context release) or ~120 s, then a short settle sleep.
#   * MAX_ATTEMPTS reached with an advancing checkpoint => FAIL-EXHAUSTED.
#   * A failed/exhausted/stuck arm does NOT abort the campaign — next arm runs.
#
# ─── Testing hooks (inert in production; all env vars unset by default) ──────
#   WM_P4_STUB_CMD  : if set, each attempt runs `eval "$WM_P4_STUB_CMD"` instead
#                     of the real python command. Lets the retry logic be tested
#                     with no orb/GPU/python. The stub's exit code drives the loop
#                     exactly like the real command's would.
#   WM_P4_GPU_WAIT  : max seconds to wait for the GPU to drain (default 120).
#   WM_P4_POLL_SEC  : GPU poll interval (default 5).
#   WM_P4_SETTLE_SEC: settle sleep after the GPU drains (default 10).
#   WM_P4_CAMPAIGN_DIR: output/log root (default runs/wm_p4). Overridden only by
#                     the test harness so it never writes under the live run.
#   GPU memory threshold for "drained" is < 2000 MiB.
# nvidia-smi / the checkpoint python one-liner can be shadowed on PATH for tests.

set -uo pipefail   # NOT -e: a failing arm must not abort the whole campaign.

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "${REPO_ROOT}"

export PYTHONPATH="${REPO_ROOT}:${REPO_ROOT}/src:${PYTHONPATH:-}"
export KMP_DUPLICATE_LIB_OK="${KMP_DUPLICATE_LIB_OK:-TRUE}"
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"

# Interpreter: the MLIP conda env. Override with PYTHON=/abs/path/to/python.
PYTHON="${PYTHON:-${HOME}/miniconda3/envs/pfpoly-gpu/bin/python}"

CAMPAIGN_DIR="${WM_P4_CAMPAIGN_DIR:-runs/wm_p4}"
mkdir -p "${CAMPAIGN_DIR}"
CAMPAIGN_LOG="${CAMPAIGN_DIR}/campaign.log"

MAX_ATTEMPTS="${MAX_ATTEMPTS:-6}"

# GPU-recovery tunables (overridable so tests run fast).
GPU_WAIT="${WM_P4_GPU_WAIT:-120}"
POLL_SEC="${WM_P4_POLL_SEC:-5}"
SETTLE_SEC="${WM_P4_SETTLE_SEC:-10}"
GPU_FREE_MIB=2000

# Shared composition / schedule (identical across arms) — copied from
# scripts/run_wm_p4_campaign.sh.
COMMON=(
    --n-acrylate 20
    --n-methacrylate 20
    --n-initiators 2
    --seed 7
    --n-cycles 15
    --biased-steps 2000
    --unbiased-steps 1500
    --equil-steps 2000
    --timestep-fs 0.25
    --temperature 333
    --backend orb
    --device cuda
)

# Per-arm distinguishing flags — copied VERBATIM from run_wm_p4_campaign.sh.
MIX_FLAGS=(--mix --mix-ps 50 --mix-settle-steps 500 --mix-platform CUDA --mix-charge-method nagl)
SOFTMAX_FLAGS=(--selection-policy softmax --selection-temperature 0.5)

# Assemble the CLI flags for a given arm name.
arm_flags() {
    local arm="$1"
    case "${arm}" in
        A1_baseline) : ;;                                   # no extra flags
        A2_mix50)    printf '%s\n' "${MIX_FLAGS[@]}" ;;
        A3_softmax)  printf '%s\n' "${SOFTMAX_FLAGS[@]}" ;;
        A4_both)     printf '%s\n' "${MIX_FLAGS[@]}" "${SOFTMAX_FLAGS[@]}" ;;
        *) echo "__UNKNOWN_ARM__" ;;
    esac
}

# Read checkpoint next_cycle for an arm; prints 0 if absent/corrupt. Never fails.
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
    # Guard against empty / non-numeric output.
    case "${val}" in
        ''|*[!0-9]*) echo 0 ;;
        *) echo "${val}" ;;
    esac
}

# Wait for the GPU to drain the dead process's CUDA context (or GPU_WAIT secs).
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

ARMS="${ARMS:-A1_baseline A2_mix50 A3_softmax A4_both}"

echo "=== WM-P4 resilient campaign start $(date '+%Y-%m-%d %H:%M:%S') ===" | tee -a "${CAMPAIGN_LOG}"
echo "  Python:       ${PYTHON}"       | tee -a "${CAMPAIGN_LOG}"
echo "  Arms:         ${ARMS}"         | tee -a "${CAMPAIGN_LOG}"
echo "  MAX_ATTEMPTS: ${MAX_ATTEMPTS}" | tee -a "${CAMPAIGN_LOG}"
[ -n "${WM_P4_STUB_CMD:-}" ] && echo "  STUB MODE:    ${WM_P4_STUB_CMD}" | tee -a "${CAMPAIGN_LOG}"

for arm in ${ARMS}; do
    mapfile -t EXTRA < <(arm_flags "${arm}")
    if [ "${EXTRA[0]:-}" = "__UNKNOWN_ARM__" ]; then
        echo "SKIP  ${arm}: unknown arm name $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "${CAMPAIGN_LOG}"
        continue
    fi
    OUT="${CAMPAIGN_DIR}/${arm}"
    ARM_LOG="${CAMPAIGN_DIR}/${arm}.log"
    CKPT="${OUT}/checkpoint.pkl"
    mkdir -p "${OUT}"

    echo "START ${arm} $(date '+%Y-%m-%d %H:%M:%S')  -> ${OUT}" | tee -a "${CAMPAIGN_LOG}"

    attempts=0
    stuck=0
    while [ "${attempts}" -lt "${MAX_ATTEMPTS}" ]; do
        prev="$(read_next_cycle "${CKPT}")"

        # Resume iff a checkpoint already exists (crash recovery or pre-existing).
        RESUME_FLAG=""
        if [ -f "${CKPT}" ]; then RESUME_FLAG="--resume"; fi

        n=$((attempts + 1))
        echo "RUN   ${arm} attempt=${n}/${MAX_ATTEMPTS} resume=$( [ -n "${RESUME_FLAG}" ] && echo yes || echo no ) from cycle=${prev} $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "${CAMPAIGN_LOG}"

        SECONDS=0
        if [ -n "${WM_P4_STUB_CMD:-}" ]; then
            # Test hook: run the stub instead of the real command. Inert in
            # production (WM_P4_STUB_CMD unset). ARM/OUT/CKPT/RESUME_FLAG are
            # exported so the stub can inspect them.
            ARM="${arm}" OUT="${OUT}" CKPT="${CKPT}" RESUME_FLAG="${RESUME_FLAG}" \
                eval "${WM_P4_STUB_CMD}" 2>&1 | tee -a "${ARM_LOG}"
            RC=${PIPESTATUS[0]}
        else
            "${PYTHON}" scripts/run_vinyl_copolymer.py \
                --output-dir "${OUT}" \
                "${COMMON[@]}" \
                "${EXTRA[@]}" \
                ${RESUME_FLAG} 2>&1 | tee -a "${ARM_LOG}"
            RC=${PIPESTATUS[0]}
        fi
        DUR=${SECONDS}

        if [ "${RC}" -eq 0 ]; then
            echo "END   ${arm} rc=0 dur=${DUR}s cycle=${prev}->$(read_next_cycle "${CKPT}") $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "${CAMPAIGN_LOG}"
            break
        fi

        now="$(read_next_cycle "${CKPT}")"
        echo "RETRY ${arm} attempt=${n} rc=${RC} cycle ${prev}->${now} dur=${DUR}s $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "${CAMPAIGN_LOG}"

        # Progress check: did we complete at least one new cycle this attempt?
        if [ "${now}" -gt "${prev}" ]; then
            stuck=0
        else
            stuck=$((stuck + 1))
            if [ "${stuck}" -ge 2 ]; then
                echo "FAIL-STUCK ${arm} (no cycle progress across 2 attempts — persistent error, not transient) $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "${CAMPAIGN_LOG}"
                break
            fi
        fi

        # GPU recovery: reap any hung process for THIS arm, then wait for the
        # GPU to release the dead CUDA context before the next attempt.
        pkill -f "output-dir ${OUT}" 2>/dev/null || true
        wait_for_gpu
        sleep "${SETTLE_SEC}"

        attempts=$((attempts + 1))
    done

    if [ "${attempts}" -ge "${MAX_ATTEMPTS}" ]; then
        echo "FAIL-EXHAUSTED ${arm} after ${MAX_ATTEMPTS} attempts $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "${CAMPAIGN_LOG}"
    fi
done

echo "=== WM-P4 resilient campaign done $(date '+%Y-%m-%d %H:%M:%S') ===" | tee -a "${CAMPAIGN_LOG}"
