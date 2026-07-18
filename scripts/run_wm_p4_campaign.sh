#!/usr/bin/env bash
# WM-P4 validation campaign — well-mixed measurement mode vs. baseline.
#
# Four arms over the SAME 20 acrylate + 20 methacrylate + 2 initiator copolymer
# (~540 atoms), seed 7, 15 cycles, orb/CUDA, run SEQUENTIALLY (never two orb
# processes on the GPU at once):
#   A1_baseline  — paper-faithful: no mixing, deterministic selection
#   A2_mix50     — + classical mixing (50 ps, 500 settle, CUDA, NAGL charges)
#   A3_softmax   — + softmax candidate selection (T=0.5), no mixing
#   A4_both      — mixing + softmax together
#
# Design refs: specs/decisions.md "2026-07-17: well-mixed 測定モード" items
#   (i) mixing (WM-P3), (iv) softmax selection (WM-P5a). Branch carries the
#   mixing nonbonded-cutoff/box-clamp fix (fix/mixing-cutoff-box-clamp).
#
# ─── RUN IT (detached, on the GPU box / WSL pfpoly-gpu) ─────────────────────
#   nohup bash scripts/run_wm_p4_campaign.sh > runs/wm_p4/campaign.out 2>&1 &
#
#   # Resume every arm from its per-cycle checkpoint after a kill:
#   RESUME=1 nohup bash scripts/run_wm_p4_campaign.sh > runs/wm_p4/campaign.out 2>&1 &
#
#   # Run only a subset (space-separated, any order — flags keyed off the name):
#   ARMS="A2_mix50 A4_both" bash scripts/run_wm_p4_campaign.sh
#
# Each arm streams to runs/wm_p4/<arm>.log (via tee) and writes a per-cycle
# checkpoint.pkl in its own output dir, so any arm is independently resumable.
# One failing arm is recorded and the campaign continues to the next.

set -uo pipefail   # NOT -e: a failing arm must not abort the whole campaign.

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "${REPO_ROOT}"

export PYTHONPATH="${REPO_ROOT}:${REPO_ROOT}/src:${PYTHONPATH:-}"
export KMP_DUPLICATE_LIB_OK="${KMP_DUPLICATE_LIB_OK:-TRUE}"
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"

# Interpreter: the MLIP conda env. Override with PYTHON=/abs/path/to/python.
PYTHON="${PYTHON:-${HOME}/miniconda3/envs/pfpoly-gpu/bin/python}"

CAMPAIGN_DIR="runs/wm_p4"
mkdir -p "${CAMPAIGN_DIR}"
CAMPAIGN_LOG="${CAMPAIGN_DIR}/campaign.log"

# RESUME=1 -> add --resume to every arm.
RESUME_FLAG=""
if [ "${RESUME:-0}" = "1" ]; then RESUME_FLAG="--resume"; fi

# Shared composition / schedule (identical across arms).
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

# Per-arm distinguishing flags.
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

ARMS="${ARMS:-A1_baseline A2_mix50 A3_softmax A4_both}"

echo "=== WM-P4 campaign start $(date '+%Y-%m-%d %H:%M:%S') ===" | tee -a "${CAMPAIGN_LOG}"
echo "  Python:  ${PYTHON}" | tee -a "${CAMPAIGN_LOG}"
echo "  Arms:    ${ARMS}"    | tee -a "${CAMPAIGN_LOG}"
echo "  Resume:  $( [ -n "${RESUME_FLAG}" ] && echo yes || echo no )" | tee -a "${CAMPAIGN_LOG}"

for arm in ${ARMS}; do
    mapfile -t EXTRA < <(arm_flags "${arm}")
    if [ "${EXTRA[0]:-}" = "__UNKNOWN_ARM__" ]; then
        echo "SKIP  ${arm}: unknown arm name $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "${CAMPAIGN_LOG}"
        continue
    fi
    OUT="${CAMPAIGN_DIR}/${arm}"
    ARM_LOG="${CAMPAIGN_DIR}/${arm}.log"
    mkdir -p "${OUT}"

    echo "START ${arm} $(date '+%Y-%m-%d %H:%M:%S')  -> ${OUT}" | tee -a "${CAMPAIGN_LOG}"
    SECONDS=0
    "${PYTHON}" scripts/run_vinyl_copolymer.py \
        --output-dir "${OUT}" \
        "${COMMON[@]}" \
        "${EXTRA[@]}" \
        ${RESUME_FLAG} 2>&1 | tee "${ARM_LOG}"
    RC=${PIPESTATUS[0]}
    DUR=${SECONDS}

    if [ "${RC}" -eq 0 ]; then
        echo "END   ${arm} rc=0 dur=${DUR}s $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "${CAMPAIGN_LOG}"
    else
        echo "FAIL  ${arm} rc=${RC} dur=${DUR}s $(date '+%Y-%m-%d %H:%M:%S') (see ${ARM_LOG}) — continuing" | tee -a "${CAMPAIGN_LOG}"
    fi
done

echo "=== WM-P4 campaign done $(date '+%Y-%m-%d %H:%M:%S') ===" | tee -a "${CAMPAIGN_LOG}"
