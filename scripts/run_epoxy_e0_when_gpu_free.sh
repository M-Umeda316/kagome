#!/usr/bin/env bash
# Track 2 / E0 — run the epoxy-amine ring-opening de-risk WHEN the local GPU is
# free. Intended to be launched later (e.g. ~6 h from now) by a scheduler; on
# firing it polls nvidia-smi and runs scripts/scan_epoxy_amine.py only if the GPU
# is idle enough, otherwise retries a few times then gives up (logged).
#
# The E0 scan is short (a few small-molecule FIRE relaxations + a ~24-point PES
# scan), so this needs the GPU free for only a few minutes.
#
# Everything is logged to ${LOG}. Tunable via env:
#   FREE_MEM_MIN_MB (default 6000)  minimum free VRAM to call the GPU "free"
#   UTIL_MAX        (default 30)    maximum GPU utilization %% to call it "free"
#   TRIES           (default 6)     number of checks
#   INTERVAL_S      (default 600)   seconds between checks (6*600 = 1 h window)
#   PYTHON          (default the pfpoly-gpu env python)  MLIP interpreter
#   EPOXY_SMILES    (default unset -> script default: propylene oxide CC1CO1)

set -uo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "${REPO_ROOT}"
export PYTHONPATH="${REPO_ROOT}:${REPO_ROOT}/src:${PYTHONPATH:-}"
export KMP_DUPLICATE_LIB_OK="${KMP_DUPLICATE_LIB_OK:-TRUE}"
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"

PYTHON="${PYTHON:-/home/shanu/miniconda3/envs/pfpoly-gpu/bin/python}"
FREE_MEM_MIN_MB="${FREE_MEM_MIN_MB:-6000}"
UTIL_MAX="${UTIL_MAX:-30}"
TRIES="${TRIES:-6}"
INTERVAL_S="${INTERVAL_S:-600}"
OUTPUT_DIR="${OUTPUT_DIR:-runs/calibration}"
LOG="${LOG:-${OUTPUT_DIR}/epoxy_e0_scheduled.log}"

mkdir -p "${OUTPUT_DIR}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG}"; }

log "=== E0 scheduled runner started (pid $$) ==="
log "thresholds: free_mem>=${FREE_MEM_MIN_MB}MB, util<=${UTIL_MAX}%, tries=${TRIES}, interval=${INTERVAL_S}s"

if ! command -v nvidia-smi >/dev/null 2>&1; then
    log "ERROR: nvidia-smi not found — cannot check GPU. Aborting."
    exit 1
fi

gpu_is_free() {
    # Reads GPU 0: free VRAM (MiB) and utilization (%). Returns 0 (free) / 1 (busy).
    local line free util
    line="$(nvidia-smi --query-gpu=memory.free,utilization.gpu \
            --format=csv,noheader,nounits 2>/dev/null | head -1)"
    free="$(echo "${line}" | awk -F',' '{gsub(/ /,"",$1); print $1}')"
    util="$(echo "${line}" | awk -F',' '{gsub(/ /,"",$2); print $2}')"
    if [ -z "${free:-}" ] || [ -z "${util:-}" ]; then
        log "  WARN: could not parse nvidia-smi output: '${line}'"
        return 1
    fi
    log "  GPU: free=${free}MB util=${util}%"
    if [ "${free}" -ge "${FREE_MEM_MIN_MB}" ] && [ "${util}" -le "${UTIL_MAX}" ]; then
        return 0
    fi
    return 1
}

run_scan() {
    log "GPU is free -> running E0 scan (scan_epoxy_amine.py --backend orb --device cuda)"
    "${PYTHON}" scripts/scan_epoxy_amine.py \
        --backend orb --device cuda \
        ${EPOXY_SMILES:+--epoxy-smiles "${EPOXY_SMILES}"} \
        --output-dir "${OUTPUT_DIR}" >> "${LOG}" 2>&1
    local rc=$?
    log "E0 scan finished with exit=${rc}. Result: ${OUTPUT_DIR}/epoxy_amine_scan.json"
    return ${rc}
}

i=1
while [ "${i}" -le "${TRIES}" ]; do
    log "check ${i}/${TRIES}"
    if gpu_is_free; then
        run_scan
        exit $?
    fi
    if [ "${i}" -lt "${TRIES}" ]; then
        log "  GPU busy — waiting ${INTERVAL_S}s before next check."
        sleep "${INTERVAL_S}"
    fi
    i=$((i + 1))
done

log "GPU stayed busy across ${TRIES} checks — E0 scan NOT run. Re-launch when free, "
log "or run manually: ${PYTHON} scripts/scan_epoxy_amine.py --backend orb --device cuda --output-dir ${OUTPUT_DIR}"
exit 2
