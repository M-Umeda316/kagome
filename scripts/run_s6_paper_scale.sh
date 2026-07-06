#!/usr/bin/env bash
# S6: Paper-scale vinyl radical polymerization reproduction
# Paper: arXiv:2511.22874, Table S1, SI S-3
# System: 200 methyl acrylate monomers + 10 AIBN initiators (≈2520 atoms)
# Required GPU VRAM: ≥24 GB (OrbMol-v2 single-step footprint ~9.5 GB at this scale)
#
# Usage:
#   # From repo root, with conda env active:
#   bash scripts/run_s6_paper_scale.sh
#
#   # Or with explicit seed/output-dir:
#   SEED=42 OUTPUT_DIR=runs/s6_seed42 bash scripts/run_s6_paper_scale.sh
#
# For 16 GB GPU (half-scale 100+5), run run_vinyl_aibn.py directly. This mirrors
# the full-scale flags below; only system size and n-cycles differ (RF22). This
# is the validated recipe (runs/s6_half_100x5: conversion 26%, T~338 K, swap-free):
#   python scripts/run_vinyl_aibn.py --seed 7 --output-dir runs/s6_half_scale \
#       --n-monomers 100 --n-initiators 5 --activation --activation-f2 0.3 \
#       --activation-f1-max 250 --activation-steps 5000 --f2 2 --friction-per-fs 0.01 \
#       --density 0.5 --temperature 333.0 --no-barostat --backend orb --device cuda \
#       --n-cycles 30 --biased-steps 2000 --unbiased-steps 1500 --equil-steps 2000 \
#       --timestep-fs 0.25 --minimize --minimize-fmax 1.0
#
# Environment: pfpoly-gpu (or equivalent clone; see docs below)
# Estimated wall-clock: 12-48 h depending on GPU and n_cycles
#
# --- Parameter rationale ---
# n_monomers=200, n_initiators=10  : Paper Table S1, Section 3 (200+10 system)
# activation                        : AIBN V^d C-N dissociation before propagation
# activation-f2=0.3, f1_max=250    : OrbMol-v2 C-N barrier ~39 kcal/mol requires these
#                                     (f2=10/f1_max=125 paper defaults insufficient for OrbMol-v2)
#                                     See specs/decisions.md 2026-06-18
# f2=2.0                           : Capture width. Paper 10 / repo prior 5 leave a dead
#                                     zone between the [3,6] candidate window and the ~2.5 A
#                                     bias-capture shell, so selected pairs feel ~0 force and
#                                     0 formations result in a melt. f2=2 (reach ~0.71 A)
#                                     bridges it. See specs/decisions.md 2026-06-26.
# friction_per_fs=0.01             : Langevin friction. Lowering f2 injects bias work that,
#                                     with addition heat, accumulates over cycles; 0.01 (vs
#                                     default 0.001) dissipates it and pins T near 333 K.
# n_cycles=50                      : Paper reports multi-hundred cycles; 50 is a feasible start
# biased_steps=2000                : Validated in S2-S3 runs
# unbiased_steps=1500              : Relaxation window; with friction 0.01 keeps T at 333 K
#                                     (decisions.md 2026-06-26). Paper uses 2000.
# timestep_fs=0.25                 : Paper value. REQUIRED for reactive multi-radical
#                                    stability — 1.0 fs numerically explodes the open-shell
#                                    melt (1e6-1e10 K). See specs/decisions.md 2026-06-25
#                                    CORRECTION and specs/validity-domain.md §2.1/§3.
# density=0.5                      : Paper SI S-3 (methyl acrylate melt)
# temperature=333.0                : Paper Table S1 (60°C)
# no-barostat                      : NVT validated (barostat is unstable for open-shell systems)
# seed=7                           : Matches S2-S3 best-result seed

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"
export KMP_DUPLICATE_LIB_OK="${KMP_DUPLICATE_LIB_OK:-TRUE}"

# Source calibrated parameters if available (from calibrate_tdbb.py).
# These set ACTIVATION_F2, ACTIVATION_F1_MAX, ACTIVATION_STEPS, F2,
# F1_MAX_FORMATION, F1_MAX_DISSOCIATION. Manual overrides below take
# precedence (env vars are set only if not already defined).
CALIB_ENV="${REPO_ROOT}/runs/calibration/tdbb_params.env"
if [ -f "${CALIB_ENV}" ]; then
    echo "Loading calibrated parameters from ${CALIB_ENV}"
    # shellcheck source=/dev/null
    source "${CALIB_ENV}"
fi

SEED="${SEED:-7}"
OUTPUT_DIR="${OUTPUT_DIR:-runs/s6_paper_scale_seed${SEED}}"
DEVICE="${DEVICE:-cuda}"
N_CYCLES="${N_CYCLES:-50}"
BIASED_STEPS="${BIASED_STEPS:-2000}"
UNBIASED_STEPS="${UNBIASED_STEPS:-1500}"
EQUIL_STEPS="${EQUIL_STEPS:-2000}"
F2="${F2:-2}"
F1_MAX_FORMATION="${F1_MAX_FORMATION:-250.0}"
F1_MAX_DISSOCIATION="${F1_MAX_DISSOCIATION:-125.0}"
FRICTION_PER_FS="${FRICTION_PER_FS:-0.01}"
ACTIVATION_F2="${ACTIVATION_F2:-0.3}"
ACTIVATION_F1_MAX="${ACTIVATION_F1_MAX:-250.0}"
ACTIVATION_STEPS="${ACTIVATION_STEPS:-5000}"
# RESUME=1 continues from ${OUTPUT_DIR}/checkpoint.pkl after a killed run
# (skips build-time activation, restarts at the saved cycle). Checkpoints are
# written every cycle by default.
RESUME="${RESUME:-0}"
RESUME_FLAG=""
if [ "${RESUME}" = "1" ]; then RESUME_FLAG="--resume"; fi

echo "=== S6 paper-scale run ==="
echo "  Seed:           ${SEED}"
echo "  Output dir:     ${OUTPUT_DIR}"
echo "  Device:         ${DEVICE}"
echo "  N cycles:       ${N_CYCLES}"
echo "  Biased steps:   ${BIASED_STEPS}"
echo "  Unbiased steps: ${UNBIASED_STEPS}"
echo "  Activation:     f2=${ACTIVATION_F2}, f1_max=${ACTIVATION_F1_MAX}, steps=${ACTIVATION_STEPS}"
echo "  Production:     f2=${F2}, f1_max_form=${F1_MAX_FORMATION}, f1_max_dissoc=${F1_MAX_DISSOCIATION}"
echo ""

# Warn if VRAM might be insufficient
if command -v nvidia-smi &>/dev/null; then
    VRAM_MB=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d ' ')
    VRAM_GB=$(echo "scale=1; ${VRAM_MB:-0} / 1024" | bc 2>/dev/null || echo "?")
    echo "  GPU VRAM:       ${VRAM_GB} GB (need ≥24 GB for 200+10 system)"
    if [ "${VRAM_MB:-0}" -lt 20000 ] 2>/dev/null; then
        echo "  WARNING: VRAM may be insufficient. Consider --n-monomers 100 --n-initiators 5 for 16 GB GPU."
    fi
fi

echo ""
echo "Starting run..."

python scripts/run_vinyl_aibn.py \
    --seed "${SEED}" \
    --output-dir "${OUTPUT_DIR}" \
    --n-monomers 200 \
    --n-initiators 10 \
    --activation \
    --activation-f2 "${ACTIVATION_F2}" \
    --activation-f1-max "${ACTIVATION_F1_MAX}" \
    --activation-steps "${ACTIVATION_STEPS}" \
    --f2 "${F2}" \
    --f1-max-formation "${F1_MAX_FORMATION}" \
    --f1-max-dissociation "${F1_MAX_DISSOCIATION}" \
    --friction-per-fs "${FRICTION_PER_FS}" \
    --density 0.5 \
    --temperature 333.0 \
    --no-barostat \
    --backend orb \
    --device "${DEVICE}" \
    --n-cycles "${N_CYCLES}" \
    --biased-steps "${BIASED_STEPS}" \
    --unbiased-steps "${UNBIASED_STEPS}" \
    --equil-steps "${EQUIL_STEPS}" \
    --timestep-fs 0.25 \
    --minimize \
    --minimize-fmax 1.0 \
    ${RESUME_FLAG}

echo ""
echo "Run complete. Generating figures..."

# n_reactive_sites is auto-read from the trajectory header (RF2: α denominator
# = n_monomers, written by PolymerizationWorkflow). No CLI override needed.
python scripts/reproduce_figures.py \
    --trajectory "${OUTPUT_DIR}/trajectory.jsonl" \
    --bonds "${OUTPUT_DIR}/bonds.jsonl" \
    --target-temperature 333.0 \
    --timestep-fs 0.25 \
    --output-dir "${OUTPUT_DIR}/figures"

echo ""
echo "=== S6 done ==="
echo "  Artifacts: ${OUTPUT_DIR}/"
echo "  Figures:   ${OUTPUT_DIR}/figures/"
echo ""
echo "Reproduction command:"
echo "  SEED=${SEED} OUTPUT_DIR=${OUTPUT_DIR} bash scripts/run_s6_paper_scale.sh"
