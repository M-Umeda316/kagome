#!/usr/bin/env bash
# Nylon-6,6 step-growth polycondensation — PAPER-SCALE reproduction.
# Paper: arXiv:2511.22874, SI S-3, Fig. 4c (Carothers DPn vs conversion).
# System: 100 hexamethylenediamine + 100 adipic acid (equimolar) ≈ 4400 atoms.
#
# Purpose: produce many (p, DPn) points for the Carothers curve vs paper Fig. 4c.
# This is the scale-up of the validated small run (runs/nylon66_f2test: 10+10,
# 4 amide bonds, p=0.20, DPn=1.25). The reaction mechanics (paper-faithful
# conjunctive reaction event: amide C-N forms AND N-H/C-OH dissociate together;
# water is bias-only) are validated; only system size and n_cycles change here.
#
# ─── RUN IT (on a machine with a ≥24 GB GPU) ───────────────────────────────
#   # From the repo root, with the MLIP conda env active (python on PATH):
#   bash scripts/run_nylon66_paper_scale.sh
#
#   # Point at a specific interpreter / tune knobs without editing the file:
#   PYTHON=/path/to/envs/pfpoly-gpu/bin/python N_CYCLES=300 \
#       OUTPUT_DIR=runs/nylon_paper_seed7 bash scripts/run_nylon66_paper_scale.sh
#
#   # Resume a killed run from its last cycle checkpoint (checkpoints are written
#   # every cycle by default; resume skips build/minimize/equilibration):
#   RESUME=1 OUTPUT_DIR=runs/nylon_paper_seed7 bash scripts/run_nylon66_paper_scale.sh
#
# Estimated wall-clock: LONG (a Carothers curve needs many cycles). Each cycle is
# ~4000 MD steps over ~4400 atoms; budget on the order of minutes/cycle on a
# modern GPU, so N_CYCLES=200 is many hours. The run is fully resumable — run it
# in chunks with RESUME=1 and raise N_CYCLES as conversion climbs.
#
# ─── Parameter rationale ───────────────────────────────────────────────────
# n_diamines=100, n_diacids=100 : Paper SI S-3 (equimolar).
# density=0.5                   : Paper SI S-3 initial density (g/mL).
# NPT (barostat on), P=1 atm    : Paper ensemble. The MC barostat accounts for the
#                                 bias energy in its acceptance (decisions.md D2
#                                 2026-07-03). Set NO_BAROSTAT=1 to fall back to
#                                 NVT — that is what the small run validated, use
#                                 it if the barostat proves unstable at scale.
# temperature=300               : Paper SI S-3 (nylon-6,6).
# f2=2.0                        : OrbMol-v2 capture width. The paper default f2=10
#                                 (for PFP) leaves a dead-zone: the amide-forming
#                                 bias force is ~0 across the [3,6] Å candidate
#                                 window on OrbMol's PES, so 0 amides form. f2=2
#                                 bridges it (decisions.md 2026-07-08), mirroring
#                                 the vinyl/MA recipe.
# biased=2000, unbiased=2000    : Paper schedule (2000+2000 steps = 500 fs each).
# minimize + equil=2000         : Pre-TDBB relaxation of the compressed 0.5 g/mL
#                                 box; without it OrbMol segfaults on the first
#                                 biased step (decisions.md 2026-07-08). Paper p.20.
# checkpoint every cycle        : Resumable long run (RESUME=1).

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PYTHONPATH="${REPO_ROOT}:${REPO_ROOT}/src:${PYTHONPATH:-}"
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"
export KMP_DUPLICATE_LIB_OK="${KMP_DUPLICATE_LIB_OK:-TRUE}"

# Interpreter: override with PYTHON=/abs/path/to/python for a specific conda env.
PYTHON="${PYTHON:-python}"

SEED="${SEED:-7}"
OUTPUT_DIR="${OUTPUT_DIR:-runs/nylon66_paper_scale_seed${SEED}}"
DEVICE="${DEVICE:-cuda}"
N_DIAMINES="${N_DIAMINES:-100}"
N_DIACIDS="${N_DIACIDS:-100}"
N_CYCLES="${N_CYCLES:-200}"
BIASED_STEPS="${BIASED_STEPS:-2000}"
UNBIASED_STEPS="${UNBIASED_STEPS:-2000}"
EQUIL_STEPS="${EQUIL_STEPS:-2000}"
F2="${F2:-2}"
DENSITY="${DENSITY:-0.5}"
TEMPERATURE="${TEMPERATURE:-300.0}"
PRESSURE="${PRESSURE:-1.0}"

# NPT by default (paper). NO_BAROSTAT=1 -> NVT (the small-run-validated ensemble).
BAROSTAT_FLAG=""
if [ "${NO_BAROSTAT:-0}" = "1" ]; then BAROSTAT_FLAG="--no-barostat"; fi

# RESUME=1 continues from ${OUTPUT_DIR}/checkpoint.pkl after a killed run.
RESUME_FLAG=""
if [ "${RESUME:-0}" = "1" ]; then RESUME_FLAG="--resume"; fi

# Reactive-site count for the figure alpha(t) denominator: each diamine has 2
# amine_N ends, each diacid 2 carboxyl_C ends -> 2*(N_DIAMINES+N_DIACIDS).
N_REACTIVE_SITES=$(( 2 * (N_DIAMINES + N_DIACIDS) ))

echo "=== Nylon-6,6 paper-scale run ==="
echo "  Python:         ${PYTHON}"
echo "  Seed:           ${SEED}"
echo "  Output dir:     ${OUTPUT_DIR}"
echo "  Device:         ${DEVICE}"
echo "  System:         ${N_DIAMINES} diamine + ${N_DIACIDS} diacid (${N_REACTIVE_SITES} reactive ends)"
echo "  Ensemble:       $( [ -n "${BAROSTAT_FLAG}" ] && echo 'NVT (--no-barostat)' || echo 'NPT, '"${PRESSURE}"' atm' )"
echo "  Schedule:       ${N_CYCLES} cycles x (${BIASED_STEPS} biased + ${UNBIASED_STEPS} unbiased), f2=${F2}, T=${TEMPERATURE} K"
echo "  Resume:         $( [ -n "${RESUME_FLAG}" ] && echo yes || echo 'no (fresh; checkpoints written each cycle)' )"
echo ""

if command -v nvidia-smi >/dev/null 2>&1; then
    VRAM_MB=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d ' ')
    echo "  GPU VRAM: ${VRAM_MB:-?} MB (this ~4400-atom system wants >=24 GB; shrink N_DIAMINES/N_DIACIDS if OOM)"
    if [ "${VRAM_MB:-0}" -lt 20000 ] 2>/dev/null; then
        echo "  WARNING: <20 GB VRAM — reduce system size, e.g. N_DIAMINES=50 N_DIACIDS=50."
    fi
    echo ""
fi

echo "Starting run..."
"${PYTHON}" scripts/run_nylon66.py \
    --seed "${SEED}" \
    --output-dir "${OUTPUT_DIR}" \
    --n-diamines "${N_DIAMINES}" \
    --n-diacids "${N_DIACIDS}" \
    --n-cycles "${N_CYCLES}" \
    --biased-steps "${BIASED_STEPS}" \
    --unbiased-steps "${UNBIASED_STEPS}" \
    --equil-steps "${EQUIL_STEPS}" \
    --f2 "${F2}" \
    --density "${DENSITY}" \
    --temperature "${TEMPERATURE}" \
    --pressure "${PRESSURE}" \
    --minimize \
    --minimize-fmax 1.0 \
    --backend orb \
    --device "${DEVICE}" \
    ${BAROSTAT_FLAG} \
    ${RESUME_FLAG}

echo ""
echo "Run complete. Generating Carothers + stability figures..."
"${PYTHON}" scripts/reproduce_figures.py \
    --trajectory "${OUTPUT_DIR}/trajectory.jsonl" \
    --bonds "${OUTPUT_DIR}/bonds.jsonl" \
    --topology "${OUTPUT_DIR}/topology.jsonl" \
    --summary "${OUTPUT_DIR}/summary.json" \
    --n-reactive-sites "${N_REACTIVE_SITES}" \
    --target-temperature "${TEMPERATURE}" \
    --output-dir "${OUTPUT_DIR}/figures"

echo ""
echo "=== Nylon paper-scale done ==="
echo "  Artifacts: ${OUTPUT_DIR}/  (summary.json has carothers_p / carothers_dpn)"
echo "  Figures:   ${OUTPUT_DIR}/figures/  (dpn_vs_conversion.png = measured DPn vs 1/(1-p))"
echo ""
echo "Resume/extend:"
echo "  RESUME=1 N_CYCLES=$(( N_CYCLES + 100 )) OUTPUT_DIR=${OUTPUT_DIR} bash scripts/run_nylon66_paper_scale.sh"
