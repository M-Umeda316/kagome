#!/usr/bin/env bash
# Epoxy-amine bulk curing — PAPER-SCALE run (Track 2 / E2 comparison feeder).
# Method paper: arXiv:2511.22874 (TDBB workflow); the epoxy-amine bulk system is
# the repo's Track 2 extension (decisions.md 2026-07-08 E0 design, 2026-07-09 E1).
# Comparison target: Provenzano 2025 (ACS Appl. Polym. Mater. 7(8):4876) via
# scripts/compare_epoxy_external.py once this run's conversion curve exists.
#
# System (default): 100 DGEBA + 50 DETA ≈ 5900 atoms (2:1, paper-Table-S3-style
# reproduction formulation). For the E2 external comparison use the Provenzano
# stoichiometric 5:2 formulation instead: N_EPOXIES=100 N_AMINES=40
# (decisions.md 2026-07-12 — E2 protocol; do NOT mix the two in one run dir).
#
# ─── RUN IT (≥24 GB GPU wanted; measured ~27 GB @5900 atoms, ~1.9 s/step) ────
#   # From the repo root, with the MLIP conda env active (python on PATH):
#   bash scripts/run_epoxy_amine_paper_scale.sh
#
#   # Specific interpreter / knobs without editing the file:
#   PYTHON=/path/to/envs/pfpoly-gpu/bin/python N_CYCLES=150 \
#       OUTPUT_DIR=runs/epoxy_paper_seed7 bash scripts/run_epoxy_amine_paper_scale.sh
#
#   # Resume a killed run from its last per-cycle checkpoint:
#   RESUME=1 OUTPUT_DIR=runs/epoxy_paper_seed7 bash scripts/run_epoxy_amine_paper_scale.sh
#
#   # E2 comparison formulation (Provenzano 5:2):
#   N_AMINES=40 OUTPUT_DIR=runs/epoxy_e2_seed7 bash scripts/run_epoxy_amine_paper_scale.sh
#
# Estimated wall-clock: LONG. ~3500 MD steps/cycle over ~5900 atoms at ~1.9 s/step
# ≈ 1.8 h/cycle eager; N_CYCLES=100 is on the order of a week. Fully resumable —
# run in chunks with RESUME=1 and raise N_CYCLES as conversion climbs. Judge
# continuation from summary.json epoxide conversion around cycle 30-50, as with
# the nylon paper-scale run.
#
# ─── Parameter rationale ───────────────────────────────────────────────────
# n_epoxies=100, n_amines=50   : Track 2 paper-scale target (decisions.md
#                                2026-07-09 E1; mid30 = 30+15 validated the
#                                mechanism at 1770 atoms, conv 0.167, zero
#                                spurious dissociations).
# density=0.5                  : Same initial packing as MA/nylon (paper SI S-3
#                                convention adopted for Track 2).
# NPT (barostat on), P=1 atm   : Paper ensemble; MC barostat accounts for bias
#                                energy (decisions.md D2 2026-07-03). The
#                                validated mid30 run was NVT — set NO_BAROSTAT=1
#                                to fall back to that ensemble if the barostat
#                                proves unstable at scale.
# temperature=333              : E1 recipe (run_epoxy_amine.py default; mid30
#                                temps healthy: mean 350 K, zero frames >500 K).
# f2=2.0                       : OrbMol-v2 capture width — E0 scan: paper f2=10
#                                leaves a dead-zone (max bias force ~0 in the
#                                [3,6] Å window); f2=2 bridges it (decisions.md
#                                2026-07-09 E0 result, user-approved).
# friction=0.01 /fs            : Exotherm-safe recipe (dE ≈ −29.5 kcal/mol ring
#                                opening ≈ MA-level; decisions.md 2026-07-09).
# biased=2000, unbiased=1500   : E1/mid30 schedule (mirrors the MA recipe).
# minimize + equil=2000        : Pre-TDBB relaxation of the compressed box
#                                (decisions.md 2026-07-08 — segfault guard).
# checkpoint every cycle       : Resumable long run (RESUME=1).
#
# NOT in git (provision separately on a new machine):
#   models/orbmol-v2-teqabfhg-20260523.ckpt  (gitignored; backends/orb_backend.py
#   resolves it relative to the repo root), the pfpoly-gpu conda env
#   (docs/installation.md), and runs/ artifacts.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PYTHONPATH="${REPO_ROOT}:${REPO_ROOT}/src:${PYTHONPATH:-}"
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"
export KMP_DUPLICATE_LIB_OK="${KMP_DUPLICATE_LIB_OK:-TRUE}"

# Interpreter: override with PYTHON=/abs/path/to/python for a specific conda env.
PYTHON="${PYTHON:-python}"

SEED="${SEED:-7}"
OUTPUT_DIR="${OUTPUT_DIR:-runs/epoxy_amine_paper_scale_seed${SEED}}"
DEVICE="${DEVICE:-cuda}"
N_EPOXIES="${N_EPOXIES:-100}"
N_AMINES="${N_AMINES:-50}"
N_CYCLES="${N_CYCLES:-100}"
BIASED_STEPS="${BIASED_STEPS:-2000}"
UNBIASED_STEPS="${UNBIASED_STEPS:-1500}"
EQUIL_STEPS="${EQUIL_STEPS:-2000}"
F2="${F2:-2}"
DENSITY="${DENSITY:-0.5}"
TEMPERATURE="${TEMPERATURE:-333.0}"
PRESSURE="${PRESSURE:-1.0}"
FRICTION_PER_FS="${FRICTION_PER_FS:-0.01}"

# NPT by default (paper ensemble). NO_BAROSTAT=1 -> NVT (the mid30-validated ensemble).
BAROSTAT_FLAG=""
if [ "${NO_BAROSTAT:-0}" = "1" ]; then BAROSTAT_FLAG="--no-barostat"; fi

# RESUME=1 continues from ${OUTPUT_DIR}/checkpoint.pkl after a killed run.
RESUME_FLAG=""
if [ "${RESUME:-0}" = "1" ]; then RESUME_FLAG="--resume"; fi

# ─── Perf flags (decisions.md 追補 2026-07-22/23, runs/scaleup_a + scaleup_matrix) ─
# NO_EMPTY_CACHE=1 (default): skip per-step torch.cuda.empty_cache(). On WSL,
#   expandable_segments absorbs fragmentation (reserved flat, 1.29-1.41x faster)
#   — the measured recommended operation for WSL runs. Set NO_EMPTY_CACHE=0 only
#   on Windows-native python, where expandable_segments is a no-op.
# COMPILE=1 (opt-in, default 0): torch.compile — 1.24x alone, 1.65x combined
#   with NO_EMPTY_CACHE=1, and ~8% less VRAM (welcome at this ~5900-atom scale).
#   Forces match eager within the TF32 noise floor. Kept opt-in until the first
#   long soak on a production workload (barostat + bond formation + resume);
#   note the NPT volume moves here vary graph sizes more than the NVT probe
#   the 1.65x was measured on (decisions.md 2026-07-23).
NO_EMPTY_CACHE="${NO_EMPTY_CACHE:-1}"
COMPILE="${COMPILE:-0}"
PERF_FLAGS=""
if [ "${NO_EMPTY_CACHE}" = "1" ]; then PERF_FLAGS="--no-empty-cache"; fi
if [ "${COMPILE}" = "1" ]; then PERF_FLAGS="${PERF_FLAGS} --compile"; fi

# Epoxide-basis conversion denominator for figures: each DGEBA carries 2 epoxide
# rings -> 2*N_EPOXIES (mid30 manifest: n_reactive_sites=60 at N_EPOXIES=30).
N_REACTIVE_SITES=$(( 2 * N_EPOXIES ))

echo "=== Epoxy-amine paper-scale run ==="
echo "  Python:         ${PYTHON}"
echo "  Seed:           ${SEED}"
echo "  Output dir:     ${OUTPUT_DIR}"
echo "  Device:         ${DEVICE}"
echo "  System:         ${N_EPOXIES} DGEBA + ${N_AMINES} DETA (${N_REACTIVE_SITES} epoxide sites)"
echo "  Ensemble:       $( [ -n "${BAROSTAT_FLAG}" ] && echo 'NVT (--no-barostat)' || echo 'NPT, '"${PRESSURE}"' atm' )"
echo "  Schedule:       ${N_CYCLES} cycles x (${BIASED_STEPS} biased + ${UNBIASED_STEPS} unbiased), f2=${F2}, T=${TEMPERATURE} K"
echo "  Resume:         $( [ -n "${RESUME_FLAG}" ] && echo yes || echo 'no (fresh; checkpoints written each cycle)' )"
echo "  Perf:           empty_cache=$( [ "${NO_EMPTY_CACHE}" = "1" ] && echo off || echo on ), compile=$( [ "${COMPILE}" = "1" ] && echo on || echo off )"
echo ""

if command -v nvidia-smi >/dev/null 2>&1; then
    VRAM_MB=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d ' ')
    echo "  GPU VRAM: ${VRAM_MB:-?} MB (measured ~27 GB @5900 atoms; shrink N_EPOXIES/N_AMINES if OOM)"
    if [ "${VRAM_MB:-0}" -lt 24000 ] 2>/dev/null; then
        echo "  WARNING: <24 GB VRAM — reduce system size, e.g. N_EPOXIES=50 N_AMINES=25."
    fi
    echo ""
fi

echo "Starting run..."
"${PYTHON}" scripts/run_epoxy_amine.py \
    --seed "${SEED}" \
    --output-dir "${OUTPUT_DIR}" \
    --n-epoxies "${N_EPOXIES}" \
    --n-amines "${N_AMINES}" \
    --n-cycles "${N_CYCLES}" \
    --biased-steps "${BIASED_STEPS}" \
    --unbiased-steps "${UNBIASED_STEPS}" \
    --equil-steps "${EQUIL_STEPS}" \
    --f2 "${F2}" \
    --density "${DENSITY}" \
    --temperature "${TEMPERATURE}" \
    --pressure "${PRESSURE}" \
    --friction-per-fs "${FRICTION_PER_FS}" \
    --minimize \
    --minimize-fmax 1.0 \
    --backend orb \
    --device "${DEVICE}" \
    ${PERF_FLAGS} \
    ${BAROSTAT_FLAG} \
    ${RESUME_FLAG}

echo ""
echo "Run complete. Generating conversion + network + stability figures..."
"${PYTHON}" scripts/reproduce_figures.py \
    --trajectory "${OUTPUT_DIR}/trajectory.jsonl" \
    --bonds "${OUTPUT_DIR}/bonds.jsonl" \
    --topology "${OUTPUT_DIR}/topology.jsonl" \
    --summary "${OUTPUT_DIR}/summary.json" \
    --n-reactive-sites "${N_REACTIVE_SITES}" \
    --target-temperature "${TEMPERATURE}" \
    --species-figures \
    --output-dir "${OUTPUT_DIR}/figures"

echo ""
echo "=== Epoxy paper-scale done ==="
echo "  Artifacts: ${OUTPUT_DIR}/  (summary.json has epoxide conversion)"
echo "  Figures:   ${OUTPUT_DIR}/figures/  (species_concentrations, gel_curve, ...)"
echo ""
echo "E2 external comparison (after this run):"
echo "  python scripts/compare_epoxy_external.py --run-dir ${OUTPUT_DIR} ..."
echo "  (fetch reference data first: python scripts/fetch_provenzano2025.py)"
echo ""
echo "Resume/extend:"
echo "  RESUME=1 N_CYCLES=$(( N_CYCLES + 50 )) OUTPUT_DIR=${OUTPUT_DIR} bash scripts/run_epoxy_amine_paper_scale.sh"
