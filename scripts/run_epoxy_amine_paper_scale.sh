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
# ─── RUN IT (≥24 GB GPU wanted; estimated ~27 GB @5900 atoms, ~1.9 s/step —
#     extrapolated from 16 GB-machine half-scale measurements; paper-scale
#     numbers not yet measured directly) ──────────────────────────────────
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

# shellcheck source=scripts/_paper_scale_common.sh
source "$(dirname "$0")/_paper_scale_common.sh"

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
build_barostat_flag npt

# RESUME=1 continues from ${OUTPUT_DIR}/checkpoint.pkl after a killed run.
build_resume_flag

# Perf flags: see _paper_scale_common.sh for NO_EMPTY_CACHE/COMPILE rationale
# (decisions.md 追補 2026-07-22/23, runs/scaleup_a + scaleup_matrix). At this
# ~5900-atom scale COMPILE=1 is kept opt-in until the first long soak on a
# production workload (barostat + bond formation + resume); note the NPT
# volume moves here vary graph sizes more than the NVT probe the 1.65x was
# measured on (decisions.md 2026-07-23).
build_perf_flags

echo "=== Epoxy-amine paper-scale run ==="
echo "  Python:         ${PYTHON}"
echo "  Seed:           ${SEED}"
echo "  Output dir:     ${OUTPUT_DIR}"
echo "  Device:         ${DEVICE}"
echo "  System:         ${N_EPOXIES} DGEBA + ${N_AMINES} DETA"
echo "  Ensemble:       $( [ -n "${BAROSTAT_FLAG}" ] && echo 'NVT (--no-barostat)' || echo 'NPT, '"${PRESSURE}"' atm' )"
echo "  Schedule:       ${N_CYCLES} cycles x (${BIASED_STEPS} biased + ${UNBIASED_STEPS} unbiased), f2=${F2}, T=${TEMPERATURE} K"
echo "  Resume:         $( [ -n "${RESUME_FLAG}" ] && echo yes || echo 'no (fresh; checkpoints written each cycle)' )"
echo "  Perf:           empty_cache=$( [ "${NO_EMPTY_CACHE}" = "1" ] && echo off || echo on ), compile=$( [ "${COMPILE}" = "1" ] && echo on || echo off )"
echo ""

check_vram 28000 "epoxy-amine ~5900 atoms; estimated ~27 GB (extrapolated from 16 GB-machine half-scale measurements, paper-scale not yet measured) — shrink N_EPOXIES/N_AMINES if OOM, e.g. N_EPOXIES=50 N_AMINES=25"

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
# n-reactive-sites is intentionally omitted: reproduce_figures.py auto-reads
# the true count from the trajectory header, and passing an explicit value
# here would overwrite that (correct) figure with a stale/manual one.
FIGURES_CMD=(
    "${PYTHON}" scripts/reproduce_figures.py
    --trajectory "${OUTPUT_DIR}/trajectory.jsonl"
    --bonds "${OUTPUT_DIR}/bonds.jsonl"
    --summary "${OUTPUT_DIR}/summary.json"
    --target-temperature "${TEMPERATURE}"
    --timestep-fs 0.25
    --output-dir "${OUTPUT_DIR}/figures"
)
if [ -f "${OUTPUT_DIR}/topology.jsonl" ]; then
    FIGURES_CMD+=(--topology "${OUTPUT_DIR}/topology.jsonl" --species-figures)
else
    echo "  (topology.jsonl not found; skipping --topology/--species-figures)"
fi
if ! "${FIGURES_CMD[@]}"; then
    echo "WARNING: figure generation failed; run manually:"
    echo "  ${PYTHON} scripts/reproduce_figures.py --trajectory ${OUTPUT_DIR}/trajectory.jsonl --bonds ${OUTPUT_DIR}/bonds.jsonl --summary ${OUTPUT_DIR}/summary.json --target-temperature ${TEMPERATURE} --timestep-fs 0.25 --output-dir ${OUTPUT_DIR}/figures"
fi

echo ""
echo "=== Epoxy paper-scale done ==="
echo "  Artifacts: ${OUTPUT_DIR}/  (summary.json has epoxide conversion)"
echo "  Figures:   ${OUTPUT_DIR}/figures/  (species_concentrations, gel_curve, ...)"
echo ""
echo "E2 external comparison (after this run):"
echo "  ${PYTHON} scripts/compare_epoxy_external.py --run-dir ${OUTPUT_DIR} ..."
echo "  (fetch reference data first: ${PYTHON} scripts/fetch_provenzano2025.py)"
echo ""
echo "Resume/extend:"
echo "  RESUME=1 N_CYCLES=$(( N_CYCLES + 50 )) OUTPUT_DIR=${OUTPUT_DIR} bash scripts/run_epoxy_amine_paper_scale.sh"
