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
#   # Well-mixed measurement mode (NOT paper-faithful) — the answer to the
#   # no-mixing paper-scale plateau at p≈12% (decisions.md 2026-08-04):
#   MIX=1 MIX_PLATFORM=CUDA OUTPUT_DIR=runs/nylon_mix_seed7 bash scripts/run_nylon66_paper_scale.sh
#   # NOTE: mix_ps=25 was derived on a 40-monomer system (WM-P5b), NOT at paper
#   # scale. Check the first few dozen cycles before leaning on it: mixing.jsonl
#   # RMS displacement should show real diffusion, and summary.json
#   # n_mixing_skipped must stay near zero (a non-trivial tally means the
#   # measurement is compromised). See decisions.md 2026-08-04 / 追補 2026-07-22.
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
# friction=0.01 /fs             : Bias-heat dissipation for the f2=2 recipe, same
#                                 rationale as the epoxy-amine/vinyl runs (dissolves
#                                 the addition/formation bias work so T stays near
#                                 target; decisions.md 2026-07-30 addendum).
# checkpoint every cycle        : Resumable long run (RESUME=1).

set -euo pipefail

# shellcheck source=scripts/_paper_scale_common.sh
source "$(dirname "$0")/_paper_scale_common.sh"

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
FRICTION_PER_FS="${FRICTION_PER_FS:-0.01}"

# NPT by default (paper). NO_BAROSTAT=1 -> NVT (the small-run-validated ensemble).
build_barostat_flag npt

# RESUME=1 continues from ${OUTPUT_DIR}/checkpoint.pkl after a killed run.
build_resume_flag

# Perf flags: see _paper_scale_common.sh for NO_EMPTY_CACHE/COMPILE rationale
# (decisions.md 追補 2026-07-22/23, runs/scaleup_a + scaleup_matrix). At this
# ~4400-atom scale COMPILE=1 is kept opt-in until the first long soak on a
# production workload (barostat + bond formation + resume); note the NPT
# volume moves here vary graph sizes more than the NVT probe the 1.65x was
# measured on (decisions.md 2026-07-23).
build_perf_flags

# MIX=1 (+ optional MIX_PS / MIX_SETTLE_STEPS / MIX_PLATFORM) enables the
# per-cycle classical mixing stage; see build_mix_flags in
# _paper_scale_common.sh and the RUN IT note above (decisions.md 2026-08-04).
build_mix_flags

echo "=== Nylon-6,6 paper-scale run ==="
echo "  Python:         ${PYTHON}"
echo "  Seed:           ${SEED}"
echo "  Output dir:     ${OUTPUT_DIR}"
echo "  Device:         ${DEVICE}"
echo "  System:         ${N_DIAMINES} diamine + ${N_DIACIDS} diacid"
echo "  Ensemble:       $( [ -n "${BAROSTAT_FLAG}" ] && echo 'NVT (--no-barostat)' || echo 'NPT, '"${PRESSURE}"' atm' )"
echo "  Schedule:       ${N_CYCLES} cycles x (${BIASED_STEPS} biased + ${UNBIASED_STEPS} unbiased), f2=${F2}, T=${TEMPERATURE} K"
echo "  Resume:         $( [ -n "${RESUME_FLAG}" ] && echo yes || echo 'no (fresh; checkpoints written each cycle)' )"
echo "  Perf:           empty_cache=$( [ "${NO_EMPTY_CACHE}" = "1" ] && echo off || echo on ), compile=$( [ "${COMPILE}" = "1" ] && echo on || echo off )"
echo "  Mixing:         $(mix_status)"
echo ""

check_vram 24000 "nylon66 ~4400 atoms; wants >=24 GB — shrink N_DIAMINES/N_DIACIDS if OOM, e.g. N_DIAMINES=50 N_DIACIDS=50"

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
    --friction-per-fs "${FRICTION_PER_FS}" \
    --minimize \
    --minimize-fmax 1.0 \
    --backend orb \
    --device "${DEVICE}" \
    ${PERF_FLAGS} \
    ${MIX_FLAGS} \
    ${BAROSTAT_FLAG} \
    ${RESUME_FLAG}

echo ""
echo "Run complete. Generating Carothers + stability figures..."
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
    FIGURES_CMD+=(--topology "${OUTPUT_DIR}/topology.jsonl")
else
    echo "  (topology.jsonl not found; skipping --topology)"
fi
if ! "${FIGURES_CMD[@]}"; then
    echo "WARNING: figure generation failed; run manually:"
    echo "  ${PYTHON} scripts/reproduce_figures.py --trajectory ${OUTPUT_DIR}/trajectory.jsonl --bonds ${OUTPUT_DIR}/bonds.jsonl --summary ${OUTPUT_DIR}/summary.json --target-temperature ${TEMPERATURE} --timestep-fs 0.25 --output-dir ${OUTPUT_DIR}/figures"
fi

echo ""
echo "=== Nylon paper-scale done ==="
echo "  Artifacts: ${OUTPUT_DIR}/  (summary.json has carothers_p / carothers_dpn)"
echo "  Figures:   ${OUTPUT_DIR}/figures/  (dpn_vs_conversion.png = measured DPn vs 1/(1-p))"
echo ""
echo "Resume/extend:"
echo "  RESUME=1 N_CYCLES=$(( N_CYCLES + 100 )) OUTPUT_DIR=${OUTPUT_DIR} bash scripts/run_nylon66_paper_scale.sh"
