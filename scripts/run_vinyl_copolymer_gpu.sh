#!/usr/bin/env bash
# Vinyl COPOLYMERIZATION (methyl acrylate + methyl methacrylate) — GPU run.
# Paper: arXiv:2511.22874, Section 2/3, Table S1 (4-group ij+ik+jl criterion).
# Design: specs/decisions.md 2026-07-16 (documented copolymer extension — both
#         species' alpha-Cs share one vinyl_alpha_C group; no reactivity-ratio
#         bias, so the sequence is left to the MLIP energetics).
#
# System (default): 50 acrylate + 50 methacrylate + 5 isobutyronitrile radicals
#         (closed-shell initiator model; ~1500 atoms). Env-overridable below.
#
# ─── RUN IT (on a GPU machine, e.g. WSL pfpoly-gpu) ─────────────────────────
#   # From the repo root, with the MLIP conda env active (python on PATH):
#   bash scripts/run_vinyl_copolymer_gpu.sh
#
#   # Point at a specific interpreter / tune knobs without editing the file:
#   PYTHON=/path/to/envs/pfpoly-gpu/bin/python \
#       N_ACRYLATE=8 N_METHACRYLATE=8 N_INITIATORS=2 N_CYCLES=10 \
#       OUTPUT_DIR=runs/copoly_seed7 bash scripts/run_vinyl_copolymer_gpu.sh
#
#   # Resume a killed run from its last cycle checkpoint (written every cycle;
#   # resume skips build/minimize/equilibration and restarts at the saved cycle):
#   RESUME=1 OUTPUT_DIR=runs/copoly_seed7 bash scripts/run_vinyl_copolymer_gpu.sh
#
# Estimated wall-clock: each cycle is ~3500 MD steps over ~1500 atoms; budget on
# the order of minutes/cycle on a modern GPU. The run is fully resumable — run it
# in chunks with RESUME=1 and raise N_CYCLES as conversion climbs.
#
# ─── Parameter rationale (all from the validated S6 vinyl recipe) ───────────
# f2=2.0             : OrbMol-v2 capture width. Paper default 10 (for PFP) leaves
#                      a dead-zone — the addition bias force is ~0 across the
#                      [3,6] Å candidate window on OrbMol's PES, so 0 formations.
#                      f2=2 bridges it (decisions.md 2026-06-26).
# f1_max 250/125     : Formation/dissociation bias caps, same recipe.
# friction=0.01      : Langevin friction. Lowering f2 injects bias work that,
#                      with addition heat, accumulates; 0.01 (vs 0.001) dissipates
#                      it and pins T near 333 K (decisions.md 2026-06-26).
# timestep=0.25 fs   : Paper value. REQUIRED for reactive multi-radical stability
#                      — 1.0 fs numerically explodes the open-shell melt
#                      (decisions.md 2026-06-25, validity-domain §2.1/§3).
# density=0.5        : Paper SI S-3 (acrylate melt). temperature=333 K (Table S1).
# NVT (no barostat)  : Validated ensemble for open-shell vinyl systems.
# biased=2000, unbiased=1500, equil=2000 : S6 validated schedule.
# minimize fmax=1.0  : Pre-TDBB relaxation of the compressed box.
#
# NOTE: this uses the CLOSED-SHELL isobutyronitrile initiator (radical C carries a
# placeholder H; the first addition over-coordinates by 1 in the raw geometry,
# idealized in the emitted topology). Copolymer trajectory bond-topology output is
# off (single-stride helpers don't fit the mixed layout); conversion / reaction
# events are unaffected. See specs/decisions.md 2026-07-16.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PYTHONPATH="${REPO_ROOT}:${REPO_ROOT}/src:${PYTHONPATH:-}"
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"
export KMP_DUPLICATE_LIB_OK="${KMP_DUPLICATE_LIB_OK:-TRUE}"

# Interpreter: override with PYTHON=/abs/path/to/python for a specific conda env.
PYTHON="${PYTHON:-python}"

SEED="${SEED:-7}"
OUTPUT_DIR="${OUTPUT_DIR:-runs/copolymer_seed${SEED}}"
DEVICE="${DEVICE:-cuda}"
# BACKEND=toy DEVICE=cpu runs a GPU-free smoke of this exact pipeline.
BACKEND="${BACKEND:-orb}"
N_ACRYLATE="${N_ACRYLATE:-50}"
N_METHACRYLATE="${N_METHACRYLATE:-50}"
N_INITIATORS="${N_INITIATORS:-5}"
N_CYCLES="${N_CYCLES:-30}"
BIASED_STEPS="${BIASED_STEPS:-2000}"
UNBIASED_STEPS="${UNBIASED_STEPS:-1500}"
EQUIL_STEPS="${EQUIL_STEPS:-2000}"
TIMESTEP_FS="${TIMESTEP_FS:-0.25}"
F2="${F2:-2}"
F1_MAX_FORMATION="${F1_MAX_FORMATION:-250.0}"
F1_MAX_DISSOCIATION="${F1_MAX_DISSOCIATION:-125.0}"
FRICTION_PER_FS="${FRICTION_PER_FS:-0.01}"
DENSITY="${DENSITY:-0.5}"
TEMPERATURE="${TEMPERATURE:-333.0}"
MINIMIZE_FMAX="${MINIMIZE_FMAX:-1.0}"

# NVT by default (validated). BAROSTAT=1 -> NPT at PRESSURE atm.
BAROSTAT_FLAG=""
if [ "${BAROSTAT:-0}" = "1" ]; then BAROSTAT_FLAG="--barostat --pressure ${PRESSURE:-1.0}"; fi

# RESUME=1 continues from ${OUTPUT_DIR}/checkpoint.pkl after a killed run.
RESUME_FLAG=""
if [ "${RESUME:-0}" = "1" ]; then RESUME_FLAG="--resume"; fi

N_MONOMERS=$(( N_ACRYLATE + N_METHACRYLATE ))

echo "=== Vinyl copolymer (acrylate + methacrylate) run ==="
echo "  Python:         ${PYTHON}"
echo "  Seed:           ${SEED}"
echo "  Output dir:     ${OUTPUT_DIR}"
echo "  Device:         ${DEVICE}"
echo "  Composition:    ${N_ACRYLATE} acrylate + ${N_METHACRYLATE} methacrylate + ${N_INITIATORS} initiator (${N_MONOMERS} monomers)"
echo "  Ensemble:       $( [ -n "${BAROSTAT_FLAG}" ] && echo 'NPT' || echo 'NVT' )"
echo "  Schedule:       ${N_CYCLES} cycles x (${BIASED_STEPS} biased + ${UNBIASED_STEPS} unbiased), dt=${TIMESTEP_FS} fs"
echo "  TDBB:           f2=${F2}, f1_max_form=${F1_MAX_FORMATION}, f1_max_dissoc=${F1_MAX_DISSOCIATION}, friction=${FRICTION_PER_FS}"
echo "  T=${TEMPERATURE} K, density=${DENSITY} g/mL"
echo "  Resume:         $( [ -n "${RESUME_FLAG}" ] && echo yes || echo 'no (fresh; checkpoints written each cycle)' )"
echo ""

if command -v nvidia-smi >/dev/null 2>&1; then
    VRAM_MB=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d ' ')
    echo "  GPU VRAM: ${VRAM_MB:-?} MB (~1500-atom default wants >=16 GB; shrink N_ACRYLATE/N_METHACRYLATE if OOM)"
    if [ "${VRAM_MB:-0}" -lt 12000 ] 2>/dev/null; then
        echo "  WARNING: <12 GB VRAM — reduce system, e.g. N_ACRYLATE=20 N_METHACRYLATE=20 N_INITIATORS=2."
    fi
    echo ""
fi

echo "Starting run..."
"${PYTHON}" scripts/run_vinyl_copolymer.py \
    --seed "${SEED}" \
    --output-dir "${OUTPUT_DIR}" \
    --n-acrylate "${N_ACRYLATE}" \
    --n-methacrylate "${N_METHACRYLATE}" \
    --n-initiators "${N_INITIATORS}" \
    --n-cycles "${N_CYCLES}" \
    --biased-steps "${BIASED_STEPS}" \
    --unbiased-steps "${UNBIASED_STEPS}" \
    --equil-steps "${EQUIL_STEPS}" \
    --timestep-fs "${TIMESTEP_FS}" \
    --f2 "${F2}" \
    --f1-max-formation "${F1_MAX_FORMATION}" \
    --f1-max-dissociation "${F1_MAX_DISSOCIATION}" \
    --friction-per-fs "${FRICTION_PER_FS}" \
    --density "${DENSITY}" \
    --temperature "${TEMPERATURE}" \
    --minimize \
    --minimize-fmax "${MINIMIZE_FMAX}" \
    --backend "${BACKEND}" \
    --device "${DEVICE}" \
    ${BAROSTAT_FLAG} \
    ${RESUME_FLAG}

echo ""
echo "Run complete. Generating figures..."
"${PYTHON}" scripts/reproduce_figures.py \
    --trajectory "${OUTPUT_DIR}/trajectory.jsonl" \
    --bonds "${OUTPUT_DIR}/bonds.jsonl" \
    --summary "${OUTPUT_DIR}/summary.json" \
    --n-reactive-sites "${N_MONOMERS}" \
    --target-temperature "${TEMPERATURE}" \
    --timestep-fs "${TIMESTEP_FS}" \
    --output-dir "${OUTPUT_DIR}/figures"

echo ""
echo "=== Copolymer run done ==="
echo "  Artifacts: ${OUTPUT_DIR}/  (summary.json has n_acrylate / n_methacrylate / confirmed_formations)"
echo "  Figures:   ${OUTPUT_DIR}/figures/"
echo ""
echo "Resume/extend:"
echo "  RESUME=1 N_CYCLES=$(( N_CYCLES + 20 )) OUTPUT_DIR=${OUTPUT_DIR} bash scripts/run_vinyl_copolymer_gpu.sh"
