#!/usr/bin/env bash
# Scale-up 4-condition matrix (empty_cache x compile) + drift bracket.
# 5 sequential paper-scale arms; see specs/decisions.md 追補 2026-07-23 (c).
set -u
cd /mnt/c/Users/shanu/Documents/Python/kagome
export PYTHONPATH=.:src
PY="$HOME/miniconda3/envs/pfpoly-gpu/bin/python"
COMMON="--device cuda --systems vinyl_methyl_acrylate --md-steps 300 --warmup-steps 30 --mem-sample-every 10"
mkdir -p runs/scaleup_matrix

run_arm() {
    local name="$1"; shift
    echo "[$(date +%H:%M:%S)] START $name"
    $PY -m scripts.profile_vram $COMMON "$@" \
        --output-dir "runs/scaleup_matrix/$name" \
        > "runs/scaleup_matrix/$name.log" 2>&1 \
        || echo "FAILED $name"
    echo "[$(date +%H:%M:%S)] END $name"
}

run_arm m1_ecOn_eager
run_arm m2_ecOff_eager --no-empty-cache
run_arm m3_ecOn_compile --compile
run_arm m4_ecOff_compile --no-empty-cache --compile
run_arm m5_ecOn_eager_bracket
echo ALL_DONE
