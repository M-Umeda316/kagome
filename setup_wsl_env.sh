#!/usr/bin/env bash
# WSL側に pfpoly-gpu conda環境を作成するスクリプト
# 使い方: wsl -e bash -l /mnt/c/Users/shanu/Documents/Python/pfpoly/setup_wsl_env.sh
set -euo pipefail

ENV_NAME="pfpoly-gpu"
PYTHON_VERSION="3.12"

# ── conda 初期化 ──────────────────────────────────────────────────────────
source ~/miniconda3/etc/profile.d/conda.sh

# ── 既存環境があれば削除確認 ──────────────────────────────────────────────
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "環境 ${ENV_NAME} は既に存在します。再作成しますか? [y/N]"
    read -r ans
    if [[ "$ans" =~ ^[Yy]$ ]]; then
        conda env remove -n "$ENV_NAME" -y
    else
        echo "中止しました。"
        exit 0
    fi
fi

# ── conda 環境作成 (rdkit + openmm は conda-forge から) ───────────────────
echo "=== conda 環境を作成中 ==="
conda create -n "$ENV_NAME" -y -c conda-forge \
    "python=${PYTHON_VERSION}" \
    rdkit \
    openmm \
    openff-toolkit \
    openff-forcefields \
    h5py

conda activate "$ENV_NAME"

# ── PyTorch (CUDA 12.6) ──────────────────────────────────────────────────
echo "=== PyTorch (cu126) をインストール中 ==="
pip install torch torchvision torchaudio \
    --index-url https://download.pytorch.org/whl/cu126

# ── ML バックエンド ───────────────────────────────────────────────────────
echo "=== ML バックエンドをインストール中 ==="
pip install mace-torch "ase>=3.22"
pip install orb-models

# ── 解析・ユーティリティ ──────────────────────────────────────────────────
echo "=== 解析ライブラリをインストール中 ==="
pip install \
    "numpy>=1.24" \
    "scipy>=1.10" \
    "matplotlib>=3.7" \
    "pandas>=2.0" \
    pyyaml \
    tqdm

# ── 開発用 ────────────────────────────────────────────────────────────────
echo "=== 開発ツールをインストール中 ==="
pip install pytest

# ── kagome 本体を editable install ────────────────────────────────────────
echo "=== kagome を editable install 中 ==="
PFPOLY_DIR="/mnt/c/Users/shanu/Documents/Python/pfpoly"
pip install -e "$PFPOLY_DIR"

# ── 検証 ──────────────────────────────────────────────────────────────────
echo ""
echo "=== インストール検証 ==="
python -c "
import torch
print(f'torch:       {torch.__version__}')
print(f'CUDA avail:  {torch.cuda.is_available()}')
if torch.cuda.is_available():
    print(f'CUDA device: {torch.cuda.get_device_name(0)}')
    print(f'CUDA ver:    {torch.version.cuda}')

import numpy; print(f'numpy:       {numpy.__version__}')
import scipy; print(f'scipy:       {scipy.__version__}')
import matplotlib; print(f'matplotlib:  {matplotlib.__version__}')
import ase; print(f'ase:         {ase.__version__}')
import rdkit; print(f'rdkit:       {rdkit.__version__}')

try:
    import openmm; print(f'openmm:      {openmm.__version__}')
except Exception as e:
    print(f'openmm:      SKIP ({e})')

try:
    from mace.calculators import mace_mp
    print('mace-torch:  OK')
except Exception as e:
    print(f'mace-torch:  SKIP ({e})')

try:
    import orb_models; print(f'orb-models:  OK')
except Exception as e:
    print(f'orb-models:  SKIP ({e})')

print()
print('全依存関係のインストールが完了しました。')
"

echo ""
echo "=== 完了 ==="
echo "使い方: conda activate ${ENV_NAME}"
echo "テスト: conda run -n ${ENV_NAME} pytest tests/ -x"
