# インストールガイド

KAGOME は **Python 3.11+** で動作します。用途に応じて以下の環境を構築してください。

## 環境の全体像

| 環境 | 用途 | OS | 主な依存 |
|---|---|---|---|
| ML 本番環境 | TDBB 重合シミュレーション | Windows / Linux | PyTorch, OrbMol-v2 or MACE |
| 古典準備環境 | OpenFF/OpenMM による構造前処理 | WSL2 / Linux | OpenMM, OpenFF, RDKit |

> **Windows + WSL2 の二段構成を推奨**: OpenFF/OpenMM は Windows ネイティブでは動作しないため、WSL2 内で古典構造準備を行い、Windows 側で GPU を使った ML 本番を実行します。Linux 環境では単一の conda 環境で両方をカバーできます。

---

## ML 本番環境 (Windows / Linux)

### 1. conda 環境の作成

```bash
conda create -n pfpoly-gpu python=3.12 -c conda-forge -y
conda activate pfpoly-gpu
```

### 2. PyTorch のインストール

```bash
# CUDA 12.8 (推奨)
pip install torch==2.10.0 --index-url https://download.pytorch.org/whl/cu128

# CPU のみ (GPU なし)
pip install torch==2.10.0 --index-url https://download.pytorch.org/whl/cpu
```

> PyTorch は **pip でインストール**してください。conda-forge 版は CUDA ランタイムが欠落する場合があります。

### 3. プロジェクトのインストール

```bash
pip install -e ".[dev]"
```

### 4. バックエンドの追加

```bash
# OrbMol-v2 (推奨 -- Apache-2.0, 有機系に高精度)
pip install -e ".[orb]"

# MACE-MP-0 (代替 -- MIT, 汎用だが有機系精度はやや劣る)
pip install -e ".[mace]"
```

### 5. プロットツール (任意)

```bash
pip install -e ".[plot]"
```

### 動作確認

```bash
# スモークテスト (Toy バックエンドで TDBB ワークフローを検証)
python scripts/run_smoke.py

# ユニットテスト
pytest -q tests/unit
```

---

## 古典構造準備環境 (WSL2 / Linux)

OpenFF/OpenMM による古典力場での構造前処理に使います。論文スケール (100+ モノマー) の実行では、密度 0.5 g/mL への圧縮・熱平衡化をこの環境で行ってから ML 本番に渡します。

### WSL2 セットアップ (Windows ユーザー)

```bash
# WSL2 内で実行
conda create -n pfpoly-prep python=3.12 -c conda-forge -y
conda activate pfpoly-prep

# OpenFF/OpenMM スタック
conda install -c conda-forge \
    openff-toolkit-base openff-interchange \
    openff-forcefields openmm rdkit numpy scipy ase -y

# プロジェクトのインストール
pip install -e .
```

### Linux セットアップ (単一環境)

Linux では ML 本番環境に OpenFF/OpenMM を追加インストールすることで、単一環境で完結できます。

```bash
conda activate pfpoly-gpu

# OpenFF/OpenMM を追加
conda install -c conda-forge \
    openff-toolkit-base openff-interchange \
    openff-forcefields openmm rdkit -y
```

---

## 論文スケール実行 (S6) のハードウェア要件

200 モノマー + 10 AIBN (2520 原子) の論文スケール実行には追加要件があります。

| 要件 | 推奨 |
|---|---|
| GPU VRAM | 24 GB 以上 (RTX 3090 / 4090 / A100) |
| CUDA | 12.1 以上 |
| RAM | 32 GB 以上 |
| ストレージ | 20 GB 以上 |

### 実行時間の目安

| システムサイズ | GPU | 1 サイクルあたり | 50 サイクル |
|---|---|---|---|
| 2520 原子 | RTX 4090 (24 GB) | 15-30 分 | 12-25 時間 |
| 2520 原子 | A100 (40/80 GB) | 5-15 分 | 4-12 時間 |
| 1260 原子 (半スケール) | RTX 4060 Ti (16 GB) | 5-10 分 | 4-8 時間 |

> 16 GB GPU では半スケール (100+5, 1260 原子) で実行可能です。`--n-monomers 100 --n-initiators 5` を使用してください。

---

## GPU メモリに関する注意

OrbMol-v2 は近傍グラフのサイズがステップごとに変動し、CUDA メモリアロケータのフラグメンテーションを起こすことがあります。以下の環境変数が自動設定されます:

```
PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
```

VRAM 不足で `CUDA out of memory` が発生する場合は:

1. `--n-monomers` を減らす (半スケール: 100)
2. `--device cpu` で CPU 実行 (低速だが VRAM 不要)
3. より大きな VRAM の GPU を使用

---

## ライセンスチェック

新しい依存を追加する前に、ライセンスの互換性を確認してください。

```bash
python scripts/check_dependency_licenses.py --approved specs/approved_dependencies.yaml
```

詳細は [ライセンスポリシー](license-policy.md) を参照してください。
