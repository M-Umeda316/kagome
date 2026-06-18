# S6 実行環境セットアップガイド

## 概要

S6 は論文スケール（200 monomer + 10 AIBN = ≈2520 atoms）の実行です。
OrbMol-v2 のシングルステップ VRAM フットプリントが約 9.5 GB のため、**24 GB 以上の GPU** が必要です。

---

## 必要スペック

| 項目 | 要件 |
|------|------|
| GPU VRAM | ≥ 24 GB（RTX 3090 / 4090 / A10G / A100 等） |
| CUDA | ≥ 12.1 |
| RAM | ≥ 32 GB 推奨 |
| ストレージ | ≥ 20 GB（conda env + run artifacts） |

---

## 環境クローン手順

### 現在の環境情報（Windows pfpoly-gpu）

主要パッケージバージョン:
- Python 3.12
- PyTorch 2.10.0+cu128（pip 経由, CUDA 12.8 対応）
- orb-models 0.7.0（pip 経由）
- numpy 2.4.4, scipy 1.17.1, matplotlib 3.11.0
- ase 3.28.0, rdkit 2026.3.3

### Linux/WSL2 での構築（推奨手順）

```bash
# 1. リポジトリをクローン
git clone <repo_url> pfpoly
cd pfpoly

# 2. conda 環境を作成（Python 3.12）
conda create -n pfpoly-gpu python=3.12 -c conda-forge -y
conda activate pfpoly-gpu

# 3. PyTorch を pip でインストール（CUDA 12.8 対応ビルド）
#    Note: CUDA 12.1 の場合は cu121 に変更
pip install torch==2.10.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu128

# 4. orb-models と科学系パッケージを pip でインストール
pip install orb-models==0.7.0
pip install ase rdkit scipy matplotlib numpy

# 5. プロジェクト依存をインストール
pip install -e ".[dev]"   # または requirements があれば: pip install -r requirements.txt

# 6. GPU が認識されていることを確認
python -c "import torch; print('CUDA available:', torch.cuda.is_available()); print('Device:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A')"
```

### conda environment.yml からの復元（代替手順）

```bash
# Windows 環境から export する場合（別ハード用なので platform 指定が必要）
conda env export -n pfpoly-gpu --no-builds | grep -v "^prefix:" > environment.yml

# Linux 側で復元する際は pip 部分のみ参考にし、
# torch は上記の pip コマンドで個別インストールすること
# (Windows ビルドの torch は Linux では使えない)
```

---

## 実行手順

```bash
# repo root から
conda activate pfpoly-gpu

# 環境変数設定（OrbMol-v2 VRAM 断片化対策）
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export KMP_DUPLICATE_LIB_OK=TRUE

# S6 実行（デフォルト: seed=7, 50 cycles）
bash scripts/run_s6_paper_scale.sh

# カスタム設定で実行する場合
SEED=42 OUTPUT_DIR=runs/s6_seed42 N_CYCLES=100 bash scripts/run_s6_paper_scale.sh

# GPU が 16 GB しかない場合（論文スケール半分: 100+5）
DEVICE=cuda python scripts/run_vinyl_aibn.py \
    --seed 7 --output-dir runs/s6_half_scale \
    --n-monomers 100 --n-initiators 5 \
    --activation --activation-f2 0.3 --activation-f1-max 250 \
    --f2 5.0 --density 0.5 --temperature 333.0 --no-barostat \
    --backend orb --device cuda \
    --n-cycles 30 --biased-steps 2000 --unbiased-steps 500 \
    --equil-steps 2000 --timestep-fs 1.0
```

---

## ポータビリティに関する注意

### pip 経由 PyTorch の理由

現在の環境では PyTorch は **pip 経由**（`torch==2.10.0+cu128`）でインストールされています。
これは conda-forge の PyTorch ビルドが CUDA バージョンの柔軟性に欠けるため、
公式 PyTorch wheel を直接使う方が確実だからです。

**conda-forge 経由の PyTorch で GPU は使えるか？**
- 原理的には可能ですが、パッケージ名と CUDA バージョンの対応が複雑
- conda-forge の `pytorch-gpu` メタパッケージは存在しますが、
  特定の CUDA バージョン（12.1 等）向けビルドが常に最新とは限らない
- **推奨**: `pip install torch --index-url https://download.pytorch.org/whl/cu128` の方が
  確実に最新 CUDA 対応ビルドを取得できる
- インストール後に必ず `torch.cuda.is_available()` で確認すること

### WSL2 での CUDA 動作確認

```bash
# nvidia-smi が WSL2 から見えるか確認
nvidia-smi

# Python から CUDA を確認
python -c "import torch; print(torch.version.cuda); print(torch.cuda.is_available())"
```

WSL2 では Windows ホスト側の NVIDIA ドライバ（≥ 470.76）が CUDA をサポートしていれば
WSL2 内から GPU が見えます。CUDA Toolkit は WSL2 側に別途インストール不要（ドライバ経由で利用）。

---

## 実行時間の目安

| システム | GPU | 予想時間/cycle | 50 cycles |
|---------|-----|----------------|-----------|
| 200+10 (2520 atoms) | RTX 4090 (24 GB) | ~15-30 min | 12-25 h |
| 200+10 (2520 atoms) | A100 (40/80 GB) | ~5-15 min | 4-12 h |
| 100+5 (1260 atoms) | RTX 4060 Ti (16 GB) | ~5-10 min | 4-8 h |

---

## 成果物の確認

実行後に以下のファイルが生成されます:

```
runs/s6_paper_scale_seed7/
├── manifest.json        # seed, git SHA, backend, timestamp
├── summary.json         # n_formations, activation_dissociations 等
├── trajectory.jsonl     # 全フレームの位置・エネルギー・温度
├── bonds.jsonl          # 全反応イベント（formation/dissociation）
└── figures/
    ├── conversion_vs_step.{png,pdf}   # α(t) + Eq.11 フィット
    ├── energy_vs_step.{png,pdf}       # biased/unbiased エネルギー分離
    ├── base_energy.{png,pdf}          # base ポテンシャル推移
    └── temperature_vs_step.{png,pdf}  # NVT 温度安定性
```

figure-comparison.md に結果を追記してください（S5 最終確認項目）。
