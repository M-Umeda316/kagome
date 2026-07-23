# S6 実行環境セットアップガイド

## 概要

S6 は論文スケール（200 monomer + 10 AIBN = ≈2520 atoms）の実行です。
OrbMol-v2 のこのスケールでの実測デバイスピークは **~9.5-9.8 GB**（WSL +
`expandable_segments`、runs/scaleup_a/c、decisions.md 追補 2026-07-22/23）で、
**WSL 上なら 16 GB GPU で足ります**。かつての「24 GB 以上必要」という指針は
Windows native（`expandable_segments` が no-op で断片化ハングした環境、
2026-06-15）由来のものです。**実行は WSL(Linux) を必須とします** — native
Windows の pfpoly 環境は CPU ビルドの torch になり orb が動きません。

---

## 必要スペック

| 項目 | 要件 |
|------|------|
| GPU VRAM | ≥ 16 GB（WSL/Linux + expandable_segments 前提。実測ピーク ~9.8 GB） |
| OS | Linux または WSL2（Windows native は不可 — 上記） |
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

### Linux / WSL2（唯一のサポート環境）

旧 Windows PowerShell 版ランチャー（`run_s6_paper_scale.ps1`）は、native
Windows で orb が動かないため 2026-07-23 に削除しました。実行は常に WSL の
`pfpoly-gpu` 環境から `run_s6_paper_scale.sh` を使います。

```bash
# repo root から
conda activate pfpoly-gpu

# S6 実行（デフォルト: seed=7, 50 cycles）
# PYTORCH_CUDA_ALLOC_CONF / KMP_DUPLICATE_LIB_OK / PYTHONPATH はスクリプト内で自動設定
bash scripts/run_s6_paper_scale.sh

# カスタム設定で実行する場合
SEED=42 OUTPUT_DIR=runs/s6_seed42 N_CYCLES=100 bash scripts/run_s6_paper_scale.sh

# 中断したランの再開（毎サイクル checkpoint が書かれる）
RESUME=1 OUTPUT_DIR=runs/s6_seed42 bash scripts/run_s6_paper_scale.sh

# 速度フラグ（decisions.md 2026-07-22/23）: 既定で --no-empty-cache が付く。
# torch.compile を試すとき（併用で ~1.65x、本番 soak は未実施の opt-in）:
COMPILE=1 bash scripts/run_s6_paper_scale.sh
```

半分スケール（100+5）を明示的に走らせたい場合の検証済みレシピ（f2=2 /
dt=0.25 fs / friction=0.01。旧記載の f2=5, dt=1.0 fs は開殻反応系で温度発散
するため使用禁止 — decisions.md 2026-06-25/26, validity-domain.md §2.1/§3）は
`scripts/run_s6_paper_scale.sh` ヘッダのコマンド例を参照。

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
