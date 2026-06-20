# KAGOME

**Kinetic Accelerated Growth Orchestrated by Molecular Engine**

論文 "Ready-to-Use Polymerization Simulations Combining Universal Machine Learning Interatomic Potential with Time-Dependent Bond Boosting for Polymer and Interface Design" (arXiv:2511.22874, Mori et al.) の手法を、商用利用可能 (commercial-safe) な形で再現実装するリポジトリです。

## 特徴

- **TDBB (Time-Dependent Bond Boosting)** ワークフローの論文忠実な実装 (Eq. 2-8)
- **差し替え可能な MLIP バックエンド**: OrbMol-v2 (推奨), MACE-MP-0, ASE アダプタ, 古典 FF
- **決定論的再現性**: seed, config, git SHA, バックエンド名を RunManifest に記録
- **商用ライセンス・ガードレール**: 依存ライセンスの明示的管理
- **二段構成パイプライン**: 古典 FF (OpenMM/OpenFF) による構造準備 + ML ポテンシャルによる TDBB 本番

## 前提条件

- Python 3.11+
- GPU: CUDA 対応 GPU (OrbMol-v2 使用時)
- Windows 11 + WSL2 (Ubuntu 24.04) の二段構成を推奨

## インストール

### ML 本番環境 (Windows / Linux)

```bash
# conda 環境の作成
conda create -n pfpoly-gpu python=3.12 -c conda-forge -y
conda activate pfpoly-gpu

# PyTorch (CUDA 12.8)
pip install torch==2.10.0 --index-url https://download.pytorch.org/whl/cu128

# プロジェクトのインストール
pip install -e ".[dev]"

# OrbMol-v2 バックエンド (推奨)
pip install -e ".[orb]"

# MACE バックエンド (代替)
pip install -e ".[mace]"

# プロットツール
pip install -e ".[plot]"
```

### 古典構造準備環境 (WSL2 / Linux)

OpenFF/OpenMM は Windows ネイティブでは動作しないため、WSL2 内に別環境を構築します。

```bash
conda create -n pfpoly-prep python=3.12 -c conda-forge -y
conda activate pfpoly-prep
conda install -c conda-forge openff-toolkit-base openff-interchange \
    openff-forcefields openmm rdkit numpy scipy ase -y
pip install -e .
```

## 使い方

### 1. スモークテスト

```bash
# Toy バックエンドでの動作確認
python scripts/run_smoke.py

# ユニットテスト
pytest -q tests/unit

# ライセンスチェック
python scripts/check_dependency_licenses.py --approved specs/approved_dependencies.yaml
```

### 2. ビニルラジカル重合 (基本)

```bash
# 小規模テスト (20 monomer + 1 initiator, CPU)
python scripts/run_vinyl_aibn.py \
    --n-monomers 20 --n-initiators 1 --seed 42 \
    --backend orb --device cpu --no-barostat \
    --n-cycles 5 --biased-steps 2000 --unbiased-steps 500 \
    --output-dir runs/test_small
```

### 3. 論文スケール再現 (二段パイプライン)

```bash
# Step 1: WSL2 で古典構造準備 (OpenFF/Sage, 0.50 g/mL)
wsl -d Ubuntu-24.04 -- bash -lc '
    cd /mnt/c/Users/<user>/Documents/Python/pfpoly
    ~/miniconda3/envs/pfpoly-prep/bin/python scripts/prep_structure.py \
        --n-monomers 100 --n-initiators 5 --seed 42 \
        --charge-method gasteiger --platform CPU \
        --output runs/prep/paper100.json
'

# Step 2: Windows GPU で TDBB 本番
python scripts/run_vinyl_aibn.py \
    --load-structure runs/prep/paper100.json \
    --n-monomers 100 --n-initiators 5 --seed 42 \
    --backend orb --device cuda --no-barostat \
    --n-cycles 30 --biased-steps 2000 --unbiased-steps 500 \
    --equil-steps 2000 --timestep-fs 1.0 \
    --output-dir runs/vinyl_aibn_paper100

# Step 3: 図の生成
python scripts/reproduce_figures.py \
    --trajectory runs/vinyl_aibn_paper100/trajectory.jsonl \
    --bonds runs/vinyl_aibn_paper100/bonds.jsonl \
    --n-reactive-sites 205 --target-temperature 333 \
    --output-dir runs/vinyl_aibn_paper100/figures
```

### 4. AIBN 活性化 + 連鎖重合

```bash
python scripts/run_vinyl_aibn.py \
    --n-monomers 20 --n-initiators 1 --seed 42 \
    --activation --activation-f2 0.3 --activation-f1-max 250 \
    --f2 5.0 --density 0.5 --temperature 333 --no-barostat \
    --backend orb --device cuda \
    --n-cycles 5 --biased-steps 2000 --unbiased-steps 500 \
    --output-dir runs/s4_activation
```

## プロジェクト構成

```
src/kagome/
  boost/          # TDBB ポテンシャル・力・total_bias (Eq. 2-5, 8)
  reactive/       # 反応グループ定義、候補選択、ボンドトラッキング (Eq. 6-7)
  workflows/      # 重合ワークフロー (biased/unbiased 交互ループ)
  backends/       # MLIP/Calculator アダプタ (OrbMol-v2, MACE, ASE, Toy)
  integrators/    # Velocity Verlet, Langevin (BAOAB), MC barostat, FIRE
  analysis/       # 転化率 (Eq. 11), 深さ分解密度, Carothers 解析
  chem/           # RDKit ベースの分子ビルダー
  io/             # JSONL トラジェクトリ I/O
  prep/           # 古典構造準備 (OpenMM/OpenFF Sage)
scripts/          # 実行スクリプト、スキャン、図の生成
configs/          # 実験設定ファイル (YAML)
specs/            # 要件定義、決定記録、ライセンスマトリクス
paper/            # 論文からの構造化ノート
```

## アーキテクチャ

### バックエンドインターフェース

すべての MLIP は `src/kagome/backends/base.py` の `Calculator` プロトコルに準拠します。

```python
class Calculator(Protocol):
    def compute(self, positions, species, cell=None) -> tuple[float, NDArray]:
        '''エネルギー [kcal/mol] と力 [kcal/(mol*Angstrom)] を返す。'''
        ...
```

### ワークフロー

`PolymerizationWorkflow.run()` が以下のサイクルを繰り返します:

1. **候補選択**: 反応グループから候補ペアを列挙・スコアリング (Eq. 7)
2. **Biased phase**: TDBB バイアス付き MD (Eq. 2-5)
3. **Unbiased phase**: バイアスなし MD で緩和
4. **結合判定**: 距離閾値で結合形成/解離を確認
5. **グループ更新**: 確認された反応に基づき反応グループを更新

### 単位系

LAMMPS 'real' スタイル: エネルギー=kcal/mol, 距離=Angstrom, 時間=fs, 質量=amu。
定数は `src/kagome/units.py` で一元管理。

### 実験トレーサビリティ

すべての実行は `RunManifest` (seed, config, git SHA, backend, dirty flag) を出力し、
トラジェクトリ (JSONL) と結合イベントログ (JSONL) から図を再生成可能です。

## 現在の進捗

| マイルストーン | 状態 |
|---|---|
| S1: 単鎖伝播デモ | 完了 (pentamer 構築, radical 移動確認) |
| S2: メルト駆動での結合形成 | 完了 (OrbMol-v2, formations > 0) |
| S3: 多ラジカル・高スピン近似 | 完了 (2-radical, alpha=27.3%) |
| S4: AIBN 活性化 + 連鎖重合 | 完了 (C-N 解離 + 結合形成実証) |
| S5: 図の再現・比較 | 一部完了 (energy, conversion, temperature) |
| S6: 論文スケール実行 (200+10) | 未着手 (24 GB+ GPU が必要) |

## ドキュメント

| ファイル | 内容 |
|---|---|
| `CLAUDE.md` | 開発ルール・非交渉要件 |
| `specs/decisions.md` | 設計判断と根拠の記録 |
| `specs/requirements.md` | 機能・非機能要件 |
| `specs/tasks.md` | タスクバックログ |
| `specs/figure-comparison.md` | 論文図 vs 再現図の比較 |
| `specs/dependency-license-matrix.md` | 依存ライセンスの管理 |
| `specs/handoff-plan-v5.md` | コードレビュー是正計画 (RF1-RF13) |
| `specs/handoff-plan-v6.md` | レビュー是正計画 (RF14-RF24) |
| `specs/s6-environment-setup.md` | S6 実行環境セットアップ |

## ライセンス

本プロジェクトの依存ライセンス状況は `specs/dependency-license-matrix.md` を参照してください。
商用利用に関するガードレールは `CLAUDE.md` に定義されています。
