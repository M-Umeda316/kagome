# クイックスタート

このガイドでは、KAGOME を使ったシミュレーションを段階的に実行する方法を説明します。

> 環境構築がまだの場合は、先に [インストールガイド](installation.md) を参照してください。

---

## Step 1: スモークテスト

Toy バックエンド (LJ ポテンシャル) で TDBB ワークフローの動作を確認します。MLIP は不要です。

```bash
python scripts/run_smoke.py
```

出力先: `runs/smoke/`

| ファイル | 内容 |
|---|---|
| `trajectory.jsonl` | 各ステップの位置・エネルギー・温度 |
| `bonds.jsonl` | 結合イベント (試行/確認) |
| `manifest.json` | 再現情報 (seed, config, git SHA, backend) |
| `summary.json` | 実行サマリー |

---

## Step 2: 小規模ビニル重合 (CPU)

OrbMol-v2 バックエンドで実際のビニルラジカル重合を実行します。

```bash
python scripts/run_vinyl_aibn.py \
    --n-monomers 20 --n-initiators 1 --seed 42 \
    --backend orb --device cpu --no-barostat \
    --n-cycles 5 --biased-steps 2000 --unbiased-steps 500 \
    --output-dir runs/test_small
```

主要パラメータ:

| パラメータ | 値 | 意味 |
|---|---|---|
| `--n-monomers` | 20 | メチルアクリレートモノマー数 |
| `--n-initiators` | 1 | 開始剤 (イソブチロニトリル) 数 |
| `--backend orb` | - | OrbMol-v2 バックエンド |
| `--no-barostat` | - | NVT アンサンブル (NPT を無効化) |
| `--n-cycles` | 5 | biased/unbiased サイクル数 |
| `--biased-steps` | 2000 | biased phase ステップ数 (論文準拠) |

> CPU 実行は低速ですが、GPU なしで動作確認できます。

---

## Step 3: GPU で高速実行

GPU が使える場合は `--device cuda` を指定します。

```bash
python scripts/run_vinyl_aibn.py \
    --n-monomers 20 --n-initiators 1 --seed 42 \
    --backend orb --device cuda --no-barostat \
    --n-cycles 15 --biased-steps 2000 --unbiased-steps 500 \
    --output-dir runs/test_gpu
```

---

## Step 4: AIBN 活性化 + 連鎖重合

論文の Table S1 に記載された AIBN の C-N 結合解離 (V^d) から開始する完全なワークフロー:

```bash
python scripts/run_vinyl_aibn.py \
    --n-monomers 20 --n-initiators 1 --seed 42 \
    --activation --activation-f2 0.3 --activation-f1-max 250 \
    --f2 5.0 --density 0.5 --temperature 333 --no-barostat \
    --backend orb --device cuda \
    --n-cycles 5 --biased-steps 2000 --unbiased-steps 500 \
    --output-dir runs/s4_activation
```

`--activation` フラグにより:
1. 完全な AIBN 分子を構築 (事前切断しない)
2. V^d バイアスで C-N アゾ結合を解離 (活性化フェーズ)
3. ラジカル生成後、V^f バイアスで連鎖重合を実行

---

## Step 5: 論文スケール再現 (二段パイプライン)

100+ モノマーの論文スケールでは、古典力場による構造準備を先に行います。

### 5a. 古典構造準備 (WSL2)

```bash
wsl -d Ubuntu-24.04 -- bash -lc '
    cd /mnt/c/Users/<user>/Documents/Python/pfpoly
    ~/miniconda3/envs/pfpoly-prep/bin/python scripts/prep_structure.py \
        --n-monomers 100 --n-initiators 5 --seed 42 \
        --charge-method gasteiger --platform CPU \
        --output runs/prep/paper100.json
'
```

OpenFF Sage 2.2 力場で:
- 希薄配置 → 0.5 g/mL まで圧縮
- NVT 熱平衡化 (333 K, 50000 ステップ)

### 5b. TDBB 本番 (Windows GPU)

```bash
python scripts/run_vinyl_aibn.py \
    --load-structure runs/prep/paper100.json \
    --n-monomers 100 --n-initiators 5 --seed 42 \
    --backend orb --device cuda --no-barostat \
    --n-cycles 30 --biased-steps 2000 --unbiased-steps 500 \
    --equil-steps 2000 --timestep-fs 1.0 \
    --output-dir runs/vinyl_aibn_paper100
```

`--load-structure` で前処理済み構造を読み込み、ML ポテンシャルで TDBB サイクルを実行します。

### 5c. 図の生成

```bash
python scripts/reproduce_figures.py \
    --trajectory runs/vinyl_aibn_paper100/trajectory.jsonl \
    --bonds runs/vinyl_aibn_paper100/bonds.jsonl \
    --n-reactive-sites 205 --target-temperature 333 \
    --output-dir runs/vinyl_aibn_paper100/figures
```

生成される図:
- `energy_vs_step.png` -- biased/unbiased エネルギー推移
- `conversion_vs_step.png` -- モノマー転化率 alpha(t) + Eq. 11 フィット
- `temperature_vs_step.png` -- 温度安定性の確認

---

## 出力ファイルの読み方

すべての実行は以下の標準出力を生成します:

### trajectory.jsonl

JSONL 形式のトラジェクトリ。1 行目はヘッダ (スキーマ, 元素種, 原子数)、以降は各フレーム。

```python
from kagome.io.readers import read_trajectory
header, frames = read_trajectory('runs/test_small/trajectory.jsonl')
print(f'原子数: {header["n_atoms"]}, フレーム数: {len(frames)}')
```

### bonds.jsonl

結合イベントログ。各行は 1 つの結合イベント (試行/確認/拒否)。

```python
from kagome.io.readers import read_bond_events
events = read_bond_events('runs/test_small/bonds.jsonl')
confirmed = [e for e in events if e.event_type == 'confirmed_formation']
print(f'確認された結合形成: {len(confirmed)}')
```

### manifest.json

実験の再現に必要なすべてのメタデータ (seed, config, git SHA, backend, model_id)。

---

## 次のステップ

- [バックエンドの選択](backends.md) -- OrbMol-v2 / MACE の使い分け
- [コンフィグリファレンス](configuration.md) -- 全パラメータの詳細
- [論文再現ガイド](paper-reproduction.md) -- S1-S6 マイルストーンの詳細
- [アーキテクチャ](architecture.md) -- コードの構造と拡張方法
