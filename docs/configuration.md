# コンフィグリファレンス

KAGOME のシミュレーションパラメータの一覧です。パラメータは CLI 引数で指定し、実行時に `manifest.json` に記録されます。

---

## TDBB パラメータ

Time-Dependent Bond Boosting (Eq. 2-5, 8) の制御パラメータ。

| パラメータ | CLI 引数 | デフォルト | 論文値 | 単位 | 説明 |
|---|---|---|---|---|---|
| f2 | `--f2` | 10.0 | 10.0 | A^-2 | ガウシアン幅パラメータ (Eq. 2-3)。小さいほど広いバイアス井戸 |
| gamma | - | 1.0 | 1.0 | kcal/(mol*step) | f1 の時間増加率 (Eq. 5) |
| f1_max (形成) | - | 250.0 | 250.0 | kcal/mol | V^f の最大振幅 |
| f1_max (解離) | - | 125.0 | 125.0 | kcal/mol | V^d の最大振幅 |
| lambda_vdw | - | 0.60 | 0.60 | - | r0 の vdW 半径スケーリング (Eq. 4) |

### OrbMol-v2 向け調整値

OrbMol-v2 のエネルギー障壁は PFP と異なるため、以下の調整が有効です:

| パラメータ | 調整値 | 理由 |
|---|---|---|
| `--f2` | 5.0 | バイアス井戸を広げ、捕獲半径を拡大 |
| `--select-rmin` | 1.5 | 候補選択の下限を縮小 |
| `--select-rmax` | 3.0 | 候補選択の上限を縮小 |

---

## AIBN 活性化パラメータ

AIBN の C-N アゾ結合を V^d バイアスで解離するためのパラメータ。

| パラメータ | CLI 引数 | デフォルト | 説明 |
|---|---|---|---|
| 有効化 | `--activation` | False | AIBN 活性化フェーズを有効にする |
| f2 (活性化) | `--activation-f2` | 0.3 | V^d のガウシアン幅 (OrbMol-v2 の C-N 障壁に合わせて調整) |
| f1_max (活性化) | `--activation-f1-max` | 250.0 | V^d の最大振幅 |
| 最大ステップ | `--activation-steps` | 3000 | 活性化フェーズの最大ステップ数 |

---

## MD パラメータ

| パラメータ | CLI 引数 | デフォルト | 論文値 | 単位 | 説明 |
|---|---|---|---|---|---|
| タイムステップ | `--timestep-fs` | 0.25 | 0.25 | fs | MD 積分のタイムステップ |
| biased ステップ | `--biased-steps` | 500 | 2000 | steps | biased phase のステップ数 |
| unbiased ステップ | `--unbiased-steps` | 500 | 2000 | steps | unbiased phase のステップ数 |
| サイクル数 | `--n-cycles` | 3 | - | - | biased/unbiased の繰り返し回数 |
| 平衡化ステップ | `--equil-steps` | 2000 | - | steps | TDBB 前のバイアスなし平衡化 |
| 温度 | `--temperature` | 333.0 | 333 | K | NVT/NPT の目標温度 |
| 圧力 | `--pressure` | 1.0 | 1.0 | atm | NPT の目標圧力 |
| barostat 無効化 | `--no-barostat` | False | - | - | NVT で実行 (barostat を無効化) |

> **注意**: デフォルトの `--biased-steps 500` と `--unbiased-steps 500` は高速テスト用です。論文準拠の実行では `2000` を使用してください。

---

## 系の構成パラメータ

| パラメータ | CLI 引数 | デフォルト | 説明 |
|---|---|---|---|
| モノマー数 | `--n-monomers` | 8 | メチルアクリレートモノマー数 |
| 開始剤数 | `--n-initiators` | 2 | 開始剤 (AIBN) の数 |
| 密度 | `--density` | 0.5 | 初期密度 [g/mL]。論文 SI S-3 は 0.5 |
| ボックスサイズ | `--box-size` | None | 直接指定 [A]。省略時は密度から計算 |
| スピン多重度 | `--spin` | 1 | スピン多重度 (2S+1)。ラジカル 1 個なら 2 |
| 開始剤 SMILES | `--initiator-smiles` | None | 開始剤の SMILES を直接指定 |

---

## 構造準備パラメータ

| パラメータ | CLI 引数 | デフォルト | 説明 |
|---|---|---|---|
| 前処理済み構造 | `--load-structure` | None | PreparedStructure JSON を読み込む |
| 圧縮バックエンド | `--compress-backend` | `classical` | `classical` (OpenMM) or `ml` |
| 圧縮プラットフォーム | `--compress-platform` | `CPU` | OpenMM プラットフォーム |
| FIRE 最小化 | `--minimize` | True | TDBB 前に FIRE 最小化を実行 |
| FIRE 閾値 | `--minimize-fmax` | 1.0 | FIRE 収束閾値 [kcal/(mol*A)] |

### prep_structure.py 固有パラメータ

| パラメータ | CLI 引数 | デフォルト | 説明 |
|---|---|---|---|
| 目標密度 | `--target-density` | 0.5 | 目標密度 [g/mL] |
| 電荷法 | `--charge-method` | `gasteiger` | `gasteiger` (軽量) or `nagl` (高精度) |
| 力場 | `--forcefield` | `openff-2.2.0.offxml` | OpenFF Sage 力場 |
| 圧縮ステージ | `--compress-stages` | 20 | 幾何圧縮の段階数 |
| NVT ステップ | `--nvt-steps` | 50000 | 熱平衡化ステップ数 |
| プラットフォーム | `--platform` | `CPU` | OpenMM プラットフォーム |

---

## バックエンドパラメータ

| パラメータ | CLI 引数 | デフォルト | 説明 |
|---|---|---|---|
| バックエンド | `--backend` | `orb` | `orb` (OrbMol-v2) or `mace` |
| デバイス | `--device` | `cpu` | `cpu` or `cuda` |
| MACE モデル | `--model` | `small` | MACE モデルサイズ (`--backend mace` 時のみ) |

---

## 候補選択パラメータ

| パラメータ | CLI 引数 | デフォルト | 論文値 | 説明 |
|---|---|---|---|---|
| 最小距離 | `--select-rmin` | None (3.0) | 3.0 A | 候補ペアの最小距離 |
| 最大距離 | `--select-rmax` | None (6.0) | 6.0 A | 候補ペアの最大距離 |

> `None` の場合、テンプレートのデフォルト値 (3.0, 6.0) が使用されます。

---

## 出力制御パラメータ

| パラメータ | CLI 引数 | デフォルト | 説明 |
|---|---|---|---|
| 出力ディレクトリ | `--output-dir` | `runs/vinyl_aibn` | 全出力ファイルの保存先 |
| シード | `--seed` | 7 | 乱数シード (再現性のため) |

---

## YAML コンフィグファイル

`configs/` 以下に論文準拠のパラメータセットが用意されています。CLI 引数による指定が優先されます。

### configs/boost/paper_faithful.yaml

```yaml
boost:
  mode: paper_faithful
  timestep_fs: 0.25          # PDF p.7
  biased_steps: 2000         # PDF p.7
  unbiased_steps: 2000       # PDF p.7
  lambda_vdw: 0.60           # Eq. 4
  f2_invA2: 10.0             # PDF p.7
  gamma: 1.0                 # PDF p.7
  f1_max_form_kcal_mol: 250.0
  f1_max_break_kcal_mol: 125.0
```

> このファイルは **文書化目的** です。実行時のパラメータは CLI 引数で指定し、`manifest.json` に記録されます。

---

## Langevin 積分器パラメータ (コード内)

CLI では直接設定できませんが、コード内で以下のデフォルト値が使用されます:

| パラメータ | デフォルト | 論文値 | 説明 |
|---|---|---|---|
| temperature_K | 300.0 | 333 | スクリプト側で上書き |
| friction_per_fs | 0.001 | 0.001 (= 1.0 ps^-1) | Langevin 摩擦係数 |

---

## MC Barostat パラメータ (コード内)

| パラメータ | デフォルト | 説明 |
|---|---|---|
| pressure_atm | 1.0 | 目標圧力 |
| frequency | 25 | バロスタット試行頻度 (ステップ) |
| max_volume_change_frac | 0.01 | 最大体積変化率 |
