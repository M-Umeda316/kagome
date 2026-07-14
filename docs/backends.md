# バックエンド

KAGOME では、MLIP (Machine-Learned Interatomic Potential) を差し替え可能なバックエンドとして管理しています。すべてのバックエンドは `Calculator` プロトコルに準拠し、エネルギー [kcal/mol] と力 [kcal/(mol*A)] を返します。

---

## 比較表

| | OrbMol-v2 | MACE-MP-0 | Classical (OpenFF) | Toy |
|---|---|---|---|---|
| **推奨度** | **推奨** | 代替 | 構造準備用 | テスト用 |
| **ライセンス** | Apache-2.0 | MIT | MIT/LGPL/CC-BY-4.0 | 内部コード |
| **商用利用** | OK | OK | OK (帰属表示必要) | - |
| **有機系精度** | 高 | 中 | 低 (力場依存) | なし |
| **スピン対応** | あり | なし | なし | なし |
| **GPU 対応** | あり | あり | あり (OpenMM) | なし |
| **周期境界** | あり (要 nvalchemiops) | あり | あり | あり |
| **用途** | 本番シミュレーション | 汎用/検証 | 圧縮・前処理 | 単体テスト |

---

## OrbMol-v2 (推奨)

Apache-2.0 ライセンスの汎用 MLIP。コードと重みの両方が商用利用可能です。有機分子系に対して高い精度を持ち、スピン多重度のサポートによりラジカル系のシミュレーションに適しています。

### インストール

```bash
pip install -e ".[orb]"
```

### 使い方

```bash
python scripts/run_vinyl_aibn.py --backend orb --device cuda
```

### パラメータ

| パラメータ | デフォルト | 説明 |
|---|---|---|
| `model` | `orbmol_v2` | モデル名 |
| `device` | `cpu` | `cpu` or `cuda` |
| `compile` | `False` | torch.compile の使用 (Linux のみ)。run スクリプトの `--compile` から指定可 |
| `charge` | `0` | 分子の全電荷 |
| `spin` | `1` | スピン多重度 (2S+1) |
| `empty_cache` | `True` | 毎ステップ `torch.cuda.empty_cache()` を呼ぶ。32 GB 以上の GPU では `--no-empty-cache` で無効化可 (CPU 時間 ~9% 削減、decisions.md 2026-07-14) |

### 注意事項

- **VRAM**: 2520 原子で約 9.5 GB。24 GB 以上の GPU を推奨
- **Windows**: `nvalchemiops` は `torch.compile` が必要なため PME 静電が利用不可。非周期計算は問題なし
- **ローカルチェックポイント**: `models/orbmol-v2-teqabfhg-20260523.ckpt` があれば優先使用。なければ自動ダウンロード
- **CUDA メモリ**: 近傍グラフサイズの変動により `empty_cache()` を毎ステップ呼び出し (16 GB GPU での断片化対策、decisions.md 2026-06-15)。VRAM に余裕がある GPU (32 GB 以上) では `empty_cache=False` (`--no-empty-cache`) で無効化でき、CPU 時間を ~9% 削減 (py-spy 実測、decisions.md 2026-07-14)
- **torch.compile**: `compile=True` (run スクリプトの `--compile`) で有効。upstream 公称 ~1.7x。従来は `TORCHDYNAMO_DISABLE=1` が無条件設定されていたため無効化されていたが、compile=False の時のみ設定するよう修正済み (decisions.md 2026-07-14)。採用前に compile 有無での数値等価性検証が必要

---

## MACE-MP-0

MIT ライセンスの汎用 MLIP。コードも重みも商用利用可能です。元素周期表の広い範囲をカバーしますが、有機ポリマー系の精度は OrbMol-v2 に劣ります。

### インストール

```bash
pip install -e ".[mace]"
```

### 使い方

```bash
python scripts/run_vinyl_aibn.py --backend mace --model small --device cuda
```

### パラメータ

| パラメータ | デフォルト | 説明 |
|---|---|---|
| `model` | `small` | `small`, `medium`, `large` |
| `device` | `cpu` | `cpu` or `cuda` |
| `default_dtype` | `float64` | 計算精度 |

### 注意事項

- **MACE-OFF はブロック**: ASL ライセンスのため商用利用不可。`mace_off` / `mace-off` / `off23` を指定するとエラー
- **スピン非対応**: ラジカル系では closed-shell 近似になる
- **ローカルチェックポイント**: `models/mace-mp-0-small.model` があれば優先使用

---

## Classical (OpenFF/OpenMM)

OpenMM + OpenFF Sage 2.2 力場による古典計算。トポロジー情報 (結合・角度・二面角) を使うため、分子構造の定義が必要です。主に構造準備 (密度圧縮) に使用します。

### 用途

TDBB 本番ではなく、以下の前処理で使用:
- 希薄配置から目標密度への圧縮
- FIRE 最小化の加速
- GPU VRAM を ML バックエンドに温存

### 使い方

`run_vinyl_aibn.py` では `--compress-backend classical` (デフォルト) で自動使用されます。

```bash
# 構造準備スクリプトで直接使用
python scripts/prep_structure.py \
    --n-monomers 100 --n-initiators 5 \
    --charge-method gasteiger --platform CPU \
    --output runs/prep/structure.json
```

### パラメータ

| パラメータ | デフォルト | 説明 |
|---|---|---|
| `charge_method` | `gasteiger` | `gasteiger` (軽量) or `nagl` (GNN AM1-BCC 代替) |
| `forcefield` | `openff-2.2.0.offxml` | OpenFF Sage 力場 |
| `platform` | `CPU` | OpenMM プラットフォーム |
| `cutoff_nm` | `0.8` | 非結合カットオフ (ボックス辺に応じて自動調整) |

---

## Toy

単純なペアワイズ LJ ポテンシャル。実際の化学的意味はなく、TDBB ワークフローのロジックテストにのみ使用します。

```bash
python scripts/run_smoke.py
```

| パラメータ | デフォルト | 説明 |
|---|---|---|
| `epsilon` | `0.1` | LJ 深さ (kcal/mol) |
| `sigma` | `1.5` | LJ サイズパラメータ (A) |

---

## バックエンドの選び方

### 論文再現が目的

→ **OrbMol-v2** を推奨。有機系の精度が高く、スピン対応あり。論文は PFP (Matlantis) を使用していますが、商用ライセンスの関係で OrbMol-v2 が最も近い代替です。

### 汎用・探索的な使用

→ **MACE-MP-0** が選択肢。広い元素カバレッジで、有機系以外 (金属酸化物等) にも対応。精度は OrbMol-v2 に劣るが、MIT ライセンスで制約なし。

### 構造準備のみ

→ **Classical (OpenFF)** で十分。ML ポテンシャルの GPU 時間を本番に温存できます。

### CI/テスト

→ **Toy** バックエンド。外部依存なしで TDBB ロジックのテストが可能。

---

## 新しいバックエンドの追加

1. `src/kagome/backends/` に新ファイルを作成
2. `Calculator` ABC を継承し、`compute` と `name` を実装
3. エネルギーは **kcal/mol**、力は **kcal/(mol*A)** で返すこと
4. `model_id` プロパティでチェックポイント/重みの識別子を返す
5. ライセンスを `specs/dependency-license-matrix.md` に記載
6. テストを `tests/unit/` に追加

```python
from kagome.backends.base import Calculator

class MyBackend(Calculator):
    def compute(self, positions, species, cell=None):
        # ... 計算 ...
        return energy_kcal_mol, forces_kcal_mol_per_A

    @property
    def name(self):
        return 'my-backend'

    @property
    def model_id(self):
        return f'my-backend:{self._checkpoint_hash}'
```
