# アーキテクチャ

KAGOME は、論文 (arXiv:2511.22874) の TDBB (Time-Dependent Bond Boosting) ワークフローを再現する反応分子動力学シミュレーションパッケージです。

---

## モジュール構成

```
src/kagome/
  boost/          TDBB バイアスポテンシャル・力 (Eq. 2-5, 8)
  reactive/       反応グループ・候補選択・結合トラッキング (Eq. 6-7)
  workflows/      重合ワークフロー・実験プロベナンス (Fig. 1)
  backends/       MLIP/Calculator アダプタ
  integrators/    Velocity Verlet, Langevin (BAOAB), MC barostat, FIRE
  analysis/       転化率, 密度プロファイル, Carothers 解析
  chem/           RDKit ベースの分子ビルダー
  io/             JSONL トラジェクトリ I/O
  prep/           古典構造準備 (OpenMM/OpenFF)
  units.py        単位定数の一元管理
  geometry.py     PBC 処理 (minimum image, wrap)
```

---

## データフロー

```
                  ┌─────────────────┐
                  │  分子ビルダー     │  chem/builders.py
                  │  SMILES → 3D 構造 │  scripts/_systems.py
                  └────────┬────────┘
                           │ positions, species, cell
                           │ template, groups
                           ▼
              ┌────────────────────────┐
              │  構造準備 (任意)         │  prep/openmm_equilibrate.py
              │  古典FF: 圧縮 + 熱平衡化 │  integrators/minimize.py
              └────────────┬───────────┘
                           │ PreparedStructure (.json)
                           ▼
    ┌──────────────────────────────────────────┐
    │         TDBB ワークフロー                  │  workflows/polymerization.py
    │                                          │
    │  ┌─ [FIRE 最小化] ─ [NPT 平衡化] ─┐      │
    │  │                                │      │
    │  ▼                                │      │
    │  ┌──────────────────────────┐     │      │
    │  │ Biased Phase (2000 steps)│◄────┘      │
    │  │  候補選択 (Eq. 7)         │            │
    │  │  TDBB バイアス (Eq. 2-5)  │            │
    │  │  MD 積分 + 結合チェック    │            │
    │  └───────────┬──────────────┘            │
    │              │                           │
    │              ▼                           │
    │  ┌──────────────────────────┐            │
    │  │ Unbiased Phase (2000 st) │            │
    │  │  バイアスなし MD           │            │
    │  │  結合確認/拒否            │            │
    │  └───────────┬──────────────┘            │
    │              │                           │
    │              ▼                           │
    │  ┌──────────────────────────┐            │
    │  │ グループ更新               │            │
    │  │  確認された反応に基づき     │            │
    │  │  反応グループを変更        │────► 次サイクル
    │  └──────────────────────────┘            │
    └──────────────────────────────────────────┘
                           │
                           │ trajectory.jsonl, bonds.jsonl
                           │ manifest.json, selection.jsonl
                           ▼
              ┌────────────────────────┐
              │  解析 + 図の生成        │  analysis/
              │  転化率, kp_eff フィット │  scripts/reproduce_figures.py
              │  密度プロファイル       │
              └────────────────────────┘
```

---

## 主要インターフェース

### Calculator プロトコル

すべての MLIP バックエンドは `Calculator` ABC に準拠します。

```python
class Calculator(ABC):
    @abstractmethod
    def compute(self, positions, species, cell=None) -> tuple[float, NDArray]:
        '''エネルギー [kcal/mol] と力 [kcal/(mol*A)] を返す。'''

    @property
    @abstractmethod
    def name(self) -> str: ...

    @property
    def model_id(self) -> str:
        '''プロベナンス用の重み/チェックポイント識別子。'''
        return self.name

    @property
    def supports_spin(self) -> bool:
        return False
```

新しいバックエンドを追加する場合は、この ABC を継承して `compute` と `name` を実装します。エネルギーと力は **kcal/mol / kcal/(mol*A)** で返す必要があります。

### Integrator プロトコル

```python
class Integrator(Protocol):
    def pre_force(self, positions, velocities, forces, masses, dt, rng, cell): ...
    def post_force(self, velocities, forces, masses, dt): ...
```

MD ステップは `pre_force → compute → [bias] → post_force → [barostat]` の順で実行されます。

### ReactionTemplate

反応系を定義するテンプレート。反応グループ (Eq. 6) とペア仕様 (距離窓、形成/解離) を保持します。

```python
template = ReactionTemplate(
    name='vinyl_radical',
    groups=['radical_C', 'vinyl_alpha_C', 'chain_C', 'vinyl_beta_C'],
    pairs=[
        PairSpec('radical_C', 'vinyl_alpha_C', is_formation=True, r_min=3.0, r_max=6.0),
        PairSpec('radical_C', 'chain_C', is_formation=True, constraint_only=True),
        PairSpec('vinyl_alpha_C', 'vinyl_beta_C', is_formation=True, constraint_only=True),
    ]
)
```

---

## 単位系

LAMMPS 'real' スタイルを全モジュールで統一使用:

| 物理量 | 単位 | 定数名 |
|---|---|---|
| エネルギー | kcal/mol | - |
| 距離 | Angstrom (A) | - |
| 時間 | femtosecond (fs) | - |
| 質量 | amu (g/mol) | - |
| 温度 | Kelvin (K) | - |
| ボルツマン定数 | kcal/(mol*K) | `KB = 0.001987204` |
| 力→加速度 | A/fs^2 | `FORCE_CONV = 4.184e-4` |
| 圧力 | kcal/(mol*A^3) | `ATM_TO_KCAL_MOL_A3 = 1.4584e-5` |
| eV → kcal/mol | - | `EV_TO_KCAL_MOL = 23.060548` |

### 単位変換の境界

バックエンドと内部コードの間で単位変換が発生します:

| 境界 | 変換元 | 変換先 |
|---|---|---|
| ASE/MACE | eV, eV/A | kcal/mol, kcal/(mol*A) |
| OrbMol-v2 | eV, eV/A | kcal/mol, kcal/(mol*A) |
| OpenMM (古典) | kJ/mol, kJ/(mol*nm) | kcal/mol, kcal/(mol*A) |

内部コード間での単位変換は不要です。すべてのモジュールは LAMMPS real 単位で統一されています。

---

## 論文方程式とコードの対応

| 方程式 | 内容 | 実装 |
|---|---|---|
| Eq. 2 | V^f: 結合形成バイアス | `boost/tdbb.py::formation_potential` |
| Eq. 3 | V^d: 結合解離バイアス | `boost/tdbb.py::dissociation_potential` |
| Eq. 4 | r0: vdW 半径からの目標距離 | `boost/tdbb.py::target_distance` |
| Eq. 5 | f1(t): 時間依存振幅 | `boost/tdbb.py::boost_amplitude` |
| Eq. 6 | G_X: 反応原子グループ | `reactive/groups.py::ReactiveGroup` |
| Eq. 7 | d_ijkl: 候補スコアリング | `reactive/selection.py::find_candidates` |
| Eq. 8 | 全バイアスポテンシャル | `boost/tdbb.py::total_bias` |
| Eq. 11 | alpha(t) 指数フィット | `analysis/conversion.py::fit_conversion_exponential` |
| PDF p.12 | rho_rxn(z) 深さ密度 | `analysis/density.py::reaction_density_profile` |

---

## 積分器

| 積分器 | ファイル | 用途 |
|---|---|---|
| Velocity Verlet | `integrators/verlet.py` | 基本 NVE 積分 |
| Langevin BAOAB | `integrators/langevin.py` | NVT 温度制御 (論文 Section 2) |
| MC Barostat | `integrators/mc_barostat.py` | NPT 圧力制御 |
| FIRE | `integrators/minimize.py` | エネルギー最小化 + ボックス圧縮 |
| Maxwell-Boltzmann | `integrators/init_velocities.py` | 初期速度生成 |

---

## ポストサイクル更新戦略

結合形成が確認された後の反応グループ更新は、反応タイプに応じて異なります:

| 戦略 | クラス | 用途 |
|---|---|---|
| デフォルト | `DefaultPostCycleUpdater` | ステップ成長 (ナイロン): 反応した原子をグループから除去 |
| ビニル連鎖伝播 | `VinylChainPropagationUpdater` | ラジカル重合: beta-C が新ラジカルに、連鎖が伸長 |

---

## PBC (周期境界条件)

現在は直方体 (orthorhombic) セルのみをサポートしています。三斜晶系セルは `ValueError` を送出します。

- `minimum_image(r_vec, cell)` -- 最小イメージ規約による変位ベクトル
- `wrap_positions(positions, cell)` -- 位置を一次セル内に折り返し

---

## 実験トレーサビリティ

すべての実行は `RunManifest` を出力し、以下を記録します:

- **seed**: 乱数シード
- **config_path**: 使用した設定ファイル
- **git_sha**: コミットハッシュ (自動取得)
- **git_dirty**: 未コミットの変更があるか
- **backend**: バックエンド名と `model_id`
- **timestamp**: UTC ISO 形式のタイムスタンプ
- **extra**: 実効パラメータ (biased_steps, f2 等)

`trajectory.jsonl` と `bonds.jsonl` から図を再生成でき、`manifest.json` から実行条件を完全に復元できます。
