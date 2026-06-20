# 作業計画書 v3 — 化学系拡充・OrbMol-v2デフォルト化・定性比較検証

最終更新: 2026-06-13
前提文書: `specs/handoff-plan-v2.md`（T-A〜T-G1 完了）、`CLAUDE.md`（絶対遵守ルール）
論文: `paper/2511.22874v1.pdf`（28ページ, 本文+Supporting Information）

---

## 0. PDFから確定した事実（v2で未確定だった項目の解消）

### P-1: γ の単位 → **kcal/(mol·step) のまま維持**

PDFの記述（p.7, Section 3）:

> "alternating biased and unbiased dynamics every 2000 steps (500 fs), and bias parameters of f_1^max = 250 kcal/mol for V^f, f_1^max = 125 kcal/mol for V^d, f2 = 10 Å^{-2}, and γ = 1.0"

γ の単位は明示されていない。Eq. 5 の `f1(t) = γt` の `t` の単位も明記なし。

**判断**: 現行コードの kcal/(mol·step) を維持する。理由:
- PDF にも単位は不明。変更する根拠がない。
- 現行実装で saturation at step 250 (62.5 fs) ← biased phase 2000 steps (500 fs) の 12.5% で飽和。
- もし t が物理時間 fs なら saturation at 250 fs = 1000 steps ← biased phase の 50% で飽和。
- どちらでも biased phase 内で飽和するため、定性動作は同一。定量差は γ スケーリングで吸収される（Fig. S4 参照）。
- `decisions.md` の PENDING CONFIRMATION を「PDF confirmed: unit not stated, maintaining kcal/(mol·step)」に更新する。

### P-2: unbiased ステップ数 → **2000 で確定**

PDFの記述（p.7）:

> "alternating biased and **unbiased** dynamics every **2000 steps** (500 fs)"

biased と unbiased の**両方**が 2000 steps であることが PDF で確認された。

**アクション**: `paper/notes.md` と `configs/boost/paper_faithful.yaml` のコメントを確定形に更新。

### P-3: NPT ensemble → 1 atm, Langevin + MC barostat で確定

PDF Supporting Information（p.20）:

> "Equilibration simulations were performed in the NPT ensemble at 300 K and **1 atm**. Temperature was controlled using a **Langevin thermostat** with a coupling constant of **1.0 ps^{-1}**, and pressure was maintained using a **Monte Carlo barostat**."
> "Production simulations using reactive acceleration MD were then carried out at **333 K and 1 atm** in the NPT ensemble."

**確定事項**:
- 目標圧力: **1 atm** ← 現行コードの `MCBarostatParams(pressure_atm=1.0)` と一致 ✅
- バロスタット種類: **Monte Carlo** ← 現行の `MCBarostat` と一致 ✅
- サーモスタット: **Langevin, coupling = 1.0 ps^{-1}** ← 現行の `friction_per_fs=0.01` = 0.01/fs = 10/ps ≠ 1.0/ps

**⚠️ 仮説**: 現行の `friction_per_fs=0.01` (= 10 ps^{-1}) は論文の 1.0 ps^{-1} の **10倍**。
`decisions.md` に記録し、`friction_per_fs=0.001` (= 1.0 ps^{-1}) をデフォルトに変更するかオーナー判断を仰ぐ。

### P-4: 生産シミュレーション温度 → **333 K**（vinyl系）、300 K（nylon/epoxy）

PDF（p.21）: "Production simulations ... at **333 K** and 1 atm" (vinyl radical polymerization)
PDF（p.22）: "production ... at **300 K** and 1 atm" (nylon-6,6)
PDF（p.24）: "production ... in the NVT ensemble at **333 K**" (epoxy, NVTで実施)

**⚠️ 注意**: 現行スクリプトは `--temperature 500.0` をデフォルトにしている → 論文値と異なる。修正が必要。

### P-5: 論文の反応パターンテーブル（Table S1-S3）→ 候補選択パラメータ確定

**Table S1 (Vinyl radical polymerization)**:
| Reaction | i-j rmin/rmax | i-k rmin/rmax | j-l rmin/rmax | V^f | V^d |
|---|---|---|---|---|---|
| Activation | 0.0 / 3.0 | 0.0 / 3.0 | 0.0 / 3.0 | — | i-k, j-l |
| Initiation | 3.0 / 6.0 | 0.0 / 3.0 | 0.0 / 3.0 | i-j | — |
| Propagation | 3.0 / 6.0 | 0.0 / 3.0 | 0.0 / 3.0 | i-j | — |

**⚠️ 重要発見**: 論文のビニル重合は **4グループ** テンプレート (i, j, k, l) を使っている。
- i-j: ラジカルC と ビニルα-C (3.0-6.0 Å, formation)
- i-k: ラジカルC と 既存結合C (0.0-3.0 Å, confirmation) ← **現行コードにない**
- j-l: ビニルα-C と ビニルβ-C (0.0-3.0 Å, confirmation) ← **現行コードにない**

これは**結合判定の追加条件**であり、単に i-j 距離だけでなく、k と l の存在確認（ラジカルが既存鎖末端にいること、ビニル二重結合が存在すること）も行っている。

**しかし**: 現行の2グループ実装は**機能的に等価**な近似と見なせる:
- Propagation で i-k (0.0-3.0) は「iがすでに鎖末端にいる（ラジカル）」ことの確認
- j-l (0.0-3.0) は「jがモノマー（二重結合あり）」であることの確認
- 初期状態で正しくグループ割り当てされていれば、これらは自動的に満たされる

**判断**: 現行の2グループ実装を維持する。理由:
1. 連鎖成長ロジック（propagation_map）が同じ機能を提供
2. i-k/j-l は結合存在の確認であり、グループ割り当てが正しければ冗長
3. 4グループ拡張は epoxy 系（T-G3）で必要だが、vinyl では不要

**Table S2 (Nylon-6,6)**:
| Reaction | i-j | i-k | j-l | V^f | V^d |
|---|---|---|---|---|---|
| Condensation | 3.0 / 6.0 | 0.0 / 3.0 | 0.0 / 3.0 | i-j, k-l | i-k, j-l |

**⚠️ 重要**: Nylon では formation (i-j = N-C amide, k-l = O-H water) **と** dissociation (i-k = N-H breaking, j-l = C-O breaking) が**同時に**適用される。現行コードでは `is_formation` は True/False の二値だが、nylon ではペアごとに混在する。

**Table S3 (Epoxy curing)**:
全3反応パターンとも 4グループ (i,j,k,l), formation: i-j/k-l, dissociation: i-k/j-l

### P-6: 系の規模

PDF（p.21-23）:
- Vinyl: **200 monomer + 10 AIBN** (solvent-free bulk)
- Nylon-6,6: **100 hexamethylenediamine + 100 adipic acid** (equimolar)
- Epoxy: **100 DGEBA + 50 DETA + CuO slab**

### P-7: AIBN分解は**反応パターンに含まれている**

Fig. S1: 反応パターンに「Activation (AIBN decomposition)」が含まれる。Table S1 にも明記。
→ 現行の「pre-formed radical」設計とは異なる。ただし、論文の Activation 行は V^d (dissociation) を i-k, j-l に適用し、AIBN の C-N 結合を切断する。この実装には:
- AIBN 全分子の SMILES (`CC(C)(C#N)N=NC(C)(C)C#N`) からの3D座標生成
- 4グループテンプレートによるActivation反応パターンの追加

**判断**: T-G1 の「pre-formed radical」近似は維持する（decisions.md に記録済み）。AIBN分解の忠実な再現は Phase 10 に延期。

---

## 1. タスク一覧と依存関係

```
P-0 (PDF確定事項の反映) ──────────────────────────────────────────────┐
   ├─ P-0a: γ decisions.md 更新                                     │
   ├─ P-0b: unbiased=2000 確定                                      │
   ├─ P-0c: friction = 1.0 ps^{-1} 変更 (Ask-first)                │
   └─ P-0d: 温度デフォルト 500→333 K 変更                           │
                                                                     │
T-G1a (vinyl+AIBN OrbMol-v2 E2E) ──┐                               │
                                    │                               │
T-G2 (nylon-6,6) ──────────────────┤── T-V (定性比較検証) ──────────┘
                                    │
T-G3 (epoxy/CuO) ──────────────────┘
                                    
T-OD (OrbMol-v2 デフォルト化) ─── 全スクリプトに反映
```

**推奨作業順序**:
1. P-0 (PDF確定事項の反映) — 0.5日
2. T-OD (OrbMol-v2 デフォルト化) — 0.5日
3. T-G1a (vinyl+AIBN OrbMol-v2 E2E) — 1日
4. T-G2 (nylon-6,6) — 2日
5. T-G3 (epoxy/CuO) — 3日（4グループテンプレートの実運用初回）
6. T-V (定性比較検証) — 1日

---

## 2. 各タスクの詳細

### P-0: PDF確定事項の反映

#### P-0a: γ 単位の decisions.md 更新
- **変更ファイル**: `specs/decisions.md`
- **内容**: 「PENDING CONFIRMATION」を以下に更新:
  > PDF (p.7) confirmed: γ=1.0, unit not stated. biased/unbiased each 2000 steps (500 fs).
  > Maintaining kcal/(mol·step). Saturation at step 250 is within biased phase.
  > Fig. S4 (p.25-26) confirms γ acts as global scaling factor — unit choice affects rate but not relative trends.
- **テスト**: なし（文書更新のみ）

#### P-0b: unbiased ステップ数確定
- **変更ファイル**: `paper/notes.md`, `configs/boost/paper_faithful.yaml`
- **内容**: `unbiased_steps: 2000` のコメントを "confirmed from PDF p.7" に更新
- **テスト**: なし

#### P-0c: Langevin friction の修正
- **現状**: `friction_per_fs=0.01` (= 10 ps^{-1})
- **論文**: 1.0 ps^{-1} = 0.001 fs^{-1}
- **変更ファイル**: 全 `run_*.py` スクリプトの `LangevinParams` デフォルト値
- **⚠️ Ask-first**: `LangevinParams` のデフォルト値変更は科学的影響あり。具体的には:
  - 現行 0.01/fs: 強い熱浴結合 → 温度安定だが動力学に干渉
  - 論文 0.001/fs: 弱い結合 → より物理的な動力学
- **対処**: `LangevinParams` の `friction_per_fs` デフォルトを `0.001` に変更。全 run スクリプトから明示的な `friction_per_fs` 指定を除去（デフォルト使用）。`decisions.md` に記録。
- **テスト**: 既存テストは `LangevinParams(temperature_K=300.0)` で friction 省略 → 新デフォルト適用。パスすることを確認。
- **受け入れ基準**: 全テストパス、`decisions.md` に変更記録

#### P-0d: 温度デフォルトの修正
- **現状**: `run_mace_pbc.py`, `run_vinyl_aibn.py` で `--temperature 500.0`
- **論文**: vinyl系は 333 K、nylon は 300 K、epoxy は 333 K
- **変更ファイル**: `scripts/run_vinyl_aibn.py` のデフォルト温度を `333.0` に変更、`scripts/run_mace_pbc.py` も `333.0` に変更
- **テスト**: なし（CLI デフォルト変更のみ）

---

### T-OD: OrbMol-v2 デフォルトバックエンド化

#### 目的
論文は PFP（proprietary）を使用。MACE-MP-0 は無機結晶向き。OrbMol-v2 は OPoly26（高分子 DFT データ）で訓練されており、有機高分子系により適合する。オーナー指示に従い OrbMol-v2 をデフォルトに変更する。

#### 変更内容

1. **`scripts/run_vinyl_aibn.py` の修正**:
   - import を `mace_backend` → `orb_backend` に変更
   - `create_mace_calculator` → `create_orb_calculator` に変更
   - `--model` CLI 引数を削除（OrbMol-v2 は model 指定不要）
   - `--backend` 引数を追加: `choices=['orb', 'mace']`, default='orb'
   - backend 引数に応じて `create_orb_calculator` or `create_mace_calculator` を呼び分け
   - PBC セルを渡す（OrbMol-v2 は PBC 対応、ただし `nvalchemiops` なしでは PME なし）
   - barostat はそのまま維持

2. **`scripts/run_mace_pbc.py` の修正**:
   - `--backend` 引数を追加（default='mace' — このスクリプトは MACE+PBC テスト用に維持）

3. **新規 `scripts/run_vinyl_aibn_orb.py` は作らない** — `run_vinyl_aibn.py` に `--backend` で統合。

4. **`specs/decisions.md` に記録**:
   > OrbMol-v2 selected as default backend for polymer systems per owner instruction.
   > Reason: OPoly26 training data includes polymer-relevant DFT (ωB97M-V/def2-TZVPD).
   > MACE-MP-0 retains compatibility via --backend mace.
   > nvalchemiops (periodic PME) status: blocked_pending_review.
   > Non-periodic or short-range-only periodic runs do not require nvalchemiops.

5. **`pyproject.toml` の修正**: `[project.optional-dependencies]` に以下を追加:
   ```
   default = [
       "orb-models>=0.7",
       "ase>=3.22",
   ]
   ```

#### 受け入れ基準
- `python scripts/run_vinyl_aibn.py --seed 7 --n-monomers 2 --n-initiators 1 --n-cycles 1 --biased-steps 10 --unbiased-steps 10 --box-size 18` が OrbMol-v2 で完走
- `--backend mace` で MACE-MP-0 にフォールバック可能
- `decisions.md` に変更記録
- 既存テスト全パス

---

### T-G1a: vinyl+AIBN 系 OrbMol-v2 での論文スケール E2E 実行

#### 目的
T-G1 で実装した vinyl+AIBN 系を OrbMol-v2 バックエンドで論文スケール（2000+2000 steps × 10 cycles）実行し、定性トレンドを確認する。

#### 前提
- P-0 (温度333K, friction 0.001/fs) 完了後
- T-OD (OrbMol-v2 デフォルト化) 完了後

#### 実行コマンド

```bash
# 小規模テスト (CPU 可, ~10分)
python scripts/run_vinyl_aibn.py \
    --seed 7 --output-dir runs/vinyl_aibn_orb_test \
    --n-monomers 4 --n-initiators 1 --n-cycles 2 \
    --biased-steps 100 --unbiased-steps 100 \
    --box-size 14.0 --temperature 333.0

# 論文スケール (GPU 推奨, CPU で ~数時間)
python scripts/run_vinyl_aibn.py \
    --seed 7 --output-dir runs/vinyl_aibn_orb_paper \
    --n-monomers 8 --n-initiators 2 --n-cycles 10 \
    --biased-steps 2000 --unbiased-steps 2000 \
    --box-size 16.0 --temperature 333.0

# 図生成
python scripts/reproduce_figures.py \
    --trajectory runs/vinyl_aibn_orb_paper/trajectory.jsonl \
    --bonds runs/vinyl_aibn_orb_paper/bonds.jsonl \
    --n-reactive-sites 18 \
    --target-temperature 333.0 \
    --output-dir runs/vinyl_aibn_orb_paper/figures
```

注: `--n-reactive-sites 18` = 8 monomers × 2 (alpha+beta) + 2 initiators = 18

#### ⚠️ 計算コスト見積もり
- 8 monomer × 12 atoms + 2 initiator × 12 atoms = **120 atoms** (H 含む)
- 2000+2000 × 10 cycles = **40,000 steps** × OrbMol-v2 compute
- CPU 見積もり: ~1-3 秒/step × 40,000 = ~11-33 時間
- GPU 見積もり: ~0.1 秒/step × 40,000 = ~1.1 時間
- **⚠️ Ask-first トリガー 7: 大規模計算の本実行前にオーナーへ確認**
- 注: 論文規模は 200 monomer + 10 AIBN だが、計算資源の制約から 8+2 で開始するのが妥当

#### 受け入れ基準
- E2E 実行が完走し `summary.json`, `trajectory.jsonl`, `bonds.jsonl` が出力される
- `confirmed_formations >= 1`（最低1つの結合形成）
- 温度が 333 K 付近で安定（trajectory の temperature_K を確認）
- `specs/figure-comparison.md` を更新: energy/conversion/temperature プロットの定性比較
- git SHA + seed + config が `manifest.json` に記録される

---

### T-G2: nylon-6,6 ステップ成長重合

#### 目的
論文の nylon-6,6 系を実装し、Carothers 方程式（DPn vs 転化率の非線形関係）を定性的に再現する。

#### 論文の仕様（PDF p.22）
- 系: hexamethylenediamine 100 mol + adipic acid 100 mol（等モル）
- 温度: 300 K, 圧力: 1 atm, NPT
- 反応: 縮合 (amine N + carboxylic acid C → amide, + H2O 生成)
- Table S2 より: i-j: 3.0-6.0 Å (formation: i-j + k-l), dissociation: i-k + j-l

#### 設計
1. **SMILES**:
   - hexamethylenediamine: `NCCCCCCN`
   - adipic acid: `OC(=O)CCCCC(=O)O`

2. **反応テンプレート（4グループ）**:
   - Gi: amine N (1° amine の N)
   - Gj: carboxylic acid C (C=O の C)
   - Gk: amine H (N-H の H) ← dissociation でN-H 結合切断
   - Gl: carboxylic acid OH (O-H の O) ← dissociation でC-OH 結合切断
   - Pair set P:
     - (i, j): formation (N-C amide bond) — r_min=3.0, r_max=6.0
     - (k, l): formation (O-H water formation) — r_min=0.0, r_max=3.0（近接確認）
     - (i, k): dissociation (N-H bond breaking) — r_min=0.0, r_max=3.0
     - (j, l): dissociation (C-OH bond breaking) — r_min=0.0, r_max=3.0

3. **`PairSpec.is_formation`**: nylon では**ペアごとに異なる** (i-j: True, k-l: True, i-k: False, j-l: False)。
   現行の `PairSpec` データクラスはこれをサポート済み（`is_formation: bool` がペアごと）。

4. **系ビルダー追加先**: `scripts/_systems.py` に `build_nylon66_system()` を追加

5. **原子特定ヘルパー**: RDKit で:
   - hexamethylenediamine: 末端 N (2つ)
   - adipic acid: カルボニル C (2つ)、ヒドロキシル O (2つ)、N-H の H、O-H の H
   - これらの正確な SMARTS パターンによる特定が必要

6. **水の生成**: 論文では水が生成される（Fig. S2）。現行コードには原子の追加/削除の機構がない。
   - **判断**: 水の生成はシミュレートしない。論文も「水除去なし」条件で実施（p.11: "without continuous water removal"）。TDBB は原子間距離のバイアスのみ — 原子種変更は行わない。反応は「N-C 距離の短縮」と「N-H/C-OH 距離の伸長」として表現され、TDBB が自然に駆動する。

7. **Carothers 比較**: `src/analysis/` に `carothers.py` を新規作成
   - DPn = 1/(1-p) (p = 転化率α) — 論文 Fig. 4c の理論曲線
   - `reproduce_figures.py` に DPn vs α プロット関数を追加

#### 変更ファイル
| ファイル | 変更内容 |
|---|---|
| `scripts/_systems.py` | `_find_amine_n()`, `_find_carboxyl_c_o()`, `build_nylon66_system()` 追加 |
| `scripts/run_nylon66.py` | 新規実行スクリプト |
| `src/analysis/carothers.py` | DPn(α) 理論曲線関数 |
| `scripts/reproduce_figures.py` | `plot_carothers_comparison()` 追加 |
| `tests/unit/test_systems.py` | nylon 系ビルダーのテスト追加 |
| `specs/decisions.md` | T-G2 設計判断記録 |

#### 受け入れ基準
- `python scripts/run_nylon66.py --seed 7 --n-cycles 1 --biased-steps 10 --unbiased-steps 10` が完走
- 系ビルダーが化学的に妥当なグループ割り当てを返す
- 4グループ+混合 formation/dissociation テンプレートのテスト
- Carothers 理論曲線との比較プロットが生成される
- `decisions.md` に設計判断記録

---

### T-G3: epoxy/CuO 界面硬化

#### 目的
論文の DGEBA+DETA on CuO 系を実装。4グループテンプレートの完全適用と界面反応の定性再現。

#### 論文の仕様（PDF p.22-24）
- 系: DGEBA 100 + DETA 50 + CuO(001) スラブ (8×8×6 supercell + 20 Å vacuum)
- 温度: 333 K, **NVT** (バロスタットなし)
- 反応パターン（3種, Table S3）:
  1. RO by primary amine: i-j/k-l formation, i-k/j-l dissociation
  2. RO by secondary amine: 同上
  3. RO by hydroxyl group: 同上
- 全パターン: i-j 3.0-6.0 Å, i-k/j-l 0.0-3.0 Å

#### 設計
1. **SMILES**:
   - DGEBA: `CC(C)(c1ccc(OCC2CO2)cc1)c1ccc(OCC2CO2)cc1` (簡略化, RDKit で構造生成)
   - DETA: `NCCNCCN`
   - CuO スラブ: ASE の `surface()` ビルダー + 表面水酸化

2. **4グループテンプレート**:
   - Gi: epoxy O
   - Gj: 1° amine N
   - Gk: 2° amine N
   - Gl: surface OH
   - P = {(i,j), (i,k), (j,l)}

3. **⚠️ 制約**: CuO スラブ生成には ASE の `ase.build.surface()` が必要。Cu, O の MACE-MP-0 / OrbMol-v2 対応を確認する必要がある。

4. **NVT 実施**: 論文は epoxy 系を NVT で実施。`--no-barostat` フラグで対応済み。

#### 変更ファイル
| ファイル | 変更内容 |
|---|---|
| `scripts/_systems.py` | `build_epoxy_cuo_system()` 追加 |
| `scripts/run_epoxy_cuo.py` | 新規実行スクリプト |
| `tests/unit/test_systems.py` | epoxy 系ビルダーのテスト追加 |
| `specs/decisions.md` | T-G3 設計判断記録 |

#### 受け入れ基準
- E2E 短縮実行が完走
- 4グループテンプレート + 3反応パターンのテスト
- 界面モデルが化学的に妥当（CuO スラブ + 有機層の配置）
- `decisions.md` に設計判断記録

#### ⚠️ Ask-first
- CuO スラブの結晶構造データの出典を確認
- OrbMol-v2 が Cu を適切に扱えるか検証（OPoly26 には金属表面は含まれていない可能性）

---

### T-V: 定性比較検証

#### 目的
実装した全化学系について、論文の定性トレンドとの比較を体系的に行い、`specs/figure-comparison.md` を更新する。

#### 比較項目

| 論文 Fig | 内容 | 対応プロット | 検証基準 |
|---|---|---|---|
| Fig. 2b | styrene α(t) (I/M比依存) | vinyl α(t) 単調増加 | α(t) が時間とともに増加 ✅ |
| Fig. 2c | Mn vs conversion (線形) | — | 分子量計算は未実装（延期可） |
| Fig. 2d | Rp vs [I][M] (線形) | — | 複数 I/M 比の実行が必要（延期可） |
| Fig. 3a | α(t) + exp fit | vinyl α(t) + Eq.11 fit | kp_eff が正の値で fit |
| Fig. 4b | nylon α(t) (60-80%で飽和) | nylon α(t) | 完全転化しない（平衡限界） |
| Fig. 4c | DPn vs α (Carothers) | DPn vs α plot | 理論曲線と定性一致 |
| Fig. 5a | epoxy c(t)/c(0) (種濃度変化) | — | 濃度追跡は実装後 |
| Fig. 5b | ρ_rxn(z) depth profile | density_profile.png | 界面付近で反応密度が高い |

#### 対処方法

1. 各系の論文スケール実行を行う（T-G1a, T-G2, T-G3 の受け入れ基準参照）
2. `scripts/reproduce_figures.py` で全プロットを生成
3. `specs/figure-comparison.md` を更新: 各図について「定性一致 / 部分一致 / 未達」を判定
4. 不一致があればその理由を分析（バックエンド差、系規模の差など）

#### 受け入れ基準
- 全化学系の論文スケール実行が完了
- Fig. 2b 相当 (α(t) 増加) が vinyl 系で確認 ← **最重要**
- Fig. 4b 相当 (α(t) 飽和) が nylon 系で確認
- Fig. 3a 相当 (Eq.11 fit) が vinyl 系で kp_eff > 0
- `specs/figure-comparison.md` が全図についてステータスを記載

---

## 3. 未着手の論文要素（延期項目）

以下は論文に含まれるが、現段階では実装しない項目。

| 項目 | 理由 |
|---|---|
| AIBN 熱分解（Activation反応パターン）| pre-formed radical 近似で代替（decisions.md 記録済み）|
| 連鎖停止・連鎖移動 | 論文も「excluded」と明記（p.21: "termination and chain-transfer reactions were excluded"）|
| 200 monomer スケール | 計算コスト制約。8 monomer で定性トレンド確認後に検討 |
| Mn (数平均分子量) 計算 | 結合トポロジー解析の追加実装が必要。Phase 10 |
| Fig. S4 感度分析 | γ スケーリングの確認は論文スケール実行後に検討 |

---

## 4. 各タスクが触るファイルの対応表

| タスク | 変更ファイル | 新規ファイル |
|---|---|---|
| P-0a | `specs/decisions.md` | — |
| P-0b | `paper/notes.md`, `configs/boost/paper_faithful.yaml` | — |
| P-0c | `src/integrators/langevin.py`, `specs/decisions.md` | — |
| P-0d | `scripts/run_vinyl_aibn.py`, `scripts/run_mace_pbc.py` | — |
| T-OD | `scripts/run_vinyl_aibn.py`, `pyproject.toml`, `specs/decisions.md` | — |
| T-G1a | — | `runs/vinyl_aibn_orb_paper/` (出力) |
| T-G2 | `scripts/_systems.py`, `scripts/reproduce_figures.py`, `tests/unit/test_systems.py`, `specs/decisions.md` | `scripts/run_nylon66.py`, `src/analysis/carothers.py` |
| T-G3 | `scripts/_systems.py`, `tests/unit/test_systems.py`, `specs/decisions.md` | `scripts/run_epoxy_cuo.py` |
| T-V | `specs/figure-comparison.md` | — |

---

## 5. 質問必須トリガー（再掲）

1. P-0c: friction 値の変更 (0.01 → 0.001) — 科学的影響あり
2. T-G1a: 論文スケール実行前のコスト確認
3. T-G3: CuO スラブの元素対応 (OrbMol-v2 + Cu)
4. 論文が複数の妥当な数学的解釈を許す場合
5. 商用利用権が不明なバックエンドを追加しようとする場合
6. 図再現でスムージング・フィルタ・平均化を論文と変える場合

---

## 6. 完了の定義

本計画の完了条件:
1. P-0 の全 PDF 確定事項が反映され、`decisions.md` に記録されている
2. OrbMol-v2 がデフォルトバックエンドとして全 run スクリプトで使用可能
3. vinyl+AIBN, nylon-6,6, epoxy/CuO の3系がE2E実行可能
4. 各系の定性比較が `specs/figure-comparison.md` に記載されている
5. 全ユニットテストがパスしている
6. 全実行の seed, config path, git SHA, backend, output_dir が `manifest.json` に記録されている
