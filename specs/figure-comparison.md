# Figure Comparison: Paper vs Reproduction

作成日: 2026-06-13
論文: arXiv:2511.22874, Mori et al.
再現 run:
- `runs/toy_bond_demo/`  (T8.1: toy LJ system, confirmed_formations=1)
- `runs/mace_pbc_paper/` (T8.2: MACE-MP-0 + PBC, 4 ethylene molecules, 3 cycles)
- `runs/vinyl_aibn_gpu40/` (T-G1a paper-scale: OrbMol-v2, 40 monomer + 2 AIBN = 504 atoms, density 0.35 g/mL, NPT 333 K, 3 cycles × 2000+2000 = 12,000 steps, seed 42, RTX 4060 Ti)
- `runs/vinyl_aibn_paper100/` (2026-06-15: OrbMol-v2, 100 monomer + 5 initiator = 1260 atoms, **classical OpenFF/Sage prep to 0.50 g/mL in WSL**, NVT 333 K, FIRE minimize + 2000 ML equil + 3×2000+2000 = 14,000 steps, seed 42, RTX 4060 Ti)

---

## 2026-06-15 追記: paper100(古典prep + NVT)結果

WSL の OpenFF/Sage で 0.50 g/mL に高密度化した構造を `--load-structure` で読み込み、Windows GPU で NVT 本番を完走。

| 項目 | 結果 | 判定 |
|---|---|---|
| E2E 完走（14,000 steps, 1260 atoms, 0.50 g/mL） | ✅ RC=0, summary/trajectory/bonds/figures 出力 | 達成 |
| 温度安定（333 K 付近） | ✅ 全体 mean 295 K / max 336 K（旧 unprepped: mean 527 / max 1364）, biased 300 / unbiased 314 K | 達成（スパイク解消） |
| VRAM（16 GB GPU） | ✅ 2.4–4.7 GB / 16 GB, util 75–87% | 達成 |
| 候補検出・選択・bias 印加 | ✅ biased 各 14/11/19 候補、5/4/5 選択、bias_E=1000–1250（f1_max=250/pair フル） | 達成 |
| confirmed_formations ≥ 1 | ❌ formations=0, dissociations=0, propagation=0 | **未達** |

**結論**: インフラ的な障壁(論文密度での GPU 完走)は**解決**。`formations=0` は密度・規模・温度のいずれでもなく、**TDBB バイアスの捕捉範囲(~2.5 Å)と候補リスト範囲(3–6 Å, Table S1)の不一致**に切り分け済み(decisions.md 2026-06-15 参照)。これは TDBB の科学的意味に関わる Ask-first 事項で、論文の力場定義の再読+承認が必要。

**注**: 200+10(2520 atoms)は OrbMol-v2 の単発フットプリント ~9.5 GB のため 16 GB では sustained 実行不可(24 GB+ GPU が必要)。

## 2026-06-14 追記: gpu40 論文スケール検証結果

| 項目 | 結果 | 判定 |
|---|---|---|
| E2E 完走（12,000 steps） | ✅ exit 0, summary/trajectory/bonds/figures 出力 | 達成 |
| 候補検出・選択・bias 印加 | ✅ 各 biased で 2/4/5 候補、2 選択、bias_E=500 kcal/mol | 達成（machinery 正常） |
| confirmed_formations >= 1 | ❌ formations=0, dissociations=0, propagation=0 | **未達** |
| 温度安定（333 K 付近） | ❌ mean 527 K, std 182 K, max 1364 K（step 51 で 1364 K に急騰） | **未達** |
| 温度の時間収束 | △ 2000step毎: 850→607→516→438→392→360 K（最終 cycle で漸近） | 部分 |

**根本原因（decisions.md 2026-06-14 参照）**: スクリプトに **エネルギー最小化・平衡化が無く**、格子配置の初期構造の近接衝突 → 初期温度爆発 → 候補ペアが結合距離(~2-2.5 Å)に収束せず → formations=0。論文（PDF p.20）は production 前に NPT 平衡化を実施。本ワークフローは未実装。

### 2026-06-14 追記2: 平衡化実装後の再検証（minimize + equil）

pre-TDBB 段階（FIRE 最小化 + 2000 step NPT 平衡化）を実装し gpu40 を再実行（`runs/vinyl_aibn_gpu40/`、baseline は `runs/vinyl_aibn_gpu40_no_equil/`）。

| 指標 | no-equil | minimize+equil | 判定 |
|---|---|---|---|
| 全体温度 mean/std/max | 527 / 182 / 1364 K | **310 / 40 / 373 K** | **解消** |
| biased phase 平均温度 | 586 K | **316 K** | 設定333Kに一致 |
| unbiased phase 平均温度 | 469 K | **330 K** | 設定333Kに一致 |
| 初期温度スパイク | 1364 K | 373 K（最大） | **解消** |
| confirmed_formations | 0 | **0（変化なし）** | 未達 |

**結論**: 平衡化実装は**温度不安定バグを解消**（real fix, 保持）。しかし formations=0 は不変で、これは温度問題ではなく **bias 捕捉範囲(~2.5 Å) と候補リスト範囲(3-6 Å, Table S1) の不一致＋系規模**が原因（V^f の井戸幅 ~0.32 Å、4-6 Å では力ほぼ0）。TDBB パラメータは論文確定値のため変更しない。

**次アクション**: formations>0 には論文規模（200 monomer + 10 AIBN）かつ複数サイクルが必要（大規模 GPU、オーナー承認待ち、Ask-first トリガー 2/7）。bias 井戸の拡幅は TDBB の科学的意味を変えるため論文再読＋承認なしには行わない（トリガー 3）。

## 図一覧と対応

| 再現図 | 論文 Fig | 比較結果 | 備考 |
|---|---|---|---|
| energy_vs_step.png | Fig. 2 (推定) | **定性一致** | biased/unbiased エネルギー分離を確認（OrbMol-v2 S2 含む） |
| base_energy.png | Fig. 2 (推定) | 定性一致 | base ポテンシャルの緩和トレンドを確認 |
| temperature_vs_step.png | -- | 参考 | NVT 333 K 安定。反応イベント時にスパイク（物理的に正） |
| conversion_vs_step.png | Fig. 3/4 (推定) | **定性一致** | α(t) 単調増加 + Eq.11 指数フィット overlay（OrbMol-v2 S2） |
| s2_diagnostics.png | -- (S2固有) | 検証用 | min_pair_dist + candidates per cycle、2 seed 比較 |
| conversion_vs_step.png (toy) | Fig. 3/4 (推定) | 定性一致 | α(t) 単調増加を確認（toy系、以前のもの） |

---

## 2026-06-18: S2 figure 生成 (OrbMol-v2, formations > 0)

対象 run: `runs/s2_sweep5_seed7/` (20+1, f2=5, [3,6] window, 15 cycles, seed 7, confirmed_formations=2)
比較用: `runs/s2_sweep4/` (同条件, 5 cycles, seed 42, confirmed_formations=1)

### conversion_vs_step.png (S2, OrbMol-v2)

- **初の MLIP 系での非ゼロ α(t)。** 2 段の階段状転化率（cycle 6, 12 で各 1 formation）
- Eq. 11 指数フィット overlay: kp_eff = 1.58e-06 fs^-1, R^2 = 0.653
- R^2 が低いのは 2 events のみの階段関数への指数フィットのため。統計的にはデータ不足だが、手法の動作実証として有効
- α_max = 9.5% (2/21 reactive sites) -- 論文の 60-80% には規模・サイクル数が不足
- Artifact: `runs/s2_sweep5_seed7/figures/conversion_vs_step.{png,pdf}`

### energy_vs_step.png (S2, OrbMol-v2)

- biased(赤)/unbiased(青) の明確な分離: 論文 Fig. 2 と定性一致
- biased phase で total energy 400-500 kcal/mol (bias 250 + base ~200)
- unbiased phase で 150-200 kcal/mol に緩和
- 15 cycle にわたりパターンが安定（TDBB machinery の安定動作を確認）
- Artifact: `runs/s2_sweep5_seed7/figures/energy_vs_step.{png,pdf}`

### temperature_vs_step.png (S2, OrbMol-v2)

- NVT 333 K に安定。平常時 300-350 K 範囲
- step ~17,600 と ~32,100 に温度スパイク（~490 K, ~600 K）-- 反応イベント（C-C 結合形成）に伴うエネルギー放出。Langevin が数百 step で収束。物理的に正しい挙動
- Artifact: `runs/s2_sweep5_seed7/figures/temperature_vs_step.{png,pdf}`

### s2_diagnostics.png (S2 固有検証図)

- 上段: min_pair_distance vs cycle（2 seed 比較）
  - 反応 cycle (★) で min_pair_dist が ~2 Å に急降下（r0=2.04 Å 以下）
  - 非反応 cycle は 3.2-5.0 Å 範囲に分布（f2=5 capture radius ~3.2 Å の上下）
  - f2=10 capture radius (2.87 Å) と f2=5 capture radius (3.2 Å) の参照線で差異を可視化
- 下段: candidates per cycle（棒グラフ）-- [3,6] window で常に 2-5 candidates
- Artifact: `runs/s2_figures/s2_diagnostics.{png,pdf}`

### 非適用図（vinyl 周期系）

| 図 | 理由 |
|---|---|
| density_profile (Eq. 12) | 界面/硬化系の深さ方向密度用。等方的な周期 melt box では z 方向に構造がない |
| Carothers DPn (Fig. 4c) | ステップ成長重合用。ビニルラジカル重合は連鎖成長であり Carothers は科学的に不適合 |

---

## 詳細比較（以前の run、参考）

### 1. energy_vs_step（エネルギー vs ステップ）

**論文が示すもの（推定）**:
- biased フェーズでエネルギーが高くなる（バイアス印加）
- unbiased フェーズで緩和してエネルギーが下がる
- サイクルを繰り返すごとにベースエネルギーが低下（重合によるエネルギー利得）

**再現図（MACE+PBC, runs/mace_pbc_paper/figures/energy_vs_step.png）**:
- ✅ biased（赤）と unbiased（青）の明確な分離を確認
- ✅ bias エネルギーが biased フェーズで上昇し unbiased で 0 に戻る
- △ ベースエネルギーの単調低下は不明瞭（3 サイクルかつ formations=0 のため）

**差異の説明**:
- 論文は PFP（Matlantis）使用、本再現は MACE-MP-0。有機系での精度差あり。
- 論文は confirmed_formations > 0 の系（エネルギー利得あり）。本再現は formations=0。
- バックエンドの違いによる定量差は `specs/decisions.md` 2026-06-12「Using MACE instead of paper's PFP backend」に記録済み。

### 2. base_energy（base ポテンシャルエネルギー）

**論文が示すもの（推定）**:
- 全体的に下降トレンド（重合反応による安定化）

**再現図**:
- △ MACE+PBC 系では下降トレンドは確認できていない（formations=0、短いサイクル数）
- ✅ toy_bond_demo ではサイクル 0 後に base エネルギーが低下（LJ well への落下）を確認

### 3. temperature_vs_step（瞬間温度）

**論文対応**: 論文の Fig. には温度プロットの直接対応不明（NVT 検証用）

**再現図（MACE+PBC）**:
- ✅ Langevin サーモスタットで目標温度 500 K 周辺に推移
- △ 序盤に大きな温度スパイクがある（初期速度ゼロからの加速）

**差異の説明**:
- 初期速度ゼロは意図的（Maxwell-Boltzmann 初期化が未実装）。
- 実験比較には初期速度を適切に設定する改善が推奨される（follow-up task）。

### 4. conversion_vs_step（転化率 α(t)）

**論文が示すもの（推定, Eq. 11-12）**:
- サイクル数とともに α(t) が単調増加
- 最終的に 60-80%（ポリマー系）または 95%+（架橋系）に到達

**再現図（toy_bond_demo）**:
- ✅ サイクル 0 終了時（step=600）に α が 0 → 50%（2 反応サイト中 1 形成）に上昇
- ✅ 単調増加のトレンドを確認（TDBB 機械の動作証明）
- 分母: n_reactive_sites=2（trajectory header から正しく読み取り）

**再現図（MACE+PBC, runs/mace_pbc_paper）**:
- △ α=0 のまま（formations=0）。単調増加なし。
- 原因: エチレン直接 C-C 結合障壁 ~40+ kcal/mol、PBC あり 500+500 ステップでも不十分
- → 論文スケール（2000+2000 ステップ）かつ適切な反応系（radical initiator あり）で改善見込み

---

## 未対応の論文図

| 論文 Fig | 内容（推定） | 再現状態 | 理由 |
|---|---|---|---|
| Fig. 5 | depth-resolved reaction density ρ_rxn(z) | 未対応 | formations=0 のため密度計算不可。T8.2 で formations が得られたら追加。 |
| Fig. 6 | per-system (nylon/epoxy) 比較 | 未対応 | nylon/epoxy 系は T10.x（オーナー承認後） |
| 感度分析 Fig. S4 | γ のスケーリング効果 | 未対応 | γ 単位確定後（T6.1 PDF 確認待ち） |

---

## 定性一致の評価

| 項目 | 状態 |
|---|---|
| biased/unbiased エネルギー分離 | ✅ 確認済み（toy + OrbMol-v2 S2） |
| Langevin 温度制御 | ✅ 確認済み（333 K 安定、反応時スパイクは物理的に正） |
| TDBB 機械による bond formation | ✅ toy + OrbMol-v2 両方で確認 |
| MLIP 系での bond formation | ✅ **達成**（S2: OrbMol-v2, f2=5, 3 formations / 20 cycles） |
| α(t) 単調増加トレンド | ✅ **達成**（S2 sweep5: 0→4.8%→9.5% 階段状増加） |
| Eq. 11 指数フィット | ✅ 実装済み（kp_eff, R^2 表示。統計的にはデータ点不足） |
| 深さ分解密度プロット | N/A vinyl 周期系では非適用（界面/硬化系用） |
| Carothers DPn | N/A vinyl 連鎖成長では非適用（ステップ成長用） |

---

## 再現コマンド

```bash
# S2 sweep5_seed7 (OrbMol-v2, formations=2, 主要 figure)
python scripts/run_vinyl_aibn.py \
    --n-monomers 20 --n-initiators 1 --initiator-smiles "C[C](C)C#N" --spin 2 \
    --density 0.5 --backend orb --device cuda --no-barostat \
    --n-cycles 15 --biased-steps 2000 --unbiased-steps 500 --equil-steps 2000 \
    --timestep-fs 1.0 --f2 5.0 --seed 7 --output-dir runs/s2_sweep5_seed7
python scripts/reproduce_figures.py \
    --trajectory runs/s2_sweep5_seed7/trajectory.jsonl \
    --bonds runs/s2_sweep5_seed7/bonds.jsonl \
    --n-reactive-sites 21 --target-temperature 333.0 --timestep-fs 1.0 \
    --output-dir runs/s2_sweep5_seed7/figures

# S2 diagnostics (sweep4 + sweep5 comparison)
python scripts/reproduce_figures.py \
    --summary runs/s2_sweep4/summary.json runs/s2_sweep5_seed7/summary.json \
    --output-dir runs/s2_figures

# T8.1 toy system (確定: formations=1)
python scripts/run_toy_bond_demo.py --seed 7 --output-dir runs/toy_bond_demo
python scripts/reproduce_figures.py \
    --trajectory runs/toy_bond_demo/trajectory.jsonl \
    --bonds runs/toy_bond_demo/bonds.jsonl \
    --output-dir runs/toy_bond_demo/figures

# T8.2 MACE + PBC (中規模: 500+500 x 3 cycles)
python scripts/run_mace_pbc.py --seed 7 --output-dir runs/mace_pbc_paper \
    --n-molecules 4 --n-cycles 3 --biased-steps 500 --unbiased-steps 500 \
    --box-size 14.0 --temperature 500.0
python scripts/reproduce_figures.py \
    --trajectory runs/mace_pbc_paper/trajectory.jsonl \
    --bonds runs/mace_pbc_paper/bonds.jsonl \
    --n-reactive-sites 8 --target-temperature 500.0 \
    --output-dir runs/mace_pbc_paper/figures

# T8.2 論文スケール (2000+2000 steps — GPU 推奨, CPU で ~数時間)
python scripts/run_mace_pbc.py --seed 7 --output-dir runs/mace_pbc_paper_full \
    --n-molecules 4 --n-cycles 10 --biased-steps 2000 --unbiased-steps 2000 \
    --box-size 14.0 --temperature 500.0
```
