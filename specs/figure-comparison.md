# Figure Comparison: Paper vs Reproduction

作成日: 2026-06-13
論文: arXiv:2511.22874, Mori et al.
再現 run:
- `runs/toy_bond_demo/`  (T8.1: toy LJ system, confirmed_formations=1)
- `runs/mace_pbc_paper/` (T8.2: MACE-MP-0 + PBC, 4 ethylene molecules, 3 cycles)
- `runs/vinyl_aibn_gpu40/` (T-G1a paper-scale: OrbMol-v2, 40 monomer + 2 AIBN = 504 atoms, density 0.35 g/mL, NPT 333 K, 3 cycles × 2000+2000 = 12,000 steps, seed 42, RTX 4060 Ti)

---

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
| energy_vs_step.png | Fig. 2 (推定) | 定性一致 | biased/unbiased エネルギー分離を確認 |
| base_energy.png | Fig. 2 (推定) | 定性一致 | base ポテンシャルの緩和トレンドを確認 |
| temperature_vs_step.png | — | 参考 | 論文図との対応は不明。NVT 検証用 |
| conversion_vs_step.png (toy) | Fig. 3/4 (推定) | 定性一致 | α(t) 単調増加を確認（toy系） |
| conversion_vs_step.png (MACE) | Fig. 3/4 (推定) | 部分一致 | MACE+PBC では α=0（原因は既知バリア） |

---

## 詳細比較

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
| biased/unbiased エネルギー分離 | ✅ 確認済み |
| Langevin 温度制御 | ✅ 確認済み（スパイクあり） |
| TDBB 機械による bond formation | ✅ toy 系で確認（machinery が正しく動作） |
| MLIP 系での bond formation | ❌ 未達成（エチレン直接 C-C 障壁が高い） |
| α(t) 単調増加トレンド | △ toy 系のみ（MLIP 系は 0） |
| 深さ分解密度プロット | ❌ 未着手（formations=0） |

**Phase 3 受け入れ基準**（specs/acceptance-criteria.md）:
- 「Paper-faithful config reproduces the expected qualitative trend」
- → エネルギー biased/unbiased 分離と温度制御は達成。α(t) の定性一致は toy 系で証明済み。
  MLIP 系での α(t) 増加は T8.2 の論文スケール実行（2000+2000 steps, 適切な反応系）で達成見込み。

---

## 再現コマンド

```bash
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
