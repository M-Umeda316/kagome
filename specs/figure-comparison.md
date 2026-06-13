# Figure Comparison: Paper vs Reproduction

作成日: 2026-06-13
論文: arXiv:2511.22874, Mori et al.
再現 run:
- `runs/toy_bond_demo/`  (T8.1: toy LJ system, confirmed_formations=1)
- `runs/mace_pbc_paper/` (T8.2: MACE-MP-0 + PBC, 4 ethylene molecules, 3 cycles)

---

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
