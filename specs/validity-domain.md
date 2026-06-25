# TDBB + MLIP 妥当域 (Validity Domain)

本ドキュメントは、商用安全な MLIP バックエンド上の TDBB ワークフローが
**どこまで物理的に信頼できるか**を、実測エビデンス付きで画定する。

プロジェクト目標(`specs/decisions.md` 2026-06-25「プロジェクト方針の確定」):
- (i) 使える商用安全 TDBB シミュレーション能力(道具)
- (ii) その妥当域の知見 ← **本ドキュメントが (ii) の中核成果物**

> **重要な訂正(2026-06-25 CORRECTION, decisions.md 参照)**: 当初本書は
> 「密な多ラジカル melt は妥当域外(OrbMol-v2 のスピンの壁)」としていたが、
> これは **誤り**だった。一連の爆発は **timestep 1fs + 平衡化不足の run 条件
> アーティファクト**であり、**論文値 0.25fs + 適切な平衡化**を使えば OrbMol-v2
> 単独で 200+10/spin=21 の多ラジカル melt が ~404K で安定に走り候補も潤沢
> (40/41/41)。よって多ラジカル melt は §1 の妥当領域に移動した。

凡例: ✅ 検証済み(妥当) / ⛔ 妥当域外(限界) / 🔶 部分的/要検証

---

## バックエンド(商用安全, 2 系統)

| バックエンド | License(コード/重み) | 特性 | 役割 |
|---|---|---|---|
| **OrbMol-v2** | Apache-2.0 / Apache-2.0 | 既統合・VRAM 実測済み(2640原子 9.5GB)。周期は nvalchemiops 要 | **主力** |
| **AIMNet2-NSE** | MIT / MIT | ラジカル/開殻特化訓練、spin-charge equilibration、周期(Ewald/PME)ネイティブ | **第2(クロス検証・ラジカル系)** |

両者で同一物理が再現されれば信頼度↑(backend-agnostic 原則, CLAUDE.md)。AIMNet2-NSE は
spike で PES 一致を確認(下記 §1.1)。OrbMol-v2 の代替でなく補完。

---

## 1. ✅ 妥当(信頼できる)領域

### 1.1 反応エネルギー論 — 単一ラジカル付加 (doublet) — 2 バックエンド一致
- メチルラジカル `[CH3]` + エチレン `C=C`、spin=2、forming C–C を 3.5→1.54 Å 拘束緩和スキャン(`scripts/scan_radical_addition.py`)。
- **OrbMol-v2**: 障壁 ~+6.1、井戸 −27.8 kcal/mol(1.54 Å)。**AIMNet2-NSE**: 障壁 ~+6.2、井戸 −28.7。教科書(障壁 ~7、発熱 ~−23)と一致。**2 バックエンドが一致**=クロス検証成立。

### 1.2 多ラジカル melt の同時シミュレーション(良条件下) ← 訂正で昇格
- **OrbMol-v2, 200+10, spin=21(高スピン和, 18解離), timestep 0.25fs + equil 2000**: 温度 **mean 404K / max 437K で安定**、qualified candidates **40/41/41 と全サイクル維持**(runs/orb_fullcond)。
- **必須条件**: timestep **0.25fs**(論文値)+ 平衡化。**1fs では数値発散して爆発**(§2.1)。
- AIMNet2-NSE も多ラジカル(40+4, spin=9, 0.25fs+equil)で安定(~875K、やや高め)+候補維持(17/15/17)。

### 1.3 閉殻系 / 段階重合 (spin=1 が物理的に正しい)
- 縮合重合(ナイロン6,6 等)は閉殻 → spin=1 が物理的に正しい。スピン問題なし。
- 🔶 ナイロンの障壁/キネティクス定量検証は未実施(オープン項目)。

### 1.4 インフラ — paper-scale が 16 GB GPU + WSL で実走可能
- 200+10 (2640 atoms, 0.5 g/mL) 持続 MD ピーク VRAM = 9.55 GB / 16 GB(WSL2 で expandable_segments 有効)。

---

## 2. ⛔ 妥当域外 / 要注意

### 2.1 反応性 MD を timestep 1fs で回すこと(数値不安定)
- **症状**: 多ラジカル melt を **1fs** で回すと温度が 1e6〜1e10 K に発散(s6_equil0/s6_fix_verify/s6_spincap1)。
- **原因**: 反応性・開殻・強バイアスの硬い系で 1fs は過大 → 数値ドリフト/発散。spin の問題ではなかった(2026-06-25 CORRECTION)。
- **対策**: 反応条件は **timestep 0.25fs(論文値)を既定**にする(§3)。

### 2.2 多ラジカル系の相対反応性(障壁差依存)
- モノマー間の相対反応性は正確なエネルギー論依存。OrbMol-v2 の絶対精度は未検証 → 定量比較は要検証。AIMNet2-NSE(ラジカル訓練)の方が有利な可能性、クロス検証で確認。

### 2.3 spin-agnostic 実行(OrbMol-v2/AIMNet2 とも不可)
- 両モデルとも charge+spin(total multiplicity)を要求。スピン非設定運用は不可。多ラジカルは高スピン和(N+1)を渡す(良条件下で安定、§1.2)。

### 2.4 エポキシ/CuO 界面(範囲外 — 既決)
- `decisions.md` 2026-06-20「再現スコープを有機系に限定」。

---

## 3. 運用ガイダンス(妥当域内で使うための既定)

| 項目 | 推奨 | 根拠 |
|---|---|---|
| **timestep(反応条件)** | **0.25 fs(論文値)を既定** | 1fs は反応系で数値発散; §1.2/§2.1, decisions 2026-06-25 CORRECTION |
| 平衡化 | 活性化後に十分な平衡化(例 equil 2000) | §1.2 |
| 活性化順序 | activate を 333K 平衡化の前に | decisions 2026-06-24(commit 165295e) |
| charge/spin | 必ず両方設定。多ラジカルは高スピン和 N+1 | §2.3 |
| 単一ラジカル/閉殻 | doublet=2 / singlet=1 | §1.1/§1.3 |
| 周期境界 | OrbMol-v2: nvalchemiops。AIMNet2: ネイティブ | backend |

---

## 4. ガードレール(ワークストリーム B で実装予定)
道具(i)として無効領域で黙って壊れないため:
- 反応条件で timestep が過大(>0.5fs 等)→ 警告
- 生産中の温度が設定値を大きく超過(例 >3×setpoint)→ 早期に異常検知して停止
- spin/charge 未設定での呼び出し → 明示エラー

---

## 5. 再現可能な validation スイート

| スクリプト/run | 検証内容 | 状態 |
|---|---|---|
| `scripts/scan_radical_addition.py --backend orb/aimnet` | 単一ラジカル付加の障壁 vs 教科書(2 backend 一致) | ✅ |
| `scripts/scan_radical_addition_2rad.py` | 高スピン近似の誤差(2 ラジカル) | ✅ |
| `scripts/diag_spin_sweep.py` | force vs spin multiplicity | ✅ |
| runs/orb_fullcond (200+10, 0.25fs, equil2000) | 多ラジカル melt の温度安定性 | ✅ |
| (新規予定) 縮合反応の障壁チェック | ナイロン amide 化の障壁 vs 参照 | 🔶 |

---

## 6. オープン項目(次タスク)
- **OrbMol-v2 + 0.25fs + equil で 200+10 を長尺実行し転化率を実測**(元の S6 目標が達成可能か; biased を 2000×0.25fs に)。
- timestep 0.25fs を反応条件の既定に(S6 スクリプト/configs を更新; 現状 1fs)。
- ナイロン(段階重合)の障壁/キネティクス検証 → §1.3 を ✅ に昇格。
- AIMNet2-NSE と OrbMol-v2 のクロス検証(同系で formation/温度/転化率を比較)。
- ガードレール実装(§4)。
- 目標スコープ(論文追試の要否)を新事実で再検討(decisions 2026-06-25 CORRECTION の follow-up c)。
