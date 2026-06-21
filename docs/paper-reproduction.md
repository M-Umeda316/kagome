# 論文再現ガイド

このドキュメントでは、論文 (arXiv:2511.22874, Mori et al.) の手法を KAGOME で再現する方法と現在の進捗を記録します。

---

## 対象論文

**"Ready-to-Use Polymerization Simulations Combining Universal Machine Learning Interatomic Potential with Time-Dependent Bond Boosting for Polymer and Interface Design"**

- arXiv: 2511.22874
- 著者: Mori et al. (2025)
- 手法: TDBB (Time-Dependent Bond Boosting) + uMLIP (Universal Machine-Learned Interatomic Potential)

---

## 再現マイルストーン

| マイルストーン | 状態 | 主な成果 |
|---|---|---|
| S1: 単鎖伝播デモ | 完了 | pentamer 構築、ラジカル移動確認 |
| S2: メルト駆動での結合形成 | 完了 | OrbMol-v2, alpha=9.5%, kp_eff=1.58e-06 fs^-1 |
| S3: 多ラジカル・高スピン近似 | 完了 | 2-radical, alpha=27.3%, R^2=0.833 |
| S4: AIBN 活性化 + 連鎖重合 | 完了 | C-N 解離 + 結合形成実証 |
| S5: 図の再現・比較 | 一部完了 | energy, conversion, temperature 図生成済み |
| S6: 論文スケール実行 (200+10) | 未着手 | 24 GB+ GPU が必要 |

### S1: 単鎖伝播デモ

ビニルラジカル重合の連鎖伝播メカニズムをユニットレベルで検証。pentamer (5 量体) の構築でラジカルの移動 (atom indices 1→12→24→36→48) を確認。

### S2: メルト駆動での結合形成

OrbMol-v2 バックエンドで実際のメルト環境での TDBB による結合形成を初めて達成。

- 系: 20 モノマー + 1 開始剤, f2=5, [3,6] A 窓
- seed 7: 15 サイクルで 2 結合形成
- seed 42: 5 サイクルで 1 結合形成
- 転化率: alpha_max = 9.5%

### S3: 多ラジカル・高スピン近似

複数ラジカル系でのスピン取り扱いを検証。高スピン近似 (spin = 2 * n_radicals) の妥当性を確認。

- 系: 20 モノマー + 2 開始剤, 15 サイクル
- 6 結合形成 (steps 5620, 6729, 12422, 15681, 16423, 17594)
- 転化率: alpha_max = 27.3%
- Eq. 11 フィット: kp_eff = 1.07e-05 fs^-1, R^2 = 0.833

### S4: AIBN 活性化 + 連鎖重合

AIBN の C-N アゾ結合を V^d バイアスで解離し、生成したラジカルで連鎖重合を実行する完全ワークフロー。

- 活性化パラメータ: f2=0.3, f1_max=250 (OrbMol-v2 の C-N 障壁 39.4 kcal/mol に調整)
- 系: 5 モノマー + 1 AIBN
- C-N 解離: step 175-176 で 2/2 結合解離
- 結合形成: step 2269 で 1 結合
- 転化率: alpha = 7.1%

### S5: 図の再現

再現された図:

| 図 | 対応する論文図 | 内容 | 再現品質 |
|---|---|---|---|
| `energy_vs_step` | Fig. 2 (推定) | biased/unbiased エネルギー分離 | 定性的一致 |
| `conversion_vs_step` | Fig. 3-4 (推定) | alpha(t) + Eq. 11 フィット | 定性的一致 (スケール差あり) |
| `temperature_vs_step` | 検証用 | NVT 333 K 安定性 | 良好 |
| `s2_diagnostics` | 検証用 | min pair distance, candidates/cycle | 参考図 |

未再現 (理由付き):

| 図 | 理由 |
|---|---|
| Fig. 5 (深さ密度 rho_rxn(z)) | 界面/硬化系が必要 (ビニルメルトでは N/A) |
| Fig. 6 (系別比較) | ナイロン/エポキシ系が未実装 |
| Fig. S4 (gamma 感度) | 延期中 |

### S6: 論文スケール (未着手)

200 モノマー + 10 AIBN (2520 原子) での完全再現。24 GB+ GPU が必要。セットアップ手順は `specs/s6-environment-setup.md` に記載済み。

---

## 論文との差異

### バックエンド

論文は PFP (Matlantis) を使用していますが、商用ライセンスの制約により KAGOME は OrbMol-v2 (Apache-2.0) を使用しています。

| 項目 | 論文 (PFP) | KAGOME (OrbMol-v2) |
|---|---|---|
| C-N 障壁 | 不明 | 39.4 kcal/mol |
| 活性化 f2 | 10.0 | 0.3 (調整済み) |
| 活性化 f1_max | 不明 | 250 kcal/mol |
| スピン対応 | 不明 | 高スピン近似 |

### スケール

現在の最大再現は 20+2 (S3) で alpha=27.3%。論文は 200+10 で alpha=60-80% を報告。スケール差が転化率に直結するため、S6 でのスケールアップが必要です。

### アンサンブル

論文は NPT (1 atm) を使用。KAGOME は現在 NVT (`--no-barostat`) で実行。MC barostat は実装済みですが、バイアス力との相互作用を避けるため無効化しています。

---

## 方程式→コード対応表

| 方程式 | 内容 | 実装ファイル | 関数/クラス |
|---|---|---|---|
| Eq. 2 | V^f 結合形成バイアス | `boost/tdbb.py` | `formation_potential` |
| Eq. 3 | V^d 結合解離バイアス | `boost/tdbb.py` | `dissociation_potential` |
| Eq. 4 | r0 = lambda * sum(r_vdw) | `boost/tdbb.py` | `target_distance` |
| Eq. 5 | f1(t) = min(gamma*t, f1_max) | `boost/tdbb.py` | `boost_amplitude` |
| Eq. 6 | G_X 反応原子グループ | `reactive/groups.py` | `ReactiveGroup` |
| Eq. 7 | d_ijkl 候補スコアリング | `reactive/selection.py` | `find_candidates`, `score_candidates` |
| Eq. 8 | 全バイアスポテンシャル | `boost/tdbb.py` | `total_bias` |
| Eq. 11 | alpha(t) = 1 - exp(-kp*t) | `analysis/conversion.py` | `fit_conversion_exponential` |
| PDF p.12 | rho_rxn(z) 深さ密度 | `analysis/density.py` | `reaction_density_profile` |

Eq. 9 (Rp = kp[M][P*]) と Eq. 10 (Rp ∝ [M][I]) は理論的比較用の式であり、シミュレーションループでは使用されません。

---

## 反応テンプレート

### ビニルラジカル重合 (Table S1)

4 グループ構成: radical_C (i), vinyl_alpha_C (j), chain_C (k), vinyl_beta_C (l)

| ペア | 型 | 距離窓 | score_pair | 用途 |
|---|---|---|---|---|
| i-j | formation | [3, 6] A | Yes | C-C 結合形成 |
| i-k | formation | - | constraint_only | 拘束 |
| j-l | formation | - | constraint_only | 拘束 |

AIBN 活性化 (Table S1 Activation 行): azo_C - azo_N の V^d 解離。

### ナイロン 6,6 (Table S2) -- 実装済み・検証延期

4 グループ: amine_N (i), carboxyl_C (j), amine_H (k), carboxyl_OH (l)

| ペア | 型 | 距離窓 | score_pair | 用途 |
|---|---|---|---|---|
| i-j | formation | [3, 6] A | Yes | N-C 結合形成 |
| i-k | dissociation | [0, 3] A | Yes | N-H 解離 |
| j-l | dissociation | [0, 3] A | Yes | C-OH 解離 |
| k-l | formation | [0, 100] A | No | H-OH 水生成 (バイアスのみ) |

---

## 再現手順

### 小規模検証 (GPU 不要)

```bash
# S2 相当: 20 モノマー + 1 開始剤
python scripts/run_vinyl_aibn.py \
    --n-monomers 20 --n-initiators 1 --seed 42 \
    --backend orb --device cpu --no-barostat \
    --n-cycles 5 --biased-steps 2000 --unbiased-steps 500 \
    --output-dir runs/reproduce_s2

# 図の生成
python scripts/reproduce_figures.py \
    --trajectory runs/reproduce_s2/trajectory.jsonl \
    --bonds runs/reproduce_s2/bonds.jsonl \
    --target-temperature 333 \
    --output-dir runs/reproduce_s2/figures
```

### 論文スケール再現

[クイックスタートの Step 5](quickstart.md#step-5-論文スケール再現-二段パイプライン) を参照してください。

---

## 残された課題

1. **S6 スケールアップ**: 24 GB+ GPU での 200+10 (2520 原子) 実行
2. **統計的再現**: 論文は 3 independent runs を報告。現在は単一 seed
3. **ナイロン/エポキシ系**: Fig. 6 と Fig. 5 の再現に必要
4. **gamma 感度分析**: Fig. S4 の再現
5. **NPT アンサンブル**: barostat とバイアス力の共存検証
