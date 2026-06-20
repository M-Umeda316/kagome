# 作業計画書 v2（外部委託向けハンドオフ — 論文対比ギャップ分析版）

最終更新: 2026-06-13
作成者: プロジェクトオーナー側（Claude Code 経由）
対象読者: 本リポジトリの開発を引き継ぐ外部委託先のエンジニア／研究者
前提文書: `CLAUDE.md`（絶対遵守ルール）、`specs/handoff-plan.md`（初版計画）

---

## 0. この文書の目的

本書は、**論文 arXiv:2511.22874 と現行コードの事実ベースの差異分析**に基づき、残りの実装タスクを優先度順に整理したものです。

初版 `specs/handoff-plan.md` のルール（絶対遵守事項・アーキテクチャ規約・質問必須トリガー・Definition of Done）は**そのまま有効**です。本書はそれらを繰り返さず、**差異の具体的な事実・影響・対処方法**に特化します。

⚠️ 「推測」は含めていません。各ギャップには**コード上の根拠**（ファイル名・行番号）と**論文上の根拠**（arXiv HTML で確認した記述）を併記します。

---

## 1. 論文 vs コード — ギャップ一覧（事実のみ）

以下は全ソースファイルを読み、arXiv HTML 版論文と突き合わせた結果です。

### G-A: アンサンブル — NPT（論文）vs NVT/NVE（コード）

| 項目 | 事実 |
|---|---|
| **論文の記述** | arXiv HTML: "initial relaxation with a classical force field, followed by refinement using the PFP uMLIP under the NPT ensemble" — NPT を明記 |
| **コードの現状** | `src/integrators/langevin.py` は Langevin BAOAB（NVT）。`src/integrators/verlet.py` は VelocityVerlet（NVE）。**バロスタット実装は存在しない**。`src/integrators/` ディレクトリに NPT 関連コードなし |
| **影響** | 圧力制御がないため、定圧条件での密度変化や体積緩和を再現できない。重合で体積が変わる系では定量的差異が大きい |
| **論文で不明な点** | バロスタット種類（Berendsen / Parrinello-Rahman 等）、目標圧力値、圧力結合定数 — いずれも arXiv HTML に記載なし |

### G-B: 結合確認閾値 — r < 0.6×ΣvdW（論文）vs r < 1.2×r0（コード）

| 項目 | 事実 |
|---|---|
| **論文の記述** | arXiv HTML: "60% of the sum of their van der Waals radii" — これが結合判定基準 |
| **コードの現状** | `src/reactive/bonds.py:33` で `threshold_fraction: float = 1.2`。`check_outcomes()` (行 74) は `r <= threshold_fraction * pair.r0` で判定。`pair.r0` は `target_distance()` = λ×ΣvdW = 0.6×ΣvdW。実効閾値 = 1.2 × 0.6 × ΣvdW = **0.72 × ΣvdW** |
| **実行スクリプトでの値** | `scripts/run_mace_pbc.py` と `scripts/run_toy_bond_demo.py` では `threshold_fraction=1.3` を渡す → 実効 = 0.78×ΣvdW |
| **差異** | 論文: r < 0.6×ΣvdW = r0。コード: r < 1.2〜1.3×r0。コードは論文より **20〜30% 緩い** |
| **影響** | 結合形成の判定が緩くなるため、実際には未結合のペアを "formed" と誤判定するリスクがある |

### G-C: 初期速度 — Maxwell-Boltzmann（標準）vs ゼロ（コード）

| 項目 | 事実 |
|---|---|
| **論文の記述** | arXiv HTML に初期速度の生成方法は**明記なし**。ただし論文は NPT 平衡化を前提としており、MD の標準手法として Maxwell-Boltzmann 初期化が暗黙の前提 |
| **コードの現状** | `scripts/run_mace_pbc.py:61` で `velocities=np.zeros((n_atoms, 3))`。`scripts/run_toy_bond_demo.py:49` も同様。**Maxwell-Boltzmann 初期化関数は存在しない** |
| **影響** | ゼロ速度から Langevin サーモスタットで加速するため、序盤に温度スパイクが発生（`specs/figure-comparison.md` の温度プロットで確認済み）。平衡到達まで無駄なステップが必要 |

### G-D: 位置ラッピング（PBC 内座標折り返し）— 未実装

| 項目 | 事実 |
|---|---|
| **論文の記述** | PBC 使用が前提（NPT アンサンブル） |
| **コードの現状** | `src/integrators/verlet.py:54` と `src/integrators/langevin.py:55-60` で `positions += dt * velocities` のみ。**modulo 演算やラッピングなし**。`src/geometry.py` の `minimum_image()` は距離計算時の補正のみで、座標自体は折り返さない |
| **影響** | 長時間シミュレーションで原子がセル外に漂流する。距離計算は MIC で補正されるが、座標の出力や解析（密度プロット等）が不正確になる |

### G-E: Eq. 11 の指数フィッティング — 未実装

| 項目 | 事実 |
|---|---|
| **論文の記述** | arXiv HTML Eq. 11: α(t) = 1 - exp(-kp_eff·t) — 転化率の指数関数フィット。論文はこのフィットから kp_eff を抽出して異なる系を比較 |
| **コードの現状** | `src/analysis/conversion.py` は生の α = N_reacted / N_total を計算するのみ（行 14-16）。**指数フィット関数は存在しない**。`conversion_timeseries()` は離散的な step-count ベースの α(t) を返す |
| **影響** | 論文 Fig. 3-4 の定量比較（kp_eff 値の抽出）ができない。定性トレンド確認は可能 |

### G-F: 深さ分解反応密度プロット — 未接続

| 項目 | 事実 |
|---|---|
| **論文の記述** | arXiv HTML Eq. 12: ρ_rxn(z) — 界面近傍の反応密度を z 方向で分解 |
| **コードの現状** | `src/analysis/density.py` に `reaction_density_profile()` は**実装済み**。しかし `scripts/reproduce_figures.py` にはこの関数を呼ぶコードが**存在しない**（行 1-203 のどこにもインポートされていない） |
| **影響** | Fig. 5 相当の深さ分解密度プロットを自動生成できない |

### G-G: 化学系 — エチレンのみ（コード）vs 6 種+2 種（論文）

| 項目 | 事実 |
|---|---|
| **論文の記述** | ラジカル重合 6 系（methyl acrylate, methacrylate, styrene, vinyl acetate, diphenylethylene, dimethyl itaconate）+ nylon-6,6 + DGEBA/DETA epoxy on CuO |
| **コードの現状** | `scripts/_systems.py` にはエチレン (`build_ethylene_box()`) のみ。`build_template_and_groups()` は 2 グループ・1 ペアのビニル重合テンプレートのみ |
| **影響** | 論文の主要結果を再現するモノマー系がない。特に AIBN 開始剤を含むラジカル重合系、ステップ成長の nylon-6,6、4 グループの epoxy 系が欠如 |

### G-H: 反応後の原子マーキング

| 項目 | 事実 |
|---|---|
| **論文の記述** | arXiv HTML: "marked as 'reacted' and removed from the reactive set" |
| **コードの現状** | `src/workflows/polymerization.py:308-318` の `_update_groups_after_cycle()` で、confirmed formation の原子を `group.remove_atom()` で除去。**機能としては実装済み** |
| **差異** | なし — 論文と一致。✅ |

### G-I: BoostState リセットタイミング

| 項目 | 事実 |
|---|---|
| **論文の記述** | 明記なし |
| **コードの現状** | `BoostState` は `_run_biased_phase()` 冒頭で `BoostState()` として新規作成（`polymerization.py:198`）。各 biased フェーズ開始時に step=0 からリセットされる |
| **差異** | `specs/decisions.md`（2026-06-11）に「t resets each biased segment」として記録済み。✅ |

### G-J: γ 単位の未確定

| 項目 | 事実 |
|---|---|
| **論文の記述** | arXiv HTML: γ=1.0 は確認済み。単位は明記なし |
| **コードの現状** | `src/boost/tdbb.py:49-51` の `boost_amplitude()` は `gamma * t` で t は**ステップ数**（整数）。`BoostState.advance()` が `self.step += 1` した後に呼ぶ |
| **差異** | t が物理時間（fs）なら `gamma * t * timestep_fs` とすべき。現在は kcal/(mol·step)。`specs/decisions.md` に PENDING CONFIRMATION として記録済み |
| **影響** | saturation タイミングが変わる（step 250 vs step 1000）。定性的動作は同一だが定量差あり |

### G-K: unbiased ステップ数の未確認

| 項目 | 事実 |
|---|---|
| **論文の記述** | arXiv HTML で biased=2000 は確認。unbiased のステップ数は**明記なし** |
| **コードの現状** | `configs/boost/paper_faithful.yaml` で unbiased_steps=2000 を採用 |
| **差異** | 2000 は初期ノートからの値。PDF からの確認が必要 |

---

## 2. 対処済み項目（コード確認で確認）

以下は初版計画時に懸念されたが、**既に解消済み**の項目です。

| 項目 | 状態 | 根拠 |
|---|---|---|
| G5: 品質ゲートスクリプト欠如 | ✅ 解消 | `scripts/validate_configs.py`, `check_dependency_licenses.py`, `check_seed_defined.py`, `check_output_path.py` が存在。`tests/unit/test_quality_gates.py` で 15 テスト合格 |
| G6: 承認台帳の矛盾 | ✅ 解消 | `specs/approved_dependencies.yaml` が rewritten（mace, ase, orb-models, matplotlib, pytorch 登録済み） |
| T7.1: PBC 経路 | ✅ 解消 | `scripts/run_mace_pbc.py` が MACE-MP-0 + PBC で動作確認済み |
| T8.1: 結合形成の機械実証 | ✅ 解消 | `scripts/run_toy_bond_demo.py` で confirmed_formations=1 達成 |
| T9.1: 図比較記録 | ✅ 部分解消 | `specs/figure-comparison.md` 作成済み。ただし formations=0 の MACE 系は未達 |
| α(t) 分母の修正 | ✅ 解消 | `reproduce_figures.py` が `--n-reactive-sites` 引数と trajectory ヘッダの `n_reactive_sites` を使用 |
| 2 グループテンプレートの妥当性 | ✅ 記録済み | `specs/decisions.md`（2026-06-13）でビニル重合には 2 グループが正しいことを確認 |

---

## 3. 残タスク — 優先度順

> すべてのタスクは `specs/handoff-plan.md` の「2. 絶対遵守ルール」「6. 質問必須トリガー」に従うこと。

### 優先度 1: 論文忠実性の修正（科学的正確性に直結）

#### T-A: 結合確認閾値を論文準拠に修正（G-B 対処）

- **目的**: 論文の "60% of the sum of their van der Waals radii" を正しく実装する
- **現状の問題**:
  - `src/reactive/bonds.py:33` — `threshold_fraction=1.2` のデフォルト
  - `src/reactive/bonds.py:74` — `r <= threshold_fraction * pair.r0` で判定
  - `pair.r0` = 0.6×ΣvdW（`src/boost/tdbb.py:56-58`）
  - 実効閾値 = 1.2 × r0 = 0.72×ΣvdW ≠ 0.6×ΣvdW（論文）
- **対処方法**:
  1. `BondTracker.__init__()` の `threshold_fraction` デフォルトを `1.0` に変更
  2. すべての呼び出し元（`scripts/run_mace_pbc.py`, `scripts/run_toy_bond_demo.py`, 他）で `threshold_fraction` の明示指定を `1.0` に修正（または省略して新デフォルトを使用）
  3. `specs/decisions.md` に変更理由を記録: "Paper states bonding criterion as r < 0.6×ΣvdW = r0. Previous threshold_fraction=1.2 was a lenient approximation."
- **受け入れ基準**:
  - `BondTracker(threshold_fraction=1.0)` がデフォルト
  - 全実行スクリプトで `threshold_fraction=1.0` を使用（または未指定）
  - `decisions.md` に変更記録
  - 既存テストがパス（テスト内の toy 系が新閾値で動作するか確認。動かない場合は toy 系の epsilon/sigma を調整して r0 内に収まるようにする）
- **⚠️ 注意**: toy_bond_demo（T8.1）は `threshold_fraction=1.3` で confirmed_formations=1 を達成した（r=1.825Å, r0=2.04Å, 閾値=2.65Å）。`threshold_fraction=1.0` に変更すると閾値が 2.04Å になり、r=1.825 < 2.04 なので**引き続き pass する**。ただし別条件の実行では結果が変わりうるため、新デフォルトで E2E を再確認すること。
- **依存**: なし
- **論文アンカー**: "60% of the sum of their van der Waals radii"（Section 2, bond confirmation criterion）

#### T-B: NPT アンサンブルの実装（G-A 対処）

- **目的**: 論文の NPT 条件を再現するバロスタットを追加する
- **現状の問題**:
  - `src/integrators/langevin.py` — NVT のみ（温度制御あり、圧力制御なし）
  - `src/integrators/verlet.py` — NVE のみ
  - バロスタット実装なし
- **対処方法**:
  1. `src/integrators/npt.py` を新規作成
  2. Langevin BAOAB + Berendsen/Parrinello-Rahman バロスタットを実装（**バロスタット種類は論文に記載なし → `specs/decisions.md` に選択理由を記録してから実装**）
  3. `Integrator` プロトコル（`src/integrators/verlet.py:12-33`）を拡張: `pre_force()` と `post_force()` のシグネチャに `cell` パラメータを追加するか、バロスタット専用のフックを設ける
  4. `SimulationState`（`src/workflows/polymerization.py:43-50`）の `cell` を可変にする（NPT ではセルサイズが変わる）
  5. `PolymerizationWorkflow._run_biased_phase()` と `_run_unbiased_phase()` でセルサイズ更新をサポート
- **受け入れ基準**:
  - `LangevinNPTIntegrator` が存在し、`Integrator` プロトコルを満たす
  - 目標圧力・圧力結合定数が config で指定可能
  - アルゴンガスの NPT テスト: 目標温度・圧力に収束することをユニットテストで確認
  - 既存 NVT/NVE テストに影響なし
  - バロスタット種類の選択理由が `decisions.md` に記録
- **⚠️ Ask-first**: バロスタットの種類と目標圧力（論文に記載なし）はオーナーに確認してから実装
- **依存**: なし（並行作業可能）
- **論文アンカー**: "refinement using the PFP uMLIP under the NPT ensemble"
- **工数見積もり**: 中〜大。Integrator プロトコル拡張 + セルサイズ更新のワークフロー結合が必要

#### T-C: Maxwell-Boltzmann 初期速度の実装（G-C 対処）

- **目的**: MD の標準手法に従い、目標温度に対応する初期速度を生成する
- **現状の問題**:
  - `scripts/run_mace_pbc.py:61` — `velocities=np.zeros((n_atoms, 3))`
  - 初期化関数なし
- **対処方法**:
  1. `src/integrators/` または `src/utils.py` に `maxwell_boltzmann_velocities(masses, temperature_K, rng)` を追加
  2. 各原子に `v_i ~ N(0, sqrt(kT/m_i))` で速度を割り当て、重心運動量をゼロに補正
  3. `scripts/run_mace_pbc.py` と他の実行スクリプトで呼び出す
  4. 単位に注意: 速度は Å/fs、`KB` は kcal/(mol·K)、`FORCE_CONV` を使って `kT[kcal/mol] × FORCE_CONV / m[amu]` → `v²[Å²/fs²]`
- **受け入れ基準**:
  - `maxwell_boltzmann_velocities()` が存在
  - 生成された速度から計算した瞬間温度が目標温度 ±10% に収まるテスト
  - 重心運動量がゼロ（`sum(m_i * v_i)` ≈ 0）のテスト
  - 実行スクリプトでゼロ速度が MB 初期化に置き換わっている
- **依存**: なし
- **論文アンカー**: MD 標準手法（論文は明示していないが NPT 平衡化を前提）

#### T-D: 位置ラッピング（PBC 座標折り返し）の実装（G-D 対処）

- **目的**: PBC シミュレーションで原子座標をセル内に保つ
- **現状の問題**:
  - `src/integrators/verlet.py:54` — `positions += dt * velocities` のみ
  - `src/integrators/langevin.py:55-60` — 同様にラッピングなし
  - `src/geometry.py:8-21` — `minimum_image()` は距離補正のみ
- **対処方法**:
  1. `src/geometry.py` に `wrap_positions(positions, cell)` を追加: `positions = positions % box` （直方体セル）
  2. 各積分器の `pre_force()` 末尾（drift 後）でラッピングを適用
  3. または `PolymerizationWorkflow` の各ステップ末尾でラッピング
  4. ラッピングの適用箇所は**積分器内が望ましい**（位置更新の直後が物理的に正しい）
- **受け入れ基準**:
  - 1000 ステップ後にすべての原子座標が `0 <= x < box` に収まるテスト
  - 非周期（cell=None）ではラッピングしないことのテスト
  - 既存テスト（非周期系）に影響なし
- **依存**: なし
- **論文アンカー**: PBC 使用の前提

### 優先度 2: 解析・可視化の完成

#### T-E: Eq. 11 指数フィッティングの実装（G-E 対処）

- **目的**: α(t) = 1 - exp(-kp_eff·t) のフィットから kp_eff を抽出する
- **現状の問題**:
  - `src/analysis/conversion.py` は生の α(t) のみ
- **対処方法**:
  1. `src/analysis/conversion.py` に `fit_conversion_exponential(steps, alpha, timestep_fs)` を追加
  2. `scipy.optimize.curve_fit` で kp_eff をフィット（scipy はオプション依存）
  3. `scripts/reproduce_figures.py` の `plot_conversion_vs_step()` にフィット曲線の重ね描きオプションを追加
  4. scipy を `pyproject.toml` の optional dependency に追加し、`specs/dependency-license-matrix.md` に登録（BSD ライセンス）
- **受け入れ基準**:
  - `fit_conversion_exponential()` が kp_eff と R² を返す
  - 合成データ（既知の kp_eff）に対してフィットが ±5% で一致するテスト
  - `reproduce_figures.py` のプロットにフィット曲線が表示される
  - scipy のライセンスが `dependency-license-matrix.md` に記録
- **依存**: なし（ただし意味のあるフィットには formations > 0 のデータが必要 → T-A 後が実用的）
- **論文アンカー**: Eq. 11: α(t) = 1 - exp(-kp_eff·t)

#### T-F: 深さ分解密度プロットの接続（G-F 対処）

- **目的**: `reaction_density_profile()` を図生成パイプラインに接続する
- **現状の問題**:
  - `src/analysis/density.py` に関数はあるが `scripts/reproduce_figures.py` が呼ばない
- **対処方法**:
  1. `scripts/reproduce_figures.py` に `plot_density_profile()` 関数を追加
  2. `--trajectory` と `--bonds` から反応イベントの位置を抽出し `reaction_density_profile()` に渡す
  3. CLI 引数 `--cell-xy-area` と `--z-bins` を追加
  4. trajectory.jsonl にはフレーム位置が含まれるため、イベント発生ステップの位置を抽出するヘルパーが必要
- **受け入れ基準**:
  - `reproduce_figures.py` に `--bonds` 指定時かつイベントが存在する場合に `density_profile.png` が出力される
  - density.py のテストが存在し、既知の入力に対して正しいビン分布を返す
- **依存**: formations > 0 のデータが必要（T-A 後）
- **論文アンカー**: Eq. 12（arXiv HTML 番号）: ρ_rxn(z)

### 優先度 3: 化学系の拡充

#### T-G: 論文の主要反応系の追加（G-G 対処）

- **目的**: 論文で扱われている反応系のうち最低 1 つを実装し、TDBB の系横断的適用性を実証する
- **現状の問題**:
  - `scripts/_systems.py` にエチレンのみ
  - 論文は 6 種のビニル系 + nylon-6,6 + epoxy を対象
- **対処方法（段階的）**:
  1. **T-G1: ビニルモノマー + AIBN 開始剤**（最小拡張）
     - `scripts/_systems.py` に `build_vinyl_monomer_box(monomer_smiles, n_molecules, n_initiators, ...)` を追加
     - 最も単純な methyl acrylate + AIBN から開始
     - SMILES → 座標変換は RDKit（BSD ライセンス、商用可）か ASE の molecule builder を使用
     - テンプレートは現行の 2 グループ（radical C + alkene C）を流用
     - `specs/dependency-license-matrix.md` に RDKit 登録（使用する場合）
  2. **T-G2: nylon-6,6（ステップ成長）**
     - hexamethylenediamine + adipic acid の系ビルダー追加
     - テンプレートは 2 グループ（amine N + carboxylic acid C）
     - Carothers 比較（Eq. 10）の検証に使用
  3. **T-G3: epoxy/CuO（4 グループテンプレート）**
     - DGEBA + DETA + CuO 表面
     - 4 グループテンプレート: Gi(epoxy O), Gj(1° amine N), Gk(2° amine N), Gl(surface OH)
     - P = {(i,j),(i,k),(j,l)}
     - **selection.py の既存 N グループ対応機構を使用**（新規ロジック不要）
     - 表面モデル生成は ASE の surface builder を使用
- **受け入れ基準（T-G1 のみ必須、T-G2/G3 はオーナー承認後）**:
  - T-G1: methyl acrylate + AIBN の系で E2E 実行が完走し、trajectory/bonds/summary が出力される
  - 系ビルダーの分子配置が化学的に妥当（結合長が不自然でない）
  - 系選定理由が `decisions.md` に記録
- **⚠️ Ask-first**: 新規依存（RDKit 等）を追加する前にライセンス確認
- **依存**: T-A（正しい閾値）、T-C（初期速度）、T-D（ラッピング）が完了していることが望ましい
- **論文アンカー**: Section 3 Methods（対象系の列挙）、Table（reaction templates per system）

### 優先度 4: 未確定パラメータの確認

#### T-H: γ 単位の確定（G-J 対処）

- **目的**: Eq. 5 の γ の単位を PDF から確定する
- **現状の問題**:
  - `specs/decisions.md` で PENDING CONFIRMATION
  - 現行コードは kcal/(mol·step)
- **対処方法**:
  1. 論文 PDF（arXiv:2511.22874）の Section 3 Methods または Supplementary Materials から γ の単位を確認
  2. kcal/(mol·fs) の場合: `src/boost/tdbb.py:49-51` の `boost_amplitude()` で `t * timestep_fs` を使うよう変更。**この変更は Ask-first トリガー 6（単位系変更）に該当 → オーナー確認必須**
  3. `specs/decisions.md` を確定形に更新
  4. `configs/boost/paper_faithful.yaml` のコメントを更新
- **受け入れ基準**:
  - `decisions.md` の γ エントリが "PENDING" ではなく確定
  - コードと config が確定した単位系を反映
- **⚠️ Ask-first**: 単位変更はオーナー確認後のみ実施
- **依存**: 論文 PDF の入手
- **論文アンカー**: Eq. 5, Section 3 Methods

#### T-I: unbiased ステップ数の確認（G-K 対処）

- **目的**: unbiased フェーズのステップ数を PDF で確認する
- **対処方法**: T-H と同時に実施。PDF で確認し、`paper/notes.md` のテーブルと config を更新。
- **依存**: 論文 PDF の入手

---

## 4. タスク依存関係図

```
T-H, T-I (PDF確認)  ─────────────────────────────────────────────┐
                                                                  │
T-A (閾値修正) ──┐                                                │
T-C (初期速度) ──┤                                                │
T-D (ラッピング) ┤── T-G1 (ビニル+AIBN) ── T-G2 (nylon) ── T-G3 (epoxy)
T-B (NPT) ───────┘         │
                            │
T-E (Eq.11 フィット) ───────┤
T-F (密度プロット) ─────────┤
                            ↓
                    論文スケール再実行 + 図再現（specs/figure-comparison.md 更新）
```

**推奨作業順序:**

1. **T-A**（閾値修正） — 1 行変更 + テスト確認、即日完了
2. **T-C**（MB 初期速度） — 関数 1 つ追加、1 日
3. **T-D**（ラッピング） — 関数 1 つ + 積分器修正、1 日
4. **T-B**（NPT） — 中規模実装、2-3 日（Ask-first 後）
5. **T-E**（Eq. 11 フィット） — 小規模、0.5 日
6. **T-F**（密度プロット接続） — 小規模、0.5 日
7. **T-G1**（ビニル+AIBN） — 系ビルダー + E2E、2 日
8. **T-H, T-I**（PDF 確認） — PDF 入手次第、0.5 日
9. **T-G2, T-G3** — オーナー承認後

---

## 5. 各タスクに共通の作業手順

`specs/handoff-plan.md` の「7. 各タスク共通の Definition of Done」に加え、以下を遵守:

1. **コード変更前に**: 変更するファイルの現在の状態を確認。特に他タスクとの競合に注意。
2. **仮定を入れる場合**: `specs/decisions.md` にテンプレート形式で記録してから実装。
3. **テスト追加**: 新規関数には必ずテスト。`tests/unit/` にファイル追加。
4. **E2E 確認**: 変更後に `scripts/run_toy_bond_demo.py --seed 7` が引き続き `confirmed_formations >= 1` を返すことを確認（リグレッション防止）。
5. **図再生成**: `scripts/reproduce_figures.py` でプロットが壊れていないことを確認。

---

## 6. 既知の制約・注意事項

### 6.1 バックエンド差異
- 論文は PFP/Matlantis（proprietary）を使用。本実装は MACE-MP-0（MIT）がデフォルト。
- **定量的な一致は期待しない**。定性トレンド（biased/unbiased エネルギー分離、α(t) 増加）の一致が Phase 3 受け入れ基準。
- この前提は `specs/decisions.md`（2026-06-12「Using MACE instead of paper's PFP backend」）に記録済み。

### 6.2 Windows 環境
- OrbMol-v2 の周期 PME には `nvalchemiops`（NVIDIA）が必要 → `blocked_pending_review`
- Windows での `torch.compile` に互換性問題あり
- MACE-MP-0 + PBC は動作確認済み（`scripts/run_mace_pbc.py`）

### 6.3 計算コスト
- 論文スケール（2000+2000 steps × 10 cycles）は CPU で数時間〜。
- 4 分子以上の系は GPU 推奨。
- **本実行前にオーナーへコスト見積もりを共有すること**（Ask-first トリガー 7）。

---

## 7. ファイル対応表（各タスクが触るファイル）

| タスク | 変更ファイル | 新規ファイル |
|---|---|---|
| T-A | `src/reactive/bonds.py`, `scripts/run_*.py` | — |
| T-B | `src/integrators/verlet.py` (Protocol), `src/workflows/polymerization.py` | `src/integrators/npt.py` |
| T-C | `scripts/run_mace_pbc.py`, `scripts/run_toy_bond_demo.py` | `src/utils.py` (or `src/integrators/init_velocities.py`) |
| T-D | `src/integrators/verlet.py`, `src/integrators/langevin.py` | — (`src/geometry.py` に `wrap_positions()` 追加) |
| T-E | `src/analysis/conversion.py`, `scripts/reproduce_figures.py` | — |
| T-F | `scripts/reproduce_figures.py` | — |
| T-G1 | `scripts/_systems.py` | `scripts/run_vinyl_aibn.py` |
| T-G2 | `scripts/_systems.py` | `scripts/run_nylon66.py` |
| T-G3 | `scripts/_systems.py` | `scripts/run_epoxy_cuo.py` |
| T-H | `src/boost/tdbb.py` (条件付き), `specs/decisions.md`, `configs/boost/paper_faithful.yaml` | — |
| T-I | `paper/notes.md`, `configs/boost/paper_faithful.yaml` | — |

---

## 8. 質問必須トリガー（再掲 — 重要）

以下に該当したら **作業を止めてオーナー（m.umeda1502@gmail.com）に確認**すること。

1. 論文が複数の妥当な数学的解釈を許す場合
2. 商用利用権が不明なバックエンドを追加しようとする場合
3. 簡略化が TDBB の科学的意味を変える場合
4. 図再現でスムージング・フィルタ・平均化を論文と変える場合
5. 論文の引用なしにデフォルトのハイパーパラメータを導入する場合
6. 単位系（`src/units.py`）を変更する場合
7. 大規模計算を本実行する前

**T-B（NPT バロスタット種類）と T-H（γ 単位）は確実にトリガーに該当します。**

---

（本書に不明点・矛盾があれば、実装を止めてオーナーに確認すること。`specs/handoff-plan.md` の絶対遵守ルール・アーキテクチャ規約はそのまま有効です。）
