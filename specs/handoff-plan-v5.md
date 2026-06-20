# Handoff Plan v5 — コードレビュー指摘事項の修正計画

Date: 2026-06-18
Paper: arXiv:2511.22874, Mori et al. (paper/2511.22874v1.pdf)
種別: **コード品質・整合性・再現性の是正**（v4 までの「科学的スコープ拡張」とは直交）。
v4 を supersede しない。v4 = 重合再現の前進、v5 = 既存実装の不整合・責務・テスト・再現性の修正。

このドキュメントは 2026-06-18 のリポジトリ全体レビューで洗い出した指摘を、別担当者が
着手できる粒度の実装タスクに落としたもの。各タスクは独立して着手可能（依存は明記）。

---

## 実装者向けグラウンドルール（CLAUDE.md 非交渉要件の再掲）

着手前に必ず守ること:

1. **すべての変更は paper artifact（claim/式/図/表/method 段落）を1つ以上参照する。**
   各タスクに "Paper/guardrail anchor" を明記済み。無ければ追記してから着手。
2. **paper で一意に定まらない解釈・仮定は、実装前に `specs/decisions.md` に追記する**
   （テンプレートは decisions.md 冒頭）。本計画で「要 PDF 確認」と書いた項目は、
   先に PDF を読み decision を残してからコードを触ること。
3. **TDBB のパラメータ・物理は変えない**: f2=10、f1_max=250/125、γ=1.0、λ=0.6、
   r0=λ·Σr_vdw、production の候補窓 [3,6]（Table S1）。すべて paper 確定値。
   本計画の修正はいずれも「記録・整合・責務・テスト」であり TDBB の科学的意味を変えない。
   RF2/RF5 は当初「要 PDF 確認」だったが、2026-06-18 に PDF 突き合わせ済み
   （下記「PDF 突き合わせ結果」）。paper に明示された定義に合わせる修正であり、
   解釈の余地は無い（owner 承認は不要、ただし候補順位/α スケールが変わるため decision 記録は必須）。

## PDF 突き合わせ結果（2026-06-18, RF2/RF5 用）

`paper/2511.22874v1.pdf`（PyMuPDF でテキスト抽出。本文 28 ページ）を直接確認した確定事実。

### F-1 候補スコア d_ijkl（RF5 の根拠） — 本文 p.4, Eq.7
- 目標ペア集合 **P = {(i, j), (i, k), (j, l)}**（本文 p.4 そのまま）。
- 非重複選択のスコアは **dijkl = rij + rik + rjl**（**3項固定**、本文 p.4）。反応種に依らず同一式。
- 「各反応種 r は独自の (GI,GJ,GK,GL)・P(r)・距離窓を持てる」とあるが、**スコア式は上記3項のまま**。

### F-2 Table S1（vinyl, SI p.21）
全反応種が **同じ3ペア (i–j, i–k, j–l) で群同定**（距離窓を持つのはこの3ペアのみ）:

| Reaction | i–j 窓 | i–k 窓 | j–l 窓 | 適用バイアス |
|---|---|---|---|---|
| Activation | 0–3 | 0–3 | 0–3 | V^d を **i–k, j–l** |
| Initiation | 3–6 | 0–3 | 0–3 | V^f を **i–j** |
| Propagation| 3–6 | 0–3 | 0–3 | V^f を **i–j** |

→ vinyl の Initiation/Propagation は i–j のみ V^f、i–k/j–l は窓拘束のみ（バイアス無し）。
   これは現行 `_systems.build_vinyl_aibn_system`（i–j 形成 [3,6]、i–k/j–l constraint_only [0,3]）と**完全一致**。
   `score_candidates` が全 pair を合算しても、vinyl は3ペアなので d_ijkl と一致する（=現状 vinyl は正しい）。

### F-3 Table S2（nylon-6,6 Condensation, SI p.22） — RF5 の核心
群同定は **i–j (3–6), i–k (0–3), j–l (0–3) の3ペアのみ**。**k–l ペアには距離窓が無い**。

| Reaction | i–j | i–k | j–l | 適用バイアス |
|---|---|---|---|---|
| Condensation | 3–6 | 0–3 | 0–3 | V^f を **i–j, k–l** / V^d を **i–k, j–l** |

割り当て（現行コードと照合）: i=amine_N, j=carboxyl_C, k=amine_H, l=carboxyl_OH。
- i–j: アミド結合形成 (V^f) … 現行コード一致
- i–k: N–H 解離 (V^d) … 一致
- j–l: C–OH 解離 (V^d) … 一致
- **k–l: H–OH（水）形成 (V^f)。群同定・距離窓・スコアには非関与、バイアスのみ** … **現行コードはこれをスコアに含めてしまっている（不整合）**

→ **結論(RF5)**: スコア d_ijkl は常に i–j, i–k, j–l の3項。nylon の k–l は**バイアス専用**でスコア・候補同定に入れてはならない。
   現行 `score_candidates`（`selection.py:116-123`）が template.pairs 全合算 → nylon で r_kl を余計に加算。要修正。

### F-4 monomer conversion α（RF2 の根拠） — 本文 p.9 Fig.2 キャプション, Eq.11
- **α = 1 − [M]/[M]₀**（[M]₀=初期モノマー濃度, [M]=瞬時モノマー濃度）。Eq.11 α(t)=1−exp(−k*_p·t) の α は「monomer conversion」。
- すなわち **分母 = 初期モノマー数 n_monomers**、分子 = 反応済みモノマー数。
  vinyl では1形成イベント(radical_C+vinyl_alpha_C)につきモノマー1個が消費されるので
  **α = confirmed_formations / n_monomers**。
- 開始剤・脱離基・constraint 群は分母に**含めない**。
- nylon（step-growth）は monomer conversion ではなく**反応進行度 p**（Carothers, DPn=1/(1−p)）で別管理
  （`src/analysis/carothers.py`）。Eq.11 の α プロット（`reproduce_figures.py`）は vinyl の monomer-conversion 用。

### 補足（RF5 と隣接、今回の修正対象外）
Table S1 の **Activation は 4群で V^d を i–k と j–l に適用**するが、現行 `build_activation_template`
は別の 2群 (azo_C, azo_N) 簡略テンプレートを使う（decisions.md 2026-06-18 で OrbMol-v2 向けに合意済の意図的簡略化）。
RF5 のスコア修正とは独立。必要なら別タスクで Table S1 準拠の4群 activation を検討。
4. **決定論**: seed を固定。リファクタは「ビット一致」を原則の受入基準にする
   （toy backend の既存テストが pos1==pos2 を要求している。`tests/unit/test_workflow.py::test_deterministic_with_seed`）。
5. 完了の定義: コード実装 + テスト緑 + 仮定の文書化 + 再現コマンド + 出力が config/seed に辿れる。
6. 各タスク完了時に `.claude/hooks/pre-commit.sh`（validate_configs → check_dependency_licenses → pytest）
   を通すこと。

---

## 優先度つきタスク一覧

| ID | 指摘 | 区分 | 依存 | 工数 | リスク | 科学的意味 |
|----|------|------|------|------|--------|-----------|
| RF1 | 設定ファイルが実行に読まれず manifest の config_path と実パラメータが乖離 | 再現性 | — | 中 | 中 | なし（記録のみ） |
| RF2 | α(t) の分母（reactive sites）が3箇所で不一致＋死蔵コード | 整合性 | — | 中 | 中 | **PDF確認済**: α=1−[M]/[M]₀ → 分母=n_monomers（p.9 Fig.2） |
| RF3 | in-bias 検出 / activation / equilibration のユニットテスト欠落 | テスト | — | 中 | 低 | なし |
| RF4 | `PolymerizationWorkflow` の MD ループ重複＋群更新の vinyl 固有直書き | 責務/設計 | RF3 | 大 | 中 | なし（不変リファクタ） |
| RF5 | `score_candidates` が nylon で k–l 項を余計に加算（d_ijkl 不一致） | 整合性 | — | 中 | 中 | **PDF確認済**: d_ijkl=r_ij+r_ik+r_jl の3項固定（Table S1/S2） |
| RF6 | notes.md / 旧 decision の陳腐化（4群化を反映せず supersede 未明示） | ドキュメント | — | 小 | 低 | なし |
| RF7 | run-until-reaction 早期終了時に `CycleLog.steps` が全ステップを記録 | 整合性 | — | 小 | 低 | なし |
| RF8 | activation の解離判定だけ絶対閾値 2.5 Å（規約不統一） | 整合性 | RF3 | 小 | 低 | 弱くあり→要記録 |
| RF9 | ヘルパ/定数の重複（_accel, EV_TO_KCAL_MOL, 温度計算） | 整理 | — | 小 | 低 | なし |
| RF10 | `Candidate.pair_distances` の死蔵 | 整理 | RF5 | 小 | 低 | なし |
| RF11 | パッケージング不整合（配布名 kagome / import は src.*、egg-info 2箇所） | 構造 | — | 中 | 中 | なし |
| RF12 | 幾何が orthorhombic 限定（epoxy/スラブ系の制約） | 文書化 | — | 小 | 低 | なし（記録のみ） |
| RF13 | NVE+barostat 時の温度フォールバック 300 K 固定 | 整合性 | — | 小 | 低 | 弱くあり |

**ID 番号と優先度の関係**: ID（RF1…RF13）は **優先度グループ順**（High=RF1–2 → Medium=RF3–6 →
Low=RF7–13）に並べてあるが、**厳密な実行順ではない**。実際の着手順は下記「推奨着手順」を正とする。

**推奨着手順**: RF6 → RF3 →（RF7, RF9, RF13）→ RF1 → RF2 → RF5 → RF4 →（RF10）→ RF11 → RF12。
RF6/RF3/RF9 は低リスクの足場固め。RF4 は RF3（テスト）で安全網を張ってから。
RF2/RF5 は PDF 突き合わせ済み（上記「PDF 突き合わせ結果」）。owner 承認は不要だが、
α スケール・候補順位が変わるため decisions.md への記録は必須。

---

# RF1 — 設定ファイルを実読込し、manifest に実効パラメータを記録

## Finding（事実）
- `configs/boost/paper_faithful.yaml`（キー `f2_invA2`, `f1_max_form_kcal_mol`,
  `f1_max_break_kcal_mol`, `lambda_vdw`, `gamma`, `timestep_fs`, `biased_steps`,
  `unbiased_steps`）は **どの実行経路でも `TDBBParams`/`PolymerizationConfig` に読み込まれていない**。
- `scripts/run_vinyl_aibn.py:302-319` は CLI 引数から config を構築し、
  `config_path='configs/boost/paper_faithful.yaml'`（同:410）は `RunManifest`
  への**文字列としてのみ**渡る（`src/workflows/manifest.py`、`polymerization.py:159-166`）。
- `scripts/validate_configs.py:23-34` は YAML キーの存在/正値を検査するだけで値を消費しない。
- 真実源が3系統に分散しドリフトしている: ① `TDBBParams` 既定（`src/boost/tdbb.py:22-26`）
  ② YAML（`f2_invA2: 10.0`）③ argparse 既定（`run_vinyl_aibn.py:109` の `--f2` 既定 10.0）。
  S4 検証は `f2=5.0` を採用（`specs/decisions.md` 2026-06-17 系記録）。

## Goal
manifest を見れば**実際に使われた TDBB パラメータが復元できる**状態にする。
（CLAUDE.md「outputs can be traced to configs and seeds」「record … config path」）

## Paper/guardrail anchor
CLAUDE.md 非交渉要件「All experiments must record seed, config path, git SHA,
backend name, and output directory」。configs/boost/paper_faithful.yaml の各値は
PDF p.7 由来（同 YAML のコメント参照）。

## 推奨アプローチ（2案、決定を decisions.md に記録すること）
- **案A（推奨・小）**: YAML を実読込せず、`RunManifest.extra` に**実効パラメータの dict**
  （f2, gamma, f1_max_formation, f1_max_dissociation, lambda_vdw, timestep_fs,
  biased_steps, unbiased_steps, n_cycles, 候補窓 r_min/r_max, backend, spin）を必ず格納する。
  config_path は「由来の目安」と位置づけ、真の記録は extra に置く。乖離が原理的に消える。
- **案B（中）**: `configs/boost/*.yaml` を読み `TDBBParams`/`PolymerizationConfig` を構築する
  ローダ（例 `src/io/config.py`）を新設し、run script はそれを土台に CLI で上書き。
  YAML キー（`f2_invA2` 等）→ dataclass フィールド（`f2` 等）の対応表を1箇所に集約。
  上書きが起きたら manifest に「YAML 値 + override 差分」を残す。

どちらでも可。案Aを既定推奨（変更面積が小さく、即座に乖離を解消できる）。

## Steps（案A）
1. `src/workflows/manifest.py`: `RunManifest` に実効パラメータを受け取り `extra` に格納する
   経路を用意（既に `extra: dict` フィールドあり）。
2. `src/workflows/polymerization.py:159-166`: `run()` で manifest を作る箇所に、
   `self.config`（TDBBParams 含む）と候補窓・backend・spin をシリアライズして渡す。
   ヘルパ `PolymerizationConfig`→dict を1つ作る（`dataclasses.asdict` で可、TDBBParams 入れ子も展開）。
3. 全 run script（`run_vinyl_aibn.py`, `run_mace.py`, `run_mace_pbc.py`, `run_orb.py`,
   `run_nylon66.py`, `run_smoke.py`）が manifest に実効値を載せることを確認（ほぼ自動）。
4. `configs/boost/paper_faithful.yaml` の冒頭コメントに「このファイルは現状 run script に
   自動適用されない。実効値は manifest.extra を見よ（または案B 採用時はローダ経由）」と明記。
   ※コメント追記も「設定変更」なので decision を残す。

## Tests
- `tests/unit/test_manifest.py`（新規）: `RunManifest` を保存→読み戻し、`extra` に
  f2/gamma/f1_max 等が含まれることを assert。
- `tests/unit/test_workflow.py`: `wf.run(state, output_dir=tmp_path, config_path=...)` 後、
  `manifest.json` を読み、TDBB 実効値が記録されていることを assert。

## Acceptance criteria
- `runs/.../manifest.json` から f2, gamma, f1_max_formation, f1_max_dissociation,
  lambda_vdw, timestep, biased/unbiased steps, n_cycles, 候補窓, backend, spin が復元可能。
- 既存の決定論テストがビット一致のまま緑（manifest 追記は数値計算に影響しない）。
- decision record（案A/B のどちらを採ったか、理由）が decisions.md にある。

## Risk
- 低。記録の追加のみで物理に触れない。`extra` に numpy 型を入れると JSON 化で失敗するため
  float/int/str に正規化すること（`asdict` 後に `float(...)` 変換）。

---

# RF2 — α(t) の分母（reactive sites）定義を一本化し死蔵コードを除去

## Finding（事実）
コンバージョン α = N_reacted / N_total の分母が3箇所で食い違う:
- `src/workflows/polymerization.py:168`（**軌跡ヘッダに採用**）:
  `n_reactive_sites = sum(len(g.atom_indices) for g in self.groups.values())`。
  vinyl 4群では radical_C + vinyl_alpha_C + **chain_C + vinyl_beta_C** まで合算し、
  constraint 専用群を分母に含めて**過大計上**。
- `scripts/run_vinyl_aibn.py:461`（図生成コマンドの提案 `--n-reactive-sites`）:
  `n_reactive = args.n_monomers * 2 + args.n_initiators`。
- `scripts/run_vinyl_aibn.py:417`: `len(radical_C)+len(vinyl_alpha_C)` を計算するが
  **summary に未使用＝死蔵**。
- `scripts/reproduce_figures.py:339-353` は CLI 値（=461 の式）をヘッダ値より優先するため、
  保存値と図の分母が別物になる。

## Goal
α の分母を **paper 定義（= 初期モノマー数 n_monomers）に一本化**し、
ヘッダ・summary・図・docstring を統一。死蔵コードを除去。

## Paper/guardrail anchor（PDF 確認済 → 上記 F-4）
- 本文 p.9 Fig.2 キャプション: **α = 1 − [M]/[M]₀**（[M]₀=初期モノマー濃度）。Eq.11 の α は monomer conversion。
- よって **分母 = n_monomers**（初期モノマー数）、分子 = 反応済みモノマー数。
  vinyl では1形成イベントにつきモノマー1個消費 → **α = confirmed_formations / n_monomers**。
- 開始剤・脱離基・constraint 群（chain_C, vinyl_beta_C）は分母に**含めない**。
- nylon は monomer conversion ではなく反応進行度 p（Carothers, `carothers.py`）で別管理。

## Steps
1. decisions.md に「α denominator = n_monomers（本文 p.9 Fig.2, α=1−[M]/[M]₀）」decision を追加。
   旧来の「全群和」「2*n_mono+n_init」が誤りであった旨も記録。
2. `src/analysis/conversion.py` に分母を返す単一関数を新設
   （例 `monomer_site_count(n_monomers) -> int` か、群から数える場合は
   **vinyl_alpha_C 群のサイズ = 初期モノマー数**を使う。constraint 群と radical_C を除外）。
   ※ vinyl では「初期 vinyl_alpha_C 群サイズ = n_monomers」なので、これを唯一の定義に採用するのが
     系非依存で安全（プロパゲーションで群が減る前の初期値を使うこと）。
3. `polymerization.py:168` の `sum(... for g in self.groups.values())` を上記定義に置換し、
   TrajectoryWriter ヘッダ `n_reactive_sites` が n_monomers になるようにする。
   ※ 初期群サイズが必要なので、`run()` 冒頭（群更新前）で算出すること。
4. `run_vinyl_aibn.py:417` の死蔵変数 `n_reactive_sites` を削除。
   図提案コマンド（同:461-470）の `--n-reactive-sites` をヘッダと同じ n_monomers に統一
   （または CLI を省略してヘッダ自動採用を促す）。
5. `reproduce_figures.py` のフォールバック順（CLI > header > n_atoms）は維持。
   提案コマンドが header と一致するので乖離が消える。conversion_timeseries の docstring を
   「分母=初期モノマー数」に更新。

## Tests
- `tests/unit/test_conversion.py`: 新分母関数が vinyl 系で n_monomers を返すこと
  （constraint 群・radical_C を除外、プロパゲーション後も初期値を保つこと）を assert。
- `tests/unit/test_workflow.py::test_trajectory_output` を拡張し、ヘッダ `n_reactive_sites`
  == n_monomers（=初期 vinyl_alpha_C 群サイズ）であることを assert。

## Acceptance criteria
- ヘッダ・summary・図提案コマンド・docstring の分母が **n_monomers** に一致。
- `run_vinyl_aibn.py` に未使用の reactive-sites 変数が無い。
- decisions.md に分母定義の根拠（本文 p.9 Fig.2, α=1−[M]/[M]₀）がある。

## Risk
- 中。α の絶対スケールが変わる（従来は分母過大で α を過小評価）。過去 run の図と数値が変わるため、
  「定義変更により α スケールが変化（旧=全群和/2*n_mono+n_init → 新=n_monomers）」を
  decisions.md と figure-comparison.md に明記。

---

# RF3 — in-bias 検出 / activation / equilibration のユニットテスト追加

## Finding（事実）
`tests/` 全走査で、以下に**ユニットテストが無い**（grep 済み）:
- `BondTracker.check_reactions_during_bias`（run-until-reaction の in-bias 検出、`src/reactive/bonds.py:45-83`）
- `PolymerizationWorkflow.run_activation` および
  `scripts/_systems.build_full_aibn_system` / `build_activation_template`
- `PolymerizationWorkflow._run_equilibration_phase`
- 形成＋解離が混在する nylon テンプレートの選択/バイアス経路

これらは S3/S4 で入った paper-fidelity の中核機能。CLAUDE.md は
「Treat TDBB equations as the most paper-critical」「Make reaction candidate
generation testable without running MD」を要求。

## Goal
最近の重要機能に MD 不要・toy backend ベースの回帰テストを張る（RF4 リファクタの安全網）。

## Paper/guardrail anchor
CLAUDE.md scientific rules（候補生成の testability、biased/unbiased ログ、
形成/解離 API の区別）。paper §2.2 step 3（run-until-reaction）。Table S1（Activation, V^d）。

## Steps（新規 `tests/unit/test_bonds.py` 拡充 or 追加ファイル）
1. **check_reactions_during_bias**:
   - formation pair で r ≤ threshold·r0 のとき1件返し `_reacted` に登録、
     2回目は重複検出しないこと。
   - dissociation pair で r > threshold·r0 のとき検出すること。
   - `check_outcomes` が in-bias 検出済みペアを二重計上しないこと（`bonds.py:114-116`）。
2. **run_activation**（toy backend で MD を短く回す or 距離判定部分を分離してテスト）:
   - `build_activation_template` が azo C-N ペアのみ（nitrile 除外）を作ること
     （`_systems.build_full_aibn_system` の `_find_aibn_azo_bonds` 経由）。
   - 解離検出（`polymerization.py:636` の `r > dissoc_threshold`）が想定どおり発火し
     `(C_idx, N_idx)` を返すこと。MD を避けたい場合は判定ロジックを純関数に切り出して検証（RF8 と連携）。
3. **_run_equilibration_phase**: toy backend + 数ステップで、frames が `phase='equilibration'`,
   `cycle=-1` で書かれること、bond_tracker が呼ばれないこと。
4. **nylon mixed bias**: `build_nylon66_system` のテンプレートで `find_candidates`→
   `_build_pair_biases` が formation/dissociation の両 PairBias を生成すること（is_formation 区別）。

## Acceptance criteria
- 上記4領域に最低1テストずつ、`pytest -q tests/unit` 緑。
- RF4 着手前にこれらが存在すること（安全網）。

## Risk
- 低。テスト追加のみ。run_activation は MD を含むため、判定部の純関数化（RF8）と
  合わせると速く決定論的にできる。

---

# RF4 — MD ループの共通化と群更新のストラテジ化（責務分離）

## Finding（事実）
- `PolymerizationWorkflow`（`src/workflows/polymerization.py`）が
  最小化・平衡化・biased・unbiased・群更新・連鎖伝播・AIBN 活性化を1クラスで担う。
- **MD ステップ骨格（pre_force → calculator.compute →（total_bias）→ post_force →
  barostat → write_frame）が4箇所で重複**:
  `_run_biased_phase`(290-411), `_run_unbiased_phase`(413-473),
  `_run_equilibration_phase`(231-288), `run_activation`(517-651)。
  `run_activation` は `_build_pair_biases` 相当（555-574）まで再実装。
- 群更新 `_update_groups_after_cycle`(475-515) が群名 `'chain_C'`/`'vinyl_beta_C'`/
  `'radical_C'` をリテラル参照＝**汎用クラスに vinyl 連鎖重合の固有規則が直書き**。
  nylon/epoxy では `groups.get(...)` が None で無害化されるだけで、機構別の差し替えができない。

## Goal
不変リファクタで重複を1つの MD ループに集約し、サイクル後の群更新を**注入可能**にする。
**数値結果はビット一致**を維持。

## Paper/guardrail anchor
CLAUDE.md「Separate numerical kernels from orchestration code」「Functions should
align with scientific concepts」。物理は不変なので paper 解釈変更なし（decision に「不変リファクタ」と明記）。

## Steps
1. RF3 のテストが緑であることを確認（安全網）。
2. **MD ステップ抽出**: 1ステップを担う内部メソッド/関数を定義。例:
   `_md_step(state, forces, dt, rng, *, bias=None, boost=None, tdbb=None, barostat_on)`
   が「pre_force → compute(base) →（bias 有れば total_bias）→ post_force → barostat 試行 →
   (energy, forces) 返却」を行う。フレーム書込みと反応検出は呼び出し側に残す
   （フェーズ固有のため）。biased/unbiased/equilibration/activation を同一ステップ関数で表現。
3. **力の構築順を厳密に保存**: 現行は「pre_force は前ステップの current_forces を使用 →
   位置更新後に base+bias を再計算 → post_force」。barostat 受理時の bias 再計算
   （`polymerization.py:357-371`）も含めて順序・丸めを変えないこと（ビット一致条件）。
4. **群更新のストラテジ化**: `PostCycleUpdater` プロトコル（例
   `update(groups, confirmed_formations, state) -> None`）を導入し、vinyl 連鎖伝播ロジックを
   `VinylChainPropagationUpdater` に移設。`PolymerizationWorkflow` は updater を DI で受け取り、
   未指定時は「形成原子を群から除去するだけ」の no-op 既定。
   既存の `propagation_map`/`chain_c_map` は updater 側へ移す。
5. `run_activation` をステップ関数の上で再構成（解離判定は RF8 で純関数化したものを使用）。
6. VDW_RADII / ATOMIC_MASSES（`polymerization.py:35-43`）を `src/reactive/` か
   `src/constants.py`（新規）へ移し、workflow から参照（責務整理。RF9 と統合可）。

## Tests
- `tests/unit/test_workflow.py::test_deterministic_with_seed` がビット一致のまま緑。
- 既存の `TestChainPropagation`（test_workflow.py:194-281）が
  `VinylChainPropagationUpdater` 経由でも同じ結果になるよう移植/再配線。
- 反応イベントで break するケースのステップ数（RF7）も確認。

## Acceptance criteria
- MD ステップ実装が1箇所に集約され、4フェーズがそれを呼ぶ。
- 群更新が DI 可能で、`PolymerizationWorkflow` 本体から vinyl 固有の群名リテラルが消える。
- 全 run script が従来どおり動作し、toy 決定論テストがビット一致。
- decision record（不変リファクタである旨、ビット一致を確認した手順）あり。

## Risk
- 中。力の評価順・barostat 受理時の再計算順を1つでも変えると数値が変わる。
  小さなコミットに分割し、各コミットで決定論テストを回すこと。

---

# RF5 — スコア d_ijkl から nylon の k–l 項を除外（PDF 確認済）

## Finding（事実）
`src/reactive/selection.py:108-126` の `score_candidates` は `constraint_only` を
区別せず **template.pairs 全部の距離を合算**する。
- vinyl（pairs (i,j),(i,k),(j,l)）: `r_ij+r_ik+r_jl` = paper の d_ijkl と一致（**現状正しい**）。
- nylon（pairs 4本、`_systems.py:569-586`）: `r_ij+r_ik+r_jl+r_kl` と**4項**になり、paper の3項定義と不一致。
  かつ (amine_H=k, carboxyl_OH=l) は `r_max=100`（同:582-585）でほぼ無拘束のため、
  **本来スコアに入らない H–OH 距離が候補ソートを支配し得る**。

## PDF 確認結果（上記 F-1〜F-3）— これにより本タスクは「確定した修正」に確定
- 本文 p.4 Eq.7: スコアは **d_ijkl = r_ij + r_ik + r_jl の3項固定**（全反応種共通）。
- Table S1/S2: 群同定（距離窓を持つペア）は **常に i–j, i–k, j–l の3ペアのみ**。
- nylon の **k–l ペアは「バイアス適用専用」**（Table S2: V^f を k–l に適用）で、
  **群同定・距離窓・スコアには非関与**。→ k–l をスコア/候補同定から外すのが paper 準拠。
- nylon の i–k, j–l は V^d バイアスを持つが、**群同定ペアなのでスコアには入る**（現状の合算で正しい）。
  ⇒ **修正対象は「k–l をスコアと候補同定から外す」一点のみ**。i–j/i–k/j–l の扱いは現状維持。

## Goal
スコア d_ijkl が全系で **r_ij + r_ik + r_jl の3項**になるよう、k–l のような「バイアス専用ペア」を
スコア・候補同定から除外する。バイアス適用（V^f on k–l）は維持する。

## Paper/guardrail anchor
本文 p.4 Eq.7（d_ijkl=r_ij+r_ik+r_jl）、SI Table S1（vinyl）, Table S2（nylon, k–l はバイアスのみ）。
Ask-first ではない（paper に明示）。ただし nylon の候補順位が変わるため decision 記録は必須。

## Steps（最小・外科的）
1. decisions.md に「scoring/identification pair set = {(i,j),(i,k),(j,l)} のみ。nylon の k–l は
   バイアス専用でスコア・候補同定に非関与（本文 p.4 Eq.7 / Table S2）」decision を追加。
2. `src/reactive/groups.py::PairSpec` に **`score_pair: bool = True`**（命名は任意。
   「群同定＝距離窓フィルタ＋d_ijkl スコアに参加するか」を表すフラグ）を追加。
   - 既存 `constraint_only`（=バイアス無しだが同定/スコアに参加）はそのまま温存。
   - 新フラグは独立軸: `score_pair=False` は「同定/スコアに非関与」を意味し、バイアス適用は
     `is_formation`/`constraint_only` 側で従来どおり制御。
3. `src/reactive/selection.py`:
   - `find_candidates` の `pair_specs` 構築（49-53行）を **`score_pair=True` のペアのみ**に限定
     （k–l の距離窓 [0,100] による疑似拘束を排除＝より忠実）。
   - `score_candidates`（116-123行）を **`score_pair=True` のペアのみ合算**に修正
     → 全系で d_ijkl = r_ij+r_ik+r_jl。
4. `scripts/_systems.build_nylon66_system` の k–l PairSpec（(amine_H, carboxyl_OH) 形成,
   `_systems.py:582-585`）に **`score_pair=False`** を設定（is_formation=True は維持＝V^f は適用され続ける）。
5. `_build_pair_biases`（`polymerization.py`）は変更不要（k–l は依然 is_formation=True で V^f が掛かる）。
   ただし k–l の atom が候補タプル内で確定していることに依存するので、k,l が i–k/j–l 経由で
   候補に含まれることを確認（含まれる: l は j–l、k は i–k で同定される）。
6. `find_candidates`/`score_candidates`/`PairSpec` の docstring に
   「距離窓フィルタ＋スコア（群同定ペア）と、バイアス適用は別概念」を明記。

## Tests
- `tests/unit/test_selection.py`:
  - vinyl: score == r_ij+r_ik+r_jl（不変）を維持。
  - nylon テンプレート（`build_nylon66_system`）で score が **3項のみ**（k–l 非加算）になることを assert。
  - `score_pair=False` のペアが候補同定（find_candidates の距離窓）に影響しないことを assert。
- `tests/unit/test_systems.py`: nylon テンプレートの k–l PairSpec が `score_pair=False` かつ
  `is_formation=True`（V^f 維持）であることを assert。

## Acceptance criteria
- 全系でスコア = r_ij+r_ik+r_jl（3項）。nylon で k–l がスコア・候補同定に入らない。
- nylon の k–l に V^f バイアスは従来どおり適用される（重合化学は不変）。
- decisions.md に PDF 参照（p.4 Eq.7 / Table S2）つきの根拠がある。

## Risk
- 中。nylon の候補順位が変わり選択＝反応経路が変わり得る（vinyl は不変＝決定論テストに影響なし）。
  既存の nylon 検証 run があれば1回再現確認。RF10（`pair_distances` 死蔵）と同時に整理可。

---

# RF6 — ドキュメント/旧 decision の陳腐化解消（4群化の反映・supersede 明示）

## Finding（事実）
- `paper/notes.md:33-37` の反応テンプレート表は
  「Vinyl … 2: Gi/Gj … {(i,j)} … Matches current 2-group implementation」と記載。
  実装は4群（`scripts/_systems.py:413-441`、git: "S4 Phase 1 … d_ijkl per Table S1"）。
- `specs/decisions.md` の「2026-06-13: 2-group template is correct for vinyl」（249-256行）は、
  後続「2026-06-18: S4 Phase1 multi-pair（4群）」（661 行〜）で事実上覆っているが、
  **旧記録に superseded 表記が無い**。
- `paper/claims.yaml` と src 内 docstring の「Eq.11/12 番号は PDF 未確認」は未解決のまま
  （`claims.yaml:40-71`、`src/analysis/conversion.py`・`density.py` の docstring）。

## Goal
ドキュメントを現実装（4群）に追従させ、覆った decision に supersede を明示。

## Paper/guardrail anchor
CLAUDE.md「Figures/notes must be reproducible … store only short notes」、
意思決定の追跡可能性。Table S1（vinyl 4群、d_ijkl）。

## Steps
1. `paper/notes.md` の表を 4群（radical_C, vinyl_alpha_C, chain_C, vinyl_beta_C、
   P={(i,j),(i,k),(j,l)}）に更新。「2-group」記述を削除し、Table S1 参照に。
2. `decisions.md` の 249-256 行 decision 冒頭に
   「**Superseded by 2026-06-18 S4 Phase1（4-group multi-pair）**」を1行追記
   （履歴は残す）。
3. Eq.11/12 番号: PDF を1回読んで確定し、`claims.yaml` と
   `src/analysis/conversion.py`・`density.py` の docstring を統一
   （確定できなければ「PDF p.X で確認、番号は HTML と一致/相違」を1行で記録）。

## Acceptance criteria
- notes.md の vinyl 行が4群。covered decision に supersede 注記。
- Eq 番号が claims.yaml と analysis docstring で一致（または相違理由が1箇所に集約）。

## Risk
- 低。文書のみ。コピー量は最小限（CLAUDE.md「short notes only」を順守、論文本文の長文転載禁止）。

---

# RF7 — run-until-reaction 早期終了時の `CycleLog.steps` を実ステップに

## Finding（事実）
`src/workflows/polymerization.py:392-402` で反応イベント検出時に `break` するが、
戻り値（404-411）は常に `steps=self.config.biased_steps`。早期終了しても全ステップ走った
体になり、`summary.json` のログ（`run_vinyl_aibn.py:441-453`）に実ステップが残らない。

## Goal
実際に回った biased ステップ数を記録。

## Steps
1. ループ内で到達 step を保持（`steps_run = step_in_phase + 1`、break 時/完走時とも更新）。
2. `CycleLog(... steps=steps_run ...)` を返す。
3. 早期終了したか（reaction で break したか）を示す bool を `CycleLog` に追加すると
   診断に有用（任意）。

## Tests
- `tests/unit/test_workflow.py`: bond_tracker 付き・近接初期配置で1サイクル走らせ、
  break した場合 `logs[0].steps < config.biased_steps` を assert。

## Acceptance criteria
- 早期終了時 `CycleLog.steps` が実ステップ＝ break 時の `step_in_phase+1`。
- 完走時は従来どおり `biased_steps`。

## Risk
- 低。ログ値のみ。物理・決定論に影響なし。

---

# RF8 — activation の解離判定規約を r0 相対へ統一（または純関数化）

## Finding（事実）
`run_activation` は固定の絶対閾値 `dissoc_threshold = 2.5`（`polymerization.py:603,636`）を使う。
一方 `BondTracker` の形成/解離判定は `threshold_fraction * r0` の相対閾値
（`src/reactive/bonds.py:69,123,135`）。同じ「解離」でモジュール間の閾値規約が不統一。

## Goal
解離判定の規約を1つに統一し、判定ロジックを MD から分離してテスト可能にする。

## Paper/guardrail anchor
Eq.3（V^d, r0 を含まない）。decisions.md「2026-06-12: Dissociation tracking uses
r0 = λ·Σr_vdW as confirmation threshold」既存記録あり → activation もこれに揃えるか、
azo C-N に固有の閾値を使う理由を decisions.md に記録。

## Steps
1. 解離判定を純関数 `is_dissociated(r, r0, threshold_fraction) -> bool` に切り出し、
   `BondTracker` と `run_activation` の両方から使う。
2. activation の 2.5 Å を、azo C-N の r0（λ·Σr_vdw(C,N)）相対に置換するか、
   絶対値を使う科学的理由（C-N 平衡長 ~1.49 Å、解離は >2.5 Å 等）を decisions.md に明記して残す。
   ※ activation パラメータ（f2=0.3, f1_max=250）は decisions.md 2026-06-18 で確定済みなので変えない。
3. RF3 のテストでこの純関数を検証。

## Acceptance criteria
- 解離判定が1つの純関数に集約され、両経路から使われる。
- activation の閾値根拠が decisions.md にある。

## Risk
- 低〜小。閾値を相対化すると activation の解離検出タイミングが微変動し得るため、
  既存 activation 検証（`runs/s4_activation_v3/`）の再現を1回確認。

---

# RF9 — 重複ヘルパ/定数の集約（DRY）

## Finding（事実）
- `_accel` が `src/integrators/verlet.py:70-79`（staticmethod）と
  `src/integrators/langevin.py:76-84`（module 関数）に同一実装で重複。
- `EV_TO_KCAL_MOL = 23.060548` が `src/backends/ase_adapter.py:18` と
  `src/backends/orb_backend.py:22` に重複定義。
- 温度計算が `_instant_temperature`（`polymerization.py:96-113`）と
  `instant_temperature_K`（`init_velocities.py:46-59`）に同式で二重実装。

## Goal
それぞれ単一ソースに集約。

## Steps
1. `_accel` を `src/integrators/_common.py`（新規）か `src/units.py` に
   `force_to_accel(forces, masses)` として置き、両 integrator が import。
2. `EV_TO_KCAL_MOL` を `src/units.py` に1つ定義し、両 backend が import
   （`src/units.py` には既に他の換算定数あり）。
3. 温度計算を1関数（masses=None を許容する版）に統一し、polymerization から import。

## Tests
- 既存 `tests/unit/test_integrators.py`・`test_init_velocities.py`・`test_backends.py`
  が緑のまま（数値不変）。

## Acceptance criteria
- 各重複が1定義に集約され、import で共有。重複定義が grep で出ない。

## Risk
- 低。値は同一なので数値不変。import 経路の循環に注意（units は末端依存に置く）。

---

# RF10 — 死蔵フィールド `Candidate.pair_distances` の整理

## Finding（事実）
`src/reactive/selection.py:75-78,100` で `Candidate.pair_distances` を埋めるが、
`score_candidates`（同:116-123）は距離を再計算し格納値を使わない。下流でほぼ未使用
（テストは `{}` を渡す: `test_selection.py:113-135`）。計算重複かつ未使用フィールド。

## Goal
死蔵を解消（除去 or 活用）。

## Steps
- RF5 と同時実施推奨。`score_candidates` を `pair_distances` を活用する形にして二重計算を消すか、
  `pair_distances` 自体を削除して `Candidate` を軽量化する（どちらかを decision で選択）。
- 削除する場合は `Candidate(atom_indices, score)` に簡約し、生成箇所/テストを更新。

## Acceptance criteria
- 距離の二重計算が無い、または未使用フィールドが無い（いずれか一貫した状態）。
- `test_selection.py` 緑。

## Risk
- 低。RF5 と絡むため、RF5 のスコア定義確定後に着手すると手戻りが無い。

---

# RF11 — パッケージング構造の整理

## Finding（事実）
- 配布名は `kagome`（`pyproject.toml:6`）だが、実体は `src/` 配下で import は全て
  `from src.xxx`（`import kagome` は成立しない）。
- `[tool.setuptools.packages.find] where = ["."]`（同:37-38）だと
  `src`/`scripts`/`tests` がトップレベル package 化する。
- `kagome.egg-info/` がリポジトリ直下と `src/` の2箇所に存在（editable install が複数箇所で走った痕跡）。
- `src/prep/openmm_equilibrate.py:86` が `from scripts._systems import _rdkit_mol` を行い、
  **ライブラリ層（src）が scripts 層に依存**する逆方向依存。

## Goal
配布名と import 名の不一致・逆方向依存を解消し、editable install を1経路に。

## Paper/guardrail anchor
N/A（実装/配布の健全性）。decision を1件残す。

## Steps（影響大のため段階的に）
1. **まず低リスク**: 不要な egg-info を1つに整理（`.gitignore` 済みか確認、リポジトリに
   コミットされているなら除外）。`pyproject.toml` の `packages.find` を `src` 配下に限定
   （例 `where=["."]` + `include=["src*"]`、または `src` レイアウトに移行）。
2. **逆方向依存の解消（RF4/RF9 と連携）**: `_rdkit_mol` など src から使う共有ロジックを
   `src/`（例 `src/chem/builders.py`）へ移し、`scripts/_systems.py` はそれを import する向きに反転。
   `src/prep/openmm_equilibrate.py` の `from scripts...` を `from src...` に。
3. **(任意・大)** import 名を `kagome` に統一する場合は `src/` → `pfpoly/` リネーム＋全 import 書換。
   影響が広い（`from src.` が全コード/テストに散在）ため、別タスク・別 PR を推奨。
   当面は「配布名と import 名が異なる」ことを README に明記する最小対応でも可。

## Acceptance criteria
- egg-info が1箇所、`pip install -e .` が1経路で再現可能。
- src → scripts の逆方向 import が無い。
- （リネームを選んだ場合）`import kagome` が成立し全テスト緑。

## Risk
- 中。import 書換は広範囲。最小対応（egg-info 整理＋逆依存解消）と、フルリネームを
  分離し、フルリネームは owner 判断で別 PR。

---

# RF12 — 幾何の orthorhombic 限定を明文化（epoxy/スラブ系の前提）

## Finding（事実）
`src/geometry.py` の `minimum_image`/`wrap_positions` は cell の対角成分のみ使用
（triclinic 非対応）。decisions.md に既記載だが、epoxy/CuO スラブ系（notes.md の系一覧）を
進める際の制約。

## Goal
制約を README/該当 docstring に明示し、非対角 cell が渡された場合に**早期に気づける**ようにする。

## Steps
1. `minimum_image`/`wrap_positions` に「orthorhombic 前提（対角のみ使用）」を docstring 強化
   （既に一部記載あり）。
2. 任意: 非対角成分が非ゼロの cell を検出したら `warnings.warn` か `ValueError`
   （静かに誤った最小像を返さないようにする）。挙動変更になるため decision を残す。

## Acceptance criteria
- 制約が README とコードで明示。非対角 cell が静かに誤計算されない。

## Risk
- 低。現行系は cubic/orthorhombic のみ。warn/raise を入れる場合は既存テストに
  非対角 cell が無いことを確認。

---

# RF13 — NVE+barostat 時の温度フォールバックの是正

## Finding（事実）
`src/workflows/polymerization.py:116-121` の `_integrator_temperature` は Langevin 以外で
常に 300 K を返す。Verlet + barostat 構成では実 KE と無関係に 300 K を仮定し
MC barostat の受理確率（`src/integrators/mc_barostat.py`）に使われる。

## Goal
barostat が使う温度を、サーモスタット温度が無い場合は**瞬時運動温度**から取るようにする。

## Paper/guardrail anchor
MC barostat 受理式は `kT` を要する（mc_barostat.py の docstring）。NPT は通常
Langevin 併用だが、Verlet+barostat の不整合を消す。

## Steps
1. Langevin のときは現行どおり target T。それ以外は
   `_instant_temperature(state.velocities, state.masses)`（RF9 で統一した関数）を渡す。
2. 渡す側（`_run_biased_phase`/`_run_unbiased_phase`/`_run_equilibration_phase` の barostat 呼び出し）で
   温度ソースを切替。

## Tests
- `tests/unit/test_mc_barostat.py` に、Verlet 構成で瞬時温度が使われる経路の smoke を追加。

## Acceptance criteria
- Verlet+barostat で 300 K 固定仮定が消え、瞬時温度が使われる。Langevin 経路は不変。

## Risk
- 低〜小。Langevin 中心運用では実害が無いが、整合性向上。瞬時温度は揺らぐため
  barostat の受理確率が変わる点を decision に1行記録。

---

## まとめ（実装者向けチェックリスト）

- [ ] RF6 ドキュメント追従＋supersede 明示（最初の足場、低リスク）
- [ ] RF3 不足テスト追加（RF4 の前提）
- [ ] RF9 重複集約 / RF7 ステップ記録 / RF13 温度フォールバック（小・独立）
- [ ] RF1 manifest に実効パラメータ記録（再現性）
- [ ] RF2 α 分母を n_monomers に一本化＋死蔵除去（PDF 確認済: α=1−[M]/[M]₀, 本文 p.9）
- [ ] RF5 スコアから nylon k–l を除外し d_ijkl=3項に（PDF 確認済: Eq.7 / Table S1/S2）
- [ ] RF4 MD ループ共通化＋群更新ストラテジ化（RF3 後、ビット一致厳守）
- [ ] RF10 pair_distances 整理（RF5 後）
- [ ] RF11 パッケージング整理（最小対応＋リネームは別 PR）
- [ ] RF12 幾何制約の明文化

各タスク完了で: paper artifact 参照を明記 → 必要なら decisions.md に decision 追記 →
テスト緑 → `.claude/hooks/pre-commit.sh` 通過 → 再現コマンドを記録。
