# Handoff Plan v6 — リポジトリ全体レビュー（2026-06-20）の是正計画

Date: 2026-06-20
Paper: arXiv:2511.22874, Mori et al. (paper/2511.22874v1.pdf)
種別: **コード品質・科学的整合・再現性・ライセンス監査の是正**（v5 の続き）。
v5 を supersede しない。v5 = 2026-06-18 レビューの RF1–RF13（実装済）、
v6 = 2026-06-20 の多面的レビューで新規に洗い出した指摘（**RF14 以降に連番**）。

このドキュメントは 2026-06-20 のリポジトリ全体レビュー（12 観点 × 敵対的検証、
確定指摘 36 件・棄却 2 件・確認済みの強み 59 件）を、別担当者が着手できる粒度の
タスクへ落としたもの。各タスクは独立して着手可能（依存は明記）。元レビューの
全 evidence と検証 reasoning は当時の構造化結果に保存されている。

---

## まず: 「正しいので触らないこと」（今回のレビューで再導出・確認した中核）

以下は手計算/再導出で**正しいことを確認済み**。是正タスクではない。
リファクタや“ついで修正”で壊さないこと。

- **TDBB 方程式（Eq.2–5, 8）**: `src/boost/tdbb.py` は `paper/claims.yaml`/`notes.md` と代数一致。
- **TDBB の力の符号**: `total_bias` は真の物理力 `F_a = -dV/dr_a = +(dV/dr)·ê_ij` を返し、
  base 力・velocity Verlet の F/m 規約・Newton 第3法則と整合（`tdbb.py:138-164`）。
- **Langevin の揺動散逸**: BAOAB の `c2 = √(kT·FORCE_CONV·(1-c1²)/m)` は平衡分散と厳密一致
  （`src/integrators/langevin.py:50`）。単位変換まで正。
- **単位定数**: KB・FORCE_CONV・EV_TO_KCAL_MOL は正値（**唯一 ATM_TO_KCAL_MOL_A3 のみ誤り → RF19**）。
- **velocity Verlet / FIRE / Maxwell-Boltzmann 初期化**: いずれも正。
- **PBC**: `minimum_image`/`wrap_positions` は orthorhombic 限定で triclinic を `ValueError`（正しい設計）。
- **manifest 中核**: seed・config_path・backend・output_dir・**git SHA**・実効パラメータを記録
  （補強点のみ RF17）。
- **ガードレール設計**: toy が安全な既定、PFP/MACE-OFF23/nvalchemiops は `blocked_pending_review`、
  proprietary 既定なし、各バックエンドは optional-dep ガード付き（運用の穴のみ RF16/RF20）。

## 棄却された指摘（再提起しないこと）

敵対的検証で**誤検知と判定**された。蒸し返さない。

- 「OrbMol-v2 を既定バックエンドにしているのはライセンス違反」→ **誤り**。
  `specs/approved_dependencies.yaml:80-85` / `dependency-license-matrix.md` で
  **Apache-2.0・approved・推奨**。CLAUDE.md が禁ずる proprietary は PFP であり、PFP は
  正しく blocked かつ既定でない。
- 「サイクル毎の候補再列挙が反応済み原子を dedupe できない」→ **誤り**。updater による
  群除去で実際には整合している。

---

## 実装者向けグラウンドルール（CLAUDE.md 非交渉要件の再掲）

着手前に必ず守ること（v5 と同一。要点のみ）:

1. すべての変更は paper artifact（claim/式/図/表/method 段落）または明示の guardrail を1つ以上参照する。
   各タスクに "Paper/guardrail anchor" を明記済み。
2. paper で一意に定まらない解釈・仮定は、実装前に `specs/decisions.md` に追記する。
3. **TDBB のパラメータ・物理は変えない**（f2=10、f1_max=250/125、γ=1.0、λ=0.6、r0=λ·Σr_vdw、
   production の候補窓 [3,6]）。本計画の修正はいずれも「記録・整合・責務・テスト・ガードレール」であり、
   TDBB の科学的意味を変えない（例外: RF15 は解離側トポロジ更新の**未実装機能の追加**で、
   既存の vinyl 経路の数値は不変に保つ）。
4. **決定論**: seed を固定。リファクタは「ビット一致」を受入基準にする
   （`tests/unit/test_workflow.py::test_deterministic_with_seed` が pos1==pos2 を要求）。
5. 完了の定義: コード実装 + テスト緑 + 仮定の文書化 + 再現コマンド + 出力が config/seed に辿れる。
6. 各タスク完了時に `.claude/hooks/pre-commit.sh`（validate_configs → check_dependency_licenses → pytest）を通す。

---

## 優先度つきタスク一覧

| ID | 指摘（統合後） | 区分 | 元レビュー重大度 | 依存 | 工数 | リスク | 科学的意味 |
|----|------|------|----|------|------|--------|-----------|
| RF14 | TDBB の力の「大きさ」が未検証（有限差分/解析値テスト追加） | テスト | High | — | 小 | 低 | なし（回帰網） |
| RF15 | 解離イベントが反応群トポロジを更新しない＋nylon に鎖延長が無い | 科学的整合 | Medium | RF14 | 大 | 中 | **あり**（縮合/硬化のトポロジ進行） |
| RF16 | ライセンスゲートがブロックリスト方式で未登録依存を検知不能＋未登録依存の登録 | ライセンス | Medium | — | 中 | 中 | なし（ガードレール） |
| RF17 | プロベナンス補強（dirty tree / モデル重み identity / 分母 / demo manifest） | 再現性 | Medium | — | 中 | 低 | なし（記録のみ） |
| RF18 | 候補ランキング・非重複棄却の判断が監査用に保存されない | 監査性 | Medium | — | 中 | 低 | なし（記録のみ） |
| RF19 | 数値・単位の小修正（ATM 定数 / barostat N+1 / 温度 DOF / Carothers クランプ） | 数値整合 | Low | — | 小 | 低 | 弱くあり→要記録 |
| RF20 | backend のスピン/周期パスの誤解防止（set_spin no-op、nvalchemiops 事前チェック） | 整合性 | Low | — | 小 | 低 | なし（明確化） |
| RF21 | テスト網羅の補強（Verlet エネルギー保存、barostat 空振り assert、step-only テスト） | テスト | Medium/Low | — | 中 | 低 | なし |
| RF22 | 設定/スクリプトの整合・文書化（unbiased_steps、nylon 分母、ps1 half-scale、t-reset、γ 単位） | 整合性/文書 | Low | — | 小 | 低 | なし（記録のみ） |
| RF23 | 再現性 seed の徹底（RDKit conformer seed、prep の charge-method 既定乖離） | 再現性 | Low | — | 小 | 低 | なし |
| RF24 | 軽微な整理（PairBias 再エクスポート、重複再計算、docstring 過約束、空 __init__） | 整理 | Info | — | 小 | 低 | なし |

**ID と優先度**: RF14–RF18 が High/Medium、RF19–RF24 が Low/Info。番号は優先度グループ順だが
厳密な実行順ではない。実際の着手順は下記「推奨着手順」を正とする。

**推奨着手順**:
RF24 →（RF19 のうち ATM 定数のみ）→ RF14 → RF21 → RF16 → RF17 → RF18 →
RF15（テスト整備後）→ RF20 → RF22 → RF23 →（RF19 の barostat/温度 DOF は最後に decision つきで）。
理由: 低リスクの足場（RF24/ATM/テスト）を先に固め、科学挙動を変える RF15 と RF19(barostat/DOF) は
回帰網が揃ってから着手する。

---

# RF14 — TDBB の力の「大きさ」を有限差分/解析値で検証する（High）

## Finding（事実）
- 最も論文中核の `formation_force_magnitude`/`dissociation_force_magnitude`（`src/boost/tdbb.py:73-107`）
  と `total_bias`（同:121-166）について、`tests/unit/test_tdbb.py:79-121,147-175` は
  **力の符号と極値（≈0）のみ**を検査し、**数値の大きさを一切 assert していない**。
  - `TestFormationForce`(82-92): `assert dv[0] > 0.0` / `< 0.0` / `approx(0.0)` のみ。
  - `TestDissociationForce`(114-121): `assert np.all(dv < 0)` のみ。
  - `TestTotalBias`(150-175): `forces[0,0]` の符号と `forces[1] == -forces[0]`（Newton 第3法則）のみ。
- 結果: 解析微分の係数 2 の脱落（`dV^f/dr = f1·2·f2·(r-r0)·exp(...)`）や f2 のスケール誤りなど、
  **符号を保ったまま大きさを壊す回帰が全テストを通過**してしまう。
- ポテンシャル本体（`formation_potential`/`dissociation_potential`）は独立な期待値で検査済み（良）。
  欠けているのは力（=MD を駆動する一次量）の大きさ。

## Goal
力の大きさを (a) 解析式の独立計算値、(b) ポテンシャルの有限差分 `-(V(r+h)-V(r-h))/(2h)` の
両方に対して許容誤差つきで照合し、係数・スケール回帰を捕捉する。

## Paper/guardrail anchor
CLAUDE.md「Treat TDBB equations as the most paper-critical component」。Eq.2/Eq.3（力 = ポテンシャルの勾配）。

## Steps
1. `tests/unit/test_tdbb.py` に、複数 r（r<r0, r=r0, r>r0; dissociation は r 小/中/大）で
   `formation_force_magnitude`/`dissociation_force_magnitude` の戻り値を、
   別途手書きした解析式の値と `np.testing.assert_allclose(..., rtol=1e-10)` で照合。
2. 有限差分照合: 同じ r で `formation_potential`/`dissociation_potential` の中心差分（h=1e-6）と
   `*_force_magnitude` を `assert_allclose(..., rtol=1e-5)`。`*_force_magnitude` は dV/dr を返す
   実装である点に注意（力ベクトルへの符号適用は `total_bias` 側、ここでは dV/dr の値を検証）。
3. `total_bias` レベル: 既知の2原子配置で `forces[idx_a]` のノルムが `|dv_dr|` と一致し、
   方向が `e_ij` と平行（formation で r>r0 のとき a→b 向き）であることを数値で assert。
   3原子以上でも `Σ forces ≈ 0`（運動量保存）を assert。
4. PBC ありの `total_bias`（cell 指定）で minimum image 経由の r が使われることを、
   セル跨ぎ配置で1ケース検証。

## Tests / Acceptance criteria
- 上記が `pytest -q tests/unit/test_tdbb.py` で緑。
- 係数 2 を意図的に外す/ f2 を 2 倍する変異を入れると**少なくとも1テストが落ちる**こと（手元で確認）。

## Risk
- 低。テスト追加のみ。物理・決定論に影響なし。

---

# RF15 — 解離イベントのトポロジ更新と nylon の鎖延長を実装する（Medium・科学的整合）

## Finding（事実）
- `DefaultPostCycleUpdater.update`（`src/workflows/polymerization.py:137-152`）は
  `tracker.confirmed_formations()` のみを処理し、`atom_a`/`atom_b` を群から除去する。
  `_update_groups_after_cycle`（同:615-616）は updater にしか委譲しない。
- `tracker.confirmed_dissociations()` は `src/reactive/bonds.py:165` に存在するが、
  **群更新からは一切呼ばれない**（参照はレポート/集計のみ。リポジトリ全走査で確認）。
- 影響1（**解離トポロジ**）: nylon 縮合テンプレート（`amine_N–amine_H`, `carboxyl_C–carboxyl_OH` が
  `is_formation=False`）で結合が切れても、遊離した H/OH が群から除去・再分類されない。
  → 同じ H/OH がサイクル横断で再選択され、縮合トポロジが正しく進行しない。
- 影響2（**鎖延長**）: vinyl は β-C をラジカルに昇格（`VinylChainPropagationUpdater`、同:155-217）するが、
  nylon は `build_nylon66_system`（`scripts/_systems.py:510-538`）が updater を持たず
  `DefaultPostCycleUpdater` を使うため、アミド結合形成後に新生鎖末端が新たな反応末端にならない。
  連続段階重合が群レベルで表現されない。
- `specs/decisions.md`（2026-06-12「Reactive group atom removal」）は「消費のみがモデル化された効果」と
  記すが、解離側の更新を**未対応の既知課題**として明示していない（RF22 で文書化も）。

## Goal
（a）確定解離イベントで遊離原子を適切に群から除去/再分類する経路、
（b）nylon 段階重合で新生末端を再活性化する updater、を実装する。
vinyl 経路の数値は**不変（ビット一致）**に保つ。

## Paper/guardrail anchor
- 論文 §2.2（biased→unbiased→トポロジ更新のサイクル）。Table S2（nylon Condensation: i–j 形成、
  i–k/j–l 解離、k–l 水形成）。
- CLAUDE.md「Distinguish bond formation bias from bond dissociation bias」
  「Record all reactive-pair selection decisions」。
- **解離後にどの群をどう更新するか**は paper に逐一明示が無いため、**実装前に decisions.md に
  解釈 decision を残すこと**（ask-first トリガ: 「a simplification changes the scientific meaning of TDBB」に近い）。

## Steps
1. decisions.md に「縮合系の解離→群更新規則」decision を追加（どの原子を除去/再分類し、
   新生末端をどの群へ昇格するか。Table S2 の割り当て i=amine_N, j=carboxyl_C, k=amine_H, l=carboxyl_OH を引用）。
2. `PostCycleUpdater` プロトコル（RF4 で導入済、`polymerization.py:113-124`）に
   解離処理経路を追加するか、`update()` 内で `confirmed_dissociations()` も走査する。
   `DefaultPostCycleUpdater` を「形成原子を消費し、解離で遊離した原子を再分類/除去する」へ拡張。
3. `NylonCondensationUpdater`（新規、`VinylChainPropagationUpdater` と同階層）を実装し、
   アミド結合形成で生じる新生末端を次の反応に使える群へ反映。`build_nylon66_system` がこれを
   `PolymerizationWorkflow(updater=...)` に注入。
4. `processed_formations` と同様に `processed_dissociations` カウンタで二重処理を防ぐ。
5. vinyl 経路は updater 既定を変えない（既存 `VinylChainPropagationUpdater` のまま）。

## Tests
- `tests/unit/test_workflow.py`: nylon テンプレートで擬似 `confirmed_dissociation` を与え、
  遊離 H/OH が群から除去/再分類されること、同じ原子が次サイクルで再選択されないことを assert。
- 新 updater 経由で段階重合の末端再活性化が起きることを toy backend + 既知配置で検証。
- vinyl の `test_deterministic_with_seed` がビット一致のまま緑（vinyl 不変の保証）。

## Acceptance criteria
- `confirmed_dissociations()` が群更新で実際に消費される。
- nylon で解離原子の再選択が起きない。新生末端が次反応に使える。
- vinyl 数値が不変。decisions.md に解離→群更新規則の根拠（Table S2）がある。

## Risk
- 中。科学挙動（nylon の反応経路）が変わる。RF14/RF21 のテスト整備後に着手。
  既存 nylon 検証 run があれば1回再現確認し、変化を figure-comparison.md に記録。

---

# RF16 — ライセンスゲートを許可リスト方式にし、未登録依存を登録する（Medium・ライセンス）

## Finding（事実）
- `scripts/check_dependency_licenses.py:85-135` は **`blocked_pending_review` 登録済みかつ import 済み**の
  時だけ exit 1。**「全 import が approved 登録されているか」は検査しない**。
- 実際にゲートを走らせると `OK`/exit 0 だが、以下は `approved_dependencies.yaml` 未掲載のまま import されている:
  - `scipy`（`src/analysis/conversion.py:96`、`pyproject.toml:30-32` の [fit]、env では scipy==1.17.1。BSD 実体だが未記録）
  - `rdkit`（`src/chem/builders.py`・`scripts/_systems.py`。YAML は「single source of truth」を謳うが欠落）
  - `openff-units` / `pint`（`src/prep/openmm_equilibrate.py` 系。一部のみ登録）
- `_INSTALL_TO_IMPORT['mace-off23'] = []`（同:43）のため、`any(... if imp)` が常に False となり、
  blocked の MACE-OFF23（ASL、商用制限）を**構造的に検知できない**（モデル文字列で選択され import が無い）。

## Goal
（a）import される全トップレベル依存が registry に存在し `approved` であることを検査して未登録/未承認で失敗させる、
（b）未登録の実在依存（scipy/rdkit/openff-units/pint）を登録、（c）MACE-OFF23 を検知可能にする。

## Paper/guardrail anchor
CLAUDE.md 商用ガードレール「New dependencies require an explicit license check before adoption」
「If license status is unclear, mark it blocked_pending_review」「single source of truth」。

## Steps
1. `check_dependency_licenses.py` に**許可リスト検査**を追加:
   import スキャン結果のうち、標準ライブラリ・自前 package（`src`/`scripts`/`kagome`）・既知の安全集合を除いた
   トップレベル module が registry の `approved` に無ければ列挙して exit 1。
   既存のブロックリスト検査は temporso維持（多層防御）。誤検知を避けるため除外集合は明示的に定義。
2. `specs/approved_dependencies.yaml` と `specs/dependency-license-matrix.md` に
   `scipy`（BSD-3）・`rdkit`（BSD-3）・`openff-units`/`pint`（要確認: 確認できなければ `blocked_pending_review`）を追加。
   実際のライセンスは upstream LICENSE を確認して evidence に記載。
3. MACE-OFF23 検知: `mace_backend` がモデル文字列/パスに `off23`/`MACE-OFF` を含む場合に
   **実行時ガード**（`blocked_pending_review` の旨で `RuntimeError`）を入れるか、ゲートで
   コード中のモデル文字列を走査。`_INSTALL_TO_IMPORT` のコメントに「import では検知不能」の理由を残す。
4. `numpy→['numpy','np']`・`matplotlib→['matplotlib','plt']` のような import 正規表現で到達不能なエイリアスは
   コメントで「到達不能（regex は実 module 名を取る）」と注記（無害だが誤解防止）。

## Tests
- `tests/unit/test_quality_gates.py`: registry 未掲載の module を import する一時ファイルを作ると
  ゲートが exit 1 になること、approved のみなら exit 0 を assert。
- MACE-OFF23 ガードが発火することを assert（実 weights は不要、文字列分岐のみ）。

## Acceptance criteria
- import される全依存が registry+approved で覆われ、未登録で失敗する。
- scipy/rdkit/openff-units/pint が registry と matrix に存在。
- MACE-OFF23 が検知可能。`.claude/hooks/pre-commit.sh` 緑。

## Risk
- 中。許可リスト化は誤検知（標準ライブラリ・相対 import）を生みやすい。除外集合を堅牢に。
  openff-units/pint のライセンス確認が取れない場合は `blocked_pending_review` で安全側に倒す。

---

# RF17 — プロベナンスの補強（dirty tree / モデル重み identity / 分母 / demo manifest）（Medium・再現性）

## Finding（事実）
- `src/workflows/manifest.py:44-52` の `_get_git_sha` は `git rev-parse HEAD` のみ。
  未コミット変更があると（現に作業ツリーは dirty）、記録 SHA は実走コードと異なる。
- `src/backends/orb_backend.py:23-27,60-74`: manifest は `backend = calc.name`（例 `orb-orbmol_v2`）のみ記録し、
  **実際に使った重みファイル identity を記録しない**。local ckpt（`models/orbmol-v2-...ckpt`）→ DL フォールバックが
  無警告で起き、同一 manifest でも別重みになり得る。uMLIP 結果では重みは一次プロベナンス。
- 変換分母 `n_monomers` は trajectory ヘッダには入るが manifest には入らない（`polymerization.py:268-282`、partially-confirmed）。
- `scripts/demo_chain_propagation.py:152-154` は TDBB を走らせるが **manifest を出力しない**。

## Goal
manifest だけで「どのコード状態・どの重み・どの分母で走ったか」を復元可能にする。

## Paper/guardrail anchor
CLAUDE.md 非交渉要件「All experiments must record seed, config path, git SHA, backend name, and output directory」
「outputs can be traced to configs and seeds」。

## Steps
1. `_get_git_sha` を `git rev-parse HEAD` + `git status --porcelain`（空でなければ `-dirty` 付与、
   または `dirty: true` を別フィールドで記録）に拡張。`RunManifest` に `git_dirty: bool` を追加。
2. バックエンドが「解決済みモデル identity」を返せるようにする（例 `Calculator.model_id` プロパティ:
   orb は重みファイルのパス + 可能なら sha256、dtype/device）。manifest.extra に格納。
   orb の local→DL フォールバックは `logger.warning` を出す。
3. `polymerization.py` の manifest 構築で `n_monomers`（α 分母）を extra に格納。
4. `demo_chain_propagation.py`（および manifest 未出力の他 demo）で `RunManifest` を保存。

## Tests
- `tests/unit/test_manifest.py`: dirty 検出のモック（status 出力あり）で `git_dirty=True` を assert。
- manifest.extra に model identity と n_monomers が含まれることを assert。

## Acceptance criteria
- dirty 作業ツリーで走らせると manifest に dirty が記録される。
- manifest から重みファイル identity と α 分母が分かる。demo も manifest を残す。

## Risk
- 低。記録の追加のみ。subprocess 失敗時は従来どおり `'unknown'`/`False` にフォールバック。

---

# RF18 — 候補ランキングと非重複棄却の判断を監査用に保存する（Medium・監査性）

## Finding（事実）
- `src/workflows/polymerization.py:484-573`（`_run_biased_phase`）は
  `find_candidates` → `score_candidates`（d_ijkl 昇順）→ `select_non_overlapping` を実行するが、
  `CycleLog`（`polymerization.py:58-67`）に残るのは **`n_candidates`/`n_selected` の整数カウントのみ**。
- bias された pair は `record_attempts` が `bonds.jsonl` に書く（`bonds.py:97-117`、良）が、
  **列挙候補・各 d_ijkl スコア・overlap で落とされた候補は破棄**される。
  「なぜ候補 X を Y より落としたか」を成果物から再構成できない。

## Goal
選択判断（スコア付き候補一覧と棄却理由）を per-cycle で永続化し、選択を後追い監査可能にする。

## Paper/guardrail anchor
CLAUDE.md「Record all reactive-pair selection decisions for debugging and auditability」
「Make reaction candidate generation testable without running MD」。Eq.7（d_ijkl）。

## Steps
1. `select_non_overlapping`（`src/reactive/selection.py:122-134`）を、棄却理由（衝突した既使用原子）を
   返せる形に拡張するか、選択ログ用の補助関数を追加（純関数のまま、MD 非依存）。
2. `output_dir` 指定時に per-cycle の `selection.jsonl`（cycle, 各候補の atom_indices・score・
   selected/rejected・rejected 理由）を書き出す軽量 writer を追加（save_interval とは独立に、
   サイクル先頭で1回）。大規模系では候補数が多いため、上位 N または selected+棄却理由のみに
   絞るオプションを設け、**絞った場合は log() でその旨を残す**（沈黙の打ち切り禁止）。
3. docstring/decisions.md に「選択監査ログの粒度（全件 or 上位 N）」decision を記録。

## Tests
- `tests/unit/test_selection.py` または `test_workflow.py`: 既知配置で selection.jsonl に
  selected と rejected（理由つき）が期待どおり並ぶことを assert。

## Acceptance criteria
- run 後に selection ログから「候補・スコア・選択/棄却・理由」が復元できる。
- 既存の決定論テストがビット一致（ログ追加は数値に影響しない）。

## Risk
- 低。記録の追加。大規模系での I/O 量に注意（粒度オプションで制御）。

---

# RF19 — 数値・単位の小修正（ATM 定数 / barostat N+1 / 温度 DOF / Carothers クランプ）（Low）

複数の小さな数値ハイジーン。**ATM 定数のみ安全な確定修正**、他3つは挙動が微変するため decision を残す。

## Finding（事実）
1. **ATM_TO_KCAL_MOL_A3 = 1.4596e-5 は転記誤り**（`src/units.py:18-21`）。
   厳密値 = 101325 / ((4184/N_A)/1e-30) = **1.4584e-5**（約 0.08% 高い。コード自身のコメント計算
   101325/6.947e9 = 1.45854e-5 とも不一致＝末尾2桁の転記ミス）。MC barostat の圧力項
   （`src/integrators/mc_barostat.py:63,114`）に入るが、既定 1 atm では P·dV が dE に比べ極小で実害は小。
2. **MC barostat の Jacobian が N（正しくは N+1）**（`mc_barostat.py:100-114`）。提案を `ln(V)` 一様で行う場合、
   詳細釣り合いの正しい項は `(N+1)·kT·ln(V'/V)`。差は O(1/N) で実質無視可だが厳密でない。
3. **診断温度が 3N 自由度**（`src/integrators/init_velocities.py:54-62`）。COM 並進を引いた `3N-3` が厳密
   （大 N で無視可）。`maxwell_boltzmann_velocities` は COM 除去するため、測定 T が僅かに低めに出る。
4. **`dpn_from_bonds` の無警告クランプ**（`src/analysis/carothers.py:36-39`）。
   `p = min(n_bonds / (n_functional_groups / 2), 0.9999)` で DPn を 10000 に上限化するが、
   完全変換到達や p>1（分母誤り）を**沈黙で隠す**。

## Goal
ATM 定数を正値に直し、barostat/温度の厳密化と Carothers クランプの可視化を（decision つきで）行う。

## Paper/guardrail anchor
- ATM/温度/MC 受容式は MD の物理。NPT 受容式は `mc_barostat.py` の docstring に明記済。
- 「a default hyperparameter / 数値を変える」ため decisions.md 記録（特に barostat/DOF）。

## Steps
1. **(安全・即)** `ATM_TO_KCAL_MOL_A3 = 1.4584e-5` に修正し、コメントの導出も一致させる。
2. barostat: `delta_H` の Jacobian を `(n_atoms + 1) * kT * delta_ln_V` に変更し、
   「ln(V) 一様提案には N+1 が正しい」を decisions.md に記録（受容率が僅かに変わる）。
   ※ あるいは「N のまま（OpenMM 互換、差は O(1/N)）」を選ぶ場合もその旨を decision に明記。
3. 温度: `instant_temperature_K` の分母を `3*n - 3`（n>1）に変更し decision を記録、
   または「3N を診断用に維持」と明記。`_integrator_temperature`（RF13 と整合）と一貫させる。
4. Carothers: クランプ作動時に `logger.warning`、p>1 入力は警告（または ValueError）で表面化。

## Tests
- `tests/unit/test_mc_barostat.py` / `test_units.py`（無ければ新規）で ATM 定数の値を assert。
- barostat 受容式の Jacobian 係数をユニットで検査（N+1 を採った場合）。
- Carothers: クランプ境界（p≥1）で警告が出ることを assert。

## Acceptance criteria
- ATM 定数が 1.4584e-5。barostat/温度 DOF の選択が decisions.md にある。
- Carothers の上限到達/異常入力が沈黙しない。

## Risk
- 低。barostat/温度の変更は受容率・診断 T を僅かに動かすため、既存 NPT/温度テストの許容幅を確認。

---

# RF20 — backend のスピン/周期パスの誤解を防ぐ（Low）

## Finding（事実）
- `ASECalculatorAdapter`（`src/backends/ase_adapter.py:19-39`）は base の `set_spin` を override しないため、
  `--backend mace` では `calc.set_spin(...)`（`scripts/run_vinyl_aibn.py:385`）が **base の no-op**
  （`src/backends/base.py:32-37`）に当たるのに、同:386 が「Spin switched: 1 -> %d」とログする。
  MACE-MP-0 はスピン非依存なので科学は妥当だが、**ログが誤解を招き監査痕として不正確**。
- `src/backends/orb_backend.py:109-120`: 周期 OrbMol-v2 は未確認 nvalchemiops（PME）に依存し
  実行時失敗しうるが、アダプタに**明示の事前チェックが無い**。

## Goal
スピン無視を可視化し、周期パスの前提を明確に失敗/警告させる。

## Paper/guardrail anchor
CLAUDE.md「auditability」「no proprietary/license-unclear backend as default」
（nvalchemiops は `blocked_pending_review`）。

## Steps
1. `Calculator.set_spin` を呼んだ際、無視するバックエンドでは `logger.warning`
   （「this backend ignores spin; requested=…」）を出す。または run script 側でバックエンドが
   spin を尊重するか確認してからログを出す。
2. orb 周期パス（cell != None）で nvalchemiops 不在/未承認を**早期に明確なエラー**にする
   （`blocked_pending_review` の旨を含む）。`run_orb` が cell を渡す経路に注意。

## Tests
- `tests/unit/test_backends.py`: ASE 経路で set_spin が warning を出す（または no-op が明示される）ことを assert。

## Acceptance criteria
- MACE 経路で「Spin switched」だけが残り実際は無視、という不整合が解消。
- 周期 OrbMol-v2 が未承認依存で静かに失敗しない。

## Risk
- 低。ログ/早期失敗の追加。既定（非周期 orb / mace）に影響なし。

---

# RF21 — テスト網羅の補強（Verlet エネルギー保存・barostat 空振り・step-only テスト）（Medium/Low）

## Finding（事実）
- `tests/unit/test_integrators.py:73-157`（`TestVelocityVerlet`）は自由粒子・定数力の運動学を
  独立期待値で検査（良）するが、**保存力下での全エネルギーの非発散/対称性（symplectic 保存）を検査するテストが無い**。
  pre_force で新しい力を使う・half-kick を二度適用する等の回帰が、定数力テストを通過し得る。
- `tests/unit/test_mc_barostat.py:87-99,138-154`: cell 変化/位置スケールの assert が `if accepted:` の
  内側にあり、**受理が一度も起きないと空振り（vacuously pass）**になり得る。
- `tests/unit/test_workflow.py:204-240`（`test_with_masses`/`test_with_langevin`）は
  `assert state.step == 10` のみで、**masses 経路 vs massless 経路、Langevin の温度効果を比較しない**
  （配線確認テスト。中核 kernel は別途検査済なので重大度は低いが網羅穴）。

## Goal
中核 kernel の「壊れたら落ちる」テストを足す。

## Paper/guardrail anchor
CLAUDE.md「Treat TDBB equations as the most paper-critical」「deterministic seeds」。
Velocity Verlet の定義的性質（symplectic エネルギー保存）。

## Steps
1. ToyCalculator（LJ）または調和ばねの保存力下で velocity Verlet を多数ステップ回し、
   全エネルギー（KE+PE）のドリフト/振幅が有界であることを assert。
2. `test_mc_barostat.py`: 受理が確実に起きる設定（大きな max_volume_change_frac や低 PES 勾配）で
   cell/位置スケールを検査し、**少なくとも1回は accepted=True を assert**（空振り防止）。
   受理率テストは別途統計で。
3. `test_workflow.py`: masses あり/なしで結果が変わること、Langevin で温度がターゲット近傍に
   収束する統計（既存 `TestIntegratorTemperature` の流儀）を追加。

## Acceptance criteria
- Verlet エネルギー保存テストが存在し、half-kick 二重適用などの変異で落ちる。
- barostat の構造 assert が空振りしない。
- workflow テストが masses/integrator の効果を実際に検査。

## Risk
- 低。テスト追加のみ。

---

# RF22 — 設定/スクリプトの整合と文書化（Low）

## Finding（事実）
- `scripts/run_vinyl_aibn.py:70-71`: 既定 `--biased-steps`/`--unbiased-steps` = 500/500 で、
  論文値（2000/2000、`configs/boost/paper_faithful.yaml:9`）と乖離。S6 本番は 2000/500。
  decisions.md（S2/S4）で根拠付き・実効値は manifest 記録だが、YAML と既定の差が読者を混乱させ得る。
- `scripts/run_nylon66.py:149-162,198-203`: `wf.run` に `n_monomers` を渡さないため trajectory ヘッダの
  `n_reactive_sites` は4群和（`polymerization.py:282`）になるが、図コマンドは `--n-reactive-sites` を
  (amine_N + carboxyl_C) で上書き → **ヘッダと図の分母が不一致**（nylon は本来 Carothers p で別管理なので
  影響は限定的だが drift）。
- `scripts/run_s6_paper_scale.ps1:11-16,85-106`: half-scale ガイドコメントが実走に含む
  `--activation-steps 5000 / --minimize / --minimize-fmax 1.0` を一部欠落（コピペで activation-steps が
  既定 3000 になる）。.sh には half-scale 例が無く2スクリプトのガイドが非対称。
- `src/workflows/polymerization.py:504,520-526`: バイアス時刻 t を**セグメント毎にリセット**する解釈
  （decisions.md 記載ありだが論文に明文なし）。
- `src/boost/tdbb.py:16-26`: `gamma` に単位注記が無い（unit-sensitive。config 側はコメントあり）。

## Goal
設定の真実源乖離・分母 drift・ガイド非対称・解釈の明文化を解消（いずれも記録/文書/軽微整合）。

## Paper/guardrail anchor
CLAUDE.md「figures reproducible from scripts」「outputs traced to configs」、
「a default hyperparameter without paper citation」回避。論文 p.7（2000 ステップ交替）/ p.9 Fig.2（α 分母）。

## Steps
1. run script の既定 steps をどうするか decision（論文値に合わせる or 既定を維持して manifest 記録に委ねる）を記録。
   少なくとも README/`--help` に「既定は探索用、本番は S6 スクリプト/明示フラグ」と明記。
2. `run_nylon66.py`: `wf.run(..., n_monomers=...)` を渡してヘッダと図コマンドの分母を一致させる
   （vinyl の `run_vinyl_aibn.py:411` に倣う）。または nylon は Carothers p のみを図化し α を出さない方針を明記。
3. `run_s6_paper_scale.ps1` の half-scale コメントに不足フラグを補完し、.sh にも同等の half-scale 例を追記（対称化）。
4. decisions.md/`tdbb.py` docstring に「t のセグメント毎リセット」「gamma の単位（kcal/(mol·step)）」を明記。

## Acceptance criteria
- steps 既定の扱いが文書化。nylon のヘッダ/図分母が一致。ps1/sh の half-scale ガイドが対称。
- t-reset と gamma 単位が docstring/decisions に明記。

## Risk
- 低。記録/文書と軽微なスクリプト整合。物理不変。

---

# RF23 — 再現性 seed の徹底（RDKit conformer / prep charge-method 既定）（Low）

## Finding（事実）
- `scripts/_systems.py:286-308,462-474,602-624` 系で **RDKit コンフォマー生成 seed が 42 にハードコード**され、
  `--seed` から繋がっていない。同一 `--seed` でも RDKit 由来の初期座標は seed に追従しない。
- `scripts/prep_structure.py:92-95`: 既定 charge-method が config dataclass の既定と乖離（prep スクリプトと
  ライブラリ層で既定が食い違う）。

## Goal
初期構造生成まで含めて `--seed` 一本で決定論を貫く。prep の既定を一致させる。

## Paper/guardrail anchor
CLAUDE.md「Use deterministic seeds everywhere possible」「record seed」。

## Steps
1. RDKit conformer 生成に渡す `randomSeed` を `--seed`（または派生 seed）から供給。
   ハードコード 42 を撤去。`_rdkit_mol`/`_rdkit_3d` の seed 引数を呼び出し元の seed に接続。
2. `prep_structure.py` の charge-method 既定を config dataclass の既定と一致させる（どちらを正とするか decision）。

## Tests
- `tests/unit/test_systems.py`/`test_prep.py`: 同一 seed で生成座標が再現、異 seed で変わることを assert。

## Acceptance criteria
- 初期構造が `--seed` に決定論的に追従。prep の既定が一意。

## Risk
- 低。seed 経路の接続。既存 run の座標が変わり得るため、固定 seed の回帰は新基準で取り直し。

---

# RF24 — 軽微な整理（Info）

## Finding（事実）
- `src/boost/__init__.py:1-11`: `total_bias` が依存する `PairBias` が boost パッケージから**再エクスポートされていない**
  （利用側は深い import が必要）。
- `src/workflows/polymerization.py:534-538,555-558`: `min_form_dist` 診断が biased 各ステップで
  pair 距離を再計算し、BondTracker 側の再計算と**二重**。
- `src/chem/builders.py:7-31`: `box_from_density` が分子数を **SMILES 文字列でキー化**し、同一 SMILES を黙ってマージ。
- `src/analysis/carothers.py:29-39`: docstring が「stoichiometric-imbalance handling / Xn 関係」を約束するが未実装（過約束）。
- `src/backends/__init__.py:1-1`: 空でパッケージ API 面を提供しない。

## Goal
読みやすさ・整合のための無害な整理。

## Paper/guardrail anchor
N/A（コード健全性）。物理不変。

## Steps
1. `boost/__init__.py` に `PairBias` を re-export（`__all__` も更新）。
2. `min_form_dist` を BondTracker の既存計算結果から得るか、診断の重複計算を1本化（RF4 の `_md_step` と整合）。
3. `box_from_density` の SMILES キー化挙動を docstring に明記、または意図しないマージを検出/警告。
4. `carothers.py` docstring を実装に合わせて縮約（未実装機能を約束しない）。
5. `backends/__init__.py` に主要シンボル（`Calculator`, 各 factory）を re-export（任意）。

## Acceptance criteria
- 上記が解消され、既存テストが緑（数値不変）。

## Risk
- 低。整理のみ。

---

## まとめ（実装者向けチェックリスト）

- [ ] RF24 軽微整理（足場・無害）
- [ ] RF19(ATM 定数のみ) 安全な確定修正
- [ ] RF14 TDBB 力の大きさテスト（High・回帰網）
- [ ] RF21 テスト網羅の補強（Verlet 保存ほか）
- [ ] RF16 ライセンスゲートを許可リスト化＋未登録依存の登録
- [ ] RF17 プロベナンス補強（dirty/重み identity/分母/demo manifest）
- [ ] RF18 選択判断の監査ログ
- [ ] RF15 解離トポロジ更新＋nylon 鎖延長（科学挙動変更、テスト整備後）
- [ ] RF20 backend スピン/周期パスの明確化
- [ ] RF22 設定/スクリプト整合・文書化
- [ ] RF23 再現性 seed の徹底
- [ ] RF19(barostat/温度 DOF) 厳密化（decision つき・最後）

各タスク完了で: paper/guardrail anchor を明記 → 必要なら decisions.md に decision 追記 →
テスト緑 → `.claude/hooks/pre-commit.sh` 通過 → 再現コマンドを記録 → 出力が config/seed に辿れることを確認。

## 元レビューとの対応（重大度マップ）

| 重大度 | 元レビュー指摘 | 対応タスク |
|---|---|---|
| High | TDBB 力の大きさ未検証 | RF14 |
| Medium | 解離→トポロジ更新なし / nylon 鎖延長なし | RF15 |
| Medium | ライセンスゲート blocklist-only / 未登録依存 | RF16 |
| Medium | git dirty 無視（+重み identity/分母/demo） | RF17 |
| Medium | 候補選択の監査ログ欠落 | RF18 |
| Medium | Verlet エネルギー保存テスト欠落 | RF21 |
| Low | ATM 定数 / barostat N+1 / 温度 DOF / Carothers クランプ | RF19 |
| Low | set_spin no-op / nvalchemiops 事前チェック | RF20 |
| Low | barostat 空振り assert / step-only テスト | RF21 |
| Low | unbiased_steps / nylon 分母 / ps1 half-scale / t-reset / gamma 単位 | RF22 |
| Low | RDKit conformer seed=42 / prep charge-method 既定 | RF23 |
| Info | PairBias 再エクスポート / min_form_dist 重複 / docstring 過約束 / 空 __init__ / SMILES マージ | RF24 |

（棄却: OrbMol-v2 既定=Apache-2.0 で正当、候補 dedup は整合済 — 本計画の対象外）
