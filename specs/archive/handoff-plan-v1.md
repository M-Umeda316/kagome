# 作業計画書（外部委託向けハンドオフ）

最終更新: 2026-06-13
作成者: プロジェクトオーナー側（Claude Code 経由）
対象読者: 本リポジトリの開発を引き継ぐ外部委託先のエンジニア／研究者

---

## 0. この文書の目的と読み方

本書は、外部委託先が **本リポジトリの現状を正確に理解し、齟齬なく残作業を実行する**ための単一の参照点です。

- まず「1. プロジェクト前提」「2. 絶対遵守ルール」を読むこと。これは交渉不可です。
- 次に「3. 現状サマリ」で、どこまで完成しているかを把握すること。
- 実作業は「5. 作業計画（Phase 6 以降）」のタスク単位で進めること。各タスクには **目的 / 対象ファイル / 受け入れ基準 / 依存関係 / 論文アンカー** が定義されています。
- 判断に迷ったら「6. 質問必須トリガー（Ask-first）」を確認し、**勝手に進めず必ずオーナーに確認**すること。

⚠️ **本書の指示は `CLAUDE.md`（リポジトリルート）と矛盾してはなりません。**矛盾を見つけた場合は実装を止めてオーナーに報告すること。

---

## 1. プロジェクト前提

### 1.1 ミッション
以下の論文の手法を、**商用利用可能（commercial-safe）**な形で再現実装する。

- 論文: "Ready-to-Use Polymerization Simulations Combining Universal Machine Learning Interatomic Potential with Time-Dependent Bond Boosting for Polymer and Interface Design"
- arXiv: 2511.22874 / Mori, Tonogai, Miyazaki, Hayashi, Takayanagi (2025-11-28)

### 1.2 中核手法
**TDBB（Time-Dependent Bond Boosting / 時間依存ボンドブースティング）** を uMLIP（universal MLIP）と組み合わせ、系ごとのパラメータ調整なしに重合・架橋シミュレーションを行う。TDBB が**本プロジェクトで最も論文忠実度が要求される部分**である。

### 1.3 運用方針
> **paper-faithful where specified, backend-agnostic where licensing is unclear**
> （論文が明記する箇所は論文忠実に。ライセンスが不明確なバックエンドには依存しない。）

- TDBB のワークフロー・数式・スケジューリング・反応選択ロジックは論文通りに再現する。
- MLIP はクリーンなバックエンドインターフェース越しに差し替え可能に保つ。
- proprietary／ライセンス不明なバックエンドを**デフォルトにしてはならない**。

---

## 2. 絶対遵守ルール（Non-negotiables）

`CLAUDE.md` より抜粋。**すべてのコミットでこれを満たすこと。**

1. **すべての実装変更は、論文アーティファクト（claim / 数式 / 図 / 表 / 手法段落）を最低1つ参照すること。** コミットメッセージかコード中のコメントに `Eq. N` 等を明記する。
2. **論文で明示的に支持されない仮定は、実装前に `specs/decisions.md` に記録すること。** （後述のテンプレートを使用）
3. **図はスクリプトから再現可能であること。**手作業での図編集は禁止。
4. **すべての実験は seed / config パス / git SHA / バックエンド名 / 出力ディレクトリを記録すること。**（`src/workflows/manifest.py` の `RunManifest` を使う）
5. **TDBB 実装を特定のモデルプロバイダから独立に保つこと。**
6. **論文の大きな著作権テキストをリポジトリに貼らないこと。**短いノートと構造化サマリのみ。

### 2.1 商用ライセンス・ガードレール
- 新規依存を追加する前に**必ずライセンス確認**を行い、`specs/dependency-license-matrix.md` を更新する。
- 新規モデルバックエンドには `specs/dependency-license-matrix.md` に商用化ステータスを文書化する。
- ソフトが許諾的でも**モデル重みが制限される場合はバックエンド全体を制限扱い**とする。
- ライセンスが不明な場合は `blocked_pending_review` とマークする。
- 論文で使われている＝商用利用可、ではない。
- PFP / Matlantis 系バックエンドは、オーナーが利用権を確認するまで**デフォルトにしない**。

### 2.2 コーディング規約
- Python 3.11+ / public インターフェースには型ヒント必須。
- **インデントは4スペース、文字列はシングルクォート優先**（オーナー個人設定）。
- config は `dataclass` / typed dict で構造化。
- 数値カーネルとオーケストレーションコードを分離。
- 決定論的 seed を常に使用。
- numpy をコア唯一の実行時依存に保つ（ASE / OpenMM / MACE / orb はオプション依存）。
- I/O・シミュレーション状態・解析を明確に分離。

### 2.3 科学的実装ルール
- TDBB の数式を最も論文忠実度の高いコンポーネントとして扱う。
- 単位をコードと config コメントに明示する（後述「単位系」参照）。
- biased / unbiased セグメントの両方をログする。
- bond formation バイアスと bond dissociation バイアスを API レベルで区別する。
- reactive-pair 選択判断をすべて記録する（監査可能性）。
- 反応候補生成は MD を回さずにテスト可能に保つ。

---

## 3. 現状サマリ（2026-06-13 時点）

### 3.1 完了済み（Phase 1〜5）
| Phase | 内容 | 状態 |
|---|---|---|
| Phase 1 | TDBB フレームワーク（Eq. 2-8）、反応グループ、候補選択、biased/unbiased ループ | ✅ DONE |
| Phase 2 | トラジェクトリ JSONL 入出力、図再生成スクリプト | ✅ DONE |
| Phase 3 準備 | 原子質量、積分器プロトコル（VelocityVerlet / Langevin）、bond イベント追跡、転化率（Eq.11-12）、深さ分解反応密度（Eq.13） | ✅ DONE |
| Phase 3 | MACE-MP-0 バックエンド、ASE アダプタ、E2E 実行 | ✅ DONE（一部未完、後述） |
| Phase 4 | OrbMol-v2 バックエンド（OMol25+OPoly26 学習、商用可） | ✅ DONE |
| Phase 5 | 配置オーバーラップ検出、温度ログ、転化率プロット結線、各種バグ修正 | ✅ DONE |

- **テスト: 96 passed**（`python -m pytest tests/ -q`）。
- 実装済み数式: Eq. 2-8（TDBB コア）、Eq. 11-12（転化率）、Eq. 13（反応密度）。

### 3.2 実装ファイル地図（参照用）
```
src/boost/tdbb.py            # Eq. 2-5, 8: TDBB ポテンシャル・力・total_bias
src/reactive/groups.py       # Eq. 6: ReactiveGroup, ReactionTemplate
src/reactive/selection.py    # Eq. 7: find_candidates / score_candidates / select_non_overlapping
src/reactive/bonds.py        # BondTracker（attempt/confirm イベント）
src/workflows/polymerization.py  # biased/unbiased 交互ループ本体
src/workflows/manifest.py    # RunManifest（seed/SHA/backend 記録）
src/integrators/verlet.py    # VelocityVerlet（pre_force/post_force 分割）
src/integrators/langevin.py  # Langevin BAOAB サーモスタット
src/backends/base.py         # Calculator プロトコル
src/backends/toy.py          # テスト用トイバックエンド
src/backends/ase_adapter.py  # ASE Calculator アダプタ
src/backends/mace_backend.py # MACE-MP-0
src/backends/orb_backend.py  # OrbMol-v2
src/geometry.py              # minimum_image（直方体セル MIC）
src/units.py                 # FORCE_CONV, KB（LAMMPS 'real' 単位系）
src/analysis/conversion.py   # Eq. 11-12
src/analysis/density.py      # Eq. 13
scripts/run_smoke.py / run_mace.py / run_orb.py / reproduce_figures.py / _systems.py
configs/boost/paper_faithful.yaml  # 論文ハイパラ（★PDF 未照合、3.4 参照）
```

### 3.3 単位系（変更厳禁・全タスク共通の前提）
`specs/decisions.md`（2026-06-12「MD unit system」）で確定済み。

- 力: kcal/(mol·Å)、質量: amu、時間: fs、エネルギー: kcal/mol、距離: Å。
- `src/units.py` の `FORCE_CONV = 4.184e-4`。全積分器は `accel = FORCE_CONV · F/m` を使う。
- 運動エネルギー: `KE[kcal/mol] = 0.5·m·v² / FORCE_CONV`。
- タイムステップ 0.25 fs。
- **この単位系を壊す変更は科学的に致命的。** 触る場合は 6. のトリガーに従い必ず確認。

### 3.4 ⚠️ 重大な未解決ギャップ（残作業の核心）

外部委託先が最優先で認識すべき6点。

**G1: 実スケールの結合形成成功例が一度もない（最重要）**
- これまでの全実行で `confirmed_formations = 0`。つまり「ポリマーが1本もできていない」。
- 原因（`specs/decisions.md` 2026-06-12「Ethylene+ethylene...」に記録済み）:
  - (a) エチレン+エチレンの直接 C–C 結合形成障壁が ~40+ kcal/mol と高い（OrbMol-v2 が正しく再現）。
  - (b) 非周期系のため、バイアス除去後に分子が自由拡散して離れる。
  - (c) unbiased ステップが 50〜200 と少なく、結合状態をトラップできない（論文は 2000）。
- **TDBB 機械自体は検証済み**（バイアス印加・attempt 記録・outcome 判定は動作）。問題は系・条件・スケール。

**G2: 実行はすべてトイスケール**
- 既存 run（`runs/orb_v2`, `runs/orb_bond_demo` 等）は biased/unbiased 各 50〜100 ステップ、分子2個。
- 論文は **biased 2000 + unbiased 2000 ステップ、より大きな系、PBC あり**。
- `specs/tasks.md` の Phase 3 未完項目「Run paper-faithful scale」がこれに該当。

**G3: 周期境界条件（PBC）が Windows でブロック中**
- OrbMol-v2 の周期 PME には `nvalchemiops`（NVIDIA）が必要だが `blocked_pending_review`（ライセンス未確認）かつ Windows で torch.compile が走り失敗。
- 現状すべて非周期（cell=None）。**G1(b) の拡散問題は PBC がないと根本解決しない。**
- MACE-MP-0 は PBC 対応可能（要検証）。

**G4: ハイパーパラメータが論文 PDF と未照合**
- `configs/boost/paper_faithful.yaml` に `notes: Values should be rechecked against the paper PDF before final release.` と明記。
- `paper/notes.md` の既定値（λ=0.60, f2=10.0, γ=1.0, f1_max_form=250, f1_max_break=125 等）も PDF 一次照合が未完。
- γ の単位（kcal/(mol·step) か kcal/(mol·fs) か）も `specs/decisions.md` で「要確認」のまま。

**G5: 品質ゲート（フック／CI）スクリプトが実在しない（リポジトリ整合性の欠陥）**
- `.claude/hooks/pre-commit.sh` と `.claude/hooks/pre-run.sh`、および `README.md` の「Minimal commands」が以下4本を呼ぶが、**`scripts/` に存在しない**:
  - `scripts/validate_configs.py`
  - `scripts/check_dependency_licenses.py`
  - `scripts/check_seed_defined.py`
  - `scripts/check_output_path.py`
- 結果として **pre-commit / pre-run フックは現状すべて失敗する**（`set -euo pipefail` で即停止）。`CLAUDE.md` が掲げる「seed/SHA/backend 記録」「ライセンス確認」「config 検証」の自動ゲートが**機能していない**。
- これは委託先が最初の `git commit` でつまずく地雷であり、優先的に解消が必要（T6.3）。

**G6: 承認依存ファイルがライセンスマトリクスと矛盾（単一の真実源が二重化・陳腐化）**
- `specs/approved_dependencies.yaml`（ライセンスチェックフックが読む承認台帳）と `specs/dependency-license-matrix.md`（人間可読のマトリクス）が**食い違っている**:
  - approved_dependencies.yaml は `openmm` / `openmm-torch` / `pytorch` のみ列挙し、**現に稼働している `mace` / `ase` / `orb-models` / `OrbMol-v2` / `matplotlib` が未登録**。
  - `pytorch` が approved_dependencies.yaml では `review_required`、マトリクスでは `approved (BSD-3)` と**判定が矛盾**。
- 承認台帳が実態（Phase 3-4 で追加した backend 群）を反映しておらず、ライセンスガードレールが空洞化している（T6.2）。

### 3.5 その他の未完・要注意項目
- `specs/tasks.md` Phase 3: 「MACE-OFF or fine-tuned model for organic polymer accuracy」未着手（※ MACE-OFF23 は ASL ライセンスで商用不可、デフォルト化禁止）。
- Eq. 9, 10（`paper/claims.yaml` の `eval.qualitative` が参照）の実体が未実装。何の量かを PDF で確認のうえ要対応（G4 と統合）。
- 図 Fig. 2-6 を**論文の図と並べて定性比較**した記録がない。
- **転化率 α(t) の分母が不正確**: `scripts/reproduce_figures.py` は `n_total_sites` を `header['n_atoms']`（全原子数）で代用している。本来の分母は**反応サイト数**（エチレンなら反応 C のみ）であり、α のスケールが過大に出る。T9.1 で修正対象。
- **反応テンプレートが論文 Eq.7 の一般形と乖離**: 実装済みの machinery（`src/reactive/selection.py`）は N グループ・任意ペア集合 P を一般的に扱えるが、現在使われている `scripts/_systems.py::build_template_and_groups` は **2 グループ・1 ペアの簡略版**（C_donor–C_acceptor）。論文 Eq.7 は 4 原子・ペア集合 P={(i,j),(i,k),(j,l)} を例示する。この簡略化は `decisions.md` 未記録。T6.1 で論文の対象反応を確認し、簡略化の妥当性を記録すること。
- `.claude/agents/`（paper-analyst, license-auditor, reproduction-auditor, experiment-runner, implementation-planner）と `.claude/skills/`（read-paper, license-check, figure-reproduction 等）が用意済み。委託先が Claude Code を使う場合は活用してよい。

---

## 4. アーキテクチャ規約（コードを書く前に必読）

### 4.1 バックエンドインターフェース
すべての MLIP は `src/backends/base.py` の `Calculator` プロトコルに準拠すること。新規バックエンドは:
- `compute(positions, species, cell) -> (energy: float, forces: ndarray)` を実装。
- エネルギーは **kcal/mol**、力は **kcal/(mol·Å)** で返す（バックエンド境界で単位変換する。例: orb_backend の `FORCE_CONV` 周辺を参照）。
- `name` プロパティを持つ（manifest に記録される）。
- オプション依存（mace-torch, orb-models 等）は遅延 import にする。コア（numpy のみ）を壊さない。

### 4.2 ワークフロー
- `PolymerizationWorkflow.run()` が biased→unbiased を `n_cycles` 回繰り返す。
- biased フェーズ冒頭で候補選択 → `_build_pair_biases` → `BondTracker.record_attempts`。
- unbiased フェーズ末尾で `BondTracker.check_outcomes` → 確定結合は次サイクルで reactive group から除去。
- `save_interval > 0` かつ `output_dir` 指定時に `trajectory.jsonl` を出力。
- **積分器は差し替え可能**（デフォルト VelocityVerlet、Langevin はオプトイン）。

### 4.3 実験トレーサビリティ
新しい実行スクリプトは必ず:
- `RunManifest` を `output_dir/manifest.json` に保存（config パス / seed / backend / git SHA）。
- `summary.json` に主要結果（confirmed_formations 等）を出力。
- 図は `scripts/reproduce_figures.py` で trajectory.jsonl / bonds.jsonl から再生成（手編集禁止）。

### 4.4 decisions.md テンプレート（仮定を入れる際は必ず使用）
```
## YYYY-MM-DD: <決定の見出し>
- Context:
- Paper anchor:
- Decision:
- Alternatives considered:
- Scientific risk:
- Licensing/commercial impact:
- Follow-up:
```

---

## 5. 作業計画（Phase 6 以降）

> 各タスクは独立に着手できるよう設計しているが、**依存関係**欄を必ず確認すること。
> 着手前に「6. 質問必須トリガー」に該当しないか確認すること。

### Phase 6: 論文一次照合 + リポジトリ整合性回復（最優先・他の全タスクの前提）

#### T6.1 論文 PDF からハイパーパラメータ・数式を一次照合する
- **目的**: G4 を解消。実装済みの値・式が論文と一致することを確認し、`paper_faithful.yaml` の「未照合」注記を外す。
- **対象**: `paper/notes.md`, `paper/claims.yaml`, `configs/boost/paper_faithful.yaml`, `specs/decisions.md`
- **やること**:
  1. 論文 PDF（arXiv:2511.22874）の本文・表・補足から、λ, f2, γ, f1_max(form/break), timestep, biased/unbiased ステップ数、温度、系サイズを抽出して `paper/notes.md` の表と突き合わせる。
  2. γ の単位（step か fs か）を PDF で確定し、`specs/decisions.md` の該当 2 エントリ（2026-06-11）を更新。
  3. Eq. 9, 10 が何の量か（おそらく転化率／反応速度／密度関連）を特定し、`paper/claims.yaml` の `equations` に追記。
  4. **論文が対象とする反応の原子構成を確認し、現在の 2 グループ簡略テンプレート（`_systems.py`）との差異を評価する。** 論文 Eq.7 の P={(i,j),(i,k),(j,l)} に相当する 4 原子反応が対象なら、簡略化の妥当性（または 4 グループテンプレートの必要性）を `decisions.md` に記録（G6 関連項目）。
  5. 値に差異があれば config を修正し、差異内容を `decisions.md` に記録。
- **受け入れ基準**:
  - `configs/boost/paper_faithful.yaml` の `notes` から「rechecked が必要」の文言が消え、PDF 参照（節・表番号）に置き換わっている。
  - γ 単位の決定が `decisions.md` に確定形で記録されている。
  - Eq. 9, 10 が `claims.yaml` に定義されている。
  - 反応テンプレートの簡略化（2 グループ vs 論文の 4 原子）の妥当性が `decisions.md` に記録されている。
- **依存**: なし（最優先）。
- **論文アンカー**: 全 Eq./ハイパラ表、Eq. 7。
- ⚠️ **論文テキストを長文転載しない**（2. の禁止事項）。短い構造化ノートのみ。

#### T6.2 承認依存台帳を実態と同期する（G6 解消）
- **目的**: `specs/approved_dependencies.yaml` と `specs/dependency-license-matrix.md` の矛盾を解消し、ライセンスガードレールを実効化する。
- **対象**: `specs/approved_dependencies.yaml`, `specs/dependency-license-matrix.md`
- **やること**:
  1. マトリクスを真実源として、現に稼働中の `mace`(MIT) / `ase`(LGPL) / `orb-models`(Apache-2.0) / `OrbMol-v2`(Apache-2.0) / `matplotlib`(PSF) を approved_dependencies.yaml に追記。
  2. `pytorch` の判定矛盾を解消（マトリクス＝BSD-3 approved に合わせるか、台帳側の `required_evidence` を満たして確定する。どちらにするかは 2.1 の方針に従い、不確かなら `blocked` 維持）。
  3. 2ファイルが同一の結論を示すことを確認。
- **受け入れ基準**:
  - 稼働中バックエンドがすべて approved_dependencies.yaml に存在し、ステータスがマトリクスと一致する。
  - `pytorch` の判定が両ファイルで一致している。
- **依存**: なし。
- **論文アンカー**: N/A（リポジトリ整合性。`CLAUDE.md` 商用ガードレールに対応）。
- ⚠️ **6. のトリガー**: ライセンス判定を `approved` に変える根拠が不確かな場合はオーナーに確認（Ask-first）。

#### T6.3 欠落している品質ゲートスクリプトを復元する（G5 解消）
- **目的**: フックと README が参照するのに実在しない4スクリプトを実装し、自動ゲートを機能させる。
- **対象**: `scripts/validate_configs.py`, `scripts/check_dependency_licenses.py`, `scripts/check_seed_defined.py`, `scripts/check_output_path.py`（いずれも新規）、必要なら `.claude/hooks/*.sh` / `README.md` の整合確認
- **やること**:
  1. `validate_configs.py`: `configs/**/*.yaml` をロードし、必須キー・型・単位整合（3.3）を検証。
  2. `check_dependency_licenses.py --approved specs/approved_dependencies.yaml`: 台帳に `blocked_pending_review` や未承認が無いか検査し、あれば非ゼロ終了。
  3. `check_seed_defined.py` / `check_output_path.py`: 実行 config／引数に seed と出力パスが定義されているかを検査（pre-run フック用）。
  4. 4本がフックの呼び出し規約（引数・終了コード）と一致することを確認し、`pytest -q tests/unit` まで含めて pre-commit フックがグリーンになることを確認。
- **受け入れ基準**:
  - `bash .claude/hooks/pre-commit.sh` がエラーなく完走する。
  - 4スクリプトに最小限のユニットテストがある（不正 config／未承認依存で非ゼロ終了することを確認）。
  - README の「Minimal commands」が実際に動く。
- **依存**: T6.2（check_dependency_licenses は同期済み台帳を前提）。
- **論文アンカー**: N/A（`CLAUDE.md`「すべての実験は seed/config/SHA/backend を記録」「ライセンス確認」に対応）。

---

### Phase 7: PBC とライセンスのブロッカー解消（G3）

#### T7.1 PBC バックエンド経路の確立
- **目的**: G1(b) の分子拡散問題を根治するため、周期境界での MD を可能にする。
- **対象**: `src/backends/mace_backend.py`（PBC 検証）、必要なら `src/backends/orb_backend.py`、`specs/dependency-license-matrix.md`
- **やること（いずれか／オーナー確認のうえ選択）**:
  - 経路A: **MACE-MP-0 + PBC** が動くか検証（MACE はライセンスクリア＝MIT、商用安全）。直方体セルで `cell` を渡し、`src/geometry.py::minimum_image` と整合することを確認。← **推奨経路（ライセンス上安全）**
  - 経路B: OrbMol-v2 の `nvalchemiops` ライセンスをオーナーが確認し、解禁された場合のみ周期 PME を有効化。Windows の torch.compile 問題は別途切り分け。
- **受け入れ基準**:
  - 周期セルでエチレン箱の E2E が完走し、`runs/` に trajectory/bonds/summary が出る。
  - 使用バックエンドの商用ステータスが `dependency-license-matrix.md` に反映されている（B 経路なら nvalchemiops の判定更新）。
- **依存**: なし（T8 系と並行可）。
- **論文アンカー**: Section 2（PBC を全シミュレーションで使用）。
- ⚠️ **6. のトリガー**: nvalchemiops を有効化する前に**必ずオーナーに利用権を確認**（Ask-first: ライセンス不明バックエンド追加）。

---

### Phase 8: 結合形成の実証（G1・G2、本プロジェクトの科学的山場）

> 目標: **少なくとも1つの系で `confirmed_formations ≥ 1` を達成し、α(t) が単調増加するトラジェクトリを得る。**

#### T8.1 低障壁テスト系での結合形成のユニット実証
- **目的**: 化学的に難しいエチレンを避け、TDBB が「正しく動けば結合が確定する」ことを最小系で示す。
- **対象**: `scripts/_systems.py`（新規系ビルダー追加）、新規 `scripts/run_radical_demo.py` 等、`tests/`
- **やること**:
  1. ラジカル鎖成長など**低い結合形成障壁**を持つテスト系、または論文が例示する系の最小版を構築する（系選定は `decisions.md` に根拠記録）。
  2. paper-scale 寄りのパラメータ（unbiased を 1000〜2000 ステップに増やす）で実行。
  3. `confirmed_formations ≥ 1` を確認。
- **受け入れ基準**:
  - `summary.json` に `confirmed_formations ≥ 1`。
  - `bonds.jsonl` に `confirmed_formation` イベントが記録され、`reproduce_figures.py` の α(t) プロットが増加を示す。
  - 系選定理由が `decisions.md` に記録されている。
- **依存**: T6.1（正しいハイパラ）、できれば T7.1（PBC）。
- **論文アンカー**: Eq. 11-12（転化率）、Section 2（TDBB 実証系）。

#### T8.2 論文スケール実行（biased 2000 + unbiased 2000）
- **目的**: G2 を解消。論文と同じ時間スケールで定性トレンドを得る。
- **対象**: 新規 config `configs/eval/paper_scale.yaml`、実行スクリプト、`runs/paper_scale/`
- **やること**:
  1. T6.1 で確定したハイパラで、biased 2000 / unbiased 2000 ステップ × 複数サイクルを実行（PBC ありが望ましい＝T7.1 依存）。
  2. 計算コストが大きいため、まず系サイズ・サイクル数を見積もり、オーナーに実行規模（CPU/GPU・所要時間）を共有してから本実行。
  3. `RunManifest` で完全トレーサビリティを確保。
- **受け入れ基準**:
  - 論文スケールのパラメータで E2E 完走。
  - エネルギー（biased で上昇、unbiased で緩和）と α(t) のトレンドが論文と**定性的に**一致。
  - 結果と逸脱が `decisions.md` に記録（バックエンドが PFP でないことによる定量差は既知）。
- **依存**: T6.1, T7.1, できれば T8.1 で機械実証済みであること。
- **論文アンカー**: Fig. 1（スケジュール）、Section 2、Fig. 2-6。
- ⚠️ **6. のトリガー**: 計算規模が大きい本実行の前にオーナーへ実行計画を共有。

---

### Phase 9: 図の再現と定性比較（評価の締め）

#### T9.1 Fig. 2-6 の再現と論文との定性比較
- **目的**: 「定性トレンド一致」という Phase 3 受け入れ基準を、図レベルで証拠化する。
- **対象**: `scripts/reproduce_figures.py`（必要なら Eq. 9/10/13 のプロット追加）、`runs/<...>/figures/`、新規 `specs/figure-comparison.md`
- **やること**:
  1. T8.2 の出力から energy / α(t) / 温度 / 反応密度（Eq.13）の図を生成。
  2. 論文 Fig. 2-6 の各図が「何を示すか」を `paper/notes.md` に整理し、再現図と並べて**定性的な一致点・相違点**を `specs/figure-comparison.md` に記述（数値完全一致は不要、PFP 非使用のため）。
  3. Eq. 9, 10 が図に必要なら `src/analysis/` に追加実装（T6.1 で定義済みの式に基づく）。
  4. **転化率 α(t) の分母を修正**（3.5 参照）: `reproduce_figures.py` が `n_total_sites` に全原子数 `n_atoms` を使っている箇所を、実際の反応サイト数に置き換える（trajectory ヘッダかフラグで反応サイト数を受け渡す）。
- **受け入れ基準**:
  - 図はすべてスクリプトから再生成可能（手編集ゼロ）。
  - α(t) の分母が反応サイト数に基づいており、値域が妥当（0–100% を超えない）。
  - `specs/figure-comparison.md` に論文図 vs 再現図の対応と差異の説明がある。
- **依存**: T6.1, T8.2。
- **論文アンカー**: Fig. 2-6, Eq. 9-13。
- ⚠️ **6. のトリガー**: スムージング・フィルタ・平均化を論文と変える場合は確認必須。

---

### Phase 10: 任意拡張（オーナー承認後のみ）

以下は**オーナーが明示承認した場合のみ**着手。先行して着手しないこと。

- **T10.1 有機特化モデルの追加**: organic 精度向上のためのモデル。ただし **MACE-OFF23 は ASL ライセンスで商用不可、デフォルト化禁止**。商用可のモデル（または fine-tuning）に限る。`dependency-license-matrix.md` 更新必須。
- **T10.2 dissociation（解離）反応の有効化**: 現状 dissociation テンプレートは未使用。架橋・硬化シミュレーションで必要になった場合に、`decisions.md`（2026-06-12「Dissociation tracking...」）の閾値仮定を MLIP の PES と照合してから有効化。
- **T10.3 三斜晶セルの MIC 対応**: 現状 `minimum_image` は直方体のみ。非立方セルが必要になった場合に拡張。

---

## 6. 質問必須トリガー（Ask-first）

以下に該当したら**作業を止めてオーナー（m.umeda1502@gmail.com）に確認**すること。勝手に進めない。

1. 論文が**複数の妥当な数学的解釈**を許す場合。
2. **商用利用権が不明なバックエンド**を追加しようとする場合（例: nvalchemiops, PFP/Matlantis, MACE-OFF）。
3. 簡略化が **TDBB の科学的意味を変える**場合。
4. 図再現で **スムージング・フィルタ・平均化を論文と変える**場合。
5. **論文の引用なしにデフォルトのハイパーパラメータを導入**する場合。
6. **単位系（`src/units.py`）を変更**する場合（3.3 参照、科学的に致命的になりうる）。
7. 大規模計算（T8.2 等）を**本実行する前**（コスト・所要時間の事前共有）。

---

## 7. 各タスク共通の Definition of Done

`CLAUDE.md` より。タスクは以下をすべて満たして初めて完了:

1. コードが実装されている。
2. テストがパスする（`python -m pytest tests/ -q`、新規ロジックには新規テスト）。
3. 仮定が `specs/decisions.md` に文書化されている。
4. 再現コマンドが書かれている（README またはスクリプト docstring）。
5. 出力が config と seed にトレースできる（`RunManifest` 経由）。

### 7.1 推奨開発フロー（CLAUDE.md「Default workflow」準拠）
1. 論文を読み `paper/claims.yaml` を更新。
2. claim/数式を `specs/tasks.md` または `specs/decisions.md` のタスクに落とす。
3. 最小スライスを実装。
4. ユニット／統合テストを追加・更新。
5. 最小 E2E 再現を実行。
6. `runs/` に成果物を保存。
7. 図をスクリプトから再生成。
8. ドキュメント更新。

### 7.2 出力の優先順位（CLAUDE.md「Output expectations」）
1. spec / decision レコードの更新
2. 最小コード変更
3. テスト
4. 再現コマンド
5. 成果物の説明

---

## 8. 環境・実行メモ

- OS: Windows 11（PowerShell 主、Git Bash 併用可）。
- `python -m pytest tests/ -q` で全テスト（現状 96 passed）。
- バックエンド依存のインストール: `pip install -e .[mace]` / `.[orb]` / `.[plot]`。
- OrbMol-v2 は upstream が Windows 動作を保証していない。検証してから依存させること。
- 既存 run の場所: `runs/smoke`, `runs/mace`, `runs/orb_v2`, `runs/orb_bond_demo` 等（すべてトイスケール、formations=0）。

---

## 9. 優先順位の要約（迷ったらこの順）

1. **T6.3 + T6.2（壊れた品質ゲート復元 + 承認台帳同期）** ← 最初の commit でつまずく地雷。着手初日に解消。
2. **T6.1（論文一次照合）** ← すべての科学的前提。
3. **T7.1（PBC 経路確立、MACE 推奨）** ← G1 の根治に必要。
4. **T8.1（低障壁系で結合形成を1件実証）** ← 「ポリマーができる」初の証拠。
5. **T8.2（論文スケール実行）** ← 定性トレンド取得。
6. **T9.1（図再現・比較、α 分母修正含む）** ← 評価の締め。
7. T10.x は**オーナー承認後のみ**。

---

（本書に不明点・矛盾があれば、実装を止めてオーナーに確認すること。）
