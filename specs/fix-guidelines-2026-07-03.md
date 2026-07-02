# 修正指針・検証結果 (2026-07-03)

ブランチ: `fix/review-findings-2026-07-03`
レビュー元: `specs/review-findings-2026-07-03.md`
検証方法: 4 並列エージェントによるコード読み取り検証 + `pytest` 全 367 テスト通過

---

## 実装済み修正一覧

### High (反応判定の科学的意味に影響)

| ID | 概要 | コミット | 検証 |
|----|------|---------|------|
| H1+L5 | バイアス中検知を暫定化、緩和後に確定 | a08938d | PASS (5/5 項目) |
| H2 | candidate_id による候補単位の原子的受理 | 8acb2e6 | PASS (5/5 項目) |

#### H1+L5 修正内容
- `check_reactions_during_bias`: `tentative_formation`/`tentative_dissociation` を発行（`confirmed_*` ではない）
- バイアス中は `_reacted` に追加しない。`_tentative` で重複防止のみ
- `record_attempts`: 新サイクル開始時に `_tentative` をクリア
- `check_outcomes`: 3-tuple キー `(*pair_key, is_formation)` で `_reacted` を管理（L5: 同一ペアの形成→解離を許容）
- 設計決定: `specs/decisions.md` D1 に記録済み

#### H2+L10 修正内容
- `PairBias`/`BondEvent` に `candidate_id: int = -1` フィールド追加
- `_build_pair_biases`: `enumerate(selected)` の index を `candidate_id` に割り当て
- `DefaultPostCycleUpdater.update`: dissociation は `confirmed_cids` に含まれる `candidate_id` の場合のみグループ編集を適用
- L10: `_apply_topology_updates` が `confirmed_dissociation` に対して `remove_bond` を実行
- `_topo_processed_dissociations` を初期化・checkpoint 保存・復元に対応
- `run_activation`: `confirmed_dissociation` イベントを `bond_tracker` に記録（`candidate_id=-1`: activation は独立した解離）
- 設計決定: `specs/decisions.md` H2 に記録済み

### Medium

| ID | 概要 | コミット | 検証 |
|----|------|---------|------|
| M1 | バロスタット受理に ΔV_bias を含める | 1cbcde9 | PASS |
| M2 | resume 時 JSONL truncation | b8a6856 | PASS |
| M3 | n_reactive_sites 分母を monomer_site_count に | c1c594d | PASS |
| M4 | activation RNG ストリーム分離 | c1c594d | PASS |
| M5 | reconstruct_topology 乖離修正 | 8bf6076 | PASS |
| M6 | check_topology_valence PBC 対応 | b8a6856 | PASS |

#### M1 修正内容
- `try_step` に `bias_energy_fn` コールバック引数追加
- `_md_step` がバイアス相で `total_bias_fast` ベースのクロージャを構築し渡す
- `dE = ΔE_base + ΔE_bias` で受理判定
- 設計決定: `specs/decisions.md` D2 に記録済み

#### M2 修正内容
- `_truncate_jsonl_after_step(path, max_step)` ユーティリティ関数追加
- resume 時に `trajectory.jsonl`, `selection.jsonl`, `topology.jsonl` を checkpoint step まで切り詰め

#### M3 修正内容
- `n_reactive_sites` フォールバック: `monomer_site_count(self.groups)` を優先
- 戻り値 0 の場合のみ全グループ合算にフォールバック + 警告ログ

#### M4 修正内容
- `run_activation` の `rng=None` 時: `np.random.default_rng([self.config.seed, 1])`（サブストリーム化）

#### M5 修正内容
- `run_activation` で解離イベントを `bond_tracker` に記録
- `reconstruct_topology.py`: `bonds.jsonl` の `confirmed_dissociation` 使用（レガシーは無条件 azo 除去にフォールバック）
- `vinyl_addition_over_coordinates` ガードを再生ループに追加
- dead code `first_step` 除去

#### M6 修正内容
- `_min_nonbonded_distance` に `box` パラメータ追加、`cKDTree(boxsize=box)` で PBC 対応
- CLI に `--box LX LY LZ` 引数追加

### Low

| ID | 概要 | コミット | 検証 |
|----|------|---------|------|
| L1 | Verlet/Langevin 質量キャッシュ: `is` 比較 | 67ef215 | PASS |
| L2 | Langevin キャッシュキーに温度・摩擦追加 | 67ef215 | PASS |
| L3 | `_get_box` で `validated_box()` 使用 | 67ef215 | PASS |
| L4 | `ReactionTemplate.__post_init__` バリデーション | 67ef215 | PASS |
| L5 | `_reacted` を 3-tuple 化 (H1 と同時) | a08938d | PASS |
| L6 | `find_candidates` に `topology` 引数追加 | 67ef215 | PASS |
| L7 | valence drop audit レコード出力 | c1c594d | PASS |
| L8 | `fit_conversion_exponential` に `production_start_step` | c1c594d | PASS |
| L9 | r≈0 ガードでエネルギー加算 | c1c594d | PASS |
| L10 | 縮合系 dissociation bond 除去 (H2 と同時) | 8acb2e6 | PASS |

### Info

| ID | 概要 | コミット | 検証 | 備考 |
|----|------|---------|------|------|
| I1 | バロスタット docstring 修正 | 180c319 | PASS | |
| I2 | FIRE 初回 uphill 分岐 | — | 未実装 | ASE 互換、影響なし |
| I3 | Eq.5 離散化解釈 | 5a8799d | PASS | decisions.md 記録のみ |
| I4 | checkpoint ラベル不一致警告 | 180c319 | PASS | |
| I5 | `_spare_hydrogen` 一意性 | — | 未実装 | 現行テンプレートで一意 |
| I6 | density.py ビン境界 | — | 未実装 | wrap 後 z∈[0,Lz) で無影響 |
| I7 | Langevin DOF=3N-3 | — | 未実装 | O(1/N)、docstring 記載済み |
| I8 | stats/フレームカウンタ非 checkpoint | — | 未実装 | 診断値のみ |

---

## 未実装の Info 項目について

I2, I5-I8 は以下の理由で実装対象外とした:
- いずれもバグではなく改善推奨
- 現行テンプレート・ワークフローで科学的影響なし
- レビュー findings にも「実害なし」と明記されている

縮合系（ナイロン）の本格運用時に I5 を再評価すること。

---

## チェックポイント互換性

- `CHECKPOINT_VERSION` を 1 → 2 に更新
- v1 checkpoint の `_reacted` (2-tuple) → v2 (3-tuple) のマイグレーションコード実装済み
- `topo_processed_dissociations` は `_ckpt.get(..., 0)` でデフォルト付き復元（後方互換）

---

## テスト追加

`tests/unit/test_bonds.py` に以下を追加:
- `test_in_bias_formation_tentative`: バイアス中は tentative イベントのみ
- `test_in_bias_dissociation_tentative`: 解離も同様
- `test_in_bias_no_duplicate_detection`: 同一ペアの重複検知防止
- `test_tentative_confirmed_after_unbiased_relaxation`: 暫定→確定フロー
- `test_tentative_not_confirmed_when_drifted_apart`: 暫定→離散→未確定

---

## 実装後の再検証で発見した問題 (2026-07-03 追記)

実装完了後、コードを再読して横断チェックを実施。以下 5 件の漏れ・新規問題を発見し、
コミット 2b7df7e で修正済み（全 374 テスト通過）。

### V1. M2 の漏れ: selection.jsonl の truncation が no-op だった【修正済み】
- **事実**: `_write_selection_audit` / `_write_valence_drop_audit` のレコードは
  `cycle` フィールドのみ持ち `step` を持たない。step ベースの
  `_truncate_jsonl_after_step` は selection.jsonl の全行を保持してしまい、
  M2 が意図した重複除去がこのファイルには効いていなかった。
- **修正**: `_truncate_jsonl_after(path, max_value, field=...)` に汎用化し、
  selection.jsonl は `field='cycle'`、閾値 `next_cycle - 1` で切り詰める。

### V2. checkpoint 互換バグ: 旧 BondEvent の candidate_id 欠落【修正済み】
- **事実**: pickle は `__init__` を経由せず `__dict__` を直接復元するため、
  candidate_id 追加**前**に保存された checkpoint の `tracker_events` を復元すると
  属性が欠落したままになる。run 終了時の `tracker.save()` 内 `asdict(ev)` が
  AttributeError で落ちる（resume した run が最後にクラッシュする）。
- **修正**: 復元時に `hasattr` チェックで `candidate_id = -1` をバックフィル。

### V3. L6 の未配線: topology ガードが呼び出し側で使われていなかった【修正済み】
- **事実**: `find_candidates(topology=...)` 引数は実装したが、
  `_run_biased_phase` の呼び出しで渡しておらず、結合済みペア除外が
  実際には一度も発火しない状態だった。
- **修正**: production 側の呼び出しに `topology=self._topology` を配線。
  activation 側は解離テンプレート（結合済みペアが対象）のため渡さない。

### V4. M1 docstring の矛盾: bias 二重計上を招く記述【修正済み】
- **事実**: `try_step` の docstring は「current_energy は bias 込みで渡すべき」と
  述べていたが、dE の内部計算は `Δbase + (new_bias − old_bias)` であり、
  bias 込みで渡すと old_bias が二重計上される。実際の呼び出し側は base のみ
  渡しており数学的に正しい（コードは正、docstring が誤り）。
- **修正**: 「BASE energy のみ。bias は bias_energy_fn が内部処理」と明記。

### V5. L4 の部分未対応: 対称ペア検証の欠落【修正済み】
- **事実**: `__post_init__` はラベル一意性とペア所属は検証していたが、
  レビュー L4 が指摘した `group_a == group_b`（対称ペア → 距離窓の黙殺）を
  検証していなかった。
- **修正**: 同一グループを参照するペアを ValueError で拒否。

### V6. H2 照合の堅牢化（潜在問題の予防）【修正済み】
- **事実**: candidate_id はサイクル毎に 0 から振り直される。updater の update()
  が複数サイクル分のイベントを一括処理した場合（現行フローでは起きないが）、
  別サイクルの同番号候補と誤照合しうる。
- **修正**: 照合キーを `(ev.cycle, ev.candidate_id)` に変更。

### 設計ノート（既知の挙動、修正不要と判断）
- `run_activation` の解離イベントは candidate_id=-1 で記録される。
  DefaultPostCycleUpdater では -1 は H2 ゲートをバイパスし無条件適用
  (RF15 のレガシー挙動)。ビニル系は VinylChainPropagationUpdater
  (formation のみ処理) を使うため影響なし。activation + DefaultPostCycleUpdater
  を組み合わせる将来の化学系では要注意（テスト
  `test_legacy_events_without_candidate_id_are_applied` が挙動を固定）。
- `_apply_topology_updates` の dissociation 除去は candidate_id ゲートなし。
  解離が check_outcomes で物理的に確定した以上、結合除去はトポロジーとして
  正しい（グループ簿記のゲートとは意味が異なる）。

### 追加テスト (tests/unit/test_workflow.py)
- `test_dissociation_skipped_when_formation_not_confirmed` — H2 ゲートの本体
- `test_dissociation_applied_when_same_candidate_formation_confirmed`
- `test_candidate_id_matching_is_per_cycle` — V6 の回帰テスト
- `test_legacy_events_without_candidate_id_are_applied` — RF15 互換
- `TestTruncateJsonl` 3 件 — step/cycle truncation と欠損ファイル no-op

---

## コミット一覧

```
2b7df7e fix: close gaps found in post-implementation verification
4eb21be docs(specs): add fix guidelines with verification results
11fb704 docs(specs): add decision record for H2 candidate_id design
8acb2e6 fix(reactive): atomic candidate acceptance and dissociation bond removal (H2/L10)
180c319 chore: barostat docstring fix and checkpoint label mismatch warning (I1/I4)
67ef215 fix: integrator cache safety, template validation, bonded pair exclusion (L1-L4/L6)
8bf6076 fix: reconstruct_topology divergence and activation event logging (M5)
b8a6856 fix: resume JSONL truncation and PBC-aware contact check (M2/M6)
c1c594d fix: Phase 3 independent small fixes (M3/M4/L7/L8/L9)
1cbcde9 fix(barostat): include bias energy change in NPT acceptance criterion (M1/D2)
a08938d fix(reactive): make in-bias detection tentative, confirm only after relaxation (H1+L5)
5a8799d docs(specs): add review findings and decision records D1/D2/I3
```
