# アルゴリズム整合性・正確性レビュー指摘事項 (2026-07-03)

対象: `src/kagome/` 全域 + `scripts/reconstruct_topology.py`, `scripts/check_topology_valence.py`
方法: コード精読 + 論文ノート (`paper/notes.md`, `paper/claims.yaml`)・`specs/decisions.md` との照合。
全指摘は実コードの該当行を確認済み。行番号は本日時点の `main` (11b355b) 基準。

**検証済みで問題なしと確認した領域**(参考): TDBB のポテンシャル/解析勾配/力の符号 (Eq.2-5, 8)、
Velocity Verlet、Langevin BAOAB(揺動散逸の厳密形)、MC バロスタットの (N+1)kT·lnV ヤコビアン、
FIRE、Maxwell-Boltzmann 初期化と 3N−3 DOF の整合、`units.py` の全定数(数値検証済み)、
直方体セルの minimum image / wrap、`total_bias` と `total_bias_fast` の代数的等価性、
チェックポイントのビット厳密 resume(物理状態)、選択パイプラインの決定性(単一 RNG + 状態保存)。

---

## 重大度: High — 反応判定の科学的意味に影響

### H1. バイアス中の `confirmed_formation` が非バイアス緩和後に再検証されない

- **箇所**: `src/kagome/reactive/bonds.py:57-97`(`check_reactions_during_bias`)、
  `src/kagome/reactive/bonds.py:129-130`(`check_outcomes` のスキップ)、
  `src/kagome/workflows/polymerization.py:870-879`(バイアス相の break)
- **事実**: バイアス相の最中に r ≤ 0.6·Σr_vdW に達したペアは、その場で
  `event_type='confirmed_formation'` として記録され `_reacted` に追加される。
  続く非バイアス相終端の `check_outcomes` は `_reacted` 入りのペアを
  「already confirmed during the biased phase」としてスキップするため、
  **f1 最大 250 kcal/mol の人工引力下で近接しただけのペアが、バイアス除去後に
  離れても確定反応として残る**。`bonds.py` のモジュール docstring
  (「confirms outcomes after unbiased relaxation」)とも矛盾する。
- **影響**: 変換率 α・DPn の過大計上。原子は `DefaultPostCycleUpdater._remove_pair`
  で反応グループから恒久的に除去され、`_apply_topology_updates` が実在しない結合を
  `topology.jsonl` に書き込む。不可逆に伝播する。
- **付随事実**: 候補窓の既定 r_min = 1.5 Å は C–C の閾値 r0 ≈ 2.04 Å(0.6·Σr_vdW)を
  下回るため、選択時点で既に r < 2.04 Å のペアはバイアス 1 ステップ目に
  動力学ゼロで「確定」しうる。
- **修正指針**:
  1. バイアス中の検知イベントは `tentative_formation` 等の**暫定イベント**とし、
     バイアス相の打ち切り(run-until-reaction)のトリガーにのみ使う。
  2. 確定 (`confirmed_*`) の発行は非バイアス緩和後の `check_outcomes` に一本化する。
     具体的には `check_reactions_during_bias` で `_reacted` に入れず、
     `_pending` に残したまま `check_outcomes` で距離を再判定する。
  3. グループ/トポロジー更新 (`updater.update`, `_apply_topology_updates`) は
     確定イベントのみを消費するので、上記変更だけで下流は正しくなる。
  4. 論文 §2.2 の解釈に関わるため、実装前に `specs/decisions.md` へ
     決定記録を追加すること(CLAUDE.md の Ask-first トリガー該当。
     「バイアス中の近接=即confirm」を意図的に採る場合も、その旨を記録する)。
- **テスト指針**: 「バイアス中に閾値を割った後、非バイアス緩和で離れるペア」を
  合成トラジェクトリで作り、confirmed にならないことを検証するユニットテストを
  `tests/unit/` に追加(現状 `BondTracker` のバイアス中/緩和後の相互作用を突く
  テストは存在しない)。

### H2. `confirmed_dissociation` が対の formation 不成立でも反応サイトを消費する

- **箇所**: `src/kagome/workflows/polymerization.py:250-276`
  (`DefaultPostCycleUpdater.update` / `_remove_pair`)
- **事実**: formation と dissociation の確定イベントは**独立に**消費され、
  `_remove_pair` は dissociation イベントの両原子を**全グループ**から除去する。
  同一候補(例: ナイロン縮合の N–C 形成 + N–H / C–OH 解離)で
  V^d 側だけが確定し V^f 側が不成立の場合、アミン N が結合を作らないまま
  `amine_N` グループから消え、反応サイトが不可逆に失われる。
- **影響**: 縮合系で変換率が停滞・過小になる。サイト簿記と実トポロジーが乖離。
- **修正指針**: 候補単位の**原子的(all-or-nothing)受理**にする。
  - `BondEvent` に候補 ID(またはサイクル内の candidate index)を持たせ、
    updater は「同一候補の formation が確定した場合のみ」その候補の
    dissociation 側のグループ編集を適用する。
  - あるいは updater を「候補ごとの結果 (formed / not formed)」を受け取る
    インターフェースに変更し、イベント列からの再構成をやめる。
  - H1 の修正(確定を check_outcomes に一本化)と同時に設計するのが効率的。

---

## 重大度: Medium

### M1. MC バロスタットの受理判定が TDBB バイアスエネルギーを無視

- **箇所**: `src/kagome/workflows/polymerization.py:654-677`(`_md_step` の呼び出し)、
  `src/kagome/integrators/mc_barostat.py:109-116`(`try_step` の dE)
- **事実**: バイアス相の NPT 体積試行で `try_step` に渡るのは `base_energy` のみ。
  dE は calculator の再計算だけで構成され、体積スケールでペア距離が変わることによる
  ΔV_bias は受理判定に入らない(バイアス力は受理**後**に再計算され、力は整合する。
  非整合なのは受理判定のみ)。
- **定量**: f2=10 Å⁻² の勾配最大点 (dr≈0.22 Å) で f1=250 のとき |dV/dr|≈670 kcal/(mol·Å)。
  最大体積移動 (Δln V=0.01) で r≈3 Å のペアは Δr≈0.01 Å → ΔV_bias は数 kcal/mol/ペア。
  kT(333 K)≈0.66 kcal/mol に対して無視できない。
- **修正指針**(いずれか):
  1. `try_step` がスケール後配置で `total_bias_fast` を再計算し、
     dE = ΔE_base + ΔE_bias で判定する(bias 計算用のコールバックか
     pairs/state/params を try_step に渡す)。
  2. バイアス相中はバロスタットを無効化する(`run_activation` は既に
     `enable_barostat=False` で同じ方針)。
  3. 「バイアスは人工的駆動ポテンシャルであり非バイアス面で受理する」を
     意図的仕様とするなら `specs/decisions.md` に決定記録を書く
     (decisions.md 該当行は「Active in both biased and unbiased phases」と
     述べるのみで、この除外を記録していない)。

### M2. resume 時、クラッシュ直前サイクルの artifact レコードが重複する

- **箇所**: `src/kagome/workflows/polymerization.py:481-503`、`src/kagome/io/trajectory.py:50`
- **事実**: `trajectory.jsonl` / `selection.jsonl` / `topology.jsonl` は resume 時
  `append=True` で開かれ、切り詰めは行われない。チェックポイントはサイクル境界で
  保存されるため、サイクル途中でクラッシュすると、途中まで書かれたレコードが残ったまま
  同サイクルが再実行され、同一 step/cycle のレコードが二重に追記される。
  (`bonds.jsonl` は復元済み tracker から全書き直しされるため安全。)
- **影響**: 物理状態はビット厳密だが、フレームを列挙する解析(密度の n_frames、
  α(t) プロット等)が再実行区間を二重計上する。
- **修正指針**: resume 時、各 JSONL をチェックポイントの step/cycle 以下のレコードまで
  切り詰めてから追記を開始する。実装は「行を読み step > ckpt_step の最初の行以降を捨てる」
  で十分(あるいはチェックポイントに各ファイルのバイトオフセットを保存する)。

### M3. `n_reactive_sites` のフォールバックが constraint-only グループも合算する

- **箇所**: `src/kagome/workflows/polymerization.py:418`
- **事実**: `n_monomers` 未指定時、`sum(len(g.atom_indices) for g in self.groups.values())`
  で**全グループ**(chain_C, vinyl_beta_C, radical_C を含む)を合算する。
  `analysis/conversion.py:16-29` は α の分母 = 初期モノマー数と明記し、
  constraint-only と radical_C の除外を要求している。ビニル系では分母が 3〜4 倍になる。
  この値は trajectory ヘッダと manifest に「α の分母」として記録される。
- **修正指針**: フォールバックを `monomer_site_count(self.groups)`
  (`kagome.analysis.conversion` に既存)に差し替える。モノマーグループ名が
  ビニル既定 (`vinyl_alpha_C`) と異なる化学系では引数で渡せるようにするか、
  `n_monomers` を必須引数に昇格して曖昧さを排除する。

### M4. `run_activation` が `config.seed` を再利用し RNG 列が本走行と相関する

- **箇所**: `src/kagome/workflows/polymerization.py:1106-1107`(`run()` 側は 505 行)
- **事実**: `rng=None` の場合 `np.random.default_rng(self.config.seed)` を新規生成。
  `run()` も同じ seed で生成するため、Langevin ノイズ列が activation と
  production 冒頭で完全に一致する(統計的相関アーティファクト)。
- **修正指針**: `np.random.default_rng([self.config.seed, 1])` のようにサブストリーム化
  するか、呼び出し側で生成した rng を必ず渡す規約にして `rng=None` 時は例外にする。
  再現性記録(seed → 派生規則)を manifest に残すこと。

### M5. `scripts/reconstruct_topology.py` の再構成がライブトポロジーと乖離しうる

- **箇所**: `scripts/reconstruct_topology.py:96-97, 111-112`
- **事実**(独立した 2 原因):
  1. azo C–N 結合を**無条件に全て**除去するが、`run_activation`
     (`polymerization.py:1084-1216`)は BondTracker を経由しない独自ループで
     解離を扱うため、解離イベントが `bonds.jsonl` に残らない。部分解離で終わった
     run では、実際には存在した結合が再構成から消える。
  2. 再生時に `apply_vinyl_addition` を無条件適用するが、ライブ側の
     `_apply_topology_updates`(`polymerization.py:1009-1020`)は過配位になる
     formation 編集を**スキップ**する。このスキップが発火した run では
     再構成結果が乖離し、過価数のトポロジーを出力しうる。
- **修正指針**: (1) activation の解離を `bonds.jsonl`(またはサイドカー)へ記録し
  再生時に読む。(2) `vinyl_addition_over_coordinates` による skip 判定を
  再生ループにも同一実装で移植する(共通関数化が望ましい)。
- **軽微**: 106 行の `first_step` は計算されるが未使用(dead code)。

### M6. `scripts/check_topology_valence.py` の近接スキャンが PBC を無視し最短非結合距離を取りこぼす

- **箇所**: `scripts/check_topology_valence.py:33-56`
- **事実**: `cKDTree(positions)` を `boxsize` なしで構築(周期境界越しの接触が不可視)。
  さらに大域最短距離 `gmin` を `tree.query(positions, k=2)`(最近接 1 個のみ)から
  導くため、真の最短非結合ペアの双方が「結合相手を最近接に持つ」場合に見逃す
  (全原子の最近接が結合相手なら gmin=inf のまま出力される)。
- **修正指針**: セルがある場合は座標が wrap 済みであることを利用して
  `cKDTree(positions, boxsize=box)` を使う。gmin は `query_pairs(r=十分大きい半径)`
  を走査して非結合ペアの最小値を取る(または非結合近傍が見つかるまで k を増やす)。

---

## 重大度: Low(潜在バグ・現行テンプレートでは未発火)

### L1. 質量キャッシュが `id(masses)` キーで stale になりうる

- **箇所**: `src/kagome/integrators/verlet.py:50-55`, `src/kagome/integrators/langevin.py:40-54`
- **事実**: 逆質量(Langevin では c2 も)を `id(masses)` でキャッシュ。
  (a) 同一配列の in-place 書き換えは検知されない。(b) CPython は解放済み
  アドレスを再利用するため、旧配列 GC 後に確保された**別の**配列が同じ id を
  持ち、前系の逆質量が適用されうる。縮合で原子が減り masses が作り直される
  ワークフローでは id 再利用の窓が現実に存在する。
- **修正指針**: キャッシュ対象の配列への**強参照**を保持し `masses is self._cached_masses`
  で比較する(参照保持により id 再利用も防げる)。

### L2. Langevin キャッシュが温度・摩擦係数をキーに含まない

- **箇所**: `src/kagome/integrators/langevin.py:40-54`
- **事実**: c1/c2 は kT と γ を含むが、キャッシュ無効化条件は `(dt, id(masses))` のみ。
  途中で `params.temperature_K` を変えるアニーリングを実装した場合、旧温度の
  サーモスタットが黙って使われ続ける(現行コードに温度変更箇所はない=潜在)。
- **修正指針**: キャッシュキーに `temperature_K` と `friction_per_fs` を追加するか、
  `LangevinParams` を frozen dataclass にする。

### L3. 積分器の `_get_box` が直方体チェックを迂回する

- **箇所**: `src/kagome/integrators/verlet.py:57-63`, `src/kagome/integrators/langevin.py:56-62`
- **事実**: `cell[0,0], cell[1,1], cell[2,2]` を直接読むため、三斜セルが渡ると
  `geometry.py` の規約(ValueError)に反して黙って直方体として wrap される。
- **修正指針**: キャッシュミス経路で `geometry.validated_box(cell)` を呼ぶ
  (キャッシュ済みなのでコストは償却される)。

### L4. 対称テンプレート・重複ラベルで候補列挙が壊れる(検証なし)

- **箇所**: `src/kagome/reactive/selection.py:50-56, 84-107`
- **事実**: `label_list.index(ps.group_a)` は同一ラベル 2 回出現時に同じ位置を返し、
  ペアキー `(i,i)` は照合ループ(`prev_depth < depth`)に決して一致せず距離窓が
  **黙って無効化**される。また `_enumerate_recursive` に同一原子重複ガードがなく、
  グループが原子を共有すると i==j 候補(後段 `add_bond` で ValueError)や
  (i,j)/(j,i) の二重列挙が起こる。現行のビニル/ナイロンテンプレートは
  ラベル・原子集合が互いに素なので未発火(潜在)。
- **修正指針**: `ReactionTemplate` 構築時にグループラベルの一意性と
  pair ラベルの相異を検証して例外にする。`_enumerate_recursive` で
  `atom_idx in chosen` をスキップ。対称反応を将来サポートするなら
  直積でなく組み合わせ (i<j) で列挙する。

### L5. `BondTracker._reacted` がペアキーのみ — 同一ペアは二度と反応できない

- **箇所**: `src/kagome/reactive/bonds.py:51, 74-75, 129-130`
- **事実**: `(min,max)` 原子ペアの集合で、イベント種別もサイクルも持たず、
  一度入ると消えない。形成→後の解離、失敗後の再形成が正当な化学系では
  2 回目のイベントが記録されない。ビニル系(各ペア一回限り)では無害。
- **修正指針**: 最低限 `(key, is_formation)` をキーにする。H1/H2 の修正で
  確定ロジックを再設計する際に、`_reacted` の意味(「このサイクルで処理済み」か
  「恒久的に反応済み」か)を明確化して docstring に書く。

### L6. `find_candidates` に結合済みペア除外がない

- **箇所**: `src/kagome/reactive/selection.py`(全体)
- **事実**: 距離窓とグループ所属のみで候補化する。r_min=1.5 Å は C–C 結合長
  1.54 Å より短く、グループ簿記が破れた場合に結合済みペアが候補になりうる。
  現状はグループ更新と `_valence_filter` が防波堤だが、後者は
  非ビニル系では no-op(`polymerization.py:911`)。
- **修正指針**: `BondTopology` を任意引数で受け取り `has_bond` のペアを除外する
  (1-3 近傍の除外も任意で)。

### L7. selection.jsonl の監査が valence filter 前に書かれる

- **箇所**: `src/kagome/workflows/polymerization.py:799-805`
- **事実**: `audited_selection` が `reason='selected'` を書いた**後**に
  `_valence_filter` が一部を落とすため、監査 artifact 上は「選択された」ままの
  候補が実際にはバイアスされない。RF18(選択理由の artifact からの再構成)に反する。
- **修正指針**: valence filter を監査書き込みの前に移すか、落とした候補に
  `reason='valence_dropped'` の追記レコードを出す。

### L8. 変換率フィットの時間軸に前処理ステップ分のオフセットが乗る

- **箇所**: `src/kagome/analysis/conversion.py:60-62, 103`
- **事実**: イベントの step は equilibration / activation を含むグローバル
  `state.step`。既定 `step_range` は 0 起点なので、α(t)=1−exp(−k t) のフィットに
  デッドタイムが混入し k_p が過小推定される。
- **修正指針**: production 開始 step を引数に取り、フィット前に
  `steps - production_start_step` で 0 起点化する(呼び出し側規約として
  docstring に明記するだけでも可だが、引数化が安全)。

### L9. `total_bias` / `total_bias_fast` の r<1e-12 ガードがエネルギーも捨てる

- **箇所**: `src/kagome/boost/tdbb.py:153-154, 213-214`
- **事実**: 原子が一致した場合に力(方向未定義なので正当)だけでなく
  エネルギー寄与も skip する。V^d(0)=f1(最大 125 kcal/mol)は定義済みなので、
  記録されるバイアスエネルギーが r→0 で不連続に 0 になる。実用上は到達不能。
- **修正指針**: `continue` 前に解析値を加算する
  (formation: `f1*(1-exp(-f2*r0**2))`、dissociation: `f1`)。

### L10. 縮合系のトポロジーが解離した脱離基結合を保持し続ける

- **箇所**: `src/kagome/workflows/polymerization.py:1025-1026`(コードコメントで既知)
- **事実**: `_apply_topology_updates` は非ビニル系では formation 結合の追加のみ行い、
  確定 dissociation に対応する N–H / C–OH 結合の除去をしない。
  縮合 run の `topology.jsonl` は stale な結合を含む。
- **修正指針**: confirmed_dissociation に対応する `remove_bond` を追加する。
  H2 の候補原子性(all-or-nothing)と併せて実装すること。

---

## 情報(バグではないが記録・改善推奨)

- **I1.** `MCBarostat.try_step` は棄却時 `(False, current_energy, None)` を返すが、
  docstring は「accepted configuration を反映」と述べる(`mc_barostat.py:134`)。
  唯一の呼び出し元はガード済みで実害なし。docstring 修正推奨。
- **I2.** FIRE は初回イテレーションで v=0 → power=0 → uphill 分岐に入り dt を
  半減させる(`minimize.py:94-104`)。ASE と同挙動で収束に影響なし。
  正準記述に合わせるなら初回のみ `power >= 0` 扱いにする。
- **I3.** `boost.advance()` がループ先頭で呼ばれるため、ループ内 1 ステップ目の
  f1=γ·1、ループ前評価(最初の half-kick 用)は f1=0。Eq.5 の離散化として
  正当だが、`specs/decisions.md` に一行記録しておくとよい。
- **I4.** checkpoint 復元時、checkpoint 側グループラベルが再構築側 `groups` に
  無い場合は黙ってスキップされる(`polymerization.py:435-437`)。ラベル集合の
  一致を assert すると呼び出し側の取り違えを早期検出できる。
- **I5.** `reactive/topology.py:107-116` の `_spare_hydrogen` は最初に見つかった
  H を除去する。現行モデル(ラジカル C にプレースホルダ H が 1 個)では一意だが、
  実 H を複数持つ 4 配位炭素では任意の H が黙って消える(docstring に前提記載あり)。
- **I6.** `analysis/density.py:42-44` は最上端ビン境界ちょうど (z==z_bins[-1]) の
  サンプルを捨てる(`np.histogram` は含める)。wrap 後の z∈[0,Lz) なら実害なし。
- **I7.** Langevin 下では COM 非保存のため DOF=3N−3 は近似(O(1/N)、
  docstring に記載済み・容認済み)。
- **I8.** `MCBarostat.stats` と TrajectoryWriter のフレームカウンタは
  checkpoint されない(診断値のみ、物理には影響なし)。

---

## 推奨着手順

1. **H1 + H2 + L5 + L10 をひとまとまり**として設計する
   (反応確定ライフサイクル: 暫定検知 → 緩和後確定 → 候補単位で原子的に
   グループ/トポロジー更新)。実装前に `specs/decisions.md` に
   「確定タイミングの論文解釈」の決定記録を書く(Ask-first トリガー該当)。
2. **M1**(バロスタット受理)は方針決定(dE に bias を含める / バイアス相で無効化 /
   意図的除外を記録)だけ先に行い、decisions.md に残す。
3. **M3, M4, L7, L9** は独立した小修正。既存テストへの追加とセットで着手可能。
4. **M2, M5, M6** は artifact/スクリプトの正確性の問題で、物理状態には影響しない。
   解析の信頼性のため H 群の次に対応。
5. **L1〜L4, L6** は潜在バグ。テンプレート検証 (L4) と参照保持キャッシュ (L1) は
   コスト小で先に潰せる。
