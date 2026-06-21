# KAGOME

**Kinetic Accelerated Growth Orchestrated by Molecular Engine**

論文 "Ready-to-Use Polymerization Simulations Combining Universal Machine Learning Interatomic Potential with Time-Dependent Bond Boosting for Polymer and Interface Design" (arXiv:2511.22874, Mori et al.) の手法を、商用利用可能 (commercial-safe) な形で再現実装するリポジトリです。

## 特徴

- **TDBB (Time-Dependent Bond Boosting)** ワークフローの論文忠実な実装 (Eq. 2-8)
- **差し替え可能な MLIP バックエンド**: OrbMol-v2 (推奨), MACE-MP-0, ASE アダプタ, 古典 FF
- **決定論的再現性**: seed, config, git SHA, バックエンド名を RunManifest に記録
- **商用ライセンス・ガードレール**: 依存ライセンスの明示的管理
- **二段構成パイプライン**: 古典 FF (OpenMM/OpenFF) による構造準備 + ML ポテンシャルによる TDBB 本番

## クイックスタート

```bash
# インストール
conda create -n pfpoly-gpu python=3.12 -c conda-forge -y
conda activate pfpoly-gpu
pip install torch==2.10.0 --index-url https://download.pytorch.org/whl/cu128
pip install -e ".[dev,orb,plot]"

# スモークテスト
python scripts/run_smoke.py

# ビニルラジカル重合 (小規模)
python scripts/run_vinyl_aibn.py \
    --n-monomers 20 --n-initiators 1 --seed 42 \
    --backend orb --device cpu --no-barostat \
    --n-cycles 5 --biased-steps 2000 --unbiased-steps 500 \
    --output-dir runs/test_small
```

## プロジェクト構成

```
src/kagome/
  boost/          # TDBB バイアスポテンシャル (Eq. 2-5, 8)
  reactive/       # 反応グループ・候補選択 (Eq. 6-7)
  workflows/      # 重合ワークフロー (Fig. 1)
  backends/       # MLIP/Calculator アダプタ
  integrators/    # Velocity Verlet, Langevin, MC barostat, FIRE
  analysis/       # 転化率, 密度, Carothers 解析
  chem/           # RDKit 分子ビルダー
  io/             # JSONL トラジェクトリ I/O
  prep/           # 古典構造準備 (OpenMM/OpenFF)
docs/             # ドキュメント
specs/            # 要件定義・設計判断・ライセンス管理
paper/            # 論文からの構造化ノート
```

## ドキュメント

| ドキュメント | 内容 |
|---|---|
| [インストールガイド](docs/installation.md) | 環境構築 (Windows GPU / WSL2 / Linux) |
| [クイックスタート](docs/quickstart.md) | 段階的チュートリアル |
| [アーキテクチャ](docs/architecture.md) | モジュール構成・データフロー・単位系 |
| [バックエンド](docs/backends.md) | OrbMol-v2 / MACE / Classical / Toy の比較と選び方 |
| [コンフィグリファレンス](docs/configuration.md) | 全パラメータの詳細 |
| [論文再現ガイド](docs/paper-reproduction.md) | S1-S6 進捗・方程式対応表 |
| [ライセンスポリシー](docs/license-policy.md) | 商用利用ポリシー・依存ライセンス |

### 開発者向け

| ファイル | 内容 |
|---|---|
| `CLAUDE.md` | 開発ルール・非交渉要件 |
| `specs/decisions.md` | 設計判断と根拠の記録 |
| `specs/tasks.md` | タスクバックログ |
| `specs/handoff-plan-v5.md` | コードレビュー是正計画 (RF1-RF13) |
| `specs/handoff-plan-v6.md` | レビュー是正計画 (RF14-RF24) |

## 現在の進捗

| マイルストーン | 状態 |
|---|---|
| S1: 単鎖伝播デモ | 完了 |
| S2: メルト駆動での結合形成 | 完了 (OrbMol-v2, alpha=9.5%) |
| S3: 多ラジカル・高スピン近似 | 完了 (alpha=27.3%) |
| S4: AIBN 活性化 + 連鎖重合 | 完了 (C-N 解離 + 結合形成) |
| S5: 図の再現・比較 | 一部完了 |
| S6: 論文スケール実行 (200+10) | 未着手 (24 GB+ GPU 必要) |

## ライセンス

依存ライセンスの詳細は [ライセンスポリシー](docs/license-policy.md) を参照してください。
