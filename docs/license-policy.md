# ライセンスポリシー

KAGOME は商用利用可能 (commercial-safe) な実装を目指しています。新しい依存の追加には明示的なライセンス確認が必要です。

---

## 基本方針

1. **ソフトウェアとモデル重みの両方** のライセンスを確認する
2. SaaS としての使用は再配布とは異なる
3. 証拠が不完全な場合は `blocked_pending_review` としてブロック
4. CC-BY-4.0 は帰属表示が必要
5. 論文での使用は商用利用の許可を意味しない
6. ライセンスが不明確なバックエンドをデフォルトにしない

---

## 依存ライセンスマトリクス

### MLIP バックエンド

| コンポーネント | カテゴリ | ライセンス | 商用利用 | 備考 |
|---|---|---|---|---|
| orb-models (コード) | MLIP | Apache-2.0 | OK | |
| OrbMol-v2 (重み) | MLIP | Apache-2.0 | OK | |
| MACE (コード) | MLIP | MIT | OK | |
| MACE-MP-0 (重み) | MLIP | MIT | OK | |
| MACE-OFF23 (重み) | MLIP | ASL | **ブロック** | 商用利用不可 |
| nvalchemiops | MLIP 補助 | Apache-2.0 | OK | Windows torch.compile 制限あり |
| PFP/Matlantis | MLIP | 不明 | **ブロック** | blocked_pending_review |

### ML フレームワーク

| コンポーネント | ライセンス | 商用利用 |
|---|---|---|
| PyTorch | BSD-3-Clause | OK |
| ASE | LGPL-2.1 | OK (import-only) |

### 古典力場スタック

| コンポーネント | ライセンス | 商用利用 | 備考 |
|---|---|---|---|
| OpenMM (コア) | MIT | OK | |
| OpenMM (GPU) | LGPL | OK | |
| openff-toolkit | MIT | OK | |
| openff-interchange | MIT | OK | |
| openff-forcefields (Sage 2.x) | CC-BY-4.0 | OK | **帰属表示必要** |
| openff-nagl (コード) | MIT | OK | |
| openff-nagl (モデル) | CC-BY-4.0 | OK | **帰属表示必要** |

### ユーティリティ

| コンポーネント | ライセンス | 商用利用 |
|---|---|---|
| RDKit | BSD-3-Clause | OK |
| SciPy | BSD-3-Clause | OK |
| NumPy | BSD-3-Clause | OK |
| matplotlib | PSF (BSD-like) | OK |
| PyYAML | MIT | OK |

---

## ブロックされているバックエンド

### PFP/Matlantis

論文で使用されている uMLIP ですが、商用ライセンスの確認が取れていません。`blocked_pending_review` ステータスです。

### MACE-OFF23

Academic Software License (ASL) のため商用利用が制限されています。KAGOME では `create_mace_calculator` が MACE-OFF 指定を明示的に拒否します。

---

## 新しい依存を追加する際のチェックリスト

1. ソフトウェアのライセンスを確認
2. モデル重みのライセンスを確認 (ML バックエンドの場合)
3. `specs/dependency-license-matrix.md` に記載
4. `specs/approved_dependencies.yaml` に追加
5. ライセンスチェックスクリプトを実行:
   ```bash
   python scripts/check_dependency_licenses.py --approved specs/approved_dependencies.yaml
   ```
6. 帰属表示が必要な場合はその旨を記録

---

## 帰属表示

以下のコンポーネントは CC-BY-4.0 ライセンスのため帰属表示が必要です:

- **OpenFF Sage 2.x 力場**: Open Force Field Initiative
- **openff-nagl モデル**: Open Force Field Initiative

詳細なライセンス状況は `specs/dependency-license-matrix.md` を参照してください。
