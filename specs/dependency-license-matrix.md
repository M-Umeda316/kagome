# Dependency and commercial-use matrix

| Component | Category | License / status | Commercial-use default? | Notes |
|---|---|---|---|---|
| MACE (code) | MLIP framework | MIT | **yes** | default uMLIP backend |
| MACE-MP-0 model | pretrained weights | MIT | **yes** | 89-element universal potential, default model |
| MACE-OFF23 model | pretrained weights | ASL (restricted) | no | organic chemistry focus, ASL restricts commercial use |
| MACE-MH/OMAT models | pretrained weights | ASL (restricted) | no | do not use as default |
| ASE | atomistic toolkit | LGPL-2.1 | yes (import-only) | used as Calculator adapter layer |
| PyTorch | ML runtime | BSD-3-Clause | **yes** | required by MACE |
| OpenMM | engine | MIT (core/CPU) + LGPL (CUDA/OpenCL) | **yes** | classical structure-prep engine; import-only, LGPL GPU plugins dynamically linked |
| OpenMM-Torch | plugin | MIT | yes | optional, not required for MACE path |
| openff-toolkit | FF application / topology | MIT | **yes** | builds OpenFF Topology + applies Sage; classical prep only |
| openff-interchange | OpenMM System export | MIT | **yes** | Interchange → OpenMM System for classical prep |
| openff-forcefields (Sage 2.x) | classical FF parameters | CC-BY-4.0 | **yes (attribution req.)** | small-molecule FF for monomer/initiator prep; cite OpenFF in distributed outputs |
| openff-nagl (code) | GNN partial charges | MIT | **yes** | AM1-BCC surrogate; avoids AmberTools(GPL)/OpenEye for charge assignment |
| openff-nagl-models | charge model weights | CC-BY-4.0 | **yes (attribution req.)** | trained NAGL weights; cite OpenFF. Fallback: RDKit Gasteiger (BSD) |
| PFP / Matlantis-related backend | model/backend | blocked_pending_review | no | do not enable by default without explicit rights confirmation |
| orb-models (code) | MLIP framework | Apache-2.0 | **yes** | OrbMol calculator; optional `[orb]` dependency |
| OrbMol-v2 model | pretrained weights | Apache-2.0 | **yes** | trained on OMol25 + OPoly26 (polymer data); recommended for organic/polymer systems |
| nvalchemiops | PME electrostatics + D3 | Apache-2.0 | **yes** | NVIDIA ALCHEMI Toolkit-Ops (`nvalchemi-toolkit-ops`); required for periodic OrbMol-v2 (PME long-range Coulomb); LICENSE verified Apache-2.0 (2026-06-20). Windows torch.compile caveat → run periodic on Linux/WSL/cloud |
| aimnet (code) | MLIP framework | MIT | **yes** | AIMNet2 calculator (`pip install aimnet`; isayevlab/aimnetcentral); 周期境界 (Ewald/PME) ネイティブ対応; spin-charge 認識 |
| AIMNet2-NSE model | pretrained weights | MIT | **yes** | 開殻/ラジカル化学(spin-charge equilibration、total spin multiplicity 入力)。HF `isayevlab/aimnet2-nse`。ラジカル系の第2バックエンド候補(spike 2026-06-25)。OrbMol-v2 を置換でなく補完(クロス検証) |
| Provenzano 2025 Zenodo dataset (xlinker) | external reference data (+ code, unused) | CC-BY-4.0 | **yes (attribution req.)** | Zenodo DOI 10.5281/zenodo.11402476; E2 external comparison baseline (LAMMPS structures + xl_trend.txt). Record-level CC-BY-4.0 is the only stated license (no per-file code license); we READ data only — xlinker.py is not executed, modified, or vendored. Cite Provenzano et al., ACS Appl. Polym. Mater. 2025, 7(8), 4876 in outputs |
| GPUMD (code) | MD engine + NEP runtime | GPL-3.0-only | yes (external-process only) | NEP89 実行エンジン。GPL のため kagome へのリンク/ベンダリング禁止 — 別プロセス実行(入出力ファイル経由)なら商用利用可。Zenodo 10.5281/zenodo.15299684 も GPL v3.0 only |
| NEP89 model (nep89_20250409) | pretrained weights | blocked_pending_review | no | GPUMD リポジトリ `potentials/nep/nep89_20250409/` に同梱、重み単体の明示ライセンスなし。リポジトリの GPL-3.0 が及ぶ解釈が自然だが weights への明示宣言なし=方針により block。著者確認 or 明示ライセンス付与を待つ |
| calorine | NEP の Python/ASE calculator | MPL-2.0 | **yes** | gitlab.com/materials-modeling/calorine; ファイル単位コピーレフト、import 利用は商用安全。NEP89 を kagome バックエンド化する場合の推奨経路 |
| Toy backend | internal test backend | internal code | yes | required for open/public CI |
| matplotlib | plotting (optional) | PSF (permissive) | yes | optional `[plot]` dependency, not required at runtime |
| RDKit | cheminformatics | BSD-3-Clause | **yes** | SMILES → 3D via EmbedMolecule + MMFF; used in `[rdkit]` extra for vinyl/AIBN system builder |
| SciPy | numerics | BSD-3-Clause | **yes** | `curve_fit` for the Eq.11 conversion fit; optional `[fit]` extra |
| PyYAML | config parsing | MIT | **yes** | parses configs and this license registry in the check scripts |
| openff-units | units (OpenFF stack) | MIT | **yes** | unit handling in classical prep; imported via the `openff` namespace; verify upstream LICENSE before release |

## Policy
- Software license and model-weight license must both be acceptable.
- SaaS access permission is not the same as model redistribution permission.
- If evidence is incomplete, block the component by default.
- CC-BY-4.0 components (openff-forcefields, openff-nagl-models) are commercial-safe but require attribution: distributed results that depend on them must cite the Open Force Field Initiative.
- License evidence verified 2026-06-14 from upstream LICENSE files: openff-toolkit/interchange/nagl = MIT; openff-forcefields + openff-nagl-models = CC-BY-4.0; OpenMM core = MIT, CUDA/OpenCL = LGPL.
- nvalchemiops verified 2026-06-20 from upstream LICENSE (github.com/NVIDIA/nvalchemi-toolkit-ops): SPDX-License-Identifier Apache-2.0, © NVIDIA CORPORATION. Unblocks periodic OrbMol-v2 PME for paper-scale runs (Linux/WSL/cloud; Windows torch.compile caveat persists).
- Provenzano 2025 Zenodo record verified 2026-07-12 (https://zenodo.org/records/15418263): record license = Creative Commons Attribution 4.0 International; archive downloaded and inspected, no per-file LICENSE inside xlinker.tgz. Data-read-only usage for E2 comparison; attribution required in distributed figures/reports.
- NEP89/GPUMD/calorine verified 2026-07-21 (NEP89 適用可能性検討, arXiv:2504.21286): GPUMD GitHub LICENSE = GPL-3.0, Zenodo record (10.5281/zenodo.15299684) = "GNU General Public License v3.0 only"; calorine GitLab LICENSE = MPL-2.0。NEP89 重み (`potentials/nep/nep89_20250409/nep89_20250409.txt`) は GPL リポジトリ同梱だが**重み自体への明示ライセンス宣言・README なし**、公式モデル配布先 gitlab.com/brucefan1983/nep-data もライセンス未確認 → 方針「evidence が不完全なら block」により `blocked_pending_review`。解除条件: (a) 著者への利用条件確認、または (b) 重みへの明示ライセンス付与の確認。訓練データは「論文出版時公開予定」段階(推論利用には不要)。GPL-3.0 の GPUMD は外部プロセス実行(ファイル I/O 経由)に限定すれば商用利用可、kagome コードへのリンク・同梱は不可。
- aimnet (code) + AIMNet2-NSE weights verified MIT 2026-06-25: github.com/isayevlab/aimnetcentral LICENSE = MIT; HF model card `isayevlab/aimnet2-nse` declares `license: mit`。コード・重みとも商用安全。ラジカル/開殻化学のための第2バックエンド候補として spike 評価中(PES 検証 + 多ラジカル高スピン安定性)。OrbMol-v2 の代替でなく補完(backend-agnostic クロス検証; specs/decisions.md 2026-06-25 参照)。
