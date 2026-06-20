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
| nvalchemiops | PME electrostatics | blocked_pending_review | no | NVIDIA; required only for periodic OrbMol-v2 with long-range Coulomb; license unconfirmed |
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
