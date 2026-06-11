# Dependency and commercial-use matrix

| Component | Category | License / status | Commercial-use default? | Notes |
|---|---|---|---|---|
| MACE (code) | MLIP framework | MIT | **yes** | default uMLIP backend |
| MACE-MP-0 model | pretrained weights | MIT | **yes** | 89-element universal potential, default model |
| MACE-OFF23 model | pretrained weights | ASL (restricted) | no | organic chemistry focus, ASL restricts commercial use |
| MACE-MH/OMAT models | pretrained weights | ASL (restricted) | no | do not use as default |
| ASE | atomistic toolkit | LGPL-2.1 | yes (import-only) | used as Calculator adapter layer |
| PyTorch | ML runtime | BSD-3-Clause | **yes** | required by MACE |
| OpenMM | engine | MIT + LGPL (GPU) | yes | optional, not required for MACE path |
| OpenMM-Torch | plugin | MIT | yes | optional, not required for MACE path |
| PFP / Matlantis-related backend | model/backend | blocked_pending_review | no | do not enable by default without explicit rights confirmation |
| Toy backend | internal test backend | internal code | yes | required for open/public CI |
| matplotlib | plotting (optional) | PSF (permissive) | yes | optional `[plot]` dependency, not required at runtime |

## Policy
- Software license and model-weight license must both be acceptable.
- SaaS access permission is not the same as model redistribution permission.
- If evidence is incomplete, block the component by default.
