# Dependency and commercial-use matrix

| Component | Category | License / status | Commercial-use default? | Notes |
|---|---|---|---|---|
| OpenMM | engine | review and pin exact upstream license before release | yes, after verification | candidate default engine |
| OpenMM-Torch | plugin | MIT | yes | good candidate for default bias implementation |
| PyTorch | runtime | verify exact distribution terms in release process | yes, after verification | runtime dependency |
| PFP / Matlantis-related backend | model/backend | blocked_pending_review | no | do not enable by default without explicit rights confirmation |
| Toy backend | internal test backend | internal code | yes | required for open/public CI |
| matplotlib | plotting (optional) | PSF (permissive) | yes | optional `[plot]` dependency, not required at runtime |

## Policy
- Software license and model-weight license must both be acceptable.
- SaaS access permission is not the same as model redistribution permission.
- If evidence is incomplete, block the component by default.
