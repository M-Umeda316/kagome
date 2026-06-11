# Paper notes

arXiv:2511.22874, Mori et al. (2025-11-28)

## Core claim
uMLIP + TDBB enables system-independent polymerization simulations without per-system parameter tuning.

## Key equations (implementation-critical)
- Eq. 2: V^f(r,t) = f1(t)(1 - exp(-f2(r-r0)^2)) — formation bias
- Eq. 3: V^d(r,t) = f1(t) exp(-f2 r^2) — dissociation bias
- Eq. 4: r0 = λ Σ r_a^vdw — target distance from vdW radii
- Eq. 5: f1(t) = min(γt, f1_max) — linear ramp with cap
- Eq. 6: G_X = {a | a ∈ X} — reactive group definition
- Eq. 7: candidate selection with distance bounds and pair set P
- Eq. 8: ΔV = Σ_groups Σ_pairs [fp V^f + (1-fp) V^d] — total bias

## Reaction-selection workflow
1. Define reactive atom groups G_I, G_J, G_K, G_L (chemically equivalent atoms).
2. Enumerate candidate tuples (i,j,k,l) satisfying distance bounds on pair set P={(i,j),(i,k),(j,l)}.
3. Score each candidate: d_ijkl = r_ij + r_ik + r_jl.
4. Sort ascending, greedily select non-overlapping (no shared atoms).

## Simulation schedule
Alternating: 2000 biased steps → 2000 unbiased steps → repeat.
Timestep: 0.25 fs. Each phase = 500 fs.

## Default hyperparameters
| Parameter | Value | Unit |
|---|---|---|
| timestep | 0.25 | fs |
| biased_steps | 2000 | steps |
| unbiased_steps | 2000 | steps |
| λ (lambda_vdw) | 0.60 | dimensionless |
| f2 | 10.0 | Å⁻² |
| γ (gamma) | 1.0 | kcal/(mol·step) |
| f1_max (formation) | 250 | kcal/mol |
| f1_max (dissociation) | 125 | kcal/mol |

## Figures to reproduce
- Fig. 2-6: qualitative trend matching (see acceptance-criteria.md)
