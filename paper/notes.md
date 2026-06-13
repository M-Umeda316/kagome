# Paper notes

arXiv:2511.22874, Mori et al. (2025-11-28)
Last verified: 2026-06-13 (arXiv HTML version; PDF needed for final confirmation)

## Core claim
uMLIP + TDBB enables system-independent polymerization simulations without per-system parameter tuning.

## Key equations (implementation-critical)
- Eq. 2: V^f(r,t) = f1(t)(1 - exp(-f2(r-r0)^2)) — formation bias
- Eq. 3: V^d(r,t) = f1(t) exp(-f2 r^2) — dissociation bias
- Eq. 4: r0 = λ Σ r_a^vdw — target distance from vdW radii
- Eq. 5: f1(t) = min(γt, f1_max) — linear ramp with cap
- Eq. 6: G_X = {a | a ∈ X} — reactive group definition (X = I, J, K, L)
- Eq. 7: candidate selection with distance bounds and pair set P
- Eq. 8: ΔV = Σ_groups Σ_pairs [fp V^f + (1-fp) V^d] — total bias

## Post-hoc analysis equations (not in simulation loop)
Numbering from arXiv HTML — needs PDF cross-check (docstrings in src/ may cite different numbers).
- Eq. 9: Rp = -d[M]/dt = kp·[M]·[P*] — polymerization rate (kinetics theory comparison)
- Eq. 10: Rp ∝ [M]·[I] — simplified rate (Carothers comparison)
- Eq. 11: α(t) = 1 - exp(-kp_eff·t) — exponential fit to conversion curve (post-hoc only)
- Eq. 12: ρ_rxn(z) = N_rxn(z) / (A·Δz·N_frames) — depth-resolved reaction density
- (Unnumbered or different number): α = N_reacted / N_total — raw conversion fraction from bond counts

## Reaction-selection workflow
1. Define reactive atom groups G_I, G_J, G_K, G_L (chemically equivalent atoms).
2. Enumerate candidate tuples (i,j,k,l) satisfying distance bounds on pair set P.
3. Score each candidate: d_ijkl = r_ij + r_ik + r_jl.
4. Sort ascending, greedily select non-overlapping (no shared atoms).

## Reaction templates per system
| System | Groups | Pair set P | Notes |
|---|---|---|---|
| Vinyl radical polymerization | 2: Gi (radical C), Gj (alkene C) | {(i,j)} | Matches current 2-group implementation |
| Nylon-6,6 (step-growth) | 2: Gi (amine N), Gj (carboxylic acid C) | {(i,j)} | |
| Epoxy curing (DGEBA+DETA on CuO) | 4: Gi (epoxy O), Gj (1° amine N), Gk (2° amine N), Gl (surface OH) | {(i,j),(i,k),(j,l)} | 4-group full template |

## Simulation schedule (Section 3 Methods — confirmed)
Alternating: **2000 biased steps** → **2000 unbiased steps** → repeat.
Timestep: 0.25 fs. Each phase = 2000 × 0.25 fs = 500 fs.
Confirmed from PDF p.7: "alternating biased and unbiased dynamics every 2000 steps (500 fs)".

## Default hyperparameters (confirmed from arXiv HTML, Section 3 Methods)
| Parameter | Value | Unit | Confirmed? |
|---|---|---|---|
| timestep | 0.25 | fs | Yes (HTML) |
| biased_steps | 2000 | steps | Yes (HTML) |
| unbiased_steps | 2000 | steps | Yes (PDF p.7) |
| λ (lambda_vdw) | 0.60 | dimensionless | Partial — HTML says "0.6 (implied)", Eq. 4 |
| f2 | 10.0 | Å⁻² | Yes (HTML, Section 3 Methods) |
| γ (gamma) | 1.0 | kcal/(mol·step) | Yes — value confirmed (PDF p.7); unit not stated, maintained as kcal/(mol·step) |
| f1_max (formation) | 250 | kcal/mol | Yes (HTML, Section 3 Methods) |
| f1_max (dissociation) | 125 | kcal/mol | Yes (HTML, Section 3 Methods) |

### γ unit — resolved
PDF p.7 confirmed γ=1.0, unit not stated. Maintaining kcal/(mol·step) — saturation at step 250.
Fig. S4 (p.25-26) confirms γ acts as global scaling factor; unit choice does not distort relative trends.
See specs/decisions.md "2026-06-11: Units convention for gamma" for full rationale.

## Systems studied in the paper
- Radical polymerization: methyl acrylate, methacrylate, styrene, vinyl acetate, diphenylethylene, dimethyl itaconate
- Initiator: AIBN (azobisisobutyronitrile)
- Nylon-6,6: hexamethylenediamine + adipic acid (step-growth, Carothers comparison)
- Epoxy curing: DGEBA + DETA hardener on hydroxylated CuO surface (interface design)

## Figures to reproduce
- Fig. 2-6: qualitative trend matching (see specs/acceptance-criteria.md)
- Fig. S4: sensitivity analysis showing gamma as global scaling factor (not system-specific)
