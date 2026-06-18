# Paper notes

arXiv:2511.22874, Mori et al. (2025-11-28)
Last verified: 2026-06-18 (PDF cross-checked for eq numbering, reaction templates, and α definition)

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
Numbering confirmed from PDF (2026-06-18). Eq.11 is the last numbered equation.
- Eq. 9 (PDF p.8): Rp = -d[M]/dt = kp·[M]·[P*] — polymerization rate (kinetics theory comparison)
- Eq. 10 (PDF p.8): Rp ∝ [M]·[I] — simplified rate (Carothers comparison)
- Eq. 11 (PDF p.9): α(t) = 1 - exp(-k*_p·t) — exponential fit to conversion curve (post-hoc only)
- Unnumbered (PDF p.12): ρ_rxn(z) = N_rxn(z) / (A·Δz·N_frames) — depth-resolved reaction density
- Unnumbered (PDF p.9 Fig.2 caption): α = 1 − [M]/[M]₀ — raw monomer conversion definition

## Reaction-selection workflow
1. Define reactive atom groups G_I, G_J, G_K, G_L (chemically equivalent atoms).
2. Enumerate candidate tuples (i,j,k,l) satisfying distance bounds on pair set P.
3. Score each candidate: d_ijkl = r_ij + r_ik + r_jl.
4. Sort ascending, greedily select non-overlapping (no shared atoms).

## Reaction templates per system (Table S1/S2, PDF p.21-22)
| System | Groups | Pair set P (scoring/identification) | Bias-only pairs | Notes |
|---|---|---|---|---|
| Vinyl radical polymerization | 4: Gi (radical_C), Gj (vinyl_alpha_C), Gk (chain_C), Gl (vinyl_beta_C) | {(i,j) [3,6], (i,k) [0,3], (j,l) [0,3]} | — | Table S1: i-j V^f, i-k/j-l constraint only (no bias) |
| Nylon-6,6 (step-growth) | 4: Gi (amine_N), Gj (carboxyl_C), Gk (amine_H), Gl (carboxyl_OH) | {(i,j) [3,6], (i,k) [0,3], (j,l) [0,3]} | k-l (H-OH water formation, V^f) | Table S2: d_ijkl = r_ij+r_ik+r_jl (3 terms); k-l is bias-only, not in score |
| Epoxy curing (DGEBA+DETA on CuO) | 4: Gi (epoxy O), Gj (1° amine N), Gk (2° amine N), Gl (surface OH) | {(i,j),(i,k),(j,l)} | — | 4-group full template |

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

## System compositions (Supporting Information, confirmed PDF S-3..S-5)
| System | Composition | Density | Ensemble | T |
|---|---|---|---|---|
| Vinyl radical | **200 monomer + 10 AIBN** (20:1) | **0.5 g/mL** initial | NPT, 1 atm | 333 K |
| Nylon-6,6 | 100 diamine + 100 diacid (equimolar) | 0.5 g/mL initial | NPT, 1 atm | 300 K |
| Epoxy/CuO | 100 DGEBA + 50 DETA on CuO(001) 8×8×6 slab | 0.2→1.0 g/cm³ (wall-compressed) | NVT | 333 K |

- All systems: 3 independent runs, results averaged. Solvent-free unless noted (e.g. styrene/toluene).
- **Key mechanism (PDF p.7, S-7): "polymerization is governed by near-contact events"** — the
  formation bias (f2=10 Å⁻², r0≈2 Å) only acts strongly within ~r0±0.5 Å. Candidates are listed at
  3–6 Å (Table S1) but bond formation occurs when thermal motion brings a listed pair to near-contact.
  This REQUIRES paper-scale density (0.5 g/mL) and many molecules (~200) so near-contact events are
  frequent. Small/dilute systems (≪200 mol, density <0.5) produce few or zero confirmed formations
  even with correct parameters. See specs/decisions.md "2026-06-13: T-G1a root-cause".

## Figures to reproduce
- Fig. 2-6: qualitative trend matching (see specs/acceptance-criteria.md)
- Fig. S4: sensitivity analysis showing gamma as global scaling factor (not system-specific)
