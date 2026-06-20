# Handoff Plan v4 — Remaining scope toward vinyl radical polymerization reproduction

Date: 2026-06-15
Paper: arXiv:2511.22874, Mori et al. (paper/2511.22874v1.pdf)
Supersedes the forward-looking parts of handoff-plan-v3. Diagnostic history is in
specs/decisions.md (2026-06-14 / 2026-06-15 entries).

## Where we are (validated this session)
- Infra blocker resolved: paper-density (0.5 g/mL) TDBB runs end-to-end on the
  RTX 4060 Ti via **WSL classical prep (OpenFF/Sage) → file handoff → NVT
  production (OrbMol-v2)**. (decisions.md "Classical prep runs in WSL", "paper100
  completes".)
- formations=0 fully diagnosed: NOT density/scale/temperature/window (the [3,6]
  window is paper-correct, Table S1). Causes were (1) detection only at end of the
  unbiased phase — fixed (Fix A: in-phase detection + run-until-reaction); (2) the
  closed-shell radical surrogate had no radical chemistry.
- OrbMol-v2 reproduces the radical addition: relaxed PES scan gives a textbook
  profile (barrier +6, product −28 kcal/mol). (scripts/scan_radical_addition.py.)
- FIRST confirmed bond: scripts/demo_radical_formation.py forms and RETAINS a
  radical-addition C–C bond (final 1.62 Å after unbiased relaxation). Full TDBB
  chain validated: select → bias → in-phase vdW detection → confirm → bond
  survives relaxation.

So the mechanism works. What remains is making it a **melt-scale chain
polymerization** and matching the paper's figures.

## Remaining scope (priority order, with dependencies)

| ID | Item | Depends on | Effort | Risk |
|----|------|-----------|--------|------|
| S1 | Single-chain propagation demo (1 radical + N monomers, doublet throughout, ≥2–3 successive additions) | — | Med | Low |
| S2 | Melt sampling: formations from the [3,6] window via diffusion (run-until-reaction length, cycles, 0.5 g/mL) | S1 | Med | Med (GPU time) |
| S3 | Multi-radical system spin for OrbMol-v2 (the main open modeling question) | — (scout early) | High | High |
| S4 | Paper fidelity: multi-pair ij+ik+jl criterion (Table S1); AIBN decomposition (Activation, V^d) | S3 | Med | Med |
| S5 | Figures/validation: conversion α(t) (Eq.11), Carothers, depth-resolved density (Eq.12) vs paper | S1–S2 | Med | Low |
| S6 | Hardware for full 200+10 (>=24 GB GPU) or stay at 100+5 | — | Small | Low |

**Critical path to "polymerization reproduced at 100+5":** S1 → S2 → S5.
**S3** is required for paper-parallel multi-chain / full scale and is the most
uncertain — scout its feasibility early, in parallel. **S4** raises fidelity.
**S6** only if full 200+10 is required.

Key constraints to respect (do NOT change without paper re-read + owner approval):
TDBB params f2=10, f1_max=250/125, γ=1.0, r0=λ·Σr_vdw (λ=0.6); production candidate
window [3,6] (Table S1). All paper-confirmed.

---

# S1 work plan — single-chain propagation demo

## Goal
Demonstrate **chain propagation** end-to-end: one open-shell radical initiates and
then adds successive methyl-acrylate monomers, so the active radical migrates along
a growing chain. Target: **propagation_events ≥ 2** (a trimer or longer), with each
new C–C bond stable after unbiased relaxation.

## Why S1 first
Along a single growing chain the unpaired electron count is **always exactly 1**
(the radical sits at the chain end). So the whole run is a **doublet (spin=2)**
throughout, which **sidesteps the multi-radical system-spin problem (S3)** while
exercising the heart of polymerization (propagation + group/topology update). It
builds directly on the validated demo_radical_formation.py.

## Design
- System: 1 radical (`C[C](C)C#N`) + N monomers (start N≈6–8). Box/density modest
  so monomers sit near the growing chain end (a near-0.5 g/mL small box, or
  pre-arrange monomers around the radical for a deterministic demo).
- Backend: OrbMol-v2, **spin=2 held constant** (1 unpaired e per chain end at all
  times — document this invariant).
- Workflow: existing PolymerizationWorkflow with in-phase detection (Fix A) and the
  propagation_map. After a formation, `_update_groups_after_cycle` removes the
  reacted α-C and adds the monomer's β-C to the `radical_C` group (the new chain-end
  radical); the next cycle should then bias (new radical, next monomer α-C).
- Selection: for a clean demonstration, reuse the demo's widened window
  (`--select-rmin 1.5`) and/or pre-position the next monomer near the chain end so
  each cycle has a selectable near-contact pair. Clearly label as a DEMO device
  (production window stays [3,6]); melt-driven selection is S2.

## Tasks
- **T-S1.1 (trace/verify bookkeeping)**: Confirm in `_update_groups_after_cycle`
  that for a confirmed formation the monomer α-C is removed from `vinyl_alpha_C`,
  the reacted radical is removed from `radical_C`, and the new β-C is appended to
  `radical_C` (chain-end migration). Confirm the NEXT cycle's `find_candidates`
  then pairs the new radical with another monomer's α-C. Add a focused unit test
  (toy backend) asserting the group/propagation update across two successive
  formations. File: tests/unit/test_workflow.py (or test_propagation.py).
- **T-S1.2 (builder/placement)**: New `scripts/demo_chain_propagation.py` (adapt
  demo_radical_formation.py): build 1 radical + N monomers; optionally pre-arrange
  monomers as a loose row near the radical lobe so successive additions are
  geometrically reachable. Keep spin=2.
- **T-S1.3 (run)**: Run multiple cycles (n_cycles ≈ N, biased run-until-reaction +
  short unbiased). After each cycle, log: which pair formed, new radical index,
  current chain length, final pair distance. Expect ≥2 successive additions.
- **T-S1.4 (spin invariant)**: Assert/verify the system stays doublet — exactly one
  radical end at all times. Decide and document that spin is held at 2 for S1
  (single chain). Note where multi-chain (S3) would diverge.
- **T-S1.5 (artifacts/figures)**: Save bonds.jsonl/trajectory; generate a simple
  conversion-vs-cycle figure for the single chain; record results in
  specs/decisions.md and a row in specs/figure-comparison.md. Commit.

## Acceptance criteria
- `propagation_events >= 2` (≥2 successive additions; a trimer or longer).
- Each new C–C bond has r < r0 (=2.04 Å) and remains bonded after the unbiased
  phase (no spring-back).
- The active radical correctly migrates to each new chain end (verified via logs +
  unit test).
- Spin remains doublet throughout; the invariant is documented.
- Reproduction command recorded.

## Risks / open questions
- **Sampling vs determinism**: successive additions may not occur within feasible
  cycles if monomers must diffuse in. For the DEMO, pre-positioning the next
  monomer near the chain end (or widened window) is acceptable and clearly labeled;
  melt-driven, undirected formation is S2.
- **Steric access** to the growing chain end after a few additions — keep N small
  first (trimer) then extend.
- **Group-update correctness** is the main code risk — covered by T-S1.1 test.
- **Spin** stays 2 only because it is a single chain; this is the explicit boundary
  with S3 (multiple simultaneous radicals).

## Reproduction command (planned)
```
# from repo root, in pfpoly-gpu (production env)
<pfpoly-gpu>/python scripts/demo_chain_propagation.py \
    --device cuda --n-monomers 8 --n-cycles 8 \
    --select-rmin 1.5 --output-dir runs/demo_chain_propagation
```
(WSL classical prep is not required for this small isolated single-chain demo;
build + place directly. spin=2 is set inside the script.)

---

# S2 work plan — melt-driven formation from the paper [3,6] window

## Goal
Demonstrate that confirmed formations occur in a real **melt** WITHOUT the S1 demo
devices — i.e. with the paper-faithful candidate window **[3,6] Å (Table S1)** and
**no directed monomer placement**. Selected pairs must diffuse from the [3,6]
window across the ~6 kcal/mol addition barrier into the bias-capture shell
(<~2.6 Å), where the bias drives them to the product. Target: **>=1 melt-driven
confirmed_formation** (ideally a few / short propagation) from [3,6].

## Scope decision — single radical, to isolate sampling from S3
A melt-driven test needs REAL radical chemistry (the closed-shell surrogate used in
paper100/paper100_fixA can never bond → those zeros were expected). But a dense melt
with MANY radicals hits the multi-radical system-spin problem (S3). So S2 uses
**ONE open-shell radical in a dense methyl-acrylate melt** → spin=2 (doublet)
throughout, S3 avoided. This cleanly isolates the question S2 exists to answer:
*does a pair selected from [3,6] reach bonding distance by melt diffusion + the
weak inward bias slope within a feasible biased-phase length / cycle count?*
(The single-bond demo + PES scan say it is physically possible once a pair reaches
the capture shell; S2 measures/achieves the RATE.)

## Physics levers (paper-allowed; do NOT touch f2/r0/λ/window)
- The bias V^f decreases monotonically from ~250 (at 3.5 Å) to 0 (at r0=2.04), so it
  adds a WEAK inward thermodynamic slope across [2.04, 3.5] that strengthens as the
  pair approaches — over a long biased phase this nudges pairs inward.
- **Biased-phase max length** (run-until-reaction already implemented; raise the cap).
- **Number of cycles** (the paper builds chains over many cycles).
- **γ** (acceleration factor): the paper explicitly tunes reaction rate via γ
  (Fig. S4); larger γ ramps f1 faster → higher rate. This is the legitimate knob.
- **Density** (0.5 g/mL paper; higher density → more near-contacts).

## Design / system preparation (two options)
- **Option A (recommended first, fast):** build 1 radical (`C[C](C)C#N`) + N monomers
  (N≈30–50) directly via the grid placer at the densest feasible density for this
  count (the placer handled 42 molecules at gpu40; try 0.5, fall back to ~0.35),
  then a short ML equilibration (existing `equil_steps`) to relax. No WSL/OpenFF
  prep → avoids classical parametrization of an open-shell radical. spin=2.
- **Option B (fidelity follow-up):** WSL classical prep at 0.5 g/mL using the
  CLOSED-SHELL initiator for packing only, then convert that initiator to the
  radical in production by removing its radical H (atom-count/index/group fix-up)
  and setting spin=2. Truer to the paper density but needs a small H-removal utility.

## Tasks
- **T-S2.1**: single-radical dense melt builder/runner (Option A): extend
  run_vinyl_aibn.py (already has --initiator-smiles/--spin) or a small script;
  1 radical + N monomers, spin=2, ML equilibration, NVT, paper window [3,6]
  (NO --select-rmin widening), NO directed placement.
- **T-S2.2**: run TDBB with run-until-reaction and a LONG biased cap (e.g. 5000–10000
  steps) over MANY cycles (e.g. 20–50); log per-cycle reaction events and the
  min radical–vinyl distance reached (to see how close pairs get).
- **T-S2.3**: if formations are too rare, sweep the paper-allowed levers — γ (e.g.
  1→2→4), biased-phase length, density — and report reaction count vs each
  (mirrors the paper's Fig. S4 sensitivity study). Record what is needed for a
  feasible melt rate on this hardware.
- **T-S2.4**: produce a conversion-like metric (cumulative formations vs cycle) and
  compare qualitatively to the paper (monotonic rise; α(t) shape — full Eq.11 fit
  is S5).
- **T-S2.5**: record results in decisions.md + figure-comparison.md; commit.

## Acceptance criteria
- **>=1 melt-driven confirmed_formation** with the paper window [3,6] and no directed
  placement / no widened window (single radical, spin=2).
- Logged evidence that selected [3,6] pairs reach <r0 during biasing (min-distance
  trace), i.e. the mechanism works from the paper window, not just from pre-placed
  contacts.
- The biased-phase length / cycles / γ needed for a feasible rate are recorded
  (quantifies the sampling cost on this hardware).
- Reproduction command recorded.

## Risks / open questions
- **Rate / GPU time**: melt-driven crossing may be rare → long runs (hours). Mitigate
  with higher density, longer biased caps, larger γ; cap total wall-clock and report
  the achieved rate even if low.
- **Open-shell radical in classical prep** (Option B only) — parametrization;
  Option A avoids it.
- If even with aggressive (but paper-allowed) γ / length the melt rate is impractical
  on this GPU, that is itself a finding (quantifies why directed placement was needed
  in S1); document and consider scale/hardware (S6).
- Stays single-radical: multi-radical melts (paper-parallel) remain S3.

## Reproduction command (planned, Option A)
```
<pfpoly-gpu>/python scripts/run_vinyl_aibn.py \
    --n-monomers 40 --n-initiators 1 --initiator-smiles "C[C](C)C#N" --spin 2 \
    --density 0.5 --backend orb --device cuda --no-barostat \
    --n-cycles 30 --biased-steps 6000 --unbiased-steps 500 --equil-steps 2000 \
    --output-dir runs/s2_melt_single_radical
```
(window stays [3,6]; no --select-rmin. If placement at 0.5 fails, lower --density or
use WSL prep per Option B.)
