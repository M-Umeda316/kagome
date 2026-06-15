# Task backlog

## Phase 1: Framework correctness (DONE)
- [x] Extract all equations and symbols from the paper into `paper/claims.yaml`
- [x] Implement `src/boost/tdbb.py` (Eq. 2-5, 8)
- [x] Implement `src/reactive/groups.py` (Eq. 6)
- [x] Implement `src/reactive/selection.py` (Eq. 7)
- [x] Implement `src/workflows/polymerization.py` (biased/unbiased loop)
- [x] Add run manifest schema (`src/workflows/manifest.py`)
- [x] Add toy backend for testing (`src/backends/toy.py`)
- [x] Add smoke test script (`scripts/run_smoke.py`)

## Phase 2: Reproduction readiness (DONE)
- [x] Add trajectory JSONL writer (`src/io/trajectory.py`)
- [x] Add trajectory reader (`src/io/readers.py`)
- [x] Integrate trajectory output into workflow
- [x] Add figure regeneration script (`scripts/reproduce_figures.py`)
- [x] Add matplotlib as optional `[plot]` dependency

## Phase 3 preparation (DONE)
- [x] Add atomic masses support to SimulationState
- [x] Extract integrator protocol + VelocityVerlet (`src/integrators/verlet.py`)
- [x] Add Langevin thermostat (`src/integrators/langevin.py`)
- [x] Add bond event tracking (`src/reactive/bonds.py`)
- [x] Add reactive group updates after confirmed reactions
- [x] Implement conversion tracking (Eq. 11-12) (`src/analysis/conversion.py`)
- [x] Implement depth-resolved reaction density (Eq. 13) (`src/analysis/density.py`)

## Phase 3: Scientific reproduction
- [x] Validate commercial status of target backends → MACE-MP-0 (MIT), ASE (LGPL), PyTorch (BSD)
- [x] Add ASE calculator adapter (`src/backends/ase_adapter.py`)
- [x] Add MACE-MP-0 backend (`src/backends/mace_backend.py`)
- [x] Run MACE+TDBB end-to-end (`scripts/run_mace.py`)
- [x] Generate figures from real MLIP output
- [x] Document backend choice and deviations in `specs/decisions.md`
- [ ] Run paper-faithful scale (2000-step phases, larger systems) for trend matching
- [ ] Add MACE-OFF or fine-tuned model for organic polymer accuracy

## Phase 4: OrbMol-v2 backend (DONE)
- [x] Validate commercial status: orb-models (Apache-2.0), weights (Apache-2.0)
- [x] Add OrbMol-v2 backend (`src/backends/orb_backend.py`)
- [x] Add orb optional dependency in pyproject.toml
- [x] Add backend unit test (`tests/unit/test_backends.py::TestOrbBackend`)
- [x] Fix unit system (FORCE_CONV), integrator pre/post split, MIC, pair key normalization
- [x] Run OrbMol-v2+TDBB end-to-end (`scripts/run_orb.py`)
- [x] Generate figures from OrbMol-v2 output → biased/unbiased energy pattern confirmed
- [x] Document OrbMol-v2 decision in `specs/decisions.md`
- [x] Update `specs/dependency-license-matrix.md`
- Note: PBC disabled (cell=None) — nvalchemiops triggers torch.compile on Windows (blocked_pending_review)

## Phase 5: Robustness and analysis wiring (DONE)
- [x] A1: Fix molecule placement overlap detection (atom-vs-atom, RuntimeError on failure)
- [x] A1: Extract shared ethylene builders to `scripts/_systems.py`
- [x] A2: Add instantaneous temperature to TrajectoryFrame and workflow output
- [x] B1: Add `read_bond_events()` to `src/io/readers.py`
- [x] B1: Wire α(t) conversion plot and temperature plot into reproduce_figures.py
- [x] B3: Document dissociation r0 assumption in specs/decisions.md
- [x] C1: Fix EV_TO_KCAL_MOL precision (23.0609 → 23.060548, NIST CODATA 2018)
- [x] C2: Remove dead TYPE_CHECKING block in orb_backend.py
- [x] C3: Add logger.warning for unknown element fallbacks in masses_from_species / _build_pair_biases
- [x] C4: De-duplicate temperature_K in run_orb / run_mace (read from LangevinParams)
- [~] B2: Bond formation demonstration — f1_max=250 reached, machinery verified, confirmed_formation=0 (expected: ethylene barrier + non-periodic diffusion). Paper-scale run (2000 steps, PBC) deferred until nvalchemiops license resolved. See decisions.md.

## Phase 6: Classical structure preparation (OpenMM/OpenFF) — decouple prep from ML production
Rationale: the all-ML pipeline spent ~1-2k expensive ML evals on packing/densification at the slowest density, making paper200 infeasible on the 16 GB GPU (throughput, not VRAM). Move structure prep to a cheap classical FF; reserve OrbMol-v2 for TDBB production. See specs/decisions.md "2026-06-14: Decouple initial-structure preparation".
- [x] P0a: License-verify + register openmm, openff-toolkit/interchange/forcefields/nagl(+models) in dependency-license-matrix.md and approved_dependencies.yaml
- [x] P0b: Record decisions D-1 (0.5 g/mL fixed prep), D-2 (nagl charges + Gasteiger fallback), D-3 (simple protocol), D-4 (self-contained OpenMM prep) in decisions.md
- [x] P0c: Correct the paper200 record (launched run did not complete; no artifacts)
- [x] P1a: Refactor `_systems._rdkit_mol(smiles, seed)` shared by `_rdkit_3d` and prep (atom-order guarantee)
- [x] P1b: Scaffold `src/prep/` (structure_io + openmm_equilibrate config/body) + `--load-structure` handoff in run_vinyl_aibn.py; Å↔nm / kcal↔kJ constants in src/units.py; tests/unit/test_prep.py (12 cases). OpenMM/OpenFF body lands but is validated in P2.
- [x] P2a: Build OpenFF Topology in builder order → Sage + charges → OpenMM System (validated in WSL)
- [x] P2b: Protocol: minimize → compress 0.25→0.5 g/mL → NVT thermalize → return (positions, cell); box_vectors set pre-interchange for PBC
- [x] P2c: scripts/prep_structure.py entry point (runs in WSL prep env) + run_vinyl_aibn.py `--load-structure` handoff
- [x] P2d: ENV PIVOT — OpenFF unusable on native Windows (MKL gemm DLL fault); prep runs in WSL conda env `pfpoly-prep`, production stays on Windows `pfpoly-gpu`, handoff via PreparedStructure JSON on /mnt/c. Validated 40+2 end-to-end. See decisions.md 2026-06-14 "Classical prep runs in WSL".
- [x] P3a: Unit tests (atom-order/species match, unit round-trip, tiny-system smoke) — tests/unit/test_prep.py (12 cases). OpenFF body validated by WSL integration run.
- [x] P3b: Integration on 40+2 (prep completes; handoff into GPU production passes species assert; FIRE E=-67.5; 6 candidates)
- [x] P4a: 40+2 prep + handoff validated end-to-end
- [x] P4b: paper-scale run — 200+10 exceeds 16 GB VRAM (OrbMol-v2 ~9.5 GB/call), so ran **100+5** instead: completed 14,000 steps NVT at 0.50 g/mL, T mean 295 K, VRAM 2.4-4.7 GB, candidates 14/11/19. confirmed_formations=0 (isolated to TDBB bias capture range — Ask-first). Figures generated; figure-comparison.md updated.
- [x] P5: ML compress_box is now only the legacy fallback (paper-scale uses WSL classical prep + --load-structure); reproduce_figures cp932 fix; docs in decisions.md/figure-comparison.md; commits per phase.

### Reproduction recipe (paper-density vinyl, this machine)
1. Prep in WSL: `wsl -d Ubuntu-24.04 -- bash -lc 'cd /mnt/c/Users/shanu/Documents/Python/pfpoly && ~/miniconda3/envs/pfpoly-prep/bin/python scripts/prep_structure.py --n-monomers 100 --n-initiators 5 --seed 42 --charge-method gasteiger --platform CPU --compress-relax-steps 0 --nvt-steps 20000 --output runs/prep/paper100.json'`
2. Produce on Windows GPU: `<pfpoly-gpu>/python scripts/run_vinyl_aibn.py --load-structure runs/prep/paper100.json --n-monomers 100 --n-initiators 5 --seed 42 --backend orb --device cuda --no-barostat --n-cycles 3 --biased-steps 2000 --unbiased-steps 2000 --equil-steps 2000 --output-dir runs/vinyl_aibn_paper100`
3. Figures: `scripts/reproduce_figures.py --trajectory ...trajectory.jsonl --bonds ...bonds.jsonl --n-reactive-sites 205 --target-temperature 333 --output-dir .../figures`

### RESOLVED (2026-06-15): formations=0 diagnosed; first bond demonstrated
- PDF Table S1 confirms the [3,6] window is paper-correct (earlier "window bug" hypothesis withdrawn).
- Fix A implemented: in-phase reaction detection during the biased phase + run-until-reaction (paper §2.2 step 3).
- OrbMol-v2 reproduces the radical addition (relaxed PES scan: barrier +6, product −28 kcal/mol; scripts/scan_radical_addition.py).
- First confirmed bond: scripts/demo_radical_formation.py forms+retains a C–C bond (final 1.62 Å). Full TDBB chain validated.
- formations=0 in melt runs was: closed-shell surrogate (no radical channel) + sampling (bias has ~0 force beyond ~2.6 Å, so pairs must diffuse inward). See specs/decisions.md 2026-06-15 entries.

## Phase 7: Toward polymerization reproduction — see specs/handoff-plan-v4.md
Remaining scope S1–S6 (single-chain propagation → melt sampling → multi-radical spin → paper fidelity → figures → hardware). Detailed roadmap and the S1 work plan are in **specs/handoff-plan-v4.md**.
- [ ] S1: single-chain propagation demo (1 radical + N monomers, doublet throughout, propagation_events>=2) — NEXT
- [ ] S2: melt-driven formations from the [3,6] window (run-until-reaction length, cycles)
- [ ] S3: multi-radical system-spin handling for OrbMol-v2 (open modeling question)
- [ ] S4: paper-fidelity (ij+ik+jl multi-pair criterion; AIBN decomposition / Activation)
- [ ] S5: figures — conversion α(t) (Eq.11), Carothers, depth-resolved density (Eq.12) vs paper
- [ ] S6: hardware for full 200+10 (>=24 GB GPU) or stay at 100+5
