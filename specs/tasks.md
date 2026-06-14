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
- [ ] P2a: Build OpenFF Topology in builder order → Sage + charges → OpenMM System
- [ ] P2b: Protocol: minimize → compress 0.25→0.5 g/mL → NVT thermalize → return (positions, cell)
- [ ] P2c: Wire `--prep {none,openmm}` (+ protocol/density/charge flags) into run_vinyl_aibn.py
- [ ] P3a: Unit tests (atom-order/species match, unit round-trip, tiny-system smoke, compression target)
- [ ] P3b: Integration on 40+2 (completes, T stable, wall-clock vs all-ML)
- [ ] P4a: Run 40+2 E2E with classical prep; save artifacts
- [ ] P4b: Run paper200 (200+10) classical prep → ML production; check formations>=1, T≈333 K, VRAM/wall-clock; figures; update figure-comparison.md
- [ ] P5: Demote ML compress_box to fallback for paper-scale; docs; final commit
