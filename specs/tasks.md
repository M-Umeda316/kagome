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
