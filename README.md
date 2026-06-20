# KAGOME

**Kinetic Accelerated Growth Orchestrated by Molecular Engine**

A reproducible, commercial-safe implementation of time-dependent bond boosting (TDBB) for polymerization and curing simulations with machine-learning interatomic potentials.

## Features
- TDBB workflow faithful to the original paper
- Swappable MLIP backends (OrbMol-v2, MACE-MP-0, ASE adapters, classical FF)
- Deterministic seeds and full run manifests for reproducibility
- Commercial-use guardrails with explicit dependency license tracking

## Installation

```bash
pip install -e ".[dev]"
```

## Quick start

```bash
# Smoke test with toy backend
python scripts/run_smoke.py

# Unit tests
pytest -q tests/unit

# License check
python scripts/check_dependency_licenses.py --approved specs/approved_dependencies.yaml
```

## Project structure

```
src/kagome/
  boost/          # Time-dependent bond boosting (TDBB)
  reactive/       # Reactive group definitions and candidate selection
  workflows/      # Polymerization and curing loops
  backends/       # MLIP/calculator adapters
  integrators/    # Velocity Verlet, Langevin, MC barostat, FIRE minimizer
  analysis/       # Conversion, density profiles, Carothers analysis
  chem/           # Molecule builders (RDKit)
  io/             # Trajectory I/O
  prep/           # Classical structure preparation (OpenMM/OpenFF)
scripts/          # Entry points, scans, figure generation
configs/          # Experiment configurations
specs/            # Requirements, decisions, dependency licenses
paper/            # Structured notes from the source paper
```

## License

See `specs/dependency-license-matrix.md` for third-party dependency status.
