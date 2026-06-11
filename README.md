# Polymerization Reproduction Starter

This repository is a Claude Code starter for reproducing a polymerization/crosslinking workflow based on a time-dependent bond boosting method with commercial-use guardrails.

## Purpose
- Reproduce the paper's TDBB-centered workflow
- Keep the MLIP backend swappable
- Avoid default reliance on license-unclear model providers
- Make every result traceable to config, code, and seed

## Suggested first steps
1. Fill `paper/claims.yaml` from the paper PDF.
2. Confirm commercial status of candidate backends in `specs/dependency-license-matrix.md`.
3. Implement `src/boost/tdbb.py` first.
4. Add candidate selection logic in `src/reactive/selection.py`.
5. Wire a toy backend for integration testing.
6. Add OpenMM integration once interfaces stabilize.

## Minimal commands
```bash
python scripts/validate_configs.py
python scripts/check_dependency_licenses.py --approved specs/approved_dependencies.yaml
pytest -q tests/unit
```
