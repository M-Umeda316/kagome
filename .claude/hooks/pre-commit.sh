#!/usr/bin/env bash
set -euo pipefail

python scripts/validate_configs.py
python scripts/check_dependency_licenses.py --approved specs/approved_dependencies.yaml
pytest -q tests/unit
