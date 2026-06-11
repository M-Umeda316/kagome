#!/usr/bin/env bash
set -euo pipefail

python scripts/check_seed_defined.py "$@"
python scripts/check_output_path.py "$@"
python scripts/check_dependency_licenses.py --approved specs/approved_dependencies.yaml
