"""Validate all YAML config files under configs/.

Checks that required keys are present and values are within expected ranges.
Exits non-zero if any config fails validation.

Usage:
    python scripts/validate_configs.py
    python scripts/validate_configs.py --config-dir configs/
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    import yaml
    _HAS_YAML = True
except ImportError:
    _HAS_YAML = False


BOOST_REQUIRED_KEYS = {
    'boost': {
        'timestep_fs': (float | int, lambda v: v > 0, 'must be > 0'),
        'biased_steps': (int, lambda v: v > 0, 'must be > 0'),
        'unbiased_steps': (int, lambda v: v > 0, 'must be > 0'),
        'lambda_vdw': (float | int, lambda v: 0.0 < v <= 1.0, 'must be in (0, 1]'),
        'f2_invA2': (float | int, lambda v: v > 0, 'must be > 0'),
        'gamma': (float | int, lambda v: v > 0, 'must be > 0'),
        'f1_max_form_kcal_mol': (float | int, lambda v: v > 0, 'must be > 0'),
        'f1_max_break_kcal_mol': (float | int, lambda v: v > 0, 'must be > 0'),
    }
}

EVAL_REQUIRED_KEYS = {
    'eval': {
        'name': (str, lambda v: len(v) > 0, 'must be non-empty'),
        'save_interval': (int, lambda v: v >= 0, 'must be >= 0'),
    }
}


def _load_yaml(path: Path) -> dict:
    text = path.read_text(encoding='utf-8')
    if _HAS_YAML:
        return yaml.safe_load(text) or {}
    raise ImportError(
        'PyYAML is required to validate configs. '
        'Install with: pip install pyyaml'
    )


def _validate_section(
    data: dict,
    section: str,
    rules: dict,
    path: Path,
    errors: list[str],
) -> None:
    if section not in data:
        errors.append(f'{path}: missing top-level key "{section}"')
        return

    section_data = data[section]
    for key, (expected_type, check_fn, check_desc) in rules.items():
        if key not in section_data:
            errors.append(f'{path}: [{section}] missing key "{key}"')
            continue
        val = section_data[key]
        if not isinstance(val, expected_type):
            errors.append(
                f'{path}: [{section}] "{key}" has type {type(val).__name__}, '
                f'expected {expected_type}'
            )
            continue
        if not check_fn(val):
            errors.append(f'{path}: [{section}] "{key}" = {val!r} — {check_desc}')


def validate_file(path: Path) -> list[str]:
    errors: list[str] = []
    try:
        data = _load_yaml(path)
    except Exception as exc:
        return [f'{path}: failed to parse — {exc}']

    parents = [p.name for p in path.parents]
    if 'boost' in parents or path.parent.name == 'boost':
        _validate_section(data, 'boost', BOOST_REQUIRED_KEYS['boost'], path, errors)
    elif 'eval' in parents or path.parent.name == 'eval':
        _validate_section(data, 'eval', EVAL_REQUIRED_KEYS['eval'], path, errors)
    else:
        # Unknown category: only check that the file is parseable YAML
        pass

    return errors


def main() -> None:
    parser = argparse.ArgumentParser(description='Validate config YAML files')
    parser.add_argument(
        '--config-dir',
        type=Path,
        default=Path('configs'),
        help='Directory containing config YAML files (default: configs/)',
    )
    args = parser.parse_args()

    if not args.config_dir.exists():
        print(f'ERROR: Config directory not found: {args.config_dir}', file=sys.stderr)
        sys.exit(1)

    yaml_files = list(args.config_dir.rglob('*.yaml')) + list(args.config_dir.rglob('*.yml'))
    if not yaml_files:
        print(f'No YAML files found in {args.config_dir}')
        sys.exit(0)

    all_errors: list[str] = []
    for f in sorted(yaml_files):
        all_errors.extend(validate_file(f))

    if all_errors:
        print('Config validation FAILED:', file=sys.stderr)
        for e in all_errors:
            print(f'  {e}', file=sys.stderr)
        sys.exit(1)

    print(f'OK: {len(yaml_files)} config file(s) validated successfully.')


if __name__ == '__main__':
    main()
