"""Check that no blocked dependencies are imported in source code.

Reads specs/approved_dependencies.yaml and scans src/ and scripts/ for
imports of any 'blocked_pending_review' entries.  Also warns about
'review_required' entries that are actively imported.

Rationale:
  approved_dependencies.yaml lists ALL evaluated dependencies, including
  blocked ones (for documentation).  The check passes if blocked items
  are only listed (not imported).  It fails if a blocked item is found
  in the source code.

Usage:
    python scripts/check_dependency_licenses.py --approved specs/approved_dependencies.yaml
    python scripts/check_dependency_licenses.py --approved specs/approved_dependencies.yaml --src-dirs src scripts
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

try:
    import yaml
    _HAS_YAML = True
except ImportError:
    _HAS_YAML = False


# Map package install names to import names (when they differ)
_INSTALL_TO_IMPORT: dict[str, list[str]] = {
    'mace-torch': ['mace'],
    'orb-models': ['orb_models', 'orb'],
    'pytorch': ['torch'],
    'ase': ['ase'],
    'numpy': ['numpy', 'np'],
    'matplotlib': ['matplotlib', 'plt'],
    'openmm': ['openmm', 'simtk'],
    'openmm-torch': ['openmmtorch'],
    'nvalchemiops': ['nvalchemiops'],
    'pfp': ['pfp', 'matlantis'],
    'mace-off23': [],  # only specific model path, hard to detect by import
}


def _load_yaml(path: Path) -> dict:
    text = path.read_text(encoding='utf-8')
    if _HAS_YAML:
        return yaml.safe_load(text) or {}
    raise ImportError(
        'PyYAML is required for YAML parsing. Install with: pip install pyyaml'
    )


def _find_imports_in_source(src_dirs: list[Path]) -> set[str]:
    """Return set of all imported module names found in Python source files."""
    imported: set[str] = set()
    import_pattern = re.compile(
        r'^\s*(?:import|from)\s+([\w.]+)', re.MULTILINE,
    )
    for src_dir in src_dirs:
        if not src_dir.exists():
            continue
        for py_file in src_dir.rglob('*.py'):
            text = py_file.read_text(encoding='utf-8', errors='replace')
            for m in import_pattern.finditer(text):
                top_level = m.group(1).split('.')[0]
                imported.add(top_level)
    return set(imported)


def _import_names_for(package: str) -> list[str]:
    return _INSTALL_TO_IMPORT.get(package, [package.replace('-', '_')])


def check(approved_path: Path, src_dirs: list[Path]) -> int:
    data = _load_yaml(approved_path)
    active_imports = _find_imports_in_source(src_dirs)

    hard_blocked_and_used: list[str] = []
    review_and_used: list[str] = []
    review_warnings: list[str] = []

    for section in ('software', 'models'):
        for entry in data.get(section, []):
            name = entry.get('name', '<unnamed>')
            status = entry.get('status', '')
            reason = entry.get('reason', entry.get('required_evidence', ''))

            import_names = _import_names_for(name)
            is_imported = any(imp in active_imports for imp in import_names if imp)

            if status == 'blocked_pending_review':
                if is_imported:
                    hard_blocked_and_used.append(
                        f'  [{section}] {name} (blocked) is imported as '
                        f'{import_names} - {reason}'
                    )
            elif status == 'review_required':
                if is_imported:
                    review_and_used.append(
                        f'  [{section}] {name} (review_required, imported) - '
                        f'missing evidence: {reason}'
                    )
                else:
                    review_warnings.append(
                        f'  [{section}] {name} - missing evidence: {reason}'
                    )

    if review_warnings:
        print('WARNING: pending license reviews (not currently imported):')
        for w in review_warnings:
            print(w)
        print()

    if review_and_used:
        print('WARNING: imported dependencies with pending license reviews:')
        for w in review_and_used:
            print(w)
        print('  These must be reviewed and approved before release.\n')

    if hard_blocked_and_used:
        print('ERROR: blocked dependencies are imported in source code:', file=sys.stderr)
        for b in hard_blocked_and_used:
            print(b, file=sys.stderr)
        print(
            '\nRemove these imports or resolve the license status in '
            'specs/approved_dependencies.yaml.',
            file=sys.stderr,
        )
        return 1

    print(f'OK: No hard-blocked dependencies are imported. ({approved_path})')
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(description='Check dependency license approval status')
    parser.add_argument(
        '--approved',
        type=Path,
        default=Path('specs/approved_dependencies.yaml'),
        help='Path to approved_dependencies.yaml',
    )
    parser.add_argument(
        '--src-dirs',
        nargs='+',
        type=Path,
        default=[Path('src'), Path('scripts')],
        help='Source directories to scan for imports (default: src scripts)',
    )
    args = parser.parse_args()

    if not args.approved.exists():
        print(f'ERROR: Approved dependencies file not found: {args.approved}', file=sys.stderr)
        sys.exit(1)

    sys.exit(check(args.approved, args.src_dirs))


if __name__ == '__main__':
    main()
