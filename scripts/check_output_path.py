"""Pre-run check: verify that an output directory argument is defined.

Called by .claude/hooks/pre-run.sh before any simulation script.
Checks that --output-dir is present so every run has a designated artifact location.

Usage (from hook):
    python scripts/check_output_path.py "$@"

Exit 0 if --output-dir is defined, 1 otherwise.
"""
from __future__ import annotations

import sys
from pathlib import Path


def check(argv: list[str]) -> int:
    if '--output-dir' not in argv:
        print(
            'ERROR: --output-dir is required for all simulation runs.\n'
            'All experiment outputs (trajectory, manifest, summary) must be saved '
            'to a designated directory for traceability. Pass --output-dir <path>.',
            file=sys.stderr,
        )
        return 1

    idx = argv.index('--output-dir')
    if idx + 1 >= len(argv):
        print('ERROR: --output-dir flag provided but no value given.', file=sys.stderr)
        return 1

    out_path = argv[idx + 1]
    if not out_path or out_path.startswith('-'):
        print(
            f'ERROR: --output-dir value looks invalid: {out_path!r}',
            file=sys.stderr,
        )
        return 1

    # Warn (not fail) if path already has content, so the user can decide
    p = Path(out_path)
    if p.exists() and any(p.iterdir()):
        print(
            f'WARNING: --output-dir {out_path!r} already exists and is non-empty. '
            'Existing files may be overwritten.',
            file=sys.stderr,
        )

    print(f'OK: output-dir={out_path!r}')
    return 0


def main() -> None:
    sys.exit(check(sys.argv[1:]))


if __name__ == '__main__':
    main()
