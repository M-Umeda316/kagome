"""Pre-run check: verify that a seed argument is defined.

Called by .claude/hooks/pre-run.sh before any simulation script.
Reads CLI args forwarded from the hook and checks that --seed is present
and is a valid non-negative integer.

Usage (from hook):
    python scripts/check_seed_defined.py "$@"

Exit 0 if seed is defined, 1 otherwise.
"""
from __future__ import annotations

import sys


def check(argv: list[str]) -> int:
    if '--seed' not in argv:
        print(
            'ERROR: --seed is required for all simulation runs.\n'
            'All experiments must be reproducible. Pass --seed <int> to your script.',
            file=sys.stderr,
        )
        return 1

    idx = argv.index('--seed')
    if idx + 1 >= len(argv):
        print('ERROR: --seed flag provided but no value given.', file=sys.stderr)
        return 1

    seed_val = argv[idx + 1]
    try:
        seed_int = int(seed_val)
        if seed_int < 0:
            raise ValueError('negative')
    except ValueError:
        print(
            f'ERROR: --seed value must be a non-negative integer, got: {seed_val!r}',
            file=sys.stderr,
        )
        return 1

    print(f'OK: seed={seed_int}')
    return 0


def main() -> None:
    # argv[0] is this script; remaining args are the forwarded CLI args
    sys.exit(check(sys.argv[1:]))


if __name__ == '__main__':
    main()
