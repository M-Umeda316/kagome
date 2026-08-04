"""Unit tests for scripts/run_vinyl_copolymer.py CLI guards.

Pure-Python (no orb/GPU/openmm): the covered guard fires in argparse before any
system build or backend creation.

Covers the symmetric CLI guard as seen from the vinyl driver: any mixing knob
given WITHOUT --mix is a usage error (parser.error -> SystemExit).

The mixing helpers themselves moved to ``scripts/_mixing_cli.py``, shared with
run_nylon66 / run_epoxy_amine (decisions.md 2026-08-04); the resume mode-switch
guard, the knob-default resolution and the same CLI guard on the two newly
wired drivers are covered in tests/unit/test_mixing_cli.py.
"""
from __future__ import annotations

import sys

import pytest

from scripts.run_vinyl_copolymer import main


@pytest.mark.parametrize('flag', [
    ['--mix-ps', '1.0'],                 # pre-existing behavior, pinned
    ['--mix-settle-steps', '100'],       # pre-existing behavior, pinned
    ['--mix-friction-per-ps', '1.0'],    # newly guarded
    ['--mix-timestep-fs', '0.5'],        # newly guarded
    ['--mix-platform', 'CPU'],           # newly guarded
    ['--mix-charge-method', 'nagl'],     # newly guarded
])
def test_mix_knob_without_mix_errors(monkeypatch, tmp_path, flag):
    """Any mixing knob supplied WITHOUT --mix is a usage error. The guard fires
    right after argparse, before any system build or backend creation."""
    monkeypatch.setattr(sys, 'argv', [
        'run_vinyl_copolymer.py',
        '--output-dir', str(tmp_path),
        '--backend', 'toy',
        *flag,
    ])
    with pytest.raises(SystemExit):
        main()
