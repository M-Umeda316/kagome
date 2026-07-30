"""Tests for src/io/readers.py."""
import json
from pathlib import Path

import pytest

from kagome.io.readers import read_bond_events


class TestReadBondEvents:

    def test_reads_confirmed_formation(self, tmp_path):
        path = tmp_path / 'bonds.jsonl'
        event = {
            'step': 100, 'cycle': 0,
            'atom_a': 0, 'atom_b': 6,
            'event_type': 'confirmed_formation',
            'distance': 1.54, 'r0': 2.0,
        }
        path.write_text(json.dumps(event) + '\n', encoding='utf-8')

        events = read_bond_events(path)
        assert len(events) == 1
        assert events[0].event_type == 'confirmed_formation'
        assert events[0].atom_a == 0
        assert events[0].atom_b == 6
        assert events[0].distance == pytest.approx(1.54)

    def test_old_format_defaults_counts_as_reaction_true(self, tmp_path):
        """A5: bonds.jsonl written before the counts_as_reaction field must load
        with counts_as_reaction=True (backward compatible)."""
        path = tmp_path / 'bonds.jsonl'
        event = {
            'step': 100, 'cycle': 0, 'atom_a': 0, 'atom_b': 6,
            'event_type': 'confirmed_formation', 'distance': 1.54, 'r0': 2.0,
        }
        path.write_text(json.dumps(event) + '\n', encoding='utf-8')
        events = read_bond_events(path)
        assert events[0].counts_as_reaction is True

    def test_reads_counts_as_reaction_false(self, tmp_path):
        path = tmp_path / 'bonds.jsonl'
        event = {
            'step': 100, 'cycle': 0, 'atom_a': 2, 'atom_b': 3,
            'event_type': 'confirmed_formation', 'distance': 0.9, 'r0': 1.0,
            'counts_as_reaction': False,
        }
        path.write_text(json.dumps(event) + '\n', encoding='utf-8')
        events = read_bond_events(path)
        assert events[0].counts_as_reaction is False

    def test_reads_multiple_events(self, tmp_path):
        path = tmp_path / 'bonds.jsonl'
        lines = [
            json.dumps({'step': 50, 'cycle': 0, 'atom_a': 0, 'atom_b': 6,
                        'event_type': 'attempted_formation', 'distance': 3.2, 'r0': 2.0}),
            json.dumps({'step': 150, 'cycle': 0, 'atom_a': 0, 'atom_b': 6,
                        'event_type': 'confirmed_formation', 'distance': 1.8, 'r0': 2.0}),
        ]
        path.write_text('\n'.join(lines) + '\n', encoding='utf-8')

        events = read_bond_events(path)
        assert len(events) == 2
        assert events[0].event_type == 'attempted_formation'
        assert events[1].event_type == 'confirmed_formation'

    def test_nonexistent_file_returns_empty(self, tmp_path):
        events = read_bond_events(tmp_path / 'nonexistent.jsonl')
        assert events == []

    def test_empty_file_returns_empty(self, tmp_path):
        path = tmp_path / 'bonds.jsonl'
        path.write_text('', encoding='utf-8')
        events = read_bond_events(path)
        assert events == []

    def test_roundtrip_with_bond_tracker(self, tmp_path):
        import numpy as np
        from kagome.boost.tdbb import PairBias
        from kagome.reactive.bonds import BondTracker
        from kagome.reactive.pairs import TrackedPair

        tracker = BondTracker(threshold_fraction=1.3)
        pairs = [TrackedPair(PairBias(idx_a=0, idx_b=1, is_formation=True, r0=2.0))]
        positions = np.array([[0.0, 0.0, 0.0], [3.5, 0.0, 0.0]])

        tracker.record_attempts(pairs, positions, step=0, cycle=0)
        tracker.check_outcomes(positions, step=100)

        bonds_path = tmp_path / 'bonds.jsonl'
        tracker.save(bonds_path)

        events = read_bond_events(bonds_path)
        assert len(events) >= 1
        assert all(hasattr(e, 'event_type') for e in events)
