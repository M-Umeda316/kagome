"""Trajectory and bond-event readers."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from kagome.io.trajectory import TrajectoryFrame
from kagome.reactive.bonds import BondEvent


def read_trajectory(path: Path) -> tuple[dict[str, Any], list[TrajectoryFrame]]:
    """Read a JSONL trajectory file.

    Returns (header_dict, list_of_frames).
    """
    header: dict[str, Any] = {}
    frames: list[TrajectoryFrame] = []

    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            record = json.loads(line)
            if record.get('_header'):
                header = record
            else:
                frames.append(TrajectoryFrame(**record))

    return header, frames


def read_topology_snapshots(
    path: Path,
) -> list[tuple[int, list[tuple[int, int, float]]]]:
    """Read a topology.jsonl written by PolymerizationWorkflow.

    Returns ``[(step, bonds), ...]`` sorted by step, where each ``bonds`` is a
    list of ``(i, j, order)``.  Empty list if the file does not exist.
    """
    if not path.exists():
        return []
    snapshots: list[tuple[int, list[tuple[int, int, float]]]] = []
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            bonds = [(int(i), int(j), float(o)) for i, j, o in rec['bonds']]
            snapshots.append((int(rec['step']), bonds))
    snapshots.sort(key=lambda s: s[0])
    return snapshots


def bonds_at_step(
    snapshots: list[tuple[int, list[tuple[int, int, float]]]],
    step: int,
) -> list[tuple[int, int, float]]:
    """Return the bond list in effect at *step* (latest snapshot with step<=).

    Snapshots are recorded at the initial state (step of the first frame) and
    after each cycle that formed a bond, so the connectivity for any frame is
    the most recent snapshot at or before that frame's step.
    """
    current: list[tuple[int, int, float]] = []
    for snap_step, bonds in snapshots:  # sorted ascending
        if snap_step <= step:
            current = bonds
        else:
            break
    return current


def read_bond_events(path: Path) -> list[BondEvent]:
    """Read a bonds.jsonl file written by BondTracker.save().

    Returns a list of BondEvent records.  Returns an empty list if the file
    does not exist (no bond events were recorded).
    """
    if not path.exists():
        return []

    events: list[BondEvent] = []
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            record = json.loads(line)
            events.append(BondEvent(**record))
    return events
