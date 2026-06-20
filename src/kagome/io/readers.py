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
