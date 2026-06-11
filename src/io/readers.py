"""Trajectory reader: loads JSONL trajectory files."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from src.io.trajectory import TrajectoryFrame


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
