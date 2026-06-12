"""JSONL trajectory writer for simulation output.

Each line is a JSON object: first line is a header with metadata,
subsequent lines are TrajectoryFrame records.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


@dataclass
class TrajectoryFrame:
    step: int
    time_fs: float
    phase: str
    cycle: int
    energy_base: float
    energy_bias: float
    energy_total: float
    positions: list[list[float]]
    n_candidates: int = 0
    n_selected: int = 0
    temperature_K: float = 0.0


class TrajectoryWriter:
    """Writes trajectory frames to a JSONL file with a metadata header."""

    def __init__(
        self,
        path: Path,
        species: list[str],
        save_interval: int = 100,
        metadata: dict[str, Any] | None = None,
    ) -> None:
        self._path = path
        self._save_interval = save_interval
        self._step_counter = 0
        path.parent.mkdir(parents=True, exist_ok=True)
        self._file = open(path, 'w', encoding='utf-8')

        header = {
            '_header': True,
            'species': species,
            'n_atoms': len(species),
            'save_interval': save_interval,
        }
        if metadata:
            header['metadata'] = metadata
        self._file.write(json.dumps(header) + '\n')

    def should_write(self, step_in_phase: int) -> bool:
        if self._save_interval <= 0:
            return False
        return step_in_phase % self._save_interval == 0

    def write_frame(self, frame: TrajectoryFrame) -> None:
        self._file.write(json.dumps(asdict(frame)) + '\n')
        self._step_counter += 1

    def close(self) -> None:
        if not self._file.closed:
            self._file.close()

    @property
    def n_frames(self) -> int:
        return self._step_counter
