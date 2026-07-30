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
    """1 タイムステップの状態スナップショット (JSONL の 1 行に対応)。"""

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
    # (3,3) simulation cell at this frame (Å), or None for non-periodic runs.
    # Recorded per frame because NPT lets the box vary between frames, so
    # depth-resolved density (analysis/density.py) needs the cell of the exact
    # event frame. Optional with a None default for backward compatibility:
    # frames written before this field existed load as cell=None.
    cell: list[list[float]] | None = None


class TrajectoryWriter:
    """Writes trajectory frames to a JSONL file with a metadata header."""

    def __init__(
        self,
        path: Path,
        species: list[str],
        save_interval: int = 100,
        metadata: dict[str, Any] | None = None,
        n_reactive_sites: int | None = None,
        append: bool = False,
        initial_bonds: list[tuple[int, int, float]] | None = None,
    ) -> None:
        self._path = path
        self._save_interval = save_interval
        self._step_counter = 0
        path.parent.mkdir(parents=True, exist_ok=True)
        # append=True (checkpoint resume): keep prior frames and the existing
        # header instead of truncating; only fresh runs write the header.
        self._file = open(path, 'a' if append else 'w', encoding='utf-8')
        if append:
            return

        header: dict[str, Any] = {
            '_header': True,
            'schema_version': 1,
            'species': species,
            'n_atoms': len(species),
            'save_interval': save_interval,
        }
        if n_reactive_sites is not None:
            # Correct denominator for alpha(t) = N_reacted / N_reactive_sites
            header['n_reactive_sites'] = n_reactive_sites
        if initial_bonds is not None:
            # Initial intramolecular connectivity so viewers use real bonds
            # instead of distance inference (specs/decisions.md 2026-07-02).
            # Time-evolving bonds live in topology.jsonl.
            header['bonds'] = [[int(i), int(j), float(o)] for i, j, o in initial_bonds]
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

    def flush(self) -> None:
        """Flush buffered frames to the OS so a SIGKILL/OOM keeps them.

        Called at cycle boundaries (after the checkpoint save) rather than every
        frame: per-frame flushing would negate the write buffer's throughput
        benefit, while a cycle-boundary flush bounds the loss to at most the
        frames of the in-progress cycle — which the resume truncation would
        rewrite anyway.
        """
        if not self._file.closed:
            self._file.flush()

    def close(self) -> None:
        if not self._file.closed:
            self._file.close()

    @property
    def n_frames(self) -> int:
        return self._step_counter
