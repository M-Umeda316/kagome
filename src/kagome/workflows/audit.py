"""JSONL audit logging for the polymerization workflow.

The workflow emits four line-per-record audit streams alongside the trajectory:

* ``selection.jsonl`` — per-cycle candidate ranking / selection decisions
  (and the valence-drop sub-records);
* ``topology.jsonl`` — connectivity snapshots after bond-changing cycles;
* ``mixing.jsonl`` — one record per WM-P3 classical mixing segment.

Each was written ad hoc with its own ``open(path, 'a')`` /
``json.dumps(record) + '\\n'`` block and its own truncate-on-init. This module
collects that single pattern behind :class:`JsonlAuditLog` so there is one place
that knows the append and truncate lifecycle. The emitted bytes are identical to
the previous inline writes: ``json.dumps`` with default separators, one record
per line, a trailing newline, UTF-8, append mode.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


class JsonlAuditLog:
    """Append-only JSONL audit stream at a fixed path.

    Thin by design: it opens the file per :meth:`append` (matching the prior
    inline behaviour, which never held the file open) so that a crash between
    cycles cannot lose a buffered audit record.
    """

    def __init__(self, path: Path) -> None:
        self._path = Path(path)

    @property
    def path(self) -> Path:
        """The backing file path (needed by the resume truncation helper)."""
        return self._path

    def truncate(self) -> None:
        """Start a fresh, empty log — discards any prior run's records.

        Mirrors the ``path.write_text('', encoding='utf-8')`` used at run start
        for a non-resuming run.
        """
        self._path.parent.mkdir(parents=True, exist_ok=True)
        self._path.write_text('', encoding='utf-8')

    def append(self, record: dict[str, Any]) -> None:
        """Append one JSON record as a single line (byte-identical to before)."""
        self._path.parent.mkdir(parents=True, exist_ok=True)
        with open(self._path, 'a', encoding='utf-8') as f:
            f.write(json.dumps(record) + '\n')
