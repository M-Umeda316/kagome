"""Portable serialization of a prepared atomic structure.

The classical structure-prep stage (``src.prep.openmm_equilibrate``) and the ML
production run (OrbMol-v2) may execute in different conda environments — OpenMM /
OpenFF for prep, OrbMol-v2 for production — so the relaxed structure is handed
over through a small JSON file rather than in-process.

Conventions (this repo's unit system):
- ``positions`` : (N, 3) in Å
- ``cell``      : (3, 3) lattice vectors in Å, or ``None`` for non-periodic
- ``species``   : element symbols, length N, in the SAME order as ``positions``

Atom ordering is preserved verbatim so the caller's group / propagation_map
indices remain valid (see specs/decisions.md decision D-4).
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

# Format version, bumped if the on-disk schema changes.
SCHEMA_VERSION = 1


@dataclass
class PreparedStructure:
    """An equilibrated structure handed from classical prep to ML production."""

    positions: NDArray[np.floating]          # (N, 3) Å
    species: list[str]                       # length N
    cell: NDArray[np.floating] | None = None  # (3, 3) Å or None
    metadata: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.positions = np.asarray(self.positions, dtype=np.float64)
        if self.positions.ndim != 2 or self.positions.shape[1] != 3:
            raise ValueError(
                f'positions must be (N, 3); got {self.positions.shape}'
            )
        if len(self.species) != self.positions.shape[0]:
            raise ValueError(
                f'species length {len(self.species)} != n_atoms '
                f'{self.positions.shape[0]}'
            )
        if self.cell is not None:
            self.cell = np.asarray(self.cell, dtype=np.float64)
            if self.cell.shape != (3, 3):
                raise ValueError(f'cell must be (3, 3); got {self.cell.shape}')

    @property
    def n_atoms(self) -> int:
        return self.positions.shape[0]

    def save(self, path: str | Path) -> None:
        payload = {
            'schema_version': SCHEMA_VERSION,
            'units': {'length': 'angstrom'},
            'species': list(self.species),
            'positions_A': self.positions.tolist(),
            'cell_A': None if self.cell is None else self.cell.tolist(),
            'metadata': self.metadata,
        }
        Path(path).write_text(json.dumps(payload), encoding='utf-8')

    @classmethod
    def load(cls, path: str | Path) -> 'PreparedStructure':
        data = json.loads(Path(path).read_text(encoding='utf-8'))
        version = data.get('schema_version')
        if version != SCHEMA_VERSION:
            raise ValueError(
                f'Unsupported structure schema_version {version!r} '
                f'(expected {SCHEMA_VERSION}) in {path}'
            )
        cell = data.get('cell_A')
        return cls(
            positions=np.array(data['positions_A'], dtype=np.float64),
            species=list(data['species']),
            cell=None if cell is None else np.array(cell, dtype=np.float64),
            metadata=data.get('metadata', {}),
        )
