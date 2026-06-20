"""Geometry utilities for periodic boundary conditions.

Orthorhombic cells only — only the diagonal elements of the (3,3) lattice
matrix are used.  Triclinic (non-zero off-diagonal) cells will raise
``ValueError``.  See specs/decisions.md "2026-06-19: RF12".
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def _check_orthorhombic(cell: NDArray[np.floating]) -> None:
    """Raise ValueError if *cell* has non-negligible off-diagonal elements."""
    off_diag = cell.copy()
    np.fill_diagonal(off_diag, 0.0)
    if np.any(np.abs(off_diag) > 1e-10):
        raise ValueError(
            'Only orthorhombic (diagonal) cells are supported. '
            f'Off-diagonal elements: {off_diag.ravel().tolist()}'
        )


def minimum_image(
    r_vec: NDArray[np.floating],
    cell: NDArray[np.floating] | None,
) -> NDArray[np.floating]:
    """Apply minimum image convention for an orthorhombic cell.

    Orthorhombic only — raises ``ValueError`` for triclinic cells.

    r_vec: displacement vector(s), shape (3,) or (N, 3)
    cell:  (3, 3) diagonal lattice matrix, or None (non-periodic)
    Returns adjusted displacement(s).
    """
    if cell is None:
        return r_vec
    _check_orthorhombic(cell)
    box = np.array([cell[0, 0], cell[1, 1], cell[2, 2]])
    return r_vec - box * np.round(r_vec / box)


def wrap_positions(
    positions: NDArray[np.floating],
    cell: NDArray[np.floating] | None,
) -> None:
    """Wrap atomic positions into the primary orthorhombic cell (in-place).

    Orthorhombic only — raises ``ValueError`` for triclinic cells.

    positions: (N, 3) array of coordinates (Å), modified in-place
    cell:      (3, 3) diagonal lattice matrix, or None (no-op)
    """
    if cell is None:
        return
    _check_orthorhombic(cell)
    box = np.array([cell[0, 0], cell[1, 1], cell[2, 2]])
    positions[:] = positions % box
