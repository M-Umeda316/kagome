"""Geometry utilities for periodic boundary conditions."""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def minimum_image(
    r_vec: NDArray[np.floating],
    cell: NDArray[np.floating] | None,
) -> NDArray[np.floating]:
    """Apply minimum image convention for an orthorhombic cell.

    r_vec: displacement vector(s), shape (3,) or (N, 3)
    cell:  (3, 3) lattice matrix (only diagonal used) or None
    Returns adjusted displacement(s).
    """
    if cell is None:
        return r_vec
    box = np.array([cell[0, 0], cell[1, 1], cell[2, 2]])
    return r_vec - box * np.round(r_vec / box)


def wrap_positions(
    positions: NDArray[np.floating],
    cell: NDArray[np.floating] | None,
) -> None:
    """Wrap atomic positions into the primary orthorhombic cell (in-place).

    positions: (N, 3) array of coordinates (Å), modified in-place
    cell:      (3, 3) lattice matrix (only diagonal used) or None (no-op)
    """
    if cell is None:
        return
    box = np.array([cell[0, 0], cell[1, 1], cell[2, 2]])
    positions[:] = positions % box
