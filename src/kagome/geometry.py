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


# ---------------------------------------------------------------------------
# Fast-path variants for inner MD loops
# ---------------------------------------------------------------------------

def validated_box(
    cell: NDArray[np.floating] | None,
) -> NDArray[np.floating] | None:
    """Validate an orthorhombic cell once and return its (3,) box diagonal.

    Call once per phase (or after barostat acceptance), then pass the
    returned ``box`` to ``minimum_image_fast`` / ``wrap_positions_fast``.
    """
    if cell is None:
        return None
    _check_orthorhombic(cell)
    return np.array([cell[0, 0], cell[1, 1], cell[2, 2]])


def validated_box_cached(
    cell: NDArray[np.floating] | None,
    cached_cell: NDArray[np.floating] | None,
    cached_box: NDArray[np.floating] | None,
) -> tuple[NDArray[np.floating] | None, NDArray[np.floating] | None]:
    """Return the (3,) box for *cell*, validating only when its contents changed.

    Shared inner-loop helper for the Velocity-Verlet / Langevin integrators.
    ``cached_cell`` is a *copy* of the cell contents last validated (or None on
    the first call) and ``cached_box`` the (3,) diagonal returned for it. The
    orthorhombic check (:func:`validated_box`) runs only when ``cell`` differs
    from ``cached_cell`` (by value), so the steady state — an unchanged cell
    passed every step — skips it entirely.

    A content comparison (not identity) is required because the MC barostat
    mutates the cell in place (``cell[:] = new_cell``): the same array object is
    handed back each step, but its values change on an accepted move.

    Returns ``(box, cell_snapshot)`` where ``box`` is the diagonal to use (None
    when ``cell`` is None) and ``cell_snapshot`` is the copy to store for the
    next call. The returned ``box`` is numerically identical to
    ``validated_box(cell)``.
    """
    if cell is None:
        return None, None
    if (cached_cell is not None and cached_box is not None
            and np.array_equal(cell, cached_cell)):
        return cached_box, cached_cell
    return validated_box(cell), np.array(cell, copy=True)


def minimum_image_fast(
    r_vec: NDArray[np.floating],
    box: NDArray[np.floating] | None,
) -> NDArray[np.floating]:
    """Minimum image without per-call orthorhombic validation.

    ``box`` must be a prevalidated (3,) diagonal from ``validated_box``. A
    lightweight shape guard catches a (3, 3) cell passed here by mistake, which
    would otherwise broadcast silently into a wrong minimum image.
    """
    if box is None:
        return r_vec
    if box.ndim != 1:
        raise ValueError(
            f'box must be a (3,) diagonal from validated_box, got ndim={box.ndim} '
            '(a full (3,3) cell?)'
        )
    return r_vec - box * np.round(r_vec / box)


def wrap_positions_fast(
    positions: NDArray[np.floating],
    box: NDArray[np.floating] | None,
) -> None:
    """In-place wrap without per-call orthorhombic validation.

    ``box`` must be a prevalidated (3,) diagonal from ``validated_box``. A
    lightweight shape guard catches a (3, 3) cell passed here by mistake.
    """
    if box is None:
        return
    if box.ndim != 1:
        raise ValueError(
            f'box must be a (3,) diagonal from validated_box, got ndim={box.ndim} '
            '(a full (3,3) cell?)'
        )
    positions[:] = positions % box
