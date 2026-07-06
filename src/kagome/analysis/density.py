"""Depth-resolved reaction density.

Paper: arXiv:2511.22874, PDF p.12 (unnumbered equation).
ρ_rxn(z) = N_rxn(z) / (A · Δz · N_frames)
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from kagome.geometry import minimum_image, wrap_positions
from kagome.reactive.bonds import BondEvent


def reaction_density_profile(
    events: list[BondEvent],
    positions_at_event: dict[int, NDArray[np.floating]],
    z_bins: NDArray[np.floating],
    area_xy: float,
    n_frames: int = 1,
    cell: NDArray[np.floating] | None = None,
    cells_at_event: dict[int, NDArray[np.floating]] | None = None,
) -> NDArray[np.floating]:
    """Depth-resolved reaction density (PDF p.12, unnumbered).

    events: confirmed bond events
    positions_at_event: {step: positions_array} for each event step
    z_bins: bin edges along z-axis (Å)
    area_xy: cross-sectional area in Å²
    n_frames: number of frames averaged over (all sampled trajectory frames in
        the analysis window, per the paper's N_frames — NOT the number of events)
    cell: (3,3) cell matrix for PBC midpoint correction (None = no PBC).
        Used as the fallback when ``cells_at_event`` has no entry for an event.
    cells_at_event: {step: (3,3) cell} giving the box at each event's frame.
        NPT lets the box vary, so the midpoint correction must use the event's
        own cell; falls back to ``cell`` when a step is absent.

    Returns density array of shape (len(z_bins)-1,).

    Note (F7): unlike ``np.histogram``, a midpoint landing exactly on the final
    bin's right edge is excluded (the ``bin_idx < len(counts)`` guard). This is
    unreachable in practice because midpoints are wrapped into [0, Lz) and never
    equal Lz exactly, so it needs no special-casing.
    """
    dz = np.diff(z_bins)
    counts = np.zeros(len(z_bins) - 1, dtype=np.float64)

    for ev in events:
        if ev.step not in positions_at_event:
            continue
        pos = positions_at_event[ev.step]
        ev_cell = cell
        if cells_at_event is not None and ev.step in cells_at_event:
            ev_cell = cells_at_event[ev.step]
        z_mid = _pbc_midpoint_z(pos[ev.atom_a], pos[ev.atom_b], ev_cell)
        bin_idx = np.searchsorted(z_bins, z_mid, side='right') - 1
        if 0 <= bin_idx < len(counts):
            counts[bin_idx] += 1

    density = counts / (area_xy * dz * n_frames)
    return density


def _pbc_midpoint_z(
    pos_a: NDArray[np.floating],
    pos_b: NDArray[np.floating],
    cell: NDArray[np.floating] | None,
) -> float:
    """Z-coordinate of the midpoint between two atoms, PBC-aware."""
    if cell is None:
        return float(0.5 * (pos_a[2] + pos_b[2]))
    delta = minimum_image(pos_b - pos_a, cell)
    mid = pos_a + 0.5 * delta
    mid_wrapped = mid.copy().reshape(1, 3)
    wrap_positions(mid_wrapped, cell)
    return float(mid_wrapped[0, 2])
