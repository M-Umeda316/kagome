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
) -> NDArray[np.floating]:
    """Depth-resolved reaction density (PDF p.12, unnumbered).

    events: confirmed bond events
    positions_at_event: {step: positions_array} for each event step
    z_bins: bin edges along z-axis (Å)
    area_xy: cross-sectional area in Å²
    n_frames: number of frames averaged over
    cell: (3,3) cell matrix for PBC midpoint correction (None = no PBC)

    Returns density array of shape (len(z_bins)-1,).
    """
    dz = np.diff(z_bins)
    counts = np.zeros(len(z_bins) - 1, dtype=np.float64)

    for ev in events:
        if ev.step not in positions_at_event:
            continue
        pos = positions_at_event[ev.step]
        z_mid = _pbc_midpoint_z(pos[ev.atom_a], pos[ev.atom_b], cell)
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
