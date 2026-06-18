"""Depth-resolved reaction density.

Paper: arXiv:2511.22874, PDF p.12 (unnumbered equation).
ρ_rxn(z) = N_rxn(z) / (A · Δz · N_frames)
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from src.reactive.bonds import BondEvent


def reaction_density_profile(
    events: list[BondEvent],
    positions_at_event: dict[int, NDArray[np.floating]],
    z_bins: NDArray[np.floating],
    area_xy: float,
    n_frames: int = 1,
) -> NDArray[np.floating]:
    """Depth-resolved reaction density (PDF p.12, unnumbered).

    events: confirmed bond events
    positions_at_event: {step: positions_array} for each event step
    z_bins: bin edges along z-axis (Å)
    area_xy: cross-sectional area in Å²
    n_frames: number of frames averaged over

    Returns density array of shape (len(z_bins)-1,).
    """
    dz = np.diff(z_bins)
    counts = np.zeros(len(z_bins) - 1, dtype=np.float64)

    for ev in events:
        if ev.step not in positions_at_event:
            continue
        pos = positions_at_event[ev.step]
        z_mid = 0.5 * (pos[ev.atom_a, 2] + pos[ev.atom_b, 2])
        bin_idx = np.searchsorted(z_bins, z_mid, side='right') - 1
        if 0 <= bin_idx < len(counts):
            counts[bin_idx] += 1

    density = counts / (area_xy * dz * n_frames)
    return density
