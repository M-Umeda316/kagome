"""Conversion tracking from bond events.

Paper: arXiv:2511.22874, Eq. 11-12.
α(t) = N_reacted(t) / N_total_reactive_sites
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from src.reactive.bonds import BondEvent


def conversion(n_reacted: int, n_total: int) -> float:
    """Eq. 11: α = N_reacted / N_total."""
    if n_total == 0:
        return 0.0
    return n_reacted / n_total


def conversion_timeseries(
    events: list[BondEvent],
    n_total_sites: int,
    step_range: NDArray[np.integer] | None = None,
) -> tuple[NDArray[np.integer], NDArray[np.floating]]:
    """Build α(step) from confirmed formation events.

    Returns (steps, alpha) arrays suitable for plotting.
    """
    formations = [e for e in events if e.event_type == 'confirmed_formation']
    if not formations:
        if step_range is not None:
            return step_range, np.zeros(len(step_range), dtype=np.float64)
        return np.array([0], dtype=np.int64), np.array([0.0], dtype=np.float64)

    formations.sort(key=lambda e: e.step)

    if step_range is None:
        max_step = formations[-1].step
        step_range = np.arange(0, max_step + 1, dtype=np.int64)

    alpha = np.zeros(len(step_range), dtype=np.float64)
    cumulative = 0
    event_idx = 0

    for i, s in enumerate(step_range):
        while event_idx < len(formations) and formations[event_idx].step <= s:
            cumulative += 1
            event_idx += 1
        alpha[i] = conversion(cumulative, n_total_sites)

    return step_range, alpha
