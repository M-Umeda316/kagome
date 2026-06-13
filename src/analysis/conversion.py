"""Conversion tracking from bond events.

Paper: arXiv:2511.22874.
- Raw conversion: α = N_reacted / N_total_reactive_sites
- Eq. 11 (arXiv HTML numbering): α(t) = 1 - exp(-kp_eff * t)  [exponential fit]
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


def fit_conversion_exponential(
    steps: NDArray[np.integer],
    alpha: NDArray[np.floating],
    timestep_fs: float = 0.25,
) -> tuple[float, float]:
    """Fit α(t) = 1 - exp(-kp_eff * t) and return (kp_eff, r_squared).

    Eq. 11 (arXiv HTML): α(t) = 1 - exp(-kp_eff * t)

    Args:
        steps:       step indices (integer), shape (N,)
        alpha:       conversion values in [0, 1], shape (N,)
        timestep_fs: MD timestep in fs (converts steps -> physical time)

    Returns:
        kp_eff:    effective polymerization rate constant (1/fs)
        r_squared: coefficient of determination for the fit quality
    """
    try:
        from scipy.optimize import curve_fit
    except ImportError as e:
        raise ImportError(
            'scipy is required for exponential fitting. '
            'Install with: pip install pfpoly[fit]'
        ) from e

    t = steps.astype(np.float64) * timestep_fs  # convert to physical time (fs)

    # Guard: need at least some non-zero alpha and non-trivial data
    if alpha.max() < 1e-9 or len(t) < 3:
        return 0.0, 0.0

    # Clip alpha to avoid log(0) in curve_fit internals
    alpha_clipped = np.clip(alpha, 0.0, 1.0 - 1e-9)

    def model(t_arr: NDArray, kp: float) -> NDArray:
        return 1.0 - np.exp(-kp * t_arr)

    try:
        popt, _ = curve_fit(model, t, alpha_clipped, p0=[1e-4], bounds=(0, np.inf), maxfev=5000)
        kp_eff = float(popt[0])
    except RuntimeError:
        return 0.0, 0.0

    residuals = alpha_clipped - model(t, kp_eff)
    ss_res = float(np.sum(residuals ** 2))
    ss_tot = float(np.sum((alpha_clipped - alpha_clipped.mean()) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else 0.0

    return kp_eff, r_squared
