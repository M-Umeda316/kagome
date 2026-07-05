"""Conversion tracking from bond events.

Paper: arXiv:2511.22874.
- α = 1 − [M]/[M]₀ (PDF p.9 Fig.2 caption). Denominator = initial monomer count.
- Eq. 11 (PDF p.9): α(t) = 1 - exp(-k*_p · t)  [exponential fit]
"""
from __future__ import annotations

import logging

import numpy as np
from numpy.typing import NDArray

from kagome.reactive.bonds import BondEvent
from kagome.reactive.groups import ReactiveGroup

logger = logging.getLogger(__name__)


def monomer_site_count(groups: dict[str, ReactiveGroup], monomer_group: str = 'vinyl_alpha_C') -> int:
    """Return the initial monomer count for α denominator.

    For vinyl systems the monomer count equals the initial size of the
    vinyl_alpha_C group (one per monomer molecule).  Constraint-only
    groups (chain_C, vinyl_beta_C) and radical_C are excluded.

    Paper anchor: PDF p.9 Fig.2 caption, α = 1 − [M]/[M]₀ where
    [M]₀ = initial monomer concentration → denominator = n_monomers.
    """
    g = groups.get(monomer_group)
    if g is None:
        return 0
    return len(g.atom_indices)


def conversion(n_reacted: int, n_total: int) -> float:
    """Raw conversion α = N_reacted / N_total.

    Denominator should be the initial monomer count (PDF p.9 Fig.2).
    """
    if n_total == 0:
        return 0.0
    return n_reacted / n_total


def conversion_timeseries(
    events: list[BondEvent],
    n_total_sites: int,
    step_range: NDArray[np.integer] | None = None,
) -> tuple[NDArray[np.integer], NDArray[np.floating]]:
    """Build α(step) from confirmed formation events.

    n_total_sites should be the initial monomer count (PDF p.9 Fig.2,
    α = 1 − [M]/[M]₀).  Returns (steps, alpha) arrays.
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
    production_start_step: int = 0,
) -> tuple[float, float]:
    """Fit α(t) = 1 - exp(-kp_eff * t) and return (kp_eff, r_squared).

    Eq. 11 (PDF p.9): α(t) = 1 - exp(-k*_p · t)

    Args:
        steps:       step indices (integer), shape (N,)
        alpha:       conversion values in [0, 1], shape (N,)
        timestep_fs: MD timestep in fs (converts steps -> physical time)
        production_start_step: global step at which production begins;
            subtracted from *steps* before converting to physical time
            so the fit starts at t=0 (L8).

    Returns:
        kp_eff:    effective polymerization rate constant (1/fs)
        r_squared: coefficient of determination for the fit quality
    """
    try:
        from scipy.optimize import curve_fit
    except ImportError as e:
        raise ImportError(
            'scipy is required for exponential fitting. '
            'Install with: pip install kagome[fit]'
        ) from e

    t = (steps.astype(np.float64) - production_start_step) * timestep_fs

    # Guard: need at least some non-zero alpha and non-trivial data
    if alpha.max() < 1e-9 or len(t) < 3:
        logger.warning(
            'fit_conversion_exponential: skipping fit (alpha.max=%.3g, n_points=%d); '
            'need alpha>0 and >=3 points. Returning (kp=0, R2=0).',
            float(alpha.max()) if alpha.size else 0.0, len(t),
        )
        return 0.0, 0.0

    # Clip alpha to avoid log(0) in curve_fit internals
    alpha_clipped = np.clip(alpha, 0.0, 1.0 - 1e-9)

    def model(t_arr: NDArray, kp: float) -> NDArray:
        return 1.0 - np.exp(-kp * t_arr)

    try:
        popt, _ = curve_fit(model, t, alpha_clipped, p0=[1e-4], bounds=(0, np.inf), maxfev=5000)
        kp_eff = float(popt[0])
    except RuntimeError as e:
        logger.warning(
            'fit_conversion_exponential: curve_fit did not converge (%s); '
            'returning (kp=0, R2=0).', e,
        )
        return 0.0, 0.0

    residuals = alpha_clipped - model(t, kp_eff)
    ss_res = float(np.sum(residuals ** 2))
    ss_tot = float(np.sum((alpha_clipped - alpha_clipped.mean()) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else 0.0

    return kp_eff, r_squared
