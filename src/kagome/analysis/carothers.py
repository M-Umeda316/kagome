"""Carothers equation for step-growth polymerization.

DPn = 1 / (1 - p)

where p is the extent of reaction (conversion) and DPn is the
number-average degree of polymerization.

Paper anchor: arXiv:2511.22874, Fig. 4c — Carothers theoretical curve
compared to TDBB-simulated DPn vs conversion.
"""
from __future__ import annotations

import logging

import numpy as np
from numpy.typing import NDArray

logger = logging.getLogger(__name__)


def dpn_carothers(p: NDArray[np.floating] | float) -> NDArray[np.floating]:
    """Carothers equation: DPn = 1 / (1 - p).

    ``p`` (extent of reaction) is clamped to 0.9999 so that p >= 1.0 yields a
    finite DPn (10000) instead of ``inf``/negative, matching ``dpn_from_bonds``.
    A clamp fires a warning because p >= 1.0 signals full conversion or a
    bond-count/denominator miscount (F9, specs/fix-plan-2026-07-06).
    """
    p_arr = np.asarray(p, dtype=np.float64)
    if np.any(p_arr > 0.9999):
        n_clamped = int(np.count_nonzero(p_arr > 0.9999))
        p_max = float(np.max(p_arr))
        logger.warning(
            'dpn_carothers: %d extent-of-reaction value(s) exceed the 0.9999 '
            'clamp (max p=%.6f); DPn capped at 10000. p>=1 would indicate full '
            'conversion or a bond-count/denominator error.',
            n_clamped, p_max,
        )
    p_clamped = np.minimum(p_arr, 0.9999)
    return 1.0 / (1.0 - p_clamped)


def dpn_from_bonds(
    n_bonds: int,
    n_functional_groups: int,
) -> float:
    """Compute DPn from bond count and initial functional group count.

    Equimolar A-A + B-B step-growth only (e.g., nylon-6,6, N_A = N_B):
      p = n_bonds / (n_functional_groups / 2)
      DPn = 1 / (1 - p)

    Stoichiometric imbalance (r != 1) and monofunctional end-capping are NOT
    modeled here. ``n_functional_groups`` is the total number of reactive end
    groups (A + B); for equimolar systems each type contributes
    ``n_functional_groups / 2``.

    Counting convention (A5): ``n_bonds`` is the number of amide bonds, i.e.
    ``confirmed_formation`` events with ``counts_as_reaction=True``. Bias-only
    water-forming k-l events must be excluded by the caller so one condensation
    counts once (specs/decisions.md 2026-07-06).
    """
    if n_functional_groups <= 0:
        return 1.0
    raw_p = n_bonds / (n_functional_groups / 2)
    if raw_p > 0.9999:
        # Surface the clamp instead of silently masking near-complete conversion
        # or a p>1 bond/denominator miscount (RF19b).
        logger.warning(
            'dpn_from_bonds: extent of reaction p=%.6f exceeds the 0.9999 clamp '
            '(n_bonds=%d, n_functional_groups=%d); DPn capped at 10000. p>1 would '
            'indicate a bond-count or denominator error.',
            raw_p, n_bonds, n_functional_groups,
        )
    p = min(raw_p, 0.9999)
    return 1.0 / (1.0 - p)
