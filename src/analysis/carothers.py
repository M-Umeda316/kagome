"""Carothers equation for step-growth polymerization.

DPn = 1 / (1 - p)

where p is the extent of reaction (conversion) and DPn is the
number-average degree of polymerization.

Paper anchor: arXiv:2511.22874, Fig. 4c — Carothers theoretical curve
compared to TDBB-simulated DPn vs conversion.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def dpn_carothers(p: NDArray[np.floating] | float) -> NDArray[np.floating]:
    """Carothers equation: DPn = 1 / (1 - p)."""
    p_arr = np.asarray(p, dtype=np.float64)
    return 1.0 / (1.0 - p_arr)


def dpn_from_bonds(
    n_bonds: int,
    n_functional_groups: int,
) -> float:
    """Compute DPn from bond count and initial functional group count.

    For equimolar A-A + B-B step-growth (e.g., nylon-6,6):
      p = n_bonds / (n_functional_groups / 2)
      DPn = 1 / (1 - p)

    n_functional_groups: total reactive sites (e.g., 2*n_diamines*2 for NH2 ends,
    but only one type limits; for equimolar, N_A = N_B = n_functional_groups / 2).
    """
    if n_functional_groups <= 0:
        return 1.0
    p = min(n_bonds / (n_functional_groups / 2), 0.9999)
    return 1.0 / (1.0 - p)
