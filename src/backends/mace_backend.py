"""MACE universal MLIP backend.

Uses MACE-MP-0 (MIT license) via ASE adapter.
Code: MIT, Model weights (MACE-MP-0): MIT.
See specs/dependency-license-matrix.md.
"""
from __future__ import annotations

from typing import Literal

from src.backends.ase_adapter import ASECalculatorAdapter


def create_mace_calculator(
    model: str = 'small',
    device: str = 'cpu',
    default_dtype: str = 'float64',
) -> ASECalculatorAdapter:
    """Create a MACE-MP-0 backed calculator.

    model: 'small', 'medium', or 'large' — MACE-MP-0 model size.
    device: 'cpu' or 'cuda'.
    default_dtype: 'float32' or 'float64'.
    """
    try:
        from mace.calculators import mace_mp
    except ImportError:
        raise ImportError(
            'mace-torch is required for this backend. '
            'Install with: pip install mace-torch'
        )

    calc = mace_mp(
        model=model,
        device=device,
        default_dtype=default_dtype,
    )
    return ASECalculatorAdapter(calc, name=f'mace-mp-0-{model}')
