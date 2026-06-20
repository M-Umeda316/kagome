"""MACE universal MLIP backend.

Uses MACE-MP-0 (MIT license) via ASE adapter.
Code: MIT, Model weights (MACE-MP-0): MIT.
See specs/dependency-license-matrix.md.
"""
from __future__ import annotations

from pathlib import Path

from src.backends.ase_adapter import ASECalculatorAdapter

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
_LOCAL_MODELS = {
    'small': _PROJECT_ROOT / 'models' / 'mace-mp-0-small.model',
}


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

    local_path = _LOCAL_MODELS.get(model)
    model_arg = str(local_path) if local_path and local_path.exists() else model

    calc = mace_mp(
        model=model_arg,
        device=device,
        default_dtype=default_dtype,
    )
    return ASECalculatorAdapter(calc, name=f'mace-mp-0-{model}')
