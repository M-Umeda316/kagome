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

    # Guard against MACE-OFF (ASL, blocked_pending_review): this factory only
    # serves MACE-MP-0. Reject any MACE-OFF model string/path so restricted
    # weights cannot be loaded through here. See specs/approved_dependencies.yaml.
    if any(tok in str(model).lower() for tok in ('mace_off', 'mace-off', 'off23')):
        raise RuntimeError(
            f'MACE-OFF weights ({model!r}) are blocked_pending_review (ASL restricts '
            'commercial use); this backend only provides MACE-MP-0. See '
            'specs/approved_dependencies.yaml.'
        )

    local_path = _LOCAL_MODELS.get(model)
    model_arg = str(local_path) if local_path and local_path.exists() else model

    calc = mace_mp(
        model=model_arg,
        device=device,
        default_dtype=default_dtype,
    )
    # Record the resolved weights (local checkpoint path or model size) for
    # provenance (RF17), distinct from the human-readable backend name.
    return ASECalculatorAdapter(
        calc, name=f'mace-mp-0-{model}', model_id=f'mace-mp-0:{model_arg}',
    )
