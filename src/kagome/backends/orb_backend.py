"""OrbMol-v2 backend.

Code: Apache-2.0, Model weights: Apache-2.0.
See specs/dependency-license-matrix.md.

OrbMol-v2 uses a conservative regressor: forces = -dE/dr (grad_forces).
charge and spin must be set on ASE Atoms via atoms.info.

Windows note: nvalchemiops PME triggers torch._inductor C++ compilation
which requires cl.exe.  TORCHDYNAMO_DISABLE=1 bypasses this.
"""
from __future__ import annotations

import logging
import os
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from kagome.backends.base import Calculator
from kagome.units import EV_TO_KCAL_MOL

logger = logging.getLogger(__name__)

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
_LOCAL_MODELS = {
    'orbmol_v2': _PROJECT_ROOT / 'models' / 'orbmol-v2-teqabfhg-20260523.ckpt',
}


def create_orb_calculator(
    model: str = 'orbmol_v2',
    device: str = 'cpu',
    compile: bool = False,
    charge: int = 0,
    spin: int = 1,
) -> 'OrbCalculatorAdapter':
    """Create an OrbMol-v2 backed calculator.

    model: 'orbmol_v2' (conservative, forces as energy gradients).
    device: 'cpu' or 'cuda'.
    compile: whether to use torch.compile (requires C++ compiler on Windows).
    charge: total molecular charge.
    spin: spin multiplicity (2S+1).
    """
    os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')
    os.environ.setdefault('TORCHDYNAMO_DISABLE', '1')
    os.environ.setdefault('PYTORCH_CUDA_ALLOC_CONF', 'expandable_segments:True')

    try:
        from orb_models.forcefield import pretrained
    except ImportError:
        raise ImportError(
            'orb-models is required for this backend. '
            'Install with: pip install orb-models'
        )

    loader = getattr(pretrained, model, None)
    if loader is None:
        raise ValueError(f'Unknown orb model: {model!r}')

    local_path = _LOCAL_MODELS.get(model)
    if local_path and local_path.exists():
        resolved_model_id = str(local_path)
        orbff, adapter = loader(
            weights_path=str(local_path), device=device, compile=compile,
        )
    else:
        # No local checkpoint: the loader downloads/uses its bundled pretrained
        # weights. Record that this happened so two runs with identical backend
        # names but different resolved weights are distinguishable (RF17).
        resolved_model_id = f'{model}:pretrained-download'
        logger.warning(
            'OrbMol local checkpoint not found for %r; falling back to the '
            'pretrained loader (downloaded weights). model_id=%s',
            model, resolved_model_id,
        )
        orbff, adapter = loader(device=device, compile=compile)
    orbff.eval()
    return OrbCalculatorAdapter(
        orbff, adapter,
        device=device,
        name=f'orb-{model}',
        charge=charge,
        spin=spin,
        model_id=resolved_model_id,
    )


class OrbCalculatorAdapter(Calculator):
    """Wraps OrbMol model+adapter to provide kagome's Calculator interface."""

    def __init__(
        self,
        model,
        atoms_adapter,
        device: str = 'cpu',
        name: str = 'orb-orbmol_v2',
        charge: int = 0,
        spin: int = 1,
        model_id: str = '',
    ) -> None:
        self._model = model
        self._adapter = atoms_adapter
        self._device = device
        self._name = name
        self._charge = charge
        self._spin = spin
        self._model_id = model_id or name
        self._pbc_checked = False

        try:
            from ase import Atoms
        except ImportError:
            raise ImportError('ASE is required for this backend.')
        self._Atoms = Atoms

    @property
    def name(self) -> str:
        return self._name

    @property
    def model_id(self) -> str:
        return self._model_id

    @property
    def supports_spin(self) -> bool:
        return True

    def set_spin(self, spin: int) -> None:
        self._spin = spin

    def _check_periodic_support(self) -> None:
        """Fail early and clearly if a periodic run needs nvalchemiops but it is
        unavailable (blocked_pending_review; torch.compile failures on Windows).

        Only fires when the dependency is genuinely missing, so a properly
        configured periodic run is unaffected (RF20)."""
        import importlib.util
        if importlib.util.find_spec('nvalchemiops') is None:
            raise RuntimeError(
                'Periodic OrbMol-v2 (cell != None) requires nvalchemiops for PME, '
                'which is blocked_pending_review (license unconfirmed) and triggers '
                'torch.compile failures on Windows. Run non-periodic (cell=None), '
                'or resolve nvalchemiops first. See specs/decisions.md and '
                'specs/approved_dependencies.yaml.'
            )

    def compute(
        self,
        positions: NDArray[np.floating],
        species: list[str],
        cell: NDArray[np.floating] | None = None,
    ) -> tuple[float, NDArray[np.floating]]:
        if cell is not None and not self._pbc_checked:
            self._pbc_checked = True
            self._check_periodic_support()
        atoms = self._Atoms(
            symbols=species,
            positions=positions,
            cell=cell,
            pbc=cell is not None,
        )
        atoms.info['charge'] = self._charge
        atoms.info['spin'] = self._spin

        batch = self._adapter.from_ase_atoms(atoms, device=self._device)
        result = self._model(batch)

        energy_ev = float(result['energy'].detach().item())
        forces_ev = result['grad_forces'].detach().cpu().numpy()

        energy = energy_ev * EV_TO_KCAL_MOL
        forces = forces_ev * EV_TO_KCAL_MOL

        # Release the per-call autograd workspace back to the allocator. The
        # neighbour-graph size varies per step, so without this the CUDA caching
        # allocator fragments and reserved VRAM creeps up until a long paper-scale
        # run exhausts memory and hangs (observed at 2520 atoms). batch/result are
        # dropped first so their tensors are freeable. See specs/decisions.md
        # 2026-06-15 VRAM record.
        if self._device == 'cuda':
            del batch, result
            import torch
            torch.cuda.empty_cache()

        return energy, forces
