"""OrbMol-v2 backend.

Code: Apache-2.0, Model weights: Apache-2.0.
See specs/dependency-license-matrix.md.

OrbMol-v2 uses a conservative regressor: forces = -dE/dr (grad_forces).
charge and spin must be set on ASE Atoms via atoms.info.

Windows note: nvalchemiops PME triggers torch._inductor C++ compilation
which requires cl.exe.  TORCHDYNAMO_DISABLE=1 bypasses this.
"""
from __future__ import annotations

import os

import numpy as np
from numpy.typing import NDArray

from src.backends.base import Calculator

# 1 eV = 23.060548 kcal/mol (NIST CODATA 2018)
EV_TO_KCAL_MOL = 23.060548


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
    # Reduce CUDA allocator fragmentation on long runs where the neighbour-graph
    # size varies per step (else reserved VRAM creeps up until the run hangs).
    # Only effective if set before torch is first imported.
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

    orbff, adapter = loader(device=device, compile=compile)
    orbff.eval()
    return OrbCalculatorAdapter(
        orbff, adapter,
        device=device,
        name=f'orb-{model}',
        charge=charge,
        spin=spin,
    )


class OrbCalculatorAdapter(Calculator):
    """Wraps OrbMol model+adapter to provide pfpoly's Calculator interface."""

    def __init__(
        self,
        model,
        atoms_adapter,
        device: str = 'cpu',
        name: str = 'orb-orbmol_v2',
        charge: int = 0,
        spin: int = 1,
    ) -> None:
        self._model = model
        self._adapter = atoms_adapter
        self._device = device
        self._name = name
        self._charge = charge
        self._spin = spin

        try:
            from ase import Atoms
        except ImportError:
            raise ImportError('ASE is required for this backend.')
        self._Atoms = Atoms

    @property
    def name(self) -> str:
        return self._name

    def set_spin(self, spin: int) -> None:
        self._spin = spin

    def compute(
        self,
        positions: NDArray[np.floating],
        species: list[str],
        cell: NDArray[np.floating] | None = None,
    ) -> tuple[float, NDArray[np.floating]]:
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
