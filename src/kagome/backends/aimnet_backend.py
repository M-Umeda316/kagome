"""AIMNet2 backend (candidate 2nd MLIP for open-shell / radical chemistry).

Code: MIT (github.com/isayevlab/aimnetcentral). Weights: MIT (HF isayevlab/*).
See specs/dependency-license-matrix.md (verified 2026-06-25).

AIMNet2-NSE is purpose-trained for open-shell radical chemistry with neural
spin-charge equilibration; it takes a total spin multiplicity (mult = 2S+1).
This backend is a COMPLEMENT to OrbMol-v2 (cross-validation, radical systems),
not a replacement — see specs/decisions.md 2026-06-25.

Spike status: under evaluation (PES radical-addition validation + multi-radical
high-spin stability). Not yet a default backend.
"""
from __future__ import annotations

import logging

import numpy as np
from numpy.typing import NDArray

from kagome.backends.base import Calculator
from kagome.units import EV_TO_KCAL_MOL

logger = logging.getLogger(__name__)


def create_aimnet_calculator(
    model: str = 'aimnet2-nse',
    device: str = 'cpu',
    charge: int = 0,
    spin: int = 2,
) -> 'AimnetCalculatorAdapter':
    """Create an AIMNet2-backed calculator.

    model: 'aimnet2-nse' (open-shell/radicals), 'aimnet2' (general), etc.
    spin: total spin multiplicity (2S+1). Default 2 (doublet, single radical).
    """
    try:
        from aimnet.calculators import AIMNet2Calculator, AIMNet2ASE
    except ImportError:
        raise ImportError(
            'aimnet is required for this backend. Install with: pip install aimnet'
        )
    base = AIMNet2Calculator(model, device=device)
    ase_calc = AIMNet2ASE(base, charge=charge, mult=spin)
    return AimnetCalculatorAdapter(ase_calc, name=f'aimnet-{model}', spin=spin, charge=charge)


class AimnetCalculatorAdapter(Calculator):
    """Wraps the AIMNet2 ASE calculator to kagome's Calculator interface."""

    def __init__(self, ase_calc, name: str = 'aimnet-aimnet2-nse',
                 spin: int = 2, charge: int = 0) -> None:
        self._ase = ase_calc
        self._name = name
        self._spin = spin
        self._charge = charge
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
        return self._name

    @property
    def supports_spin(self) -> bool:
        return True

    def set_spin(self, spin: int) -> None:
        self._spin = spin
        self._ase.set_mult(spin)

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
        atoms.calc = self._ase
        self._ase.set_mult(self._spin)
        self._ase.set_charge(self._charge)
        energy_ev = float(atoms.get_potential_energy())
        forces_ev = np.asarray(atoms.get_forces(), dtype=np.float64)
        return energy_ev * EV_TO_KCAL_MOL, forces_ev * EV_TO_KCAL_MOL

    def magnetic_moments(self) -> NDArray[np.floating] | None:
        """Per-atom spin (magnetic) moments from the last compute, if available.

        Useful for inspecting where the unpaired spin localises in a multi-radical
        system (NSE spin-charge equilibration)."""
        try:
            return np.asarray(self._ase.get_magnetic_moments(), dtype=np.float64)
        except (NotImplementedError, AttributeError):
            # ASE raises PropertyNotImplementedError (a NotImplementedError
            # subclass) when the calculator has no magnetic moments; a missing
            # method raises AttributeError. Narrow the catch so genuine bugs
            # (e.g. shape/dtype errors) are not silently swallowed (B7).
            return None
