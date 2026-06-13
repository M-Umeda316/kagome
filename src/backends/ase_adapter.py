"""ASE Calculator adapter: wraps any ASE Calculator for use with pfpoly.

ASE is LGPL-2.1. This module uses import-only access (no modification).
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from src.backends.base import Calculator

if TYPE_CHECKING:
    from ase.calculators.calculator import Calculator as ASECalc

# eV -> kcal/mol, eV/Å -> kcal/(mol·Å)  — NIST CODATA 2018
EV_TO_KCAL_MOL = 23.060548


class ASECalculatorAdapter(Calculator):
    """Wraps an ASE Calculator to provide pfpoly's Calculator interface.

    Handles unit conversion: ASE uses eV/Å, pfpoly uses kcal/mol/Å.
    """

    def __init__(self, ase_calc: ASECalc, name: str = 'ase') -> None:
        try:
            from ase import Atoms
        except ImportError:
            raise ImportError(
                'ASE is required for this backend. '
                'Install with: pip install ase'
            )
        self._calc = ase_calc
        self._name = name
        self._Atoms = Atoms

    @property
    def name(self) -> str:
        return self._name

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
        atoms.calc = self._calc

        energy_ev = atoms.get_potential_energy()
        forces_ev_a = atoms.get_forces()

        energy = float(energy_ev) * EV_TO_KCAL_MOL
        forces = forces_ev_a * EV_TO_KCAL_MOL

        return energy, forces
