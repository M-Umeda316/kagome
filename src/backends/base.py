"""Abstract calculator interface for MLIP backends."""
from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import NDArray


class Calculator(ABC):
    """Backend-agnostic energy/force calculator."""

    @abstractmethod
    def compute(
        self,
        positions: NDArray[np.floating],
        species: list[str],
        cell: NDArray[np.floating] | None = None,
    ) -> tuple[float, NDArray[np.floating]]:
        """Return (energy_kcal_mol, forces_kcal_mol_per_A) for the given configuration.

        positions: (n_atoms, 3) in Å
        species: element symbols, length n_atoms
        cell: (3, 3) lattice vectors in Å, or None for non-periodic
        """

    @property
    @abstractmethod
    def name(self) -> str:
        """Backend identifier for logging."""

    def set_spin(self, spin: int) -> None:
        """Update the system spin multiplicity (2S+1).

        No-op by default; override in backends that support spin (e.g. OrbMol-v2).
        """
