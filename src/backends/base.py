"""Abstract calculator interface for MLIP backends."""
from __future__ import annotations

import logging
from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import NDArray

logger = logging.getLogger(__name__)


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

    @property
    def model_id(self) -> str:
        """Resolved model/weights identity for provenance (RF17).

        For uMLIP backends the exact weights are a first-order provenance item;
        override to return the resolved checkpoint path/hash. Defaults to ``name``.
        """
        return self.name

    @property
    def supports_spin(self) -> bool:
        """Whether this backend actually applies spin multiplicity (RF20).

        Default False (spin-agnostic, e.g. MACE-MP-0). Callers can branch on this
        instead of assuming ``set_spin`` took effect.
        """
        return False

    def set_spin(self, spin: int) -> None:
        """Update the system spin multiplicity (2S+1).

        Default: no-op that WARNS, so a caller expecting spin to take effect has
        an audit trail rather than a silently ignored request (RF20). Override in
        spin-aware backends (e.g. OrbMol-v2) — and set ``supports_spin = True``.
        """
        logger.warning(
            '%s ignores spin multiplicity (supports_spin=False); requested '
            'spin=%d has no effect.', type(self).__name__, spin,
        )
