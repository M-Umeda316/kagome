"""Isotropic Monte Carlo barostat for NPT ensemble.

Paper anchor: arXiv:2511.22874 Section 2 — "NPT ensemble".
Barostat type (MC) and target pressure (1 atm) are not specified in the paper;
see specs/decisions.md "2026-06-13: NPT Monte Carlo barostat".

Acceptance criterion (isotropic NPT, N atoms):
    ΔH = ΔE + P_ext·ΔV - N·kT·ln(V_new/V_old)
    acc = min(1, exp(-ΔH / kT))

Volume proposal:
    ln(V_new/V_old) ~ Uniform(-max_ln_dV, max_ln_dV)
    => cell scales isotropically: new_cell = old_cell * (V_new/V_old)^(1/3)
    => positions scale with cell:  new_pos  = old_pos  * (V_new/V_old)^(1/3)
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from kagome.units import ATM_TO_KCAL_MOL_A3, KB

logger = logging.getLogger(__name__)


@dataclass
class MCBarostatParams:
    """Parameters for the isotropic Monte Carlo barostat.

    pressure_atm:          target pressure (atm). Default 1 atm.
    frequency:             MD steps between barostat attempts. Default 25.
    max_volume_change_frac: maximum |ln(V'/V)| per attempt. Default 0.01.
    """
    pressure_atm: float = 1.0
    frequency: int = 25
    max_volume_change_frac: float = 0.01


@dataclass
class MCBarostatStats:
    """Running acceptance statistics."""
    attempts: int = 0
    accepted: int = 0

    @property
    def acceptance_rate(self) -> float:
        return self.accepted / self.attempts if self.attempts > 0 else 0.0


class MCBarostat:
    """Isotropic NPT Monte Carlo barostat.

    Usage inside a simulation loop:
        if barostat.should_attempt(step):
            barostat.try_step(state, calculator, rng, temperature_K)
    """

    def __init__(self, params: MCBarostatParams) -> None:
        self.params = params
        self._pressure_kcal = params.pressure_atm * ATM_TO_KCAL_MOL_A3
        self.stats = MCBarostatStats()

    def should_attempt(self, step_in_phase: int) -> bool:
        return step_in_phase > 0 and step_in_phase % self.params.frequency == 0

    def try_step(
        self,
        positions: NDArray[np.floating],
        species: list[str],
        cell: NDArray[np.floating],
        current_energy: float,
        calculator: object,  # src.backends.base.Calculator
        rng: np.random.Generator,
        temperature_K: float,
    ) -> tuple[bool, float, NDArray[np.floating]]:
        """Attempt an isotropic volume change.

        Args:
            positions:      (N, 3) array, modified in-place on acceptance
            species:        atom species list
            cell:           (3, 3) diagonal cell matrix, modified in-place on acceptance
            current_energy: potential energy at current configuration (kcal/mol)
            calculator:     Calculator with compute(positions, species, cell) -> (E, F)
            rng:            random generator
            temperature_K:  current thermostat temperature (K)

        Returns:
            (accepted, new_energy, new_forces)
            new_energy and new_forces reflect the accepted configuration.
        """
        kT = KB * temperature_K
        box = np.array([cell[0, 0], cell[1, 1], cell[2, 2]])
        V_old = float(np.prod(box))
        n_atoms = positions.shape[0]

        # Volume proposal: ln(V'/V) uniform in [-max, max]
        delta_ln_V = rng.uniform(-self.params.max_volume_change_frac,
                                  self.params.max_volume_change_frac)
        scale = np.exp(delta_ln_V / 3.0)  # linear scale factor (per dimension)
        V_new = V_old * np.exp(delta_ln_V)

        # Scale positions and cell
        new_positions = positions * scale
        new_cell = cell * scale

        new_energy, new_forces = calculator.compute(new_positions, species, new_cell)

        dE = new_energy - current_energy
        dV = V_new - V_old
        # Jacobian term for a uniform-in-ln(V) proposal is -(N+1)*kT*ln(V'/V):
        # the +1 over the uniform-in-V form is the d(lnV) measure (RF19b). The
        # difference vs N is O(1/N); see specs/decisions.md 2026-06-20 RF19b.
        delta_H = dE + self._pressure_kcal * dV - (n_atoms + 1) * kT * delta_ln_V

        self.stats.attempts += 1
        if delta_H <= 0.0 or rng.random() < np.exp(-delta_H / kT):
            # Accept: apply changes in-place
            positions[:] = new_positions
            cell[:] = new_cell
            self.stats.accepted += 1
            logger.debug(
                'MCBarostat accepted: delta_H=%.3f kcal/mol, V %.3f->%.3f A^3, scale=%.5f',
                delta_H, V_old, V_new, scale,
            )
            return True, new_energy, new_forces
        else:
            logger.debug(
                'MCBarostat rejected: delta_H=%.3f kcal/mol, acc_prob=%.4f',
                delta_H, np.exp(-delta_H / kT),
            )
            return False, current_energy, None
