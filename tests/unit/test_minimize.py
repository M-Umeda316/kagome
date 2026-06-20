"""Unit tests for the FIRE energy minimizer.

Paper anchor: PDF p.20 — equilibration precedes reactive production; the
minimizer removes initial close contacts before dynamics.
"""
from __future__ import annotations

import numpy as np

from kagome.backends.toy import ToyCalculator
from kagome.integrators.minimize import FireParams, compress_box, fire_minimize


def test_fire_relaxes_clashing_pair_to_lj_minimum() -> None:
    """A clashing LJ pair relaxes toward r = sigma with reduced force."""
    calc = ToyCalculator(epsilon=0.1, sigma=1.5)
    # Two atoms much closer than sigma -> huge repulsive force (a "clash").
    positions = np.array([[0.0, 0.0, 0.0], [0.9, 0.0, 0.0]], dtype=np.float64)
    species = ['C', 'C']

    _, f0 = calc.compute(positions, species, None)
    fmax0 = float(np.sqrt((f0 ** 2).sum(axis=1).max()))

    result = fire_minimize(
        positions, species, None, calc,
        FireParams(fmax_kcal_mol_A=1e-3, max_steps=2000, maxstep_A=0.05),
    )

    assert result.fmax < fmax0           # force reduced
    assert result.converged              # reached threshold
    # LJ minimum of V = eps((s/r)^12 - 2(s/r)^6) is at r = sigma.
    r_final = float(np.linalg.norm(result.positions[1] - result.positions[0]))
    assert abs(r_final - 1.5) < 0.05


def test_fire_energy_does_not_increase() -> None:
    """Minimization never raises the energy of the configuration."""
    calc = ToyCalculator()
    rng = np.random.default_rng(0)
    positions = rng.uniform(0.0, 3.0, size=(5, 3))
    species = ['C'] * 5

    e0, _ = calc.compute(positions, species, None)
    result = fire_minimize(
        positions, species, None, calc,
        FireParams(fmax_kcal_mol_A=0.01, max_steps=1000),
    )
    assert result.energy <= e0 + 1e-9


def test_fire_already_converged_is_noop() -> None:
    """A well-separated pair (near-zero force) returns immediately."""
    calc = ToyCalculator(epsilon=0.1, sigma=1.5)
    positions = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]], dtype=np.float64)
    species = ['C', 'C']
    result = fire_minimize(positions, species, None, calc, FireParams())
    assert result.n_steps == 0
    assert result.converged


def test_compress_box_reaches_target_edge() -> None:
    """compress_box shrinks the cell to the target edge and scales positions."""
    calc = ToyCalculator(epsilon=0.1, sigma=1.5)
    rng = np.random.default_rng(1)
    edge0 = 20.0
    positions = rng.uniform(2.0, edge0 - 2.0, size=(8, 3))
    species = ['C'] * 8
    cell = np.diag([edge0, edge0, edge0])

    result = compress_box(
        positions, cell, target_edge_A=12.0, species=species, calculator=calc,
        n_stages=10,
    )
    assert abs(float(result.cell[0, 0]) - 12.0) < 1e-6
    assert result.n_stages == 10
    # Atoms must stay inside the (now smaller) box and remain finite.
    assert np.all(np.isfinite(result.positions))


def test_compress_box_noop_when_target_larger() -> None:
    """Requesting a larger box is a no-op (compression only)."""
    calc = ToyCalculator()
    positions = np.array([[1.0, 1.0, 1.0], [3.0, 3.0, 3.0]], dtype=np.float64)
    species = ['C', 'C']
    cell = np.diag([10.0, 10.0, 10.0])
    result = compress_box(
        positions, cell, target_edge_A=15.0, species=species, calculator=calc,
    )
    assert result.n_stages == 0
    assert float(result.cell[0, 0]) == 10.0
