"""Unit tests for the FIRE energy minimizer.

Paper anchor: PDF p.20 — equilibration precedes reactive production; the
minimizer removes initial close contacts before dynamics.
"""
from __future__ import annotations

import numpy as np
import pytest

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


# ---------------------------------------------------------------------------
# Canonical-FIRE-order verification (I1/I2/I3)
# ---------------------------------------------------------------------------
#
# A hand-written reference that is an ASE-independent literal transcription of
# ``ase.optimize.fire.FIRE.step`` (Bitzek et al., PRL 97, 170201 (2006)):
#
#   step 0:   v is None -> v = 0, skip the power/mixing branch (no dt/alpha
#             change), then kick  v += dt*F.
#   step k>0: power = F.v; if power > 0: mix on the *pre-kick* v and |v|,
#             then (check N_min, then count) grow dt/shrink alpha; else reset
#             v=0, alpha=alpha_start, dt*=f_dec.  Kick  v += dt*F  afterwards.
#   displacement limited by the global norm of dr (matches ASE); a gentle
#   two-atom system is used so the clamp never fires and per-atom vs global
#   limiting is irrelevant.


def _canonical_fire_positions(
    positions: np.ndarray,
    species: list[str],
    calc: ToyCalculator,
    p: FireParams,
    n_steps: int,
) -> list[np.ndarray]:
    """Reference FIRE trajectory (position after each of ``n_steps`` steps)."""
    pos = np.array(positions, dtype=np.float64, copy=True)
    vel: np.ndarray | None = None
    dt = p.dt_start
    alpha = p.alpha_start
    n_positive = 0

    _, forces = calc.compute(pos, species, None)
    out: list[np.ndarray] = []
    for _ in range(n_steps):
        if vel is None:
            vel = np.zeros_like(pos)
        else:
            power = float(np.vdot(forces, vel))
            if power > 0.0:
                fnorm = float(np.sqrt(np.vdot(forces, forces)))
                vnorm = float(np.sqrt(np.vdot(vel, vel)))
                vel = (1.0 - alpha) * vel + alpha * forces / fnorm * vnorm
                if n_positive > p.n_min:
                    dt = min(dt * p.f_inc, p.dt_max)
                    alpha *= p.f_alpha
                n_positive += 1
            else:
                vel[:] = 0.0
                alpha = p.alpha_start
                dt *= p.f_dec
                n_positive = 0

        vel = vel + dt * forces
        dr = dt * vel
        normdr = float(np.sqrt(np.vdot(dr, dr)))
        if normdr > p.maxstep_A:
            dr = p.maxstep_A * dr / normdr
        pos = pos + dr
        out.append(pos.copy())
        _, forces = calc.compute(pos, species, None)
    return out


def test_fire_matches_canonical_reference_step_by_step() -> None:
    """fire_minimize reproduces the canonical (ASE-order) FIRE trajectory."""
    calc = ToyCalculator(epsilon=0.1, sigma=1.5)
    # Gently compressed pair: forces stay small so the maxstep clamp never
    # fires and the dynamics exercises several positive-power steps.
    positions = np.array([[0.0, 0.0, 0.0], [1.30, 0.0, 0.0]], dtype=np.float64)
    species = ['C', 'C']
    n_steps = 12
    # fmax threshold ~0 so the run does exactly n_steps steps (no early stop).
    p = FireParams(fmax_kcal_mol_A=1e-12, max_steps=n_steps, maxstep_A=0.2)

    ref = _canonical_fire_positions(positions, species, calc, p, n_steps)

    captured: list[np.ndarray] = []

    def _on_step(step, pos, energy, fmax):  # noqa: ANN001, ARG001
        if step >= 1:                       # step 0 is the initial config
            captured.append(np.array(pos, copy=True))

    fire_minimize(positions, species, None, calc, p, on_step=_on_step)

    assert len(captured) == n_steps
    for k, (a, b) in enumerate(zip(captured, ref)):
        assert np.allclose(a, b, atol=1e-12, rtol=0.0), f'step {k + 1} mismatch'


def test_fire_mixing_is_applied_before_the_kick() -> None:
    """The first mixing step must use the pre-kick velocity (I1 regression).

    On step 1 v is zero, so the (correct) pre-kick mixing is a no-op and the
    displacement is exactly dt0**2 * F0.  The old kick-then-mix order would
    instead blend the just-kicked velocity toward F, changing |dr|.  For a
    collinear pair the direction is identical, so we check the magnitude.
    """
    calc = ToyCalculator(epsilon=0.1, sigma=1.5)
    positions = np.array([[0.0, 0.0, 0.0], [1.30, 0.0, 0.0]], dtype=np.float64)
    species = ['C', 'C']
    p = FireParams(fmax_kcal_mol_A=1e-12, max_steps=1, maxstep_A=10.0)

    _, f0 = calc.compute(positions, species, None)

    first: list[np.ndarray] = []

    def _on_step(step, pos, energy, fmax):  # noqa: ANN001, ARG001
        if step == 1:
            first.append(np.array(pos, copy=True))

    fire_minimize(positions, species, None, calc, p, on_step=_on_step)

    expected = positions + p.dt_start * (p.dt_start * f0)   # dr = dt*(dt*F)
    assert np.allclose(first[0], expected, atol=1e-12)


@pytest.mark.skipif(
    __import__('importlib').util.find_spec('ase') is None,
    reason='ase not installed',
)
def test_fire_matches_real_ase_fire() -> None:
    """fire_minimize matches ase.optimize.FIRE on the same numeric potential.

    Guards against a shared bug between fire_minimize and the hand-written
    reference above.  A gentle pair keeps the displacement under maxstep so
    ASE's global-norm clamp and our per-atom clamp never diverge.
    """
    from ase import Atoms
    from ase.calculators.calculator import Calculator as ASECalculator
    from ase.calculators.calculator import all_changes
    from ase.optimize import FIRE as ASEFIRE

    calc = ToyCalculator(epsilon=0.1, sigma=1.5)
    positions = np.array([[0.0, 0.0, 0.0], [1.30, 0.0, 0.0]], dtype=np.float64)
    species = ['C', 'C']
    n_steps = 12
    p = FireParams(
        fmax_kcal_mol_A=1e-12, max_steps=n_steps, maxstep_A=0.2,
        dt_start=0.1, dt_max=1.0, n_min=5, f_inc=1.1, f_dec=0.5,
        alpha_start=0.1, f_alpha=0.99,
    )

    class _ToyASE(ASECalculator):
        implemented_properties = ['energy', 'forces']

        def calculate(self, atoms=None, properties=('energy',),
                      system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            e, f = calc.compute(
                atoms.get_positions(),
                list(atoms.get_chemical_symbols()),
                None,
            )
            self.results['energy'] = float(e)
            self.results['forces'] = np.asarray(f, dtype=np.float64)

    atoms = Atoms(species, positions=positions.copy())
    atoms.calc = _ToyASE()
    opt = ASEFIRE(
        atoms, dt=p.dt_start, maxstep=p.maxstep_A, dtmax=p.dt_max,
        Nmin=p.n_min, finc=p.f_inc, fdec=p.f_dec, astart=p.alpha_start,
        fa=p.f_alpha, a=p.alpha_start, logfile=None,
    )
    ase_positions: list[np.ndarray] = []
    for _ in range(n_steps):
        opt.step()
        ase_positions.append(atoms.get_positions().copy())

    captured: list[np.ndarray] = []

    def _on_step(step, pos, energy, fmax):  # noqa: ANN001, ARG001
        if step >= 1:
            captured.append(np.array(pos, copy=True))

    fire_minimize(positions, species, None, calc, p, on_step=_on_step)

    assert len(captured) == n_steps
    for k, (a, b) in enumerate(zip(captured, ase_positions)):
        assert np.allclose(a, b, atol=1e-9, rtol=0.0), f'ASE step {k + 1} diff'
