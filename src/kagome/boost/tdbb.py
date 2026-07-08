"""Time-Dependent Bond Boosting (TDBB) potentials and forces.

Paper: arXiv:2511.22874, Mori et al.
Equations 2-5, 8.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from kagome.geometry import minimum_image


@dataclass
class TDBBParams:
    """Parameters for TDBB bias potentials.

    Units: distances in Å, energies in kcal/mol, f2 in Å⁻², gamma in
    kcal/(mol·step) (Eq. 5 ramp slope; the paper states gamma=1.0 without a unit,
    see specs/decisions.md 2026-06-11 'Units convention for gamma').
    """
    f2: float = 10.0
    gamma: float = 1.0
    f1_max_formation: float = 250.0
    f1_max_dissociation: float = 125.0
    lambda_vdw: float = 0.60


@dataclass
class BoostState:
    """Mutable state for a single biased segment.

    The time variable ``step`` (t in Eq. 5) is measured from the start of each
    biased segment and resets to 0 when a new biased segment begins — the
    workflow constructs a fresh BoostState per biased phase. This is a paper
    interpretation; see specs/decisions.md 2026-06-11 'Time variable t ... resets
    each biased segment'.
    """
    step: int = 0
    f1_formation: float = 0.0
    f1_dissociation: float = 0.0

    def advance(self, gamma: float, f1_max_form: float, f1_max_dissoc: float) -> None:
        self.step += 1
        self.f1_formation = boost_amplitude(self.step, gamma, f1_max_form)
        self.f1_dissociation = boost_amplitude(self.step, gamma, f1_max_dissoc)

    def reset(self) -> None:
        self.step = 0
        self.f1_formation = 0.0
        self.f1_dissociation = 0.0


# --- Eq. 5 ---

def boost_amplitude(t: int, gamma: float, f1_max: float) -> float:
    """f1(t) = min(γ·t, f1_max).  Eq. 5."""
    return min(gamma * t, f1_max)


# --- Eq. 4 ---

def target_distance(vdw_radii: NDArray[np.floating], lambda_vdw: float) -> float:
    """r0 = λ · Σ r_a^vdw.  Eq. 4."""
    return lambda_vdw * float(np.sum(vdw_radii))


# --- Eq. 2 ---

def formation_potential(
    r: NDArray[np.floating],
    r0: NDArray[np.floating] | float,
    f1: float,
    f2: float,
) -> NDArray[np.floating]:
    """V^f(r,t) = f1·(1 - exp(-f2·(r - r0)²)).  Eq. 2."""
    return f1 * (1.0 - np.exp(-f2 * (r - r0) ** 2))


def formation_force_magnitude(
    r: NDArray[np.floating],
    r0: NDArray[np.floating] | float,
    f1: float,
    f2: float,
) -> NDArray[np.floating]:
    """dV^f/dr = f1·2·f2·(r - r0)·exp(-f2·(r - r0)²).

    Positive when r > r0 → attractive force toward r0.
    """
    dr = r - r0
    return f1 * 2.0 * f2 * dr * np.exp(-f2 * dr ** 2)


# --- Eq. 3 ---

def dissociation_potential(
    r: NDArray[np.floating],
    f1: float,
    f2: float,
) -> NDArray[np.floating]:
    """V^d(r,t) = f1·exp(-f2·r²).  Eq. 3."""
    return f1 * np.exp(-f2 * r ** 2)


def dissociation_force_magnitude(
    r: NDArray[np.floating],
    f1: float,
    f2: float,
) -> NDArray[np.floating]:
    """dV^d/dr = -2·f1·f2·r·exp(-f2·r²).

    Always negative → repulsive (pushes atoms apart).
    """
    return -2.0 * f1 * f2 * r * np.exp(-f2 * r ** 2)


# --- Eq. 8 ---

@dataclass
class PairBias:
    """One pair in a reaction group with its bias mode."""
    idx_a: int
    idx_b: int
    is_formation: bool
    r0: float = 0.0
    candidate_id: int = -1
    # If False, a confirmed formation of this pair is excluded from reaction
    # counting (alpha/Carothers p); bias/topology behaviour is unchanged.
    # Used for nylon water-forming k-l pairs (specs/decisions.md 2026-07-06 A5).
    counts_as_reaction: bool = True
    # Whether this pair belongs to the paper's identification set P (§2.2 step
    # 3-4): the reaction event fires / is confirmed only when ALL trigger pairs
    # of a candidate simultaneously satisfy their bonding condition. Set from
    # PairSpec.score_pair — the bias-only nylon k-l water pair (score_pair=False)
    # is is_trigger=False and never participates in the trigger conjunction.
    # Default True keeps existing callers/tests unchanged.
    is_trigger: bool = True


def total_bias(
    pairs: list[PairBias],
    positions: NDArray[np.floating],
    state: BoostState,
    params: TDBBParams,
    cell: NDArray[np.floating] | None = None,
) -> tuple[float, NDArray[np.floating]]:
    """Compute total bias energy and forces for all active pairs.

    Eq. 8: ΔV = Σ_pairs [fp·V^f(rp,t) + (1-fp)·V^d(rp,t)]

    Returns (energy, forces) where forces has shape (n_atoms, 3).
    """
    n_atoms = positions.shape[0]
    forces = np.zeros((n_atoms, 3), dtype=np.float64)
    energy = 0.0

    for pair in pairs:
        r_vec = minimum_image(
            positions[pair.idx_b] - positions[pair.idx_a], cell,
        )
        r = np.linalg.norm(r_vec)
        if r < 1e-12:
            if pair.is_formation:
                energy += state.f1_formation * (1.0 - float(np.exp(-params.f2 * pair.r0 ** 2)))
            else:
                energy += state.f1_dissociation
            continue
        e_ij = r_vec / r

        r_arr = np.array([r])

        if pair.is_formation:
            f1 = state.f1_formation
            v = formation_potential(r_arr, pair.r0, f1, params.f2)
            dv_dr = formation_force_magnitude(r_arr, pair.r0, f1, params.f2)
        else:
            f1 = state.f1_dissociation
            v = dissociation_potential(r_arr, f1, params.f2)
            dv_dr = dissociation_force_magnitude(r_arr, f1, params.f2)

        energy += float(v[0])

        # F_a = (dV/dr) · e_ij  (toward b when dV/dr > 0)
        # F_b = -(dV/dr) · e_ij
        f_mag = float(dv_dr[0])
        forces[pair.idx_a] += f_mag * e_ij
        forces[pair.idx_b] -= f_mag * e_ij

    return energy, forces


def total_bias_fast(
    pairs: list[PairBias],
    positions: NDArray[np.floating],
    state: BoostState,
    params: TDBBParams,
    box: NDArray[np.floating] | None = None,
) -> tuple[float, NDArray[np.floating], list[float]]:
    """Scalar-math fast path for ``total_bias``.

    Returns ``(energy, forces, distances)`` where ``distances[i]`` is
    the interatomic distance for ``pairs[i]``.
    ``box`` must be a prevalidated (3,) diagonal from ``validated_box``.
    """
    n_atoms = positions.shape[0]
    forces = np.zeros((n_atoms, 3), dtype=np.float64)
    energy = 0.0
    distances: list[float] = []
    f2 = params.f2
    f1_form = state.f1_formation
    f1_dissoc = state.f1_dissociation

    for pair in pairs:
        dx = positions[pair.idx_b, 0] - positions[pair.idx_a, 0]
        dy = positions[pair.idx_b, 1] - positions[pair.idx_a, 1]
        dz = positions[pair.idx_b, 2] - positions[pair.idx_a, 2]

        if box is not None:
            bx, by, bz = float(box[0]), float(box[1]), float(box[2])
            dx -= bx * round(dx / bx)
            dy -= by * round(dy / by)
            dz -= bz * round(dz / bz)

        r = math.sqrt(dx * dx + dy * dy + dz * dz)
        distances.append(r)
        if r < 1e-12:
            if pair.is_formation:
                energy += f1_form * (1.0 - math.exp(-f2 * pair.r0 * pair.r0))
            else:
                energy += f1_dissoc
            continue
        inv_r = 1.0 / r

        if pair.is_formation:
            dr = r - pair.r0
            exp_val = math.exp(-f2 * dr * dr)
            energy += f1_form * (1.0 - exp_val)
            dv_dr = f1_form * 2.0 * f2 * dr * exp_val
        else:
            exp_val = math.exp(-f2 * r * r)
            energy += f1_dissoc * exp_val
            dv_dr = -2.0 * f1_dissoc * f2 * r * exp_val

        fx = dv_dr * dx * inv_r
        fy = dv_dr * dy * inv_r
        fz = dv_dr * dz * inv_r
        forces[pair.idx_a, 0] += fx
        forces[pair.idx_a, 1] += fy
        forces[pair.idx_a, 2] += fz
        forces[pair.idx_b, 0] -= fx
        forces[pair.idx_b, 1] -= fy
        forces[pair.idx_b, 2] -= fz

    return energy, forces, distances
