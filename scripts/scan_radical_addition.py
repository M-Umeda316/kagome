"""PES scan of a radical addition with the MLIP backend, to test whether the
potential reproduces an attractive radical-addition channel (decisive for whether
TDBB bond formation is achievable on this MLIP at all).

Default system: methyl radical (.CH3) approaching one carbon of ethylene (C2H4)
along the pi-face normal — the textbook radical addition (barrier ~7, exothermic
~-23 kcal/mol). The forming C-C distance r is scanned; at each r a single-point
energy is computed with the backend (rigid scan: fragment internal geometries are
held, so the bonded end sits higher than a relaxed scan, but the TREND reveals an
attractive channel if one exists). Total spin = doublet (one unpaired electron).

Usage:
    python scripts/scan_radical_addition.py --device cuda
"""
from __future__ import annotations

import argparse
import os

os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')

import numpy as np

from scripts._systems import _rdkit_mol


def _unit(v):
    return v / np.linalg.norm(v)


def _build_fragments():
    """Return (radical_pos, radical_sp, radical_C, radical_lobe_dir),
    (eth_pos, eth_sp, eth_C0, eth_plane_normal)."""
    # methyl radical
    rad = _rdkit_mol('[CH3]', seed=1)
    rp = np.array(rad.GetConformer().GetPositions(), dtype=np.float64)
    rsp = [a.GetSymbol() for a in rad.GetAtoms()]
    rc = next(i for i, s in enumerate(rsp) if s == 'C')
    h_idx = [i for i, s in enumerate(rsp) if s == 'H']
    lobe = _unit(rp[rc] - rp[h_idx].mean(axis=0))  # away from the H umbrella

    # ethylene
    eth = _rdkit_mol('C=C', seed=2)
    ep = np.array(eth.GetConformer().GetPositions(), dtype=np.float64)
    esp = [a.GetSymbol() for a in eth.GetAtoms()]
    c_idx = [i for i, s in enumerate(esp) if s == 'C']
    c0 = c_idx[0]
    # plane normal from C0, the other C, and one H on C0
    other_c = c_idx[1]
    h_on_c0 = None
    for b in eth.GetBonds():
        a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if a1 == c0 and esp[a2] == 'H':
            h_on_c0 = a2
        elif a2 == c0 and esp[a1] == 'H':
            h_on_c0 = a1
    normal = _unit(np.cross(ep[other_c] - ep[c0], ep[h_on_c0] - ep[c0]))
    return (rp, rsp, rc, lobe), (ep, esp, c0, normal)


def _apply_constraint(pos, ia, ib, r):
    d_vec = pos[ib] - pos[ia]
    d = np.linalg.norm(d_vec)
    if d < 1e-9:
        d_vec = np.array([0.0, 0.0, 1.0]); d = 1.0
    corr = (r - d) / 2.0 * (d_vec / d)
    pos[ib] += corr; pos[ia] -= corr
    return pos


def constrained_relax(pos, species, calc, ia, ib, r_target, n_steps=250, maxstep=0.1):
    """FIRE relaxation of all atoms while holding |pos[ia]-pos[ib]| = r_target.

    The along-bond force component on the two constrained atoms is projected out
    so the optimizer does not fight the constraint, and the exact distance is
    re-imposed after each step. Allows rehybridization (sp2->sp3) that a rigid
    scan misses.
    """
    pos = _apply_constraint(np.array(pos, dtype=np.float64), ia, ib, r_target)
    vel = np.zeros_like(pos)
    dt, alpha, n_pos = 0.1, 0.1, 0
    energy = None
    for _ in range(n_steps):
        energy, forces = calc.compute(pos, species, cell=None)
        axis = _unit(pos[ib] - pos[ia])
        forces[ia] -= np.dot(forces[ia], axis) * axis
        forces[ib] -= np.dot(forces[ib], axis) * axis
        fmax = float(np.sqrt((forces ** 2).sum(axis=1).max()))
        if fmax < 1.0:
            break
        power = float(np.sum(forces * vel))
        if power > 0.0:
            n_pos += 1
            if n_pos > 5:
                dt = min(dt * 1.1, 0.5); alpha *= 0.99
        else:
            n_pos = 0; dt *= 0.5; alpha = 0.1; vel[:] = 0.0
        vel += dt * forces
        vn = float(np.linalg.norm(vel)); fn = float(np.linalg.norm(forces))
        if fn > 1e-12:
            vel = (1.0 - alpha) * vel + alpha * (forces / fn) * vn
        dx = dt * vel
        dxm = float(np.sqrt((dx ** 2).sum(axis=1).max()))
        if dxm > maxstep:
            dx *= maxstep / dxm
        pos += dx
        pos = _apply_constraint(pos, ia, ib, r_target)
    return energy, pos


def _rotation_align(a, b):
    """Rotation matrix aligning unit vector a onto unit vector b."""
    a = _unit(a); b = _unit(b)
    v = np.cross(a, b); c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-8:
        return np.eye(3) if c > 0 else -np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * (1.0 / (1.0 + c))


def main() -> None:
    ap = argparse.ArgumentParser(description='Radical-addition PES scan')
    ap.add_argument('--backend', default='orb', choices=['orb', 'mace'])
    ap.add_argument('--device', default='cpu')
    ap.add_argument('--spin', type=int, default=2, help='doublet for one radical')
    args = ap.parse_args()

    if args.backend == 'orb':
        from src.backends.orb_backend import create_orb_calculator
        calc = create_orb_calculator(device=args.device, spin=args.spin)
    else:
        from src.backends.mace_backend import create_mace_calculator
        calc = create_mace_calculator(device=args.device)

    (rp, rsp, rc, lobe), (ep, esp, c0, normal) = _build_fragments()
    # Orient the radical so its lobe points along -normal (toward the ethylene C0).
    R = _rotation_align(lobe, -normal)
    rp_centered = rp - rp[rc]
    rp_oriented = rp_centered @ R.T

    species = esp + rsp
    rad_C_global = len(esp) + rc
    radius_list = [3.5, 3.2, 3.0, 2.8, 2.6, 2.4, 2.2, 2.0, 1.8, 1.6, 1.54]
    print(f'backend={calc.name} spin={args.spin}  (constrained-relaxed scan)')
    print(f'{"r_CC(A)":>8} {"E(kcal/mol)":>14} {"dE_vs_3.5":>12}')
    e_ref = None
    rows = []
    for r in radius_list:
        rad_pos = rp_oriented + (ep[c0] + r * normal)  # radical C at C0 + r*normal
        pos0 = np.vstack([ep, rad_pos])
        energy, _ = constrained_relax(pos0, species, calc, c0, rad_C_global, r)
        if e_ref is None:
            e_ref = energy
        rows.append((r, energy, energy - e_ref))
        print(f'{r:8.2f} {energy:14.2f} {energy - e_ref:12.2f}', flush=True)

    de = [row[2] for row in rows]
    rmin = rows[int(np.argmin(de))][0]
    print(f'\nmin dE at r={rmin:.2f} A, dE_min={min(de):.2f} kcal/mol')
    if min(de) < -5.0 and rmin <= 2.0:
        print('VERDICT: attractive radical-addition channel present '
              '(relaxed energy drops toward bonding distance) -> OrbMol-v2 can do '
              'this reaction; formations are a sampling/cycles matter.')
    else:
        print('VERDICT: NO attractive channel even in the relaxed scan -> OrbMol-v2 '
              'does not reproduce this radical addition; a different MLIP (e.g. PFP) '
              'would be required. No TDBB tuning can fix this.')


if __name__ == '__main__':
    main()
