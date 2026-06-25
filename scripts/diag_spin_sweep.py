"""Empirical OrbMol-v2 spin sweep: how energy/forces respond to the global
spin_multiplicity feature, and where high multiplicity destabilizes.

Builds N isolated methyl radicals (each contributes one unpaired electron) and
evaluates the same fixed geometry at a range of claimed spin multiplicities.
If max|F| stays bounded for low mult (1-3) but explodes toward the high-spin sum
(N+1) and beyond, that confirms OrbMol-v2 conditions on spin as a global scalar
and extrapolates catastrophically out of its training range — the mechanism
behind the spin=21 thermal blow-up at 20 radicals.

Usage:
    python scripts/diag_spin_sweep.py --device cuda --n 8
"""
from __future__ import annotations

import argparse
import itertools

import numpy as np

from scripts._systems import _rdkit_mol
from kagome.backends.orb_backend import create_orb_calculator


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--device', default='cuda')
    ap.add_argument('--n', type=int, default=8, help='number of methyl radicals')
    ap.add_argument('--spacing', type=float, default=4.5)
    args = ap.parse_args()

    rad = _rdkit_mol('[CH3]', seed=1)
    rp = np.array(rad.GetConformer().GetPositions())
    rsp = [a.GetSymbol() for a in rad.GetAtoms()]

    grid = list(itertools.product(range(4), repeat=3))
    pos_list, species = [], []
    for k in range(args.n):
        off = np.array(grid[k], dtype=float) * args.spacing
        pos_list.append(rp + off)
        species += rsp
    pos = np.vstack(pos_list)

    print(f'system: {args.n} methyl radicals, {len(species)} atoms, '
          f'high-spin-sum multiplicity = N+1 = {args.n + 1}')
    print(f'{"mult":>5} {"E(kcal/mol)":>16} {"max|F|(kcal/mol/A)":>20}')
    print('-' * 44)

    calc = create_orb_calculator(device=args.device, spin=1)
    mults = sorted({1, 2, 3, 4, 6, 9, args.n + 1, 15, 21, 31})
    for mult in mults:
        calc.set_spin(mult)
        E, F = calc.compute(pos, species, cell=None)
        fmax = float(np.sqrt((F ** 2).sum(axis=1).max()))
        flag = '  <-- N+1 (high-spin sum)' if mult == args.n + 1 else ''
        print(f'{mult:5d} {E:16.1f} {fmax:20.1f}{flag}')


if __name__ == '__main__':
    main()
