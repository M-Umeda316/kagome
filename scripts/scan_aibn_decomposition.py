"""PES scan of AIBN C-N bond homolysis with OrbMol-v2.

Scans one C-N azo bond distance (C1-N5) from equilibrium (~1.47 Å) to
dissociation (~3.5 Å). At each point, all other DOF are relaxed under a
distance constraint (constrained FIRE).

Two spin states are tested:
  - spin=1 (singlet, intact molecule)
  - spin=3 (triplet, after homolysis — 2 unpaired electrons)

The PES crossing point and barrier height determine whether V^d
(f1_max=125 kcal/mol) can drive the dissociation.

Paper anchor: Table S1 — Activation row uses V^d on the C-N azo bond.

Usage:
    conda run -n kagome-gpu python scripts/scan_aibn_decomposition.py --device cuda
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib

os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')

import numpy as np

from scripts._systems import _rdkit_3d
from scripts.scan_radical_addition import constrained_relax


_AIBN_SMILES = 'CC(C)(C#N)N=NC(C)(C)C#N'

# Atom indices from RDKit AddHs ordering (verified by _aibn_analyze.py):
# C1 (central left) = index 1, N5 (azo left) = index 5
_C_AZO = 1
_N_AZO = 5


def run_scan(calc, positions, species, c_idx, n_idx, r_list, label):
    """Scan C-N distance, relaxing all other DOF."""
    rows = []
    e_ref = None
    print(f'\n=== {label} ===')
    print(f'{"r_CN(A)":>8} {"E(kcal/mol)":>14} {"dE":>12}')
    for r in r_list:
        pos0 = np.array(positions, dtype=np.float64)
        energy, pos_out = constrained_relax(
            pos0, species, calc, c_idx, n_idx, r, n_steps=500, maxstep=0.05,
        )
        if e_ref is None:
            e_ref = energy
        de = energy - e_ref
        rows.append({'r': r, 'E': energy, 'dE': de})
        print(f'{r:8.2f} {energy:14.2f} {de:12.2f}', flush=True)
    return rows, e_ref


def plot_results(rows_s1, rows_s3, output_dir):
    """Plot spin=1 and spin=3 PES curves."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib not available, skipping plot')
        return

    r1 = [row['r'] for row in rows_s1]
    de1 = [row['dE'] for row in rows_s1]
    r3 = [row['r'] for row in rows_s3]
    de3 = [row['dE'] for row in rows_s3]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(r1, de1, 'bo-', label='spin=1 (singlet)', markersize=5)
    ax.plot(r3, de3, 'rs-', label='spin=3 (triplet)', markersize=5)
    ax.axhline(125, color='gray', linestyle='--', alpha=0.5, label='f1_max(dissoc)=125 kcal/mol')
    ax.set_xlabel('C-N distance (Å)')
    ax.set_ylabel('ΔE (kcal/mol)')
    ax.set_title('AIBN C-N homolysis PES (OrbMol-v2)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    for fmt in ('png', 'pdf'):
        path = output_dir / f'aibn_pes.{fmt}'
        fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Plot saved to {output_dir}/aibn_pes.png')


def main() -> None:
    ap = argparse.ArgumentParser(description='AIBN C-N homolysis PES scan')
    ap.add_argument('--device', default='cpu')
    ap.add_argument('--output-dir', default='runs/s4_aibn_pes')
    args = ap.parse_args()

    output_dir = pathlib.Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    positions, species = _rdkit_3d(_AIBN_SMILES)
    positions = np.array(positions, dtype=np.float64)

    r_eq = float(np.linalg.norm(positions[_C_AZO] - positions[_N_AZO]))
    print(f'AIBN: {len(species)} atoms, C-N equilibrium distance: {r_eq:.3f} Å')
    print(f'Scanning C{_C_AZO}-N{_N_AZO} bond')

    r_list = [round(r, 2) for r in np.arange(1.4, 3.6, 0.2)]

    from kagome.backends.orb_backend import create_orb_calculator

    # Spin=1 (singlet, intact AIBN)
    calc_s1 = create_orb_calculator(device=args.device, spin=1)
    rows_s1, eref_s1 = run_scan(calc_s1, positions, species, _C_AZO, _N_AZO, r_list, 'spin=1 (singlet)')

    # Spin=3 (triplet, after homolysis)
    calc_s3 = create_orb_calculator(device=args.device, spin=3)
    rows_s3, eref_s3 = run_scan(calc_s3, positions, species, _C_AZO, _N_AZO, r_list, 'spin=3 (triplet)')

    # Find crossing point and barrier
    de_s1 = [row['dE'] for row in rows_s1]
    de_s3_abs = [row['E'] for row in rows_s3]
    de_s1_abs = [row['E'] for row in rows_s1]

    # The dissociation barrier on the singlet surface
    barrier_s1 = max(de_s1)
    # Energy difference at large r (spin=3 vs spin=1 at equilibrium)
    e_dissoc = de_s3_abs[-1] - de_s1_abs[0]

    print(f'\n=== Summary ===')
    print(f'Singlet barrier (max dE on spin=1 curve): {barrier_s1:.1f} kcal/mol')
    print(f'Dissociation energy (spin=3 at r=3.4 vs spin=1 at eq): {e_dissoc:.1f} kcal/mol')
    print(f'V^d f1_max (paper): 125 kcal/mol')

    if barrier_s1 < 100 and e_dissoc < 125:
        verdict = 'PASS — V^d (125 kcal/mol) can drive AIBN C-N homolysis'
    else:
        verdict = 'FAIL — barrier too high for V^d to drive dissociation'
    print(f'GATE verdict: {verdict}')

    result = {
        'smiles': _AIBN_SMILES,
        'c_idx': _C_AZO,
        'n_idx': _N_AZO,
        'r_eq': r_eq,
        'scan_spin1': rows_s1,
        'scan_spin3': rows_s3,
        'barrier_s1': barrier_s1,
        'dissoc_energy': e_dissoc,
        'verdict': verdict,
    }
    with open(output_dir / 'aibn_pes.json', 'w') as f:
        json.dump(result, f, indent=2)

    plot_results(rows_s1, rows_s3, output_dir)


if __name__ == '__main__':
    main()
