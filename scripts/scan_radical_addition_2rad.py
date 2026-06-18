"""PES scan comparing high-spin vs low-spin for a 2-radical system.

S3 validation: verify that the high-spin approximation (spin=3, triplet) gives
the same radical-addition PES as spin=1 (singlet) when a spectator radical is
present 7 A away. If the energy difference is < 1 kcal/mol, the high-spin
approximation is justified for multi-radical polymerization with OrbMol-v2.

System: 2 methyl radicals + 1 ethylene.
  - Radical 1 (scanning): approaches ethylene C0 along the pi-face normal
  - Radical 2 (spectator): placed 7 A from ethylene center, held fixed

Usage:
    python scripts/scan_radical_addition_2rad.py --device cuda
    python scripts/scan_radical_addition_2rad.py --device cuda --output-dir runs/s3_pes
"""
from __future__ import annotations

import argparse
import json
import os

os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')

from pathlib import Path

import numpy as np

from scripts._systems import _rdkit_mol
from scripts.scan_radical_addition import (
    _build_fragments,
    _rotation_align,
    _unit,
    constrained_relax,
)


def _build_system_2rad(r_scan: float):
    """Build 2-radical + ethylene system with scanning radical at distance r_scan.

    Returns (positions, species, eth_C0_index, scanning_radical_C_index).
    """
    (rp, rsp, rc, lobe), (ep, esp, c0, normal) = _build_fragments()

    R = _rotation_align(lobe, -normal)
    rp_centered = rp - rp[rc]
    rp_oriented = rp_centered @ R.T

    rad1_pos = rp_oriented + (ep[c0] + r_scan * normal)

    rad2 = _rdkit_mol('[CH3]', seed=3)
    rp2 = np.array(rad2.GetConformer().GetPositions(), dtype=np.float64)
    rsp2 = [a.GetSymbol() for a in rad2.GetAtoms()]
    rc2 = next(i for i, s in enumerate(rsp2) if s == 'C')
    rp2_centered = rp2 - rp2[rc2]

    eth_center = ep[:len([s for s in esp if s == 'C'])].mean(axis=0)
    spectator_dir = _unit(np.array([1.0, 1.0, 0.0]))
    spectator_offset = eth_center + 7.0 * spectator_dir
    rad2_pos = rp2_centered + spectator_offset

    species = esp + rsp + rsp2
    positions = np.vstack([ep, rad1_pos, rad2_pos])

    rad1_C_global = len(esp) + rc
    return positions, species, c0, rad1_C_global


def run_scan(calc, spin_label: str):
    """Run PES scan and return list of (r, energy, dE) tuples."""
    radius_list = [3.5, 3.2, 3.0, 2.8, 2.6, 2.4, 2.2, 2.0, 1.8, 1.6, 1.54]
    print(f'\n--- spin={spin_label} (constrained-relaxed scan) ---')
    print(f'{"r_CC(A)":>8} {"E(kcal/mol)":>14} {"dE_vs_3.5":>12}')

    e_ref = None
    rows = []
    for r in radius_list:
        pos, species, c0, rad1_C = _build_system_2rad(r)
        energy, _ = constrained_relax(pos, species, calc, c0, rad1_C, r)
        if e_ref is None:
            e_ref = energy
        rows.append((r, energy, energy - e_ref))
        print(f'{r:8.2f} {energy:14.2f} {energy - e_ref:12.2f}', flush=True)

    return rows


def plot_comparison(rows_s1, rows_s3, output_dir: Path, rows_1rad=None):
    """Plot overlay of PES curves and save."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib not available -- skipping plot.')
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    r_vals = [row[0] for row in rows_s3]
    de_s1 = [row[2] for row in rows_s1]
    de_s3 = [row[2] for row in rows_s3]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True,
                                    gridspec_kw={'height_ratios': [3, 1]})

    ax1.plot(r_vals, de_s1, 'o-', color='tab:blue', label='2rad spin=1 (singlet)',
             linewidth=1.5, alpha=0.5)
    ax1.plot(r_vals, de_s3, 's-', color='tab:red', label='2rad spin=3 (triplet)',
             linewidth=1.5)
    if rows_1rad is not None:
        de_1rad = [row[2] for row in rows_1rad]
        ax1.plot(r_vals, de_1rad, '^--', color='tab:green',
                 label='1rad spin=2 (doublet)', linewidth=1.5)
    ax1.axhline(0, color='gray', linestyle=':', linewidth=0.5)
    ax1.set_ylabel('dE vs r=3.5 A (kcal/mol)')
    ax1.set_title('S3 PES validation: high-spin approximation\n'
                  '2 CH3 radicals + C2H4 (spectator at 7 A)')
    ax1.legend(fontsize=9)
    ax1.invert_xaxis()

    if rows_1rad is not None:
        de_1rad = [row[2] for row in rows_1rad]
        diff_correct = [d1 - d3 for d1, d3 in zip(de_1rad, de_s3)]
        ax2.bar(r_vals, diff_correct, width=0.12, color='tab:green', alpha=0.7,
                label='1rad/s=2 - 2rad/s=3')
    else:
        diff_correct = [s1 - s3 for s1, s3 in zip(de_s1, de_s3)]
        ax2.bar(r_vals, diff_correct, width=0.12, color='tab:green', alpha=0.7)
    ax2.axhline(0, color='gray', linestyle=':', linewidth=0.5)
    ax2.axhline(1.0, color='red', linestyle='--', linewidth=0.8, label='1 kcal/mol')
    ax2.axhline(-1.0, color='red', linestyle='--', linewidth=0.8)
    ax2.set_xlabel('r_CC (A)')
    ax2.set_ylabel('dE difference\n(kcal/mol)')
    ax2.legend(fontsize=8)
    ax2.invert_xaxis()

    fig.tight_layout()
    for fmt in ('png', 'pdf'):
        fig.savefig(output_dir / f's3_pes_comparison.{fmt}', dpi=150)
    plt.close(fig)
    print(f'\nSaved s3_pes_comparison.png/.pdf to {output_dir}')


def run_scan_1rad(calc, spin_label: str):
    """Run 1-radical PES scan (baseline) using scan_radical_addition logic."""
    radius_list = [3.5, 3.2, 3.0, 2.8, 2.6, 2.4, 2.2, 2.0, 1.8, 1.6, 1.54]
    print(f'\n--- 1-radical spin={spin_label} (constrained-relaxed scan) ---')
    print(f'{"r_CC(A)":>8} {"E(kcal/mol)":>14} {"dE_vs_3.5":>12}')

    (rp, rsp, rc, lobe), (ep, esp, c0, normal) = _build_fragments()
    R = _rotation_align(lobe, -normal)
    rp_centered = rp - rp[rc]
    rp_oriented = rp_centered @ R.T
    species_1rad = esp + rsp
    rad_C_global = len(esp) + rc

    e_ref = None
    rows = []
    for r in radius_list:
        rad_pos = rp_oriented + (ep[c0] + r * normal)
        pos0 = np.vstack([ep, rad_pos])
        energy, _ = constrained_relax(pos0, species_1rad, calc, c0, rad_C_global, r)
        if e_ref is None:
            e_ref = energy
        rows.append((r, energy, energy - e_ref))
        print(f'{r:8.2f} {energy:14.2f} {energy - e_ref:12.2f}', flush=True)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(
        description='S3: Compare high-spin vs low-spin PES for 2-radical system')
    ap.add_argument('--backend', default='orb', choices=['orb', 'mace'])
    ap.add_argument('--device', default='cpu')
    ap.add_argument('--output-dir', type=Path, default=Path('runs/s3_pes'))
    args = ap.parse_args()

    if args.backend == 'orb':
        from src.backends.orb_backend import create_orb_calculator

    print('=== S3 PES validation: high-spin approximation ===')
    print('System: 2 x CH3 radical + C2H4')
    print('Radical 1: scanning (r = 3.5 -> 1.54 A)')
    print('Radical 2: spectator at 7 A from ethylene center')

    if args.backend == 'orb':
        calc_s1 = create_orb_calculator(device=args.device, spin=1)
        calc_s3 = create_orb_calculator(device=args.device, spin=3)
        calc_1rad = create_orb_calculator(device=args.device, spin=2)
    else:
        from src.backends.mace_backend import create_mace_calculator
        calc_s1 = create_mace_calculator(device=args.device)
        calc_s3 = calc_s1
        calc_1rad = calc_s1
        print('WARNING: MACE backend ignores spin -- comparison is trivial.')

    rows_1rad = run_scan_1rad(calc_1rad, '2')
    rows_s1 = run_scan(calc_s1, '1')
    rows_s3 = run_scan(calc_s3, '3')

    print('\n=== Comparison: 1rad/s=2 vs 2rad/s=3 (correct validation) ===')
    print(f'{"r_CC(A)":>8} {"1rad/s=2":>12} {"2rad/s=3":>12} {"diff":>10}')
    max_diff_correct = 0.0
    for (r, _, de_1r), (_, _, de_s3) in zip(rows_1rad, rows_s3):
        diff = de_1r - de_s3
        max_diff_correct = max(max_diff_correct, abs(diff))
        print(f'{r:8.2f} {de_1r:12.2f} {de_s3:12.2f} {diff:10.2f}')

    de_1r_arr = [row[2] for row in rows_1rad]
    de_s3_arr = [row[2] for row in rows_s3]
    barrier_1rad = max(de_1r_arr) - min(de_1r_arr)
    barrier_s3 = max(de_s3_arr) - min(de_s3_arr)
    barrier_diff = abs(barrier_1rad - barrier_s3)
    min_1rad = min(de_1r_arr)
    min_s3 = min(de_s3_arr)
    product_diff = abs(min_1rad - min_s3)

    print(f'\nBarrier height:   1rad/s=2 {barrier_1rad:.2f}, 2rad/s=3 {barrier_s3:.2f}, '
          f'diff={barrier_diff:.2f} kcal/mol')
    print(f'Product minimum:  1rad/s=2 {min_1rad:.2f}, 2rad/s=3 {min_s3:.2f}, '
          f'diff={product_diff:.2f} kcal/mol')
    print(f'Max pointwise diff: {max_diff_correct:.2f} kcal/mol')

    if barrier_diff < 1.0 and max_diff_correct < 3.0:
        verdict = 'PASS'
        print(f'\nVERDICT: {verdict} -- high-spin approximation is valid for TDBB')
        print(f'  Barrier diff {barrier_diff:.2f} < 1 kcal/mol (kinetics preserved)')
        print(f'  Max pointwise diff {max_diff_correct:.2f} < 3 kcal/mol')
    elif max_diff_correct < 3.0:
        verdict = 'MARGINAL'
        print(f'\nVERDICT: {verdict} -- high-spin approximation is usable with caveat')
    else:
        verdict = 'FAIL'
        print(f'\nVERDICT: {verdict} -- high-spin approximation is NOT valid')

    print('\n=== Reference: 2rad spin=1 vs spin=3 ===')
    print('(spin=1 is closed-shell singlet, NOT physical for 2 radicals in DFT)')
    print(f'{"r_CC(A)":>8} {"2rad/s=1":>12} {"2rad/s=3":>12} {"diff":>10}')
    for (r, _, de1), (_, _, de3) in zip(rows_s1, rows_s3):
        print(f'{r:8.2f} {de1:12.2f} {de3:12.2f} {de1 - de3:10.2f}')

    plot_comparison(rows_s1, rows_s3, args.output_dir, rows_1rad=rows_1rad)

    results = {
        'verdict': verdict,
        'correct_comparison': '1rad/spin=2 vs 2rad/spin=3',
        'barrier_diff_kcal': round(barrier_diff, 4),
        'product_diff_kcal': round(product_diff, 4),
        'max_pointwise_diff_kcal': round(max_diff_correct, 4),
        'barrier_1rad_kcal': round(barrier_1rad, 4),
        'barrier_2rad_s3_kcal': round(barrier_s3, 4),
        'product_1rad_kcal': round(min_1rad, 4),
        'product_2rad_s3_kcal': round(min_s3, 4),
        'scan_1rad_s2': [(r, round(e, 4), round(de, 4)) for r, e, de in rows_1rad],
        'scan_2rad_s1': [(r, round(e, 4), round(de, 4)) for r, e, de in rows_s1],
        'scan_2rad_s3': [(r, round(e, 4), round(de, 4)) for r, e, de in rows_s3],
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with open(args.output_dir / 'pes_comparison.json', 'w') as f:
        json.dump(results, f, indent=2)
    print(f'Saved pes_comparison.json to {args.output_dir}')


if __name__ == '__main__':
    main()
