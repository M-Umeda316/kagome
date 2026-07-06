"""Calibrate TDBB parameters for the active OrbMol-v2 checkpoint.

Runs two PES scans:
  1. AIBN C-N homolysis (activation barrier)
  2. Radical addition .CH3 + C2H4 (propagation barrier)

From the measured barriers, recommends f1_max, f2, and activation-steps
values. Outputs a JSON file that run_s6_paper_scale.sh can source.

Usage:
    python scripts/calibrate_tdbb.py --device cuda
    python scripts/calibrate_tdbb.py --device cuda --output-dir runs/calibration
"""
from __future__ import annotations

import argparse
import json
import math
import os
import pathlib
import sys

os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')
os.environ.setdefault('PYTORCH_CUDA_ALLOC_CONF', 'expandable_segments:True')

import numpy as np

from scripts._systems import _rdkit_3d
from scripts.scan_radical_addition import constrained_relax, _build_fragments, _rotation_align

import logging
logging.basicConfig(level=logging.INFO, format='%(name)s | %(message)s')
logger = logging.getLogger(__name__)


_AIBN_SMILES = 'CC(C)(C#N)N=NC(C)(C)C#N'
_C_AZO = 1
_N_AZO = 5


def _scan_cn(calc, positions, species, r_list):
    """Scan C-N distance, return list of {r, E, dE}."""
    rows = []
    e_ref = None
    for r in r_list:
        pos0 = np.array(positions, dtype=np.float64)
        energy, _ = constrained_relax(
            pos0, species, calc, _C_AZO, _N_AZO, r, n_steps=500, maxstep=0.05,
        )
        if e_ref is None:
            e_ref = energy
        rows.append({'r': r, 'E': energy, 'dE': energy - e_ref})
    return rows


def _scan_radical_addition(calc, r_list):
    """Scan radical-vinyl C-C distance, return list of {r, E, dE}."""
    (rp, rsp, rc, lobe), (ep, esp, c0, normal) = _build_fragments()
    R = _rotation_align(lobe, -normal)
    rp_centered = rp - rp[rc]
    rp_oriented = rp_centered @ R.T
    species = esp + rsp
    rad_C_global = len(esp) + rc

    rows = []
    e_ref = None
    for r in r_list:
        rad_pos = rp_oriented + (ep[c0] + r * normal)
        pos0 = np.vstack([ep, rad_pos])
        energy, _ = constrained_relax(pos0, species, calc, c0, rad_C_global, r)
        if e_ref is None:
            e_ref = energy
        rows.append({'r': r, 'E': energy, 'dE': energy - e_ref})
    return rows


def _recommend_params(cn_barrier, rad_barrier):
    """Recommend TDBB parameters from measured barriers."""

    # --- Activation (C-N dissociation) ---
    # f1_max must exceed the singlet barrier with margin.
    # f2=0.3 puts the bias-force extremum at |r-r0| = 1/sqrt(2*0.3) = 1.29 A
    # beyond the current bond length, wide enough to keep pushing through
    # the ~3.0 A barrier region.
    # Safety factor 2.5x on barrier for dwell-time margin.
    act_f1_max = max(250.0, math.ceil(cn_barrier * 2.5 / 25) * 25)
    act_f2 = 0.3

    # Steps: f1 saturates at f1_max/gamma steps. Beyond that is dwell time.
    # More barrier → need more dwell.  Base: 5000 for 39 kcal/mol barrier
    # (the validated value, decisions.md 2026-06-26); scale linearly for
    # different barriers, with 5000 floor.
    act_steps = max(5000, int(5000 * cn_barrier / 39.0))
    # Round to nearest 1000
    act_steps = ((act_steps + 500) // 1000) * 1000

    # --- Production (radical addition / formation) ---
    # The radical addition barrier is typically ~6 kcal/mol for OrbMol-v2.
    # f1_max_formation must exceed this; we use 250 (same as paper).
    # The real lever is f2 (capture width) — currently 2.0 for OrbMol-v2.
    prod_f1_max_formation = 250.0
    prod_f1_max_dissociation = 125.0
    prod_f2 = 2.0

    # If the radical addition barrier is significantly higher than ~6 kcal/mol,
    # the capture width might need adjustment.
    if rad_barrier > 15.0:
        # Higher barrier: keep f2=2 but may need more biased steps
        logger.warning(
            'Radical addition barrier %.1f kcal/mol is high (expected ~6). '
            'Production may need more biased-steps or lower f2.',
            rad_barrier,
        )

    return {
        'activation_f2': act_f2,
        'activation_f1_max': act_f1_max,
        'activation_steps': act_steps,
        'f2': prod_f2,
        'f1_max_formation': prod_f1_max_formation,
        'f1_max_dissociation': prod_f1_max_dissociation,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description='Calibrate TDBB parameters for active checkpoint')
    ap.add_argument('--device', default='cpu')
    ap.add_argument('--output-dir', default='runs/calibration')
    args = ap.parse_args()

    output_dir = pathlib.Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    from kagome.backends.orb_backend import create_orb_calculator

    # --- Phase 1: AIBN C-N homolysis ---
    print('=' * 60)
    print('Phase 1: AIBN C-N homolysis PES scan')
    print('=' * 60)

    calc_s1 = create_orb_calculator(device=args.device, spin=1)
    print(f'Model loaded: {calc_s1.model_id}')

    positions, species = _rdkit_3d(_AIBN_SMILES)
    positions = np.array(positions, dtype=np.float64)
    r_eq = float(np.linalg.norm(positions[_C_AZO] - positions[_N_AZO]))
    print(f'C-N equilibrium distance: {r_eq:.3f} A')

    cn_r_list = [round(r, 2) for r in np.arange(1.4, 3.6, 0.2)]

    print(f'\n{"r_CN(A)":>8} {"E(kcal/mol)":>14} {"dE":>12}')
    cn_rows = _scan_cn(calc_s1, positions, species, cn_r_list)
    for row in cn_rows:
        print(f'{row["r"]:8.2f} {row["E"]:14.2f} {row["dE"]:12.2f}')

    cn_barrier = max(row['dE'] for row in cn_rows)
    cn_barrier_r = next(row['r'] for row in cn_rows if row['dE'] == cn_barrier)
    print(f'\nSinglet barrier: {cn_barrier:.1f} kcal/mol at r={cn_barrier_r:.2f} A')

    model_id = calc_s1.model_id
    del calc_s1

    # --- Phase 2: Radical addition ---
    print('\n' + '=' * 60)
    print('Phase 2: Radical addition PES scan')
    print('=' * 60)

    calc_s2 = create_orb_calculator(device=args.device, spin=2)
    print(f'Model loaded: {calc_s2.model_id}')

    rad_r_list = [3.5, 3.2, 3.0, 2.8, 2.6, 2.4, 2.2, 2.0, 1.8, 1.6, 1.54]

    print(f'\n{"r_CC(A)":>8} {"E(kcal/mol)":>14} {"dE":>12}')
    rad_rows = _scan_radical_addition(calc_s2, rad_r_list)
    for row in rad_rows:
        print(f'{row["r"]:8.2f} {row["E"]:14.2f} {row["dE"]:12.2f}')

    de_list = [row['dE'] for row in rad_rows]
    rad_barrier = max(de_list)
    rad_barrier_r = rad_rows[de_list.index(rad_barrier)]['r']
    rad_well = min(de_list)
    rad_well_r = rad_rows[de_list.index(rad_well)]['r']
    print(f'\nBarrier: {rad_barrier:.1f} kcal/mol at r={rad_barrier_r:.2f} A')
    print(f'Well:    {rad_well:.1f} kcal/mol at r={rad_well_r:.2f} A')

    del calc_s2

    # --- Phase 3: Recommend parameters ---
    print('\n' + '=' * 60)
    print('Parameter recommendations')
    print('=' * 60)

    rec = _recommend_params(cn_barrier, rad_barrier)

    print(f'\n  Activation:')
    print(f'    --activation-f2        {rec["activation_f2"]}')
    print(f'    --activation-f1-max    {rec["activation_f1_max"]}')
    print(f'    --activation-steps     {rec["activation_steps"]}')
    print(f'\n  Production:')
    print(f'    --f2                   {rec["f2"]}')
    print(f'    --f1-max-formation     {rec["f1_max_formation"]}')
    print(f'    --f1-max-dissociation  {rec["f1_max_dissociation"]}')

    result = {
        'model_id': model_id,
        'cn_barrier_kcal_mol': cn_barrier,
        'cn_barrier_r_angstrom': cn_barrier_r,
        'radical_addition_barrier_kcal_mol': rad_barrier,
        'radical_addition_barrier_r_angstrom': rad_barrier_r,
        'radical_addition_well_kcal_mol': rad_well,
        'radical_addition_well_r_angstrom': rad_well_r,
        'scan_cn': cn_rows,
        'scan_radical_addition': rad_rows,
        'recommended': rec,
    }

    out_path = output_dir / 'calibration.json'
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(f'\nResults saved to {out_path}')

    # Also write a shell-sourceable snippet
    snippet_path = output_dir / 'tdbb_params.env'
    with open(snippet_path, 'w') as f:
        f.write(f'# Auto-generated by calibrate_tdbb.py\n')
        f.write(f'# C-N barrier: {cn_barrier:.1f} kcal/mol\n')
        f.write(f'# Radical addition barrier: {rad_barrier:.1f} kcal/mol\n')
        f.write(f'# Uses conditional assignment so manual env overrides take precedence.\n')
        f.write(f'ACTIVATION_F2="${{ACTIVATION_F2:-{rec["activation_f2"]}}}"\n')
        f.write(f'ACTIVATION_F1_MAX="${{ACTIVATION_F1_MAX:-{rec["activation_f1_max"]}}}"\n')
        f.write(f'ACTIVATION_STEPS="${{ACTIVATION_STEPS:-{rec["activation_steps"]}}}"\n')
        f.write(f'F2="${{F2:-{rec["f2"]}}}"\n')
        f.write(f'F1_MAX_FORMATION="${{F1_MAX_FORMATION:-{rec["f1_max_formation"]}}}"\n')
        f.write(f'F1_MAX_DISSOCIATION="${{F1_MAX_DISSOCIATION:-{rec["f1_max_dissociation"]}}}"\n')
    print(f'Shell snippet saved to {snippet_path}')
    print(f'\nTo apply: source {snippet_path} && bash scripts/run_s6_paper_scale.sh')


if __name__ == '__main__':
    main()
