"""PES + TDBB-reach scan for the nylon-6,6 amide-forming bond, to de-risk
whether OrbMol-v2 + TDBB can actually close the amine-N -- carboxyl-C bond from
the paper's [3,6] Å candidate window BEFORE committing to long nylon runs.

This is the nylon analogue of the methyl-acrylate (MA) "capture-shell / dead-zone"
diagnosis (specs/decisions.md 2026-06-26): the TDBB formation bias
    V^f(r,t) = f1·(1 - exp(-f2·(r - r0)²))   (Eq. 2)
has appreciable pulling force only inside a Gaussian capture shell of width
~1/√f2 around r0.  For MA (C-C, r0=2.04 Å) with the OrbMol-tuned f2=5 the force
was already ~0 at 3.3-3.6 Å, so pairs *selected* in [3,6] Å could never be pulled
into the shell; lowering f2 10->2 bridged the dead-zone.  Nylon uses the SAME
[3,6] Å window (see build_nylon66_system in scripts/_systems.py and run_nylon66.py
which runs the paper default f2=10) with a shorter target r0 = λ·(r_vdw^N + r_vdw^C)
= 0.60·(1.55+1.70) = 1.95 Å — so the dead-zone risk is structurally identical and
must be measured, not assumed.

What this script measures
-------------------------
1. PES: energy vs the forming N···C distance r, from a bonded distance (~1.4 Å)
   out through the capture shell and across the whole [3,6] candidate window
   (default 1.4 -> 6.0 Å).  At each fixed r the REST of the system is relaxed with
   the same constrained-FIRE technique as scan_radical_addition.py (along-bond
   force projected out on the two constrained atoms, exact distance re-imposed
   each step), so sp2->sp3 rehybridisation of the carboxyl carbon is captured.
   Reported: barrier height (max ΔE vs reactants) and the bonded-well depth
   (ΔE at the shortest r) — i.e. does the MLIP PES even offer an accessible
   channel toward a bonded C-N adduct.
2. TDBB reach: the analytic formation-bias force magnitude |dV^f/dr| (Eq. 2
   derivative) at representative distances in the capture shell vs the [3,6] Å
   selection window, for a sweep of f2 (default 10, 5, 2), with the saturated
   amplitude f1 = f1_max_formation.  This shows directly whether the bias reaches
   the pairs the selector picks, and whether a lower f2 bridges the gap — exactly
   mirroring the MA finding.  This part is analytic (no MLIP) and is always
   computed, even with --force-only, so the reach verdict is available without GPU.

Chemistry note (IMPORTANT for interpretation)
---------------------------------------------
The minimal fragments form the C-N bond of the *tetrahedral / carbinolamine-type
addition intermediate* (amine N attacks the carbonyl C); the leaving group (water
from the -OH + one amine H) is NOT eliminated in this rigid pair scan.  That is the
correct object for the de-risk question, which is only "can OrbMol+TDBB drag N and
C together to a bonded distance at all", not "what is the full condensation
thermochemistry".  Full amide+H2O energetics is out of scope here.

Default fragments: methylamine (CN) + acetic acid (CC(=O)O) — the cleanest closed-
shell (spin=1, NO multi-radical spin complication) representatives of the HMDA
amine end and the adipic-acid carboxyl end.  --monomer-fragments switches to
fragments cut from the real monomers (see _build_amide_fragments).

Usage:
    python scripts/scan_amide_formation.py --device cuda
    python scripts/scan_amide_formation.py --device cpu --backend toy   # cheap dry run
    python scripts/scan_amide_formation.py --force-only                 # analytic reach only
"""
from __future__ import annotations

import argparse
import json
import os

os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')

from pathlib import Path

import numpy as np

from scripts._systems import (
    _DIACID_SMILES,
    _DIAMINE_SMILES,
    _find_carboxyl_c_and_oh,
    _find_terminal_amine_n,
    _rdkit_mol,
)
from scripts.scan_radical_addition import (
    _rotation_align,
    _unit,
    constrained_relax,
)
from kagome.boost.tdbb import (
    TDBBParams,
    formation_force_magnitude,
    formation_potential,
    target_distance,
)
from kagome.workflows.polymerization import VDW_RADII

# Default small closed-shell fragments (see module docstring).
_AMINE_SMILES = 'CN'         # methylamine  (primary amine end, HMDA analogue)
_ACID_SMILES = 'CC(=O)O'     # acetic acid  (carboxyl end, adipic-acid analogue)


def _build_amide_fragments(
    amine_smiles: str = _AMINE_SMILES,
    acid_smiles: str = _ACID_SMILES,
) -> tuple[tuple[np.ndarray, list[str], int, np.ndarray],
           tuple[np.ndarray, list[str], int, np.ndarray]]:
    """Build the amine and carboxylic-acid fragments for the N···C scan.

    Returns
    -------
    (amine_pos, amine_species, amine_N_local, amine_lobe_dir),
    (acid_pos,  acid_species,  acid_C_local,  acid_face_normal)

    ``amine_lobe_dir`` points from the mean of N's bonded neighbours toward N —
    i.e. along the nitrogen lone pair (the nucleophile's attack direction).
    ``acid_face_normal`` is the unit normal of the sp2 carboxyl plane (the pi
    face the nucleophile approaches).  The two fragments are later oriented so the
    amine lone pair points at the carboxyl carbon along that face normal, giving a
    chemically reasonable starting guess; the constrained relaxation then relaxes
    everything else at each fixed r.
    """
    # amine
    amine = _rdkit_mol(amine_smiles, seed=1)
    ap = np.array(amine.GetConformer().GetPositions(), dtype=np.float64)
    asp = [a.GetSymbol() for a in amine.GetAtoms()]
    n_candidates = _find_terminal_amine_n(amine_smiles)
    if not n_candidates:
        raise ValueError(f'No terminal primary amine N in {amine_smiles!r}')
    n_local = n_candidates[0]
    n_atom = amine.GetAtomWithIdx(n_local)
    nbr_idx = [nb.GetIdx() for nb in n_atom.GetNeighbors()]
    lobe = _unit(ap[n_local] - ap[nbr_idx].mean(axis=0))  # along the N lone pair

    # acid
    acid = _rdkit_mol(acid_smiles, seed=2)
    cp = np.array(acid.GetConformer().GetPositions(), dtype=np.float64)
    csp = [a.GetSymbol() for a in acid.GetAtoms()]
    cooh = _find_carboxyl_c_and_oh(acid_smiles)
    if not cooh:
        raise ValueError(f'No -C(=O)OH group in {acid_smiles!r}')
    c_local = cooh[0][0]
    c_atom = acid.GetAtomWithIdx(c_local)
    c_nbrs = [nb.GetIdx() for nb in c_atom.GetNeighbors()]
    # sp2 plane normal from two distinct neighbour bond vectors.
    v1 = cp[c_nbrs[0]] - cp[c_local]
    v2 = cp[c_nbrs[1]] - cp[c_local]
    normal = _unit(np.cross(v1, v2))
    return (ap, asp, n_local, lobe), (cp, csp, c_local, normal)


def _assemble_system(
    amine_frag, acid_frag, r: float,
) -> tuple[np.ndarray, list[str], int, int]:
    """Place the amine at N···C distance r along the carboxyl face normal.

    Returns (positions, species, acid_C_global, amine_N_global).
    """
    (ap, asp, n_local, lobe) = amine_frag
    (cp, csp, c_local, normal) = acid_frag

    # Orient the amine so its lone pair points along -normal (toward the acid C).
    R = _rotation_align(lobe, -normal)
    ap_centered = ap - ap[n_local]
    ap_oriented = ap_centered @ R.T
    amine_pos = ap_oriented + (cp[c_local] + r * normal)  # amine N at acid_C + r·normal

    species = csp + asp
    positions = np.vstack([cp, amine_pos])
    acid_C_global = c_local
    amine_N_global = len(csp) + n_local
    return positions, species, acid_C_global, amine_N_global


# --- TDBB reach analysis (analytic, no MLIP) --------------------------------

def _reach_table(
    r0: float,
    f2_values: list[float],
    f1: float,
    shell_dists: list[float],
    window_dists: list[float],
) -> dict:
    """Analytic |dV^f/dr| and V^f at representative distances for each f2.

    ``shell_dists`` are inside/near the capture shell (~r0); ``window_dists`` span
    the [3,6] Å candidate-selection window.  Force magnitude is Eq. 2's derivative
    at the saturated amplitude f1 (= f1_max_formation): the strongest the bias can
    ever pull.  If it is ~0 across the window, no amount of biased time reaches the
    selected pairs (the MA dead-zone signature).
    """
    all_dists = sorted(set(shell_dists) | set(window_dists))
    out: dict = {'r0_A': round(r0, 4), 'f1_amplitude': f1, 'by_f2': {}}
    for f2 in f2_values:
        r_arr = np.array(all_dists, dtype=np.float64)
        force = np.abs(formation_force_magnitude(r_arr, r0, f1, f2))
        pot = formation_potential(r_arr, r0, f1, f2)
        rows = [
            {
                'r_A': round(float(r), 3),
                'in_window': bool(3.0 <= r <= 6.0),
                'force_kcal_mol_A': round(float(fo), 4),
                'Vf_kcal_mol': round(float(v), 4),
            }
            for r, fo, v in zip(all_dists, force, pot)
        ]
        # dead-zone metric: max pulling force anywhere in the [3,6] window.
        win_force = [row['force_kcal_mol_A'] for row in rows if row['in_window']]
        out['by_f2'][str(f2)] = {
            'max_force_in_window_kcal_mol_A': round(max(win_force), 4) if win_force else 0.0,
            'table': rows,
        }
    return out


# --- PES scan ---------------------------------------------------------------

def _run_pes_scan(calc, amine_frag, acid_frag, radius_list: list[float]):
    """Constrained-relaxed energy vs N···C distance.  Reference = largest r."""
    print(f'backend={calc.name}  (constrained-relaxed amide N...C scan)')
    print(f'{"r_NC(A)":>8} {"E(kcal/mol)":>16} {"dE_vs_react":>14}')
    # scan from far (reactants) to bonded so e_ref is the separated reactants.
    ordered = sorted(radius_list, reverse=True)
    e_ref = None
    rows = []
    for r in ordered:
        pos, species, c_glob, n_glob = _assemble_system(amine_frag, acid_frag, r)
        energy, _ = constrained_relax(pos, species, calc, c_glob, n_glob, r)
        if e_ref is None:
            e_ref = energy
        rows.append((r, energy, energy - e_ref))
        print(f'{r:8.2f} {energy:16.2f} {energy - e_ref:14.2f}', flush=True)
    # return sorted by increasing r for stable downstream/plots
    rows.sort(key=lambda x: x[0])
    return rows


def _plot(rows, reach: dict, r0: float, output_dir: Path, fragments_label: str):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib not available -- skipping plot.')
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 9))

    # --- top: PES ---
    if rows:
        r_vals = [row[0] for row in rows]
        de = [row[2] for row in rows]
        ax1.plot(r_vals, de, 'o-', color='tab:blue', linewidth=1.5,
                 label='OrbMol PES (constrained-relaxed)')
    ax1.axvspan(3.0, 6.0, color='tab:orange', alpha=0.12,
                label='candidate window [3,6] Å')
    ax1.axvline(r0, color='tab:green', linestyle='--', linewidth=1.0,
                label=f'r0 = {r0:.2f} Å (TDBB target)')
    ax1.axhline(0, color='gray', linestyle=':', linewidth=0.5)
    ax1.set_xlabel('r(amine N -- carboxyl C)  (Å)')
    ax1.set_ylabel('dE vs separated reactants (kcal/mol)')
    ax1.set_title(f'Amide-formation PES  ({fragments_label})')
    ax1.legend(fontsize=8)

    # --- bottom: TDBB reach vs f2 ---
    for f2_str, data in reach['by_f2'].items():
        rr = [row['r_A'] for row in data['table']]
        ff = [row['force_kcal_mol_A'] for row in data['table']]
        ax2.plot(rr, ff, 'o-', linewidth=1.3, markersize=4, label=f'f2={f2_str}')
    ax2.axvspan(3.0, 6.0, color='tab:orange', alpha=0.12)
    ax2.axvline(r0, color='tab:green', linestyle='--', linewidth=1.0)
    ax2.set_yscale('symlog', linthresh=1.0)
    ax2.set_xlabel('r(amine N -- carboxyl C)  (Å)')
    ax2.set_ylabel('|dV^f/dr|  (kcal/mol/Å)  [symlog]')
    ax2.set_title('TDBB formation-bias reach (f1 saturated) vs capture-shell/window')
    ax2.legend(fontsize=8)

    fig.tight_layout()
    for fmt in ('png', 'pdf'):
        fig.savefig(output_dir / f'amide_formation_scan.{fmt}', dpi=150)
    plt.close(fig)
    print(f'Saved amide_formation_scan.png/.pdf to {output_dir}')


def main() -> None:
    ap = argparse.ArgumentParser(description='Nylon-6,6 amide-formation PES + TDBB-reach scan')
    ap.add_argument('--backend', default='orb', choices=['orb', 'mace', 'toy'])
    ap.add_argument('--device', default='cuda', help='cuda (default) or cpu')
    ap.add_argument('--spin', type=int, default=1,
                    help='closed-shell singlet (default 1) — no radical spin complication')
    ap.add_argument('--monomer-fragments', action='store_true',
                    help='use fragments cut from the real HMDA/adipic-acid monomers '
                         'instead of methylamine/acetic acid')
    ap.add_argument('--r-min', type=float, default=1.4)
    ap.add_argument('--r-max', type=float, default=6.0)
    ap.add_argument('--n-points', type=int, default=24,
                    help='number of scan distances between r-min and r-max')
    ap.add_argument('--f2-sweep', type=float, nargs='+', default=[10.0, 5.0, 2.0],
                    help='f2 values for the analytic TDBB-reach table (paper default 10)')
    ap.add_argument('--lambda-vdw', type=float, default=0.60)
    ap.add_argument('--force-only', action='store_true',
                    help='skip the MLIP PES scan; emit only the analytic reach table')
    ap.add_argument('--output-dir', type=Path, default=Path('runs/calibration'))
    args = ap.parse_args()

    # ── forming pair r0 = λ·(r_vdw^N + r_vdw^C)  (Eq. 4) ──
    vdw_N = VDW_RADII['N']
    vdw_C = VDW_RADII['C']
    r0 = target_distance(np.array([vdw_N, vdw_C]), args.lambda_vdw)

    # ── analytic TDBB-reach table (always; no MLIP needed) ──
    f1 = TDBBParams().f1_max_formation
    shell_dists = [1.5, 1.95, 2.0, 2.5, 2.8]
    window_dists = [3.0, 3.2, 3.5, 4.0, 5.0, 6.0]
    reach = _reach_table(r0, args.f2_sweep, f1, shell_dists, window_dists)

    print('=== Nylon-6,6 amide formation de-risk scan ===')
    frag_label = 'HMDA/adipic-acid fragments' if args.monomer_fragments \
        else 'methylamine + acetic acid'
    print(f'Fragments: {frag_label}   (closed-shell, spin={args.spin})')
    print(f'Forming pair: amine N -- carboxyl C   r0 = {r0:.3f} Å '
          f'(λ={args.lambda_vdw}, r_vdw N={vdw_N}, C={vdw_C})')
    print('\nTDBB formation-bias reach (|dV^f/dr|, kcal/mol/Å; f1 saturated):')
    print(f'{"r(A)":>7} ' + ' '.join(f'f2={f2:g}'.rjust(12) for f2 in args.f2_sweep))
    all_dists = sorted(set(shell_dists) | set(window_dists))
    for i, r in enumerate(all_dists):
        cells = []
        for f2 in args.f2_sweep:
            cells.append(f'{reach["by_f2"][str(f2)]["table"][i]["force_kcal_mol_A"]:12.3f}')
        win = ' *' if 3.0 <= r <= 6.0 else ''
        print(f'{r:7.2f} ' + ' '.join(cells) + win)
    print('  (* = inside the [3,6] Å candidate window)')
    for f2 in args.f2_sweep:
        mx = reach['by_f2'][str(f2)]['max_force_in_window_kcal_mol_A']
        print(f'  max |force| in [3,6] window at f2={f2:g}: {mx:.3f} kcal/mol/Å')

    rows = []
    if not args.force_only:
        if args.monomer_fragments:
            amine_frag, acid_frag = _build_amide_fragments(_DIAMINE_SMILES, _DIACID_SMILES)
        else:
            amine_frag, acid_frag = _build_amide_fragments()

        if args.backend == 'orb':
            from kagome.backends.orb_backend import create_orb_calculator
            calc = create_orb_calculator(device=args.device, spin=args.spin)
        elif args.backend == 'mace':
            from kagome.backends.mace_backend import create_mace_calculator
            calc = create_mace_calculator(device=args.device)
        else:
            from kagome.backends.toy import ToyCalculator
            calc = ToyCalculator()

        radius_list = list(np.linspace(args.r_min, args.r_max, args.n_points))
        print()
        rows = _run_pes_scan(calc, amine_frag, acid_frag, radius_list)

    # ── derive PES numbers ──
    pes_summary: dict = {}
    if rows:
        de = [row[2] for row in rows]
        r_at_min = rows[int(np.argmin(de))][0]
        barrier = max(de)
        r_at_barrier = rows[int(np.argmax(de))][0]
        well_depth = min(de)
        # is there an accessible bonded adduct? (a minimum at short r below barrier)
        short = [row for row in rows if row[0] <= 2.2]
        has_bonded_min = bool(short and min(row[2] for row in short) < barrier)
        pes_summary = {
            'barrier_kcal_mol': round(barrier, 4),
            'r_at_barrier_A': round(r_at_barrier, 3),
            'well_depth_kcal_mol': round(well_depth, 4),
            'r_at_min_A': round(r_at_min, 3),
            'has_accessible_bonded_min': has_bonded_min,
        }
        print(f'\nPES: barrier {barrier:.2f} kcal/mol at r={r_at_barrier:.2f} Å; '
              f'min dE {well_depth:.2f} at r={r_at_min:.2f} Å; '
              f'bonded min accessible={has_bonded_min}')

    # ── verdict ──
    f2_paper = str(10.0) if 10.0 in args.f2_sweep else str(args.f2_sweep[0])
    max_win_paper = reach['by_f2'].get(f2_paper, {}).get(
        'max_force_in_window_kcal_mol_A', 0.0)
    dead_zone = max_win_paper < 1.0
    verdict_lines = []
    if dead_zone:
        verdict_lines.append(
            f'DEAD-ZONE at paper f2={f2_paper}: max bias force in [3,6] window is '
            f'{max_win_paper:.3f} kcal/mol/Å (~0) -> selected pairs are NOT pulled '
            'into the capture shell. Same failure mode as MA (decisions 2026-06-26). '
            'Lower f2 (see sweep) to bridge the dead-zone before long nylon runs.')
    else:
        verdict_lines.append(
            f'REACH OK at paper f2={f2_paper}: max bias force in [3,6] window is '
            f'{max_win_paper:.3f} kcal/mol/Å -> bias reaches selected pairs.')
    # which f2 first gives a usable (>~10 kcal/mol/Å) reach at r=3.5?
    for f2 in sorted(args.f2_sweep):
        tbl = {row['r_A']: row['force_kcal_mol_A'] for row in reach['by_f2'][str(f2)]['table']}
        if tbl.get(3.5, 0.0) >= 10.0:
            verdict_lines.append(f'First f2 with |force|>=10 at r=3.5 Å: f2={f2:g}.')
            break
    for line in verdict_lines:
        print('VERDICT:', line)

    # ── write JSON ──
    results = {
        'system': {
            'fragments': frag_label,
            'amine_smiles': _DIAMINE_SMILES if args.monomer_fragments else _AMINE_SMILES,
            'acid_smiles': _DIACID_SMILES if args.monomer_fragments else _ACID_SMILES,
            'forming_pair': 'amine_N -- carboxyl_C',
            'spin': args.spin,
            'chemistry_note': 'forms the tetrahedral C-N addition intermediate; '
                              'water/leaving-group elimination NOT modelled (de-risk scope)',
        },
        'tdbb': {
            'r0_A': round(r0, 4),
            'lambda_vdw': args.lambda_vdw,
            'vdw_N': vdw_N,
            'vdw_C': vdw_C,
            'f1_max_formation': f1,
            'f2_sweep': args.f2_sweep,
            'candidate_window_A': [3.0, 6.0],
        },
        'reach': reach,
        'pes': pes_summary,
        'pes_scan': [(round(r, 4), round(e, 4), round(d, 4)) for r, e, d in rows],
        'backend': None if args.force_only else args.backend,
        'device': None if args.force_only else args.device,
        'verdict': verdict_lines,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    out_json = args.output_dir / 'amide_formation_scan.json'
    with open(out_json, 'w') as f:
        json.dump(results, f, indent=2)
    print(f'\nSaved {out_json}')

    if not args.force_only:
        _plot(rows, reach, r0, args.output_dir, frag_label)


if __name__ == '__main__':
    main()
