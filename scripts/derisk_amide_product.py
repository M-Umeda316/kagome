"""PRODUCT-side de-risk for nylon amide condensation under OrbMol-v2.

Why this exists (distinct from ``scripts/scan_amide_formation.py``)
-------------------------------------------------------------------
The TDBB workflow never *confirms* a nylon amide bond because it only edits a
side-car bond graph and never removes the leaving-group atoms from the MLIP atom
set: OrbMol therefore always sees a crowded pre-reactive complex (amine still
protonated, carboxyl still bearing -OH) and relaxes the forming C-N bond back
apart.  Before committing to the large "atom-level reaction execution" change
(option A: at commit, delete the amine H + carboxyl OH, place a water molecule,
present the REAL amide product to OrbMol, then relax) we must de-risk its core
premise: *does OrbMol-v2 actually hold a genuine amide product as a stable
minimum?*

``scan_amide_formation.py`` does NOT answer this — it scans the reactant approach
/ tetrahedral intermediate with the leaving groups still present and explicitly
does not model elimination (see its docstring).  This script answers the
PRODUCT-side question with three checks:

Check 1 — PRODUCT STABILITY (the premise of option A):
    Build the real condensation product N-methylacetamide (CC(=O)NC) plus a
    separate water molecule a few Å away, then FREE-relax (unconstrained FIRE,
    no distance constraint).  PASS if the amide C-N stays bonded (~1.30-1.45 Å)
    and does not blow apart.

Check 2 — REACTANT-COMPLEX REVERSION (the workflow failure, minimally reproduced):
    Build the un-reacted complex acetic acid (CC(=O)O) + methylamine (CN) with
    the amine N at the carboxyl C at a bonded-ish addition distance and the
    leaving groups STILL present — exactly what the workflow relaxes — then
    FREE-relax.  EXPECTED: N···C reverts (final >> initial), confirming that
    without removing the leaving groups the amide C-N does not persist.

Check 3 — REACTION ENERGETICS (sanity):
    dE_reaction = E(product + water, far apart, relaxed)
                - E(reactants, far apart, relaxed).
    Mass-balanced (acetic acid C2H4O2 + methylamine CH5N  ->  NMA C3H7NO + H2O),
    so the number is meaningful.  Gas-phase amide formation is roughly
    thermoneutral-to-mildly-endothermic; we report the number and do NOT hard-fail
    on its sign.

UNITS (verified):
    * ``Calculator.compute`` returns energy in kcal/mol and forces in kcal/mol/Å
      (src/kagome/backends/base.py:23; the orb backend converts eV->kcal/mol via
      EV_TO_KCAL_MOL, orb_backend.py:199-200).
    * FIRE (``kagome.integrators.minimize.fire_minimize``) works entirely in those
      same units; its convergence threshold is ``FireParams.fmax_kcal_mol_A``
      (minimize.py:38, default 1.0 kcal/mol/Å ≈ 0.043 eV/Å ≈ ASE's 0.05 eV/Å).
    Every reported quantity below is therefore kcal/mol (energy) or kcal/mol/Å
    (force); distances are Å (plain Euclidean — small gas-phase clusters, cell=None,
    no PBC).

Usage:
    # real chemistry (orchestrator, GPU/WSL):
    python scripts/derisk_amide_product.py --backend orb --device cuda
    # cheap plumbing dry-run (no GPU, meaningless chemistry):
    python scripts/derisk_amide_product.py --backend toy --device cpu --fire-max-steps 20
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')

# Å / ··· appear in printed output; force UTF-8 so a cp932 (Windows) console
# does not raise UnicodeEncodeError. No-op where stdout is already UTF-8 (WSL).
try:
    sys.stdout.reconfigure(encoding='utf-8')  # type: ignore[union-attr]
except Exception:  # pragma: no cover - stream may not support reconfigure
    pass

# ── self-contained import bootstrap (repo root for ``scripts``/``kagome``) ──
_REPO_ROOT = Path(__file__).resolve().parents[1]
for _p in (_REPO_ROOT, _REPO_ROOT / 'src'):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import numpy as np

from scripts._systems import _rdkit_3d
from scripts.scan_amide_formation import _assemble_system, _build_amide_fragments
from kagome.backends.base import Calculator
from kagome.chem.builders import _rdkit_mol
from kagome.integrators.minimize import FireParams, fire_minimize

# Product/reactant SMILES (mass-balanced condensation, closed-shell singlets).
_NMA_SMILES = 'CC(=O)NC'   # N-methylacetamide  (the real condensation product)
_WATER_SMILES = 'O'        # the eliminated water
_ACID_SMILES = 'CC(=O)O'   # acetic acid   (carboxyl end)
_AMINE_SMILES = 'CN'       # methylamine   (amine end)

# Amide C-N is bonded if its length is in this window (Å); ~1.32-1.35 Å real.
_BONDED_CN_MIN_A = 1.20
_BONDED_CN_MAX_A = 1.60
# Check-2 reversion threshold on the N···C distance (Å).
_REVERT_NC_A = 2.5


# ── small geometry / RDKit helpers ──────────────────────────────────────────

def _dist(positions: np.ndarray, i: int, j: int) -> float:
    """Plain Euclidean distance (Å) between atoms i and j (no PBC)."""
    return float(np.linalg.norm(positions[i] - positions[j]))


def _find_amide_cn(mol) -> tuple[int, int]:
    """Return ``(carbonyl_C_idx, amide_N_idx)`` for the amide in an RDKit mol.

    Derivation (no hardcoded indices): the amide carbonyl carbon is the unique C
    that simultaneously (a) has a C=O double bond and (b) is bonded to an N; that
    N is the amide nitrogen.  Asserts exactly one such carbon so a wrong molecule
    or atom-ordering surprise fails loudly rather than silently mislabelling the
    bond.  Indices are valid for the positions/species obtained from the SAME
    ``_rdkit_mol`` object (AddHs ordering is deterministic per SMILES).
    """
    matches: list[tuple[int, int]] = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'C':
            continue
        has_carbonyl = False
        n_neighbor: int | None = None
        for nbr in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if nbr.GetSymbol() == 'O' and bond.GetBondTypeAsDouble() == 2.0:
                has_carbonyl = True
            elif nbr.GetSymbol() == 'N':
                n_neighbor = nbr.GetIdx()
        if has_carbonyl and n_neighbor is not None:
            matches.append((atom.GetIdx(), n_neighbor))
    assert len(matches) == 1, (
        f'expected exactly one amide (C=O)-N carbon, found {matches}'
    )
    return matches[0]


def _combine_far(
    frag_a: tuple[np.ndarray, list[str]],
    frag_b: tuple[np.ndarray, list[str]],
    gap_A: float = 4.0,
) -> tuple[np.ndarray, list[str]]:
    """Concatenate two (positions, species) fragments separated along +x.

    Fragment B is translated so its minimum-x atom sits ``gap_A`` beyond
    fragment A's maximum-x atom, guaranteeing every A-B interatomic distance is
    at least ``gap_A`` (non-overlapping).  A keeps its original indices; B's
    follow A's.
    """
    pa, sa = frag_a
    pb, sb = frag_b
    pa = np.asarray(pa, dtype=np.float64)
    pb = np.asarray(pb, dtype=np.float64)
    shift = pa[:, 0].max() - pb[:, 0].min() + gap_A
    pb_shifted = pb + np.array([shift, 0.0, 0.0])
    return np.vstack([pa, pb_shifted]), list(sa) + list(sb)


def _build_product_system(
    seed: int, gap_A: float = 4.0,
) -> tuple[np.ndarray, list[str], int, int]:
    """Build N-methylacetamide + a well-separated water molecule.

    Returns ``(positions, species, amide_C_global, amide_N_global)``.  The NMA is
    placed first so its amide C/N indices (derived from the same RDKit mol) are
    global indices unchanged by concatenation.
    """
    mol = _rdkit_mol(_NMA_SMILES, seed=seed)
    nma_pos = np.array(mol.GetConformer().GetPositions(), dtype=np.float64)
    nma_sp = [a.GetSymbol() for a in mol.GetAtoms()]
    c_idx, n_idx = _find_amide_cn(mol)

    water_pos, water_sp = _rdkit_3d(_WATER_SMILES, seed=seed + 1)
    positions, species = _combine_far((nma_pos, nma_sp), (water_pos, water_sp), gap_A)
    return positions, species, c_idx, n_idx


def _free_relax(
    positions: np.ndarray,
    species: list[str],
    calc: Calculator,
    fmax_kcal_mol_A: float,
    max_steps: int,
):
    """Unconstrained FIRE relaxation (cell=None, no distance constraint)."""
    params = FireParams(fmax_kcal_mol_A=fmax_kcal_mol_A, max_steps=max_steps)
    return fire_minimize(positions, species, None, calc, params)


# ── the three checks ────────────────────────────────────────────────────────

def check1_product_stability(
    calc: Calculator, seed: int, fmax: float, max_steps: int,
) -> dict:
    """Check 1: does OrbMol hold the real amide product as a stable minimum?

    Free-relaxes N-methylacetamide + water and measures the amide C-N bond.
    PASS iff the relaxed C-N is inside the bonded window and did not blow apart.
    Also returns the relaxed product energy (kcal/mol) for reuse in Check 3.
    """
    positions, species, c_idx, n_idx = _build_product_system(seed)
    cn_initial = _dist(positions, c_idx, n_idx)
    res = _free_relax(positions, species, calc, fmax, max_steps)
    cn_final = _dist(res.positions, c_idx, n_idx)
    stayed_bonded = _BONDED_CN_MIN_A <= cn_final <= _BONDED_CN_MAX_A
    return {
        'amide_C_idx': int(c_idx),
        'amide_N_idx': int(n_idx),
        'cn_initial_A': round(cn_initial, 4),
        'cn_final_A': round(cn_final, 4),
        'bonded_window_A': [_BONDED_CN_MIN_A, _BONDED_CN_MAX_A],
        'final_max_force_kcal_mol_A': round(res.fmax, 4),
        'n_steps': res.n_steps,
        'fire_converged': bool(res.converged),
        'product_energy_kcal_mol': round(res.energy, 4),
        'passed': bool(stayed_bonded),
    }


def check2_reactant_reversion(
    calc: Calculator, seed: int, reactant_r: float, fmax: float, max_steps: int,
) -> dict:
    """Check 2: does the leaving-groups-present complex revert on free relax?

    Reuses ``_build_amide_fragments`` / ``_assemble_system`` (from
    scan_amide_formation.py) to place methylamine's N at acetic acid's carboxyl C
    at ``reactant_r`` Å with the leaving groups still present — exactly what the
    TDBB workflow relaxes — then free-relaxes.  EXPECTED reverts=True: the amide
    C-N does not persist without removing the leaving groups.
    """
    amine_frag, acid_frag = _build_amide_fragments(_AMINE_SMILES, _ACID_SMILES)
    positions, species, c_glob, n_glob = _assemble_system(
        amine_frag, acid_frag, reactant_r,
    )
    nc_initial = _dist(positions, c_glob, n_glob)
    res = _free_relax(positions, species, calc, fmax, max_steps)
    nc_final = _dist(res.positions, c_glob, n_glob)
    reverts = nc_final > _REVERT_NC_A
    return {
        'carboxyl_C_idx': int(c_glob),
        'amine_N_idx': int(n_glob),
        'reactant_r_A': round(float(reactant_r), 4),
        'nc_initial_A': round(nc_initial, 4),
        'nc_final_A': round(nc_final, 4),
        'revert_threshold_A': _REVERT_NC_A,
        'final_max_force_kcal_mol_A': round(res.fmax, 4),
        'n_steps': res.n_steps,
        'fire_converged': bool(res.converged),
        'reverts': bool(reverts),
    }


def check3_reaction_energetics(
    calc: Calculator, seed: int, product_energy_kcal_mol: float,
    fmax: float, max_steps: int,
) -> dict:
    """Check 3: dE_reaction = E(product+water) - E(reactants), both far-apart/relaxed.

    Reuses Check 1's relaxed product energy for the product term (identical
    system) and relaxes the far-apart reactants here.  Sanity number only — not a
    hard-fail — since gas-phase amide formation is roughly thermoneutral-to-mildly
    -endothermic.
    """
    acid_pos, acid_sp = _rdkit_3d(_ACID_SMILES, seed=seed)
    amine_pos, amine_sp = _rdkit_3d(_AMINE_SMILES, seed=seed + 1)
    react_pos, react_sp = _combine_far(
        (acid_pos, acid_sp), (amine_pos, amine_sp), gap_A=6.0,
    )
    res = _free_relax(react_pos, react_sp, calc, fmax, max_steps)
    d_e = product_energy_kcal_mol - res.energy
    return {
        'reactants_energy_kcal_mol': round(res.energy, 4),
        'product_energy_kcal_mol': round(product_energy_kcal_mol, 4),
        'dE_reaction_kcal_mol': round(d_e, 4),
        'reactants_fire_converged': bool(res.converged),
        'note': 'gas-phase amide formation is ~thermoneutral-to-mildly-endothermic; '
                'sign is reported for sanity, not hard-failed',
    }


def _make_calculator(backend: str, device: str) -> Calculator:
    """Instantiate the requested backend (closed-shell singlet for orb)."""
    if backend == 'orb':
        from kagome.backends.orb_backend import create_orb_calculator
        return create_orb_calculator(device=device, spin=1)
    if backend == 'toy':
        from kagome.backends.toy import ToyCalculator
        return ToyCalculator()
    raise ValueError(f'unknown backend {backend!r}')


def main() -> None:
    ap = argparse.ArgumentParser(
        description='Product-side de-risk: does OrbMol-v2 hold a real amide product?',
    )
    ap.add_argument('--backend', default='orb', choices=['orb', 'toy'])
    ap.add_argument('--device', default='cuda', help='cuda (default) or cpu')
    ap.add_argument(
        '--fire-fmax', type=float, default=1.0,
        help='FIRE convergence threshold on the largest per-atom force, in '
             'kcal/mol/Å (the repo/FIRE native unit; default 1.0 ≈ 0.043 eV/Å '
             '≈ ASE 0.05 eV/Å — see minimize.py:34-38)',
    )
    ap.add_argument('--fire-max-steps', type=int, default=500)
    ap.add_argument('--reactant-r', type=float, default=1.55,
                    help='initial amine-N···carboxyl-C distance for Check 2 (Å)')
    ap.add_argument('--output-dir', type=Path, default=Path('runs/calibration'))
    ap.add_argument('--seed', type=int, default=20260708,
                    help='RDKit embedding seed (deterministic geometry)')
    args = ap.parse_args()

    print('=== Amide PRODUCT-side de-risk (option A premise) ===')
    print(f'backend={args.backend}  device={args.device}  seed={args.seed}')
    print(f'FIRE: fmax={args.fire_fmax} kcal/mol/Å, max_steps={args.fire_max_steps}')

    calc = _make_calculator(args.backend, args.device)

    # Check 1 (also yields the relaxed product energy reused by Check 3).
    c1 = check1_product_stability(calc, args.seed, args.fire_fmax, args.fire_max_steps)
    print('\n[Check 1] PRODUCT STABILITY (N-methylacetamide + water, free relax)')
    print(f'  amide C idx={c1["amide_C_idx"]}, N idx={c1["amide_N_idx"]}')
    print(f'  C-N: {c1["cn_initial_A"]:.3f} Å -> {c1["cn_final_A"]:.3f} Å  '
          f'(bonded window {_BONDED_CN_MIN_A}-{_BONDED_CN_MAX_A} Å)')
    print(f'  final max force = {c1["final_max_force_kcal_mol_A"]:.3f} kcal/mol/Å '
          f'({c1["n_steps"]} steps, converged={c1["fire_converged"]})')
    print(f'  -> Check 1 {"PASS" if c1["passed"] else "FAIL"}')

    # Check 2.
    c2 = check2_reactant_reversion(
        calc, args.seed, args.reactant_r, args.fire_fmax, args.fire_max_steps,
    )
    print('\n[Check 2] REACTANT-COMPLEX REVERSION (acetic acid + methylamine, '
          'leaving groups present)')
    print(f'  N···C: {c2["nc_initial_A"]:.3f} Å -> {c2["nc_final_A"]:.3f} Å  '
          f'(revert threshold {_REVERT_NC_A} Å)')
    print(f'  final max force = {c2["final_max_force_kcal_mol_A"]:.3f} kcal/mol/Å '
          f'({c2["n_steps"]} steps, converged={c2["fire_converged"]})')
    print(f'  -> reverts={c2["reverts"]} (True = reproduces the workflow failure)')

    # Check 3.
    c3 = check3_reaction_energetics(
        calc, args.seed, c1['product_energy_kcal_mol'],
        args.fire_fmax, args.fire_max_steps,
    )
    print('\n[Check 3] REACTION ENERGETICS (mass-balanced, far-apart relaxed)')
    print(f'  E(reactants) = {c3["reactants_energy_kcal_mol"]:.3f} kcal/mol')
    print(f'  E(product+water) = {c3["product_energy_kcal_mol"]:.3f} kcal/mol')
    print(f'  dE_reaction = {c3["dE_reaction_kcal_mol"]:.3f} kcal/mol')

    # ── overall verdict ──
    verdict: list[str] = []
    if c1['passed']:
        verdict.append(
            f'Check 1 PASS: OrbMol-{args.backend} HOLDS the amide product '
            f'(C-N {c1["cn_final_A"]:.3f} Å after free relax). Option A\'s premise '
            'HOLDS — presenting the real amide product (leaving groups removed, '
            'water placed) gives a stable minimum for OrbMol to relax into.')
    else:
        verdict.append(
            f'Check 1 FAIL: OrbMol-{args.backend} did NOT hold the amide product '
            f'(C-N went {c1["cn_initial_A"]:.3f} -> {c1["cn_final_A"]:.3f} Å). '
            'Option A\'s premise is NOT supported by this diagnostic — removing '
            'the leaving groups alone would not yield a stable amide; investigate '
            'before the architectural change.')
    if c2['reverts']:
        verdict.append(
            f'Check 2: the leaving-groups-present complex REVERTS '
            f'(N···C {c2["nc_initial_A"]:.3f} -> {c2["nc_final_A"]:.3f} Å), '
            'reproducing the workflow failure in a minimal system and confirming '
            'the diagnosis: the crowded pre-reactive complex, not OrbMol itself, '
            'is why amide bonds never confirm.')
    else:
        verdict.append(
            f'Check 2: the leaving-groups-present complex did NOT revert '
            f'(N···C -> {c2["nc_final_A"]:.3f} Å) — unexpected; the minimal system '
            'does not reproduce the workflow failure, re-examine assumptions.')
    verdict.append(
        f'Check 3: dE_reaction = {c3["dE_reaction_kcal_mol"]:.3f} kcal/mol '
        '(sanity only; ~thermoneutral-to-mildly-endothermic expected in gas phase).')
    overall_premise_holds = bool(c1['passed'])
    for line in verdict:
        print('\nVERDICT:', line)

    results = {
        'purpose': 'product-side de-risk: does OrbMol-v2 hold a real amide product '
                   'as a stable minimum (premise of atom-level reaction execution)',
        'backend': args.backend,
        'device': args.device,
        'seed': args.seed,
        'fire': {
            'fmax_kcal_mol_A': args.fire_fmax,
            'max_steps': args.fire_max_steps,
            'cell': None,
            'unit_note': 'energy kcal/mol, force kcal/mol/Å, distance Å (no PBC)',
        },
        'systems': {
            'product': f'{_NMA_SMILES} (N-methylacetamide) + {_WATER_SMILES} (water)',
            'reactant_complex': f'{_ACID_SMILES} (acetic acid) + {_AMINE_SMILES} '
                                '(methylamine), leaving groups present',
            'mass_balance': 'CC(=O)O + CN  ->  CC(=O)NC + H2O  (C3H9NO2 both sides)',
        },
        'check1_product_stability': c1,
        'check2_reactant_reversion': c2,
        'check3_reaction_energetics': c3,
        'overall_premise_holds': overall_premise_holds,
        'verdict': verdict,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    out_json = args.output_dir / 'amide_product_derisk.json'
    with open(out_json, 'w') as f:
        json.dump(results, f, indent=2)
    print(f'\nSaved {out_json}')
    print(f'OVERALL: option A premise (OrbMol stabilizes the amide product) '
          f'{"HOLDS" if overall_premise_holds else "DOES NOT HOLD"}.')


if __name__ == '__main__':
    main()
