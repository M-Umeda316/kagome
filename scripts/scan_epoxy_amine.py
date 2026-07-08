"""Epoxy-amine ring-opening de-risk for OrbMol-v2 + TDBB (Track 2 / E0 gate).

Why this exists
---------------
Before building a bulk epoxy-amine curing system we must de-risk, exactly as we
did for nylon, that OrbMol-v2 + TDBB can (a) drive the epoxide ring-opening +
primary-amine addition, and (b) HOLD the ring-opened product as a stable minimum
(no spontaneous ring re-closure).  The nylon de-risk had TWO parts that both
mattered:

  * ``scripts/scan_amide_formation.py``  — REACTANT-side PES + TDBB-reach /
    dead-zone scan (does the bias even reach the [3,6] Å selected pairs?).
  * ``scripts/derisk_amide_product.py``  — PRODUCT-side stability de-risk (the
    DECISIVE check for nylon: does OrbMol hold the real product as a minimum?).

This script combines both ideas for the epoxy-amine reaction into one diagnostic
with three sections + one JSON result + a printed VERDICT.

The reaction chemistry (bulk epoxy-amine ring-opening addition; ORGANIC only)
----------------------------------------------------------------------------
A primary amine N attacks the terminal (less-substituted) carbon of an epoxide
(oxirane) ring.  A new N-C bond forms; the ring C-O bond (from the attacked
carbon to the ring oxygen) breaks; one amine N-H breaks and that H transfers to
the ring oxygen to make a beta-hydroxyl.  NO small molecule leaves — the oxygen
stays in the molecule as an -OH on the adjacent ring carbon (INTRAMOLECULAR
ring-opening).  Example: methylamine (CN) + propylene oxide (CC1CO1) ->
1-(methylamino)propan-2-ol (CNCC(O)C).

Forming/breaking pair r0 (Eq. 4, lambda=0.60, VDW_RADII from
kagome.workflows.polymerization):
    N-C      = 0.60*(1.55+1.70) = 1.95 A   (the FORMING bond; TDBB formation bias)
    ring C-O = 0.60*(1.70+1.52) = 1.93 A   (breaks; dissociation bias)
    H-O      = 0.60*(1.10+1.52) = 1.57 A   (H transfer)
The TDBB reach / dead-zone analysis below is for the FORMING N-C bond, which uses
the SAME [3,6] A candidate window as nylon (build_nylon66_system) — so the
dead-zone risk is structurally identical and is measured, not assumed.

UNITS (verified, matching derisk_amide_product.py)
--------------------------------------------------
    * ``Calculator.compute`` returns energy in kcal/mol, forces in kcal/mol/A
      (src/kagome/backends/base.py:23).
    * FIRE (kagome.integrators.minimize) and constrained_relax work in those
      same units; distances are Å (plain Euclidean, cell=None, NON-periodic —
      small gas-phase clusters, no PBC).

Usage:
    # real chemistry (orchestrator, GPU/WSL):
    python scripts/scan_epoxy_amine.py --backend orb --device cuda
    # analytic reach only (no MLIP, fast):
    python scripts/scan_epoxy_amine.py --backend toy --device cpu --force-only
    # cheap plumbing dry-run (no GPU, meaningless chemistry):
    python scripts/scan_epoxy_amine.py --backend toy --device cpu \
        --fire-max-steps 20 --n-points 4
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')

# Å / ··· appear in printed output; force UTF-8 so a cp932 (Windows) console does
# not raise UnicodeEncodeError.  No-op where stdout is already UTF-8 (WSL).
try:
    sys.stdout.reconfigure(encoding='utf-8')  # type: ignore[union-attr]
except Exception:  # pragma: no cover - stream may not support reconfigure
    pass

# ── self-contained import bootstrap (repo root + /src for scripts/kagome) ──
_REPO_ROOT = Path(__file__).resolve().parents[1]
for _p in (_REPO_ROOT, _REPO_ROOT / 'src'):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import numpy as np

from scripts._systems import _find_amine_h, _find_terminal_amine_n, _rdkit_3d
from scripts.derisk_amide_product import _combine_far
from scripts.scan_amide_formation import _reach_table
from scripts.scan_radical_addition import _rotation_align, _unit, constrained_relax
from kagome.backends.base import Calculator
from kagome.boost.tdbb import (
    TDBBParams,
    formation_force_magnitude,  # noqa: F401 — reused indirectly via _reach_table
    formation_potential,        # noqa: F401 — reused indirectly via _reach_table
    target_distance,
)
from kagome.chem.builders import _rdkit_mol
from kagome.integrators.minimize import FireParams, fire_minimize
from kagome.workflows.polymerization import VDW_RADII

# Default fragments/product (methylamine + propylene oxide, closed-shell singlets).
_EPOXY_SMILES = 'CC1CO1'        # propylene oxide (oxirane); 'COCC1CO1' = glycidyl
                                # methyl ether is the DGEBA-glycidyl-ether-closer alt.
_AMINE_SMILES = 'CN'           # methylamine (primary amine)
_PRODUCT_SMILES = 'CNCC(O)C'   # 1-(methylamino)propan-2-ol (ring-opened product)

# Product N-C bonded window (Å); a real C-N single bond is ~1.47 Å.
_BONDED_NC_MIN_A = 1.35
_BONDED_NC_MAX_A = 1.60
# Ring re-closure threshold on the attacked-C···former-ring-O distance (Å).
# A C-O single bond is ~1.43 Å; in the OPEN product that C and O are 1,3 (through
# the beta-carbon), so ~2.4 Å apart.  Below this threshold the oxirane re-formed.
_RECLOSE_CO_A = 1.8
# Bond-presence threshold used for the reactant-complex control diagnostics (Å).
_BONDED_GENERIC_A = 1.8


# ── RDKit atom identification (derive, never hardcode; assert unambiguous) ────

def _formula(mol) -> dict[str, int]:
    """Element -> atom count for an explicit-H RDKit mol (from ``_rdkit_mol``)."""
    counts: dict[str, int] = {}
    for atom in mol.GetAtoms():
        counts[atom.GetSymbol()] = counts.get(atom.GetSymbol(), 0) + 1
    return counts


def _n_hydrogens(atom) -> int:
    """Number of explicit-H neighbours of an atom (mol built with AddHs)."""
    return sum(1 for n in atom.GetNeighbors() if n.GetSymbol() == 'H')


def _is_hydroxyl_o(atom) -> bool:
    """True iff ``atom`` is an -OH oxygen: O with exactly one C and one H neighbour.

    This is precisely 'a hydroxyl that does NOT bridge two carbons', so a True
    result confirms the epoxide ring is OPEN at this oxygen (not a C-O-C ether).
    """
    if atom.GetSymbol() != 'O':
        return False
    nbrs = list(atom.GetNeighbors())
    if len(nbrs) != 2:
        return False
    return sorted(n.GetSymbol() for n in nbrs) == ['C', 'H']


def _find_epoxide(mol) -> tuple[int, int]:
    """Return ``(Gj_terminal_C_idx, Gl_ring_O_idx)`` for the oxirane ring.

    Derivation (no hardcoded indices):
      * The oxirane is the unique 3-membered ring with exactly two C and one O.
        Asserts exactly one such ring so a wrong molecule fails loudly.
      * ``Gl`` = that ring O.
      * ``Gj`` = the LESS-substituted ring carbon = the ring carbon with MORE
        hydrogens (the terminal CH2 that the amine attacks).
      * Tiebreak (symmetric epoxide, e.g. ethylene oxide C1CO1 where both ring
        carbons are CH2): pick the lower atom index (deterministic).
    """
    ring_info = mol.GetRingInfo()
    oxiranes: list[tuple[int, ...]] = []
    for ring in ring_info.AtomRings():
        if len(ring) != 3:
            continue
        syms = [mol.GetAtomWithIdx(i).GetSymbol() for i in ring]
        if syms.count('C') == 2 and syms.count('O') == 1:
            oxiranes.append(ring)
    assert len(oxiranes) == 1, (
        f'expected exactly one oxirane (3-ring with 2 C + 1 O), found {oxiranes}'
    )
    ring = oxiranes[0]
    o_idx = next(i for i in ring if mol.GetAtomWithIdx(i).GetSymbol() == 'O')
    c_idxs = [i for i in ring if mol.GetAtomWithIdx(i).GetSymbol() == 'C']
    h0 = _n_hydrogens(mol.GetAtomWithIdx(c_idxs[0]))
    h1 = _n_hydrogens(mol.GetAtomWithIdx(c_idxs[1]))
    if h0 > h1:
        gj = c_idxs[0]
    elif h1 > h0:
        gj = c_idxs[1]
    else:
        gj = min(c_idxs)  # tiebreak: lower index (symmetric oxirane)
    return int(gj), int(o_idx)


def _find_product_nc(mol) -> tuple[int, int, int, int]:
    """Return ``(N_idx, attacked_C_idx, beta_C_idx, hydroxyl_O_idx)`` for the product.

    The formed N-C bond is the amine N bonded to a carbon (the former epoxy
    terminal carbon Gj) that carries the beta-hydroxyl motif: that carbon is
    bonded to ANOTHER carbon (the former ring carbon) which bears a hydroxyl O.
    Also confirms the former ring O is now a genuine hydroxyl (O bonded to exactly
    one C and one H — NOT bridging two carbons), i.e. the ring stayed open.

    Asserts exactly one unambiguous match.
    """
    matches: list[tuple[int, int, int, int]] = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'N':
            continue
        for c_nbr in atom.GetNeighbors():
            if c_nbr.GetSymbol() != 'C':
                continue
            for c2 in c_nbr.GetNeighbors():
                if c2.GetIdx() == atom.GetIdx() or c2.GetSymbol() != 'C':
                    continue
                oh = [o for o in c2.GetNeighbors() if _is_hydroxyl_o(o)]
                if oh:
                    matches.append(
                        (atom.GetIdx(), c_nbr.GetIdx(), c2.GetIdx(), oh[0].GetIdx())
                    )
    assert len(matches) == 1, (
        f'expected exactly one amine-N -- (beta-hydroxyl C) N-C bond, found {matches}'
    )
    return tuple(int(x) for x in matches[0])  # type: ignore[return-value]


# ── fragment / system builders ───────────────────────────────────────────────

def _build_epoxy_fragment(
    epoxy_smiles: str, seed: int,
) -> tuple[np.ndarray, list[str], int, int, np.ndarray]:
    """Return ``(pos, species, Gj_local, Gl_local, approach_dir)`` for the epoxide.

    ``approach_dir`` is the unit vector from the ring O (Gl) through the terminal
    carbon (Gj) — i.e. the amine is placed on the side of Gj AWAY from the ring
    oxygen (a backside-of-the-C-O-bond starting guess for nucleophilic attack).
    The constrained relaxation then relaxes everything else at each fixed r.
    """
    mol = _rdkit_mol(epoxy_smiles, seed)
    pos = np.array(mol.GetConformer().GetPositions(), dtype=np.float64)
    species = [a.GetSymbol() for a in mol.GetAtoms()]
    gj, gl = _find_epoxide(mol)
    approach = _unit(pos[gj] - pos[gl])
    return pos, species, gj, gl, approach


def _build_amine_fragment(
    amine_smiles: str, seed: int,
) -> tuple[np.ndarray, list[str], int, int, np.ndarray]:
    """Return ``(pos, species, N_local, H_local, lobe_dir)`` for the primary amine.

    ``lobe_dir`` points from the mean of N's bonded neighbours toward N (along the
    lone pair — the nucleophile's attack direction), mirroring
    scan_amide_formation._build_amide_fragments.
    """
    mol = _rdkit_mol(amine_smiles, seed)
    pos = np.array(mol.GetConformer().GetPositions(), dtype=np.float64)
    species = [a.GetSymbol() for a in mol.GetAtoms()]
    n_candidates = _find_terminal_amine_n(amine_smiles)
    if not n_candidates:
        raise ValueError(f'No primary amine N (>=2 H) in {amine_smiles!r}')
    n_local = n_candidates[0]
    h_local = _find_amine_h(amine_smiles, n_local)
    n_atom = mol.GetAtomWithIdx(n_local)
    nbr_idx = [nb.GetIdx() for nb in n_atom.GetNeighbors()]
    lobe = _unit(pos[n_local] - pos[nbr_idx].mean(axis=0))
    return pos, species, n_local, h_local, lobe


def _assemble_epoxy_system(
    epoxy_frag, amine_frag, r: float,
) -> tuple[np.ndarray, list[str], int, int]:
    """Place the amine at N···Gj distance ``r`` along the epoxide approach vector.

    Epoxide atoms come first, so its Gj/Gl indices are unchanged global indices.
    Returns ``(positions, species, Gj_global, N_global)``.
    """
    (ep, esp, gj_local, _gl_local, approach) = epoxy_frag
    (ap, asp, n_local, _h_local, lobe) = amine_frag
    # Orient the amine lone pair along -approach (pointing at the epoxy Gj carbon).
    rot = _rotation_align(lobe, -approach)
    ap_centered = ap - ap[n_local]
    ap_oriented = ap_centered @ rot.T
    amine_pos = ap_oriented + (ep[gj_local] + r * approach)
    species = list(esp) + list(asp)
    positions = np.vstack([ep, amine_pos])
    return positions, species, int(gj_local), int(len(esp) + n_local)


def _dist(positions: np.ndarray, i: int, j: int) -> float:
    """Plain Euclidean distance (Å) between atoms i and j (no PBC)."""
    return float(np.linalg.norm(positions[i] - positions[j]))


def _free_relax(positions, species, calc: Calculator, fmax: float, max_steps: int):
    """Unconstrained FIRE relaxation (cell=None, no distance constraint)."""
    params = FireParams(fmax_kcal_mol_A=fmax, max_steps=max_steps)
    return fire_minimize(positions, species, None, calc, params)


def _make_calculator(backend: str, device: str) -> Calculator:
    """Instantiate the requested backend (closed-shell singlet for orb)."""
    if backend == 'orb':
        from kagome.backends.orb_backend import create_orb_calculator
        return create_orb_calculator(device=device, spin=1)
    if backend == 'toy':
        from kagome.backends.toy import ToyCalculator
        return ToyCalculator()
    raise ValueError(f'unknown backend {backend!r}')


# ── Section 1: REACTANT-side PES scan (mirror scan_amide_formation) ──────────

def section1_pes_scan(
    calc: Calculator, epoxy_frag, amine_frag, radius_list: list[float],
) -> tuple[list[tuple[float, float, float]], dict]:
    """Constrained-relaxed E vs the forming N···Gj distance.  Reference = largest r.

    At each fixed r the REST of the system is relaxed (constrained_relax projects
    out the along-bond force so the optimizer does not fight the constraint).  The
    epoxide ring and the amine N-H stay intact in this rigid pair scan — like the
    amide tetrahedral-intermediate scan, this only asks 'can OrbMol+TDBB drag N
    and the epoxy terminal C to a bonded distance at all', not the full ring-open
    thermochemistry (that is Sections 2-3).

    Returns ``(rows, summary)`` where rows = ``[(r, E, dE_vs_reactants), ...]``.
    """
    print('  backend={}  (constrained-relaxed epoxy N...C scan)'.format(calc.name))
    print(f'  {"r_NC(A)":>8} {"E(kcal/mol)":>16} {"dE_vs_react":>14}')
    ordered = sorted(radius_list, reverse=True)  # far -> bonded; e_ref = separated
    e_ref: float | None = None
    rows: list[tuple[float, float, float]] = []
    for r in ordered:
        pos, species, gj_glob, n_glob = _assemble_epoxy_system(epoxy_frag, amine_frag, r)
        energy, _ = constrained_relax(pos, species, calc, gj_glob, n_glob, r)
        if e_ref is None:
            e_ref = energy
        rows.append((r, energy, energy - e_ref))
        print(f'  {r:8.2f} {energy:16.2f} {energy - e_ref:14.2f}', flush=True)
    rows.sort(key=lambda x: x[0])

    de = [row[2] for row in rows]
    barrier = max(de)
    r_at_barrier = rows[int(np.argmax(de))][0]
    well_depth = min(de)
    r_at_min = rows[int(np.argmin(de))][0]
    short = [row for row in rows if row[0] <= 2.2]
    has_bonded_min = bool(short and min(row[2] for row in short) < barrier)
    summary = {
        'barrier_kcal_mol': round(barrier, 4),
        'r_at_barrier_A': round(r_at_barrier, 3),
        'well_depth_kcal_mol': round(well_depth, 4),
        'r_at_min_A': round(r_at_min, 3),
        'has_accessible_bonded_min': has_bonded_min,
    }
    print(f'  PES: barrier {barrier:.2f} kcal/mol at r={r_at_barrier:.2f} A; '
          f'min dE {well_depth:.2f} at r={r_at_min:.2f} A; '
          f'bonded N-C min accessible={has_bonded_min}')
    return rows, summary


# ── Section 2: PRODUCT-side stability de-risk (mirror derisk_amide_product) ──

def section2_product_stability(
    calc: Calculator, product_smiles: str, seed: int, fmax: float, max_steps: int,
) -> dict:
    """Does OrbMol HOLD the ring-opened product without re-closing the epoxide?

    Free-relaxes the product; PASS iff the N-C stays bonded (1.35-1.60 Å) AND the
    former ring O remains a hydroxyl that does NOT re-close the ring (the attacked
    carbon stays > ``_RECLOSE_CO_A`` from that O — no C-C-O 3-ring reforms).
    Also returns the relaxed product energy (kcal/mol) for reuse in Section 3.
    """
    mol = _rdkit_mol(product_smiles, seed)
    pos = np.array(mol.GetConformer().GetPositions(), dtype=np.float64)
    species = [a.GetSymbol() for a in mol.GetAtoms()]
    n_idx, c_idx, beta_idx, o_idx = _find_product_nc(mol)

    nc_initial = _dist(pos, n_idx, c_idx)
    co_initial = _dist(pos, c_idx, o_idx)
    res = _free_relax(pos, species, calc, fmax, max_steps)
    nc_final = _dist(res.positions, n_idx, c_idx)
    co_final = _dist(res.positions, c_idx, o_idx)

    stayed_bonded = _BONDED_NC_MIN_A <= nc_final <= _BONDED_NC_MAX_A
    no_reclosure = co_final > _RECLOSE_CO_A
    passed = bool(stayed_bonded and no_reclosure)
    return {
        'N_idx': n_idx, 'attacked_C_idx': c_idx,
        'beta_C_idx': beta_idx, 'hydroxyl_O_idx': o_idx,
        'nc_initial_A': round(nc_initial, 4), 'nc_final_A': round(nc_final, 4),
        'nc_bonded_window_A': [_BONDED_NC_MIN_A, _BONDED_NC_MAX_A],
        'nc_stayed_bonded': bool(stayed_bonded),
        'attackedC_O_initial_A': round(co_initial, 4),
        'attackedC_O_final_A': round(co_final, 4),
        'reclosure_threshold_A': _RECLOSE_CO_A,
        'no_ring_reclosure': bool(no_reclosure),
        'final_max_force_kcal_mol_A': round(res.fmax, 4),
        'n_steps': res.n_steps,
        'fire_converged': bool(res.converged),
        'product_energy_kcal_mol': round(res.energy, 4),
        'passed': passed,
    }


def section2_reactant_control(
    calc: Calculator, epoxy_frag, amine_frag, reactant_r: float,
    fmax: float, max_steps: int,
) -> dict:
    """Diagnostic control (not a hard pass/fail): free-relax the epoxide + amine
    with N placed at the terminal carbon (~reactant_r Å), ring INTACT.

    Reports whether it reverts (amine leaves, N···C grows) and/or spontaneously
    ring-opens (the epoxide ring C-O, Gj-Gl, stretches past a bonded threshold).
    """
    (_ep, _esp, gj_local, gl_local, _ap) = epoxy_frag
    positions, species, gj_glob, n_glob = _assemble_epoxy_system(
        epoxy_frag, amine_frag, reactant_r,
    )
    gl_glob = gl_local  # epoxide placed first -> unchanged global index
    nc_initial = _dist(positions, n_glob, gj_glob)
    ring_co_initial = _dist(positions, gj_glob, gl_glob)
    res = _free_relax(positions, species, calc, fmax, max_steps)
    nc_final = _dist(res.positions, n_glob, gj_glob)
    ring_co_final = _dist(res.positions, gj_glob, gl_glob)
    return {
        'reactant_r_A': round(float(reactant_r), 4),
        'nc_initial_A': round(nc_initial, 4), 'nc_final_A': round(nc_final, 4),
        'ring_CO_initial_A': round(ring_co_initial, 4),
        'ring_CO_final_A': round(ring_co_final, 4),
        'amine_reverted': bool(nc_final > _BONDED_GENERIC_A),
        'ring_opened': bool(ring_co_final > _BONDED_GENERIC_A),
        'final_max_force_kcal_mol_A': round(res.fmax, 4),
        'n_steps': res.n_steps,
        'fire_converged': bool(res.converged),
        'note': 'informative only; ring intact at start (no bond edits in a plain relax)',
    }


# ── Section 3: reaction energetics ───────────────────────────────────────────

def section3_energetics(
    calc: Calculator, epoxy_smiles: str, amine_smiles: str,
    product_energy_kcal_mol: float, seed: int, fmax: float, max_steps: int,
) -> dict:
    """dE_reaction = E(product, relaxed) - E(epoxide + amine, far apart, relaxed).

    Mass-balanced (epoxide + amine -> ring-opened product; no atom added/removed).
    Epoxide ring-opening by an amine is exothermic (ring-strain release); expect
    roughly -15 to -30 kcal/mol.  Sign is reported for sanity, NOT hard-failed.
    """
    epoxy_pos, epoxy_sp = _rdkit_3d(epoxy_smiles, seed=seed)
    amine_pos, amine_sp = _rdkit_3d(amine_smiles, seed=seed + 1)
    react_pos, react_sp = _combine_far(
        (epoxy_pos, epoxy_sp), (amine_pos, amine_sp), gap_A=6.0,
    )
    res = _free_relax(react_pos, react_sp, calc, fmax, max_steps)
    d_e = product_energy_kcal_mol - res.energy
    return {
        'reactants_energy_kcal_mol': round(res.energy, 4),
        'product_energy_kcal_mol': round(product_energy_kcal_mol, 4),
        'dE_reaction_kcal_mol': round(d_e, 4),
        'reactants_fire_converged': bool(res.converged),
        'exothermic': bool(d_e < 0.0),
        'note': 'epoxy-amine ring-opening is exothermic (ring-strain release); '
                'expect ~ -15 to -30 kcal/mol. Sign reported for sanity, not hard-failed',
    }


# ── atom/mass balance validation ─────────────────────────────────────────────

def _validate_atom_balance(
    epoxy_mol, amine_mol, product_mol,
    epoxy_smiles: str, amine_smiles: str, product_smiles: str,
) -> dict:
    """Assert epoxide + amine formula == product formula (intramolecular ring-open).

    Errors clearly if the caller changed the epoxy/amine but not --product-smiles
    (or vice versa) so an inconsistent product is caught before any MLIP work.
    """
    lhs = _formula(epoxy_mol)
    for el, n in _formula(amine_mol).items():
        lhs[el] = lhs.get(el, 0) + n
    rhs = _formula(product_mol)
    if lhs != rhs:
        raise ValueError(
            'ATOM/MASS BALANCE FAILED: epoxide + amine != product '
            '(epoxy-amine ring-opening adds/removes NO atoms). '
            f'epoxide({epoxy_smiles!r}) + amine({amine_smiles!r}) = {dict(sorted(lhs.items()))} '
            f'but product({product_smiles!r}) = {dict(sorted(rhs.items()))}. '
            'If you changed --epoxy-smiles/--amine-smiles you MUST pass a '
            'consistent --product-smiles.'
        )
    return {'combined_reactant_formula': dict(sorted(lhs.items())),
            'product_formula': dict(sorted(rhs.items())), 'balanced': True}


def main() -> None:
    ap = argparse.ArgumentParser(
        description='Epoxy-amine ring-opening de-risk (Track 2 / E0 gate): '
                    'reactant PES + TDBB reach, product stability, reaction energetics',
    )
    ap.add_argument('--backend', default='orb', choices=['orb', 'toy'])
    ap.add_argument('--device', default='cuda', help='cuda (default) or cpu')
    ap.add_argument('--epoxy-smiles', default=_EPOXY_SMILES,
                    help="epoxide (default 'CC1CO1' propylene oxide; "
                         "'COCC1CO1' glycidyl methyl ether = DGEBA-closer alternative)")
    ap.add_argument('--amine-smiles', default=_AMINE_SMILES,
                    help="primary amine (default 'CN' methylamine)")
    ap.add_argument('--product-smiles', default=_PRODUCT_SMILES,
                    help="ring-opened product (default 'CNCC(O)C'); MUST change "
                         'consistently if epoxy/amine change (atom balance validated)')
    ap.add_argument('--fire-fmax', type=float, default=0.5,
                    help='FIRE convergence threshold on the largest per-atom force '
                         '(kcal/mol/A; repo/FIRE native unit)')
    ap.add_argument('--fire-max-steps', type=int, default=500)
    ap.add_argument('--r-min', type=float, default=1.4)
    ap.add_argument('--r-max', type=float, default=6.0)
    ap.add_argument('--n-points', type=int, default=24,
                    help='number of PES scan distances between r-min and r-max')
    ap.add_argument('--f2-sweep', type=float, nargs='+', default=[10.0, 5.0, 2.0],
                    help='f2 values for the analytic TDBB-reach table (paper default 10)')
    ap.add_argument('--lambda-vdw', type=float, default=0.60)
    ap.add_argument('--reactant-r', type=float, default=1.55,
                    help='initial amine-N···epoxy-terminal-C distance for the control (Å)')
    ap.add_argument('--force-only', action='store_true',
                    help='skip the MLIP parts; emit only the analytic reach table + RDKit IDs')
    ap.add_argument('--output-dir', type=Path, default=Path('runs/calibration'))
    ap.add_argument('--seed', type=int, default=20260708,
                    help='RDKit embedding seed (deterministic geometry)')
    args = ap.parse_args()

    print('=== Epoxy-amine ring-opening de-risk (Track 2 / E0 gate) ===')
    print(f'backend={args.backend}  device={args.device}  seed={args.seed}')
    print(f'epoxide={args.epoxy_smiles!r}  amine={args.amine_smiles!r}  '
          f'product={args.product_smiles!r}')

    # ── RDKit atom identification + atom balance (always; pure, no MLIP) ──
    epoxy_mol = _rdkit_mol(args.epoxy_smiles, args.seed)
    amine_mol = _rdkit_mol(args.amine_smiles, args.seed + 1)
    product_mol = _rdkit_mol(args.product_smiles, args.seed + 2)
    gj, gl = _find_epoxide(epoxy_mol)
    amine_ns = _find_terminal_amine_n(args.amine_smiles)
    if not amine_ns:
        raise ValueError(f'No primary amine N (>=2 H) in {args.amine_smiles!r}')
    gi = amine_ns[0]
    gk = _find_amine_h(args.amine_smiles, gi)
    p_n, p_c, p_beta, p_o = _find_product_nc(product_mol)
    balance = _validate_atom_balance(
        epoxy_mol, amine_mol, product_mol,
        args.epoxy_smiles, args.amine_smiles, args.product_smiles,
    )
    print('\nRDKit atom IDs:')
    print(f'  epoxy: Gj(terminal C)={gj}, Gl(ring O)={gl}')
    print(f'  amine: Gi(N)={gi}, Gk(N-H)={gk}')
    print(f'  product N-C: N={p_n}, attacked C={p_c}, beta C={p_beta}, hydroxyl O={p_o}')
    print(f'  atom balance: {balance["combined_reactant_formula"]} == '
          f'{balance["product_formula"]}  (balanced)')

    # ── forming/breaking r0 (Eq. 4) ──
    r0_nc = target_distance(np.array([VDW_RADII['N'], VDW_RADII['C']]), args.lambda_vdw)
    r0_co = target_distance(np.array([VDW_RADII['C'], VDW_RADII['O']]), args.lambda_vdw)
    r0_ho = target_distance(np.array([VDW_RADII['H'], VDW_RADII['O']]), args.lambda_vdw)

    # ── analytic TDBB-reach table for the FORMING N-C bond (always; no MLIP) ──
    f1 = TDBBParams().f1_max_formation
    shell_dists = [1.5, 1.95, 2.0, 2.5, 2.8]
    window_dists = [3.0, 3.2, 3.5, 4.0, 5.0, 6.0]
    reach = _reach_table(r0_nc, args.f2_sweep, f1, shell_dists, window_dists)

    print(f'\nForming pair: amine N -- epoxy terminal C   r0(N-C)={r0_nc:.3f} A '
          f'(lambda={args.lambda_vdw}); breaking ring C-O r0={r0_co:.3f} A, '
          f'H-O r0={r0_ho:.3f} A')
    print('TDBB formation-bias reach (|dV^f/dr|, kcal/mol/A; f1 saturated):')
    all_dists = sorted(set(shell_dists) | set(window_dists))
    print(f'  {"r(A)":>7} ' + ' '.join(f'f2={f2:g}'.rjust(12) for f2 in args.f2_sweep))
    for i, r in enumerate(all_dists):
        cells = [f'{reach["by_f2"][str(f2)]["table"][i]["force_kcal_mol_A"]:12.3f}'
                 for f2 in args.f2_sweep]
        win = ' *' if 3.0 <= r <= 6.0 else ''
        print(f'  {r:7.2f} ' + ' '.join(cells) + win)
    print('  (* = inside the [3,6] A candidate window)')
    for f2 in args.f2_sweep:
        mx = reach['by_f2'][str(f2)]['max_force_in_window_kcal_mol_A']
        print(f'  max |force| in [3,6] window at f2={f2:g}: {mx:.3f} kcal/mol/A')

    # dead-zone verdict + recommended f2
    f2_paper = str(10.0) if 10.0 in args.f2_sweep else str(args.f2_sweep[0])
    max_win_paper = reach['by_f2'].get(f2_paper, {}).get(
        'max_force_in_window_kcal_mol_A', 0.0)
    dead_zone = max_win_paper < 1.0
    first_usable_f2: float | None = None
    for f2 in sorted(args.f2_sweep):
        tbl = {row['r_A']: row['force_kcal_mol_A']
               for row in reach['by_f2'][str(f2)]['table']}
        if tbl.get(3.5, 0.0) >= 10.0:
            first_usable_f2 = f2
            break

    reach_verdict: list[str] = []
    if dead_zone:
        reach_verdict.append(
            f'DEAD-ZONE at paper f2={f2_paper}: max bias force in [3,6] window is '
            f'{max_win_paper:.3f} kcal/mol/A (~0) -> selected pairs are NOT pulled '
            'into the capture shell (same failure mode as MA/nylon). Lower f2.')
    else:
        reach_verdict.append(
            f'REACH OK at paper f2={f2_paper}: max bias force in [3,6] window is '
            f'{max_win_paper:.3f} kcal/mol/A -> bias reaches selected pairs.')
    if first_usable_f2 is not None:
        reach_verdict.append(
            f'First f2 with |force|>=10 at r=3.5 A: f2={first_usable_f2:g} '
            '(recommended f2 to bridge the dead-zone).')
    for line in reach_verdict:
        print('VERDICT(reach):', line)

    # ── MLIP sections (skipped under --force-only) ──
    s1: dict = {}
    s2: dict = {}
    s2_control: dict = {}
    s3: dict = {}
    pes_rows: list[tuple[float, float, float]] = []
    if not args.force_only:
        calc = _make_calculator(args.backend, args.device)
        print(f'\nFIRE: fmax={args.fire_fmax} kcal/mol/A, max_steps={args.fire_max_steps}')

        # Section 1
        epoxy_frag = _build_epoxy_fragment(args.epoxy_smiles, args.seed)
        amine_frag = _build_amine_fragment(args.amine_smiles, args.seed + 1)
        radius_list = list(np.linspace(args.r_min, args.r_max, args.n_points))
        print('\n[Section 1] REACTANT-side PES (constrained-relaxed N...C scan)')
        pes_rows, s1 = section1_pes_scan(calc, epoxy_frag, amine_frag, radius_list)

        # Section 2
        print('\n[Section 2] PRODUCT-side stability de-risk (free relax)')
        s2 = section2_product_stability(
            calc, args.product_smiles, args.seed, args.fire_fmax, args.fire_max_steps)
        print(f'  N-C: {s2["nc_initial_A"]:.3f} -> {s2["nc_final_A"]:.3f} A '
              f'(bonded window {_BONDED_NC_MIN_A}-{_BONDED_NC_MAX_A} A) '
              f'stayed_bonded={s2["nc_stayed_bonded"]}')
        print(f'  attacked-C...ring-O: {s2["attackedC_O_initial_A"]:.3f} -> '
              f'{s2["attackedC_O_final_A"]:.3f} A (reclose thresh {_RECLOSE_CO_A} A) '
              f'no_reclosure={s2["no_ring_reclosure"]}')
        print(f'  final max force = {s2["final_max_force_kcal_mol_A"]:.3f} kcal/mol/A '
              f'({s2["n_steps"]} steps, converged={s2["fire_converged"]})')
        print(f'  -> Section 2 {"PASS" if s2["passed"] else "FAIL"}')

        print('\n[Section 2b] REACTANT-COMPLEX control (ring intact, free relax)')
        s2_control = section2_reactant_control(
            calc, epoxy_frag, amine_frag, args.reactant_r,
            args.fire_fmax, args.fire_max_steps)
        print(f'  N...C: {s2_control["nc_initial_A"]:.3f} -> '
              f'{s2_control["nc_final_A"]:.3f} A  amine_reverted={s2_control["amine_reverted"]}')
        print(f'  ring C-O: {s2_control["ring_CO_initial_A"]:.3f} -> '
              f'{s2_control["ring_CO_final_A"]:.3f} A  ring_opened={s2_control["ring_opened"]}')

        # Section 3
        print('\n[Section 3] REACTION ENERGETICS (mass-balanced, far-apart relaxed)')
        s3 = section3_energetics(
            calc, args.epoxy_smiles, args.amine_smiles,
            s2['product_energy_kcal_mol'], args.seed, args.fire_fmax, args.fire_max_steps)
        print(f'  E(reactants) = {s3["reactants_energy_kcal_mol"]:.3f} kcal/mol')
        print(f'  E(product)   = {s3["product_energy_kcal_mol"]:.3f} kcal/mol')
        print(f'  dE_reaction  = {s3["dE_reaction_kcal_mol"]:.3f} kcal/mol '
              f'(exothermic={s3["exothermic"]})')

    # ── overall verdict + E0 gate ──
    overall: list[str] = list(reach_verdict)
    channel_ok: bool | None = None
    product_ok: bool | None = None
    exo_ok: bool | None = None
    gate: str
    if args.force_only:
        gate = 'INCOMPLETE (--force-only): reach table only; run with an MLIP ' \
               'backend to evaluate the product-stability + energetics gate.'
        overall.append(gate)
    else:
        channel_ok = bool(s1.get('has_accessible_bonded_min'))
        product_ok = bool(s2.get('passed'))
        exo_ok = bool(s3.get('exothermic'))
        overall.append(
            f'Section 1: accessible N-C channel exists = {channel_ok} '
            f'(barrier {s1.get("barrier_kcal_mol")} kcal/mol, '
            f'well {s1.get("well_depth_kcal_mol")} kcal/mol).')
        if product_ok:
            overall.append(
                f'Section 2: OrbMol HOLDS the ring-opened product WITHOUT re-closure '
                f'(N-C {s2["nc_final_A"]:.3f} A, attacked-C...O {s2["attackedC_O_final_A"]:.3f} A) '
                '-> the E0 gate premise holds.')
        else:
            overall.append(
                f'Section 2: OrbMol did NOT hold the ring-opened product '
                f'(N-C {s2["nc_final_A"]:.3f} A, attacked-C...O {s2["attackedC_O_final_A"]:.3f} A, '
                f'stayed_bonded={s2["nc_stayed_bonded"]}, '
                f'no_reclosure={s2["no_ring_reclosure"]}) -> E0 premise NOT supported.')
        overall.append(
            f'Section 3: dE_reaction = {s3["dE_reaction_kcal_mol"]:.3f} kcal/mol '
            f'(exothermic={exo_ok}).')
        proceed = bool(channel_ok and product_ok and exo_ok)
        if proceed:
            gate = 'PROCEED to E1: accessible N-C channel AND product held without ' \
                   're-closure AND exothermic dE.'
        else:
            failing = []
            if not channel_ok:
                failing.append('no accessible N-C channel (Section 1)')
            if not product_ok:
                failing.append('product not held / ring re-closed (Section 2)')
            if not exo_ok:
                failing.append('dE not exothermic (Section 3)')
            gate = 'DO NOT PROCEED to E1 — failing: ' + '; '.join(failing) + '.'
        overall.append(gate)
    for line in overall[len(reach_verdict):]:
        print('VERDICT:', line)

    # ── write JSON ──
    results = {
        'purpose': 'epoxy-amine ring-opening de-risk (Track 2 / E0 gate): can '
                   'OrbMol-v2 + TDBB drive the ring-opening + amine addition and '
                   'HOLD the ring-opened product without re-closure',
        'backend': None if args.force_only else args.backend,
        'device': None if args.force_only else args.device,
        'seed': args.seed,
        'system': {
            'epoxy_smiles': args.epoxy_smiles,
            'amine_smiles': args.amine_smiles,
            'product_smiles': args.product_smiles,
            'reaction': 'epoxide + primary amine -> beta-hydroxyl amine '
                        '(intramolecular ring-opening; NO small molecule leaves)',
            'atom_ids': {
                'epoxy_Gj_terminal_C': gj, 'epoxy_Gl_ring_O': gl,
                'amine_Gi_N': gi, 'amine_Gk_H': gk,
                'product_N': p_n, 'product_attacked_C': p_c,
                'product_beta_C': p_beta, 'product_hydroxyl_O': p_o,
            },
            'atom_balance': balance,
        },
        'tdbb': {
            'lambda_vdw': args.lambda_vdw,
            'r0_NC_A': round(r0_nc, 4), 'r0_ring_CO_A': round(r0_co, 4),
            'r0_HO_A': round(r0_ho, 4),
            'vdw_N': VDW_RADII['N'], 'vdw_C': VDW_RADII['C'],
            'vdw_O': VDW_RADII['O'], 'vdw_H': VDW_RADII['H'],
            'f1_max_formation': f1, 'f2_sweep': args.f2_sweep,
            'candidate_window_A': [3.0, 6.0],
        },
        'fire': {
            'fmax_kcal_mol_A': args.fire_fmax, 'max_steps': args.fire_max_steps,
            'cell': None,
            'unit_note': 'energy kcal/mol, force kcal/mol/A, distance A (no PBC)',
        },
        'reach': reach,
        'reach_verdict': reach_verdict,
        'dead_zone_at_paper_f2': dead_zone,
        'recommended_f2': first_usable_f2,
        'section1_reactant_pes': s1,
        'pes_scan': [(round(r, 4), round(e, 4), round(d, 4)) for r, e, d in pes_rows],
        'section2_product_stability': s2,
        'section2b_reactant_control': s2_control,
        'section3_energetics': s3,
        'gate': {
            'accessible_NC_channel': channel_ok,
            'product_held_no_reclosure': product_ok,
            'dE_exothermic': exo_ok,
            'decision': gate,
        },
        'verdict': overall,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    out_json = args.output_dir / 'epoxy_amine_scan.json'
    with open(out_json, 'w') as f:
        json.dump(results, f, indent=2)
    print(f'\nSaved {out_json}')


if __name__ == '__main__':
    main()
