"""Relative reactivity of the two monomers in a vinyl copolymer run.

Reads a copolymer run directory and reports how many confirmed formations
incorporated each monomer species, normalized by that species' availability, so
acrylate vs methacrylate can be compared *relatively* (no absolute rate
constants). See specs/decisions.md 2026-07-16.

Mechanism: every ``confirmed_formation`` event is a radical_C–vinyl_alpha_C bond
(constraint-only pairs are never tracked). The alpha-C belongs to one monomer;
mapping it back to its species (via copolymer_alpha_species) tells us which
monomer was consumed. The relative reactivity index is

    R_species = (incorporations of species) / (initial count of species)

and the ratio R_acrylate / R_methacrylate says which monomer reacts faster per
available molecule. >1 means acrylate is preferentially incorporated.

Cross-propagation: the radical_C endpoint of each confirmed formation belongs
to a molecule block (initiator, or the monomer whose beta-C currently carries
the radical after chain propagation), so every event is also classified as
(terminal species -> incorporated species). The 2x2 monomer part of that table
estimates the Mayo-Lewis reactivity ratios

    r_acrylate      = (N_AA / N_AM) * (avail_M / avail_A)
    r_methacrylate  = (N_MM / N_MA) * (avail_A / avail_M)

using INITIAL availabilities (a low-conversion approximation; at equimolar feed
the availability factor is 1). Initiator-terminal events are tallied separately
and excluded from the ratios.

Usage:
    python scripts/analyze_copolymer_reactivity.py --run-dir runs/copoly_smoke_orb
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from scripts._systems import (
    _INITIATOR_SMILES,
    _METHACRYLATE_SMILES,
    _MONOMER_SMILES,
    copolymer_alpha_species,
    copolymer_atom_species,
)

_SMILES_NAME = {
    _MONOMER_SMILES: 'acrylate',
    _METHACRYLATE_SMILES: 'methacrylate',
}
_TERMINAL_NAME = {**_SMILES_NAME, _INITIATOR_SMILES: 'initiator'}


def analyze(run_dir: Path) -> dict:
    summary = json.loads((run_dir / 'summary.json').read_text(encoding='utf-8'))
    n_acr = summary['n_acrylate']
    n_mac = summary['n_methacrylate']
    n_init = summary['n_initiators']

    monomer_specs = [
        (_MONOMER_SMILES, n_acr),
        (_METHACRYLATE_SMILES, n_mac),
    ]
    alpha_species = copolymer_alpha_species(monomer_specs, n_init)
    atom_species = copolymer_atom_species(monomer_specs, n_init)
    avail = {_MONOMER_SMILES: n_acr, _METHACRYLATE_SMILES: n_mac}

    bonds_path = run_dir / 'bonds.jsonl'
    incorporated: Counter[str] = Counter()
    per_cycle: dict[int, Counter[str]] = {}
    cross: Counter[tuple[str, str]] = Counter()  # (terminal, incorporated)
    n_unmapped = 0
    for line in bonds_path.read_text(encoding='utf-8').splitlines():
        if not line.strip():
            continue
        ev = json.loads(line)
        if ev.get('event_type') != 'confirmed_formation':
            continue
        if not ev.get('counts_as_reaction', True):
            continue
        # The alpha-C is whichever endpoint is in the alpha->species map;
        # the other endpoint is the radical_C. Orientation-independent.
        a, b = ev['atom_a'], ev['atom_b']
        if a in alpha_species:
            alpha_idx, radical_idx = a, b
        elif b in alpha_species:
            alpha_idx, radical_idx = b, a
        else:
            n_unmapped += 1
            continue
        smi = alpha_species[alpha_idx]
        incorporated[smi] += 1
        per_cycle.setdefault(ev['cycle'], Counter())[smi] += 1
        # The radical's block is the chain's terminal unit at addition time
        # (initiator block before the first addition, else the monomer whose
        # beta-C received the radical in the last chain propagation).
        term_name = _TERMINAL_NAME.get(atom_species.get(radical_idx, ''), 'unknown')
        cross[(term_name, _SMILES_NAME[smi])] += 1

    total = sum(incorporated.values())
    result = {
        'n_acrylate': n_acr,
        'n_methacrylate': n_mac,
        'total_incorporations': total,
        'unmapped_events': n_unmapped,
        'by_species': {},
    }
    rel = {}
    for smi, name in _SMILES_NAME.items():
        c = incorporated.get(smi, 0)
        a = avail.get(smi, 0)
        r = (c / a) if a else float('nan')
        rel[name] = r
        result['by_species'][name] = {
            'incorporated': c,
            'available': a,
            'fraction_of_all_additions': (c / total) if total else float('nan'),
            'consumed_fraction': r,  # incorporations / initial availability
        }
    result['relative_reactivity'] = rel
    if rel.get('methacrylate'):
        result['acrylate_over_methacrylate'] = (
            rel['acrylate'] / rel['methacrylate'] if rel['methacrylate'] else float('inf')
        )
    result['per_cycle'] = {
        str(cyc): dict(cnt) for cyc, cnt in sorted(per_cycle.items())
    }

    terminals = ('initiator', 'acrylate', 'methacrylate')
    result['cross_propagation'] = {
        term: {
            mono: cross.get((term, mono), 0)
            for mono in ('acrylate', 'methacrylate')
        }
        for term in terminals
    }
    # Mayo-Lewis ratios from monomer-terminal events only, availability-
    # corrected with INITIAL counts (low-conversion approximation).
    n_aa = cross.get(('acrylate', 'acrylate'), 0)
    n_am = cross.get(('acrylate', 'methacrylate'), 0)
    n_mm = cross.get(('methacrylate', 'methacrylate'), 0)
    n_ma = cross.get(('methacrylate', 'acrylate'), 0)
    result['reactivity_ratio_estimates'] = {
        'r_acrylate': (n_aa / n_am) * (n_mac / n_acr) if n_am and n_acr else None,
        'r_methacrylate': (n_mm / n_ma) * (n_acr / n_mac) if n_ma and n_mac else None,
        'monomer_terminal_events': n_aa + n_am + n_mm + n_ma,
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, required=True)
    parser.add_argument('--json', action='store_true', help='emit raw JSON only.')
    args = parser.parse_args()

    result = analyze(args.run_dir)

    if args.json:
        print(json.dumps(result, indent=2))
        return

    print(f"=== Copolymer relative reactivity: {args.run_dir} ===")
    print(f"Composition: {result['n_acrylate']} acrylate + "
          f"{result['n_methacrylate']} methacrylate")
    print(f"Total confirmed incorporations: {result['total_incorporations']}")
    if result['unmapped_events']:
        print(f"  (warning: {result['unmapped_events']} formation events had no "
              f"alpha-C in the species map — unexpected)")
    print()
    print(f"{'species':<14}{'incorp':>8}{'avail':>8}{'%of adds':>10}{'consumed':>10}")
    for name in ('acrylate', 'methacrylate'):
        s = result['by_species'][name]
        pct = s['fraction_of_all_additions']
        cons = s['consumed_fraction']
        pct_s = '-' if pct != pct else f"{pct * 100:.0f}%"
        cons_s = '-' if cons != cons else f"{cons * 100:.0f}%"
        print(f"{name:<14}{s['incorporated']:>8}{s['available']:>8}{pct_s:>10}{cons_s:>10}")
    print()
    aom = result.get('acrylate_over_methacrylate')
    if aom is not None and result['total_incorporations'] > 0:
        if aom != aom:  # nan
            verdict = 'insufficient data'
        elif aom > 1.15:
            verdict = 'acrylate reacts MORE per available molecule'
        elif aom < 0.87:
            verdict = 'methacrylate reacts MORE per available molecule'
        else:
            verdict = 'comparable (within noise)'
        aom_s = 'inf' if aom == float('inf') else f"{aom:.2f}"
        print(f"Relative reactivity (acrylate / methacrylate, availability-normalized): "
              f"{aom_s}")
        print(f"  -> {verdict}")
        n = result['total_incorporations']
        if n < 10:
            print(f"  NOTE: only {n} events -- treat as a weak trend, not a "
                  f"statistically firm ratio. Rerun larger / multi-seed to confirm.")
    else:
        print("No confirmed incorporations yet — cannot compare reactivity.")

    if result['total_incorporations'] > 0:
        print()
        print('Cross-propagation (terminal unit -> incorporated monomer):')
        print(f"{'terminal':<14}{'-> acrylate':>12}{'-> methacrylate':>17}")
        for term in ('initiator', 'acrylate', 'methacrylate'):
            row = result['cross_propagation'][term]
            print(f"{term:<14}{row['acrylate']:>12}{row['methacrylate']:>17}")
        est = result['reactivity_ratio_estimates']
        n_term = est['monomer_terminal_events']
        print()
        r_a = est['r_acrylate']
        r_m = est['r_methacrylate']
        r_a_s = '-' if r_a is None else f'{r_a:.2f}'
        r_m_s = '-' if r_m is None else f'{r_m:.2f}'
        print(f"Mayo-Lewis estimates (initial-availability corrected, "
              f"{n_term} monomer-terminal events):")
        print(f"  r_acrylate = {r_a_s}   r_methacrylate = {r_m_s}")
        if n_term < 20:
            print(f"  NOTE: {n_term} monomer-terminal events -- ratios need "
                  f"O(10+) events PER TERMINAL to be meaningful.")


if __name__ == '__main__':
    main()
