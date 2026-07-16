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

Usage:
    python scripts/analyze_copolymer_reactivity.py --run-dir runs/copoly_smoke_orb
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from scripts._systems import (
    _METHACRYLATE_SMILES,
    _MONOMER_SMILES,
    copolymer_alpha_species,
)

_SMILES_NAME = {
    _MONOMER_SMILES: 'acrylate',
    _METHACRYLATE_SMILES: 'methacrylate',
}


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
    avail = {_MONOMER_SMILES: n_acr, _METHACRYLATE_SMILES: n_mac}

    bonds_path = run_dir / 'bonds.jsonl'
    incorporated: Counter[str] = Counter()
    per_cycle: dict[int, Counter[str]] = {}
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
        smi = alpha_species.get(ev['atom_a']) or alpha_species.get(ev['atom_b'])
        if smi is None:
            n_unmapped += 1
            continue
        incorporated[smi] += 1
        per_cycle.setdefault(ev['cycle'], Counter())[smi] += 1

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


if __name__ == '__main__':
    main()
