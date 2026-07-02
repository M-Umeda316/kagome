"""Reconstruct a bond-topology history for a run made BEFORE topology tracking.

Pre-2026-07-02 runs wrote coordinates + bonds.jsonl (confirmed formations) but
no topology.jsonl.  Because the system build is deterministic and every reaction
is recorded, the connectivity can be replayed exactly: rebuild the initial
topology, remove the activation-dissociated azo C-N bonds, then apply each
confirmed formation (radical + vinyl_alpha_C -> C-C, C=C opening, placeholder-H
shed) in step order.  Writes topology.jsonl so `export_xyz.py --format mol2` and
`check_topology_valence.py` work on the existing run — no re-run needed.

Usage:
    python scripts/reconstruct_topology.py runs/s6_half_50c
    python scripts/reconstruct_topology.py runs/s6_half_50c --n-monomers 100
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from kagome.reactive.topology import (
    BondTopology, apply_vinyl_addition, over_coordinated_atoms,
)


def _read_species_and_meta(run: Path) -> tuple[list[str], dict]:
    with open(run / 'trajectory.jsonl', encoding='utf-8') as f:
        header = json.loads(f.readline())
    meta = {}
    mf = run / 'manifest.json'
    if mf.exists():
        meta = json.loads(mf.read_text(encoding='utf-8'))
    return header['species'], meta


def _rebuild(path_kind: str, n_monomers: int, n_units: int, seed: int):
    from scripts._systems import (
        build_full_aibn_system, build_vinyl_aibn_system,
        full_aibn_initial_bonds, vinyl_initial_bonds,
    )
    rng = np.random.default_rng(seed)
    if path_kind == 'activation':
        _, sp, azo, _, groups, pmap, _ = build_full_aibn_system(
            n_monomers=n_monomers, n_aibn=n_units, box_size=40.0, rng=rng,
            rdkit_seed=seed)
        bonds = full_aibn_initial_bonds(n_monomers=n_monomers, n_aibn=n_units,
                                        rdkit_seed=seed)
        return sp, azo, groups, pmap, bonds
    _, sp, _, groups, pmap, _ = build_vinyl_aibn_system(
        n_monomers=n_monomers, n_initiators=n_units, box_size=40.0, rng=rng,
        rdkit_seed=seed)
    bonds = vinyl_initial_bonds(n_monomers=n_monomers, n_initiators=n_units,
                                rdkit_seed=seed)
    return sp, [], groups, pmap, bonds


def _detect_layout(species: list[str], n_monomers: int, seed: int):
    """Return (path_kind, n_units, azo, pmap, initial_bonds) whose rebuilt
    species matches *species* exactly."""
    n = len(species)
    candidates = []
    # activation: AIBN = 24 atoms/unit; non-activation: IBN = 12 atoms/unit.
    if (n - 12 * n_monomers) % 24 == 0 and (n - 12 * n_monomers) // 24 > 0:
        candidates.append(('activation', (n - 12 * n_monomers) // 24))
    if n % 12 == 0 and (n // 12 - n_monomers) > 0:
        candidates.append(('non-activation', n // 12 - n_monomers))
    for kind, n_units in candidates:
        sp, azo, groups, pmap, bonds = _rebuild(kind, n_monomers, n_units, seed)
        if sp == species:
            return kind, n_units, azo, pmap, bonds
    raise SystemExit(
        f'could not reconstruct layout for {n} atoms / {n_monomers} monomers; '
        'pass --n-monomers / --n-units explicitly.')


def main() -> None:
    p = argparse.ArgumentParser(description='Reconstruct topology.jsonl for a run')
    p.add_argument('run', type=Path, help='run directory')
    p.add_argument('--n-monomers', type=int, default=None,
                   help='default: manifest n_reactive_sites')
    p.add_argument('--seed', type=int, default=None,
                   help='default: manifest seed')
    args = p.parse_args()

    species, meta = _read_species_and_meta(args.run)
    extra = meta.get('extra', {})
    n_monomers = args.n_monomers or extra.get('n_reactive_sites')
    seed = args.seed if args.seed is not None else meta.get('seed', 7)
    if n_monomers is None:
        p.error('n_monomers unknown (no manifest n_reactive_sites); pass --n-monomers.')

    kind, n_units, azo, pmap, init_bonds = _detect_layout(species, n_monomers, seed)
    print(f'layout: {kind}, {n_units} initiator units, {n_monomers} monomers, '
          f'{len(species)} atoms, seed {seed}')

    topo = BondTopology.from_bonds(init_bonds)
    for c, nn in azo:  # activation dissociated the azo C-N bonds
        topo.remove_bond(c, nn)

    events = [json.loads(l) for l in open(args.run / 'bonds.jsonl', encoding='utf-8')]
    conf = sorted((e for e in events if e['event_type'] == 'confirmed_formation'),
                  key=lambda e: e['step'])

    out = args.run / 'topology.jsonl'
    with open(out, 'w', encoding='utf-8') as f:
        # initial snapshot (cycle -1) at the first recorded step, or 0.
        first_step = conf[0]['step'] if conf else 0
        f.write(json.dumps({
            'step': 0, 'cycle': -1, 'n_bonds': len(topo),
            'bonds': [[i, j, o] for i, j, o in topo.bonds()],
        }) + '\n')
        for e in conf:
            apply_vinyl_addition(topo, e['atom_a'], e['atom_b'], pmap, species)
            f.write(json.dumps({
                'step': e['step'], 'cycle': e['cycle'], 'n_bonds': len(topo),
                'bonds': [[i, j, o] for i, j, o in topo.bonds()],
            }) + '\n')

    over = over_coordinated_atoms(topo, species)
    print(f'replayed {len(conf)} confirmed formations -> {out}')
    print('valence: ' + ('OK - no over-coordinated atoms' if not over
                         else f'OVER-COORDINATED: {over}'))


if __name__ == '__main__':
    main()
