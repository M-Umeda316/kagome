"""Quantitative valence / geometry check for a run's emitted topology.

Companion to `export_xyz.py --format mol2` visual inspection: reports, per
topology snapshot, any over-coordinated atoms (should be 0 with the Layer 1/2
topology + valence guard), and, per trajectory frame, the closest NON-bonded
contact as a proxy for geometric strain (whole-molecule drag under bias).  Use
this to decide whether post-formation relaxation (Layer 3(A)) is needed.

Usage:
    python scripts/check_topology_valence.py runs/foo
    python scripts/check_topology_valence.py runs/foo/trajectory.jsonl \
        --topology runs/foo/topology.jsonl --contact-cutoff 1.2
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from kagome.io.readers import bonds_at_step, read_topology_snapshots, read_trajectory
from kagome.reactive.topology import BondTopology, over_coordinated_atoms

_MAXC = {'H': 1, 'C': 4, 'N': 3, 'O': 2, 'F': 1, 'S': 2, 'Cl': 1}


def _resolve_paths(target: Path) -> tuple[Path, Path]:
    if target.is_dir():
        return target / 'trajectory.jsonl', target / 'topology.jsonl'
    return target, target.with_name('topology.jsonl')


def _min_nonbonded_distance(
    positions: np.ndarray, bonds: list[tuple[int, int, float]], cutoff: float,
    box: np.ndarray | None = None,
) -> tuple[float, int]:
    """Return (min non-bonded distance, count below cutoff) using a KD-tree.

    When *box* is a (3,) array of cell lengths, PBC-aware distances are used
    via cKDTree(boxsize=...) (M6 fix).
    """
    try:
        from scipy.spatial import cKDTree
    except ImportError:
        return float('nan'), -1
    bonded = {(min(i, j), max(i, j)) for i, j, _ in bonds}
    tree_kw: dict = {}
    if box is not None:
        tree_kw['boxsize'] = box.astype(float)
    tree = cKDTree(positions, **tree_kw)
    pairs = tree.query_pairs(r=cutoff, output_type='ndarray')
    below = 0
    gmin = float('inf')
    dists, idx = tree.query(positions, k=2)
    for a, (d, b) in enumerate(zip(dists[:, 1], idx[:, 1])):
        key = (min(a, int(b)), max(a, int(b)))
        if key not in bonded and d < gmin:
            gmin = d
    for a, b in pairs:
        key = (min(int(a), int(b)), max(int(a), int(b)))
        if key not in bonded:
            below += 1
    return gmin, below


def main() -> None:
    p = argparse.ArgumentParser(description='Valence / close-contact check for a run')
    p.add_argument('target', type=Path, help='run dir or trajectory.jsonl')
    p.add_argument('--topology', type=Path, default=None)
    p.add_argument('--contact-cutoff', type=float, default=1.2,
                   help='report non-bonded pairs closer than this (A). Default 1.2.')
    p.add_argument('--box', type=float, nargs=3, default=None, metavar=('LX', 'LY', 'LZ'),
                   help='PBC box lengths in A for periodic distance calc (M6).')
    args = p.parse_args()

    traj_path, topo_default = _resolve_paths(args.target)
    topo_path = args.topology or topo_default

    box = np.array(args.box) if args.box else None
    header, frames = read_trajectory(traj_path)
    species = header['species']
    snapshots = read_topology_snapshots(topo_path)
    if not snapshots:
        p.error(f'no topology snapshots at {topo_path}; re-run with the '
                'topology-aware workflow.')

    print(f'run: {traj_path.parent}')
    print(f'atoms: {len(species)} | frames: {len(frames)} | '
          f'topology snapshots: {len(snapshots)}')

    # 1) Valence: over-coordination per snapshot.
    total_over = 0
    for step, bonds in snapshots:
        topo = BondTopology.from_bonds(bonds)
        over = over_coordinated_atoms(topo, species)
        if over:
            total_over += len(over)
            print(f'  [step {step}] OVER-COORDINATED: '
                  f'{[(a, species[a], topo.coordination_number(a)) for a in over]}')
    print('valence: ' + ('OK - no over-coordinated atoms' if total_over == 0
                          else f'{total_over} over-coordinated instances'))

    # 2) Geometry: closest non-bonded contact per frame.
    worst = float('inf')
    worst_step = None
    n_strained = 0
    for frame in frames:
        bonds = bonds_at_step(snapshots, frame.step)
        pos = np.asarray(frame.positions, dtype=float)
        gmin, below = _min_nonbonded_distance(pos, bonds, args.contact_cutoff, box=box)
        if below == -1:
            print('geometry: scipy not available - skipping close-contact scan')
            break
        if below > 0:
            n_strained += 1
        if gmin < worst:
            worst, worst_step = gmin, frame.step
    else:
        print(f'geometry: closest non-bonded contact = {worst:.2f} A '
              f'(at step {worst_step}); '
              f'{n_strained}/{len(frames)} frames have contacts < '
              f'{args.contact_cutoff} A')
        if worst < 1.0:
            print('  -> non-bonded atoms interpenetrate (<1.0 A): post-formation '
                  'relaxation (Layer 3(A)) likely worthwhile.')


if __name__ == '__main__':
    main()
