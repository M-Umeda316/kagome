"""Convert a kagome JSONL trajectory to XYZ or MOL2 for Winmostar / OVITO / VMD.

XYZ carries coordinates only, so viewers infer bonds by distance — under TDBB
bias whole molecules are dragged together and the viewer draws spurious bonds at
non-reactive sites, and the never-opened C=C makes reacted carbons look
over-valent (specs/decisions.md 2026-07-02).  MOL2 carries the *explicit* bond
topology tracked during the run (topology.jsonl), with per-frame connectivity
and bond orders, so the viewer shows real chemistry.

Usage:
    # coordinates only (legacy behaviour)
    python scripts/export_xyz.py runs/foo/trajectory.jsonl --cell 14.0

    # with explicit bonds (recommended): reads runs/foo/topology.jsonl
    python scripts/export_xyz.py runs/foo/trajectory.jsonl --format mol2
    python scripts/export_xyz.py runs/foo/trajectory.jsonl --format mol2 \
        --topology runs/foo/topology.jsonl -o out.mol2
"""
from __future__ import annotations

import argparse
from pathlib import Path

from kagome.io.readers import bonds_at_step, read_topology_snapshots, read_trajectory

# RDKit/Sybyl bond-order tokens for MOL2 @<TRIPOS>BOND.
_MOL2_BOND_TYPE = {1.0: '1', 2.0: '2', 3.0: '3', 1.5: 'ar'}


def _write_xyz(out_path, frames, species, n_atoms, cell) -> None:
    lattice_str = ''
    pbc_str = ''
    if cell:
        lx, ly, lz = cell
        lattice_str = (
            f' Lattice="{lx:.6f} 0.0 0.0'
            f' 0.0 {ly:.6f} 0.0'
            f' 0.0 0.0 {lz:.6f}"'
        )
        pbc_str = ' pbc="T T T"'

    with open(out_path, 'w', encoding='utf-8') as f:
        for frame in frames:
            f.write(f'{n_atoms}\n')
            f.write(
                f'step={frame.step} time_fs={frame.time_fs:.2f}'
                f' energy={frame.energy_total:.4f} phase={frame.phase}'
                f' cycle={frame.cycle}{lattice_str}{pbc_str}'
                f' Properties=species:S:1:pos:R:3\n'
            )
            for i, (x, y, z) in enumerate(frame.positions):
                f.write(f'{species[i]:2s} {x:16.8f} {y:16.8f} {z:16.8f}\n')


def _write_mol2(out_path, frames, species, n_atoms, snapshots, header_bonds) -> None:
    """Multi-frame MOL2 with per-frame explicit bonds.

    Each frame is a @<TRIPOS>MOLECULE block; the bond list is the topology in
    effect at that frame's step (falls back to the trajectory header's initial
    bonds when no topology.jsonl is available).
    """
    with open(out_path, 'w', encoding='utf-8') as f:
        for frame in frames:
            if snapshots:
                bonds = bonds_at_step(snapshots, frame.step)
            else:
                bonds = header_bonds
            f.write('@<TRIPOS>MOLECULE\n')
            f.write(f'step{frame.step}_cycle{frame.cycle}_{frame.phase}\n')
            f.write(f'{n_atoms} {len(bonds)} 0 0 0\n')
            f.write('SMALL\nNO_CHARGES\n')
            f.write('@<TRIPOS>ATOM\n')
            for i, (x, y, z) in enumerate(frame.positions):
                el = species[i]
                # Sybyl type: plain element is accepted by the common viewers.
                f.write(
                    f'{i + 1:7d} {el}{i + 1:<7d} {x:12.6f} {y:12.6f} {z:12.6f} '
                    f'{el:<5s} 1 RES 0.0000\n'
                )
            f.write('@<TRIPOS>BOND\n')
            for bid, (a, b, order) in enumerate(bonds, start=1):
                btype = _MOL2_BOND_TYPE.get(round(float(order), 1), '1')
                f.write(f'{bid:7d} {a + 1:7d} {b + 1:7d} {btype}\n')


def main() -> None:
    parser = argparse.ArgumentParser(description='JSONL trajectory → XYZ / MOL2')
    parser.add_argument('trajectory', type=Path)
    parser.add_argument('-o', '--output', type=Path, default=None)
    parser.add_argument(
        '--format', choices=('xyz', 'mol2'), default='xyz',
        help='xyz = coordinates only (viewer guesses bonds by distance); '
             'mol2 = explicit bonds from topology.jsonl (recommended).',
    )
    parser.add_argument(
        '--topology', type=Path, default=None,
        help='topology.jsonl path (default: sibling of the trajectory).',
    )
    parser.add_argument(
        '--cell', type=float, nargs='+', default=None,
        help='XYZ only — box lengths: single value for cubic, or Lx Ly Lz',
    )
    parser.add_argument(
        '--stride', type=int, default=1,
        help='keep every Nth frame (default 1 = all). Use for large trajectories.',
    )
    args = parser.parse_args()

    header, frames = read_trajectory(args.trajectory)
    species = header['species']
    n_atoms = header['n_atoms']
    if args.stride > 1:
        frames = frames[::args.stride]

    cell = None
    if args.cell:
        if len(args.cell) == 1:
            cell = (args.cell[0],) * 3
        elif len(args.cell) == 3:
            cell = tuple(args.cell)
        else:
            parser.error('--cell expects 1 or 3 values')

    if args.format == 'mol2':
        topo_path = args.topology or args.trajectory.with_name('topology.jsonl')
        snapshots = read_topology_snapshots(topo_path)
        header_bonds = [tuple(b) for b in header.get('bonds', [])]
        if not snapshots and not header_bonds:
            parser.error(
                f'no explicit bonds found (looked for {topo_path} and a "bonds" '
                'header). Re-run with the topology-aware workflow, or use '
                '--format xyz.'
            )
        out_path = args.output or args.trajectory.with_suffix('.mol2')
        _write_mol2(out_path, frames, species, n_atoms, snapshots, header_bonds)
        src = 'topology.jsonl' if snapshots else 'header bonds'
        print(f'{len(frames)} frames → {out_path} (explicit bonds from {src})')
    else:
        out_path = args.output or args.trajectory.with_suffix('.xyz')
        _write_xyz(out_path, frames, species, n_atoms, cell)
        print(f'{len(frames)} frames → {out_path}')


if __name__ == '__main__':
    main()
