"""Convert kagome JSONL trajectory to extended XYZ for Winmostar / OVITO / VMD.

Usage:
    python scripts/export_xyz.py runs/wsl_mace_gpu/trajectory.jsonl
    python scripts/export_xyz.py runs/wsl_mace_gpu/trajectory.jsonl --cell 14.0
    python scripts/export_xyz.py runs/wsl_mace_gpu/trajectory.jsonl -o out.xyz --cell 14.0 14.0 14.0
"""
from __future__ import annotations

import argparse
from pathlib import Path

from kagome.io.readers import read_trajectory


def main() -> None:
    parser = argparse.ArgumentParser(description='JSONL trajectory → XYZ')
    parser.add_argument('trajectory', type=Path)
    parser.add_argument('-o', '--output', type=Path, default=None)
    parser.add_argument(
        '--cell', type=float, nargs='+', default=None,
        help='Box lengths: single value for cubic, or Lx Ly Lz',
    )
    args = parser.parse_args()

    header, frames = read_trajectory(args.trajectory)
    species = header['species']
    n_atoms = header['n_atoms']

    out_path = args.output or args.trajectory.with_suffix('.xyz')

    lattice_str = ''
    pbc_str = ''
    if args.cell:
        if len(args.cell) == 1:
            lx = ly = lz = args.cell[0]
        elif len(args.cell) == 3:
            lx, ly, lz = args.cell
        else:
            parser.error('--cell expects 1 or 3 values')
        lattice_str = (
            f' Lattice="{lx:.6f} 0.0 0.0'
            f' 0.0 {ly:.6f} 0.0'
            f' 0.0 0.0 {lz:.6f}"'
        )
        pbc_str = ' pbc="T T T"'

    with open(out_path, 'w', encoding='utf-8') as f:
        for frame in frames:
            f.write(f'{n_atoms}\n')
            comment = (
                f'step={frame.step}'
                f' time_fs={frame.time_fs:.2f}'
                f' energy={frame.energy_total:.4f}'
                f' phase={frame.phase}'
                f' cycle={frame.cycle}'
                f'{lattice_str}{pbc_str}'
                f' Properties=species:S:1:pos:R:3'
            )
            f.write(comment + '\n')
            for i, (x, y, z) in enumerate(frame.positions):
                f.write(f'{species[i]:2s} {x:16.8f} {y:16.8f} {z:16.8f}\n')

    print(f'{len(frames)} frames → {out_path}')


if __name__ == '__main__':
    main()
