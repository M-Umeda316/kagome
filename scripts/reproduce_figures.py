"""Reproduce publication-style figures from trajectory data.

Usage:
    python scripts/reproduce_figures.py --trajectory runs/smoke/trajectory.jsonl --output-dir runs/smoke/figures
    python scripts/reproduce_figures.py --trajectory ... --bonds ... --output-dir ...

Paper: arXiv:2511.22874, Figs. 2-6 require time-series of energy and conversion.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    raise ImportError(
        'matplotlib is required for figure generation. '
        'Install with: pip install pfpoly[plot]'
    )

from src.analysis.conversion import conversion_timeseries
from src.analysis.density import reaction_density_profile
from src.io.readers import read_bond_events, read_trajectory


def plot_energy_vs_step(
    trajectory_path: Path,
    output_dir: Path,
) -> None:
    header, frames = read_trajectory(trajectory_path)

    if not frames:
        print('No frames in trajectory, nothing to plot.')
        return

    steps = np.array([f.step for f in frames])
    e_base = np.array([f.energy_base for f in frames])
    e_bias = np.array([f.energy_bias for f in frames])
    e_total = np.array([f.energy_total for f in frames])
    phases = [f.phase for f in frames]

    biased_mask = np.array([p == 'biased' for p in phases])
    unbiased_mask = ~biased_mask

    output_dir.mkdir(parents=True, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

    if np.any(biased_mask):
        ax1.scatter(steps[biased_mask], e_total[biased_mask],
                    s=4, c='tab:red', label='biased', alpha=0.7)
    if np.any(unbiased_mask):
        ax1.scatter(steps[unbiased_mask], e_total[unbiased_mask],
                    s=4, c='tab:blue', label='unbiased', alpha=0.7)
    ax1.set_ylabel('Total energy (kcal/mol)')
    ax1.legend(markerscale=3)
    ax1.set_title('Energy vs. simulation step')

    ax2.scatter(steps[biased_mask] if np.any(biased_mask) else [],
                e_bias[biased_mask] if np.any(biased_mask) else [],
                s=4, c='tab:orange', alpha=0.7)
    ax2.set_ylabel('Bias energy (kcal/mol)')
    ax2.set_xlabel('Step')

    fig.tight_layout()
    for fmt in ('png', 'pdf'):
        fig.savefig(output_dir / f'energy_vs_step.{fmt}', dpi=150)
    plt.close(fig)
    print(f'Saved energy_vs_step.png/.pdf to {output_dir}')

    fig2, ax = plt.subplots(figsize=(10, 4))
    ax.plot(steps, e_base, linewidth=0.5, color='tab:green')
    ax.set_xlabel('Step')
    ax.set_ylabel('Base energy (kcal/mol)')
    ax.set_title('Base (unbiased) potential energy')
    fig2.tight_layout()
    for fmt in ('png', 'pdf'):
        fig2.savefig(output_dir / f'base_energy.{fmt}', dpi=150)
    plt.close(fig2)
    print(f'Saved base_energy.png/.pdf to {output_dir}')


def plot_temperature_vs_step(
    trajectory_path: Path,
    output_dir: Path,
    target_temperature_K: float | None = None,
) -> None:
    """Plot instantaneous kinetic temperature vs. step (NVT validation)."""
    _, frames = read_trajectory(trajectory_path)

    temps = [f.temperature_K for f in frames]
    if not any(t > 0.0 for t in temps):
        print('No temperature data in trajectory — skipping temperature plot.')
        return

    steps = np.array([f.step for f in frames])
    temp_arr = np.array(temps)

    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(steps, temp_arr, linewidth=0.5, color='tab:purple', label='T_inst')
    if target_temperature_K is not None:
        ax.axhline(target_temperature_K, color='k', linestyle='--',
                   linewidth=0.8, label=f'T_target={target_temperature_K:.0f} K')
    ax.set_xlabel('Step')
    ax.set_ylabel('Temperature (K)')
    ax.set_title('Instantaneous kinetic temperature')
    ax.legend()
    fig.tight_layout()
    for fmt in ('png', 'pdf'):
        fig.savefig(output_dir / f'temperature_vs_step.{fmt}', dpi=150)
    plt.close(fig)
    print(f'Saved temperature_vs_step.png/.pdf to {output_dir}')


def plot_conversion_vs_step(
    bonds_path: Path,
    n_total_sites: int,
    output_dir: Path,
) -> None:
    """Plot conversion α(t) from bonds.jsonl.  Eq. 11-12."""
    events = read_bond_events(bonds_path)

    if not events:
        print(f'No bond events in {bonds_path} — skipping conversion plot.')
        return

    step_range, alpha = conversion_timeseries(events, n_total_sites)

    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(step_range, alpha * 100.0, linewidth=1.0, color='tab:brown')
    ax.set_xlabel('Step')
    ax.set_ylabel('Conversion α (%)')
    ax.set_title('Degree of polymerization vs. simulation step (Eq. 11-12)')
    ax.set_ylim(0, 105)
    fig.tight_layout()
    for fmt in ('png', 'pdf'):
        fig.savefig(output_dir / f'conversion_vs_step.{fmt}', dpi=150)
    plt.close(fig)
    print(f'Saved conversion_vs_step.png/.pdf to {output_dir}')


def plot_density_profile(
    bonds_path: Path,
    trajectory_path: Path,
    output_dir: Path,
    n_z_bins: int = 20,
    cell_xy_area: float | None = None,
) -> None:
    """Plot depth-resolved reaction density ρ_rxn(z).  Eq. 12 (arXiv HTML).

    Positions at each bond-event step are extracted from the trajectory.
    Requires --bonds with confirmed_formation events and position data.
    """
    events = read_bond_events(bonds_path)
    formations = [e for e in events if e.event_type == 'confirmed_formation']
    if not formations:
        print(f'No confirmed_formation events in {bonds_path} — skipping density plot.')
        return

    _, frames = read_trajectory(trajectory_path)
    if not frames:
        print('No frames in trajectory — skipping density plot.')
        return

    event_steps = {e.step for e in formations}
    positions_at_event: dict[int, np.ndarray] = {
        f.step: np.array(f.positions)
        for f in frames
        if f.step in event_steps and f.positions
    }

    if not positions_at_event:
        print('No trajectory frames match bond event steps — skipping density plot.')
        return

    all_z = np.concatenate([pos[:, 2] for pos in positions_at_event.values()])
    z_min, z_max = float(all_z.min()), float(all_z.max())
    z_bins = np.linspace(z_min, z_max, n_z_bins + 1)

    area_xy = cell_xy_area if cell_xy_area is not None else (z_max - z_min) ** 2
    density = reaction_density_profile(
        formations, positions_at_event, z_bins, area_xy,
        n_frames=len(positions_at_event),
    )

    z_centers = 0.5 * (z_bins[:-1] + z_bins[1:])

    output_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 8))
    ax.barh(z_centers, density, height=np.diff(z_bins) * 0.9, color='tab:cyan')
    ax.set_xlabel('Reaction density ρ_rxn (1/Å³)')
    ax.set_ylabel('z position (Å)')
    ax.set_title('Depth-resolved reaction density (Eq. 12)')
    fig.tight_layout()
    for fmt in ('png', 'pdf'):
        fig.savefig(output_dir / f'density_profile.{fmt}', dpi=150)
    plt.close(fig)
    print(f'Saved density_profile.png/.pdf to {output_dir}')


def main() -> None:
    parser = argparse.ArgumentParser(description='Reproduce figures from trajectory')
    parser.add_argument('--trajectory', type=Path, required=True)
    parser.add_argument('--bonds', type=Path, default=None,
                        help='bonds.jsonl from BondTracker (optional, for conversion/density plots)')
    parser.add_argument('--n-reactive-sites', type=int, default=None,
                        dest='n_reactive_sites',
                        help=(
                            'Number of reactive sites (denominator for alpha(t)). '
                            'Auto-read from trajectory header if omitted. '
                            'Must NOT use n_atoms — only count atoms that can react.'
                        ))
    # Backward-compat alias
    parser.add_argument('--n-total-sites', type=int, default=None,
                        dest='n_total_sites_compat',
                        help='Deprecated alias for --n-reactive-sites.')
    parser.add_argument('--target-temperature', type=float, default=None,
                        help='Target temperature in K for reference line in temperature plot')
    parser.add_argument('--cell-xy-area', type=float, default=None,
                        help='Cross-sectional area (Å²) for density normalization. '
                             'Defaults to (z_range)² if omitted.')
    parser.add_argument('--output-dir', type=Path, default=Path('runs/smoke/figures'))
    args = parser.parse_args()

    if not args.trajectory.exists():
        print(f'Trajectory file not found: {args.trajectory}')
        return

    plot_energy_vs_step(args.trajectory, args.output_dir)
    plot_temperature_vs_step(args.trajectory, args.output_dir, args.target_temperature)

    if args.bonds is not None:
        header, _ = read_trajectory(args.trajectory)

        # Priority: explicit CLI arg > trajectory header > deprecated fallback
        n_sites = args.n_reactive_sites
        if n_sites is None and args.n_total_sites_compat is not None:
            print('WARNING: --n-total-sites is deprecated; use --n-reactive-sites instead.')
            n_sites = args.n_total_sites_compat
        if n_sites is None:
            n_sites = header.get('n_reactive_sites')
        if n_sites is None:
            # Last resort: use n_atoms but warn loudly — this will over-estimate alpha
            n_atoms = header.get('n_atoms', 0)
            print(
                f'WARNING: n_reactive_sites not found in header or CLI args. '
                f'Falling back to n_atoms={n_atoms}, which over-estimates alpha(t). '
                f'Pass --n-reactive-sites to fix.'
            )
            n_sites = n_atoms

        plot_conversion_vs_step(args.bonds, n_sites, args.output_dir)
        plot_density_profile(
            args.bonds, args.trajectory, args.output_dir,
            cell_xy_area=args.cell_xy_area,
        )


if __name__ == '__main__':
    main()
