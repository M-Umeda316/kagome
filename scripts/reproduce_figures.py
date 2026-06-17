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

from src.analysis.conversion import conversion_timeseries, fit_conversion_exponential
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
        print('No temperature data in trajectory --skipping temperature plot.')
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
    timestep_fs: float = 1.0,
) -> None:
    """Plot conversion α(t) with Eq. 11 exponential fit overlay."""
    events = read_bond_events(bonds_path)

    if not events:
        print(f'No bond events in {bonds_path} --skipping conversion plot.')
        return

    step_range, alpha = conversion_timeseries(events, n_total_sites)

    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 4))
    time_ps = step_range * timestep_fs / 1000.0
    ax.plot(time_ps, alpha * 100.0, linewidth=1.0, color='tab:brown',
            label='Simulation')

    kp_eff, r_sq = fit_conversion_exponential(step_range, alpha, timestep_fs)
    if kp_eff > 0:
        t_fit = np.linspace(0, float(time_ps[-1]), 500)
        alpha_fit = (1.0 - np.exp(-kp_eff * t_fit * 1000.0)) * 100.0
        ax.plot(t_fit, alpha_fit, '--', color='tab:red', linewidth=1.2,
                label=f'Eq. 11 fit: $k_{{p,eff}}$={kp_eff:.2e} fs$^{{-1}}$')
        ax.text(0.98, 0.60, f'$R^2$ = {r_sq:.3f}',
                transform=ax.transAxes, ha='right', fontsize=9,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Conversion α (%)')
    ax.set_title('Conversion vs. time (Eq. 11)')
    ax.legend(loc='upper left')
    alpha_max = float(alpha.max()) * 100.0
    ax.set_ylim(0, max(alpha_max * 2, 5))
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
        print(f'No confirmed_formation events in {bonds_path} -- skipping density plot.')
        return

    _, frames = read_trajectory(trajectory_path)
    if not frames:
        print('No frames in trajectory --skipping density plot.')
        return

    event_steps = {e.step for e in formations}
    positions_at_event: dict[int, np.ndarray] = {
        f.step: np.array(f.positions)
        for f in frames
        if f.step in event_steps and f.positions
    }

    if not positions_at_event:
        print('No trajectory frames match bond event steps -- skipping density plot.')
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


def plot_s2_diagnostics(
    summary_paths: list[Path],
    output_dir: Path,
) -> None:
    """Plot min_pair_distance and candidates per cycle from summary.json files."""
    import json

    output_dir.mkdir(parents=True, exist_ok=True)

    n_runs = len(summary_paths)
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=False)
    ax_dist, ax_cand = axes

    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']

    for idx, spath in enumerate(summary_paths):
        with open(spath, encoding='utf-8') as f:
            summary = json.load(f)

        label = spath.parent.name
        color = colors[idx % len(colors)]

        biased_logs = [log for log in summary['logs'] if log['phase'] == 'biased']
        cycles = [log['cycle'] for log in biased_logs]
        min_dists = []
        candidates = []
        reaction_cycles = []

        for log in biased_logs:
            d = log.get('min_pair_distance')
            min_dists.append(d if d is not None else float('nan'))
            candidates.append(log['n_candidates'])
            if log.get('bias_energy', 250) < 200:
                reaction_cycles.append(log['cycle'])

        offset = idx * 0.15
        ax_dist.plot(cycles, min_dists, 'o-', color=color, markersize=5,
                     label=label, alpha=0.8)
        for rc in reaction_cycles:
            ri = cycles.index(rc)
            ax_dist.plot(rc, min_dists[ri], '*', color=color, markersize=14,
                         zorder=5)

        bar_width = 0.35
        bar_x = np.array(cycles) + offset - 0.15 * (n_runs - 1) / 2
        ax_cand.bar(bar_x, candidates, width=bar_width, color=color,
                    label=label, alpha=0.7)

    ax_dist.axhline(2.04, color='gray', linestyle=':', linewidth=0.8,
                    label='$r_0$ = 2.04 Å')
    ax_dist.axhline(2.87, color='salmon', linestyle='--', linewidth=0.8,
                    label='Capture radius (f2=10)')
    ax_dist.axhline(3.2, color='lightsalmon', linestyle='--', linewidth=0.8,
                    label='Capture radius (f2=5)')
    ax_dist.set_ylabel('min pair distance (Å)')
    ax_dist.set_title('S2 diagnostics: TDBB capture mechanism')
    ax_dist.legend(fontsize=7, ncol=2)
    ax_dist.set_ylim(0, None)

    ax_cand.set_xlabel('Cycle')
    ax_cand.set_ylabel('Candidates')
    ax_cand.legend(fontsize=8)

    fig.tight_layout()
    for fmt in ('png', 'pdf'):
        fig.savefig(output_dir / f's2_diagnostics.{fmt}', dpi=150)
    plt.close(fig)
    print(f'Saved s2_diagnostics.png/.pdf to {output_dir}')


def main() -> None:
    parser = argparse.ArgumentParser(description='Reproduce figures from trajectory')
    parser.add_argument('--trajectory', type=Path, default=None)
    parser.add_argument('--bonds', type=Path, default=None,
                        help='bonds.jsonl from BondTracker (optional, for conversion/density plots)')
    parser.add_argument('--n-reactive-sites', type=int, default=None,
                        dest='n_reactive_sites',
                        help=(
                            'Number of reactive sites (denominator for alpha(t)). '
                            'Auto-read from trajectory header if omitted. '
                            'Must NOT use n_atoms --only count atoms that can react.'
                        ))
    # Backward-compat alias
    parser.add_argument('--n-total-sites', type=int, default=None,
                        dest='n_total_sites_compat',
                        help='Deprecated alias for --n-reactive-sites.')
    parser.add_argument('--target-temperature', type=float, default=None,
                        help='Target temperature in K for reference line in temperature plot')
    parser.add_argument('--timestep-fs', type=float, default=1.0,
                        help='MD timestep in fs (for Eq. 11 exponential fit). Default 1.0.')
    parser.add_argument('--cell-xy-area', type=float, default=None,
                        help='Cross-sectional area (Å²) for density normalization. '
                             'Defaults to (z_range)² if omitted.')
    parser.add_argument('--summary', type=Path, nargs='+', default=None,
                        help='summary.json file(s) for S2 diagnostics plot.')
    parser.add_argument('--output-dir', type=Path, default=Path('runs/smoke/figures'))
    args = parser.parse_args()

    if args.trajectory is None and args.summary is None:
        parser.error('At least one of --trajectory or --summary is required.')

    if args.trajectory is not None and not args.trajectory.exists():
        print(f'Trajectory file not found: {args.trajectory}')
        return

    if args.trajectory is not None:
        plot_energy_vs_step(args.trajectory, args.output_dir)
        plot_temperature_vs_step(args.trajectory, args.output_dir, args.target_temperature)

    if args.trajectory is not None and args.bonds is not None:
        header, _ = read_trajectory(args.trajectory)

        # Priority: explicit CLI arg > trajectory header > deprecated fallback
        n_sites = args.n_reactive_sites
        if n_sites is None and args.n_total_sites_compat is not None:
            print('WARNING: --n-total-sites is deprecated; use --n-reactive-sites instead.')
            n_sites = args.n_total_sites_compat
        if n_sites is None:
            n_sites = header.get('n_reactive_sites')
        if n_sites is None:
            # Last resort: use n_atoms but warn loudly --this will over-estimate alpha
            n_atoms = header.get('n_atoms', 0)
            print(
                f'WARNING: n_reactive_sites not found in header or CLI args. '
                f'Falling back to n_atoms={n_atoms}, which over-estimates alpha(t). '
                f'Pass --n-reactive-sites to fix.'
            )
            n_sites = n_atoms

        plot_conversion_vs_step(args.bonds, n_sites, args.output_dir,
                               timestep_fs=args.timestep_fs)
        plot_density_profile(
            args.bonds, args.trajectory, args.output_dir,
            cell_xy_area=args.cell_xy_area,
        )

    if args.summary:
        plot_s2_diagnostics(args.summary, args.output_dir)


if __name__ == '__main__':
    main()
