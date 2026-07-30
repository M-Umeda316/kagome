"""Reproduce publication-style figures from trajectory data.

Usage:
    python scripts/reproduce_figures.py --trajectory runs/smoke/trajectory.jsonl --output-dir runs/smoke/figures
    python scripts/reproduce_figures.py --trajectory ... --bonds ... --output-dir ...

Paper: arXiv:2511.22874, Figs. 2-6 require time-series of energy and conversion.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    raise ImportError(
        'matplotlib is required for figure generation. '
        'Install with: pip install kagome[plot]'
    )

from kagome.analysis.carothers import (
    dpn_carothers,
    dpn_measured_from_topology,
    monomer_sets_from_bonds,
)
from kagome.analysis.conversion import conversion_timeseries, fit_conversion_exponential
from kagome.analysis.density import reaction_density_profile
from kagome.analysis.network import (
    find_epoxide_rings,
    gel_point_flory_stockmayer,
    largest_component_fraction,
    load_topology_snapshots,
    max_inter_monomer_degree,
    species_series,
)
from kagome.io.readers import (
    read_bond_events,
    read_topology_snapshots,
    read_trajectory,
)


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
    # WM-P3 (decisions.md 追補 2026-07-18 (k)): keep the mixing-mode settle
    # transient (and any classical-stage frames) out of the "unbiased" series.
    # Runs without mixing carry neither phase, so their figures are unchanged.
    unbiased_mask = ~biased_mask & np.array(
        [p not in ('mixing', 'mix_settle') for p in phases])

    output_dir.mkdir(parents=True, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

    if np.any(biased_mask):
        ax1.scatter(steps[biased_mask], e_total[biased_mask],
                    s=4, c='tab:red', label='biased', alpha=0.7)
    if np.any(unbiased_mask):
        ax1.scatter(steps[unbiased_mask], e_total[unbiased_mask],
                    s=4, c='tab:blue', label='unbiased', alpha=0.7)
    # OrbMol's compute() returns the interaction energy (model forward), NOT the
    # absolute energy with the per-element reference offset that predict() adds
    # (decisions.md 追補 2026-07-22 item 5). The absolute value is therefore not
    # physical; within a fixed-composition segment the offset is constant, so
    # trends/conservation are faithful. Label it truthfully rather than "Total".
    ax1.set_ylabel('Interaction energy + bias (kcal/mol)')
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

    # WM-P3 (decisions.md 追補 2026-07-18 (k)): the base-energy series is the
    # MLIP unbiased potential; exclude the classical->MLIP settle transient so
    # it does not show a spurious spike at each cycle boundary. Unchanged for
    # runs without mixing (no such frames).
    fig2, ax = plt.subplots(figsize=(10, 4))
    ax.plot(steps[unbiased_mask], e_base[unbiased_mask],
            linewidth=0.5, color='tab:green')
    ax.set_xlabel('Step')
    # Interaction energy (no per-element reference offset — see plot_energy_vs_step
    # comment and decisions.md 追補 2026-07-22 item 5), unbiased segment only.
    ax.set_ylabel('Interaction energy, unbiased (kcal/mol)')
    ax.set_title('Base (unbiased) interaction energy')
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
    production_start_step: int = 0,
    step_range: np.ndarray | None = None,
) -> None:
    """Plot conversion α(t) with Eq. 11 exponential fit overlay.

    production_start_step shifts the fit's t=0 to the start of production so
    k_p is not under-estimated by the equilibration/activation lead-in (L8/A4).
    step_range, when given, samples α on that grid (e.g. trajectory frame steps)
    instead of every MD step, avoiding a several-hundred-thousand-iteration
    Python loop (A6).
    """
    events = read_bond_events(bonds_path)

    if not events:
        print(f'No bond events in {bonds_path} --skipping conversion plot.')
        return

    # A5: exclude bias-only water-forming events (nylon k-l) from the reaction
    # count so one condensation = one amide bond. Old bonds.jsonl lack the field
    # and load as counts_as_reaction=True (unchanged behaviour for vinyl).
    events = [e for e in events if e.counts_as_reaction]

    step_range, alpha = conversion_timeseries(
        events, n_total_sites, step_range=step_range,
    )

    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 4))
    time_ps = step_range * timestep_fs / 1000.0
    ax.plot(time_ps, alpha * 100.0, linewidth=1.0, color='tab:brown',
            label='Simulation')

    kp_eff, r_sq = fit_conversion_exponential(
        step_range, alpha, timestep_fs,
        production_start_step=production_start_step,
    )
    if kp_eff > 0:
        # The fit is defined in production-relative time, so overlay it starting
        # at the production onset on the raw-time axis.
        t_prod_ps = production_start_step * timestep_fs / 1000.0
        t_fit = np.linspace(t_prod_ps, float(time_ps[-1]), 500)
        alpha_fit = (1.0 - np.exp(-kp_eff * (t_fit - t_prod_ps) * 1000.0)) * 100.0
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


def plot_dpn_vs_conversion(
    bonds_path: Path,
    topology_path: Path,
    n_reactive_sites: int,
    output_dir: Path,
) -> None:
    """Fig. 4c: MEASURED number-average DPn vs conversion p, over Carothers theory.

    Measured DPn at each recorded topology snapshot (cycle) is the monomer count
    divided by the number of connected molecules in the monomer graph
    (:func:`dpn_measured_from_topology`); the x-coordinate p is the counted
    amide-bond fraction from bonds.jsonl at that snapshot's step. The theoretical
    Carothers curve DPn = 1/(1-p) is overlaid. Step-growth only — the caller must
    gate this to nylon-style runs (chain-growth does not obey 1/(1-p)).

    Paper anchor: arXiv:2511.22874, Fig. 4c.
    """
    snapshots = read_topology_snapshots(topology_path)
    if not snapshots:
        print(f'No topology snapshots in {topology_path} -- skipping DPn plot.')
        return
    if n_reactive_sites <= 0:
        print('n_reactive_sites <= 0 -- skipping DPn plot.')
        return

    # Monomer atom membership = connected components of the INITIAL topology
    # (first snapshot, cycle -1): every monomer atom is intramolecularly bonded.
    monomer_atom_sets = monomer_sets_from_bonds(snapshots[0][1])
    if not monomer_atom_sets:
        print('No monomers recovered from initial topology -- skipping DPn plot.')
        return

    # Counted amide formations (A5: exclude bias-only water-forming k-l events so
    # one condensation = one amide bond), for the extent of reaction p.
    events = read_bond_events(bonds_path)
    counted_steps = sorted(
        e.step for e in events
        if e.event_type == 'confirmed_formation' and e.counts_as_reaction
    )
    counted_arr = np.array(counted_steps, dtype=np.int64)
    n_groups_each = n_reactive_sites / 2  # equimolar A-A + B-B: amine == carboxyl

    p_series: list[float] = []
    dpn_series: list[float] = []
    for step, bonds in snapshots:
        n_counted = int(np.searchsorted(counted_arr, step, side='right'))
        p = n_counted / n_groups_each if n_groups_each > 0 else 0.0
        p_series.append(min(p, 1.0))
        dpn_series.append(dpn_measured_from_topology(bonds, monomer_atom_sets))

    p_arr = np.array(p_series)
    dpn_arr = np.array(dpn_series)

    output_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 5))

    p_max = float(p_arr.max()) if p_arr.size else 0.0
    p_curve = np.linspace(0.0, min(max(p_max, 0.05), 0.99), 200)
    ax.plot(p_curve, dpn_carothers(p_curve), '-', color='tab:gray',
            linewidth=1.2, label=r'Carothers $1/(1-p)$')
    ax.plot(p_arr, dpn_arr, 'o-', color='tab:blue', markersize=6,
            linewidth=1.0, label='Simulation (measured)')

    ax.set_xlabel('Conversion p (extent of reaction)')
    ax.set_ylabel(r'Number-average DPn')
    ax.set_title('Carothers: measured DPn vs conversion (Fig. 4c)')
    ax.legend(loc='upper left')
    ax.set_xlim(0, max(p_max * 1.1, 0.05))
    ax.set_ylim(1.0, None)
    fig.tight_layout()
    for fmt in ('png', 'pdf'):
        fig.savefig(output_dir / f'dpn_vs_conversion.{fmt}', dpi=150)
    plt.close(fig)
    print(f'Saved dpn_vs_conversion.png/.pdf to {output_dir}')


def _analysis_frames(frames: list) -> list:
    """Frames in the MLIP analysis window, excluding classical-mixing frames.

    WM-P2 (decisions.md 追補 2026-07-18 (k) extension): mirrors the
    energy-plot ``unbiased_mask`` exclusion (~line 68-69) so that N_frames in
    the density-profile denominator does not count ``mix_settle`` frames (real
    MLIP frames written by ``_run_plain_md`` during a classical-mixing run's
    settle stage) or ``mixing`` frames (defensive; classical-stage frames, if
    ever written to the trajectory). No reaction can fire during either phase,
    so including them would systematically depress rho_rxn(z) for mixing runs.
    Runs without mixing carry neither phase on any frame, so this is a no-op
    (returns all frames unchanged) for non-mixing trajectories.
    """
    return [f for f in frames if f.phase not in ('mixing', 'mix_settle')]


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

    N_frames (below) counts only frames in the MLIP analysis window: it
    excludes classical mix_settle (and mixing) frames from a mixing-mode run,
    the same exclusion the energy plot applies (see ``_analysis_frames``).
    """
    events = read_bond_events(bonds_path)
    # A5: rho_rxn is the density of *reactions*. Exclude water-forming (and other
    # counts_as_reaction=False) formation events so a condensation reaction is
    # placed once, at its primary bond. Missing field -> True (vinyl unaffected).
    formations = [
        e for e in events
        if e.event_type == 'confirmed_formation' and e.counts_as_reaction
    ]
    if not formations:
        print(f'No confirmed_formation events in {bonds_path} -- skipping density plot.')
        return

    _, frames = read_trajectory(trajectory_path)
    if not frames:
        print('No frames in trajectory --skipping density plot.')
        return

    # WM-P2 (decisions.md 追補 2026-07-18 (k) extension): exclude classical
    # mix_settle/mixing frames from the analysis window, mirroring the
    # energy-plot's unbiased_mask (~line 68-69). Event steps only ever land in
    # biased/unbiased segments, so positions_at_event/cells_at_event would be
    # unaffected either way; analysis_frames is used everywhere for
    # consistency, and N_frames below MUST use it (that is the actual bug fix).
    analysis_frames = _analysis_frames(frames)

    event_steps = {e.step for e in formations}
    positions_at_event: dict[int, np.ndarray] = {
        f.step: np.array(f.positions)
        for f in analysis_frames
        if f.step in event_steps and f.positions
    }
    # Per-event cell for the PBC midpoint correction (NPT: box varies per frame).
    cells_at_event: dict[int, np.ndarray] = {
        f.step: np.array(f.cell)
        for f in analysis_frames
        if f.step in event_steps and f.cell is not None
    }

    if not positions_at_event:
        print('No trajectory frames match bond event steps -- skipping density plot.')
        return

    # Representative cell (mean over all sampled frames that recorded one). Used
    # for the default cross-sectional area and z-bin range so the profile follows
    # the paper definition rho_rxn(z) = N_rxn(z)/(A*dz*N_frames) instead of the
    # data min/max span. See specs/decisions.md 2026-07-06 (A1/A2/A3).
    frame_cells = [np.array(f.cell) for f in analysis_frames if f.cell is not None]

    if cell_xy_area is not None:
        area_xy = cell_xy_area
    elif frame_cells:
        mean_cell = np.mean(frame_cells, axis=0)
        area_xy = float(abs(mean_cell[0, 0] * mean_cell[1, 1]))
    else:
        raise ValueError(
            'Density profile needs a cross-sectional area: no cell recorded in '
            'the trajectory frames and --cell-xy-area not given. Re-run with a '
            'periodic trajectory (frames carry a cell) or pass --cell-xy-area. '
            'The old (z_max-z_min)**2 fallback was not a valid xy area.'
        )

    if frame_cells:
        # Wrapped positions live in [0, Lz); bin over the full box height.
        lz = float(np.mean([c[2, 2] for c in frame_cells]))
        z_bins = np.linspace(0.0, lz, n_z_bins + 1)
    else:
        all_z = np.concatenate([pos[:, 2] for pos in positions_at_event.values()])
        z_bins = np.linspace(float(all_z.min()), float(all_z.max()), n_z_bins + 1)

    # N_frames is every sampled trajectory frame in the analysis window (paper
    # definition), NOT the number of event steps (A2) -- and NOT including
    # classical mix_settle/mixing frames from a mixing-mode run (WM-P2), since
    # no reaction can fire during those and counting them would depress the
    # density.
    density = reaction_density_profile(
        formations, positions_at_event, z_bins, area_xy,
        n_frames=len(analysis_frames),
        cells_at_event=cells_at_event or None,
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


def plot_species_concentrations(
    topology_path: Path,
    species: list[str],
    output_dir: Path,
    timestep_fs: float = 1.0,
) -> None:
    """Fig. 5-style species concentration traces for epoxy-amine curing.

    Raw per-snapshot counts from topology.jsonl — NO smoothing, filtering, or
    averaging (paper Fig. 5 plots c(t)/c(0)).  Convention: species with a
    nonzero initial count (epoxide, 1° amine, 2° amine) are plotted as
    c(t)/c(0) on the left axis; produced species whose initial count is zero
    (3° amine, hydroxyl — c/c0 undefined) are plotted as raw counts on a twin
    right axis.  The x-axis is time in ps (step * timestep_fs / 1000); with
    the default timestep_fs=1.0 it is numerically the step count in fs.

    Paper anchor: arXiv:2511.22874 Fig. 5; design: specs/decisions.md
    2026-07-09 Track 2 / E2 network-analysis entry.
    """
    snapshots = load_topology_snapshots(topology_path)
    if not snapshots:
        print(f'No topology snapshots in {topology_path} -- skipping species plot.')
        return

    series = species_series(snapshots, species)
    steps = np.array([row['step'] for row in series], dtype=np.float64)
    time_ps = steps * timestep_fs / 1000.0

    traces = [
        ('n_epoxide', 'epoxide', 'tab:blue'),
        ('n_amine_primary', '1° amine', 'tab:green'),
        ('n_amine_secondary', '2° amine', 'tab:orange'),
        ('n_amine_tertiary', '3° amine', 'tab:red'),
        ('n_hydroxyl', 'hydroxyl', 'tab:purple'),
    ]

    output_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax_counts = ax.twinx()

    ratio_max = 1.0
    for key, label, color in traces:
        counts = np.array([row[key] for row in series], dtype=np.float64)
        c0 = counts[0]
        if c0 > 0:
            ratios = counts / c0
            ratio_max = max(ratio_max, float(ratios.max()))
            ax.plot(time_ps, ratios, 'o-', color=color, markersize=5,
                    linewidth=1.0, label=f'{label} c/c₀')
        else:
            ax_counts.plot(time_ps, counts, 's--', color=color, markersize=5,
                           linewidth=1.0, label=f'{label} (count)')

    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('c(t) / c(0) (consumed species)')
    ax_counts.set_ylabel('Count (produced species)')
    # 2° amine is consumed AND produced (1° -> 2° -> 3°), so its c/c0 can
    # exceed 1 while the amine pool shifts; do not clip it.
    ax.set_ylim(0, ratio_max * 1.05)
    ax_counts.set_ylim(bottom=0)
    ax.set_title('Species concentrations during curing (Fig. 5)')

    handles1, labels1 = ax.get_legend_handles_labels()
    handles2, labels2 = ax_counts.get_legend_handles_labels()
    ax.legend(handles1 + handles2, labels1 + labels2, loc='center left',
              fontsize=8)

    fig.tight_layout()
    for fmt in ('png', 'pdf'):
        fig.savefig(output_dir / f'species_concentrations.{fmt}', dpi=150)
    plt.close(fig)
    print(f'Saved species_concentrations.png/.pdf to {output_dir}')


def plot_gel_curve(
    topology_path: Path,
    species: list[str],
    output_dir: Path,
    f: int = 2,
    g: int = 5,
    r: float = 1.0,
) -> None:
    """Gelation curve: largest monomer-component fraction vs epoxide conversion.

    x = epoxide conversion (1 - n_epoxide / n_epoxide_initial) per topology
    snapshot; y = fraction of monomers in the largest connected monomer-level
    component (percolation indicator).  A vertical dashed line marks the
    Flory-Stockmayer theoretical gel point alpha_gel = 1/sqrt(r(f-1)(g-1))
    (textbook baseline, not from the paper; defaults f=2 DGEBA, g=5 DETA,
    r=1 -> 0.5).  Design: specs/decisions.md 2026-07-09 Track 2 / E2 entry.
    """
    snapshots = load_topology_snapshots(topology_path)
    if not snapshots:
        print(f'No topology snapshots in {topology_path} -- skipping gel plot.')
        return

    # Monomer membership from the initial (cycle -1) snapshot, as in
    # plot_dpn_vs_conversion.
    monomer_sets = monomer_sets_from_bonds(snapshots[0][2])
    if not monomer_sets:
        print('No monomers recovered from initial topology -- skipping gel plot.')
        return

    n_epoxide_initial = len(find_epoxide_rings(snapshots[0][2], species))
    if n_epoxide_initial == 0:
        print('No epoxide rings in initial topology -- skipping gel plot.')
        return

    conversions: list[float] = []
    fractions: list[float] = []
    for _step, _cycle, bonds in snapshots:
        n_epoxide = len(find_epoxide_rings(bonds, species))
        conversions.append(1.0 - n_epoxide / n_epoxide_initial)
        fractions.append(largest_component_fraction(bonds, monomer_sets))

    alpha_gel = gel_point_flory_stockmayer(f, g, r)

    output_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(conversions, fractions, 'o-', color='tab:blue', markersize=6,
            linewidth=1.0, label='Largest component (measured)')
    ax.axvline(alpha_gel, color='k', linestyle='--', linewidth=1.0,
               label=f'Flory-Stockmayer α_gel = {alpha_gel:.3f}')

    ax.set_xlabel('Epoxide conversion α')
    ax.set_ylabel('Largest monomer-component fraction')
    ax.set_title('Gelation: network percolation vs conversion')
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.05)
    ax.legend(loc='upper left', fontsize=8)
    fig.tight_layout()
    for fmt in ('png', 'pdf'):
        fig.savefig(output_dir / f'gel_curve.{fmt}', dpi=150)
    plt.close(fig)
    print(f'Saved gel_curve.png/.pdf to {output_dir}')


def _species_from_summary(summary_path: Path) -> list[str] | None:
    """Rebuild the run's species list from summary.json (primary path).

    run_epoxy_amine.py records epoxy_smiles/amine_smiles/n_epoxies/n_amines in
    summary.json; the layout order is epoxies first then amines with RDKit
    seeds (seed, seed + 1), matching build_epoxy_amine_system.  The seed is
    read from summary.json if present, else from manifest.json next to it
    (run_epoxy_amine writes it there); element ordering from SMILES + AddHs is
    actually seed-independent, so a missing seed only degrades determinism
    bookkeeping, not the rebuilt species — we fall back to 42 with a warning.
    Requires RDKit (returns None with a message if unavailable).
    """
    try:
        summary = json.loads(summary_path.read_text(encoding='utf-8'))
    except (OSError, json.JSONDecodeError) as exc:
        print(f'Could not read summary {summary_path}: {exc}')
        return None

    required = ('epoxy_smiles', 'amine_smiles', 'n_epoxies', 'n_amines')
    if any(key not in summary for key in required):
        print(f'{summary_path} lacks {required} -- cannot rebuild species.')
        return None

    seed = summary.get('seed')
    if seed is None:
        manifest_path = summary_path.parent / 'manifest.json'
        if manifest_path.exists():
            try:
                manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
                seed = manifest.get('seed', (manifest.get('extra') or {}).get('seed'))
            except (OSError, json.JSONDecodeError):
                seed = None
    if seed is None:
        print('WARNING: no seed in summary/manifest; using 42 (species '
              'ordering is seed-independent, see _species_from_summary).')
        seed = 42

    try:
        from scripts._systems import layout_species
    except ImportError:
        from _systems import layout_species  # script-dir sys.path fallback

    try:
        return layout_species([
            (summary['epoxy_smiles'], int(summary['n_epoxies']), int(seed)),
            (summary['amine_smiles'], int(summary['n_amines']), int(seed) + 1),
        ])
    except ImportError as exc:
        print(f'RDKit unavailable -- cannot rebuild species from summary: {exc}')
        return None


def _resolve_species(
    species_json: Path | None,
    summary_path: Path | None,
    trajectory_path: Path | None,
) -> list[str] | None:
    """Species list for the network-analysis figures.

    Priority: explicit --species-json > rebuild from --summary (primary path,
    needs RDKit) > 'species' recorded in the trajectory.jsonl header (written
    by run_epoxy_amine trajectory output; no RDKit needed).
    """
    if species_json is not None:
        return list(json.loads(species_json.read_text(encoding='utf-8')))

    if summary_path is not None and summary_path.exists():
        species = _species_from_summary(summary_path)
        if species is not None:
            return species

    if trajectory_path is not None and trajectory_path.exists():
        with open(trajectory_path, 'r', encoding='utf-8') as fh:
            first = json.loads(fh.readline())
        if first.get('_header') and first.get('species'):
            print(f'Species taken from trajectory header: {trajectory_path}')
            return list(first['species'])

    print('No species source available (--species-json / --summary rebuild / '
          'trajectory header) -- skipping species figures.')
    return None


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


def _infer_production_start_step(manifest_path: Path) -> int | None:
    """Best-effort production-onset step from a run's manifest.json.

    Priority (A4): ``extra['production_start_step']`` — the exact post-equilibration
    step the workflow recorded — is used first. Otherwise falls back to summing the
    pre-production step counts recorded in ``extra`` (equilibration and, when
    present, activation). Returns None if none of these are available so the caller
    can warn and fall back to 0. The --production-start-step CLI arg is the reliable
    override when equilibration was run outside the workflow manifest.
    """
    import json
    if not manifest_path.exists():
        return None
    try:
        data = json.loads(manifest_path.read_text(encoding='utf-8'))
    except (OSError, json.JSONDecodeError):
        return None
    extra = data.get('extra') or {}
    # The workflow records the true production onset directly; prefer it.
    recorded = extra.get('production_start_step')
    if recorded is not None:
        return int(recorded)
    equil = extra.get('equil_steps')
    activation = extra.get('activation_steps')
    if equil is None and activation is None:
        return None
    return int(equil or 0) + int(activation or 0)


def _infer_timestep_fs(manifest_path: Path) -> float | None:
    """Best-effort MD timestep (fs) from a run's manifest.json.

    Reads ``extra['timestep_fs']`` — the workflow serializes PolymerizationConfig
    (which carries ``timestep_fs``) into the manifest's ``extra`` block. Returns
    None when the manifest is missing or lacks the field, so the caller can warn
    and fall back to a default. Mirrors :func:`_infer_production_start_step`.
    """
    import json
    if not manifest_path.exists():
        return None
    try:
        data = json.loads(manifest_path.read_text(encoding='utf-8'))
    except (OSError, json.JSONDecodeError):
        return None
    extra = data.get('extra') or {}
    ts = extra.get('timestep_fs')
    if ts is None:
        return None
    return float(ts)


def _resolve_timestep_fs(
    cli_timestep_fs: float | None,
    trajectory: Path | None,
    bonds: Path | None,
    topology: Path | None,
) -> float:
    """Single source of truth for the MD timestep used by the time-axis figures.

    Priority (approved 2026-07-30): explicit ``--timestep-fs`` > manifest.json
    ``extra['timestep_fs']`` (looked up next to bonds/trajectory/topology) >
    1.0 fs. The last case emits a warning because the ps time axis (Eq. 11 fit,
    species traces) would be inaccurate.
    """
    if cli_timestep_fs is not None:
        return cli_timestep_fs
    manifest_dir = None
    for candidate in (bonds, trajectory, topology):
        if candidate is not None:
            manifest_dir = candidate.parent
            break
    inferred = (_infer_timestep_fs(manifest_dir / 'manifest.json')
                if manifest_dir is not None else None)
    if inferred is None:
        print(
            'WARNING: --timestep-fs not given and could not be inferred from '
            'manifest.json; using 1.0 fs. The time axis (ps) may be inaccurate '
            '-- pass --timestep-fs with the run\'s actual value.'
        )
        return 1.0
    return inferred


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
    parser.add_argument('--timestep-fs', type=float, default=None,
                        help='MD timestep in fs (for the Eq. 11 fit and species '
                             'time axis). If omitted, inferred from manifest.json '
                             "(extra['timestep_fs']); if that is unavailable, "
                             'defaults to 1.0 with a warning.')
    parser.add_argument('--production-start-step', type=int, default=None,
                        dest='production_start_step',
                        help='Global step where production (the cycle loop) begins; '
                             'the Eq. 11 fit measures k_p from t=0 there (L8/A4). '
                             'If omitted, inferred from manifest.json '
                             '(equil_steps + activation_steps), else 0.')
    parser.add_argument('--cell-xy-area', type=float, default=None,
                        help='Cross-sectional area (Å²) for density normalization. '
                             'If omitted, computed from the trajectory cell '
                             '(mean Lx*Ly); errors if no cell was recorded.')
    parser.add_argument('--topology', type=Path, default=None,
                        help='topology.jsonl for the Carothers Fig. 4c measured-DPn '
                             'plot (step-growth/nylon only). If omitted, looked up '
                             'next to --bonds.')
    parser.add_argument('--summary', type=Path, nargs='+', default=None,
                        help='summary.json file(s) for S2 diagnostics plot.')
    parser.add_argument('--species-figures', action='store_true',
                        dest='species_figures',
                        help='Generate network-analysis figures (species '
                             'concentrations, gel curve) from --topology. '
                             'Species come from --species-json, or are rebuilt '
                             'from the first --summary (epoxy-amine runs; '
                             'needs RDKit), or read from the trajectory header.')
    parser.add_argument('--species-json', type=Path, default=None,
                        dest='species_json',
                        help='JSON file containing the species list (element '
                             'symbols per atom) for --species-figures. '
                             'Overrides the summary-based rebuild.')
    parser.add_argument('--output-dir', type=Path, default=Path('runs/smoke/figures'))
    args = parser.parse_args()

    if args.trajectory is None and args.summary is None:
        parser.error('At least one of --trajectory or --summary is required.')

    if args.trajectory is not None and not args.trajectory.exists():
        print(f'Trajectory file not found: {args.trajectory}')
        return

    # Resolve the MD timestep once (CLI > manifest > 1.0) so the Eq. 11 fit and
    # the species time axis share one authoritative value (approved 2026-07-30).
    timestep_fs = _resolve_timestep_fs(
        args.timestep_fs, args.trajectory, args.bonds, args.topology,
    )

    if args.trajectory is not None:
        plot_energy_vs_step(args.trajectory, args.output_dir)
        plot_temperature_vs_step(args.trajectory, args.output_dir, args.target_temperature)

    if args.trajectory is not None and args.bonds is not None:
        header, frames = read_trajectory(args.trajectory)

        # Sample alpha(t) on the recorded frame steps rather than every MD step
        # (A6): turns a several-hundred-thousand-iteration loop into one pass over
        # the (few thousand) sampled frames.
        step_range = None
        if frames:
            step_range = np.array(sorted({f.step for f in frames}), dtype=np.int64)

        # Production onset for the Eq. 11 fit (L8/A4): CLI arg wins; else infer
        # from the run manifest; else 0 with a warning.
        production_start = args.production_start_step
        if production_start is None:
            inferred = _infer_production_start_step(args.bonds.parent / 'manifest.json')
            if inferred is None:
                print(
                    'WARNING: could not infer --production-start-step from '
                    'manifest.json (no equil_steps/activation_steps); using 0. '
                    'k_p may be under-estimated if the run had an '
                    'equilibration/activation lead-in.'
                )
                production_start = 0
            else:
                production_start = inferred

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
                               timestep_fs=timestep_fs,
                               production_start_step=production_start,
                               step_range=step_range)
        plot_density_profile(
            args.bonds, args.trajectory, args.output_dir,
            cell_xy_area=args.cell_xy_area,
        )

        # Fig. 4c (Carothers measured DPn vs conversion): step-growth only. Gate
        # on the run summary's carothers_p marker so chain-growth (vinyl) runs —
        # which also emit topology.jsonl but do not obey DPn = 1/(1-p) — no-op.
        topo_path = args.topology or (args.bonds.parent / 'topology.jsonl')
        summary_path = args.bonds.parent / 'summary.json'
        is_step_growth = False
        if summary_path.exists():
            try:
                is_step_growth = 'carothers_p' in json.loads(
                    summary_path.read_text(encoding='utf-8'))
            except (OSError, json.JSONDecodeError):
                is_step_growth = False
        want_dpn = (is_step_growth or args.topology is not None) and topo_path.exists()
        if want_dpn:
            # Carothers DPn = 1/(1-p) assumes LINEAR step-growth. A branched
            # network (epoxy-amine, functionality f>2) has monomers bonded to
            # more than two neighbours; the theory does not apply, so skip the
            # plot with a clear warning rather than emit a misleading figure
            # (approved 2026-07-30). Linear systems (nylon) are unaffected.
            snaps = load_topology_snapshots(topo_path)
            max_deg = (
                max_inter_monomer_degree(
                    snaps[-1][2], monomer_sets_from_bonds(snaps[0][2]))
                if snaps else 0
            )
            if max_deg > 2:
                print(
                    f'WARNING: branching detected (max inter-monomer bonds per '
                    f'monomer = {max_deg}); Carothers DPn plot assumes linear '
                    f'step-growth -- skipped.'
                )
            else:
                plot_dpn_vs_conversion(
                    args.bonds, topo_path, n_sites, args.output_dir)

    if args.species_figures:
        if args.topology is None or not args.topology.exists():
            parser.error('--species-figures requires --topology '
                         '(topology.jsonl of an epoxy-amine run).')
        species = _resolve_species(
            args.species_json,
            args.summary[0] if args.summary else None,
            args.trajectory,
        )
        if species is not None:
            plot_species_concentrations(
                args.topology, species, args.output_dir,
                timestep_fs=timestep_fs,
            )
            plot_gel_curve(args.topology, species, args.output_dir)

    if args.summary:
        # S2 diagnostics needs per-cycle 'logs'; epoxy-amine summaries carry
        # them too, so --summary serves both plots.
        plot_s2_diagnostics(args.summary, args.output_dir)


if __name__ == '__main__':
    main()
