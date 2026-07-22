"""Scale-up task (b): edge budget vs memory/speed/force-accuracy on OrbMol-v2.

Context (specs/decisions.md 追補 2026-07-22, item 4(b))
------------------------------------------------------
Activation memory for the autograd-based force pass scales with the edge
count, which is set by the graph parameters radius=6.0 Å /
max_num_neighbors=120 (orb_models pretrained.py defaults, kagome does not
override). This spike measures, on the paper-scale methyl-acrylate system:

  1. The *actual* per-atom neighbour-count distribution at 6.0 Å — if the
     120 cap never binds, lowering the cap frees zero memory (edges are a
     real COO list, not padded; see forcefield_adapter.py).
  2. A sweep over (radius, max_num_neighbors): edge count, per-call peak
     activation memory, sec/call, and force deviation vs the (6.0, 120)
     baseline. Radius is a TRAINING-time parameter (adapter
     is_compatible_with asserts on it), so any reduction is a measured
     accuracy trade-off, never a silent default change.

This script MEASURES ONLY. Production defaults are unchanged; adopting any
reduced setting requires a decisions.md entry + user approval (CLAUDE.md
ask-first: simplifications that change scientific meaning).

Usage (WSL, pfpoly-gpu, repo root):
    PYTHONPATH=.:src python -m scripts.spike_edge_budget \
        --device cuda --output-dir runs/scaleup_b
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import time
from pathlib import Path

os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')
os.environ.setdefault('TORCHDYNAMO_DISABLE', '1')
os.environ.setdefault('PYTORCH_CUDA_ALLOC_CONF', 'expandable_segments:True')

import numpy as np

from scripts.profile_vram import PAPER_SYSTEMS, _build_system
from kagome.integrators.init_velocities import maxwell_boltzmann_velocities
from kagome.integrators.langevin import LangevinIntegrator, LangevinParams
from kagome.workflows.polymerization import masses_from_species

logging.basicConfig(level=logging.INFO, format='%(name)s | %(message)s')
logger = logging.getLogger('spike_edge_budget')

_BYTES_PER_GB = 1024 ** 3

# (radius Å, max_num_neighbors). First entry is the baseline all deltas are
# measured against — the orb_models orbmol_v2 defaults.
SWEEP: list[tuple[float, int]] = [
    (6.0, 120),
    (6.0, 80),
    (6.0, 60),
    (6.0, 40),
    (5.5, 120),
    (5.0, 120),
    (4.5, 120),
    (4.0, 120),
]


def _neighbor_stats(calc, positions, species, cell) -> dict:
    """Per-atom neighbour-count distribution under the CURRENT calc adapter."""
    import ase
    import torch

    atoms = ase.Atoms(symbols=species, positions=positions, cell=cell,
                      pbc=cell is not None)
    atoms.info['charge'] = calc._charge
    atoms.info['spin'] = calc._spin
    batch = calc._adapter.from_ase_atoms(atoms, device=calc._device)
    receivers = batch.receivers.detach().cpu()
    senders = batch.senders.detach().cpu()
    n_atoms = len(species)
    # The cap's direction (in- vs out-degree) is an upstream implementation
    # detail (knn_alchemi), so count both; 'counts' = the tighter-capped side.
    in_deg = torch.bincount(receivers, minlength=n_atoms).numpy()
    out_deg = torch.bincount(senders, minlength=n_atoms).numpy()
    cap = int(calc._adapter.max_num_neighbors)
    counts = in_deg if in_deg.max() >= out_deg.max() else out_deg
    stats = {
        'n_edges': int(receivers.numel()),
        'in_deg_max': int(in_deg.max()),
        'out_deg_max': int(out_deg.max()),
        'neighbors_mean': round(float(counts.mean()), 2),
        'neighbors_median': int(np.median(counts)),
        'neighbors_p95': int(np.percentile(counts, 95)),
        'neighbors_max': int(counts.max()),
        'cap': cap,
        'atoms_at_cap': int((counts >= cap).sum()),
        'atoms_at_cap_pct': round(100.0 * float((counts >= cap).mean()), 2),
    }
    del batch
    return stats


def _measure_config(calc, positions, species, cell, repeats: int) -> dict:
    """Timed force calls + per-call peak activation memory for one adapter."""
    import torch

    # Warm-up (kernel compilation, allocator growth) excluded from stats.
    calc.compute(positions, species, cell=cell)
    torch.cuda.synchronize()
    torch.cuda.reset_peak_memory_stats()
    t0 = time.perf_counter()
    for _ in range(repeats):
        energy, forces = calc.compute(positions, species, cell=cell)
    torch.cuda.synchronize()
    elapsed = (time.perf_counter() - t0) / repeats
    return {
        'sec_per_call': round(elapsed, 4),
        'peak_alloc_gb': round(torch.cuda.max_memory_allocated() / _BYTES_PER_GB, 3),
        'peak_reserved_gb': round(torch.cuda.max_memory_reserved() / _BYTES_PER_GB, 3),
        'energy_kcal_mol': energy,
        'forces': forces,
    }


def _force_deltas(f_base: np.ndarray, f_test: np.ndarray) -> dict:
    """Force deviation vs baseline in kcal/mol/Å, absolute and relative."""
    diff = f_test - f_base
    rmse = float(np.sqrt(np.mean(diff ** 2)))
    f_rms = float(np.sqrt(np.mean(f_base ** 2)))
    return {
        'force_rmse_kcal_mol_A': round(rmse, 4),
        'force_maxabs_kcal_mol_A': round(float(np.abs(diff).max()), 4),
        'force_rmse_rel_pct': round(100.0 * rmse / f_rms, 3),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--device', default='cuda')
    parser.add_argument('--output-dir', type=Path, default=Path('runs/scaleup_b'))
    parser.add_argument('--system', default='vinyl_methyl_acrylate',
                        choices=[s.name for s in PAPER_SYSTEMS])
    parser.add_argument('--scale', type=float, default=1.0)
    parser.add_argument('--repeats', type=int, default=10,
                        help='timed force calls per config (after 1 warm-up)')
    parser.add_argument('--md-warm-steps', type=int, default=50,
                        help='Langevin steps (baseline adapter) to produce a second, '
                             'thermally-perturbed evaluation snapshot (0 disables)')
    parser.add_argument('--minimize-steps', type=int, default=50)
    parser.add_argument('--minimize-fmax', type=float, default=10.0)
    parser.add_argument('--no-compress', action='store_true')
    parser.add_argument('--compress-backend', default='classical',
                        choices=['classical', 'ml'])
    parser.add_argument('--compress-platform', default='CPU',
                        choices=['CPU', 'CUDA', 'OpenCL', 'Reference'])
    parser.add_argument('--timestep-fs', type=float, default=1.0)
    parser.add_argument('--seed', type=int, default=7)
    args = parser.parse_args()

    import torch
    if args.device == 'cuda' and not torch.cuda.is_available():
        raise SystemExit('CUDA not available; run on the GPU workstation.')

    from orb_models.forcefield.forcefield_adapter import ForcefieldAtomsAdapter
    from kagome.backends.orb_backend import create_orb_calculator

    logger.info('Creating OrbMol-v2 backend on %s (empty_cache off: measuring '
                'allocator behaviour per config)...', args.device)
    calc = create_orb_calculator(device=args.device, empty_cache=False)
    base_adapter = calc._adapter
    spec = next(s for s in PAPER_SYSTEMS if s.name == args.system)

    positions, species, cell = _build_system(spec, args, calc)
    n_atoms = len(species)
    logger.info('%s: %d atoms', spec.name, n_atoms)

    if args.minimize_steps > 0:
        from kagome.integrators.minimize import FireParams, fire_minimize
        min_res = fire_minimize(
            positions, species, cell, calc,
            params=FireParams(fmax_kcal_mol_A=args.minimize_fmax,
                              max_steps=args.minimize_steps),
        )
        positions = min_res.positions

    # Evaluation snapshots: relaxed, plus (optionally) one after a short
    # baseline-adapter MD so deltas are not specific to a minimized geometry.
    snapshots: list[tuple[str, np.ndarray]] = [('relaxed', positions.copy())]
    if args.md_warm_steps > 0:
        rng = np.random.default_rng(args.seed)
        masses = masses_from_species(species)
        velocities = maxwell_boltzmann_velocities(masses, spec.temperature_K, rng)
        integrator = LangevinIntegrator(LangevinParams(temperature_K=spec.temperature_K))
        pos_md = positions.copy()
        _, forces = calc.compute(pos_md, species, cell=cell)
        for _ in range(args.md_warm_steps):
            integrator.pre_force(pos_md, velocities, forces, masses,
                                 args.timestep_fs, rng, cell)
            _, forces = calc.compute(pos_md, species, cell=cell)
            integrator.post_force(velocities, forces, masses, args.timestep_fs)
        snapshots.append((f'md{args.md_warm_steps}', pos_md))

    base_radius = float(base_adapter.radius)
    base_cap = int(base_adapter.max_num_neighbors)
    logger.info('Baseline adapter: radius=%.1f Å, max_num_neighbors=%d',
                base_radius, base_cap)
    neigh = _neighbor_stats(calc, snapshots[-1][1], species, cell)
    logger.info('Neighbour stats (baseline, %s snapshot): %s',
                snapshots[-1][0], neigh)

    configs = []
    baseline_forces: dict[str, np.ndarray] = {}
    for radius, cap in SWEEP:
        entry: dict = {'radius_A': radius, 'max_num_neighbors': cap,
                       'snapshots': {}}
        try:
            if (radius, cap) == (base_radius, base_cap):
                calc._adapter = base_adapter
            else:
                calc._adapter = ForcefieldAtomsAdapter(
                    radius=radius, max_num_neighbors=cap)
            for snap_name, snap_pos in snapshots:
                torch.cuda.empty_cache()
                m = _measure_config(calc, snap_pos, species, cell, args.repeats)
                forces = m.pop('forces')
                nstats = _neighbor_stats(calc, snap_pos, species, cell)
                m['n_edges'] = nstats['n_edges']
                m['atoms_at_cap_pct'] = nstats['atoms_at_cap_pct']
                if (radius, cap) == (base_radius, base_cap):
                    baseline_forces[snap_name] = forces
                    m.update({'force_rmse_kcal_mol_A': 0.0,
                              'force_maxabs_kcal_mol_A': 0.0,
                              'force_rmse_rel_pct': 0.0})
                else:
                    m.update(_force_deltas(baseline_forces[snap_name], forces))
                entry['snapshots'][snap_name] = m
            entry['status'] = 'ok'
            logger.info('radius=%.1f cap=%d: %s', radius, cap,
                        json.dumps(entry['snapshots'], default=str)[:200])
        except Exception as exc:  # record and continue — an assert inside the
            # model on graph params is itself a useful result.
            entry['status'] = 'error'
            entry['error'] = str(exc)[:300]
            logger.warning('radius=%.1f cap=%d FAILED: %s', radius, cap,
                           str(exc)[:160])
        configs.append(entry)
    calc._adapter = base_adapter

    report = {
        'system': spec.name,
        'n_atoms': n_atoms,
        'scale': args.scale,
        'seed': args.seed,
        'repeats': args.repeats,
        'snapshots': [name for name, _ in snapshots],
        'gpu_name': torch.cuda.get_device_name(0) if args.device == 'cuda' else 'cpu',
        'baseline': {'radius_A': base_radius, 'max_num_neighbors': base_cap},
        'neighbor_stats_baseline': neigh,
        'env': {k: os.environ.get(k) for k in
                ('PYTORCH_CUDA_ALLOC_CONF', 'TORCHDYNAMO_DISABLE')},
        'sweep': configs,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.output_dir / 'edge_budget.json'
    out_path.write_text(json.dumps(report, indent=2), encoding='utf-8')

    print('\n=== edge budget sweep (vs radius=6.0, cap=120) ===')
    print(f'{args.system}  {n_atoms} atoms  |  neighbours: mean '
          f"{neigh['neighbors_mean']} p95 {neigh['neighbors_p95']} "
          f"max {neigh['neighbors_max']} at-cap {neigh['atoms_at_cap_pct']}%")
    hdr = (f'{"radius":>6} {"cap":>4} {"snap":>8} {"edges":>9} {"GB":>6} '
           f'{"s/call":>7} {"F-RMSE":>8} {"F-max":>8} {"rel%":>6}')
    print(hdr)
    print('-' * len(hdr))
    for c in configs:
        if c['status'] != 'ok':
            print(f"{c['radius_A']:>6} {c['max_num_neighbors']:>4} "
                  f"{'ERROR':>8}  {c['error'][:60]}")
            continue
        for snap_name, m in c['snapshots'].items():
            print(f"{c['radius_A']:>6} {c['max_num_neighbors']:>4} "
                  f"{snap_name:>8} {m['n_edges']:>9} {m['peak_alloc_gb']:>6} "
                  f"{m['sec_per_call']:>7} {m['force_rmse_kcal_mol_A']:>8} "
                  f"{m['force_maxabs_kcal_mol_A']:>8} {m['force_rmse_rel_pct']:>6}")
    print(f'\nReport: {out_path}')


if __name__ == '__main__':
    main()
