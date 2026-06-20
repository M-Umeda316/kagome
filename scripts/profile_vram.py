"""Peak-VRAM profiler for paper-scale organic systems on OrbMol-v2 (CUDA).

Purpose
-------
Measure the *sustained-MD* peak GPU VRAM of each organic reproduction target so
the required card size (e.g. RTX 5000 Ada = 32 GB) can be confirmed by
measurement rather than extrapolation. This addresses the open item in
specs/decisions.md "2026-06-20: 再現スコープを有機系に限定" (follow-up b) and the
VRAM-ceiling analysis in specs/decisions.md "2026-06-15: 200+10 exceeds 16 GB".

What it measures
----------------
For each system it (1) builds at paper density (direct placement, or dilute
placement + ML compression like the production run scripts), (2) optionally
relaxes close contacts with a short FIRE minimize, then (3) runs a short Langevin
MD on CUDA with the *real* OrbMol-v2 backend. The sustained-MD regime is the one
that caused the historical 16 GB hang (per-step neighbour-graph size variation
fragments the allocator), so it is what the verdict is based on.

It records two peaks:
  * torch reserved peak   — torch.cuda.max_memory_reserved (allocator-internal)
  * device peak           — nvidia-smi polled device memory (incl. CUDA context,
                            cuBLAS, nvalchemiops PME kernels) when available.

Standalone & one-way friendly
-----------------------------
Self-contained: needs only this repo + the kagome-gpu env. Writes a JSON report
and prints a table, so results can be read on the workstation without sending
anything back.

OOM is a valid result: if a system exceeds VRAM the run is caught, recorded as
'OOM', memory is freed, and the next system still runs.

Usage (on the RTX 5000 Ada workstation, inside kagome-gpu, from repo root)
--------------------------------------------------------------------------
NOTE: run as a module (-m) so the ``scripts`` package resolves under the
editable install; ``python scripts/profile_vram.py`` will NOT find scripts._systems.

    python -m scripts.profile_vram --device cuda --output-dir runs/vram_profile

    # quick smoke test at small counts first (cheap, validates the script):
    python -m scripts.profile_vram --device cuda --scale 0.25 --md-steps 30 \
        --output-dir runs/vram_profile_smoke

    # only the worst case (VRAM ceiling):
    python -m scripts.profile_vram --systems nylon66

Paper anchors: SI S-3/S-4 (vinyl 0.5 g/mL, 333 K; nylon 0.5 g/mL, 300 K),
Table S1/S2. See specs/decisions.md for the density/temperature rationale.
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import subprocess
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path

# Must be set before torch/orb import so the profiled conditions match
# production (run_vinyl_aibn.py / orb_backend.py set the same vars). Expandable
# segments defragments the per-step neighbour-graph allocations; this is a no-op
# on Windows-native CUDA but DOES work on Linux/WSL. See decisions.md 2026-06-15.
os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')
os.environ.setdefault('TORCHDYNAMO_DISABLE', '1')
os.environ.setdefault('PYTORCH_CUDA_ALLOC_CONF', 'expandable_segments:True')

import numpy as np

from scripts._systems import (
    _DIACID_SMILES,
    _DIAMINE_SMILES,
    _INITIATOR_SMILES,
    _METHACRYLATE_SMILES,
    _MONOMER_SMILES,
    box_from_density,
    build_nylon66_system,
    build_vinyl_aibn_system,
)
from kagome.integrators.init_velocities import maxwell_boltzmann_velocities
from kagome.integrators.langevin import LangevinIntegrator, LangevinParams
from kagome.workflows.polymerization import masses_from_species

logging.basicConfig(level=logging.INFO, format='%(name)s | %(message)s')
logger = logging.getLogger('profile_vram')

_BYTES_PER_GB = 1024 ** 3


@dataclass
class SystemSpec:
    """One paper-scale organic system to profile."""

    name: str
    kind: str                       # 'vinyl' | 'nylon'
    counts: dict[str, int]          # builder counts (scaled by --scale)
    density_g_per_ml: float
    temperature_K: float
    monomer_smiles: str | None = None   # vinyl only; None -> methyl acrylate


# Paper-scale organic targets (epoxy/CuO is out of scope — OrbMol-v2 domain is
# organic only; see decisions.md 2026-06-20). The vinyl monomer is varied to
# bracket the atom-count range: methyl acrylate (smallest) → styrene (largest
# buildable vinyl). Since 2026-06-20 the builder also handles 1,1-disubstituted
# vinyls (methacrylate, diphenylethylene, dimethyl itaconate) whose propagating
# radical is tertiary; methacrylate is included as a representative tertiary
# case. The largest organic target is still nylon-6,6 (~4400 atoms), which sets
# the VRAM ceiling.
PAPER_SYSTEMS: list[SystemSpec] = [
    SystemSpec(
        name='vinyl_methyl_acrylate',
        kind='vinyl',
        counts={'n_monomers': 200, 'n_initiators': 10},
        density_g_per_ml=0.5,
        temperature_K=333.0,
        monomer_smiles=_MONOMER_SMILES,           # 'C=CC(=O)OC' (~2520 atoms)
    ),
    SystemSpec(
        name='vinyl_methacrylate',
        kind='vinyl',
        counts={'n_monomers': 200, 'n_initiators': 10},
        density_g_per_ml=0.5,
        temperature_K=333.0,
        monomer_smiles=_METHACRYLATE_SMILES,      # methyl methacrylate (tertiary radical)
    ),
    SystemSpec(
        name='vinyl_styrene',
        kind='vinyl',
        counts={'n_monomers': 200, 'n_initiators': 10},
        density_g_per_ml=0.5,
        temperature_K=333.0,
        monomer_smiles='C=Cc1ccccc1',             # styrene (~3320 atoms)
    ),
    SystemSpec(
        name='nylon66',
        kind='nylon',
        counts={'n_diamines': 100, 'n_diacids': 100},
        density_g_per_ml=0.5,
        temperature_K=300.0,                      # ~4400 atoms (VRAM ceiling)
    ),
]


class _DevicePoller:
    """Background thread sampling device VRAM via nvidia-smi (peak, in GB).

    Cross-checks the torch allocator number with the true device footprint
    (CUDA context + cuBLAS + nvalchemiops kernels are not in torch reserved).
    Degrades to None if nvidia-smi is unavailable (e.g. some WSL2 setups)."""

    def __init__(self, gpu_index: int, interval_s: float = 0.25) -> None:
        self._gpu_index = gpu_index
        self._interval = interval_s
        self._stop = threading.Event()
        self._peak_mb = 0.0
        self._ok = True
        self._thread: threading.Thread | None = None

    def _query(self) -> float | None:
        try:
            out = subprocess.run(
                ['nvidia-smi', f'--id={self._gpu_index}',
                 '--query-gpu=memory.used', '--format=csv,noheader,nounits'],
                capture_output=True, text=True, timeout=5,
            )
            if out.returncode != 0:
                return None
            return float(out.stdout.strip().splitlines()[0])
        except (OSError, ValueError, IndexError, subprocess.TimeoutExpired):
            return None

    def _run(self) -> None:
        while not self._stop.is_set():
            val = self._query()
            if val is None:
                self._ok = False
                return
            self._peak_mb = max(self._peak_mb, val)
            self._stop.wait(self._interval)

    def start(self) -> None:
        if self._query() is None:
            self._ok = False
            return
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def stop_and_peak_gb(self) -> float | None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=2)
        if not self._ok:
            return None
        return self._peak_mb / 1024.0


def _build_system(spec: SystemSpec, args, calc):
    """Build positions/species/cell at paper density, mirroring the production
    run scripts: try direct placement, else dilute placement + compression.

    Compression uses the classical FF by default (--compress-backend classical),
    so densification does not consume MLIP GPU time and the measured GPU peak
    isolates the MD regime (decision 2026-06-20). Pass --compress-backend ml to
    measure the all-ML upper bound instead."""
    scale, seed, no_compress = args.scale, args.seed, args.no_compress
    rng = np.random.default_rng(seed)

    if spec.kind == 'vinyl':
        n_mono = max(1, int(round(spec.counts['n_monomers'] * scale)))
        n_init = max(1, int(round(spec.counts['n_initiators'] * scale)))
        counts_for_density = {spec.monomer_smiles: n_mono, _INITIATOR_SMILES: n_init}
        target_edge = box_from_density(counts_for_density, spec.density_g_per_ml)
        # Placement order: initiators (seed) then monomers (seed+1) — must match
        # the classical-FF topology built for compression.
        molecule_specs = (
            (_INITIATOR_SMILES, n_init, seed),
            (spec.monomer_smiles, n_mono, seed + 1),
        )

        def _build(edge: float, gen):
            pos, sp, *_ = build_vinyl_aibn_system(
                n_monomers=n_mono, n_initiators=n_init, box_size=edge, rng=gen,
                monomer_smiles=spec.monomer_smiles, initiator_smiles=_INITIATOR_SMILES,
                rdkit_seed=seed,
            )
            return pos, sp
    else:  # nylon
        n_dia = max(1, int(round(spec.counts['n_diamines'] * scale)))
        n_acid = max(1, int(round(spec.counts['n_diacids'] * scale)))
        counts_for_density = {_DIAMINE_SMILES: n_dia, _DIACID_SMILES: n_acid}
        target_edge = box_from_density(counts_for_density, spec.density_g_per_ml)
        molecule_specs = (
            (_DIAMINE_SMILES, n_dia, seed),
            (_DIACID_SMILES, n_acid, seed + 1),
        )

        def _build(edge: float, gen):
            pos, sp, *_ = build_nylon66_system(
                n_diamines=n_dia, n_diacids=n_acid, box_size=edge, rng=gen,
                rdkit_seed=seed,
            )
            return pos, sp

    # Try direct placement at paper density first.
    try:
        positions, species = _build(target_edge, rng)
        cell = np.diag([target_edge, target_edge, target_edge]).astype(float)
        logger.info('  built directly at %.2f g/mL (box %.1f Å)',
                    spec.density_g_per_ml, target_edge)
        return positions, species, cell
    except RuntimeError:
        logger.info('  direct placement failed at %.2f g/mL; placing dilute then compressing',
                    spec.density_g_per_ml)

    # Fallback: dilute placement, then compress to target (ML-driven, GPU).
    place_edge = None
    for place_density in (0.25, 0.20, 0.15, 0.10):
        edge = box_from_density(counts_for_density, place_density)
        if edge <= target_edge:
            continue
        try:
            positions, species = _build(edge, np.random.default_rng(seed))
            place_edge = edge
            break
        except RuntimeError:
            continue
    if place_edge is None:
        raise RuntimeError(f'{spec.name}: could not place even at 0.10 g/mL')

    if no_compress:
        achieved = box_from_density(counts_for_density, spec.density_g_per_ml)
        logger.warning(
            '  --no-compress: profiling at dilute box %.1f Å (NOT paper density %.2f g/mL); '
            'VRAM is an UNDER-estimate', place_edge, spec.density_g_per_ml)
        cell = np.diag([place_edge, place_edge, place_edge]).astype(float)
        return positions, species, cell

    from kagome.backends.classical_backend import make_compress_calculator
    from kagome.integrators.minimize import compress_box
    from kagome.prep.openmm_equilibrate import MoleculeSpec

    specs = [MoleculeSpec(smi, n, rdkit_seed=s) for smi, n, s in molecule_specs]
    compress_calc = make_compress_calculator(
        args.compress_backend, specs, calc, platform=args.compress_platform,
        target_edge_A=target_edge,
    )
    place_cell = np.diag([place_edge, place_edge, place_edge]).astype(float)
    logger.info('  compressing box %.1f -> %.1f Å on %s (compress-backend=%s)...',
                place_edge, target_edge, compress_calc.name, args.compress_backend)
    result = compress_box(positions, place_cell, target_edge, species, compress_calc)
    return result.positions, species, result.cell


def _run_md_probe(positions, velocities, species, cell, masses, calc,
                  integrator, md_steps, dt_fs, rng):
    """Run md_steps of Langevin MD; one OrbMol-v2 call per step (production
    pattern). Atoms move, so the neighbour graph varies per step — the regime
    that drove the historical VRAM creep."""
    _, forces = calc.compute(positions, species, cell=cell)
    t0 = time.perf_counter()
    for _ in range(md_steps):
        integrator.pre_force(positions, velocities, forces, masses, dt_fs, rng, cell)
        _, forces = calc.compute(positions, species, cell=cell)
        integrator.post_force(velocities, forces, masses, dt_fs)
    elapsed = time.perf_counter() - t0
    return elapsed / max(1, md_steps)


def _profile_one(spec, args, calc, integrator, total_vram_gb):
    import torch

    result: dict = {
        'name': spec.name,
        'kind': spec.kind,
        'density_g_per_ml': spec.density_g_per_ml,
        'temperature_K': spec.temperature_K,
    }
    poller = _DevicePoller(args.gpu_index)
    poller.start()
    try:
        torch.cuda.empty_cache()
        torch.cuda.reset_peak_memory_stats()

        positions, species, cell = _build_system(spec, args, calc)
        n_atoms = len(species)
        result['n_atoms'] = n_atoms
        logger.info('  %s: %d atoms', spec.name, n_atoms)

        build_peak_gb = torch.cuda.max_memory_reserved() / _BYTES_PER_GB

        # Optional short relax of close contacts (production minimizes pre-TDBB).
        # fire_minimize does NOT mutate in place — it returns relaxed positions.
        if args.minimize_steps > 0:
            from kagome.integrators.minimize import FireParams, fire_minimize
            min_res = fire_minimize(
                positions, species, cell, calc,
                params=FireParams(fmax_kcal_mol_A=args.minimize_fmax,
                                  max_steps=args.minimize_steps),
            )
            positions = min_res.positions

        # Reset so the reported md peak isolates the sustained-MD regime.
        torch.cuda.reset_peak_memory_stats()
        rng = np.random.default_rng(args.seed)
        masses = masses_from_species(species)
        velocities = maxwell_boltzmann_velocities(masses, spec.temperature_K, rng)
        s_per_step = _run_md_probe(
            positions, velocities, species, cell, masses, calc, integrator,
            args.md_steps, args.timestep_fs, rng)

        md_reserved_gb = torch.cuda.max_memory_reserved() / _BYTES_PER_GB
        md_alloc_gb = torch.cuda.max_memory_allocated() / _BYTES_PER_GB
        device_peak_gb = poller.stop_and_peak_gb()

        # Verdict basis: device peak if available (it includes CUDA context +
        # PME kernels), else torch reserved + a context margin.
        if device_peak_gb is not None:
            peak_gb = max(device_peak_gb, md_reserved_gb, build_peak_gb)
            peak_source = 'nvidia-smi device'
        else:
            peak_gb = max(md_reserved_gb, build_peak_gb) + args.context_margin_gb
            peak_source = f'torch reserved + {args.context_margin_gb:.1f} GB margin'

        budget = args.vram_budget_gb if args.vram_budget_gb > 0 else total_vram_gb
        headroom = budget - peak_gb
        fits = headroom >= args.headroom_gb

        result.update({
            'status': 'ok',
            'torch_md_reserved_gb': round(md_reserved_gb, 2),
            'torch_md_allocated_gb': round(md_alloc_gb, 2),
            'torch_build_reserved_gb': round(build_peak_gb, 2),
            'device_peak_gb': None if device_peak_gb is None else round(device_peak_gb, 2),
            'peak_gb': round(peak_gb, 2),
            'peak_source': peak_source,
            'budget_gb': round(budget, 2),
            'headroom_gb': round(headroom, 2),
            'fits': fits,
            'sec_per_step': round(s_per_step, 3),
        })
        logger.info('  -> peak %.1f GB (%s); budget %.0f GB; headroom %.1f GB; fits=%s',
                    peak_gb, peak_source, budget, headroom, fits)
        return result

    except (torch.cuda.OutOfMemoryError, RuntimeError, MemoryError) as exc:
        poller.stop_and_peak_gb()
        msg = str(exc)
        is_oom = isinstance(exc, (torch.cuda.OutOfMemoryError, MemoryError)) \
            or 'out of memory' in msg.lower()
        result.update({
            'status': 'oom' if is_oom else 'error',
            'fits': False,
            'error': msg[:300],
        })
        logger.warning('  -> %s: %s', 'OOM' if is_oom else 'ERROR', msg[:160])
        return result
    finally:
        try:
            torch.cuda.empty_cache()
        except Exception:
            pass


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--device', default='cuda', help='cuda (required for VRAM) or cpu')
    parser.add_argument('--gpu-index', type=int, default=0,
                        help='nvidia-smi GPU index for device-peak polling')
    parser.add_argument('--output-dir', type=Path, default=Path('runs/vram_profile'))
    parser.add_argument('--systems', nargs='*', default=None,
                        help='subset of system names (default: all). '
                             f'Choices: {[s.name for s in PAPER_SYSTEMS]}')
    parser.add_argument('--scale', type=float, default=1.0,
                        help='scale molecule counts (e.g. 0.5 for a cheap smoke test)')
    parser.add_argument('--md-steps', type=int, default=100,
                        help='Langevin MD steps to probe the sustained regime')
    parser.add_argument('--timestep-fs', type=float, default=1.0,
                        help='MD timestep (fs); 1.0 matches S6 production')
    parser.add_argument('--minimize-steps', type=int, default=50,
                        help='capped FIRE relax of close contacts before MD (0 disables)')
    parser.add_argument('--minimize-fmax', type=float, default=10.0)
    parser.add_argument('--no-compress', action='store_true',
                        help='skip compression; profile at dilute density '
                             '(faster but UNDER-estimates VRAM)')
    parser.add_argument('--compress-backend', default='classical',
                        choices=['classical', 'ml'],
                        help='Calculator for box compression to paper density. '
                             '"classical" (default) uses OpenMM/OpenFF Sage so the GPU peak '
                             'isolates the MD regime (production uses classical prep too); '
                             '"ml" measures the all-ML upper bound. The MD probe always uses '
                             'OrbMol-v2.')
    parser.add_argument('--compress-platform', default='CPU',
                        choices=['CPU', 'CUDA', 'OpenCL', 'Reference'],
                        help='OpenMM platform for --compress-backend classical (default CPU, '
                             'keeps the GPU free so the device-peak poll reflects MD only).')
    parser.add_argument('--vram-budget-gb', type=float, default=0.0,
                        help='target card VRAM (GB); 0 = auto-detect device total. '
                             'RTX 5000 Ada = 32.')
    parser.add_argument('--headroom-gb', type=float, default=3.0,
                        help='required free headroom for a "fits" verdict')
    parser.add_argument('--context-margin-gb', type=float, default=1.5,
                        help='CUDA-context margin added to torch reserved when '
                             'nvidia-smi device peak is unavailable')
    parser.add_argument('--seed', type=int, default=7)
    args = parser.parse_args()

    if args.device != 'cuda':
        logger.warning('device=%s: VRAM numbers are only meaningful on cuda.', args.device)

    import torch
    if args.device == 'cuda' and not torch.cuda.is_available():
        raise SystemExit('CUDA not available; run on the GPU workstation.')

    total_vram_gb = 0.0
    if args.device == 'cuda':
        free_b, total_b = torch.cuda.mem_get_info()
        total_vram_gb = total_b / _BYTES_PER_GB
        logger.info('GPU: %s | total VRAM %.1f GB | free %.1f GB',
                    torch.cuda.get_device_name(0), total_vram_gb, free_b / _BYTES_PER_GB)

    selected = PAPER_SYSTEMS
    if args.systems:
        wanted = set(args.systems)
        selected = [s for s in PAPER_SYSTEMS if s.name in wanted]
        unknown = wanted - {s.name for s in PAPER_SYSTEMS}
        if unknown:
            raise SystemExit(f'unknown system(s): {sorted(unknown)}')

    logger.info('Creating OrbMol-v2 backend on %s...', args.device)
    from kagome.backends.orb_backend import create_orb_calculator
    calc = create_orb_calculator(device=args.device)
    integrator = LangevinIntegrator(LangevinParams(temperature_K=300.0))

    results = []
    for spec in selected:
        logger.info('=== %s (scale=%.2f) ===', spec.name, args.scale)
        results.append(_profile_one(spec, args, calc, integrator, total_vram_gb))

    report = {
        'gpu_name': (torch.cuda.get_device_name(0) if args.device == 'cuda' else 'cpu'),
        'total_vram_gb': round(total_vram_gb, 2),
        'vram_budget_gb': args.vram_budget_gb or round(total_vram_gb, 2),
        'scale': args.scale,
        'md_steps': args.md_steps,
        'timestep_fs': args.timestep_fs,
        'no_compress': args.no_compress,
        'compress_backend': args.compress_backend,
        'seed': args.seed,
        'env': {k: os.environ.get(k) for k in
                ('PYTORCH_CUDA_ALLOC_CONF', 'TORCHDYNAMO_DISABLE')},
        'systems': results,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.output_dir / 'vram_profile.json'
    out_path.write_text(json.dumps(report, indent=2), encoding='utf-8')

    # Console summary table.
    print('\n=== VRAM profile summary ===')
    print(f'GPU: {report["gpu_name"]}  total {total_vram_gb:.1f} GB  '
          f'budget {report["vram_budget_gb"]} GB  scale {args.scale}')
    header = f'{"system":26} {"atoms":>6} {"peak GB":>8} {"headroom":>9} {"fits":>6}  source'
    print(header)
    print('-' * len(header))
    for r in results:
        if r.get('status') == 'ok':
            print(f'{r["name"]:26} {r["n_atoms"]:>6} {r["peak_gb"]:>8.1f} '
                  f'{r["headroom_gb"]:>9.1f} {str(r["fits"]):>6}  {r["peak_source"]}')
        else:
            atoms = r.get('n_atoms', '?')
            print(f'{r["name"]:26} {str(atoms):>6} {"-":>8} {"-":>9} '
                  f'{r.get("status", "?").upper():>6}  {r.get("error", "")[:50]}')
    print(f'\nReport: {out_path}')


if __name__ == '__main__':
    main()
