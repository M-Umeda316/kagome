"""Radical vinyl COPOLYMERIZATION: methyl acrylate + methyl methacrylate + AIBN.

System:  n_acrylate     methyl acrylate      (C=CC(=O)OC,     12 atoms)
       + n_methacrylate methyl methacrylate  (C=C(C)C(=O)OC,  15 atoms)
       + n_initiators   isobutyronitrile radical model (CC(C)C#N)
Backend: toy (cheap LJ smoke, GPU-free) | orb (OrbMol-v2) | mace | aimnet
Ensemble: NVT (default) or NPT with --barostat

Design:  specs/decisions.md — "2026-07-16: 共重合ビニル系ビルダー" (documented
         extension of the paper's single-monomer setup: both species' alpha-Cs
         share one vinyl_alpha_C group, no reactivity-ratio bias — sequence is
         left to the MLIP energetics).
Paper:   arXiv:2511.22874, Section 2/3, Table S1 (4-group ij+ik+jl criterion).

Usage:
    # cheap smoke (no GPU): toy backend, tiny system
    python scripts/run_vinyl_copolymer.py --backend toy --n-acrylate 3 \
        --n-methacrylate 3 --n-initiators 1 --n-cycles 2 --biased-steps 20 \
        --unbiased-steps 20 --output-dir runs/copoly_smoke

    # real MLIP, NVT
    python scripts/run_vinyl_copolymer.py --backend orb --seed 7 \
        --n-acrylate 8 --n-methacrylate 8 --n-initiators 2 \
        --output-dir runs/copoly_orb
"""
from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path

os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')
os.environ.setdefault('PYTORCH_CUDA_ALLOC_CONF', 'expandable_segments:True')

import numpy as np

from scripts._systems import (
    _INITIATOR_SMILES,
    _METHACRYLATE_SMILES,
    _MONOMER_SMILES,
    box_from_density,
    build_vinyl_copolymer_system,
    copolymer_initial_bonds,
)
from kagome.backends.base import Calculator
from kagome.boost.tdbb import TDBBParams
from kagome.integrators.init_velocities import maxwell_boltzmann_velocities
from kagome.integrators.langevin import LangevinIntegrator, LangevinParams
from kagome.integrators.mc_barostat import MCBarostat, MCBarostatParams
from kagome.reactive.bonds import BondTracker
from kagome.workflows.polymerization import (
    PolymerizationConfig,
    PolymerizationWorkflow,
    SimulationState,
    load_checkpoint,
    masses_from_species,
)

logging.basicConfig(level=logging.INFO, format='%(name)s | %(message)s')
logger = logging.getLogger(__name__)


def _create_backend(backend: str, device: str, model: str, spin: int = 1,
                    compile_model: bool = False,
                    empty_cache: bool = True) -> Calculator:
    if backend == 'toy':
        from kagome.backends.toy import ToyCalculator
        return ToyCalculator()
    if backend == 'orb':
        from kagome.backends.orb_backend import create_orb_calculator
        return create_orb_calculator(device=device, spin=spin,
                                     compile=compile_model,
                                     empty_cache=empty_cache)
    if backend == 'aimnet':
        from kagome.backends.aimnet_backend import create_aimnet_calculator
        aimnet_model = model if model and model.startswith('aimnet') else 'aimnet2-nse'
        return create_aimnet_calculator(model=aimnet_model, device=device, spin=spin)
    from kagome.backends.mace_backend import create_mace_calculator
    return create_mace_calculator(model=model, device=device)


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Vinyl copolymerization (acrylate + methacrylate) via TDBB.',
    )
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--seed', type=int, default=7)
    # composition
    parser.add_argument('--n-acrylate', type=int, default=8,
                        help='methyl acrylate monomers (C=CC(=O)OC).')
    parser.add_argument('--n-methacrylate', type=int, default=8,
                        help='methyl methacrylate monomers (C=C(C)C(=O)OC).')
    parser.add_argument('--n-initiators', type=int, default=2,
                        help='isobutyronitrile radical models (CC(C)C#N).')
    parser.add_argument('--n-cycles', type=int, default=10)
    # MD steps
    parser.add_argument('--biased-steps', type=int, default=2000)
    parser.add_argument('--unbiased-steps', type=int, default=1500)
    parser.add_argument('--equil-steps', type=int, default=2000)
    # 0.25 fs is REQUIRED for reactive multi-radical stability — 1.0 fs numerically
    # explodes the open-shell melt (decisions.md 2026-06-25, validity-domain §2.1).
    parser.add_argument('--timestep-fs', type=float, default=0.25)
    parser.add_argument('--minimize', dest='minimize', action='store_true',
                        default=True, help='FIRE-minimize before dynamics (default on).')
    parser.add_argument('--no-minimize', dest='minimize', action='store_false')
    parser.add_argument('--minimize-fmax', type=float, default=1.0)
    # box / density
    parser.add_argument('--box-size', type=float, default=None,
                        help='cubic box edge in Angstrom. Overrides --density.')
    parser.add_argument('--density', type=float, default=0.5,
                        help='target density g/mL (paper SI S-3).')
    # thermostat / ensemble
    parser.add_argument('--temperature', type=float, default=333.0)
    parser.add_argument('--friction-per-fs', type=float, default=0.01)
    parser.add_argument('--barostat', action='store_true',
                        help='enable NPT MC barostat (default NVT).')
    parser.add_argument('--pressure', type=float, default=1.0)
    # TDBB — defaults mirror the validated S6 vinyl recipe (decisions.md
    # 2026-06-26): f2=2 bridges the OrbMol capture dead-zone (paper's 10 leaves
    # ~0 bias force in the [3,6] Å window); f1_max 250/125 from the same recipe.
    parser.add_argument('--f2', type=float, default=2.0)
    parser.add_argument('--f1-max-formation', type=float, default=250.0)
    parser.add_argument('--f1-max-dissociation', type=float, default=125.0)
    parser.add_argument('--select-rmin', type=float, default=None)
    parser.add_argument('--select-rmax', type=float, default=None)
    # backend
    parser.add_argument('--backend', choices=['toy', 'orb', 'mace', 'aimnet'],
                        default='orb')
    parser.add_argument('--device', default='cuda')
    parser.add_argument('--model', default='')
    parser.add_argument('--compile', action='store_true')
    parser.add_argument('--no-empty-cache', action='store_false', dest='empty_cache')
    parser.add_argument('--spin', type=int, default=1)
    # checkpoint / resume
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--no-checkpoint', action='store_true')
    args = parser.parse_args()

    n_monomers = args.n_acrylate + args.n_methacrylate
    if n_monomers == 0:
        parser.error('Need at least one monomer (--n-acrylate / --n-methacrylate).')

    ckpt_file = args.output_dir / 'checkpoint.pkl'
    resuming = bool(args.resume and ckpt_file.exists())
    run_checkpoint_path = None if (args.no_checkpoint and not args.resume) else ckpt_file
    if args.resume and not ckpt_file.exists():
        logger.warning('--resume given but %s not found; starting fresh.', ckpt_file)

    rng = np.random.default_rng(args.seed)

    # monomer_specs order == placement order; keep density counts aggregated per
    # SMILES so box_from_density sums molar mass correctly.
    monomer_specs = [
        (_MONOMER_SMILES, args.n_acrylate),
        (_METHACRYLATE_SMILES, args.n_methacrylate),
    ]
    counts = {
        _MONOMER_SMILES: args.n_acrylate,
        _METHACRYLATE_SMILES: args.n_methacrylate,
        _INITIATOR_SMILES: args.n_initiators,
    }

    if args.box_size is not None:
        target_edge = args.box_size
    else:
        target_edge = box_from_density(counts, args.density)
        logger.info('Box edge from density %.2f g/mL: %.2f Å (paper SI S-3)',
                    args.density, target_edge)

    calc = _create_backend(args.backend, args.device, args.model, spin=args.spin,
                           compile_model=args.compile, empty_cache=args.empty_cache)
    logger.info('Backend: %s', calc.name)

    def _build(edge: float, gen: np.random.Generator):
        return build_vinyl_copolymer_system(
            monomer_specs=monomer_specs,
            n_initiators=args.n_initiators,
            box_size=edge,
            rng=gen,
            rdkit_seed=args.seed,
        )

    logger.info(
        'Building copolymer system: %d acrylate + %d methacrylate + %d initiator, '
        'target box %.1f Å...',
        args.n_acrylate, args.n_methacrylate, args.n_initiators, target_edge,
    )
    try:
        positions, species, template, groups, propagation_map, chain_c_map = _build(target_edge, rng)
        cell = np.diag([target_edge, target_edge, target_edge])
    except RuntimeError:
        logger.warning('Direct placement at %.2f Å failed — placing dilute then '
                       'compressing.', target_edge)
        place_edge = None
        for place_density in (0.25, 0.20, 0.15, 0.10):
            edge = box_from_density(counts, place_density)
            if edge <= target_edge:
                continue
            try:
                positions, species, template, groups, propagation_map, chain_c_map = _build(edge, rng)
                place_edge = edge
                break
            except RuntimeError:
                continue
        if place_edge is None:
            raise RuntimeError('Could not place the system even at 0.10 g/mL.')
        from kagome.integrators.minimize import compress_box
        from kagome.backends.classical_backend import make_compress_calculator
        from kagome.prep.openmm_equilibrate import MoleculeSpec
        # MoleculeSpec order/seeds must match the placement order in
        # build_vinyl_copolymer_system: initiators (rdkit_seed) first, then each
        # monomer spec (rdkit_seed + 1 + k). See specs/decisions.md 2026-07-16.
        specs = [MoleculeSpec(_INITIATOR_SMILES, args.n_initiators, rdkit_seed=args.seed)]
        for k, (smi, cnt) in enumerate(monomer_specs):
            specs.append(MoleculeSpec(smi, cnt, rdkit_seed=args.seed + 1 + k))
        # 'classical' uses OpenMM/OpenFF so densification skips MLIP GPU time;
        # the toy backend has no classical params, so reuse it via 'ml'.
        compress_calc = make_compress_calculator(
            'ml' if args.backend == 'toy' else 'classical', specs, calc,
            target_edge_A=target_edge,
        )
        place_cell = np.diag([place_edge, place_edge, place_edge])
        result = compress_box(positions, place_cell, target_edge, species, compress_calc)
        positions, cell = result.positions, result.cell

    if args.select_rmin is not None or args.select_rmax is not None:
        for ps in template.pairs:
            if args.select_rmin is not None:
                ps.r_min = args.select_rmin
            if args.select_rmax is not None:
                ps.r_max = args.select_rmax

    logger.info(
        'System: %d atoms  (%d radical_C, %d vinyl_alpha_C sites), box %.2f Å; '
        'propagation map %d entries',
        len(species), len(groups['radical_C'].atom_indices),
        len(groups['vinyl_alpha_C'].atom_indices), float(cell[0, 0]),
        len(propagation_map),
    )

    langevin_params = LangevinParams(
        temperature_K=args.temperature, friction_per_fs=args.friction_per_fs)
    config = PolymerizationConfig(
        timestep_fs=args.timestep_fs,
        biased_steps=args.biased_steps,
        unbiased_steps=args.unbiased_steps,
        n_cycles=args.n_cycles,
        tdbb=TDBBParams(
            f2=args.f2,
            gamma=1.0,
            f1_max_formation=args.f1_max_formation,
            f1_max_dissociation=args.f1_max_dissociation,
            lambda_vdw=0.60,
        ),
        seed=args.seed,
        save_interval=50,
        minimize=args.minimize,
        minimize_fmax=args.minimize_fmax,
        equil_steps=args.equil_steps,
    )

    integrator = LangevinIntegrator(langevin_params)
    tracker = BondTracker()
    barostat = MCBarostat(MCBarostatParams(pressure_atm=args.pressure, frequency=25)) \
        if args.barostat else None
    logger.info('Ensemble: %s', 'NPT (P=%.2f atm)' % args.pressure if barostat else 'NVT')

    masses = masses_from_species(species)
    velocities = maxwell_boltzmann_velocities(masses, args.temperature, rng)
    state = SimulationState(
        positions=positions, velocities=velocities, species=species,
        cell=cell, masses=masses,
    )

    # Trajectory bond-topology tracking: copolymer_initial_bonds builds the full
    # initial bond list with the SAME running-offset convention as the builder
    # (initiators first, then each monomer_specs entry in placement order), so
    # topology.jsonl is now emitted for the heterogeneous copolymer layout too
    # (specs/decisions.md 2026-07-17 "well-mixed 測定モード" 前提工事; previously
    # skipped here — 2026-07-16 known limitation). Retroactive reconstruction
    # (scripts/reconstruct_topology.py) is unrelated and remains single-monomer-
    # only; it is not needed for new runs since tracking is on from the start.
    initial_bonds = copolymer_initial_bonds(monomer_specs, args.n_initiators)
    wf = PolymerizationWorkflow(
        config, calc, template, groups,
        integrator=integrator,
        bond_tracker=tracker,
        barostat=barostat,
        propagation_map=propagation_map,
        propagation_target_group='radical_C',
        chain_c_map=chain_c_map,
        initial_bonds=initial_bonds,
    )

    if resuming:
        _ck = load_checkpoint(ckpt_file)
        _extra = _ck.get('extra', {}) or {}
        _spin = _extra.get('spin')
        if _spin is not None and getattr(calc, 'supports_spin', False):
            calc.set_spin(int(_spin))
        config = PolymerizationConfig(
            timestep_fs=config.timestep_fs,
            biased_steps=config.biased_steps,
            unbiased_steps=config.unbiased_steps,
            n_cycles=config.n_cycles,
            tdbb=config.tdbb,
            seed=config.seed,
            save_interval=config.save_interval,
            minimize=False,
            minimize_fmax=config.minimize_fmax,
            equil_steps=0,
        )
        wf.config = config

    logger.info('Starting TDBB: %d cycles × (%d biased + %d unbiased steps), T=%.0f K',
                config.n_cycles, config.biased_steps, config.unbiased_steps,
                args.temperature)

    logs = wf.run(
        state,
        output_dir=args.output_dir,
        config_path='configs/boost/paper_faithful.yaml',
        n_monomers=n_monomers,
        checkpoint_path=run_checkpoint_path,
        resume=resuming,
        checkpoint_extra={'spin': getattr(calc, '_spin', None)},
    )

    n_form = len(tracker.confirmed_formations())
    n_dissoc = len(tracker.confirmed_dissociations())
    logger.info('Confirmed formations: %d, dissociations: %d', n_form, n_dissoc)

    summary = {
        'total_steps': state.step,
        'n_acrylate': args.n_acrylate,
        'n_methacrylate': args.n_methacrylate,
        'n_monomers': n_monomers,
        'n_initiators': args.n_initiators,
        'n_atoms': len(species),
        'box_size_A': float(cell[0, 0]),
        'cell_periodic': True,
        'backend': calc.name,
        'temperature_K': langevin_params.temperature_K,
        'friction_per_fs': args.friction_per_fs,
        'ensemble': 'NPT' if barostat else 'NVT',
        'biased_steps': args.biased_steps,
        'unbiased_steps': args.unbiased_steps,
        'n_cycles': args.n_cycles,
        'minimize': args.minimize,
        'equil_steps': args.equil_steps,
        'confirmed_formations': n_form,
        'confirmed_dissociations': n_dissoc,
        'propagation_events': n_form,
        'logs': [
            {
                'cycle': log.cycle,
                'phase': log.phase,
                'steps': log.steps,
                'n_candidates': log.n_candidates,
                'n_selected': log.n_selected,
                'bias_energy': log.bias_energy,
                'min_pair_distance': (None if log.min_pair_distance == float('inf')
                                      else log.min_pair_distance),
            }
            for log in logs
        ],
    }
    out_path = args.output_dir / 'summary.json'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(summary, indent=2), encoding='utf-8')
    logger.info('Done. Results in %s', args.output_dir)

    print('\nTo generate figures:')
    print(
        f'  python scripts/reproduce_figures.py '
        f'--trajectory {args.output_dir}/trajectory.jsonl '
        f'--bonds {args.output_dir}/bonds.jsonl '
        f'--n-reactive-sites {n_monomers} '
        f'--target-temperature {args.temperature} '
        f'--output-dir {args.output_dir}/figures'
    )


if __name__ == '__main__':
    main()
