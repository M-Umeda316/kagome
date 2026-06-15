"""Radical vinyl polymerization with methyl acrylate + AIBN initiator.

System:  n_monomers methyl acrylate (C=CC(=O)OC) +
         n_initiators isobutyronitrile radical model (CC(C)C#N)
Backend: OrbMol-v2 (default) or MACE-MP-0 via --backend mace
Ensemble: NPT (Langevin + MC barostat) or NVT with --no-barostat

Design:  specs/decisions.md — "T-G1: vinyl radical polymerization system"
Paper:   arXiv:2511.22874, Section 3, Table S1

Usage:
    python scripts/run_vinyl_aibn.py --seed 7 --output-dir runs/vinyl_aibn
    python scripts/run_vinyl_aibn.py --seed 7 --backend mace --output-dir runs/vinyl_aibn_mace
    python scripts/run_vinyl_aibn.py --seed 7 --no-barostat --output-dir runs/vinyl_aibn_nvt
"""
from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path

os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')
# Long paper-scale runs creep up in VRAM because the OrbMol-v2 neighbour graph
# size varies per step (atoms move), fragmenting the CUDA caching allocator until
# it exhausts memory and the run hangs. Expandable segments defragment this at no
# per-step cost. Must be set before torch is imported (orb import happens lazily
# in _create_backend). See specs/decisions.md 2026-06-15 VRAM record.
os.environ.setdefault('PYTORCH_CUDA_ALLOC_CONF', 'expandable_segments:True')

import numpy as np

from scripts._systems import build_vinyl_aibn_system
from src.backends.base import Calculator
from src.boost.tdbb import TDBBParams
from src.integrators.init_velocities import maxwell_boltzmann_velocities
from src.integrators.langevin import LangevinIntegrator, LangevinParams
from src.integrators.mc_barostat import MCBarostat, MCBarostatParams
from src.reactive.bonds import BondTracker
from src.workflows.polymerization import (
    PolymerizationConfig,
    PolymerizationWorkflow,
    SimulationState,
    masses_from_species,
)

logging.basicConfig(level=logging.INFO, format='%(name)s | %(message)s')
logger = logging.getLogger(__name__)


def _create_backend(backend: str, device: str, model: str, spin: int = 1) -> Calculator:
    if backend == 'orb':
        from src.backends.orb_backend import create_orb_calculator
        return create_orb_calculator(device=device, spin=spin)
    else:
        from src.backends.mace_backend import create_mace_calculator
        return create_mace_calculator(model=model, device=device)


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Radical vinyl polymerization: methyl acrylate + AIBN'
    )
    parser.add_argument('--seed', type=int, default=7)
    parser.add_argument('--output-dir', type=Path, default=Path('runs/vinyl_aibn'))
    parser.add_argument('--n-monomers', type=int, default=8)
    parser.add_argument('--n-initiators', type=int, default=2)
    parser.add_argument('--n-cycles', type=int, default=3)
    parser.add_argument('--biased-steps', type=int, default=500)
    parser.add_argument('--unbiased-steps', type=int, default=500)
    parser.add_argument('--box-size', type=float, default=None,
                        help='Box edge (Å). If omitted, computed from --density.')
    parser.add_argument('--density', type=float, default=0.5,
                        help='Initial density (g/mL). Paper SI S-3 uses 0.5 for vinyl. '
                             'Used only when --box-size is omitted.')
    parser.add_argument('--temperature', type=float, default=333.0)
    parser.add_argument('--pressure', type=float, default=1.0,
                        help='Target pressure (atm). Default 1.0 (assumed, not stated in paper).')
    parser.add_argument('--no-barostat', action='store_true',
                        help='Disable NPT barostat and run NVT instead.')
    parser.add_argument('--backend', type=str, default='orb',
                        choices=['orb', 'mace'],
                        help='MLIP backend (default: orb = OrbMol-v2)')
    parser.add_argument('--device', type=str, default='cpu')
    parser.add_argument('--model', type=str, default='small',
                        help='MACE model size (only used with --backend mace)')
    parser.add_argument('--initiator-smiles', type=str, default=None,
                        help='Override the initiator SMILES (e.g. "C[C](C)C#N" for the '
                             'real open-shell 2-cyanoprop-2-yl radical). Default: closed-shell model.')
    parser.add_argument('--spin', type=int, default=1,
                        help='Total spin multiplicity (2S+1) passed to OrbMol-v2. '
                             'Use 2 (doublet) for a single radical. Default 1 (singlet).')
    parser.add_argument('--minimize', dest='minimize', action='store_true', default=True,
                        help='FIRE energy minimization before TDBB (default: on). '
                             'Relaxes initial close contacts (paper anchor PDF p.20).')
    parser.add_argument('--no-minimize', dest='minimize', action='store_false',
                        help='Skip pre-TDBB energy minimization.')
    parser.add_argument('--minimize-fmax', type=float, default=1.0,
                        help='FIRE convergence threshold (kcal/mol/Å). Default 1.0.')
    parser.add_argument('--equil-steps', type=int, default=2000,
                        help='Unbiased NPT equilibration steps before TDBB '
                             '(paper anchor PDF p.20; length not specified, default 2000 '
                             '= 500 fs matching a TDBB block). 0 disables.')
    parser.add_argument('--timestep-fs', type=float, default=0.25,
                        help='MD timestep (fs). Default 0.25 fs (conservative, validated '
                             'for FIRE densification + ML NVT). 1.0 fs is standard for '
                             'organic ML MD and gives 4x speed for the same physical time.')
    parser.add_argument('--load-structure', type=Path, default=None,
                        help='Load a classically pre-equilibrated structure (JSON from '
                             'scripts/prep_structure.py) and skip build/place/compress. '
                             'positions+cell come from the file; the short ML re-equil '
                             '(--minimize/--equil-steps) still runs. See decision D-4.')
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    from scripts._systems import (
        _INITIATOR_SMILES,
        _MONOMER_SMILES,
        box_from_density,
    )
    counts = {_MONOMER_SMILES: args.n_monomers, _INITIATOR_SMILES: args.n_initiators}

    if args.box_size is not None:
        target_edge = args.box_size
    else:
        target_edge = box_from_density(counts, args.density)
        logger.info(
            'Box edge from density %.2f g/mL: %.2f Å (paper SI S-3)',
            args.density, target_edge,
        )

    # Backend is created before the system so it can drive FIRE densification.
    calc = _create_backend(args.backend, args.device, args.model, spin=args.spin)
    logger.info('Backend: %s (spin=%d)', calc.name, args.spin)

    _init_smiles = args.initiator_smiles or _INITIATOR_SMILES

    def _build(edge: float, gen: np.random.Generator):
        return build_vinyl_aibn_system(
            n_monomers=args.n_monomers,
            n_initiators=args.n_initiators,
            box_size=edge,
            rng=gen,
            initiator_smiles=_init_smiles,
        )

    if args.load_structure is not None:
        # Consume a classically pre-equilibrated structure (positions + cell).
        # template/groups/propagation_map are composition-derived (not position-
        # dependent), so rebuild them at a dilute box (always places) and overlay
        # the loaded coordinates. Atom order must match 1:1 (decision D-4),
        # asserted via the species list.
        from src.prep.structure_io import PreparedStructure

        prepared = PreparedStructure.load(args.load_structure)
        meta_edge = box_from_density(counts, 0.10)
        _, ref_species, template, groups, propagation_map = _build(meta_edge, rng)
        if list(ref_species) != list(prepared.species):
            raise ValueError(
                'Loaded structure species do not match the builder '
                f'(loaded N={len(prepared.species)}, builder N={len(ref_species)}). '
                'Ensure --n-monomers/--n-initiators match the prepped structure.'
            )
        positions = prepared.positions
        species = prepared.species
        cell = (prepared.cell if prepared.cell is not None
                else np.diag([target_edge, target_edge, target_edge]))
        logger.info(
            'Loaded pre-equilibrated structure from %s (%d atoms, box %.2f Å).',
            args.load_structure, len(species), float(cell[0, 0]),
        )
    else:
        logger.info(
            'Building vinyl/AIBN system: %d monomers + %d initiators, target box %.1f Å...',
            args.n_monomers, args.n_initiators, target_edge,
        )
        try:
            # Direct placement at the target box (small systems take this path,
            # bit-for-bit identical to prior behaviour).
            positions, species, template, groups, propagation_map = _build(target_edge, rng)
            cell = np.diag([target_edge, target_edge, target_edge])
        except RuntimeError:
            # Greedy placer stalls at high molecule counts even when the density is
            # physically feasible.  Place dilute, then FIRE-compress to the target
            # (see specs/decisions.md 2026-06-14 densification record).
            logger.warning(
                'Direct placement at %.2f Å failed — placing dilute then compressing.',
                target_edge,
            )
            place_edge = None
            for place_density in (0.25, 0.20, 0.15, 0.10):
                edge = box_from_density(counts, place_density)
                if edge <= target_edge:
                    continue
                try:
                    positions, species, template, groups, propagation_map = _build(edge, rng)
                    place_edge = edge
                    logger.info(
                        'Placed at dilute density %.2f g/mL (box %.2f Å); compressing to %.2f Å.',
                        place_density, edge, target_edge,
                    )
                    break
                except RuntimeError:
                    continue
            if place_edge is None:
                raise RuntimeError(
                    'Could not place the system even at dilute density 0.10 g/mL.'
                )
            from src.integrators.minimize import compress_box
            place_cell = np.diag([place_edge, place_edge, place_edge])
            result = compress_box(positions, place_cell, target_edge, species, calc)
            positions, cell = result.positions, result.cell

    logger.info(
        'System: %d atoms total  (%d radical_C, %d vinyl_alpha_C sites), box %.2f Å',
        len(species),
        len(groups['radical_C'].atom_indices),
        len(groups['vinyl_alpha_C'].atom_indices),
        float(cell[0, 0]),
    )
    logger.info('Propagation map: %d entries', len(propagation_map))

    langevin_params = LangevinParams(temperature_K=args.temperature)
    config = PolymerizationConfig(
        timestep_fs=args.timestep_fs,
        biased_steps=args.biased_steps,
        unbiased_steps=args.unbiased_steps,
        n_cycles=args.n_cycles,
        tdbb=TDBBParams(
            f2=10.0,
            gamma=1.0,
            f1_max_formation=250.0,
            f1_max_dissociation=125.0,
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
    barostat = None if args.no_barostat else MCBarostat(
        MCBarostatParams(pressure_atm=args.pressure, frequency=25)
    )
    if barostat:
        logger.info('NPT barostat enabled: P=%.2f atm', args.pressure)
    else:
        logger.info('Barostat disabled — running NVT.')

    masses = masses_from_species(species)
    velocities = maxwell_boltzmann_velocities(masses, args.temperature, rng)

    state = SimulationState(
        positions=positions,
        velocities=velocities,
        species=species,
        cell=cell,
        masses=masses,
    )

    logger.info(
        'Pre-TDBB: minimize=%s (fmax=%.2f), equilibration=%d steps',
        config.minimize, config.minimize_fmax, config.equil_steps,
    )
    logger.info(
        'Starting TDBB: %d cycles × (%d biased + %d unbiased steps), T=%.0f K',
        config.n_cycles, config.biased_steps, config.unbiased_steps, args.temperature,
    )

    wf = PolymerizationWorkflow(
        config, calc, template, groups,
        integrator=integrator,
        bond_tracker=tracker,
        barostat=barostat,
        propagation_map=propagation_map,
        propagation_target_group='radical_C',
    )
    logs = wf.run(
        state,
        output_dir=args.output_dir,
        config_path='configs/boost/paper_faithful.yaml',
    )

    n_form = len(tracker.confirmed_formations())
    n_dissoc = len(tracker.confirmed_dissociations())
    logger.info('Confirmed formations: %d, dissociations: %d', n_form, n_dissoc)

    n_reactive_sites = len(groups['radical_C'].atom_indices) + len(groups['vinyl_alpha_C'].atom_indices)
    summary = {
        'total_steps': state.step,
        'n_monomers': args.n_monomers,
        'n_initiators': args.n_initiators,
        'n_atoms': len(species),
        'box_size_A': float(cell[0, 0]),
        'cell_periodic': True,
        'backend': calc.name,
        'temperature_K': langevin_params.temperature_K,
        'biased_steps': args.biased_steps,
        'unbiased_steps': args.unbiased_steps,
        'n_cycles': args.n_cycles,
        'minimize': args.minimize,
        'minimize_fmax': args.minimize_fmax,
        'equil_steps': args.equil_steps,
        'confirmed_formations': n_form,
        'confirmed_dissociations': n_dissoc,
        'propagation_events': n_form,  # each formation triggers one propagation
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

    n_reactive = args.n_monomers * 2 + args.n_initiators
    print('\nTo generate figures:')
    print(
        f'  python scripts/reproduce_figures.py '
        f'--trajectory {args.output_dir}/trajectory.jsonl '
        f'--bonds {args.output_dir}/bonds.jsonl '
        f'--n-reactive-sites {n_reactive} '
        f'--target-temperature {args.temperature} '
        f'--output-dir {args.output_dir}/figures'
    )


if __name__ == '__main__':
    main()
