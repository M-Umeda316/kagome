"""Classical (OpenMM/OpenFF) structure preparation — standalone entry point.

Runs in the dedicated ``pfpoly-prep`` conda environment (OpenFF + OpenMM), which
is kept separate from the ML production environment (``pfpoly-gpu``, OrbMol-v2)
because OpenFF pulls openff-nagl → a second PyTorch that would collide with the
production torch. This script writes the relaxed structure to a JSON file that
``scripts/run_vinyl_aibn.py --load-structure`` consumes. See specs/decisions.md
2026-06-14 "Decouple initial-structure preparation" (decision D-4).

Pipeline: dilute grid placement → classical minimize → compress to the target
(0.5 g/mL) density → NVT thermalization → save PreparedStructure.

Usage (from the repo root, in the prep env):
    conda run -n pfpoly-prep python scripts/prep_structure.py \
        --n-monomers 200 --n-initiators 10 --seed 42 \
        --output runs/prep/paper200_structure.json
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

# RDKit + OpenMM + MKL can each ship their own libiomp5md, and loading two
# triggers a hard abort ("OMP: Error #15") on Windows. Allow the duplicate
# (same workaround as scripts/run_vinyl_aibn.py). Must be set before numpy/
# rdkit/openmm import.
os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')

# Make the repo importable when run as a bare script in an env that does not
# have the project installed (the fresh prep env).
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import numpy as np

from scripts._systems import (
    _INITIATOR_SMILES,
    _MONOMER_SMILES,
    box_from_density,
    build_vinyl_aibn_system,
)
from src.prep.openmm_equilibrate import (
    ClassicalPrepConfig,
    MoleculeSpec,
    equilibrate_structure,
)

logging.basicConfig(level=logging.INFO, format='%(name)s | %(message)s')
logger = logging.getLogger(__name__)

# rdkit_seed offsets used by build_vinyl_aibn_system (initiator=seed, monomer=seed+1).
_RDKIT_SEED = 42


def _place_dilute(n_monomers, n_initiators, counts, target_edge, rng):
    """Place the system at the densest feasible dilute density above the target.

    Returns (positions, species, cell, place_edge). The greedy placer is reliable
    at low density; the classical stage then compresses to the target edge.
    """
    for place_density in (0.25, 0.20, 0.15, 0.10):
        edge = box_from_density(counts, place_density)
        if edge <= target_edge:
            continue
        try:
            positions, species, _, _, _ = build_vinyl_aibn_system(
                n_monomers=n_monomers, n_initiators=n_initiators,
                box_size=edge, rng=rng,
            )
            logger.info(
                'Placed at dilute density %.2f g/mL (box %.2f Å).',
                place_density, edge,
            )
            return positions, species, np.diag([edge, edge, edge]), edge
        except RuntimeError:
            continue
    raise RuntimeError('Could not place the system even at 0.10 g/mL.')


def main() -> None:
    parser = argparse.ArgumentParser(description='Classical structure prep (OpenMM/OpenFF)')
    parser.add_argument('--n-monomers', type=int, default=200)
    parser.add_argument('--n-initiators', type=int, default=10)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--target-density', type=float, default=0.5,
                        help='Target initial density (g/mL). Paper SI S-3: 0.5.')
    parser.add_argument('--temperature', type=float, default=333.0)
    parser.add_argument('--charge-method', choices=['gasteiger', 'nagl'],
                        default='gasteiger',
                        help='Partial-charge method. Default gasteiger (lightweight; '
                             'the structure is overwritten by the ML re-equil, see D-2).')
    parser.add_argument('--forcefield', type=str, default='openff-2.2.0.offxml')
    parser.add_argument('--compress-stages', type=int, default=20)
    parser.add_argument('--compress-relax-steps', type=int, default=200)
    parser.add_argument('--nvt-steps', type=int, default=50000)
    parser.add_argument('--platform', type=str, default='CPU',
                        choices=['CPU', 'CUDA', 'OpenCL', 'Reference'])
    parser.add_argument('--output', type=Path, required=True,
                        help='Output JSON path (PreparedStructure).')
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    counts = {_MONOMER_SMILES: args.n_monomers, _INITIATOR_SMILES: args.n_initiators}
    target_edge = box_from_density(counts, args.target_density)
    logger.info(
        'Target: %d monomers + %d initiators at %.2f g/mL (box %.2f Å).',
        args.n_monomers, args.n_initiators, args.target_density, target_edge,
    )

    positions, species, cell, place_edge = _place_dilute(
        args.n_monomers, args.n_initiators, counts, target_edge, rng,
    )

    # Placement order is initiators-first then monomers (build_vinyl_aibn_system),
    # so the MoleculeSpec list must follow that order for the OpenFF topology.
    molecule_specs = [
        MoleculeSpec(_INITIATOR_SMILES, args.n_initiators, rdkit_seed=_RDKIT_SEED),
        MoleculeSpec(_MONOMER_SMILES, args.n_monomers, rdkit_seed=_RDKIT_SEED + 1),
    ]

    cfg = ClassicalPrepConfig(
        target_density_g_per_ml=args.target_density,
        temperature_K=args.temperature,
        charge_method=args.charge_method,
        forcefield=args.forcefield,
        compress_stages=args.compress_stages,
        compress_relax_steps=args.compress_relax_steps,
        nvt_steps=args.nvt_steps,
        platform=args.platform,
        seed=args.seed,
        metadata={
            'n_monomers': args.n_monomers,
            'n_initiators': args.n_initiators,
            'place_edge_A': place_edge,
        },
    )

    prepared = equilibrate_structure(positions, species, cell, molecule_specs, cfg)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    prepared.save(args.output)
    logger.info(
        'Saved prepared structure: %d atoms, box %.2f Å -> %s',
        prepared.n_atoms, float(prepared.cell[0, 0]), args.output,
    )


if __name__ == '__main__':
    main()
