"""Demonstrate CHAIN PROPAGATION end-to-end (handoff-plan-v4 S1).

One open-shell radical (C[C](C)C#N) adds successive methyl-acrylate monomers; the
active radical migrates along the growing chain. Along a single chain the unpaired
electron count is always 1, so the run is a doublet (spin=2) throughout — this
sidesteps the multi-radical system-spin problem (S3) while exercising propagation.

To make propagation deterministic (the melt-driven version is S2), the script runs
one TDBB cycle at a time and, before each cycle, places the next unreacted monomer
productively next to the CURRENT chain-end radical (vinyl pi-face on the radical
lobe, alpha-C at ~2.5 A — just at the addition barrier from the PES scan). The
workflow then forms the bond (in-phase detection), the radical migrates to the
monomer's beta-C, and the next monomer is positioned against the new chain end.

DEMO devices (clearly non-paper-faithful, for mechanism demonstration only): the
candidate window is widened (--select-rmin) so the near-contact pair is selectable,
and monomers are directed into place. Production window stays [3,6]; melt sampling
is S2. TDBB params (f2, f1_max, r0) are unchanged.

Usage:
    python scripts/demo_chain_propagation.py --device cuda --n-monomers 4
"""
from __future__ import annotations

import argparse
import logging
import os
from dataclasses import asdict
from pathlib import Path

os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')

import numpy as np

from scripts._systems import build_vinyl_aibn_system
from scripts.scan_radical_addition import _rotation_align, _unit
from kagome.boost.tdbb import TDBBParams
from kagome.integrators.init_velocities import maxwell_boltzmann_velocities
from kagome.integrators.langevin import LangevinIntegrator, LangevinParams
from kagome.reactive.bonds import BondTracker
from kagome.workflows.polymerization import (
    PolymerizationConfig, PolymerizationWorkflow, SimulationState, masses_from_species,
)

logging.basicConfig(level=logging.INFO, format='%(name)s | %(message)s')
logger = logging.getLogger(__name__)

RADICAL_SMILES = 'C[C](C)C#N'
N_INIT_ATOMS = 11
_HEAVY = {'C', 'N', 'O'}


def _radical_lobe(pos, species, rc):
    """Approach normal at the chain-end radical: normal to the plane of its
    nearest (<=1.8 A) heavy neighbours."""
    nbrs = [i for i in range(len(species))
            if i != rc and species[i] in _HEAVY
            and np.linalg.norm(pos[i] - pos[rc]) < 1.8]
    nbrs = sorted(nbrs, key=lambda i: np.linalg.norm(pos[i] - pos[rc]))[:3]
    if len(nbrs) >= 2:
        n = np.cross(pos[nbrs[1]] - pos[rc], pos[nbrs[0]] - pos[rc])
        if np.linalg.norm(n) > 1e-6:
            return _unit(n)
    # fallback: away from neighbour centroid
    if nbrs:
        return _unit(pos[rc] - pos[np.array(nbrs)].mean(axis=0))
    return np.array([0.0, 0.0, 1.0])


def _place_monomer_at_radical(pos, species, rc, al, be, mono_atoms, approach):
    """Position one monomer so its vinyl alpha-C sits `approach` A from the chain-end
    radical on its lobe, pi-face toward the radical. Mutates pos in place."""
    lobe = _radical_lobe(pos, species, rc)
    # vinyl pi-face normal from alpha, beta and one H on alpha
    h_on_al = min((i for i in mono_atoms if species[i] == 'H'),
                  key=lambda i: np.linalg.norm(pos[i] - pos[al]))
    vnormal = _unit(np.cross(pos[be] - pos[al], pos[h_on_al] - pos[al]))
    idx = np.array(sorted(mono_atoms))
    R = _rotation_align(vnormal, -lobe)            # face the radical
    block = pos[idx] - pos[al]
    pos[idx] = block @ R.T + (pos[rc] + approach * lobe)
    return float(np.linalg.norm(pos[al] - pos[rc]))


def main() -> None:
    ap = argparse.ArgumentParser(description='TDBB chain-propagation demo')
    ap.add_argument('--device', default='cpu')
    ap.add_argument('--n-monomers', type=int, default=4)
    ap.add_argument('--approach', type=float, default=2.5)
    ap.add_argument('--select-rmin', type=float, default=1.5)
    ap.add_argument('--biased-steps', type=int, default=2000)
    ap.add_argument('--unbiased-steps', type=int, default=300)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--output-dir', default='runs/demo_chain_propagation')
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    pos, sp, template, groups, pmap, _ = build_vinyl_aibn_system(
        n_monomers=args.n_monomers, n_initiators=1, box_size=60.0, rng=rng,
        initiator_smiles=RADICAL_SMILES,
    )
    pos = np.array(pos, dtype=np.float64)
    for ps in template.pairs:
        if ps.is_formation:
            ps.r_min = args.select_rmin
    n_per_mono = (len(sp) - N_INIT_ATOMS) // args.n_monomers
    # alpha -> (its monomer's atom indices), in placement order (monomers after init)
    mono_atoms_by_alpha = {}
    for m in range(args.n_monomers):
        start = N_INIT_ATOMS + m * n_per_mono
        atoms = list(range(start, start + n_per_mono))
        al = next(a for a in groups['vinyl_alpha_C'].atom_indices if a in atoms)
        mono_atoms_by_alpha[al] = atoms

    from kagome.backends.orb_backend import create_orb_calculator
    calc = create_orb_calculator(device=args.device, spin=2)  # single chain -> doublet

    masses = masses_from_species(sp)
    vel = maxwell_boltzmann_velocities(masses, 333.0, rng)
    state = SimulationState(positions=pos, velocities=vel, species=sp,
                            cell=None, masses=masses)
    tracker = BondTracker()
    cfg = PolymerizationConfig(
        timestep_fs=0.25, biased_steps=args.biased_steps,
        unbiased_steps=args.unbiased_steps, n_cycles=1, seed=args.seed,
        save_interval=0, minimize=False, equil_steps=0,
        tdbb=TDBBParams(f2=10.0, gamma=1.0, f1_max_formation=250.0,
                        f1_max_dissociation=125.0, lambda_vdw=0.60),
    )
    wf = PolymerizationWorkflow(
        cfg, calc, template, groups,
        integrator=LangevinIntegrator(LangevinParams(temperature_K=333.0)),
        bond_tracker=tracker, barostat=None,
        propagation_map=pmap, propagation_target_group='radical_C',
    )

    reacted_alphas: set[int] = set()
    for step_i in range(args.n_monomers):
        radlist = groups['radical_C'].atom_indices
        assert len(radlist) == 1, f'single-chain invariant broken: radical_C={radlist}'
        rc = radlist[0]
        remaining = [a for a in groups['vinyl_alpha_C'].atom_indices
                     if a not in reacted_alphas and a in mono_atoms_by_alpha]
        if not remaining:
            break
        # pick the monomer whose alpha is currently closest (any works since we place it)
        al = remaining[0]
        be = pmap[al]
        d = _place_monomer_at_radical(state.positions, sp, rc, al, be,
                                      mono_atoms_by_alpha[al], args.approach)
        logger.info('Step %d: chain-end radical=%d, placing monomer alpha=%d at %.2f A',
                    step_i + 1, rc, al, d)
        wf.config.seed = args.seed + step_i  # vary thermal noise per segment
        before = len(tracker.confirmed_formations())
        wf.run(state)
        after = len(tracker.confirmed_formations())
        reacted_alphas.add(al)
        if after > before:
            logger.info('  -> formed (cumulative formations=%d)', after)
        else:
            logger.info('  -> no bond this segment; stopping chain.')
            break

    n_form = len(tracker.confirmed_formations())
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    tracker.save(out / 'bonds.jsonl')

    # Provenance (RF17): this demo drives wf.run() per segment without output_dir,
    # so emit a manifest here so the run is traceable to seed/backend/weights/config.
    from kagome.workflows.manifest import RunManifest, _normalize_value
    manifest_extra = _normalize_value(asdict(cfg))
    manifest_extra.update(
        backend=calc.name,
        model_id=getattr(calc, 'model_id', calc.name),
        n_monomers=args.n_monomers,
        n_initiators=1,
        demo='chain_propagation',
        select_rmin=args.select_rmin,
        approach=args.approach,
    )
    RunManifest(
        config_path='(demo: scripts/demo_chain_propagation.py)',
        seed=args.seed,
        backend=calc.name,
        output_dir=str(out),
        extra=manifest_extra,
    ).save(out / 'manifest.json')
    print(f'\nRESULT: confirmed_formations={n_form}  propagation_events={n_form}  '
          f'final chain-end radical={groups["radical_C"].atom_indices}')
    if n_form >= 2:
        print('SUCCESS: chain propagation demonstrated (>=2 successive additions).')
    else:
        print('Propagation not reached (<2 additions); try --approach 2.3 / more steps.')


if __name__ == '__main__':
    main()
