"""Demonstrate confirmed_formations >= 1 end-to-end through the TDBB workflow.

Places ONE open-shell 2-cyanoprop-2-yl radical (C[C](C)C#N) productively next to
one methyl acrylate monomer — radical carbon on the vinyl pi-face normal at ~2.5 A
(just at the addition barrier found by scripts/scan_radical_addition.py: TS ~2.2 A,
exothermic product ~1.54 A), radical lobe pointed at the terminal vinyl carbon —
then runs a biased TDBB segment. With the correct geometry, spin=doublet, and the
in-phase reaction detection, OrbMol-v2's radical-addition channel + the bias drive
the pair over the barrier to a confirmed C-C bond.

This isolates the chemistry/mechanism from the (separate) melt-sampling problem.

Usage:
    python scripts/demo_radical_formation.py --device cuda
"""
from __future__ import annotations

import argparse
import logging
import os

os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')

import numpy as np

from scripts._systems import _INITIATOR_SMILES, _MONOMER_SMILES, build_vinyl_aibn_system
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

RADICAL_SMILES = 'C[C](C)C#N'   # real open-shell 2-cyanoprop-2-yl radical
N_INIT_ATOMS = 11               # atoms in the radical block (placed first)


def main() -> None:
    ap = argparse.ArgumentParser(description='Demonstrate a TDBB radical-addition bond')
    ap.add_argument('--device', default='cpu')
    ap.add_argument('--approach', type=float, default=2.5,
                    help='initial radical_C--vinyl_alpha_C distance (A)')
    ap.add_argument('--biased-steps', type=int, default=2000)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--select-rmin', type=float, default=1.5,
                    help='DEMO-ONLY: widen the candidate window r_min so a near-contact '
                         'pair is selectable. The paper window is [3,6] (Table S1), which '
                         'excludes the <2.6 A bias-capture shell; the paper relies on '
                         'diffusion to deliver pairs inward during a long biased phase. '
                         'Here we select the pre-positioned near-contact pair directly to '
                         'demonstrate the workflow forms the bond. Set 3.0 for paper-faithful.')
    ap.add_argument('--output-dir', default='runs/demo_radical_formation')
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    pos, sp, template, groups, pmap, _ = build_vinyl_aibn_system(
        n_monomers=1, n_initiators=1, box_size=30.0, rng=rng,
        initiator_smiles=RADICAL_SMILES,
    )
    # DEMO-ONLY: widen the formation-pair selection window so the pre-positioned
    # near-contact pair is selected (see --select-rmin help).
    for ps in template.pairs:
        if ps.is_formation:
            ps.r_min = args.select_rmin
    logger.info('DEMO: formation candidate window set to [%.1f, %.1f] A (paper: [3,6]).',
                args.select_rmin, template.pairs[0].r_max)
    pos = np.array(pos, dtype=np.float64)
    rc = groups['radical_C'].atom_indices[0]
    al = groups['vinyl_alpha_C'].atom_indices[0]
    be = pmap[al]                       # beta carbon (other end of the C=C)
    mono = slice(N_INIT_ATOMS, len(sp))
    init = slice(0, N_INIT_ATOMS)

    # --- vinyl pi-face normal: from alpha-C, beta-C and one H on alpha-C ---
    h_mono = [i for i in range(*mono.indices(len(sp))) if sp[i] == 'H']
    h_on_al = min(h_mono, key=lambda i: np.linalg.norm(pos[i] - pos[al]))
    normal = _unit(np.cross(pos[be] - pos[al], pos[h_on_al] - pos[al]))

    # --- radical lobe: normal to the plane of its three carbon neighbours ---
    c_init = [i for i in range(*init.indices(len(sp))) if sp[i] == 'C' and i != rc]
    c_nbr = sorted(c_init, key=lambda i: np.linalg.norm(pos[i] - pos[rc]))[:3]
    lobe = _unit(np.cross(pos[c_nbr[1]] - pos[c_nbr[0]], pos[c_nbr[2]] - pos[c_nbr[0]]))

    # --- orient the radical so its lobe points at alpha-C (-normal), then place
    #     radical_C at alpha-C + approach*normal (pi-face approach) ---
    R = _rotation_align(lobe, -normal)
    block = pos[init] - pos[rc]
    pos[init] = block @ R.T            # rotate about radical_C (now at origin)
    target = pos[al] + args.approach * normal
    pos[init] += target                # radical_C -> target

    d0 = float(np.linalg.norm(pos[rc] - pos[al]))
    logger.info('Placed radical_C(%d) at %.2f A from vinyl_alpha_C(%d) on the pi-face.',
                rc, d0, al)

    from kagome.backends.orb_backend import create_orb_calculator
    calc = create_orb_calculator(device=args.device, spin=2)  # one radical -> doublet

    masses = masses_from_species(sp)
    vel = maxwell_boltzmann_velocities(masses, 333.0, rng)
    state = SimulationState(positions=pos, velocities=vel, species=sp,
                            cell=None, masses=masses)
    cfg = PolymerizationConfig(
        timestep_fs=0.25, biased_steps=args.biased_steps, unbiased_steps=500,
        n_cycles=1, seed=args.seed, save_interval=0, minimize=False, equil_steps=0,
        tdbb=TDBBParams(f2=10.0, gamma=1.0, f1_max_formation=250.0,
                        f1_max_dissociation=125.0, lambda_vdw=0.60),
    )
    tracker = BondTracker()
    wf = PolymerizationWorkflow(
        cfg, calc, template, groups,
        integrator=LangevinIntegrator(LangevinParams(temperature_K=333.0)),
        bond_tracker=tracker, barostat=None,
        propagation_map=pmap, propagation_target_group='radical_C',
    )
    from pathlib import Path
    wf.run(state, output_dir=Path(args.output_dir))

    d_final = float(np.linalg.norm(state.positions[rc] - state.positions[al]))
    n_form = len(tracker.confirmed_formations())
    print(f'\nRESULT: confirmed_formations={n_form}  '
          f'final r(radical_C,alpha)={d_final:.2f} A  (r0=2.04, bond if <2.04)')
    if n_form >= 1:
        print('SUCCESS: TDBB formed a radical-addition C-C bond end-to-end.')
    else:
        print('No bond this run; try --approach 2.3 or a longer --biased-steps.')


if __name__ == '__main__':
    main()
