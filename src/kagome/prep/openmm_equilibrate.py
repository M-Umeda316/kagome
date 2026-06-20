"""Classical (OpenMM + OpenFF) initial-structure preparation.

Reserves the ML potential (OrbMol-v2) for the TDBB reactive production by moving
the cheap structure-preparation work — packing relaxation, densification to the
paper's 0.5 g/mL initial density, and thermalization — onto a classical force
field (OpenFF Sage) run through OpenMM's optimized kernels.

Pipeline (decision D-3 "simple" protocol):
    minimize  →  compress dilute box to the target density  →  NVT thermalize
The box is held fixed at the paper density (decision D-1): we do NOT run a
classical barostat to the FF's own equilibrium density. Density evolution is left
to the ML NPT production loop.

Atom-order invariant (decision D-4): the OpenFF Topology is built from the SAME
RDKit Mol objects the system builder uses (``kagome.chem.builders._rdkit_mol``), in
the SAME molecule order, so the returned positions map 1:1 onto the caller's
``species`` / ``groups`` / ``propagation_map``.

All OpenMM / OpenFF imports are deferred to call time so this module imports
cleanly in environments that only run ML production. See specs/decisions.md
2026-06-14 "Decouple initial-structure preparation".
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from kagome.prep.structure_io import PreparedStructure
from kagome.units import ANGSTROM_PER_NM, NM_PER_ANGSTROM

logger = logging.getLogger(__name__)


@dataclass
class MoleculeSpec:
    """One distinct molecule species and how many copies the box contains.

    ``smiles`` and ``rdkit_seed`` must match what the system builder passed to
    ``_rdkit_mol`` so the OpenFF connectivity/ordering aligns with the placed
    coordinates. The specs MUST be listed in the same order the builder placed
    the molecules (for the vinyl system: initiators first, then monomers).
    """

    smiles: str
    count: int
    rdkit_seed: int = 42


@dataclass
class ClassicalPrepConfig:
    """Parameters for the classical structure-prep stage.

    Lengths/counts not fixed by the paper are engineering choices (decision D-3).
    """

    target_density_g_per_ml: float = 0.5      # paper initial density (SI S-3)
    temperature_K: float = 333.0              # production setpoint (PDF p.7)
    protocol: str = 'simple'                  # only 'simple' is implemented
    charge_method: str = 'nagl'              # 'nagl' | 'gasteiger'
    forcefield: str = 'openff-2.2.0.offxml'  # Sage 2.2
    nagl_model: str = 'openff-gnn-am1bcc-0.1.0-rc.3.pt'
    minimize_tolerance_kj_mol_nm: float = 10.0
    compress_stages: int = 20
    compress_minimize_iters: int = 200        # max minimizer iters per shrink stage
    compress_relax_steps: int = 200           # MD steps relaxing each shrink (after minimize)
    nvt_steps: int = 50_000                    # thermalization (classical, cheap)
    timestep_fs: float = 0.5                    # conservative: all-atom (H), dense box
    friction_per_ps: float = 1.0
    platform: str = 'CPU'                      # 'CUDA'|'OpenCL'|'CPU'|'Reference'
    seed: int = 42
    metadata: dict = field(default_factory=dict)


def _build_openff_topology(species: list[str], specs: list[MoleculeSpec], cfg):
    """Build an OpenFF Topology in builder order and verify the element order.

    Returns ``(topology, unique_offmols)``. Charges are assigned per unique
    molecule (nagl, or RDKit Gasteiger fallback) before replication so the
    downstream System carries partial charges without AmberTools/OpenEye.
    """
    from openff.toolkit import Molecule, Topology

    from kagome.chem.builders import _rdkit_mol

    unique: list = []
    ordered: list = []
    for spec in specs:
        rdmol = _rdkit_mol(spec.smiles, spec.rdkit_seed)
        offmol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)
        _assign_charges(offmol, rdmol, cfg)
        unique.append(offmol)
        ordered.extend([offmol] * spec.count)

    topology = Topology.from_molecules(ordered)

    top_symbols = [atom.symbol for atom in topology.atoms]
    if top_symbols != list(species):
        raise ValueError(
            'OpenFF topology atom order does not match the builder species '
            f'(topology N={len(top_symbols)}, species N={len(species)}). '
            'Check that MoleculeSpec order/SMILES/seed match the builder.'
        )
    return topology, unique


def _assign_charges(offmol, rdmol, cfg) -> None:
    """Assign partial charges to an OpenFF Molecule (nagl, else Gasteiger)."""
    if cfg.charge_method == 'nagl':
        try:
            from openff.nagl_models import load_nagl_model_specs  # noqa: F401
            from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper

            offmol.assign_partial_charges(
                cfg.nagl_model, toolkit_registry=NAGLToolkitWrapper(),
            )
            return
        except Exception as exc:  # noqa: BLE001 - fall back, don't fail prep
            logger.warning(
                'NAGL charge assignment failed (%s); falling back to Gasteiger.',
                exc,
            )

    from openff.units import unit
    from rdkit.Chem import AllChem

    AllChem.ComputeGasteigerCharges(rdmol)
    charges = np.array(
        [float(a.GetDoubleProp('_GasteigerCharge')) for a in rdmol.GetAtoms()],
        dtype=np.float64,
    )
    charges = np.nan_to_num(charges, nan=0.0, posinf=0.0, neginf=0.0)
    offmol.partial_charges = charges * unit.elementary_charge


def _make_system(topology, unique_offmols, cfg):
    """Create an OpenMM System from the topology with the Sage force field."""
    from openff.toolkit import ForceField

    ff = ForceField(cfg.forcefield)
    interchange = ff.create_interchange(
        topology, charge_from_molecules=unique_offmols,
    )
    return interchange.to_openmm(combine_nonbonded_forces=True)


def equilibrate_structure(
    positions_A: NDArray[np.floating],
    species: list[str],
    cell_A: NDArray[np.floating],
    molecule_specs: list[MoleculeSpec],
    cfg: ClassicalPrepConfig | None = None,
) -> PreparedStructure:
    """Relax + densify + thermalize a dilute placement with a classical FF.

    Parameters
    ----------
    positions_A : (N, 3) Å dilute placement from the system builder.
    species     : element symbols, length N, builder order.
    cell_A      : (3, 3) Å cubic cell of the dilute placement.
    molecule_specs : distinct molecules in placement order (init, then monomer).
    cfg         : prep parameters; defaults if None.

    Returns a PreparedStructure (positions + cell in Å) at the target density,
    atom order preserved. Implementation lands in Phase 2 (P2a/P2b).
    """
    cfg = cfg or ClassicalPrepConfig()
    if cfg.protocol != 'simple':
        raise NotImplementedError(f'protocol {cfg.protocol!r} not implemented')

    import openmm
    from openmm import unit as ommunit

    cell_A = np.asarray(cell_A, dtype=np.float64)
    start_edge_A = float(cell_A[0, 0])
    target_edge_A = _target_edge_A(species, molecule_specs, cfg)

    topology, unique = _build_openff_topology(species, molecule_specs, cfg)

    # The box must be set on the topology BEFORE create_interchange so the System
    # is built periodic (PME + cutoffs); otherwise to_openmm() yields a
    # non-periodic (NoCutoff) System and densification is meaningless.
    from openff.units import unit as offunit

    topology.box_vectors = (
        np.eye(3) * start_edge_A * NM_PER_ANGSTROM
    ) * offunit.nanometer

    system = _make_system(topology, unique, cfg)

    integrator = openmm.LangevinMiddleIntegrator(
        cfg.temperature_K * ommunit.kelvin,
        cfg.friction_per_ps / ommunit.picosecond,
        cfg.timestep_fs * ommunit.femtosecond,
    )
    integrator.setRandomNumberSeed(cfg.seed)
    platform = openmm.Platform.getPlatformByName(cfg.platform)
    context = openmm.Context(system, integrator, platform)

    pos_nm = positions_A * NM_PER_ANGSTROM
    _set_cubic_box(context, system, start_edge_A * NM_PER_ANGSTROM)
    context.setPositions(pos_nm * ommunit.nanometer)
    context.setVelocitiesToTemperature(
        cfg.temperature_K * ommunit.kelvin, cfg.seed,
    )

    # 1) minimize close contacts of the dilute placement
    openmm.LocalEnergyMinimizer.minimize(
        context, cfg.minimize_tolerance_kj_mol_nm, 0,
    )
    logger.info('Classical prep: minimized dilute placement.')

    # 2) deterministic compression start_edge -> target_edge
    _compress(context, system, start_edge_A, target_edge_A, cfg)

    # 3) final full minimization at the target box, then fresh velocities, so the
    #    NVT below starts from a low-force state (otherwise residual close
    #    contacts at liquid density blow up the integrator: "coordinate is NaN").
    openmm.LocalEnergyMinimizer.minimize(
        context, cfg.minimize_tolerance_kj_mol_nm, 0,
    )
    context.setVelocitiesToTemperature(cfg.temperature_K * ommunit.kelvin, cfg.seed)
    logger.info('Classical prep: final minimization at target density done.')

    # 4) NVT thermalization at the (fixed) target box. Run in chunks and bail out
    #    early with a clear error if the integrator produces a non-finite state.
    chunk = max(1, cfg.nvt_steps // 20)
    done = 0
    while done < cfg.nvt_steps:
        integrator.step(min(chunk, cfg.nvt_steps - done))
        done += chunk
        e = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            ommunit.kilojoule_per_mole
        )
        if not np.isfinite(e):
            raise RuntimeError(
                f'Classical prep: NVT became non-finite after {done} steps '
                f'(timestep {cfg.timestep_fs} fs may be too large for the dense box).'
            )
    logger.info('Classical prep: NVT thermalization %d steps done.', cfg.nvt_steps)

    state = context.getState(getPositions=True)
    out_pos_nm = np.array(
        state.getPositions().value_in_unit(ommunit.nanometer), dtype=np.float64,
    )
    out_positions_A = out_pos_nm * ANGSTROM_PER_NM
    out_cell_A = np.diag([target_edge_A, target_edge_A, target_edge_A])

    meta = {
        'prep': 'openmm-openff',
        'forcefield': cfg.forcefield,
        'charge_method': cfg.charge_method,
        'target_density_g_per_ml': cfg.target_density_g_per_ml,
        'temperature_K': cfg.temperature_K,
        'start_edge_A': start_edge_A,
        'target_edge_A': target_edge_A,
        'nvt_steps': cfg.nvt_steps,
        'compress_stages': cfg.compress_stages,
        'seed': cfg.seed,
        **cfg.metadata,
    }
    return PreparedStructure(
        positions=out_positions_A, species=list(species),
        cell=out_cell_A, metadata=meta,
    )


def _target_edge_A(species, specs, cfg) -> float:
    """Cubic edge (Å) for the target density given the molecule counts."""
    from kagome.chem.builders import box_from_density

    counts = {spec.smiles: spec.count for spec in specs}
    return box_from_density(counts, cfg.target_density_g_per_ml)


def _set_cubic_box(context, system, edge_nm: float) -> None:
    """Set an isotropic cubic periodic box on both the system and context."""
    import openmm
    from openmm import unit as ommunit

    vectors = (
        openmm.Vec3(edge_nm, 0.0, 0.0) * ommunit.nanometer,
        openmm.Vec3(0.0, edge_nm, 0.0) * ommunit.nanometer,
        openmm.Vec3(0.0, 0.0, edge_nm) * ommunit.nanometer,
    )
    system.setDefaultPeriodicBoxVectors(*vectors)
    context.setPeriodicBoxVectors(*vectors)


def _compress(context, system, start_edge_A, target_edge_A, cfg) -> None:
    """Shrink the cubic box to the target edge in geometric stages.

    Each stage scales box + positions affinely by a constant linear ratio, then
    runs a short MD relaxation of the close contacts the shrink introduces.
    No-op if the box is already at/below the target edge.
    """
    from openmm import unit as ommunit

    if target_edge_A >= start_edge_A:
        logger.info(
            'Classical prep: start edge %.2f Å <= target %.2f Å — no compression.',
            start_edge_A, target_edge_A,
        )
        return

    ratio = (target_edge_A / start_edge_A) ** (1.0 / cfg.compress_stages)
    cur_edge_A = start_edge_A
    logger.info(
        'Classical prep: compressing %.2f Å -> %.2f Å in %d stages.',
        start_edge_A, target_edge_A, cfg.compress_stages,
    )
    import openmm

    for stage in range(1, cfg.compress_stages + 1):
        state = context.getState(getPositions=True)
        pos_nm = np.array(
            state.getPositions().value_in_unit(ommunit.nanometer),
            dtype=np.float64,
        )
        pos_nm *= ratio
        cur_edge_A *= ratio
        _set_cubic_box(context, system, cur_edge_A * NM_PER_ANGSTROM)
        context.setPositions(pos_nm * ommunit.nanometer)
        # Affine scaling introduces close contacts; energy-minimize (stable, no
        # timestep) to remove them BEFORE any MD, otherwise the huge overlap
        # forces blow up the integrator ("Particle coordinate is NaN").
        openmm.LocalEnergyMinimizer.minimize(
            context, cfg.minimize_tolerance_kj_mol_nm, cfg.compress_minimize_iters,
        )
        if cfg.compress_relax_steps > 0:
            context.getIntegrator().step(cfg.compress_relax_steps)
    logger.info('Classical prep: compression complete (edge %.2f Å).', cur_edge_A)
