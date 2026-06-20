"""Classical force-field backend (OpenMM + OpenFF Sage).

Exposes a classical FF through pfpoly's ``Calculator`` interface so it can drive
``src.integrators.minimize.compress_box`` (and any other Calculator consumer)
in-process. The intended use is the cheap structure-preparation work —
densifying a dilute placement to the paper's 0.5 g/mL initial density — without
spending OrbMol-v2 GPU time on it. The MD production still runs on the MLIP.

Why topology-aware: unlike an MLIP (which infers bonding from geometry), a
classical FF needs the molecular graph to assign parameters. So this calculator
is constructed from ``MoleculeSpec`` list (the same specs the system builder and
``src/prep/openmm_equilibrate`` use), builds the OpenMM System once, then
``compute`` only updates coordinates and the periodic box per call.

Atom-order invariant: the OpenFF Topology is built in builder order via
``src.prep.openmm_equilibrate._build_openff_topology``, which validates the
element order against the caller's ``species``. So returned forces map 1:1 onto
the caller's coordinates.

Decision: specs/decisions.md 2026-06-20 "WSL 単一env で古典FFを Calculator 化し
compress の既定バックエンドにする". Licensing: OpenMM (MIT core / LGPL GPU,
import-only), OpenFF stack (MIT), Sage (CC-BY-4.0, attribution), RDKit Gasteiger
(BSD) — all approved (specs/approved_dependencies.yaml). No new dependency.

OpenMM imports are deferred to call time so this module imports cleanly in
environments that only run ML production.
"""
from __future__ import annotations

import logging

import numpy as np
from numpy.typing import NDArray

from src.backends.base import Calculator
from src.prep.openmm_equilibrate import ClassicalPrepConfig, MoleculeSpec
from src.units import KCAL_PER_KJ, NM_PER_ANGSTROM

logger = logging.getLogger(__name__)

# kJ/mol/nm -> kcal/mol/Å : (1/4.184) kcal per kJ, and per-nm -> per-Å is ×0.1.
_FORCE_KJNM_TO_KCALA = KCAL_PER_KJ * NM_PER_ANGSTROM


def _safe_cutoff_nm(cutoff_nm: float | None, min_box_edge_A: float | None) -> float | None:
    """Clamp the nonbonded cutoff so it stays below half the smallest box edge.

    OpenMM rejects a periodic box edge < 2×cutoff. The densest box a compression
    reaches is ``min_box_edge_A`` (the target edge); the cutoff is clamped to a
    1 Å margin below half of it, floored at 0.4 nm so it stays physically sane.
    ``cutoff_nm=None`` means 'use the force field default' (no clamp)."""
    if cutoff_nm is None or min_box_edge_A is None:
        return cutoff_nm
    safe = (min_box_edge_A / 2.0 - 1.0) * NM_PER_ANGSTROM
    return max(0.4, min(cutoff_nm, safe))


def create_classical_calculator(
    molecule_specs: list[MoleculeSpec],
    charge_method: str = 'gasteiger',
    forcefield: str = 'openff-2.2.0.offxml',
    platform: str = 'CPU',
    cutoff_nm: float | None = 0.8,
    min_box_edge_A: float | None = None,
) -> 'ClassicalCalculator':
    """Create a classical (OpenMM/OpenFF Sage) Calculator.

    molecule_specs : distinct molecules in builder placement order (e.g. for the
        vinyl system: initiators first, then monomers). Their SMILES/seed must
        match what the system builder passed so the topology aligns with coords.
    charge_method  : 'gasteiger' (RDKit, default — lightweight, no model
        download; the prep structure is overwritten by the ML re-equilibration
        anyway, decision D-2) or 'nagl'.
    forcefield     : OpenFF offxml (Sage 2.2 by default).
    platform       : OpenMM platform name ('CPU' | 'CUDA' | 'OpenCL' | 'Reference').
    cutoff_nm      : nonbonded cutoff ceiling (nm). Defaults to 0.8 nm so the
        cutoff stays below half the dense paper box edge (~18-21 Å), otherwise
        OpenMM rejects a periodic box smaller than twice the cutoff during
        compression. Pass None to keep the force field's own default (0.9 nm).
    min_box_edge_A : the smallest box edge compression will reach (the target
        edge). When given, the cutoff is clamped to a 1 Å margin below half of
        it, so small (sub-16 Å) boxes are accepted too.

    The OpenMM System is built lazily on the first ``compute`` call (when the
    ``species`` order can be validated against the topology).
    """
    cfg = ClassicalPrepConfig(
        charge_method=charge_method,
        forcefield=forcefield,
        platform=platform,
    )
    return ClassicalCalculator(
        molecule_specs, cfg, platform=platform,
        cutoff_nm=_safe_cutoff_nm(cutoff_nm, min_box_edge_A),
    )


def make_compress_calculator(
    compress_backend: str,
    molecule_specs: list[MoleculeSpec],
    ml_calculator: Calculator,
    platform: str = 'CPU',
    target_edge_A: float | None = None,
) -> Calculator:
    """Return the Calculator used for box compression to paper density.

    'classical' (the run-script default) builds an OpenMM/OpenFF Sage Calculator
    from ``molecule_specs`` so densification does not consume MLIP GPU time
    (decision 2026-06-20). 'ml' reuses the production MLIP calculator. The MD
    production itself always runs on the MLIP regardless of this choice.

    target_edge_A : the densest (smallest) box the compression reaches, so the
        classical cutoff can be clamped below half of it (OpenMM box constraint).
    """
    if compress_backend == 'classical':
        return create_classical_calculator(
            molecule_specs, platform=platform, min_box_edge_A=target_edge_A,
        )
    if compress_backend == 'ml':
        return ml_calculator
    raise ValueError(f'unknown compress_backend {compress_backend!r}')


class ClassicalCalculator(Calculator):
    """OpenMM+OpenFF Sage energies/forces behind the pfpoly Calculator API."""

    def __init__(
        self,
        molecule_specs: list[MoleculeSpec],
        cfg: ClassicalPrepConfig,
        platform: str = 'CPU',
        cutoff_nm: float | None = 0.8,
    ) -> None:
        self._specs = list(molecule_specs)
        self._cfg = cfg
        self._platform_name = platform
        self._cutoff_nm = cutoff_nm
        self._name = f'classical-{cfg.forcefield}'
        self._context = None       # built lazily on first compute
        self._system = None
        self._n_atoms = sum(s.count for s in molecule_specs) if molecule_specs else 0

    @property
    def name(self) -> str:
        return self._name

    @property
    def model_id(self) -> str:
        return f'{self._name}/{self._cfg.charge_method}'

    def _build(self, species: list[str], cell: NDArray[np.floating] | None) -> None:
        """Build the OpenMM System+Context once (validates atom order)."""
        import openmm
        from openff.units import unit as offunit

        from src.prep.openmm_equilibrate import _build_openff_topology, _make_system

        topology, unique = _build_openff_topology(species, self._specs, self._cfg)

        periodic = cell is not None
        if periodic:
            # The box must be set BEFORE create_interchange so the System is built
            # periodic (PME); otherwise to_openmm() yields a NoCutoff System and
            # densification is meaningless (mirrors openmm_equilibrate).
            edge_nm = float(np.asarray(cell, dtype=np.float64)[0, 0]) * NM_PER_ANGSTROM
            topology.box_vectors = (np.eye(3) * edge_nm) * offunit.nanometer

        system = _make_system(topology, unique, self._cfg)

        if periodic and self._cutoff_nm is not None:
            self._apply_cutoff(system, self._cutoff_nm)

        integrator = openmm.VerletIntegrator(1.0 * openmm.unit.femtosecond)
        platform = openmm.Platform.getPlatformByName(self._platform_name)
        self._system = system
        self._context = openmm.Context(system, integrator, platform)
        logger.info(
            'ClassicalCalculator built: %d atoms, FF=%s, charges=%s, platform=%s, '
            'periodic=%s, cutoff_nm=%s.',
            len(species), self._cfg.forcefield, self._cfg.charge_method,
            self._platform_name, periodic, self._cutoff_nm,
        )

    @staticmethod
    def _apply_cutoff(system, cutoff_nm: float) -> None:
        """Clamp every nonbonded cutoff so a dense periodic box is accepted.

        OpenMM rejects a periodic box edge smaller than twice the cutoff. The
        paper boxes (~18-21 Å) are barely above twice the FF default (0.9 nm), so
        we shrink the cutoff to keep compression stable. This is a preparation
        device only; the structure is re-equilibrated by the MLIP afterward.

        Sage applies a switching function (default 0.8 nm switch under a 0.9 nm
        cutoff). Shrinking the cutoff to the switch distance would violate
        ``switch < cutoff``, so the switch distance is lowered in step, keeping a
        0.1 nm switching window.
        """
        import openmm

        nm = openmm.unit.nanometer
        switch_nm = max(0.0, cutoff_nm - 0.1)
        for force in system.getForces():
            if isinstance(force, openmm.NonbondedForce):
                if force.getUseSwitchingFunction():
                    force.setSwitchingDistance(switch_nm * nm)
                force.setCutoffDistance(cutoff_nm * nm)

    def compute(
        self,
        positions: NDArray[np.floating],
        species: list[str],
        cell: NDArray[np.floating] | None = None,
    ) -> tuple[float, NDArray[np.floating]]:
        import openmm
        from openmm import unit as ommunit

        if self._context is None:
            self._build(species, cell)

        if cell is not None:
            cell = np.asarray(cell, dtype=np.float64)
            vectors = (
                openmm.Vec3(*(cell[0] * NM_PER_ANGSTROM)),
                openmm.Vec3(*(cell[1] * NM_PER_ANGSTROM)),
                openmm.Vec3(*(cell[2] * NM_PER_ANGSTROM)),
            )
            self._context.setPeriodicBoxVectors(*[v * ommunit.nanometer for v in vectors])

        pos_nm = np.asarray(positions, dtype=np.float64) * NM_PER_ANGSTROM
        self._context.setPositions(pos_nm * ommunit.nanometer)

        state = self._context.getState(getEnergy=True, getForces=True)
        energy_kj = state.getPotentialEnergy().value_in_unit(ommunit.kilojoule_per_mole)
        forces_kjnm = np.asarray(
            state.getForces(asNumpy=True).value_in_unit(
                ommunit.kilojoule_per_mole / ommunit.nanometer
            ),
            dtype=np.float64,
        )

        energy = float(energy_kj) * KCAL_PER_KJ
        forces = forces_kjnm * _FORCE_KJNM_TO_KCALA
        return energy, forces
