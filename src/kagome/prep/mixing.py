"""Bond-graph -> classical (OpenFF/OpenMM) translator for well-mixed mixing.

Part of the "well-mixed 測定モード" (specs/decisions.md 2026-07-17). Between TDBB
cycles the reactivity-measurement protocol runs a short classical MD *mixing*
stage so encounter statistics stop being dominated by the frozen initial packing.
This module is the translator (Phase P2): it turns the *evolving* reactive bond
graph (:class:`kagome.reactive.topology.BondTopology`) plus the element list and
coordinates into an OpenMM ``System`` parameterized with OpenFF Sage 2.2 + NAGL
charges (Gasteiger fallback), and hands back a mapping so the mixed coordinates
can be written straight back onto the MLIP atom order.

It shares the charge-assignment implementation with the prep stage
(:mod:`kagome.prep.charges` — NAGL first, Gasteiger fallback) and mirrors the
rest of :mod:`kagome.prep.openmm_equilibrate`; the difference is that the connectivity
comes from a live ``BondTopology`` (molecules span non-contiguous global indices
once chains grow) rather than from static per-SMILES ``MoleculeSpec`` blocks, so
this module keeps an explicit ``omm <-> mlip`` index map instead of relying on a
contiguous layout.

Radical / placeholder handling (specs/decisions.md 2026-07-17, decision (iii)):

* A growth-end radical centre (a heavy atom whose summed bond order is below its
  neutral valence — typically a 3-coordinate carbon after the vinyl opened) is
  open-shell, which OpenFF cannot parameterize. We inject a *cap H* per missing
  valence unit that exists only in the classical world and is dropped on
  write-back.
* A placeholder H shed from an initiator radical (``apply_vinyl_addition``) is
  left disconnected from the graph (coordination 0). OpenFF rejects a bare
  ``[H]`` (RadicalsNotSupportedError), so such atoms are appended to the OpenMM
  ``System`` as nonbonded-only, charge-neutral particles carrying a Sage-derived
  aliphatic-H Lennard-Jones size (no bonded terms), then mapped back unchanged.

Scope (P2): translator + tests only. Workflow integration (the mixing phase,
``MixConfig``, velocity re-draw, ``mix_settle``) is Phase P3 and lives elsewhere;
the public API here (:func:`build_classical_mix`, :class:`ClassicalMix`,
:class:`FragmentParamCache`, :class:`MixTranslatorConfig`) is designed to be the
clean seam P3 builds on. No TDBB code is touched.

All OpenMM / OpenFF / RDKit imports are deferred to call time so this module
imports cleanly in ML-only environments.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from kagome.analysis.carothers import monomer_sets_from_bonds
from kagome.prep.charges import CHARGE_METHODS, assign_charges
from kagome.reactive.topology import BondTopology
from kagome.units import NM_PER_ANGSTROM

logger = logging.getLogger(__name__)

# Neutral (uncharged, closed-shell) valence per element: the summed bond order a
# saturated atom carries. A heavy atom whose summed order falls below this is an
# open-shell radical centre in this reactive world (the only source of a deficit
# here is an opened vinyl / shed placeholder H — see module docstring and
# specs/decisions.md 2026-07-17). Mirrors topology.MAX_COORDINATION but keyed on
# *bond order* (so a C=O carbonyl carbon, order sum 4, is correctly NOT a radical
# even though its coordination number is 3).
NEUTRAL_VALENCE: dict[str, int] = {
    'H': 1, 'C': 4, 'N': 3, 'O': 2, 'F': 1, 'S': 2, 'Cl': 1,
}

# Cap-H bond length placed away from the radical centre's existing neighbours (Å).
# A generic C-H single-bond length; only the starting geometry, relaxed by the
# classical minimizer P3 runs before mixing.
_CAP_H_BOND_A = 1.09

# Tilt (degrees) of each cap-H direction off the vacant-valence axis when a
# centre needs MORE than one cap: half the tetrahedral angle, so 2-4 caps spread
# azimuthally like the vacant sp3 lobes instead of stacking on one point (which
# would give a zero H-C-H angle and NaN angle forces). Deterministic geometry,
# no rng. See specs/decisions.md 2026-07-17 addendum.
_MULTI_CAP_TILT_DEG = 54.7356

# Cached Sage-derived aliphatic-H (sigma_nm, epsilon_kJ/mol), keyed by forcefield
# name, used for placeholder-H nonbonded params. Derived from the force field
# itself (methane's H) rather than hardcoded, so it tracks the chosen FF.
_REF_H_VDW: dict[str, tuple[float, float]] = {}


@dataclass
class MixTranslatorConfig:
    """Force-field / charge options for the bond-graph -> OpenMM translation.

    Defaults mirror the classical prep stage (Sage 2.2 + NAGL, Gasteiger
    fallback) so the mixing PES matches the prep PES. ``charge_method`` is
    ``'nagl'`` | ``'gasteiger'``; NAGL falls back to Gasteiger automatically on
    any failure (e.g. no network for the model download).

    ``combine_nonbonded_forces`` must stay ``True`` (validated): placeholder H
    are appended to the single combined ``NonbondedForce`` only, and a split
    vdW/electrostatics layout would leave them missing from the vdW
    ``CustomNonbondedForce`` (particle-count mismatch at Context creation).
    The field is kept to document the constraint and reserve the option.
    """

    charge_method: str = 'nagl'                     # 'nagl' | 'gasteiger'
    forcefield: str = 'openff-2.2.0.offxml'         # Sage 2.2
    nagl_model: str = 'openff-gnn-am1bcc-0.1.0-rc.3.pt'
    combine_nonbonded_forces: bool = True           # must remain True (see above)

    def __post_init__(self) -> None:
        if self.charge_method not in CHARGE_METHODS:
            raise ValueError(
                f'charge_method must be one of {CHARGE_METHODS}; '
                f'got {self.charge_method!r}'
            )
        if not self.combine_nonbonded_forces:
            raise ValueError(
                'combine_nonbonded_forces=False is not supported: placeholder '
                'H particles are added to the combined NonbondedForce only, so '
                'a split vdW force would have a particle-count mismatch '
                '(specs/decisions.md 2026-07-17 addendum).'
            )


@dataclass
class FragmentParamCache:
    """Caches OpenFF-parameterized (charge-assigned) fragment templates.

    Keyed by :func:`fragment_cache_key` — the fragment's canonical graph (RDKit
    canonical SMILES of the *capped* molecule) prefixed with the charge/FF
    options (``charge_method``, ``nagl_model``, ``forcefield``), so a cache
    shared across differently-configured builds can never serve a template
    charged under another method. Free monomers and equal-length chains are
    graph-isomorphic, so the expensive per-fragment charge assignment (NAGL
    model inference) runs once per distinct graph instead of once per molecule
    per cycle. Cache values are ``(template, method_used)`` pairs — the charged
    OpenFF ``Molecule`` plus the charge method that actually produced it
    (``'gasteiger'`` when NAGL fell back); interchange maps each stored template
    onto every isomorphic copy in the topology by graph isomorphism, so
    per-instance atom order does not matter.

    ``hits`` / ``misses`` are exposed for test assertions and audit logging.
    """

    _templates: dict = field(default_factory=dict)
    hits: int = 0
    misses: int = 0

    def get_or_assign(self, key: str, build_template):
        """Return the cached ``(template, method_used)`` for ``key`` or build it.

        ``build_template`` is a zero-arg callable returning a
        ``(charged OpenFF Molecule, method_used)`` pair; it is only invoked on a
        miss.
        """
        if key in self._templates:
            self.hits += 1
            return self._templates[key]
        self.misses += 1
        value = build_template()
        self._templates[key] = value
        return value


def fragment_cache_key(canonical_smiles: str, cfg: MixTranslatorConfig) -> str:
    """Cache key for a capped fragment graph under the given charge/FF options.

    Includes the requested charge method, NAGL model, and force field so
    templates charged under one configuration are never reused by another
    (review 2026-07-17: charge-template contamination).
    """
    return (
        f'{cfg.charge_method}|{cfg.nagl_model}|{cfg.forcefield}|'
        f'{canonical_smiles}'
    )


@dataclass
class ClassicalMix:
    """A classical OpenMM system built from a reactive bond graph.

    .. warning::
        ``positions_A`` are a *starting guess*, not an equilibrated state: cap H
        sit at idealized directions and a shed placeholder H still sits ~1 Å
        from its former parent (severe LJ overlap; energy is finite but huge).
        The consumer (P3) MUST energy-minimize before any dynamics —
        ``metadata['requires_minimization']`` records this contract.

    Attributes
    ----------
    system : openmm.System
        Periodic Sage-parameterized system. Atom order is the OpenMM order
        (molecule-by-molecule, cap H appended within each molecule, placeholder
        H appended last), NOT the MLIP order.
    positions_A : (M, 3) float ndarray
        Starting coordinates in Å, in OpenMM atom order (``M`` = number of
        OpenMM particles = real atoms + cap H + placeholder H).
    box_vectors_A : (3, 3) float ndarray
        Periodic cell in Å.
    omm_to_mlip : list[int | None]
        Length ``M``; ``omm_to_mlip[k]`` is the MLIP global index of OpenMM
        particle ``k``, or ``None`` if it is a classical-only cap H to discard.
    mlip_to_omm : dict[int, int]
        Inverse map over every MLIP atom (real + placeholder H). Complete: every
        index ``0..N-1`` appears exactly once.
    n_mlip_atoms : int
        ``N`` — the MLIP atom count (write-back output length).
    metadata : dict
        Audit fields (counts of cap/placeholder H, charge method, etc.).
    """

    system: object
    positions_A: NDArray[np.floating]
    box_vectors_A: NDArray[np.floating]
    omm_to_mlip: list[int | None]
    mlip_to_omm: dict[int, int]
    n_mlip_atoms: int
    metadata: dict = field(default_factory=dict)

    def write_back(
        self, omm_positions_A: NDArray[np.floating],
    ) -> NDArray[np.floating]:
        """Map mixed OpenMM coordinates (Å, OpenMM order) back to MLIP order.

        Cap H are dropped; real atoms and placeholder H are restored to their
        original global indices. Returns an ``(N, 3)`` Å array aligned with the
        caller's ``species`` / ``groups`` / ``propagation_map``.
        """
        omm_positions_A = np.asarray(omm_positions_A, dtype=np.float64)
        if omm_positions_A.shape != (len(self.omm_to_mlip), 3):
            raise ValueError(
                f'omm_positions must be ({len(self.omm_to_mlip)}, 3); '
                f'got {omm_positions_A.shape}'
            )
        if not np.isfinite(omm_positions_A).all():
            bad = sorted(
                np.where(~np.isfinite(omm_positions_A).all(axis=1))[0].tolist()
            )
            raise ValueError(
                'write_back: input coordinates are non-finite at OpenMM '
                f'particle indices {bad[:10]}'
                f'{"..." if len(bad) > 10 else ""} — the classical MD likely '
                'diverged (NaN); fix the mixing run, not the mapping.'
            )
        out = np.full((self.n_mlip_atoms, 3), np.nan, dtype=np.float64)
        for omm_idx, mlip_idx in enumerate(self.omm_to_mlip):
            if mlip_idx is not None:
                out[mlip_idx] = omm_positions_A[omm_idx]
        if np.isnan(out).any():
            missing = sorted(np.where(np.isnan(out).any(axis=1))[0].tolist())
            raise RuntimeError(
                f'write_back: no OpenMM coordinate for MLIP atoms {missing}'
            )
        return out


# ── graph decomposition (pure Python, MD-free / no heavy imports) ─────────────

def connected_components(
    topology: BondTopology, n_atoms: int,
) -> list[list[int]]:
    """Partition atoms ``0..n_atoms-1`` into bonded components (sorted).

    Each component is returned as a sorted list of global indices; the list of
    components is sorted by each component's minimum index for determinism.
    Isolated atoms (no bonds) come back as singletons.

    Reuses the carothers component finder instead of a second union-find (the
    ``network.largest_component_fraction`` precedent): one self-edge per atom
    keeps isolated atoms (shed placeholder H) in the component list.
    """
    edges: list[tuple[int, int]] = [(a, a) for a in range(n_atoms)]
    edges.extend(
        (i, j) for i, j, _order in topology.bonds()
        if i < n_atoms and j < n_atoms
    )
    return monomer_sets_from_bonds(edges)


def radical_cap_counts(
    topology: BondTopology, species: list[str], atoms: list[int],
) -> dict[int, int]:
    """Number of cap H each heavy atom in ``atoms`` needs to close its shell.

    ``deficit = round(neutral_valence[element] - summed_bond_order)``; a positive
    deficit marks an open-shell radical centre (an opened vinyl growth end). H
    atoms are never capped (that would build H2) — a lone/deficient H is a
    placeholder handled separately.

    An element OUTSIDE :data:`NEUTRAL_VALENCE` that is itself valence-deficient
    (per RDKit's default valence) raises ``ValueError`` instead of being
    silently skipped: an uncapped open-shell atom would otherwise reach OpenFF
    (or grow implicit H) and corrupt the ``omm <-> mlip`` index map (review
    2026-07-17). Saturated unknown elements pass through uncapped.
    """
    counts: dict[int, int] = {}
    for atom in atoms:
        el = species[atom]
        if el == 'H':
            continue
        order_sum = sum(topology.order(atom, nbr)
                        for nbr in topology.neighbors(atom))
        if el not in NEUTRAL_VALENCE:
            _check_unknown_element_saturated(el, atom, order_sum)
            continue
        deficit = round(NEUTRAL_VALENCE[el] - order_sum)
        if deficit > 0:
            counts[atom] = deficit
    return counts


def _check_unknown_element_saturated(el: str, atom: int, order_sum: float) -> None:
    """Raise if an element outside NEUTRAL_VALENCE is valence-deficient.

    Uses RDKit's periodic-table default valence as the reference. Elements whose
    default valence RDKit does not define (returns < 0, e.g. transition metals)
    are rejected outright — the translator cannot decide whether they need caps.
    """
    from rdkit import Chem

    default = Chem.GetPeriodicTable().GetDefaultValence(el)
    if default < 0:
        raise ValueError(
            f'element {el!r} (atom {atom}) is outside NEUTRAL_VALENCE and has '
            'no RDKit default valence; the translator cannot determine its cap '
            'count. Extend NEUTRAL_VALENCE to support it.'
        )
    if round(default - order_sum) > 0:
        raise ValueError(
            f'atom {atom} (element {el!r}) is valence-deficient (bond-order sum '
            f'{order_sum} < default valence {default}) but {el!r} is outside '
            'NEUTRAL_VALENCE — refusing to translate an uncapped open-shell '
            'atom. Extend NEUTRAL_VALENCE to support it.'
        )


def is_placeholder_h(
    topology: BondTopology, species: list[str], atom: int,
) -> bool:
    """True if ``atom`` is a disconnected placeholder H (H with no bonds)."""
    return species[atom] == 'H' and topology.coordination_number(atom) == 0


# ── RDKit fragment construction ──────────────────────────────────────────────

_RDKIT_BOND_TYPE = {1.0: 'SINGLE', 1.5: 'AROMATIC', 2.0: 'DOUBLE', 3.0: 'TRIPLE'}


def _build_fragment_rdkit(
    topology: BondTopology,
    species: list[str],
    atoms: list[int],
    cap_counts: dict[int, int],
):
    """Build a sanitized RDKit Mol for one component, with cap H appended.

    Returns ``(mol, local_to_global, cap_parents)`` where ``local_to_global[k]``
    is the MLIP index of local atom ``k`` for the ``len(atoms)`` real atoms, or
    ``None`` for each appended cap H, and ``cap_parents`` lists the global index
    of the radical parent for each cap H in append order. Real atoms are added in
    ascending global-index order so the OpenMM order is deterministic.
    """
    from rdkit import Chem

    rw = Chem.RWMol()
    global_to_local: dict[int, int] = {}
    local_to_global: list[int | None] = []
    for g in atoms:
        local_to_global.append(g)
        atom = Chem.Atom(species[g])
        # The bond graph is the complete connectivity: never let RDKit grow
        # implicit hydrogens (they would add OpenMM particles that are absent
        # from local_to_global and shift the omm<->mlip map).
        atom.SetNoImplicit(True)
        global_to_local[g] = rw.AddAtom(atom)

    bt = Chem.BondType
    type_map = {
        'SINGLE': bt.SINGLE, 'AROMATIC': bt.AROMATIC,
        'DOUBLE': bt.DOUBLE, 'TRIPLE': bt.TRIPLE,
    }
    atom_set = set(atoms)
    for i, j, order in topology.bonds():
        if i in atom_set and j in atom_set:
            name = _RDKIT_BOND_TYPE.get(float(order))
            if name is None:
                raise ValueError(f'unsupported bond order {order} in fragment')
            rw.AddBond(global_to_local[i], global_to_local[j], type_map[name])
            if name == 'AROMATIC':
                rw.GetAtomWithIdx(global_to_local[i]).SetIsAromatic(True)
                rw.GetAtomWithIdx(global_to_local[j]).SetIsAromatic(True)

    # Cap each radical centre with H (classical-only; discarded on write-back).
    cap_parents: list[int] = []
    for g in sorted(cap_counts):
        if g not in atom_set:
            continue
        for _ in range(cap_counts[g]):
            cap_atom = Chem.Atom('H')
            cap_atom.SetNoImplicit(True)
            cap_local = rw.AddAtom(cap_atom)
            rw.AddBond(global_to_local[g], cap_local, bt.SINGLE)
            local_to_global.append(None)
            cap_parents.append(g)

    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    return mol, local_to_global, cap_parents


def _reference_h_vdw(cfg: MixTranslatorConfig) -> tuple[float, float]:
    """Aliphatic-H (sigma_nm, epsilon_kJ/mol) from the FF itself (methane H)."""
    if cfg.forcefield in _REF_H_VDW:
        return _REF_H_VDW[cfg.forcefield]
    import openmm
    from openmm import unit as ommunit
    from openff.toolkit import ForceField, Molecule, Topology

    methane = Molecule.from_smiles('C')
    methane.assign_partial_charges('gasteiger')
    top = Topology.from_molecules([methane])
    ff = ForceField(cfg.forcefield)
    inter = ff.create_interchange(top, charge_from_molecules=[methane])
    system = inter.to_openmm(combine_nonbonded_forces=True)
    nb = next(f for f in system.getForces()
              if isinstance(f, openmm.NonbondedForce))
    _q, sigma, eps = nb.getParticleParameters(1)  # atom 1 is an H of methane
    vdw = (sigma.value_in_unit(ommunit.nanometer),
           eps.value_in_unit(ommunit.kilojoule_per_mole))
    _REF_H_VDW[cfg.forcefield] = vdw
    return vdw


def _clamp_nonbonded_cutoff(system, box_nm) -> None:
    """Shrink the nonbonded cutoff so it never exceeds half the periodic box.

    Sage's default cutoff is 0.9 nm; a reactive run's box can be smaller than
    1.8 nm (dense / small systems), and OpenMM refuses to create a Context when
    the cutoff exceeds half the box edge ('cutoff distance cannot be greater
    than half the periodic box size'). The mixing translator, unlike the prep
    stage, has no control over the caller's box, so clamp the cutoff (and keep
    the switching distance just inside it) to the box actually handed in. A
    no-op for the usual case where the box is comfortably larger than 2×cutoff.
    """
    import openmm
    from openmm import unit as ommunit

    min_edge_nm = float(np.min(np.diag(np.asarray(box_nm, dtype=np.float64))))
    # Safety margin below the hard half-box limit (OpenMM uses the in-plane
    # box size for triclinic; our boxes are cubic, so 0.49 is a clean guard).
    max_cutoff_nm = 0.49 * min_edge_nm
    for force in system.getForces():
        if not isinstance(force, openmm.NonbondedForce):
            continue
        cutoff_nm = force.getCutoffDistance().value_in_unit(ommunit.nanometer)
        if cutoff_nm <= max_cutoff_nm:
            continue
        force.setCutoffDistance(max_cutoff_nm * ommunit.nanometer)
        if force.getUseSwitchingFunction():
            # Keep the switch a little inside the new cutoff (mirror Sage's
            # ~0.1 nm switch width, but never let it go non-positive).
            switch_nm = max(max_cutoff_nm - 0.1, 0.5 * max_cutoff_nm)
            force.setSwitchingDistance(switch_nm * ommunit.nanometer)
        logger.warning(
            'Mixing box edge %.2f nm is small; clamped nonbonded cutoff '
            '%.3f -> %.3f nm to satisfy the half-box limit.',
            min_edge_nm, cutoff_nm, max_cutoff_nm,
        )


# ── public entry point ───────────────────────────────────────────────────────

def build_classical_mix(
    topology: BondTopology,
    species: list[str],
    positions_A: NDArray[np.floating],
    cell_A: NDArray[np.floating],
    cfg: MixTranslatorConfig | None = None,
    cache: FragmentParamCache | None = None,
) -> ClassicalMix:
    """Translate a reactive bond graph into a classical OpenMM system.

    Parameters
    ----------
    topology : BondTopology
        The current (evolving) reactive bond graph. Bond orders are honoured;
        heavy-atom valence deficits are capped with classical-only H.
    species : list[str]
        Element symbols, length ``N``, MLIP order (indices match ``topology``).
    positions_A : (N, 3) float
        Coordinates in Å, MLIP order.
    cell_A : (3, 3) float
        Periodic cell in Å (required — the mixing MD is periodic/NVT).
    cfg : MixTranslatorConfig, optional
        Force-field / charge options; defaults to Sage 2.2 + NAGL.
    cache : FragmentParamCache, optional
        Reused across cycles to skip re-charging isomorphic fragments. A fresh
        cache is created if ``None`` (no cross-call reuse).

    Returns
    -------
    ClassicalMix
        The system, starting positions (OpenMM order), and the ``omm <-> mlip``
        maps for coordinate write-back.
    """
    cfg = cfg or MixTranslatorConfig()
    cache = cache if cache is not None else FragmentParamCache()

    positions_A = np.asarray(positions_A, dtype=np.float64)
    n_atoms = len(species)
    if positions_A.shape != (n_atoms, 3):
        raise ValueError(
            f'positions must be ({n_atoms}, 3); got {positions_A.shape}'
        )
    cell_A = np.asarray(cell_A, dtype=np.float64)
    if cell_A.shape != (3, 3):
        raise ValueError(f'cell must be (3, 3); got {cell_A.shape}')

    import openmm
    from openmm import unit as ommunit
    from openff.toolkit import Molecule, Topology
    from openff.units import unit as offunit

    components = connected_components(topology, n_atoms)

    instance_offmols: list = []          # one per fragment, Topology order
    unique_templates: list = []          # charged templates, charge_from_molecules
    unique_keys: set[str] = set()        # cache keys already in unique_templates
    omm_to_mlip: list[int | None] = []   # OpenMM particle -> MLIP idx (or None)
    frag_positions: list[NDArray] = []   # per-fragment (n_local, 3) Å
    placeholder_atoms: list[int] = []    # lone H global indices
    charge_methods_used: set[str] = set()
    n_cap_h = 0

    for atoms in components:
        # A lone hydrogen with no bonds is a shed placeholder -> nonbonded ghost.
        if len(atoms) == 1 and is_placeholder_h(topology, species, atoms[0]):
            placeholder_atoms.append(atoms[0])
            continue

        cap_counts = radical_cap_counts(topology, species, atoms)
        mol, local_to_global, cap_parents = _build_fragment_rdkit(
            topology, species, atoms, cap_counts,
        )
        n_cap_h += sum(1 for g in local_to_global if g is None)

        offmol = Molecule.from_rdkit(mol, allow_undefined_stereo=True)
        if offmol.n_atoms != len(local_to_global):
            raise ValueError(
                f'fragment atom count changed in RDKit->OpenFF conversion '
                f'({offmol.n_atoms} != {len(local_to_global)}) for the fragment '
                f'containing atoms {atoms[:6]}... — the omm<->mlip map would '
                'be corrupt (implicit H or dropped atoms?).'
            )
        key = fragment_cache_key(_fragment_key(mol), cfg)

        def _make_template(m=mol):
            template = Molecule.from_rdkit(m, allow_undefined_stereo=True)
            method = assign_charges(template, cfg.charge_method, cfg.nagl_model)
            return template, method

        template, method_used = cache.get_or_assign(key, _make_template)
        charge_methods_used.add(method_used)
        if key not in unique_keys:
            unique_keys.add(key)
            unique_templates.append(template)

        instance_offmols.append(offmol)
        omm_to_mlip.extend(local_to_global)
        frag_positions.append(
            _fragment_positions(
                local_to_global, cap_parents, positions_A, topology, n_atoms,
            )
        )

    # Assemble the OpenFF topology (all fragment instances) and Sage system.
    off_top = Topology.from_molecules(instance_offmols)
    _verify_topology_elements(off_top, omm_to_mlip, species)
    off_top.box_vectors = (
        cell_A * NM_PER_ANGSTROM
    ) * offunit.nanometer

    from openff.toolkit import ForceField

    ff = ForceField(cfg.forcefield)
    interchange = ff.create_interchange(
        off_top, charge_from_molecules=unique_templates,
    )
    system = interchange.to_openmm(
        combine_nonbonded_forces=cfg.combine_nonbonded_forces,
    )

    positions_list = (
        [p for frag in frag_positions for p in frag] if frag_positions else []
    )

    # Append placeholder H as nonbonded-only, charge-neutral particles.
    if placeholder_atoms:
        sigma_nm, eps_kj = _reference_h_vdw(cfg)
        nb = next(f for f in system.getForces()
                  if isinstance(f, openmm.NonbondedForce))
        for g in placeholder_atoms:
            if species[g] != 'H':
                raise ValueError(
                    f'placeholder particle {g} is {species[g]!r}, expected H'
                )
            sys_idx = system.addParticle(1.008 * ommunit.amu)
            nb_idx = nb.addParticle(
                0.0 * ommunit.elementary_charge,
                sigma_nm * ommunit.nanometer,
                eps_kj * ommunit.kilojoule_per_mole,
            )
            assert sys_idx == nb_idx == system.getNumParticles() - 1
            omm_to_mlip.append(g)
            positions_list.append(positions_A[g])

    # Ensure the system box matches (interchange sets it from off_top, but keep
    # System and returned box explicitly consistent for P3's context setup).
    box_nm = cell_A * NM_PER_ANGSTROM
    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(*box_nm[0]) * ommunit.nanometer,
        openmm.Vec3(*box_nm[1]) * ommunit.nanometer,
        openmm.Vec3(*box_nm[2]) * ommunit.nanometer,
    )
    _clamp_nonbonded_cutoff(system, box_nm)

    positions_out = (
        np.asarray(positions_list, dtype=np.float64)
        if positions_list else np.zeros((0, 3), dtype=np.float64)
    )

    mlip_to_omm = {
        mlip_idx: omm_idx
        for omm_idx, mlip_idx in enumerate(omm_to_mlip)
        if mlip_idx is not None
    }
    if len(mlip_to_omm) != n_atoms:
        missing = sorted(set(range(n_atoms)) - set(mlip_to_omm))
        raise RuntimeError(
            f'index map incomplete: {len(mlip_to_omm)}/{n_atoms} MLIP atoms '
            f'mapped (missing {missing[:10]}{"..." if len(missing) > 10 else ""})'
        )

    methods = sorted(charge_methods_used) if charge_methods_used else []
    meta = {
        'translator': 'graph->openff-openmm',
        'forcefield': cfg.forcefield,
        # The method(s) that actually produced the charges — NOT the request —
        # so a NAGL->Gasteiger fallback is recorded truthfully (review
        # 2026-07-17). 'mixed' can occur when a cache carries templates from an
        # earlier build in which NAGL still worked.
        'charge_method': methods[0] if len(methods) == 1 else 'mixed',
        'charge_method_requested': cfg.charge_method,
        'charge_methods_used': methods,
        'nagl_fallback': (
            cfg.charge_method == 'nagl' and 'gasteiger' in charge_methods_used
        ),
        # positions_A are a starting guess (cap H idealized, placeholder H
        # overlapping its former parent): P3 must minimize before dynamics.
        'requires_minimization': True,
        'n_mlip_atoms': n_atoms,
        'n_components': len(components),
        'n_fragments': len(instance_offmols),
        'n_cap_h': n_cap_h,
        'n_placeholder_h': len(placeholder_atoms),
        'n_omm_particles': system.getNumParticles(),
        'cache_hits': cache.hits,
        'cache_misses': cache.misses,
    }

    return ClassicalMix(
        system=system,
        positions_A=positions_out,
        box_vectors_A=cell_A.copy(),
        omm_to_mlip=omm_to_mlip,
        mlip_to_omm=mlip_to_omm,
        n_mlip_atoms=n_atoms,
        metadata=meta,
    )


def _fragment_key(mol) -> str:
    """Canonical graph key for the fragment (RDKit canonical SMILES)."""
    from rdkit import Chem

    return Chem.MolToSmiles(mol)


def _verify_topology_elements(
    off_top, omm_to_mlip: list[int | None], species: list[str],
) -> None:
    """Verify every OpenFF topology atom's element against the index map.

    The prep-stage analogue of decision D-4's atom-order check: particle ``k``
    of the assembled topology must be ``species[omm_to_mlip[k]]`` (or ``'H'``
    for a classical-only cap, ``omm_to_mlip[k] is None``). Any mismatch means
    the fragment assembly reordered/changed atoms and the coordinate write-back
    would silently scramble the system — fail loudly instead (review
    2026-07-17).
    """
    symbols = [atom.symbol for atom in off_top.atoms]
    if len(symbols) != len(omm_to_mlip):
        raise ValueError(
            f'OpenFF topology has {len(symbols)} atoms but the omm<->mlip map '
            f'covers {len(omm_to_mlip)} — fragment assembly changed the atom '
            'count.'
        )
    for k, (symbol, mlip_idx) in enumerate(zip(symbols, omm_to_mlip)):
        expected = 'H' if mlip_idx is None else species[mlip_idx]
        if symbol != expected:
            raise ValueError(
                f'element mismatch at OpenMM particle {k}: topology has '
                f'{symbol!r}, expected {expected!r} '
                f'(mlip index {mlip_idx}) — atom ordering was not preserved.'
            )


def _cap_positions(
    parent: int,
    count: int,
    positions_A: NDArray[np.floating],
    topology: BondTopology,
    n_atoms: int,
) -> NDArray[np.floating]:
    """Deterministic starting coordinates for ``count`` cap H on one centre.

    All caps sit ``_CAP_H_BOND_A`` (Å) from the parent. The vacant-valence axis
    is the negated sum of unit vectors to the parent's existing real neighbours
    (``[0, 0, 1]`` for an isolated heavy atom). One cap goes straight along the
    axis; two or more are tilted ``_MULTI_CAP_TILT_DEG`` off it and spread
    azimuthally (``2*pi/count`` apart) so no two caps coincide — coincident caps
    would give zero cap-parent-cap angles and NaN angle forces (review
    2026-07-17). Pure geometry, no rng: bitwise reproducible.

    Returns a ``(count, 3)`` array in Å.
    """
    p = positions_A[parent]
    dirs = np.zeros(3, dtype=np.float64)
    for nbr in topology.neighbors(parent):
        if nbr < n_atoms:
            d = positions_A[nbr] - p
            norm = np.linalg.norm(d)
            if norm > 1e-9:
                dirs += d / norm
    norm = np.linalg.norm(dirs)
    axis = -dirs / norm if norm > 1e-9 else np.array([0.0, 0.0, 1.0])

    if count == 1:
        return (p + _CAP_H_BOND_A * axis)[None, :]

    # Orthonormal frame (axis, v, w) for the azimuthal spread.
    ref = np.array([1.0, 0.0, 0.0])
    if abs(float(axis @ ref)) > 0.9:
        ref = np.array([0.0, 1.0, 0.0])
    v = ref - float(ref @ axis) * axis
    v /= np.linalg.norm(v)
    w = np.cross(axis, v)

    tilt = np.deg2rad(_MULTI_CAP_TILT_DEG)
    out = np.empty((count, 3), dtype=np.float64)
    for k in range(count):
        az = 2.0 * np.pi * k / count
        direction = (
            np.cos(tilt) * axis
            + np.sin(tilt) * (np.cos(az) * v + np.sin(az) * w)
        )
        out[k] = p + _CAP_H_BOND_A * direction
    return out


def _fragment_positions(
    local_to_global: list[int | None],
    cap_parents: list[int],
    positions_A: NDArray[np.floating],
    topology: BondTopology,
    n_atoms: int,
) -> NDArray[np.floating]:
    """Coordinates (Å) for a fragment's OpenMM atoms, cap H included.

    Real atoms take their MLIP coordinate; cap H are placed near their radical
    parent (see :func:`_cap_positions` — a centre's caps are spread so they
    never coincide).
    """
    out = np.zeros((len(local_to_global), 3), dtype=np.float64)

    # Group consecutive caps by parent so multi-cap centres get spread
    # positions (cap_parents is in append order: same parent's caps adjacent).
    cap_slots = [k for k, g in enumerate(local_to_global) if g is None]
    assert len(cap_slots) == len(cap_parents)
    per_parent: dict[int, list[int]] = {}
    for slot, parent in zip(cap_slots, cap_parents):
        per_parent.setdefault(parent, []).append(slot)

    for k, g in enumerate(local_to_global):
        if g is not None:
            out[k] = positions_A[g]
    for parent, slots in per_parent.items():
        coords = _cap_positions(
            parent, len(slots), positions_A, topology, n_atoms,
        )
        for slot, coord in zip(slots, coords):
            out[slot] = coord
    return out
