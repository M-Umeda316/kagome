"""Shared system builders for polymerization scripts."""
from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from kagome.reactive.groups import PairSpec, ReactiveGroup, ReactionTemplate

# ── SMILES constants ─────────────────────────────────────────────────────────

_MONOMER_SMILES = 'C=CC(=O)OC'   # methyl acrylate
# 1,1-disubstituted vinyls (CH2=CR2): propagating radical is tertiary. Buildable
# since _find_vinyl_alpha_beta accepts a 0-H beta carbon (2026-06-20).
_METHACRYLATE_SMILES = 'C=C(C)C(=O)OC'  # methyl methacrylate (tertiary radical)
_INITIATOR_SMILES = 'CC(C)C#N'   # isobutyronitrile (closed-shell IBN radical model)
_AIBN_SMILES = 'CC(C)(C#N)N=NC(C)(C)C#N'  # azobisisobutyronitrile (intact, for activation)
_DIAMINE_SMILES = 'NCCCCCCN'     # hexamethylenediamine
_DIACID_SMILES = 'OC(=O)CCCCC(=O)O'  # adipic acid
# Bulk epoxy-amine curing (Track 2 E1, decisions.md 2026-07-09). Organic-only:
# the paper's CuO surface is excluded (2026-06-20, OrbMol domain), monomers
# match the paper's resin (Table S3: DGEBA + DETA) and ref#2 (Provenzano 2025).
_DGEBA_SMILES = 'CC(C)(c1ccc(OCC2CO2)cc1)c1ccc(OCC2CO2)cc1'  # bisphenol A diglycidyl ether
_DETA_SMILES = 'NCCNCCN'         # diethylenetriamine (2 primary + 1 secondary N)

from kagome.chem.builders import box_from_density, _rdkit_mol  # noqa: F401 — re-export


def build_ethylene_box(
    n_molecules: int,
    box_size: float,
    rng: np.random.Generator,
    min_sep: float = 3.0,
) -> tuple[np.ndarray, list[str]]:
    """Place n ethylene (C2H4) molecules with atom-level overlap prevention.

    Uses a for-else pattern: raises RuntimeError if valid placement is not
    found in 200 attempts per molecule (instead of silently accepting overlap).
    Overlap check is atom-vs-atom (not center-vs-atom) so min_sep is exact.
    """
    from ase.build import molecule

    ethylene = molecule('C2H4')
    template_pos = ethylene.get_positions()
    template_symbols = list(ethylene.get_chemical_symbols())

    placed_mols: list[np.ndarray] = []
    all_species: list[str] = []
    centered = template_pos - template_pos.mean(axis=0)

    for mol_idx in range(n_molecules):
        for _attempt in range(200):
            offset = rng.uniform(2.0, box_size - 2.0, size=3)
            angle = rng.uniform(0.0, 2.0 * np.pi)
            c, s = np.cos(angle), np.sin(angle)
            rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
            candidate = centered @ rot.T + offset

            ok = True
            for prev in placed_mols:
                diffs = candidate[:, np.newaxis, :] - prev[np.newaxis, :, :]
                if np.min(np.linalg.norm(diffs, axis=2)) < min_sep:
                    ok = False
                    break
            if ok:
                placed_mols.append(candidate)
                all_species.extend(template_symbols)
                break
        else:
            raise RuntimeError(
                f'Could not place molecule {mol_idx + 1}/{n_molecules} without overlap '
                f'in box_size={box_size:.1f} A after 200 attempts. '
                f'Try increasing --box-size.'
            )

    return np.vstack(placed_mols), all_species


def build_template_and_groups(
    n_molecules: int,
) -> tuple[ReactionTemplate, dict[str, ReactiveGroup]]:
    """Define reactive groups for ethylene C=C bond activation."""
    atoms_per_mol = 6  # C2H4: 2C + 4H

    group_a_indices = [i * atoms_per_mol + 0 for i in range(n_molecules)]
    group_b_indices = [i * atoms_per_mol + 1 for i in range(n_molecules)]

    template = ReactionTemplate(
        name='vinyl_polymerization',
        groups=['C_donor', 'C_acceptor'],
        pairs=[
            PairSpec(
                group_a='C_donor',
                group_b='C_acceptor',
                is_formation=True,
                r_min=3.0,
                r_max=6.0,
            ),
        ],
    )

    groups = {
        'C_donor': ReactiveGroup('C_donor', group_a_indices),
        'C_acceptor': ReactiveGroup('C_acceptor', group_b_indices),
    }

    return template, groups


# ── RDKit-based helpers ──────────────────────────────────────────────────────
# _rdkit_mol and box_from_density live in src/chem/builders.py (canonical
# location) and are re-exported at the top of this file for backward compat.


def _rdkit_3d(smiles: str, seed: int = 42) -> tuple[np.ndarray, list[str]]:
    """SMILES → (positions Å, species list) via RDKit EmbedMolecule + MMFF.

    Thin wrapper over `_rdkit_mol` so the builder and the classical prep stage
    share one atom-ordering source (see `_rdkit_mol`).
    """
    mol = _rdkit_mol(smiles, seed)
    conf = mol.GetConformer()
    positions = np.array(conf.GetPositions(), dtype=np.float64)
    species = [atom.GetSymbol() for atom in mol.GetAtoms()]
    return positions, species


def _rdkit_local_bonds(smiles: str, seed: int = 42) -> list[tuple[int, int, float]]:
    """Return the molecule's intramolecular bonds as ``(i, j, order)`` in the
    SAME atom ordering as ``_rdkit_3d`` (AddHs + embed).

    Connectivity is embedding-independent, but we build from the same
    ``_rdkit_mol`` so local indices align 1:1 with the builder's positions and
    the group index arithmetic.  Order is RDKit's ``GetBondTypeAsDouble``
    (1.0 single, 1.5 aromatic, 2.0 double, 3.0 triple).
    """
    mol = _rdkit_mol(smiles, seed)
    bonds: list[tuple[int, int, float]] = []
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bonds.append((int(i), int(j), float(b.GetBondTypeAsDouble())))
    return bonds


def layout_bonds(
    specs: list[tuple[str, int, int]],
) -> list[tuple[int, int, float]]:
    """Global intramolecular bond topology for a contiguous fragment layout.

    ``specs`` is the ordered ``[(smiles, count, rdkit_seed), ...]`` of the
    placement — the SAME order the builders use in ``_place_fragments_in_box``
    (fragments laid out contiguously, each occupying ``len(AddHs(mol))`` atoms).
    Returns ``(i, j, order)`` bonds in global indices that line up 1:1 with the
    builder's ``positions``/``groups``.  Non-breaking companion for trajectory
    topology output (specs/decisions.md 2026-07-02); a consistency test asserts
    the layout matches the builders.
    """
    bonds: list[tuple[int, int, float]] = []
    offset = 0
    for smiles, count, seed in specs:
        local = _rdkit_local_bonds(smiles, seed=seed)
        n_per = len(_rdkit_mol(smiles, seed).GetAtoms())
        for k in range(count):
            base = offset + k * n_per
            bonds.extend((base + a, base + b, o) for a, b, o in local)
        offset += count * n_per
    return bonds


def layout_species(specs: list[tuple[str, int, int]]) -> list[str]:
    """Element symbols for a contiguous fragment layout (same order as
    :func:`layout_bonds`).  Placement-free — used to validate a retroactive
    reconstruction against a run's recorded ``species`` without re-placing atoms."""
    out: list[str] = []
    for smiles, count, seed in specs:
        sp = [a.GetSymbol() for a in _rdkit_mol(smiles, seed).GetAtoms()]
        out.extend(sp * count)
    return out


def vinyl_initial_bonds(
    n_monomers: int,
    n_initiators: int,
    monomer_smiles: str = _MONOMER_SMILES,
    initiator_smiles: str = _INITIATOR_SMILES,
    rdkit_seed: int = 42,
) -> list[tuple[int, int, float]]:
    """Initial topology for :func:`build_vinyl_aibn_system` (pre-formed radicals:
    initiators first, then monomers)."""
    return layout_bonds([
        (initiator_smiles, n_initiators, rdkit_seed),
        (monomer_smiles, n_monomers, rdkit_seed + 1),
    ])


def vinyl_reaction_maps(
    n_monomers: int,
    n_initiators: int,
    monomer_smiles: str = _MONOMER_SMILES,
    initiator_smiles: str = _INITIATOR_SMILES,
    rdkit_seed: int = 42,
) -> tuple[dict[int, int], list[tuple[int, int]]]:
    """Placement-free (propagation_map alpha->beta, azo_bonds) for the
    :func:`build_vinyl_aibn_system` layout.  Pre-formed radicals have no azo
    bonds, so the second element is empty.  For retroactive topology
    reconstruction without re-placing atoms in a box (specs/decisions.md
    2026-07-02)."""
    local_alpha, local_beta = _find_vinyl_alpha_beta(monomer_smiles)
    n_per_init = len(_rdkit_mol(initiator_smiles, rdkit_seed).GetAtoms())
    n_per_mono = len(_rdkit_mol(monomer_smiles, rdkit_seed + 1).GetAtoms())
    mono_offset = n_initiators * n_per_init
    pmap = {
        mono_offset + j * n_per_mono + local_alpha:
        mono_offset + j * n_per_mono + local_beta
        for j in range(n_monomers)
    }
    return pmap, []


def full_aibn_reaction_maps(
    n_monomers: int,
    n_aibn: int,
    monomer_smiles: str = _MONOMER_SMILES,
    aibn_smiles: str = _AIBN_SMILES,
    rdkit_seed: int = 42,
) -> tuple[dict[int, int], list[tuple[int, int]]]:
    """Placement-free (propagation_map alpha->beta, global azo C-N bonds) for the
    :func:`build_full_aibn_system` layout (AIBN first, then monomers).  For
    retroactive topology reconstruction without re-placing atoms in a box."""
    local_alpha, local_beta = _find_vinyl_alpha_beta(monomer_smiles)
    local_azo = _find_aibn_azo_bonds(aibn_smiles)
    n_per_aibn = len(_rdkit_mol(aibn_smiles, rdkit_seed).GetAtoms())
    n_per_mono = len(_rdkit_mol(monomer_smiles, rdkit_seed + 1).GetAtoms())
    azo = [
        (i * n_per_aibn + c, i * n_per_aibn + nn)
        for i in range(n_aibn) for c, nn in local_azo
    ]
    mono_offset = n_aibn * n_per_aibn
    pmap = {
        mono_offset + j * n_per_mono + local_alpha:
        mono_offset + j * n_per_mono + local_beta
        for j in range(n_monomers)
    }
    return pmap, azo


def full_aibn_initial_bonds(
    n_monomers: int,
    n_aibn: int,
    monomer_smiles: str = _MONOMER_SMILES,
    aibn_smiles: str = _AIBN_SMILES,
    rdkit_seed: int = 42,
) -> list[tuple[int, int, float]]:
    """Initial (intact-AIBN) topology for :func:`build_full_aibn_system` (AIBN
    first, then monomers).  The two azo C-N bonds per AIBN are removed from the
    tracked topology after activation dissociates them (run_vinyl_aibn)."""
    return layout_bonds([
        (aibn_smiles, n_aibn, rdkit_seed),
        (monomer_smiles, n_monomers, rdkit_seed + 1),
    ])


def _find_vinyl_alpha_beta(smiles: str) -> tuple[int, int]:
    """Return (alpha_idx, beta_idx) local atom indices for the CH2=CR- motif.

    alpha-C: =CH2 terminus (2 H) → Group J; the radical attacks here.
    beta-C:  the other vinyl carbon → propagation target (becomes the new
             radical). 1 H for mono-substituted vinyls (CH2=CH-R, e.g. methyl
             acrylate, styrene); 0 H for 1,1-disubstituted vinyls (CH2=CR2,
             e.g. methacrylate, diphenylethylene, dimethyl itaconate), where the
             propagating radical is tertiary.

    Regiochemistry (head-to-tail): the radical always adds to the terminal =CH2
    (alpha) and the unpaired electron localizes on the substituted carbon (beta).
    Tertiary-radical stability for 1,1-disubstituted monomers is left to the MLIP
    energetics (decision 2026-06-20 "atom-typing only"); no extra bias is added.
    Aromatic ring bonds report bond order 1.5, so only the genuine vinyl double
    bond matches — styrene/diphenylethylene phenyl rings are not mistaken for it.
    """
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() != 2.0:
            continue
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if a.GetSymbol() != 'C' or b.GetSymbol() != 'C':
            continue
        h_a = sum(1 for n in a.GetNeighbors() if n.GetSymbol() == 'H')
        h_b = sum(1 for n in b.GetNeighbors() if n.GetSymbol() == 'H')
        # alpha must be the =CH2 terminus (2 H); beta is the substituted carbon
        # (1 H mono-substituted, 0 H for 1,1-disubstituted).
        if h_a == 2 and h_b <= 1:
            return a.GetIdx(), b.GetIdx()
        if h_b == 2 and h_a <= 1:
            return b.GetIdx(), a.GetIdx()
    raise ValueError(f'No CH2=CR- vinyl terminus found in {smiles!r}')


def _find_ibn_radical_c(smiles: str) -> int:
    """Return local index of the radical-bearing C in the isobutyronitrile model.

    In the closed-shell representation CC(C)C#N (isobutyronitrile), the radical-
    centre C is bonded to 3 C-neighbours (2 methyls + nitrile C) and 1 H.
    This is the unique C with exactly 3 C-neighbours in the molecule.
    """
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'C':
            continue
        nbrs = atom.GetNeighbors()
        n_c = sum(1 for n in nbrs if n.GetSymbol() == 'C')
        if n_c == 3:
            return atom.GetIdx()
    raise ValueError(f'No C bonded to exactly 3 C-neighbours found in {smiles!r}')


def _find_chain_c_neighbor(smiles: str, radical_idx: int) -> int:
    """Return local index of a C neighbor of the radical C (group k, Table S1).

    Picks the first non-nitrile C neighbor (a methyl C). The nitrile C is bonded
    to N via triple bond, so we exclude it.
    """
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    radical = mol.GetAtomWithIdx(radical_idx)
    for nbr in radical.GetNeighbors():
        if nbr.GetSymbol() != 'C':
            continue
        has_triple = any(
            b.GetBondTypeAsDouble() == 3.0 for b in nbr.GetBonds()
        )
        if not has_triple:
            return nbr.GetIdx()
    raise ValueError(
        f'No non-nitrile C neighbor of radical C (idx {radical_idx}) in {smiles!r}'
    )


def _rotation_matrix(axis: np.ndarray, angle: float) -> np.ndarray:
    """Rodrigues rotation matrix for arbitrary axis/angle."""
    c, s = np.cos(angle), np.sin(angle)
    x, y, z = axis / np.linalg.norm(axis)
    return np.array([
        [c + x * x * (1 - c),     x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
        [y * x * (1 - c) + z * s, c + y * y * (1 - c),     y * z * (1 - c) - x * s],
        [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)    ],
    ], dtype=np.float64)


def _place_fragments_in_box(
    fragments: list[tuple[np.ndarray, list[str]]],
    box_size: float,
    rng: np.random.Generator,
    min_sep: float = 2.5,
    max_attempts: int = 500,
) -> tuple[np.ndarray, list[str]]:
    """Place a list of (positions, species) molecules in a periodic box.

    Uses grid-guided placement: each molecule is seeded near a distinct grid
    cell centre, then given a random 3-D rotation + small jitter and accepted on
    the first non-overlapping pose. Localising the search this way makes
    placement reliable at paper density (0.5 g/mL), where global rejection
    sampling stalls. The grid is purely an initial-configuration device and does
    not bias the subsequent dynamics. See specs/decisions.md
    "2026-06-13: grid-guided initial placement at paper density".

    Raises RuntimeError if any molecule cannot be placed without overlap.
    """
    placed: list[np.ndarray] = []
    all_species: list[str] = []

    n_mols = len(fragments)
    # Grid sized for ~70% cell occupancy (headroom factor 0.7) so each molecule's
    # assigned cell is likely to have empty neighbours. The earlier ncells =
    # ceil(N^(1/3)) gave only ~N cells: at high counts (e.g. 210 molecules -> 216
    # cells, 97% occupancy) late molecules become boxed in by filled neighbours
    # and placement fails despite the density being physically feasible. The 0.7
    # factor leaves ncells unchanged for small systems (e.g. 42 -> 4) and only
    # adds cells where crowding would otherwise stall the search. Placement is a
    # non-physical initial device and does not bias the dynamics (see
    # specs/decisions.md 2026-06-13 grid-placement record).
    ncells = int(np.ceil((n_mols / 0.7) ** (1.0 / 3.0)))
    cell = box_size / ncells
    cell_centers = np.array(
        [
            ((i + 0.5) * cell, (j + 0.5) * cell, (k + 0.5) * cell)
            for i in range(ncells)
            for j in range(ncells)
            for k in range(ncells)
        ]
    )
    rng.shuffle(cell_centers)

    box = np.array([box_size, box_size, box_size])

    for mol_idx, (template_pos, mol_species) in enumerate(fragments):
        centered = template_pos - template_pos.mean(axis=0)
        base = cell_centers[mol_idx]
        for _attempt in range(max_attempts):
            axis = rng.standard_normal(3)
            if np.linalg.norm(axis) < 1e-10:
                axis = np.array([0.0, 0.0, 1.0])
            R = _rotation_matrix(axis, rng.uniform(0.0, 2.0 * np.pi))
            # Jitter shrinks with attempt count so later tries hug the cell centre.
            jitter = rng.uniform(-0.5, 0.5, 3) * cell * (1.0 - _attempt / max_attempts)
            offset = (base + jitter) % box_size
            candidate = centered @ R.T + offset

            ok = True
            for prev in placed:
                diffs = candidate[:, np.newaxis, :] - prev[np.newaxis, :, :]
                diffs = diffs - box * np.round(diffs / box)
                if np.min(np.linalg.norm(diffs, axis=2)) < min_sep:
                    ok = False
                    break
            if ok:
                placed.append(candidate)
                all_species.extend(mol_species)
                break
        else:
            raise RuntimeError(
                f'Could not place molecule {mol_idx + 1}/{len(fragments)} '
                f'in box_size={box_size:.1f} A after {max_attempts} attempts. '
                f'Try increasing --box-size or lowering --density.'
            )

    return np.vstack(placed), all_species


def build_vinyl_aibn_system(
    n_monomers: int,
    n_initiators: int,
    box_size: float,
    rng: np.random.Generator,
    monomer_smiles: str = _MONOMER_SMILES,
    initiator_smiles: str = _INITIATOR_SMILES,
    min_sep: float = 2.5,
    rdkit_seed: int = 42,
) -> tuple[np.ndarray, list[str], ReactionTemplate, dict[str, ReactiveGroup],
           dict[int, int], dict[int, int]]:
    """Build a vinyl + AIBN-initiator system for radical polymerization.

    Placement order: all initiators first, then all monomers.
    This keeps global index arithmetic straightforward.

    Returns
    -------
    positions       : (N, 3) Å array
    species         : list of element symbols
    template        : ReactionTemplate (4-group: radical_C, vinyl_alpha_C, chain_C, vinyl_beta_C)
    groups          : dict of ReactiveGroup (4 groups)
    propagation_map : {alpha_C_global_idx: beta_C_global_idx} for each monomer
    chain_c_map     : {radical_C_global_idx: chain_C_global_idx} for each radical

    Paper anchor: Fig. 1, Section 2, Table S1 (ij+ik+jl criterion).
    Design: specs/decisions.md — "T-G1: vinyl radical polymerization system".

    NOTE (closed-shell approximation, specs/decisions.md 2026-07-02): the default
    initiator is isobutyronitrile (``_INITIATOR_SMILES``), a *closed-shell* model
    whose radical carbon carries a placeholder H standing in for the unpaired
    electron.  On the first addition that carbon would become 5-coordinate in the
    raw geometry; the emitted topology idealizes this by shedding the placeholder
    H (see reactive/topology.apply_vinyl_addition), but the H atom remains in the
    coordinates.  For valence-faithful initiator chemistry (genuine 3-coordinate
    radicals, no placeholder H) use the ``--activation`` path
    (:func:`build_full_aibn_system`), which decomposes intact AIBN.
    """
    # ── generate single-molecule 3D templates ──
    init_pos, init_sp = _rdkit_3d(initiator_smiles, seed=rdkit_seed)
    mono_pos, mono_sp = _rdkit_3d(monomer_smiles, seed=rdkit_seed + 1)

    local_radical = _find_ibn_radical_c(initiator_smiles)
    local_chain_c = _find_chain_c_neighbor(initiator_smiles, local_radical)
    local_alpha, local_beta = _find_vinyl_alpha_beta(monomer_smiles)

    n_per_init = len(init_sp)
    n_per_mono = len(mono_sp)

    # ── place in box: initiators first, then monomers ──
    fragments = (
        [(init_pos, init_sp)] * n_initiators
        + [(mono_pos, mono_sp)] * n_monomers
    )
    positions, species = _place_fragments_in_box(
        fragments, box_size, rng, min_sep=min_sep,
    )

    # ── global index arithmetic ──
    init_offset = 0
    mono_offset = n_initiators * n_per_init

    radical_C_indices = [
        init_offset + i * n_per_init + local_radical
        for i in range(n_initiators)
    ]
    chain_C_indices = [
        init_offset + i * n_per_init + local_chain_c
        for i in range(n_initiators)
    ]
    alpha_C_indices = [
        mono_offset + j * n_per_mono + local_alpha
        for j in range(n_monomers)
    ]
    beta_C_indices = [
        mono_offset + j * n_per_mono + local_beta
        for j in range(n_monomers)
    ]
    propagation_map: dict[int, int] = {
        mono_offset + j * n_per_mono + local_alpha:
        mono_offset + j * n_per_mono + local_beta
        for j in range(n_monomers)
    }
    chain_c_map: dict[int, int] = {
        radical_C_indices[i]: chain_C_indices[i]
        for i in range(n_initiators)
    }

    # ── reaction template (Table S1: 4-group ij+ik+jl criterion) ──
    template = ReactionTemplate(
        name='radical_vinyl_polymerization',
        groups=['radical_C', 'vinyl_alpha_C', 'chain_C', 'vinyl_beta_C'],
        pairs=[
            PairSpec(
                group_a='radical_C',
                group_b='vinyl_alpha_C',
                is_formation=True,
                r_min=3.0,
                r_max=6.0,
            ),
            PairSpec(
                group_a='radical_C',
                group_b='chain_C',
                is_formation=True,
                r_min=0.0,
                r_max=3.0,
                constraint_only=True,
            ),
            PairSpec(
                group_a='vinyl_alpha_C',
                group_b='vinyl_beta_C',
                is_formation=True,
                r_min=0.0,
                r_max=3.0,
                constraint_only=True,
            ),
        ],
    )
    groups = {
        'radical_C':      ReactiveGroup('radical_C',      radical_C_indices),
        'vinyl_alpha_C':  ReactiveGroup('vinyl_alpha_C',  alpha_C_indices),
        'chain_C':        ReactiveGroup('chain_C',        chain_C_indices),
        'vinyl_beta_C':   ReactiveGroup('vinyl_beta_C',   beta_C_indices),
    }

    return positions, species, template, groups, propagation_map, chain_c_map


def build_vinyl_copolymer_system(
    monomer_specs: list[tuple[str, int]],
    n_initiators: int,
    box_size: float,
    rng: np.random.Generator,
    initiator_smiles: str = _INITIATOR_SMILES,
    min_sep: float = 2.5,
    rdkit_seed: int = 42,
) -> tuple[np.ndarray, list[str], ReactionTemplate, dict[str, ReactiveGroup],
           dict[int, int], dict[int, int]]:
    """Build a MIXED-monomer vinyl + initiator system (copolymerization).

    Generalizes :func:`build_vinyl_aibn_system` to several distinct monomer
    species in one box — e.g. methyl acrylate + methyl methacrylate. Unlike the
    single-species builder, atom counts differ between species, so global index
    arithmetic uses a *running offset accumulated in placement order* rather than
    a single fixed per-monomer stride (specs/decisions.md 2026-07-16).

    Both species' alpha carbons are registered into the SAME ``vinyl_alpha_C``
    group and both beta carbons into the SAME ``vinyl_beta_C`` group, so the
    radical adds to either monomer with no extra selection bias — the
    copolymer sequence is left to the MLIP energetics (paper-faithful:
    §2/Table S1 add no per-addition bias; no reactivity-ratio control).

    Parameters
    ----------
    monomer_specs : list of (smiles, count)
        Each distinct monomer species and how many copies to place. Order is
        the placement order (after the initiators).
    n_initiators : int
        Number of isobutyronitrile radical models (closed-shell approximation,
        see :func:`build_vinyl_aibn_system` NOTE).

    Placement order: all initiators first, then each monomer spec in the given
    order (spec[0] × count[0], spec[1] × count[1], …).

    Returns the SAME tuple shape as :func:`build_vinyl_aibn_system`:
    ``(positions, species, template, groups, propagation_map, chain_c_map)``.

    Paper anchor: Fig. 1, Section 2, Table S1 (ij+ik+jl criterion). The mixed
    system is a documented extension of the paper's per-monomer setup.
    """
    if not monomer_specs:
        raise ValueError('monomer_specs must contain at least one (smiles, count).')
    if any(count < 0 for _, count in monomer_specs):
        raise ValueError('monomer counts must be non-negative.')

    # ── single-molecule 3D templates ──
    init_pos, init_sp = _rdkit_3d(initiator_smiles, seed=rdkit_seed)
    local_radical = _find_ibn_radical_c(initiator_smiles)
    local_chain_c = _find_chain_c_neighbor(initiator_smiles, local_radical)

    # One template + local alpha/beta per DISTINCT species (geometry seed varies
    # per species; atom ordering is deterministic per SMILES).
    mono_templates: list[tuple[np.ndarray, list[str], int, int]] = []
    for k, (smiles, _count) in enumerate(monomer_specs):
        mpos, msp = _rdkit_3d(smiles, seed=rdkit_seed + 1 + k)
        local_alpha, local_beta = _find_vinyl_alpha_beta(smiles)
        mono_templates.append((mpos, msp, local_alpha, local_beta))

    # ── build the placement list AND parallel per-fragment metadata ──
    # frag_meta[i] = ('init', None, None) or ('mono', local_alpha, local_beta),
    # walked in the SAME order the placer concatenates species.
    fragments: list[tuple[np.ndarray, list[str]]] = [(init_pos, init_sp)] * n_initiators
    frag_meta: list[tuple[str, int | None, int | None]] = [('init', None, None)] * n_initiators
    for k, (_smiles, count) in enumerate(monomer_specs):
        mpos, msp, local_alpha, local_beta = mono_templates[k]
        fragments += [(mpos, msp)] * count
        frag_meta += [('mono', local_alpha, local_beta)] * count

    positions, species = _place_fragments_in_box(
        fragments, box_size, rng, min_sep=min_sep,
    )

    # ── global index arithmetic via running offset (heterogeneous atom counts) ──
    radical_C_indices: list[int] = []
    chain_C_indices: list[int] = []
    alpha_C_indices: list[int] = []
    beta_C_indices: list[int] = []
    propagation_map: dict[int, int] = {}
    chain_c_map: dict[int, int] = {}

    offset = 0
    for (_fpos, fsp), (kind, local_alpha, local_beta) in zip(fragments, frag_meta):
        if kind == 'init':
            radical = offset + local_radical
            chain = offset + local_chain_c
            radical_C_indices.append(radical)
            chain_C_indices.append(chain)
            chain_c_map[radical] = chain
        else:
            alpha = offset + local_alpha
            beta = offset + local_beta
            alpha_C_indices.append(alpha)
            beta_C_indices.append(beta)
            propagation_map[alpha] = beta
        offset += len(fsp)

    # ── reaction template (identical to the single-monomer 4-group template) ──
    template = ReactionTemplate(
        name='radical_vinyl_copolymerization',
        groups=['radical_C', 'vinyl_alpha_C', 'chain_C', 'vinyl_beta_C'],
        pairs=[
            PairSpec(
                group_a='radical_C',
                group_b='vinyl_alpha_C',
                is_formation=True,
                r_min=3.0,
                r_max=6.0,
            ),
            PairSpec(
                group_a='radical_C',
                group_b='chain_C',
                is_formation=True,
                r_min=0.0,
                r_max=3.0,
                constraint_only=True,
            ),
            PairSpec(
                group_a='vinyl_alpha_C',
                group_b='vinyl_beta_C',
                is_formation=True,
                r_min=0.0,
                r_max=3.0,
                constraint_only=True,
            ),
        ],
    )
    groups = {
        'radical_C':      ReactiveGroup('radical_C',      radical_C_indices),
        'vinyl_alpha_C':  ReactiveGroup('vinyl_alpha_C',  alpha_C_indices),
        'chain_C':        ReactiveGroup('chain_C',        chain_C_indices),
        'vinyl_beta_C':   ReactiveGroup('vinyl_beta_C',   beta_C_indices),
    }

    return positions, species, template, groups, propagation_map, chain_c_map


def copolymer_alpha_species(
    monomer_specs: list[tuple[str, int]],
    n_initiators: int,
    initiator_smiles: str = _INITIATOR_SMILES,
) -> dict[int, str]:
    """Map each monomer's alpha-C global index to its species SMILES.

    Reconstructs the SAME global indices as :func:`build_vinyl_copolymer_system`
    without needing the rng or placement — the offsets depend only on the atom
    counts and the placement order (initiators first, then each monomer spec in
    order), both deterministic per SMILES. Used by the reactivity analysis to
    label which species each confirmed formation incorporated
    (specs/decisions.md 2026-07-16).

    Returns ``{alpha_C_global_idx: monomer_smiles}``.
    """
    n_per_init = len(_rdkit_3d(initiator_smiles)[1])

    alpha_species: dict[int, str] = {}
    offset = n_initiators * n_per_init
    for smiles, count in monomer_specs:
        _pos, sp = _rdkit_3d(smiles)
        n_per_mono = len(sp)
        local_alpha, _local_beta = _find_vinyl_alpha_beta(smiles)
        for _ in range(count):
            alpha_species[offset + local_alpha] = smiles
            offset += n_per_mono
    return alpha_species


def copolymer_atom_species(
    monomer_specs: list[tuple[str, int]],
    n_initiators: int,
    initiator_smiles: str = _INITIATOR_SMILES,
) -> dict[int, str]:
    """Map EVERY atom's global index to the SMILES of its molecule block.

    Same deterministic offsets as :func:`build_vinyl_copolymer_system` and
    :func:`copolymer_alpha_species` (initiators first, then each monomer spec
    in order). Initiator atoms map to ``initiator_smiles``. Used to classify
    the radical_C endpoint of a confirmed formation: the block holding the
    radical is the chain's terminal unit, giving the cross-propagation table
    (which terminal species added which monomer; specs/decisions.md 2026-07-16).

    Returns ``{atom_global_idx: block_smiles}`` covering all atoms.
    """
    atom_species: dict[int, str] = {}
    offset = 0
    n_per_init = len(_rdkit_3d(initiator_smiles)[1])
    for _ in range(n_initiators):
        for i in range(n_per_init):
            atom_species[offset + i] = initiator_smiles
        offset += n_per_init
    for smiles, count in monomer_specs:
        n_per_mono = len(_rdkit_3d(smiles)[1])
        for _ in range(count):
            for i in range(n_per_mono):
                atom_species[offset + i] = smiles
            offset += n_per_mono
    return atom_species


def copolymer_initial_bonds(
    monomer_specs: list[tuple[str, int]],
    n_initiators: int,
    initiator_smiles: str = _INITIATOR_SMILES,
) -> list[tuple[int, int, float]]:
    """Full initial intramolecular bond topology for
    :func:`build_vinyl_copolymer_system`.

    Same deterministic running-offset convention as :func:`copolymer_alpha_species`
    / :func:`copolymer_atom_species`: initiators first (``n_initiators`` copies),
    then each ``monomer_specs`` entry in placement order, each repeated ``count``
    times. Delegates the per-fragment RDKit bond extraction to :func:`layout_bonds`
    (already running-offset-correct across heterogeneous specs) rather than
    duplicating it. The seed passed to ``layout_bonds`` only affects the discarded
    3-D conformer, not atom indexing/connectivity (see ``_rdkit_mol``), so a fixed
    seed reproduces the same bonds. Enables trajectory topology tracking for the
    copolymer driver (specs/decisions.md 2026-07-17 "well-mixed 測定モード" 前提工事).
    """
    specs = [(initiator_smiles, n_initiators, 42)]
    specs += [(smiles, count, 42) for smiles, count in monomer_specs]
    return layout_bonds(specs)


# ── nylon-6,6 step-growth system ────────────────────────────────────────────

def _find_terminal_amine_n(smiles: str) -> list[int]:
    """Return local indices of terminal primary amine N in the molecule."""
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    indices = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'N':
            continue
        h_count = sum(1 for n in atom.GetNeighbors() if n.GetSymbol() == 'H')
        if h_count >= 2:
            indices.append(atom.GetIdx())
    return indices


def _find_amine_h(smiles: str, n_idx: int) -> int:
    """Return local index of one H bonded to amine N at n_idx."""
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    atom = mol.GetAtomWithIdx(n_idx)
    for nbr in atom.GetNeighbors():
        if nbr.GetSymbol() == 'H':
            return nbr.GetIdx()
    raise ValueError(f'No H bonded to N at index {n_idx} in {smiles!r}')


def _find_carboxyl_c_and_oh(smiles: str) -> list[tuple[int, int]]:
    """Return [(carboxyl_C_idx, hydroxyl_O_idx), ...] for each -C(=O)OH."""
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    results = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'C':
            continue
        o_double = None
        o_single = None
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() != 'O':
                continue
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if bond.GetBondTypeAsDouble() == 2.0:
                o_double = nbr.GetIdx()
            elif bond.GetBondTypeAsDouble() == 1.0:
                h_on_o = sum(1 for n in nbr.GetNeighbors() if n.GetSymbol() == 'H')
                if h_on_o >= 1:
                    o_single = nbr.GetIdx()
        if o_double is not None and o_single is not None:
            results.append((atom.GetIdx(), o_single))
    return results


def build_nylon66_system(
    n_diamines: int,
    n_diacids: int,
    box_size: float,
    rng: np.random.Generator,
    diamine_smiles: str = _DIAMINE_SMILES,
    diacid_smiles: str = _DIACID_SMILES,
    min_sep: float = 2.5,
    rdkit_seed: int = 42,
) -> tuple[np.ndarray, list[str], ReactionTemplate, dict[str, ReactiveGroup]]:
    """Build a nylon-6,6 step-growth polycondensation system.

    4-group template matching Table S2 of arXiv:2511.22874.
    Groups: amine_N (i), carboxyl_C (j), amine_H (k), carboxyl_OH (l)

    Returns (positions, species, template, groups).

    Paper anchor: PDF p.22, Table S2, Fig. S2.
    """
    diamine_pos, diamine_sp = _rdkit_3d(diamine_smiles, seed=rdkit_seed)
    diacid_pos, diacid_sp = _rdkit_3d(diacid_smiles, seed=rdkit_seed + 1)

    local_amine_ns = _find_terminal_amine_n(diamine_smiles)
    local_carboxyl_pairs = _find_carboxyl_c_and_oh(diacid_smiles)

    n_per_diamine = len(diamine_sp)
    n_per_diacid = len(diacid_sp)

    fragments = (
        [(diamine_pos, diamine_sp)] * n_diamines
        + [(diacid_pos, diacid_sp)] * n_diacids
    )
    positions, species = _place_fragments_in_box(
        fragments, box_size, rng, min_sep=min_sep,
    )

    diamine_offset = 0
    diacid_offset = n_diamines * n_per_diamine

    amine_N_indices: list[int] = []
    amine_H_indices: list[int] = []
    for i in range(n_diamines):
        base = diamine_offset + i * n_per_diamine
        for n_idx in local_amine_ns:
            amine_N_indices.append(base + n_idx)
            h_idx = _find_amine_h(diamine_smiles, n_idx)
            amine_H_indices.append(base + h_idx)

    carboxyl_C_indices: list[int] = []
    carboxyl_OH_indices: list[int] = []
    for j in range(n_diacids):
        base = diacid_offset + j * n_per_diacid
        for c_idx, oh_idx in local_carboxyl_pairs:
            carboxyl_C_indices.append(base + c_idx)
            carboxyl_OH_indices.append(base + oh_idx)

    template = ReactionTemplate(
        name='nylon66_condensation',
        groups=['amine_N', 'carboxyl_C', 'amine_H', 'carboxyl_OH'],
        pairs=[
            PairSpec(
                group_a='amine_N', group_b='carboxyl_C',
                is_formation=True, r_min=3.0, r_max=6.0,
            ),
            PairSpec(
                group_a='amine_N', group_b='amine_H',
                is_formation=False, r_min=0.0, r_max=3.0,
            ),
            PairSpec(
                group_a='carboxyl_C', group_b='carboxyl_OH',
                is_formation=False, r_min=0.0, r_max=3.0,
            ),
            PairSpec(
                group_a='amine_H', group_b='carboxyl_OH',
                is_formation=True, r_min=0.0, r_max=100.0,
                score_pair=False, count_as_reaction=False,
            ),
        ],
    )
    groups = {
        'amine_N':      ReactiveGroup('amine_N',      amine_N_indices),
        'carboxyl_C':   ReactiveGroup('carboxyl_C',   carboxyl_C_indices),
        'amine_H':      ReactiveGroup('amine_H',      amine_H_indices),
        'carboxyl_OH':  ReactiveGroup('carboxyl_OH',  carboxyl_OH_indices),
    }

    return positions, species, template, groups


# ── bulk epoxy-amine curing system (Track 2 E1) ─────────────────────────────

def _find_epoxide_c_o(smiles: str) -> list[tuple[int, int]]:
    """Return [(terminal_C_idx, ring_O_idx), ...] for each epoxide ring.

    terminal_C is the less-substituted (more-H) ring carbon — the carbon a
    primary/secondary amine attacks in the ring-opening addition (decisions.md
    2026-07-08 E0 design; E0 scan verified this channel on OrbMol 2026-07-09).
    """
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    patt = Chem.MolFromSmarts('C1CO1')
    results = []
    for match in mol.GetSubstructMatches(patt):
        c_atoms = [i for i in match if mol.GetAtomWithIdx(i).GetSymbol() == 'C']
        o_atoms = [i for i in match if mol.GetAtomWithIdx(i).GetSymbol() == 'O']
        if len(c_atoms) != 2 or len(o_atoms) != 1:
            continue

        def _h_count(idx: int) -> int:
            return sum(1 for n in mol.GetAtomWithIdx(idx).GetNeighbors()
                       if n.GetSymbol() == 'H')

        terminal = max(c_atoms, key=_h_count)
        results.append((terminal, o_atoms[0]))
    return results


def _find_amine_n_h(smiles: str) -> list[tuple[int, list[int]]]:
    """Return [(N_idx, [H_idx, ...]), ...] for every sp3 amine N with >=1 H.

    Both primary (2 H) and secondary (1 H) amines are reactive toward epoxide
    ring opening (1° -> 2° -> 3° multi-addition), so all N-H hydrogens are
    registered — the group updater retires the N only when its last H is gone.
    """
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    results = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'N' or atom.GetIsAromatic():
            continue
        h_indices = [n.GetIdx() for n in atom.GetNeighbors()
                     if n.GetSymbol() == 'H']
        if h_indices:
            results.append((atom.GetIdx(), h_indices))
    return results


def build_epoxy_amine_system(
    n_epoxies: int,
    n_amines: int,
    box_size: float,
    rng: np.random.Generator,
    epoxy_smiles: str = _DGEBA_SMILES,
    amine_smiles: str = _DETA_SMILES,
    min_sep: float = 2.5,
    rdkit_seed: int = 42,
) -> tuple[np.ndarray, list[str], ReactionTemplate, dict[str, ReactiveGroup],
           dict[int, list[int]]]:
    """Build a bulk epoxy-amine curing system (ring-opening addition).

    4-group template structurally identical to nylon condensation (Table S2
    machinery; decisions.md 2026-07-08 E0 design / 2026-07-09 E0 PROCEED):
    Groups: amine_N (i), epoxy_C (j), amine_H (k), ring_O (l)
    Pairs:  (i,j) N-C formation / (i,k) N-H dissociation /
            (j,l) ring C-O dissociation / (k,l) hydroxyl O-H, bias-only.
    Unlike nylon, no atom leaves the molecule: the ring O stays bonded to the
    beta carbon and becomes the beta-hydroxyl.

    Returns (positions, species, template, groups, amine_h_map) where
    amine_h_map maps each global amine-N index to the global indices of ALL its
    N-H hydrogens (consumed one per addition; drives 1° -> 2° -> 3° group
    reassignment in EpoxyAmineAdditionUpdater).
    """
    epoxy_pos, epoxy_sp = _rdkit_3d(epoxy_smiles, seed=rdkit_seed)
    amine_pos, amine_sp = _rdkit_3d(amine_smiles, seed=rdkit_seed + 1)

    local_epoxide_pairs = _find_epoxide_c_o(epoxy_smiles)
    local_amine_nh = _find_amine_n_h(amine_smiles)
    if not local_epoxide_pairs:
        raise ValueError(f'No epoxide ring found in {epoxy_smiles!r}')
    if not local_amine_nh:
        raise ValueError(f'No amine N-H found in {amine_smiles!r}')

    n_per_epoxy = len(epoxy_sp)
    n_per_amine = len(amine_sp)

    fragments = (
        [(epoxy_pos, epoxy_sp)] * n_epoxies
        + [(amine_pos, amine_sp)] * n_amines
    )
    positions, species = _place_fragments_in_box(
        fragments, box_size, rng, min_sep=min_sep,
    )

    epoxy_offset = 0
    amine_offset = n_epoxies * n_per_epoxy

    epoxy_C_indices: list[int] = []
    ring_O_indices: list[int] = []
    for i in range(n_epoxies):
        base = epoxy_offset + i * n_per_epoxy
        for c_idx, o_idx in local_epoxide_pairs:
            epoxy_C_indices.append(base + c_idx)
            ring_O_indices.append(base + o_idx)

    amine_N_indices: list[int] = []
    amine_H_indices: list[int] = []
    amine_h_map: dict[int, list[int]] = {}
    for j in range(n_amines):
        base = amine_offset + j * n_per_amine
        for n_idx, h_indices in local_amine_nh:
            global_n = base + n_idx
            global_hs = [base + h for h in h_indices]
            amine_N_indices.append(global_n)
            amine_H_indices.extend(global_hs)
            amine_h_map[global_n] = global_hs

    template = ReactionTemplate(
        name='epoxy_amine_ring_opening',
        groups=['amine_N', 'epoxy_C', 'amine_H', 'ring_O'],
        pairs=[
            PairSpec(
                group_a='amine_N', group_b='epoxy_C',
                is_formation=True, r_min=3.0, r_max=6.0,
            ),
            PairSpec(
                group_a='amine_N', group_b='amine_H',
                is_formation=False, r_min=0.0, r_max=3.0,
            ),
            PairSpec(
                group_a='epoxy_C', group_b='ring_O',
                is_formation=False, r_min=0.0, r_max=3.0,
            ),
            PairSpec(
                group_a='amine_H', group_b='ring_O',
                is_formation=True, r_min=0.0, r_max=100.0,
                score_pair=False, count_as_reaction=False,
            ),
        ],
    )
    groups = {
        'amine_N': ReactiveGroup('amine_N', amine_N_indices),
        'epoxy_C': ReactiveGroup('epoxy_C', epoxy_C_indices),
        'amine_H': ReactiveGroup('amine_H', amine_H_indices),
        'ring_O':  ReactiveGroup('ring_O',  ring_O_indices),
    }

    return positions, species, template, groups, amine_h_map


# ── AIBN decomposition (Activation) ────────────────────────────────────────
# _AIBN_SMILES is defined in the top SMILES-constants block (needed earlier by
# full_aibn_initial_bonds' default arg).


def _find_aibn_azo_bonds(smiles: str = _AIBN_SMILES) -> list[tuple[int, int]]:
    """Return [(C_idx, N_idx), ...] for azo C-N single bonds in AIBN.

    Excludes nitrile C#N bonds. AIBN has two: C1-N5 and C7-N6.
    """
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    pairs = []
    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
        s1, s2 = a1.GetSymbol(), a2.GetSymbol()
        if bond.GetBondTypeAsDouble() != 1.0:
            continue
        if not (('C' in (s1, s2)) and ('N' in (s1, s2))):
            continue
        c_atom = a1 if s1 == 'C' else a2
        is_nitrile = any(b.GetBondTypeAsDouble() == 3.0 for b in c_atom.GetBonds())
        if not is_nitrile:
            c_idx = c_atom.GetIdx()
            n_idx = (a2 if s1 == 'C' else a1).GetIdx()
            pairs.append((c_idx, n_idx))
    return pairs


def _find_all_radical_centers(smiles: str) -> list[int]:
    """Return all C atoms bonded to exactly 3 C neighbors in the molecule.

    In AIBN (CC(C)(C#N)N=NC(C)(C)C#N), both central carbons (C1, C7) have
    3 C neighbors (2 methyls + nitrile C) and 1 N neighbor.
    """
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    result = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'C':
            continue
        n_c = sum(1 for n in atom.GetNeighbors() if n.GetSymbol() == 'C')
        if n_c == 3:
            result.append(atom.GetIdx())
    return result


def build_full_aibn_system(
    n_monomers: int,
    n_aibn: int,
    box_size: float,
    rng: np.random.Generator,
    monomer_smiles: str = _MONOMER_SMILES,
    aibn_smiles: str = _AIBN_SMILES,
    min_sep: float = 2.5,
    rdkit_seed: int = 42,
) -> tuple[np.ndarray, list[str],
           list[tuple[int, int]],
           ReactionTemplate, dict[str, ReactiveGroup],
           dict[int, int], dict[int, int]]:
    """Build a system with full AIBN molecules + vinyl monomers.

    Unlike build_vinyl_aibn_system (pre-formed radicals), this uses intact AIBN
    molecules that will be decomposed via V^d activation. Each AIBN contributes
    2 radical centers after decomposition.

    Returns
    -------
    positions       : (N, 3) array
    species         : element symbols
    aibn_azo_bonds  : [(C_global, N_global), ...] for activation V^d targets
    template        : propagation ReactionTemplate (4-group, for post-activation use)
    groups          : propagation groups (radical_C initially empty, filled after activation)
    propagation_map : {alpha_C → beta_C} for each monomer
    chain_c_map     : {radical_C → chain_C} (filled after activation)
    """
    aibn_pos, aibn_sp = _rdkit_3d(aibn_smiles, seed=rdkit_seed)
    mono_pos, mono_sp = _rdkit_3d(monomer_smiles, seed=rdkit_seed + 1)

    local_radical_centers = _find_all_radical_centers(aibn_smiles)
    local_azo_bonds = _find_aibn_azo_bonds(aibn_smiles)
    local_alpha, local_beta = _find_vinyl_alpha_beta(monomer_smiles)

    n_per_aibn = len(aibn_sp)
    n_per_mono = len(mono_sp)

    fragments = (
        [(aibn_pos, aibn_sp)] * n_aibn
        + [(mono_pos, mono_sp)] * n_monomers
    )
    positions, species = _place_fragments_in_box(
        fragments, box_size, rng, min_sep=min_sep,
    )

    aibn_offset = 0
    mono_offset = n_aibn * n_per_aibn

    aibn_azo_bonds_global: list[tuple[int, int]] = []
    for i in range(n_aibn):
        base = aibn_offset + i * n_per_aibn
        for c_local, n_local in local_azo_bonds:
            aibn_azo_bonds_global.append((base + c_local, base + n_local))

    radical_C_indices: list[int] = []
    chain_C_indices: list[int] = []
    for i in range(n_aibn):
        base = aibn_offset + i * n_per_aibn
        for rc_local in local_radical_centers:
            radical_C_indices.append(base + rc_local)
            chain_c_local = _find_chain_c_neighbor(aibn_smiles, rc_local)
            chain_C_indices.append(base + chain_c_local)

    alpha_C_indices = [
        mono_offset + j * n_per_mono + local_alpha
        for j in range(n_monomers)
    ]
    beta_C_indices = [
        mono_offset + j * n_per_mono + local_beta
        for j in range(n_monomers)
    ]
    propagation_map: dict[int, int] = {
        mono_offset + j * n_per_mono + local_alpha:
        mono_offset + j * n_per_mono + local_beta
        for j in range(n_monomers)
    }
    chain_c_map: dict[int, int] = {
        radical_C_indices[k]: chain_C_indices[k]
        for k in range(len(radical_C_indices))
    }

    template = ReactionTemplate(
        name='radical_vinyl_polymerization',
        groups=['radical_C', 'vinyl_alpha_C', 'chain_C', 'vinyl_beta_C'],
        pairs=[
            PairSpec('radical_C', 'vinyl_alpha_C', is_formation=True,
                     r_min=3.0, r_max=6.0),
            PairSpec('radical_C', 'chain_C', is_formation=True,
                     r_min=0.0, r_max=3.0, constraint_only=True),
            PairSpec('vinyl_alpha_C', 'vinyl_beta_C', is_formation=True,
                     r_min=0.0, r_max=3.0, constraint_only=True),
        ],
    )
    groups = {
        'radical_C':     ReactiveGroup('radical_C',     radical_C_indices),
        'vinyl_alpha_C': ReactiveGroup('vinyl_alpha_C', alpha_C_indices),
        'chain_C':       ReactiveGroup('chain_C',       chain_C_indices),
        'vinyl_beta_C':  ReactiveGroup('vinyl_beta_C',  beta_C_indices),
    }

    return (positions, species, aibn_azo_bonds_global,
            template, groups, propagation_map, chain_c_map)


def prune_undissociated_centers(
    groups: dict[str, ReactiveGroup],
    chain_c_map: dict[int, int],
    dissociated_c_indices: Iterable[int],
) -> tuple[dict[str, ReactiveGroup], dict[int, int], int]:
    """Drop radical centres that did not actually dissociate during activation.

    ``build_full_aibn_system`` registers *every* AIBN radical centre in
    ``radical_C``/``chain_C``/``chain_c_map`` unconditionally (by design — the
    dissociation result is unknown at build time).  When the V^d activation phase
    fails to cleave all azo C-N bonds, the still-intact centres are 4-coordinate
    (3 C + 1 azo-N, no H): they cannot add a monomer without over-coordinating,
    so the Layer-2 valence guard drops them *every cycle* (wasting candidate
    slots) and ``production_spin = n_radicals + 1`` over-counts the multiplicity.
    Un-cleaved AIBN has not formed a radical and is chemically inert, so those
    centres are removed here to reconcile the reactive groups with the actual
    activation outcome (specs/decisions.md 2026-07-06).

    Pure function: the input ``groups``/``chain_c_map`` are not mutated; new
    objects are returned.  ``radical_C``/``chain_C`` retain only centres whose
    global index is in ``dissociated_c_indices``.

    Parameters
    ----------
    groups
        Reactive groups; must contain ``'radical_C'`` (and usually ``'chain_C'``).
    chain_c_map
        ``{radical_C_global_idx: chain_C_global_idx}`` for each radical centre.
    dissociated_c_indices
        Global indices of the azo carbons whose C-N bond actually dissociated.

    Returns
    -------
    (new_groups, new_chain_c_map, n_pruned)
        Reconciled group dict, reconciled chain_c_map, and the number of radical
        centres removed.
    """
    dissociated = {int(c) for c in dissociated_c_indices}

    new_groups = {
        label: ReactiveGroup(g.label, list(g.atom_indices))
        for label, g in groups.items()
    }
    new_chain_c_map = {
        rc: cc for rc, cc in chain_c_map.items() if rc in dissociated
    }

    radical_group = new_groups.get('radical_C')
    if radical_group is None:
        return new_groups, new_chain_c_map, 0

    pruned_radicals = [
        rc for rc in radical_group.atom_indices if rc not in dissociated
    ]
    radical_group.atom_indices = [
        rc for rc in radical_group.atom_indices if rc in dissociated
    ]

    chain_group = new_groups.get('chain_C')
    if chain_group is not None:
        pruned_chain_cs = {
            chain_c_map[rc] for rc in pruned_radicals if rc in chain_c_map
        }
        chain_group.atom_indices = [
            cc for cc in chain_group.atom_indices if cc not in pruned_chain_cs
        ]

    return new_groups, new_chain_c_map, len(pruned_radicals)


def production_spin(n_radicals: int, cap: int | None = None) -> int:
    """High-spin total multiplicity (2S+1) for ``n_radicals`` unpaired electrons.

    Each surviving radical centre contributes one unpaired electron, so the
    high-spin coupling gives multiplicity ``n_radicals + 1`` (paper anchor:
    Table S1 activation → open-shell radicals; specs/decisions.md 2026-07-02
    Layer 3).  ``cap`` (from ``--production-spin-cap``) optionally clamps the
    value for diagnostics.
    """
    spin = n_radicals + 1
    if cap is not None:
        spin = min(spin, cap)
    return spin


def build_activation_template(
    aibn_azo_bonds: list[tuple[int, int]],
) -> tuple[ReactionTemplate, dict[str, ReactiveGroup]]:
    """Build activation (AIBN decomposition) template for V^d on C-N bonds.

    Paper anchor: Table S1 — Activation row, V^d applied to azo C-N bonds.
    Each AIBN molecule has 2 C-N bonds; both are targeted for dissociation.

    Returns (template, groups) for the activation phase only.
    """
    azo_C_indices = [c for c, _ in aibn_azo_bonds]
    azo_N_indices = [n for _, n in aibn_azo_bonds]

    template = ReactionTemplate(
        name='aibn_activation',
        groups=['azo_C', 'azo_N'],
        pairs=[
            PairSpec(
                group_a='azo_C',
                group_b='azo_N',
                is_formation=False,
                r_min=0.0,
                r_max=3.0,
            ),
        ],
    )
    groups = {
        'azo_C': ReactiveGroup('azo_C', azo_C_indices),
        'azo_N': ReactiveGroup('azo_N', azo_N_indices),
    }
    return template, groups
