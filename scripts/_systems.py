"""Shared system builders for polymerization scripts."""
from __future__ import annotations

import numpy as np

from src.reactive.groups import PairSpec, ReactiveGroup, ReactionTemplate

# ── SMILES constants ─────────────────────────────────────────────────────────

_MONOMER_SMILES = 'C=CC(=O)OC'   # methyl acrylate
_INITIATOR_SMILES = 'CC(C)C#N'   # isobutyronitrile (closed-shell IBN radical model)
_DIAMINE_SMILES = 'NCCCCCCN'     # hexamethylenediamine
_DIACID_SMILES = 'OC(=O)CCCCC(=O)O'  # adipic acid

_AVOGADRO = 6.02214076e23  # mol^-1


def box_from_density(
    counts: dict[str, int],
    density_g_per_ml: float = 0.5,
) -> float:
    """Cubic box edge length (Å) for given molecule counts at a target density.

    Paper anchor: Supporting Information S-3..S-4 — vinyl and nylon initial
    configurations are generated at 0.5 g/mL. The formation bias is governed by
    near-contact events (PDF p.7, S-7), so reproducing the paper density is a
    prerequisite for confirmed bond formation. See specs/decisions.md
    "2026-06-13: T-G1a root-cause".
    """
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    total_mw = 0.0  # g/mol
    for smiles, n in counts.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f'Invalid SMILES for density calc: {smiles!r}')
        total_mw += Descriptors.MolWt(mol) * n
    mass_g = total_mw / _AVOGADRO
    vol_cm3 = mass_g / density_g_per_ml
    vol_a3 = vol_cm3 * 1.0e24  # 1 cm^3 = 1e24 A^3
    return float(vol_a3 ** (1.0 / 3.0))


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

def _rdkit_mol(smiles: str, seed: int = 42):
    """SMILES → 3D-embedded, MMFF-optimized RDKit Mol with explicit hydrogens.

    Single source of atom ordering for both the system builder (`_rdkit_3d`) and
    the classical structure-prep stage. The OpenFF Topology used by prep is built
    from this exact Mol (via Molecule.from_rdkit) so that AddHs atom ordering
    matches the builder's `species`/`positions` order; otherwise the global
    `groups` / `propagation_map` indices would be invalidated when prep hands
    relaxed coordinates back. See specs/decisions.md 2026-06-14 decision D-4.

    The conformer seed only affects 3-D coordinates, not atom indexing, so the
    returned connectivity/ordering is deterministic for a given SMILES.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as exc:
        raise ImportError(
            'RDKit is required for vinyl/AIBN system builder. '
            'Install with: pip install "pfpoly[rdkit]"'
        ) from exc

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f'RDKit could not parse SMILES: {smiles!r}')
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    result = AllChem.EmbedMolecule(mol, params)
    if result == -1:
        raise RuntimeError(f'EmbedMolecule failed for {smiles!r}')
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


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


def _find_vinyl_alpha_beta(smiles: str) -> tuple[int, int]:
    """Return (alpha_idx, beta_idx) local atom indices for the CH2=CH- motif.

    alpha-C: =CH2 end (2 H, vinyl terminus) → Group J
    beta-C:  =CH- inner  (1 H) → propagation target
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
        if h_a == 2 and h_b == 1:
            return a.GetIdx(), b.GetIdx()
        if h_b == 2 and h_a == 1:
            return b.GetIdx(), a.GetIdx()
    raise ValueError(f'No CH2=CH- pattern found in {smiles!r}')


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
            offset = np.clip(base + jitter, 2.0, box_size - 2.0)
            candidate = centered @ R.T + offset

            ok = True
            for prev in placed:
                diffs = candidate[:, np.newaxis, :] - prev[np.newaxis, :, :]
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
                score_pair=False,
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


# ── AIBN decomposition (Activation) ────────────────────────────────────────

_AIBN_SMILES = 'CC(C)(C#N)N=NC(C)(C)C#N'


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
