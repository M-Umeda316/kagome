"""Shared system builders for polymerization scripts."""
from __future__ import annotations

import numpy as np

from src.reactive.groups import PairSpec, ReactiveGroup, ReactionTemplate

# ── vinyl / AIBN system ──────────────────────────────────────────────────────

_MONOMER_SMILES = 'C=CC(=O)OC'   # methyl acrylate
_INITIATOR_SMILES = 'CC(C)C#N'   # isobutyronitrile (closed-shell IBN radical model)


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
                r_min=1.6,
                r_max=4.5,
            ),
        ],
    )

    groups = {
        'C_donor': ReactiveGroup('C_donor', group_a_indices),
        'C_acceptor': ReactiveGroup('C_acceptor', group_b_indices),
    }

    return template, groups


# ── RDKit-based helpers ──────────────────────────────────────────────────────

def _rdkit_3d(smiles: str, seed: int = 42) -> tuple[np.ndarray, list[str]]:
    """SMILES → (positions Å, species list) via RDKit EmbedMolecule + MMFF."""
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

    Applies a random 3-D rotation + random translation per molecule.
    Raises RuntimeError if any molecule cannot be placed without overlap.
    """
    placed: list[np.ndarray] = []
    all_species: list[str] = []

    for mol_idx, (template_pos, mol_species) in enumerate(fragments):
        centered = template_pos - template_pos.mean(axis=0)
        for _attempt in range(max_attempts):
            axis = rng.standard_normal(3)
            if np.linalg.norm(axis) < 1e-10:
                axis = np.array([0.0, 0.0, 1.0])
            R = _rotation_matrix(axis, rng.uniform(0.0, 2.0 * np.pi))
            offset = rng.uniform(2.0, box_size - 2.0, 3)
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
                f'Try increasing --box-size.'
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
) -> tuple[np.ndarray, list[str], ReactionTemplate, dict[str, ReactiveGroup], dict[int, int]]:
    """Build a vinyl + AIBN-initiator system for radical polymerization.

    Placement order: all initiators first, then all monomers.
    This keeps global index arithmetic straightforward.

    Returns
    -------
    positions      : (N, 3) Å array
    species        : list of element symbols
    template       : ReactionTemplate (radical_C + vinyl_alpha_C → bond)
    groups         : {'radical_C': Group I, 'vinyl_alpha_C': Group J}
    propagation_map: {alpha_C_global_idx: beta_C_global_idx} for each monomer

    Paper anchor: Fig. 1, Section 2 (radical chain-growth mechanism).
    Design: specs/decisions.md — "T-G1: vinyl radical polymerization system".
    """
    # ── generate single-molecule 3D templates ──
    init_pos, init_sp = _rdkit_3d(initiator_smiles, seed=rdkit_seed)
    mono_pos, mono_sp = _rdkit_3d(monomer_smiles, seed=rdkit_seed + 1)

    local_radical = _find_ibn_radical_c(initiator_smiles)
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
    alpha_C_indices = [
        mono_offset + j * n_per_mono + local_alpha
        for j in range(n_monomers)
    ]
    propagation_map: dict[int, int] = {
        mono_offset + j * n_per_mono + local_alpha:
        mono_offset + j * n_per_mono + local_beta
        for j in range(n_monomers)
    }

    # ── reaction template ──
    template = ReactionTemplate(
        name='radical_vinyl_polymerization',
        groups=['radical_C', 'vinyl_alpha_C'],
        pairs=[
            PairSpec(
                group_a='radical_C',
                group_b='vinyl_alpha_C',
                is_formation=True,
                r_min=1.6,
                r_max=4.5,
            ),
        ],
    )
    groups = {
        'radical_C':    ReactiveGroup('radical_C',    radical_C_indices),
        'vinyl_alpha_C': ReactiveGroup('vinyl_alpha_C', alpha_C_indices),
    }

    return positions, species, template, groups, propagation_map
