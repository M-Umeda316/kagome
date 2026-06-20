"""Shared chemistry helpers: RDKit molecule building and density-based box sizing."""
from __future__ import annotations

_AVOGADRO = 6.02214076e23  # mol^-1


def box_from_density(
    counts: dict[str, int],
    density_g_per_ml: float = 0.5,
) -> float:
    """Cubic box edge length (Å) for given molecule counts at a target density.

    ``counts`` is keyed by SMILES; because it is a dict, two entries with the
    same SMILES string are necessarily the same key and their counts collapse
    into one. This is intended — identical SMILES denote the same molecule, so
    the molar mass contribution is summed via the ``* n`` count. Callers that
    need to distinguish chemically identical but separately-tracked species must
    aggregate the count before calling (a single key with the total ``n``).

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
            'Install with: pip install "kagome[rdkit]"'
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
