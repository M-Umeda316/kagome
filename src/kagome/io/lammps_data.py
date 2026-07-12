"""Minimal read-only parser for LAMMPS data files (external comparison, E2).

Purpose: read the topology (element per atom + bond list) out of the
Provenzano 2025 crosslinker structures (``data.relaxed00`` /
``data.xlinker``, Zenodo DOI 10.5281/zenodo.11402476, CC-BY-4.0) so that
:mod:`kagome.analysis.network` can compute the same species / gel metrics on
their classical-MD structures as on our TDBB runs (specs/decisions.md
2026-07-12 Track 2 / E2 design).

Scope is deliberately narrow — this is NOT a general LAMMPS reader:

- ``Masses``  : atom type -> mass (amu); mapped to element symbols by nearest
                standard atomic mass among {H, C, N, O} within a tolerance.
- ``Atoms``   : only the ``full`` style is supported
                (``id mol type q x y z [ix iy iz]``); coordinates and charges
                are ignored, only the id -> type mapping is kept.
- ``Bonds``   : ``id type atom1 atom2``; returned as 0-indexed pairs
                (LAMMPS atom ids are 1-indexed).
- everything else (Velocities, Angles, Dihedrals, Impropers, all Coeffs
  sections, ...) is skipped.

Units: masses in amu (LAMMPS ``units real``); box bounds in the file's
length unit (Angstrom for ``units real``) — both passed through unchanged.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)

# Standard atomic masses (amu) for the only elements present in the
# epoxy-amine (DGEBA + DETA) all-atom system.  Class II force fields carry
# slightly different values per atom type (e.g. 12.01115 for C), hence the
# nearest-mass match with a tolerance instead of exact lookup.
_ELEMENT_MASSES: dict[str, float] = {
    'H': 1.008,
    'C': 12.011,
    'N': 14.007,
    'O': 15.999,
}

#: Maximum |mass - standard mass| (amu) accepted when inferring an element.
MASS_TOLERANCE_AMU: float = 0.5


@dataclass
class LammpsData:
    """Topology-only view of a LAMMPS data file.

    - ``species``: element symbol per atom, indexed by 0-based atom index
      (LAMMPS atom id - 1) — the same convention as the ``species`` argument
      of :mod:`kagome.analysis.network`.
    - ``bonds``: ``(i, j)`` 0-indexed atom pairs.
    - ``n_atoms``: atom count from the header.
    - ``box_bounds``: ``((xlo, xhi), (ylo, yhi), (zlo, zhi))`` in the file's
      length unit (Angstrom for ``units real``).
    """

    species: list[str] = field(default_factory=list)
    bonds: list[tuple[int, int]] = field(default_factory=list)
    n_atoms: int = 0
    box_bounds: tuple[tuple[float, float], ...] = ()


def element_from_mass(mass: float) -> str:
    """Element symbol whose standard atomic mass is nearest to ``mass`` (amu).

    Only {H, C, N, O} are considered (the epoxy-amine system contains nothing
    else).  Raises ``ValueError`` if no standard mass lies within
    ``MASS_TOLERANCE_AMU`` of ``mass`` — an unknown element must fail loudly
    rather than corrupt the species list.
    """
    best_symbol, best_diff = '', float('inf')
    for symbol, standard in _ELEMENT_MASSES.items():
        diff = abs(mass - standard)
        if diff < best_diff:
            best_symbol, best_diff = symbol, diff
    if best_diff > MASS_TOLERANCE_AMU:
        raise ValueError(
            f'mass {mass} amu matches no element in '
            f'{sorted(_ELEMENT_MASSES)} within {MASS_TOLERANCE_AMU} amu'
        )
    return best_symbol


def _strip_comment(line: str) -> str:
    """Drop a trailing ``# ...`` comment (e.g. section style hints like
    ``Atoms # full``) and surrounding whitespace."""
    return line.split('#', 1)[0].strip()


def read_lammps_data(path: Path | str) -> LammpsData:
    """Parse a LAMMPS data file into a :class:`LammpsData` topology.

    Layout assumptions (matched against the actual Provenzano 2025 files):
    the first line is a free-text title; header lines (counts, box bounds)
    follow until the first section keyword; every section header line starts
    with a letter while data rows start with a number, so any
    letter-initial line switches the current section.  Atom ids may appear
    in any order (the real files are unordered).  Only ``Masses``, ``Atoms``
    (style ``full``) and ``Bonds`` are consumed; other sections are skipped.
    """
    path = Path(path)
    n_atoms = 0
    n_bonds_declared = 0
    bounds: dict[str, tuple[float, float]] = {}
    masses: dict[int, float] = {}
    type_by_atom: dict[int, int] = {}
    bonds: list[tuple[int, int]] = []

    section: str | None = None
    with open(path, 'r', encoding='utf-8') as fh:
        first_line = True
        for raw in fh:
            if first_line:  # free-text title line
                first_line = False
                continue
            line = _strip_comment(raw)
            if not line:
                continue

            if line[0].isalpha():  # section header (data rows start numeric)
                section = line
                continue

            tokens = line.split()
            if section is None:  # still in the header block
                if len(tokens) == 2 and tokens[1] == 'atoms':
                    n_atoms = int(tokens[0])
                elif len(tokens) == 2 and tokens[1] == 'bonds':
                    n_bonds_declared = int(tokens[0])
                elif len(tokens) == 4 and tokens[2].endswith('lo'):
                    # e.g. '-39.86 40.01 xlo xhi'
                    bounds[tokens[2][0]] = (float(tokens[0]), float(tokens[1]))
                continue

            if section == 'Masses':
                masses[int(tokens[0])] = float(tokens[1])
            elif section == 'Atoms':
                # full style: id mol type q x y z [ix iy iz]
                if len(tokens) < 7:
                    raise ValueError(
                        f'{path}: Atoms row has {len(tokens)} columns, '
                        f'expected >= 7 (style "full"): {raw.strip()!r}'
                    )
                type_by_atom[int(tokens[0])] = int(tokens[2])
            elif section == 'Bonds':
                # id type atom1 atom2 (atom ids 1-indexed -> 0-indexed)
                bonds.append((int(tokens[2]) - 1, int(tokens[3]) - 1))
            # other sections: skipped

    if len(type_by_atom) != n_atoms:
        raise ValueError(
            f'{path}: header declares {n_atoms} atoms but the Atoms section '
            f'has {len(type_by_atom)} rows'
        )
    if n_bonds_declared and len(bonds) != n_bonds_declared:
        raise ValueError(
            f'{path}: header declares {n_bonds_declared} bonds but the Bonds '
            f'section has {len(bonds)} rows'
        )

    element_by_type = {t: element_from_mass(m) for t, m in masses.items()}

    species: list[str] = []
    for atom_id in range(1, n_atoms + 1):
        if atom_id not in type_by_atom:
            raise ValueError(f'{path}: atom id {atom_id} missing from Atoms')
        atom_type = type_by_atom[atom_id]
        if atom_type not in element_by_type:
            raise ValueError(
                f'{path}: atom id {atom_id} has type {atom_type} with no '
                f'Masses entry'
            )
        species.append(element_by_type[atom_type])

    box = tuple(bounds.get(axis, (0.0, 0.0)) for axis in ('x', 'y', 'z'))
    logger.info(
        'read_lammps_data(%s): %d atoms, %d bonds, %d atom types',
        path, n_atoms, len(bonds), len(masses),
    )
    return LammpsData(
        species=species, bonds=bonds, n_atoms=n_atoms, box_bounds=box,
    )
