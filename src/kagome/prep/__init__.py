"""Classical structure-preparation package.

Decouples initial-structure preparation (packing / densification /
thermalization with a classical force field via OpenMM + OpenFF) from the ML
TDBB production run. The two stages may run in different conda environments, so
the relaxed structure is handed over through a small JSON file
(:class:`kagome.prep.structure_io.PreparedStructure`).

See specs/decisions.md 2026-06-14 "Decouple initial-structure preparation".
"""
