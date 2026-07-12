# Test fixture: Provenzano 2025 `xl_trend.txt`

`xl_trend.txt` is redistributed unmodified from the Provenzano 2025
crosslinker dataset as a small unit-test fixture for
`scripts/compare_epoxy_external.py` (E2 external comparison,
specs/decisions.md 2026-07-12).

## Attribution (CC-BY-4.0)

- G. Provenzano et al., *ACS Applied Polymer Materials* **2025**, 7(8), 4876.
  DOI: 10.1021/acsapm.4c04208
- Dataset: Zenodo, DOI [10.5281/zenodo.11402476](https://doi.org/10.5281/zenodo.11402476)
  (record https://zenodo.org/records/15418263, file `xlinker.tgz`,
  member `xlinker/xl_trend.txt`)
- License: Creative Commons Attribution 4.0 International (CC-BY-4.0)

The file records the crosslinking-degree progression of the published 45%
example run: 3 whitespace-separated columns (cutoff radius `Radi` in
Angstrom, iteration `Iter` at that radius, crosslink degree `%` = reacted
amine H / initial reactive sites x 100), 50 data rows after a 3-line header.

The full dataset is NOT committed to this repository; fetch it with
`python scripts/fetch_provenzano2025.py` (see
`specs/dependency-license-matrix.md` for the license record).
