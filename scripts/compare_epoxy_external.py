"""E2 external comparison: our TDBB epoxy-amine run vs Provenzano 2025 (ref#2).

Cross-method comparison against the classical-MD distance-based crosslinking
protocol of Provenzano et al., ACS Appl. Polym. Mater. 2025, 7(8), 4876
(dataset: Zenodo DOI 10.5281/zenodo.11402476, CC-BY-4.0; fetch with
``scripts/fetch_provenzano2025.py``).  Protocol: specs/decisions.md
2026-07-12 【Track 2 / E2 設計】.

Comparison axis is CONVERSION, never time — ref#2's progression axis
(cutoff radius x iteration) has no physical-time mapping, and our TDBB
kinetics are biased (f2/friction deviations), so only conversion-resolved
network structure is compared.  Two conversion bases are reported:

- epoxide basis : 1 - (closed epoxide rings) / (initial epoxide rings)
                  — the repo-wide basis (dimensionless, 0..1).
- amine-H basis : 1 - (current N-H bonds) / (initial N-H bonds)
                  — ref#2's own definition ("reacted amine H / initial
                  reactive sites"); equals the epoxide basis at 1:1
                  stoichiometry (their 800:320 DGEBA:DETA system).

Monomer-count-dependent metrics are normalized to fractions so a
45,600-atom system and a 1,770-atom system are comparable:

- largest component fraction : already a fraction of monomers.
- tertiary amine fraction    : tertiary N / all amine N.
- fully-ring-opened DGEBA fraction : fully reacted epoxy monomers /
                                     epoxy monomers.

The ref#2 gel-curve point is COMPUTED BY US from their published
``data.xlinker`` structure (the 45%-conversion example) — it is not a value
claimed in their paper; every figure carries that caveat.  No smoothing or
filtering is applied anywhere (repo rule).

Usage:
    python scripts/compare_epoxy_external.py \
        --run-dir runs/epoxy_amine_mid30 \
        --external-dir data/external/provenzano2025/xlinker
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from kagome.analysis.carothers import monomer_sets_from_bonds
from kagome.analysis.network import (
    amine_ranks,
    crosslink_counts,
    find_epoxide_rings,
    gel_point_flory_stockmayer,
    largest_component_fraction,
    load_topology_snapshots,
)
from kagome.io.lammps_data import read_lammps_data

logger = logging.getLogger('compare_epoxy_external')

ZENODO_DOI = '10.5281/zenodo.11402476'
REF2_CITATION = ('Provenzano et al., ACS Appl. Polym. Mater. 2025, '
                 '7(8), 4876')
REF2_CAVEAT = (f'ref#2 point computed by this repo from the published '
               f'structure (Zenodo {ZENODO_DOI}), not a value claimed by '
               f'{REF2_CITATION}.')


# ── Parsers ──────────────────────────────────────────────────────────────────

def parse_xl_trend(path: Path | str) -> list[tuple[float, int, float]]:
    """Parse ref#2's ``xl_trend.txt`` crosslink-degree progression.

    Format (verified against the Zenodo file): a ``#`` comment line, a column
    header line (``Radi    Iter    %``), a dashes separator, then 50
    whitespace-separated data rows.  Returns ``(radius_A, iteration,
    percent)`` tuples: cutoff radius in Angstrom, iteration number at that
    radius (1-based), and crosslink degree in percent (ref#2's amine-H
    basis: reacted amine H / initial reactive sites x 100).
    """
    rows: list[tuple[float, int, float]] = []
    with open(path, 'r', encoding='utf-8') as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith('#') or line.startswith('-'):
                continue
            tokens = line.split()
            if len(tokens) != 3:
                raise ValueError(f'{path}: expected 3 columns, got {line!r}')
            try:
                rows.append((float(tokens[0]), int(tokens[1]), float(tokens[2])))
            except ValueError:
                if tokens[0] == 'Radi':  # column header line
                    continue
                raise ValueError(f'{path}: unparseable row {line!r}') from None
    return rows


# ── Network metrics (shared between our runs and ref#2 structures) ──────────

def network_metrics(
    bonds_initial: list,
    bonds_final: list,
    species: list[str],
) -> dict:
    """Conversion + normalized network metrics from two bond topologies.

    ``bonds_initial`` defines the monomer sets, the initial closed-epoxide
    rings and the initial N-H bond count; ``bonds_final`` is the state to
    measure.  Bond items are ``(i, j)`` or ``(i, j, order)``.  All
    conversions and fractions are dimensionless (0..1).  Works identically
    for our ``topology.jsonl`` snapshots and for ref#2 LAMMPS structures
    parsed by :func:`kagome.io.lammps_data.read_lammps_data`.
    """
    monomer_sets = monomer_sets_from_bonds(bonds_initial)
    rings_initial = find_epoxide_rings(bonds_initial, species)
    n_epoxide_initial = len(rings_initial)
    n_epoxide_final = len(find_epoxide_rings(bonds_final, species))

    amine_h_initial = sum(amine_ranks(bonds_initial, species).values())
    ranks_final = amine_ranks(bonds_final, species)
    amine_h_final = sum(ranks_final.values())
    n_amine_n = len(ranks_final)

    counts = crosslink_counts(
        bonds_final, species, rings_initial, monomer_sets)

    # Epoxy monomers = monomers owning at least one initial epoxide ring
    # (the denominator for the fully-ring-opened DGEBA fraction).
    atom2mono: dict[int, int] = {}
    for m, atoms in enumerate(monomer_sets):
        for a in atoms:
            atom2mono[int(a)] = m
    epoxy_monomers = {atom2mono[o] for _c1, _c2, o in rings_initial
                      if o in atom2mono}
    n_epoxy_monomers = len(epoxy_monomers)

    return {
        'n_monomers': len(monomer_sets),
        'n_epoxy_monomers': n_epoxy_monomers,
        'n_amine_n': n_amine_n,
        'n_epoxide_initial': n_epoxide_initial,
        'n_epoxide_final': n_epoxide_final,
        'epoxide_conversion': (
            1.0 - n_epoxide_final / n_epoxide_initial
            if n_epoxide_initial else 0.0),
        'amine_h_initial': amine_h_initial,
        'amine_h_final': amine_h_final,
        'amine_h_conversion': (
            1.0 - amine_h_final / amine_h_initial
            if amine_h_initial else 0.0),
        'largest_component_fraction': largest_component_fraction(
            bonds_final, monomer_sets),
        'n_tertiary_amine': counts['n_tertiary_amine'],
        'tertiary_amine_fraction': (
            counts['n_tertiary_amine'] / n_amine_n if n_amine_n else 0.0),
        'n_fully_reacted_epoxy_monomers':
            counts['n_fully_reacted_epoxy_monomers'],
        'fully_reacted_epoxy_fraction': (
            counts['n_fully_reacted_epoxy_monomers'] / n_epoxy_monomers
            if n_epoxy_monomers else 0.0),
    }


def run_series(
    snapshots: list[tuple[int, int, list]],
    species: list[str],
) -> list[dict]:
    """Per-topology-snapshot metrics for our TDBB run.

    ``snapshots`` from :func:`kagome.analysis.network.load_topology_snapshots`
    (first record = initial topology, cycle -1).  Each row: step, cycle,
    epoxide/amine-H conversion (bases documented in the module docstring),
    largest component fraction and normalized crosslink fractions — all
    measured against the INITIAL snapshot's rings/monomers.
    """
    if not snapshots:
        return []
    bonds_initial = snapshots[0][2]
    rows: list[dict] = []
    for step, cycle, bonds in snapshots:
        metrics = network_metrics(bonds_initial, bonds, species)
        rows.append({
            'step': int(step),
            'cycle': int(cycle),
            'epoxide_conversion': metrics['epoxide_conversion'],
            'amine_h_conversion': metrics['amine_h_conversion'],
            'largest_component_fraction':
                metrics['largest_component_fraction'],
            'tertiary_amine_fraction': metrics['tertiary_amine_fraction'],
            'fully_reacted_epoxy_fraction':
                metrics['fully_reacted_epoxy_fraction'],
        })
    return rows


def external_metrics(initial_path: Path, final_path: Path) -> dict:
    """Ref#2 network metrics from their LAMMPS data files.

    ``initial_path`` = ``data.relaxed00`` (uncrosslinked), ``final_path`` =
    ``data.xlinker`` (the published 45%-conversion example).  Both parsed
    with :func:`kagome.io.lammps_data.read_lammps_data`; the species lists
    must agree atom-by-atom (their protocol changes atom TYPES on reaction
    but never elements) — a mismatch raises ``ValueError`` because every
    ring/monomer mapping would be invalid.
    """
    initial = read_lammps_data(initial_path)
    final = read_lammps_data(final_path)
    if initial.species != final.species:
        raise ValueError(
            f'element lists differ between {initial_path} and {final_path}; '
            f'cannot map rings/monomers across the two structures')
    metrics = network_metrics(initial.bonds, final.bonds, initial.species)
    metrics['n_atoms'] = initial.n_atoms
    return metrics


# ── Species resolution (same precedence as reproduce_figures) ────────────────

def _resolve_run_species(
    run_dir: Path,
    species_json: Path | None,
) -> list[str] | None:
    """--species-json > summary.json rebuild (RDKit) > trajectory header.

    Delegates to ``scripts.reproduce_figures._resolve_species`` so the
    precedence stays identical across analysis scripts.
    """
    try:
        from scripts.reproduce_figures import _resolve_species
    except ImportError:
        from reproduce_figures import _resolve_species  # script-dir fallback
    return _resolve_species(
        species_json,
        run_dir / 'summary.json',
        run_dir / 'trajectory.jsonl',
    )


# ── Figures ──────────────────────────────────────────────────────────────────

def plot_gel_curve_comparison(
    series: list[dict],
    theirs: dict | None,
    alpha_gel: float,
    output_dir: Path,
) -> None:
    """Gel curve (largest component fraction vs epoxide conversion), ours vs
    the single ref#2 point computed from ``data.xlinker``.

    x axis: epoxide conversion (1 - closed rings / initial rings) — stated on
    the label because ref#2's own reporting uses the amine-H basis (identical
    only at 1:1 stoichiometry).  No smoothing.  Saves
    ``gel_curve_comparison.png/.pdf`` (dpi 150).
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    output_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.5, 5.5))

    ax.plot(
        [row['epoxide_conversion'] for row in series],
        [row['largest_component_fraction'] for row in series],
        'o-', color='tab:blue', markersize=6, linewidth=1.0,
        label='ours: TDBB run (largest component)',
    )
    ax.axvline(alpha_gel, color='k', linestyle='--', linewidth=1.0,
               label=f'Flory-Stockmayer α_gel = {alpha_gel:.3f}')
    if theirs is not None:
        ax.plot(
            theirs['epoxide_conversion'],
            theirs['largest_component_fraction'],
            marker='*', color='tab:red', markersize=16, linestyle='none',
            label='ref#2 data.xlinker (computed by us)',
        )

    ax.set_xlabel('Epoxide conversion α (1 − closed epoxide rings / initial rings)')
    ax.set_ylabel('Largest monomer-component fraction')
    ax.set_title('E2 gel curve: TDBB (ours) vs Provenzano 2025 (ref#2)')
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.05)
    ax.legend(loc='upper left', fontsize=8)
    fig.text(
        0.01, 0.01,
        'Note: ref#2 point computed by this repo from the published '
        f'structure (Zenodo {ZENODO_DOI});\nnot a value claimed by '
        'Provenzano et al., ACS Appl. Polym. Mater. 2025, 7(8), 4876.',
        fontsize=6.5, color='dimgray',
    )
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    for fmt in ('png', 'pdf'):
        fig.savefig(output_dir / f'gel_curve_comparison.{fmt}', dpi=150)
    plt.close(fig)
    print(f'Saved gel_curve_comparison.png/.pdf to {output_dir}')


def plot_progression_comparison(
    series: list[dict],
    xl_trend: list[tuple[float, int, float]] | None,
    output_dir: Path,
) -> None:
    """Side-by-side progression panels — EXPLICITLY non-comparable axes.

    Left: our epoxide conversion vs TDBB cycle (cycle -1 = initial
    topology).  Right: ref#2 crosslink degree (%, their amine-H basis) vs
    cumulative crosslinking iteration from ``xl_trend.txt``, with the
    cutoff-radius schedule annotated.  The panels share no axis: TDBB cycles
    and cutoff-radius iterations have no mutual (or physical-time) mapping,
    so only the progression SHAPE may be eyeballed.  No smoothing.  Saves
    ``progression_comparison.png/.pdf`` (dpi 150).
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    output_dir.mkdir(parents=True, exist_ok=True)
    fig, (ax_ours, ax_ref) = plt.subplots(1, 2, figsize=(11, 4.8))

    ax_ours.plot(
        [row['cycle'] for row in series],
        [row['epoxide_conversion'] for row in series],
        'o-', color='tab:blue', markersize=5, linewidth=1.0,
    )
    ax_ours.set_xlabel('TDBB cycle (cycle −1 = initial topology)')
    ax_ours.set_ylabel('Epoxide conversion α')
    ax_ours.set_title('ours: TDBB progression')
    ax_ours.set_ylim(bottom=0.0)

    if xl_trend:
        cumulative = list(range(1, len(xl_trend) + 1))
        ax_ref.plot(
            cumulative, [pct for _r, _i, pct in xl_trend],
            's-', color='tab:red', markersize=4, linewidth=1.0,
        )
        # Annotate where the cutoff radius steps up (their schedule axis).
        for pos, (radius, iteration, _pct) in zip(cumulative, xl_trend):
            if iteration == 1:
                ax_ref.axvline(pos, color='gray', linestyle=':',
                               linewidth=0.6)
                ax_ref.text(pos, ax_ref.get_ylim()[1] * 0.02,
                            f'{radius:g} Å', fontsize=6.5, color='gray',
                            rotation=90, va='bottom')
        ax_ref.set_xlabel('Cumulative iteration (cutoff radius × 10 iterations)')
        ax_ref.set_ylabel('Crosslink degree (%) — reacted amine H / initial sites')
        ax_ref.set_title('ref#2: distance-based crosslinking (xl_trend.txt)')
        ax_ref.set_ylim(bottom=0.0)
    else:
        ax_ref.text(0.5, 0.5, 'external xl_trend.txt unavailable',
                    ha='center', va='center', transform=ax_ref.transAxes)
        ax_ref.set_axis_off()

    fig.suptitle(
        'E2 progression shapes — axes NOT comparable (TDBB cycles vs '
        'cutoff-radius iterations); shape only',
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    for fmt in ('png', 'pdf'):
        fig.savefig(output_dir / f'progression_comparison.{fmt}', dpi=150)
    plt.close(fig)
    print(f'Saved progression_comparison.png/.pdf to {output_dir}')


# ── Report ───────────────────────────────────────────────────────────────────

_TABLE_ROWS: list[tuple[str, str]] = [
    ('epoxide_conversion', 'Epoxide conversion (ring basis)'),
    ('amine_h_conversion', 'Amine-H conversion (ref#2 basis)'),
    ('largest_component_fraction', 'Largest component fraction'),
    ('tertiary_amine_fraction', 'Tertiary amine fraction (3° N / all N)'),
    ('fully_reacted_epoxy_fraction', 'Fully-ring-opened DGEBA fraction'),
    ('n_monomers', 'Monomers'),
    ('n_epoxide_initial', 'Initial epoxide rings'),
    ('amine_h_initial', 'Initial amine N-H'),
]


def _fmt(value: float | int | None) -> str:
    if value is None:
        return 'n/a'
    if isinstance(value, int):
        return str(value)
    return f'{value:.4f}'


def write_report(
    ours_final: dict,
    theirs: dict | None,
    series: list[dict],
    xl_trend: list[tuple[float, int, float]] | None,
    alpha_gel: float,
    run_dir: Path,
    external_dir: Path | None,
    output_dir: Path,
) -> None:
    """Write ``e2_comparison.json`` (full numbers) and ``e2_comparison.md``
    (matched-metric table).  Fractions are dimensionless (0..1); the two
    conversion bases are labeled per row."""
    output_dir.mkdir(parents=True, exist_ok=True)

    payload = {
        'protocol': 'specs/decisions.md 2026-07-12 Track 2 / E2',
        'reference': REF2_CITATION,
        'reference_dataset': f'Zenodo DOI {ZENODO_DOI} (CC-BY-4.0)',
        'caveat': REF2_CAVEAT,
        'run_dir': str(run_dir),
        'external_dir': str(external_dir) if external_dir else None,
        'flory_stockmayer_alpha_gel': alpha_gel,
        'ours_final': ours_final,
        'ours_series': series,
        'theirs_final': theirs,
        'theirs_xl_trend_rows': len(xl_trend) if xl_trend else 0,
        'theirs_xl_trend_final_percent':
            xl_trend[-1][2] if xl_trend else None,
    }
    json_path = output_dir / 'e2_comparison.json'
    json_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + '\n',
        encoding='utf-8')
    print(f'Saved {json_path}')

    lines = [
        '# E2 external comparison: TDBB (ours) vs Provenzano 2025 (ref#2)',
        '',
        f'- ours: `{run_dir}` (final topology snapshot)',
        f'- ref#2: `{external_dir}` — `data.relaxed00` (initial) vs '
        f'`data.xlinker` (published 45% example)' if external_dir
        else '- ref#2: external data unavailable (metrics skipped)',
        f'- reference: {REF2_CITATION}; dataset Zenodo DOI {ZENODO_DOI} '
        f'(CC-BY-4.0)',
        f'- Flory-Stockmayer gel point (f=2, g=5, r=1): '
        f'α_gel = {alpha_gel:.3f} (epoxide basis)',
        f'- **Caveat**: {REF2_CAVEAT}',
        '- Conversion bases: epoxide = 1 − closed rings / initial rings; '
        'amine-H = 1 − current N−H / initial N−H (ref#2 definition). '
        'Identical only at 1:1 stoichiometry.',
        '- All fractions dimensionless (0..1); no smoothing/filtering '
        'applied.',
        '',
        '| Metric | ours (TDBB) | ref#2 (computed by us) |',
        '|---|---|---|',
    ]
    for key, label in _TABLE_ROWS:
        lines.append(
            f'| {label} | {_fmt(ours_final.get(key))} | '
            f'{_fmt(theirs.get(key)) if theirs else "n/a"} |'
        )
    if xl_trend:
        lines += [
            '',
            f'ref#2 `xl_trend.txt`: {len(xl_trend)} rows, final crosslink '
            f'degree {xl_trend[-1][2]:g}% (their amine-H basis) — '
            f'progression axis is cutoff radius × iteration, not time.',
        ]
    md_path = output_dir / 'e2_comparison.md'
    md_path.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print(f'Saved {md_path}')


# ── CLI ──────────────────────────────────────────────────────────────────────

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--run-dir', type=Path, required=True,
                        help='our TDBB run directory (topology.jsonl, ...)')
    parser.add_argument('--external-dir', type=Path, default=None,
                        help='extracted ref#2 xlinker directory '
                             '(data.relaxed00, data.xlinker, xl_trend.txt); '
                             'external overlays skipped if omitted/missing')
    parser.add_argument('--output-dir', type=Path, default=None,
                        help='default: <run-dir>/figures/e2_comparison')
    parser.add_argument('--species-json', type=Path, default=None,
                        help='JSON list of element symbols overriding '
                             'summary/trajectory species resolution')
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')

    run_dir: Path = args.run_dir
    output_dir: Path = (args.output_dir
                        or run_dir / 'figures' / 'e2_comparison')

    snapshots = load_topology_snapshots(run_dir / 'topology.jsonl')
    if not snapshots:
        print(f'ERROR: no topology snapshots in {run_dir / "topology.jsonl"}')
        return 1
    species = _resolve_run_species(run_dir, args.species_json)
    if species is None:
        print('ERROR: could not resolve species for the run '
              '(--species-json / summary rebuild / trajectory header)')
        return 1

    series = run_series(snapshots, species)
    ours_final = network_metrics(snapshots[0][2], snapshots[-1][2], species)

    theirs: dict | None = None
    xl_trend: list[tuple[float, int, float]] | None = None
    external_dir: Path | None = args.external_dir
    if external_dir is not None:
        initial_path = external_dir / 'data.relaxed00'
        final_path = external_dir / 'data.xlinker'
        if initial_path.exists() and final_path.exists():
            theirs = external_metrics(initial_path, final_path)
        else:
            logger.warning(
                'external structures missing under %s (need data.relaxed00 '
                'and data.xlinker) -- ref#2 metrics skipped', external_dir)
        trend_path = external_dir / 'xl_trend.txt'
        if trend_path.exists():
            xl_trend = parse_xl_trend(trend_path)
        else:
            logger.warning('%s missing -- ref#2 progression panel skipped',
                           trend_path)
    else:
        logger.warning('--external-dir not given -- external overlays '
                       'skipped (run scripts/fetch_provenzano2025.py first)')

    alpha_gel = gel_point_flory_stockmayer(f=2, g=5, r=1.0)

    plot_gel_curve_comparison(series, theirs, alpha_gel, output_dir)
    plot_progression_comparison(series, xl_trend, output_dir)
    write_report(ours_final, theirs, series, xl_trend, alpha_gel,
                 run_dir, external_dir, output_dir)
    return 0


if __name__ == '__main__':
    sys.exit(main())
