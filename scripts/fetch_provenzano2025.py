"""Fetch the Provenzano 2025 crosslinker dataset (E2 external comparison).

Downloads ``xlinker.tgz`` from Zenodo (DOI 10.5281/zenodo.11402476, record
15418263, license CC-BY-4.0) into ``data/external/provenzano2025/`` and
extracts it there, yielding ``data/external/provenzano2025/xlinker/`` with
the LAMMPS structures (``data.relaxed00``, ``data.xlinker``) and the
crosslink-trend example (``xl_trend.txt``) used by
``scripts/compare_epoxy_external.py``.

The archive is NOT committed to the repo (``data/`` is gitignored); this
script makes the download reproducible instead (specs/decisions.md
2026-07-12 Track 2 / E2 design, item (d)).  An ``ATTRIBUTION.md`` with the
CC-BY-4.0 citation is written next to the extracted tree.

The code inside the archive (``xlinker.py``, ``lib/``) is reference material
only — it is never executed, modified, or vendored by this repo.

Usage:
    python scripts/fetch_provenzano2025.py                         # download
    python scripts/fetch_provenzano2025.py --from-archive x.tgz    # local tgz
    python scripts/fetch_provenzano2025.py --from-dir xlinker/     # local dir
"""
from __future__ import annotations

import argparse
import shutil
import sys
import tarfile
import urllib.request
from pathlib import Path

ZENODO_URL = 'https://zenodo.org/records/15418263/files/xlinker.tgz?download=1'
DEFAULT_DEST = Path('data/external/provenzano2025')

#: Files that must exist under <dest>/xlinker/ after extraction.
EXPECTED_FILES = (
    'xlinker.py',
    'xlpars.txt',
    'reaxpar.txt',
    'xl_trend.txt',
    'data.relaxed00',
    'data.xlinker',
    'lib',
)

ATTRIBUTION = """\
# Attribution — Provenzano 2025 crosslinker dataset

The contents of `xlinker/` were obtained from:

- G. Provenzano et al., *ACS Applied Polymer Materials* **2025**, 7(8), 4876.
  DOI: 10.1021/acsapm.4c04208
- Dataset: Zenodo, DOI [10.5281/zenodo.11402476](https://doi.org/10.5281/zenodo.11402476)
  (record https://zenodo.org/records/15418263, file `xlinker.tgz`)
- License: Creative Commons Attribution 4.0 International (CC-BY-4.0)

Usage in this repository (kagome): the LAMMPS data structures
(`data.relaxed00`, `data.xlinker`) and the crosslink-trend example
(`xl_trend.txt`) are READ as external reference data for the E2
cross-method comparison (specs/decisions.md 2026-07-12).  The crosslinker
code (`xlinker.py`, `lib/`) is not executed, modified, or redistributed.

Retrieved via `scripts/fetch_provenzano2025.py`.
"""


def _safe_extract(archive: Path, dest: Path) -> None:
    """Extract ``archive`` into ``dest``, refusing path-traversal members.

    Rejects absolute member names, ``..`` components, and links pointing
    outside ``dest`` (equivalent to Python 3.12's ``filter='data'`` intent,
    implemented explicitly for portability).
    """
    dest = dest.resolve()
    with tarfile.open(archive, 'r:*') as tar:
        for member in tar.getmembers():
            name = Path(member.name)
            if name.is_absolute() or '..' in name.parts:
                raise ValueError(
                    f'refusing unsafe archive member: {member.name!r}'
                )
            target = (dest / name).resolve()
            if not target.is_relative_to(dest):
                raise ValueError(
                    f'archive member escapes destination: {member.name!r}'
                )
            if member.islnk() or member.issym():
                raise ValueError(
                    f'refusing link member in archive: {member.name!r}'
                )
        try:
            tar.extractall(dest, filter='data')
        except TypeError:  # Python < 3.11.4 without the filter argument
            tar.extractall(dest)


def _download(url: str, target: Path) -> None:
    print(f'Downloading {url}\n  -> {target}')
    target.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as response, open(target, 'wb') as out:
        shutil.copyfileobj(response, out)
    print(f'Downloaded {target.stat().st_size} bytes.')


def _verify(xlinker_dir: Path) -> list[str]:
    """Names from EXPECTED_FILES missing under ``xlinker_dir``."""
    return [n for n in EXPECTED_FILES if not (xlinker_dir / n).exists()]


def _ensure_gitignored(repo_root: Path) -> None:
    """Add ``data/`` to .gitignore if it is not already ignored."""
    gitignore = repo_root / '.gitignore'
    lines = (
        gitignore.read_text(encoding='utf-8').splitlines()
        if gitignore.exists() else []
    )
    if any(line.strip() in ('data/', 'data', '/data/', '/data') for line in lines):
        return
    with open(gitignore, 'a', encoding='utf-8') as fh:
        fh.write('/data/\n')
    print(f'Added data/ to {gitignore}')


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    source = parser.add_mutually_exclusive_group()
    source.add_argument(
        '--from-archive', type=Path, metavar='PATH',
        help='use an already-downloaded xlinker.tgz instead of the network',
    )
    source.add_argument(
        '--from-dir', type=Path, metavar='PATH',
        help='copy an already-extracted xlinker/ directory tree',
    )
    parser.add_argument(
        '--dest', type=Path, default=DEFAULT_DEST,
        help=f'destination directory (default: {DEFAULT_DEST})',
    )
    args = parser.parse_args(argv)

    dest: Path = args.dest
    dest.mkdir(parents=True, exist_ok=True)

    if args.from_dir is not None:
        src_dir: Path = args.from_dir
        if not src_dir.is_dir():
            print(f'ERROR: --from-dir {src_dir} is not a directory')
            return 1
        target_dir = dest / 'xlinker'
        if target_dir.exists():
            shutil.rmtree(target_dir)
        shutil.copytree(src_dir, target_dir)
        print(f'Copied {src_dir} -> {target_dir}')
    else:
        archive = args.from_archive
        if archive is not None:
            if not archive.is_file():
                print(f'ERROR: --from-archive {archive} is not a file')
                return 1
            local = dest / 'xlinker.tgz'
            if archive.resolve() != local.resolve():
                shutil.copy2(archive, local)
                print(f'Copied {archive} -> {local}')
        else:
            local = dest / 'xlinker.tgz'
            _download(ZENODO_URL, local)
        _safe_extract(local, dest)
        print(f'Extracted {local} -> {dest}')

    xlinker_dir = dest / 'xlinker'
    missing = _verify(xlinker_dir)
    if missing:
        print(f'ERROR: expected files missing under {xlinker_dir}: {missing}')
        return 1

    (dest / 'ATTRIBUTION.md').write_text(ATTRIBUTION, encoding='utf-8')
    print(f'Wrote {dest / "ATTRIBUTION.md"}')

    repo_root = Path(__file__).resolve().parent.parent
    _ensure_gitignored(repo_root)

    print(f'OK: {xlinker_dir} populated ({len(EXPECTED_FILES)} expected entries present).')
    return 0


if __name__ == '__main__':
    sys.exit(main())
