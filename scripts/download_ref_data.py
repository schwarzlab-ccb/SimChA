#!/usr/bin/env python3
"""Download reference genome FASTA files for SimChA.

Downloads the GRCh37 (hg19) and/or GRCh38 (hg38) reference genomes from the
Google ``genomics-public-data`` bucket, normalizes the chromosome headers to
UCSC style (``>chr1`` ...), and places the result at ``data/<assembly>/genome.fa``.

These files are required when running SimChA with ``-v``/``--variants`` or
``-f``/``--fasta`` for the selected assembly.

Run from the SimChA repository root:

    python scripts/download_ref_data.py            # both assemblies
    python scripts/download_ref_data.py --assembly hg19
"""

from __future__ import annotations

import argparse
import gzip
import re
import shutil
import sys
import tempfile
import urllib.error
import urllib.request
from pathlib import Path

HG19_BASE_API = "https://storage.googleapis.com/genomics-public-data/references/GRCh37/"
HG19_CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]
HG38_URL = (
    "https://storage.googleapis.com/genomics-public-data/references/"
    "GRCh38_Verily/GRCh38_Verily_v1.genome.fa"
)

# hg19 is stored per-chromosome with NCBI-style FASTA headers such as
# ">gi|...|gb|...| Homo sapiens chromosome N, ...". Convert those to ">chrN".
_NCBI_HEADER = re.compile(
    r"^>gi\|\d+\|gb\|[^|]+\| Homo sapiens chromosome ([0-9XY]+),.*"
)
_BARE_HEADER = re.compile(r"^>([0-9XY]+)\b")


def _default_data_dir() -> Path:
    """The repository ``data/`` folder, relative to this script."""
    return Path(__file__).resolve().parent.parent / "data"


def _download(url: str, dest: Path) -> None:
    """Download ``url`` to ``dest``, streaming to avoid loading it in memory."""
    print(f"Downloading {url}", flush=True)
    with urllib.request.urlopen(url) as response, dest.open("wb") as out:
        shutil.copyfileobj(response, out)


def _normalize_header(line: str) -> str:
    """Rewrite an NCBI/bare FASTA header line to UCSC ``>chrN`` style."""
    match = _NCBI_HEADER.match(line)
    if match:
        return f">chr{match.group(1)}\n"
    match = _BARE_HEADER.match(line)
    if match:
        return f">chr{match.group(1)}\n"
    return line


def download_hg19(data_dir: Path) -> None:
    """Assemble ``data/hg19/genome.fa`` from the per-chromosome hg19 files."""
    out_dir = data_dir / "hg19"
    out_dir.mkdir(parents=True, exist_ok=True)
    genome_fa = out_dir / "genome.fa"

    with genome_fa.open("w") as genome:
        for chrom in HG19_CHROMS:
            url = f"{HG19_BASE_API}{chrom}.fa.gz"
            with tempfile.NamedTemporaryFile(suffix=".fa.gz", delete=False) as tmp:
                tmp_path = Path(tmp.name)
            try:
                _download(url, tmp_path)
                # Decompress and append, fixing headers as we go.
                with gzip.open(tmp_path, "rt") as fa:
                    for line in fa:
                        genome.write(
                            _normalize_header(line) if line.startswith(">") else line
                        )
            finally:
                tmp_path.unlink(missing_ok=True)
    print(f"Wrote {genome_fa}", flush=True)


def download_hg38(data_dir: Path) -> None:
    """Download ``data/hg38/genome.fa`` (single pre-assembled FASTA)."""
    out_dir = data_dir / "hg38"
    out_dir.mkdir(parents=True, exist_ok=True)
    genome_fa = out_dir / "genome.fa"
    _download(HG38_URL, genome_fa)
    print(f"Wrote {genome_fa}", flush=True)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-a",
        "--assembly",
        choices=["hg19", "hg38", "both"],
        default="both",
        help="Which reference assembly to download (default: both).",
    )
    parser.add_argument(
        "-d",
        "--data-dir",
        type=Path,
        default=_default_data_dir(),
        help="Target data directory (default: the repository 'data/' folder).",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if args.assembly in ("hg19", "both"):
            download_hg19(args.data_dir)
        if args.assembly in ("hg38", "both"):
            download_hg38(args.data_dir)
    except (urllib.error.URLError, OSError) as err:
        print(f"Error: {err}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
