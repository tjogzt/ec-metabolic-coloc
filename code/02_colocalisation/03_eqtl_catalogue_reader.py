"""
03_eqtl_catalogue_reader.py
===========================

Purpose
-------
Read cis-eQTL summary statistics from an eQTL Catalogue dataset. The file holds
every variant tested in each gene's cis window, not only the significant pairs,
which is exactly what the ABF colocalisation test requires.

Two access paths are provided:

    qtl_region()     tabix query for one interval (local file or remote URL)
    stream_regions() one sequential pass over a local file, retaining only the
                     intervals of interest

The streaming path exists because four of the five tissue files were shipped
without a tabix index; on those, ``tabix`` returns an empty result and exits
successfully, which would silently produce zero colocalisation tests instead of
an error.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
path      eQTL Catalogue *.{dataset}.all.tsv.gz
regions   dict of chromosome -> list of (start, end) half-open intervals
tissue    GTEx tissue label; resolved to a dataset id through config.TISSUES

Outputs
-------
DataFrame with the columns gene_id, chromosome, position, ref, alt, variant,
maf, pvalue, beta, se, median_tpm, rsid (numeric columns coerced), or an empty
DataFrame when nothing falls in the requested window.

Notes
-----
Column positions are resolved from the header line rather than being hard
coded, so a change in the upstream column order cannot shift the meaning of
``beta`` and ``se`` unnoticed.
"""
from __future__ import annotations

import gzip
import pathlib
import subprocess
from bisect import bisect_right
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

# Columns retained from the upstream file, in output order.
WANTED = ["molecular_trait_id", "chromosome", "position", "ref", "alt", "variant",
          "maf", "pvalue", "beta", "se", "median_tpm", "rsid"]
OUT_COLS = ["gene_id", "chromosome", "position", "ref", "alt", "variant",
            "maf", "pvalue", "beta", "se", "median_tpm", "rsid"]
NUMERIC = ["position", "maf", "pvalue", "beta", "se", "median_tpm"]


def read_header(path: "str | pathlib.Path") -> List[str]:
    """Return the column names of an eQTL Catalogue summary file."""
    with gzip.open(path, "rt") as fh:
        return fh.readline().rstrip("\n").split("\t")


def column_positions(path) -> List[int]:
    """Map :data:`WANTED` onto the file's actual column indices.

    Parameters
    ----------
    path : eQTL Catalogue summary statistics.

    Returns
    -------
    List of 0-based column indices, aligned with :data:`WANTED`.

    Raises
    ------
    KeyError
        If a required column is missing, so that the failure happens at read
        time rather than producing a silently truncated table.
    """
    header = read_header(path)
    return [header.index(name) for name in WANTED]


def _to_frame(rows: List[List[str]], want: Sequence[int]) -> pd.DataFrame:
    """Build a typed DataFrame from raw tab-delimited rows."""
    if not rows:
        return pd.DataFrame(columns=OUT_COLS)
    df = pd.DataFrame(rows).iloc[:, list(want)]
    df.columns = OUT_COLS
    for col in NUMERIC:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def qtl_region(chrom: str, start: int, end: int,
               src: "str | pathlib.Path",
               want: Optional[Sequence[int]] = None) -> pd.DataFrame:
    """Fetch one interval with tabix.

    Parameters
    ----------
    chrom : chromosome, without the "chr" prefix.
    start, end : 1-based inclusive interval.
    src : local path or remote URL accepted by tabix.
    want : column indices from :func:`column_positions`; computed if omitted.

    Returns
    -------
    DataFrame over :data:`OUT_COLS`, empty if the query returns nothing.

    Notes
    -----
    A non-zero tabix exit status and an empty result are treated identically:
    both yield an empty frame. The streaming alternative is the correct path
    when a file lacks an index.
    """
    if want is None:
        want = column_positions(src)
    proc = subprocess.run(["tabix", str(src), f"{chrom}:{start}-{end}"],
                          capture_output=True, text=True)
    if proc.returncode != 0 or not proc.stdout.strip():
        return pd.DataFrame(columns=OUT_COLS)
    rows = [line.split("\t") for line in proc.stdout.strip().split("\n")]
    return _to_frame(rows, want)


def stream_regions(path: "str | pathlib.Path",
                   regions: Dict[str, List[Tuple[int, int]]],
                   out_path: Union[str, pathlib.Path],
                   want: Optional[Sequence[int]] = None) -> int:
    """Single sequential pass, keeping rows that fall inside any interval.

    Parameters
    ----------
    path : local eQTL Catalogue summary file.
    regions : chromosome -> list of (start, end) intervals. Lists must be
        sorted; they are sorted here defensively.
    out_path : where to write the retained subset, as a header-led TSV.
    want : column indices, computed from ``path`` if omitted.

    Returns
    -------
    Number of rows written.

    Notes
    -----
    Consecutive lookups use a bisect over interval starts rather than a linear
    scan, so the cost is one pass over the file regardless of how many loci are
    requested. The temporary subset is written once and reused, which is why the
    caller may skip a repeated scan.
    """
    if want is None:
        want = column_positions(path)
    starts = {c: sorted(r[0] for r in v) for c, v in regions.items()}
    for c, v in regions.items():
        regions[c] = sorted(v)

    kept = 0
    with gzip.open(path, "rt") as fh, open(out_path, "w") as out:
        out.write("\t".join(OUT_COLS) + "\n")
        fh.readline()                      # discard the upstream header
        for line in fh:
            parts = line.split("\t", 4)
            if len(parts) < 5:
                continue
            chrom = parts[1]
            spans = regions.get(chrom)
            if spans is None:
                continue
            try:
                pos = int(parts[2])
            except ValueError:
                continue
            i = bisect_right(starts[chrom], pos) - 1
            if i >= 0 and pos <= spans[i][1]:
                fields = line.rstrip("\n").split("\t")
                out.write("\t".join(fields[j] for j in want) + "\n")
                kept += 1
    return kept


def build_gene_index(path: "str | pathlib.Path",
                     chunksize: int = 2_000_000) -> pd.DataFrame:
    """Derive each gene's cis window from the file itself.

    Parameters
    ----------
    path : local eQTL Catalogue summary file.
    chunksize : rows per read, to bound memory.

    Returns
    -------
    DataFrame with gene_id, chrom, lo, hi: the minimum and maximum tested
    position for each gene, which is the cis window the upstream pipeline used.
    """
    acc: Dict[str, Tuple[str, int, int]] = {}
    for chunk in pd.read_csv(path, sep="\t",
                             usecols=["molecular_trait_id", "chromosome", "position"],
                             dtype={"chromosome": str}, chunksize=chunksize):
        for gene, sub in chunk.groupby("molecular_trait_id"):
            chrom = sub["chromosome"].iloc[0]
            lo, hi = int(sub["position"].min()), int(sub["position"].max())
            if gene in acc:
                _, lo0, hi0 = acc[gene]
                acc[gene] = (chrom, min(lo, lo0), max(hi, hi0))
            else:
                acc[gene] = (chrom, lo, hi)
    return pd.DataFrame([(g, v[0], v[1], v[2]) for g, v in acc.items()],
                        columns=["gene_id", "chrom", "lo", "hi"])


def clean_eqtl(df: pd.DataFrame, min_snps: int, maf_lo: float,
               maf_hi: float) -> pd.DataFrame:
    """Apply the shared variant filters used before every colocalisation test.

    Parameters
    ----------
    df : output of :func:`qtl_region` or a streamed subset.
    min_snps : minimum number of usable variants for the gene to be tested.
    maf_lo, maf_hi : minor-allele-frequency window.

    Returns
    -------
    The filtered frame, possibly empty.

    Notes
    -----
    Rows without a finite beta or se, with non-positive se, or with a MAF
    outside the window are removed here so that the ABF is never fed a variant
    whose standard error would invert or divide by zero.
    """
    df = df.dropna(subset=["position", "beta", "se"])
    df = df[(df["maf"] >= maf_lo) & (df["maf"] <= maf_hi) & (df["se"] > 0)]
    return df if len(df) >= min_snps else df.iloc[0:0]
