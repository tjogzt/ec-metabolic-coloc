"""
02_gwas_reader.py
=================

Purpose
-------
Load a genome-wide association summary-statistics file into memory, discover
independent risk loci from it, and look up individual variants by position and
allele. Three public file layouts are handled:

    ``gwas_ssf``    GWAS Catalog standard format      (variant_id first)
    ``gwas_hm``     harmonised GCST release           (variant_id last)
    ``finngen``     FinnGen endpoint release          (#chrom, pos, ref, alt)

Column names are resolved through a layout table rather than hard-coded, so a
new release adds one entry instead of editing the analysis code.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
path        gzip-compressed summary statistics
layout      one of the keys of LAYOUTS (default: auto-detected from the header)
chroms      chromosomes to keep (default: 1..22)
seed        not used here; retained in the signature for API symmetry

Outputs
-------
load_gwas()  dict: chromosome -> DataFrame[pos, ref, alt, beta, se, pval],
             sorted by position and de-duplicated
peaks()      DataFrame of loci: chrom, start, end, lead, p
lookup()     (beta, se, pval) for one variant, allele-orientation aware

Notes
-----
Only single-nucleotide variants are retained: the colocalisation test joins on
exact position, so indels and multi-allelic records would silently fail to
match and shrink the test set without raising an error.
"""
from __future__ import annotations

import gzip
import pathlib
from typing import Dict, Iterable, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

# Column layout per source. Keys are the names used internally; values are the
# column headings found in the file.
LAYOUTS: Dict[str, Dict[str, str]] = {
    "gwas_ssf": {
        "chrom": "chromosome",
        "pos": "base_pair_location",
        "ref": "other_allele",
        "alt": "effect_allele",
        "beta": "beta",
        "se": "standard_error",
        "pval": "p_value",
    },
    "gwas_hm": {          # harmonised GCST: same names, variant_id moved to the end
        "chrom": "chromosome",
        "pos": "base_pair_location",
        "ref": "other_allele",
        "alt": "effect_allele",
        "beta": "beta",
        "se": "standard_error",
        "pval": "p_value",
    },
    "finngen": {
        "chrom": "#chrom",
        "pos": "pos",
        "ref": "ref",
        "alt": "alt",
        "beta": "beta",
        "se": "sebeta",
        "pval": "pval",
    },
}

_KEEP = ["pos", "ref", "alt", "beta", "se", "pval"]


def detect_layout(path: "str | pathlib.Path") -> str:
    """Identify the file layout from its header line.

    Parameters
    ----------
    path : gzip-compressed summary statistics.

    Returns
    -------
    Key into :data:`LAYOUTS`.

    Raises
    ------
    ValueError
        When the header matches no known layout, so that an unseen format fails
        loudly instead of being parsed under the wrong column names.
    """
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    if header and header[0].startswith("#chrom"):
        return "finngen"
    if "base_pair_location" in header and "variant_id" in header:
        return "gwas_hm" if header[-1] == "variant_id" else "gwas_ssf"
    raise ValueError(f"unrecognised GWAS layout, header starts with {header[:5]}")


def load_gwas(path, layout: Optional[str] = None,
              chunksize: int = 2_000_000) -> Dict[str, pd.DataFrame]:
    """Read a GWAS file into a per-chromosome dictionary.

    Parameters
    ----------
    path : gzip-compressed summary statistics.
    layout : layout key, or None to auto-detect from the header.
    chunksize : rows per read; keeps peak memory bounded on files with tens of
        millions of variants.

    Returns
    -------
    dict mapping chromosome (as str, no "chr" prefix) to a DataFrame with the
    columns listed in ``_KEEP``, sorted by position and de-duplicated on
    (pos, ref, alt).

    Notes
    -----
    Non-finite betas, non-positive standard errors and any variant whose
    alleles are not single bases are dropped here, once, so downstream code can
    assume clean input.
    """
    path = pathlib.Path(path)
    if layout is None:
        layout = detect_layout(path)
    colmap = LAYOUTS[layout]
    usecols = list(dict.fromkeys(colmap.values()))

    store: Dict[str, pd.DataFrame] = {}
    for chunk in pd.read_csv(path, sep="\t", comment="#", low_memory=False,
                             usecols=usecols, chunksize=chunksize):
        chunk = chunk.rename(columns={v: k for k, v in colmap.items()})
        chunk["chrom"] = chunk["chrom"].astype(str).str.replace("chr", "", regex=False)

        # Keep SNVs only, and only rows that can enter the ABF at all.
        chunk = chunk[chunk["ref"].astype(str).str.len().eq(1)
                      & chunk["alt"].astype(str).str.len().eq(1)]
        chunk = chunk[np.isfinite(chunk["beta"]) & np.isfinite(chunk["se"])
                      & (chunk["se"] > 0)]

        for chrom, sub in chunk.groupby("chrom"):
            sub = sub.sort_values("pos")[_KEEP]
            prior = store.get(chrom)
            store[chrom] = (sub if prior is None
                            else pd.concat([prior, sub]).sort_values("pos"))

    for chrom, d in store.items():
        store[chrom] = d.drop_duplicates(subset=["pos", "ref", "alt"],
                                         keep="first").reset_index(drop=True)
    return store


def peaks(store: Dict[str, pd.DataFrame], p_thresh: float,
          gap: int) -> pd.DataFrame:
    """Cluster genome-wide significant variants into loci.

    Parameters
    ----------
    store : output of :func:`load_gwas`.
    p_thresh : significance threshold for a variant to seed or extend a locus.
    gap : maximum distance in base pairs between consecutive significant
        variants of the same locus; beyond this a new locus starts.

    Returns
    -------
    DataFrame with one row per locus: chrom, start, end, lead, p. Sorted by p,
    so the strongest locus comes first.

    Notes
    -----
    The lead variant is the most significant variant in the window, not the
    first one encountered; this matters when the strongest signal sits at the
    edge of a block.
    """
    rows = []
    for chrom, d in store.items():
        d = d[d["pval"] < p_thresh].sort_values("pos")
        if d.empty:
            continue
        current = None
        for x in d.itertuples(index=False):
            if current is None or x.pos - current["end"] > gap:
                if current:
                    rows.append(current)
                current = {"chrom": chrom, "start": x.pos, "end": x.pos,
                           "lead": x.pos, "p": x.pval}
            else:
                current["end"] = x.pos
                if x.pval < current["p"]:
                    current["p"], current["lead"] = x.pval, x.pos
        if current:
            rows.append(current)
    return pd.DataFrame(rows).sort_values("p").reset_index(drop=True)


def lookup(store: Dict[str, pd.DataFrame], chrom: str, pos: int,
           ref: str, alt: str) -> Optional[Tuple[float, float, float]]:
    """Fetch one variant by position, resolving strand orientation automatically.

    Parameters
    ----------
    store : output of :func:`load_gwas`.
    chrom : chromosome, with or without the "chr" prefix.
    pos : 1-based position.
    ref, alt : alleles as given by the *query* dataset.

    Returns
    -------
    (beta, se, pval) oriented to the query's ref/alt, or None when the variant
    is absent.

    Notes
    -----
    If the stored alleles are (alt, ref) relative to the query, beta is negated
    so that the returned effect always refers to the query's alt allele. This is
    the step that makes cross-build joins safe; skipping it inverts the sign of
    every swapped variant without any error being raised.
    """
    d = store.get(str(chrom).replace("chr", ""))
    if d is None:
        return None
    i = int(np.searchsorted(d["pos"].values, pos, side="left"))
    best = None
    while i < len(d) and d["pos"].iat[i] == pos:
        stored_ref = str(d["ref"].iat[i]).upper()
        stored_alt = str(d["alt"].iat[i]).upper()
        if (stored_ref, stored_alt) == (ref, alt):
            best = (d["beta"].iat[i], d["se"].iat[i], d["pval"].iat[i])
        elif (stored_ref, stored_alt) == (alt, ref):
            best = (-d["beta"].iat[i], d["se"].iat[i], d["pval"].iat[i])
        if best is not None:
            break
        i += 1
    return best
