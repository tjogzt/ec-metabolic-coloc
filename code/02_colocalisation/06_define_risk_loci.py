"""
06_define_risk_loci.py
======================

Purpose
-------
Turn the outcome GWAS into a set of independent risk loci, and pair each locus
with the genes that have a significant cis-eQTL in a given GTEx tissue.

The gene list comes from the GTEx significant variant-gene pairs archive, which
ships as a single tar of per-tissue tables. Only the one member that is needed
is read; the archive is never unpacked.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
tar_path   GTEx_Analysis_v8_eQTL.tar (or an archive with the same layout)
tissue     GTEx tissue label, e.g. "Uterus"
loci       DataFrame from 02_gwas_reader.peaks(), columns chrom/start/end
pad        half-width of the window placed around each locus, in base pairs

Outputs
-------
genes_near_loci() returns dict locus_index -> sorted list of Ensembl gene ids.

Notes
-----
A gene is assigned to a locus when any of its significant eQTL variants falls
within the padded window. The same gene may therefore appear at two loci when
two signals sit close together, which is intended: the locus, not the gene, is
the unit of the colocalisation test.

The archive must stay open while the member stream is being consumed, because
the returned stream is a view onto the tar's file object. The member is
therefore located and read inside a single ``with tarfile.open(...)`` block.
"""
from __future__ import annotations

import gzip
import pathlib
import tarfile
from typing import Dict, List, Optional, Tuple

import pandas as pd

VARIANT_COL = "variant_id"
GENE_COL = "gene_id"


def _find_member(tar: tarfile.TarFile, tissue: str) -> Optional[tarfile.TarInfo]:
    """Locate the significant-pairs member for one tissue.

    Parameters
    ----------
    tar : open tar handle; the caller owns its lifetime.
    tissue : GTEx tissue label without spaces.

    Returns
    -------
    The member, or None when the tissue is not present in the archive.
    """
    for member in tar.getmembers():
        if "signif_variant_gene_pairs" not in member.name:
            continue
        if pathlib.Path(member.name).name.startswith(tissue + "."):
            return member
    return None


def genes_near_loci(tar_path, tissue: str, loci: pd.DataFrame,
                    pad: int) -> Dict[int, List[str]]:
    """Collect the eGenes that fall inside a padded window around each locus.

    Parameters
    ----------
    tar_path : GTEx significant-pairs archive.
    tissue : GTEx tissue label.
    loci : DataFrame with columns chrom, start, end.
    pad : window half-width in base pairs.

    Returns
    -------
    dict mapping locus row index to a sorted list of unique gene ids.

    Raises
    ------
    SystemExit
        If the tissue is absent. Returning an empty mapping instead would look
        like "no signal" rather than "wrong tissue label", which is a silent
        failure mode this package avoids throughout.
    """
    # Windows are computed once; the member has millions of rows and the
    # comparison runs per row.
    windows: List[Tuple[int, str, int, int]] = [
        (i, str(row["chrom"]).replace("chr", ""),
         int(row["start"]) - pad, int(row["end"]) + pad)
        for i, row in loci.iterrows()
    ]
    out: Dict[int, List[str]] = {i: [] for i in range(len(loci))}

    with tarfile.open(tar_path, "r") as tar:
        member = _find_member(tar, tissue)
        if member is None:
            raise SystemExit(f"tissue {tissue} not found in {pathlib.Path(tar_path).name}")
        stream = gzip.GzipFile(fileobj=tar.extractfile(member))
        for line in stream:
            fields = line.decode().split("\t")
            if fields[0] == VARIANT_COL:
                continue
            # variant_id is chr<pos>_<ref>_<alt>_b38 in GTEx v8, but the suffix
            # is not guaranteed across releases; take the first two fields by
            # index rather than unpacking a fixed arity.
            parts = fields[0].split("_")
            # GTEx writes "chr17_48127316_G_A_b38"; the locus table carries bare
            # chromosome names, so the prefix is stripped on this side too.
            chrom, pos = parts[0].replace("chr", ""), int(parts[1])
            for i, win_chrom, lo, hi in windows:
                if chrom == win_chrom and lo <= pos <= hi:
                    out[i].append(fields[1].split(".")[0])   # drop Ensembl version

    return {k: sorted(set(v)) for k, v in out.items()}
