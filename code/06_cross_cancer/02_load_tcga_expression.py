#!/usr/bin/env python3
"""
02_load_tcga_expression.py
==========================

Purpose
-------
Load the TCGA pan-cancer gene-expression matrix, keep the samples that are
primary tumours or solid-tissue normals, drop genes that are not expressed, and
cache the result as a compressed array for the scoring stages.

The matrix is streamed rather than parsed as a whole: the full file is roughly
20,000 genes by 11,000 samples, which pandas would materialise densely in
memory. Streaming keeps peak memory proportional to the retained columns.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
data/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
    UCSC Xena pan-cancer adjusted gene expression (Illumina HiSeq RNASeqV2).
    The first field of each row is a HUGO gene symbol; the header carries the
    sample ids.
data/TCGA_phenotype_denseDataOnlyDownload.tsv.gz
    Xena phenotype table with the columns ``sample``, ``sample_type`` and
    ``_primary_disease``. The file name may be overridden through
    ``config.config.PAN_CANCER_PHENO_FILE`` when a local mirror is named
    differently.

Outputs
-------
<PATHS.cache>/expr.npz
    ``genes``    object array of gene symbols, in matrix row order
    ``samples``  object array of sample ids, in matrix column order
    ``X``        float32 array, genes x samples
<PATHS.cache>/meta.tsv
    one row per retained sample; columns sample_id / sample_type / cancer_type

Usage
-----
    python3 code/06_cross_cancer/02_load_tcga_expression.py

Notes
-----
* Tumours and normals come from the same matrix, hence the same processing
  pipeline and platform. This is deliberately not a GTEx comparison: a
  cross-resource contrast would confound tissue processing with tissue type,
  which is exactly the comparison the normal-tissue control is meant to keep
  clean.
* The low-expression filter keeps a gene whose mean across the retained samples
  is at least 0.5 (log2 scale). The threshold describes the matrix rather than
  the analysis, so it is fixed here and not exposed in config.PARAMS.
* Rows with too few fields, or with a non-numeric value in a retained column,
  are skipped rather than aborting the load: Xena ships a few annotation-like
  rows that would otherwise fail the cast.
"""
from __future__ import annotations

import gzip
import pathlib
import sys
import time
from typing import IO, List, Tuple

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from config.config import PATHS  # noqa: E402
import config.config as _config  # noqa: E402

# Matrix file name as documented in data/README.md.
EXPR_FILE = "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"

# Phenotype file name; overridable from config so the script carries no literal
# path of its own.
PHENO_FILE = getattr(_config, "PAN_CANCER_PHENO_FILE",
                     "TCGA_phenotype_denseDataOnlyDownload.tsv.gz")

# Sample types entering the analysis: tumours for the main scan, solid-tissue
# normals for the same-platform control.
SAMPLE_TYPES = ("Primary Tumor", "Solid Tissue Normal")

# Mean expression below which a gene counts as unexpressed (log2 scale).
MIN_MEAN_EXPRESSION = 0.5


def read_phenotype(path: pathlib.Path) -> pd.DataFrame:
    """Read the Xena phenotype table and normalise its column names.

    Releases disagree on the column titles, so the two that matter are renamed
    to a fixed vocabulary and everything downstream can ignore the difference.

    Parameters
    ----------
    path : path to the phenotype table (gzipped TSV).

    Returns
    -------
    DataFrame restricted to the two sample types of interest.
    """
    ph = pd.read_csv(path, sep="\t")
    ph = ph.rename(columns={"sample": "sample_id", "_primary_disease": "cancer_type"})
    missing = {"sample_id", "sample_type", "cancer_type"} - set(ph.columns)
    if missing:
        raise SystemExit(f"phenotype table lacks column(s): {sorted(missing)}")
    return ph[ph.sample_type.isin(SAMPLE_TYPES)]


def open_matrix(path: pathlib.Path) -> Tuple[IO[str], List[str]]:
    """Open the matrix and return the handle plus its header field names.

    The caller must close the handle. Splitting the header off here mirrors the
    original two-phase read: column positions are resolved once, then every row
    is a pure index selection.
    """
    fh = gzip.open(path, "rt", errors="replace")
    header = fh.readline().rstrip("\n").split("\t")
    return fh, header[1:]


def row_values(line: str, n_cols: int, idx_keep: List[int]) -> Tuple[str, np.ndarray] | None:
    """Extract (symbol, values) from one line, or None when it is unusable.

    ``None`` covers both a truncated row and a row whose retained fields are not
    all numeric; the caller simply skips it.
    """
    p = line.rstrip("\n").split("\t")
    if len(p) < n_cols + 1:
        return None
    try:
        vals = np.array([float(p[1 + i]) for i in idx_keep], dtype=np.float32)
    except ValueError:
        return None
    if not np.isfinite(vals).all():
        # Missing values are treated as absent signal. The convention is kept
        # from the source pipeline so that scores stay comparable.
        vals = np.nan_to_num(vals, nan=0.0)
    return p[0], vals


def main() -> None:
    PATHS.ensure()
    t0 = time.time()

    expr_path = PATHS.data / EXPR_FILE
    pheno_path = PATHS.data / PHENO_FILE
    for p in (expr_path, pheno_path):
        if not p.exists():
            raise SystemExit(f"missing input: {p}\nsee data/README.md for how to obtain it")

    # ---- phenotype ------------------------------------------------------
    ph = read_phenotype(pheno_path)
    if ph.empty:
        raise SystemExit("phenotype table holds no primary tumour or normal sample")

    keep_samples = set(ph["sample_id"])
    n_tumour = int((ph.sample_type == "Primary Tumor").sum())
    n_normal = int((ph.sample_type == "Solid Tissue Normal").sum())
    print(f"Phenotype: {len(ph)} samples (tumour {n_tumour} / normal {n_normal}), "
          f"{ph.cancer_type.nunique()} cancer types", flush=True)
    ph.to_csv(PATHS.cache / "meta.tsv", sep="\t", index=False)

    # ---- expression matrix, streamed ------------------------------------
    fh, cols = open_matrix(expr_path)
    # Resolve every wanted sample onto its column position once, so that the
    # per-row work is a plain index selection.
    idx_keep = [i for i, c in enumerate(cols) if c in keep_samples]
    keep_names = [cols[i] for i in idx_keep]
    print(f"Matrix: {len(cols)} columns, {len(idx_keep)} retained", flush=True)

    genes: List[str] = []
    blocks: List[np.ndarray] = []
    n_seen = 0
    try:
        for line in fh:
            n_seen += 1
            parsed = row_values(line, len(cols), idx_keep)
            if parsed is None:
                continue
            symbol, vals = parsed
            # Low-expression filter: unexpressed genes would only enlarge the
            # decile-matched null pool without carrying signal.
            if vals.mean() < MIN_MEAN_EXPRESSION:
                continue
            genes.append(symbol)
            blocks.append(vals)
            if len(genes) % 2000 == 0:
                print(f"  {len(genes)} genes kept of {n_seen} seen ... "
                      f"{time.time() - t0:.0f}s", flush=True)
    finally:
        fh.close()

    if not blocks:
        raise SystemExit("no gene passed the expression filter; check the matrix file")

    # Stack once after the loop: vstack over a list of small rows is cheap,
    # whereas growing one array row by row would reallocate every iteration.
    X = np.vstack(blocks)
    print(f"Retained {X.shape[0]} genes x {X.shape[1]} samples "
          f"({time.time() - t0:.0f}s)", flush=True)

    np.savez_compressed(
        PATHS.cache / "expr.npz",
        genes=np.array(genes, dtype=object),
        samples=np.array(keep_names, dtype=object),
        X=X,
    )
    print(f"Cached {PATHS.cache / 'expr.npz'} ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
