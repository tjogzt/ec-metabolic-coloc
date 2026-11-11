"""
04_finngen_scan.py
==================

Purpose
-------
Scan a FinnGen endpoint release and report its strongest association, so that a
trait can be excluded as an outcome on objective grounds rather than by
assertion.

This is the check that excluded the FinnGen R12 endometrial cancer endpoint from
the colocalisation analysis: across roughly 21 million variants no single
variant reached genome-wide significance, so any colocalisation against that
endpoint would have measured power rather than shared causal architecture.

The statistic produced here is computed from the released file by this package.
It is not a value reported in the FinnGen resource paper, and is never cited as
one.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
path       FinnGen endpoint file (columns #chrom, pos, ref, alt, pval, ...)
p_thresh   genome-wide significance threshold from config.PARAMS
chunksize  rows per read, to bound memory

Outputs
-------
scan() returns a dict:

    n_variants      variants read after filtering to SNVs
    min_p           strongest association observed
    n_significant   variants below p_thresh
    lead            dict of chrom/pos/ref/alt at the strongest variant

Notes
-----
The file is read with the shared GWAS layout table, so a change in the FinnGen
column order is a one-line edit in 02_gwas_reader.LAYOUTS rather than a change
here.
"""
from __future__ import annotations

import pathlib
import sys
from typing import Dict

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from importlib import import_module

_gwas = import_module("02_gwas_reader")


def scan(path, p_thresh: float, chunksize: int = 2_000_000) -> Dict[str, object]:
    """Find the strongest association in a FinnGen endpoint file.

    Parameters
    ----------
    path : FinnGen endpoint release, gzip-compressed.
    p_thresh : threshold used to count genome-wide significant variants.
    chunksize : rows per read.

    Returns
    -------
    dict with n_variants, min_p, n_significant and the lead variant's alleles.

    Notes
    -----
    Only SNVs with a finite p-value are counted, matching the filtering applied
    when the file is used as a colocalisation outcome; otherwise a rare indel
    could decide whether the endpoint is declared usable.
    """
    import pandas as pd

    store = _gwas.load_gwas(path, layout="finngen", chunksize=chunksize)
    n = sum(len(v) for v in store.values())
    if n == 0:
        return {"n_variants": 0, "min_p": float("nan"), "n_significant": 0, "lead": None}

    best = None
    n_sig = 0
    for chrom, d in store.items():
        p = d["pval"].to_numpy(dtype=float)
        p = p[np.isfinite(p)]
        if p.size == 0:
            continue
        n_sig += int((p < p_thresh).sum())
        i = int(np.argmin(p))
        if best is None or p[i] < best["p"]:
            row = d.iloc[i]
            best = {"chrom": chrom, "pos": int(row["pos"]), "ref": str(row["ref"]),
                    "alt": str(row["alt"]), "p": float(p[i])}

    return {"n_variants": n, "min_p": best["p"], "n_significant": n_sig,
            "lead": {k: best[k] for k in ("chrom", "pos", "ref", "alt")}}


def main() -> None:
    """CLI entry point: scan a file named in config and print the summary."""
    from config.config import PATHS, PARAMS, GWAS_FINNGEN_FILE

    path = PATHS.data / GWAS_FINNGEN_FILE
    if not path.exists():
        raise SystemExit(f"missing FinnGen file: {path}")
    res = scan(path, PARAMS.gwas_p_thresh)
    print(f"variants scanned : {res['n_variants']:,}")
    print(f"strongest p      : {res['min_p']:.3e}")
    print(f"significant      : {res['n_significant']}")
    if res["n_significant"] == 0:
        print("verdict          : excluded as an outcome (no genome-wide signal)")


if __name__ == "__main__":
    main()
