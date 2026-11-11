"""
01_check_inputs.py
==================

Purpose
-------
Verify that every input the pipeline expects is present, readable and the
expected size class, before any long-running stage starts.

A missing or truncated input is the most common reason a colocalisation run
returns no tests instead of an error: tabix reports success on an unindexed or
partially downloaded file, and the stage then writes an empty table. Checking
first turns that silent failure into a clear message.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
The data directory resolved by config.PATHS.data.

Outputs
-------
A report on stdout. Exit status 0 when every required input is present, 1 when
anything is missing, so that a shell chain can stop on failure.

Notes
-----
Sizes below are the values observed for the published run and are used only as a
sanity band (0.5x to 4x). A file far outside its band is much more likely to be a
truncated download than a legitimate alternative release.
"""
from __future__ import annotations

import gzip
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))

from config.config import (PATHS, TISSUES, GWAS_PRIMARY_FILE, GWAS_SECONDARY_FILE,
                           GWAS_FINNGEN_FILE, COHORTS)

# name -> approximate size in bytes, used as a ±4x sanity band.
REQUIRED = {
    GWAS_PRIMARY_FILE: 373 * 10**6,
    GWAS_SECONDARY_FILE: 247 * 10**6,
    "GTEx_Analysis_v8_eQTL.tar": 1.5 * 10**9,
}

OPTIONAL = {
    GWAS_FINNGEN_FILE: 773 * 10**6,
}

# UCSC over.chain files vary by two orders of magnitude depending on how the
# chains were collapsed, so they are checked for presence and integrity only,
# not against a size band.
CHAIN_FILES = ("hg19ToHg38.over.chain.gz", "hg38ToHg19.over.chain.gz")


def _sample_gzip(path: pathlib.Path, n: int = 4096) -> bool:
    """Return True if the first ``n`` bytes decompress without error.

    Parameters
    ----------
    path : gzip-compressed file.
    n : bytes to decompress for the probe.

    Returns
    -------
    True when the header and first block decode, which catches the common
    "HTML error page saved as .gz" failure.
    """
    try:
        with gzip.open(path, "rb") as fh:
            fh.read(n)
        return True
    except Exception:
        return False


def check_file(name: str, expected: float) -> bool:
    """Report one input: presence, size band and gzip integrity."""
    path = PATHS.data / name
    if not path.exists():
        print(f"  MISSING   {name}")
        return False
    size = path.stat().st_size
    lo, hi = expected * 0.5, expected * 4.0
    flag = "ok      " if lo <= size <= hi else "SIZE    "
    print(f"  {flag} {name:<45s} {size / 1e6:10.1f} MB")
    if path.suffix == ".gz" and not _sample_gzip(path):
        print(f"            -> does not decompress; likely a truncated download")
        return False
    return lo <= size <= hi


def main() -> int:
    print(f"data directory: {PATHS.data}")
    if not PATHS.data.exists():
        print("  directory does not exist; set ECIR_DATA or create it")
        return 1

    ok = True
    print("\nrequired:")
    for name, size in REQUIRED.items():
        ok &= check_file(name, size)

    print("\neQTL datasets (one per tissue):")
    for tissue, dataset in TISSUES.items():
        path = PATHS.data / f"{dataset}.all.tsv.gz"
        if path.exists():
            print(f"  ok       {tissue:<28s} {dataset} "
                  f"{path.stat().st_size / 1e9:5.2f} GB")
        else:
            # Adipose is expected to be absent in this package and is fetched
            # remotely by stage 02 with --remote.
            print(f"  absent   {tissue:<28s} {dataset} "
                  f"(stage 02 must be run with --remote)")

    print("\noptional:")
    for name, size in OPTIONAL.items():
        check_file(name, size)
    for name in CHAIN_FILES:
        path = PATHS.data / name
        if path.exists() and _sample_gzip(path):
            print(f"  ok       {name:<45s} {path.stat().st_size / 1e3:10.0f} KB")
        else:
            print(f"  MISSING   {name}")

    print("\nbulk cohorts:")
    for cohort in COHORTS:
        hits = list(PATHS.data.glob(f"{cohort.gse}*"))
        print(f"  {'ok' if hits else 'absent':<8} {cohort.gse:<12s} "
              f"{cohort.tissue:<28s} ({cohort.phenotype})")

    print("\nresult:", "all required inputs present" if ok else "inputs missing")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
