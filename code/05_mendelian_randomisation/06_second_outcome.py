#!/usr/bin/env python3
"""
06_second_outcome.py
====================

Purpose
-------
Repeat the univariable fasting insulin MR against a second endometrial cancer
dataset, to establish whether the association depends on one outcome study.

The second dataset is GWAS_SECONDARY_FILE (GCST006465, O'Mara et al. 2018,
GRCh37). Because it is on the older assembly while the instrument set is on
GRCh38, the instruments are lifted back to GRCh37 here and matched on
coordinates; no rsID join is used, since the secondary file's identifier
convention differs from the primary one.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
config.PATHS.tables / instruments_fi.tsv      fasting insulin instrument set
config.PATHS.data / GWAS_SECONDARY_FILE       GCST006465, GRCh37, tabix-style TSV
config.PATHS.data / CHAIN_HG38_TO_HG19        chain file for the reverse lift
config.PARAMS.mr_palindromic_lo               ambiguous frequency band, lower
config.PARAMS.mr_palindromic_hi               ambiguous frequency band, upper

Outputs
-------
config.PATHS.tables / harmonised_second.tsv, columns chr, pos, fi_beta, fi_se,
ec2_beta, ec2_se, align: one row per matched variant, on the secondary
dataset's effect allele.

config.PATHS.tables / mr_second_outcome.tsv, columns outcome, method, beta, se,
p, OR, n: IVW, MR-Egger and weighted median for the secondary outcome.

Notes
-----
Interpretation. The two outcome datasets are not independent: both contain the
Endometrial Cancer Association Consortium case series. Agreement between them is
therefore consistency across a change of dataset, not replication, and the
manuscript states it that way. This script exists to make that check
reproducible, not to claim an independent validation.

Why only three estimators. MR-PRESSO and the leave-one-out analysis are not
repeated here: the second dataset is a consistency check on the direction and
magnitude of the estimate, and its instrument count is too small for a residual
screen to name a variant reliably.

Below five usable variants the secondary analysis is not reported at all; the
estimators would still return numbers, but the Egger fit has at most two residual
degrees of freedom and the weighted median bootstrap is degenerate.
"""
from __future__ import annotations

import gzip
import importlib
import pathlib
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# Repository root, so that ``config`` resolves when the file is run directly.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

from config.config import (  # noqa: E402  (import after sys.path edit)
    PATHS,
    PARAMS,
    SEED,
    CHAIN_HG38_TO_HG19,
    GWAS_SECONDARY_FILE,
)

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# The chain reader and the estimator library have numeric file names, which
# cannot appear in an import statement; importlib resolves them from the
# directory added above.
_select = importlib.import_module("02_select_instruments")
mr = importlib.import_module("01_mr_estimators")
Chain = _select.Chain
INSTRUMENTS_FI = _select.INSTRUMENTS_FI

SECOND_OUTCOME_NAME = "GCST006465"
HARMONISED_SECOND = "harmonised_second.tsv"
RESULTS_SECOND = "mr_second_outcome.tsv"

# Minimum number of matched variants for the secondary analysis to be reported.
MIN_VARIANTS = 5

PALINDROMIC = {frozenset("AT"), frozenset("CG")}
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def harmonise_second(ea1: str, oa1: str, b1: float, eaf1: float,
                     ea2: str, oa2: str) -> Optional[Tuple[float, str]]:
    """Realign the instrument effect onto the secondary outcome's effect allele.

    Parameters
    ----------
    ea1, oa1 : str
        Effect and other allele in the instrument (exposure) table.
    b1 : float
        Exposure effect size.
    eaf1 : float
        Allele frequency used for the palindromic guard.
    ea2, oa2 : str
        Effect and other allele in the secondary outcome file.

    Returns
    -------
    (beta, label) for an alignable variant, or None when the allele pair is
    incompatible or palindromic with an ambiguous frequency.

    Notes
    -----
    Identical decision logic to ``harmonise`` in 03_harmonise_datasets.py, but
    without the frequency bookkeeping: the instrument table carries no frequency
    and none of the estimates on this path use one, so only the sign of the
    effect has to be settled.

    The caller supplies the *outcome* frequency for the palindromic guard,
    because that is the only frequency available for this variant pair. Both
    frequencies describe the same variant, so the ambiguity test is unaffected;
    the value is only ever compared against the band, never used arithmetically.
    """
    if ea1 == ea2 and oa1 == oa2:
        return b1, "same"

    if ea1 == oa2 and oa1 == ea2:
        return -b1, "flip"

    if COMPLEMENT.get(ea1) == ea2 and COMPLEMENT.get(oa1) == oa2:
        if frozenset((ea1, oa1)) in PALINDROMIC and \
                PARAMS.mr_palindromic_lo <= eaf1 <= PARAMS.mr_palindromic_hi:
            return None
        return b1, "strand"

    if COMPLEMENT.get(ea1) == oa2 and COMPLEMENT.get(oa1) == ea2:
        if frozenset((ea1, oa1)) in PALINDROMIC and \
                PARAMS.mr_palindromic_lo <= eaf1 <= PARAMS.mr_palindromic_hi:
            return None
        return -b1, "strand_flip"

    return None     # incompatible allele pair


def lift_instruments(instruments: pd.DataFrame, chain: Chain
                     ) -> Dict[Tuple[str, int], Tuple[str, int, str, str, float, float]]:
    """Map the GRCh38 instrument positions back to GRCh37.

    Parameters
    ----------
    instruments : DataFrame
        Fasting insulin instrument table (chr, pos in GRCh38).
    chain : Chain
        Reader opened on the hg38ToHg19 chain file.

    Returns
    -------
    dict keyed by GRCh37 (chrom, pos) holding (GRCh38 chrom, GRCh38 pos,
    EA, OA, beta, se). Unmappable variants are dropped.

    Notes
    -----
    The key is the GRCh37 position because that is what the secondary outcome
    file is indexed on; the GRCh38 coordinates are carried through so that the
    output table keeps the build the rest of the pipeline uses.
    """
    needed: Dict[Tuple[str, int], Tuple[str, int, str, str, float, float]] = {}
    for row in instruments.itertuples(index=False):
        mapped = chain.lift(str(row.chr), int(row.pos))
        if mapped is None:
            continue
        needed[(str(mapped[0]), int(mapped[1]))] = (
            str(row.chr), int(row.pos), row.EA, row.OA, row.beta, row.se)
    return needed


def read_second_outcome(path: pathlib.Path,
                        needed: Dict[Tuple[str, int], Tuple[str, int, str, str, float, float]]
                        ) -> Dict[Tuple[str, int], Tuple[str, str, float, float, float]]:
    """Collect effects from the secondary outcome at the required positions.

    Parameters
    ----------
    path : GRCh37 summary statistics, gzip TSV.
    needed : mapping GRCh37 position -> instrument record.

    Returns
    -------
    dict keyed by GRCh37 position, holding (EA, OA, beta, se, eaf).
    """
    found: Dict[Tuple[str, int], Tuple[str, str, float, float, float]] = {}
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {name: i for i, name in enumerate(header)}
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            # This file labels chromosomes with the "chr" prefix, unlike the
            # primary outcome; the prefix is stripped to match the chain reader.
            key = (cols[ix["chromosome"]].replace("chr", ""),
                   int(cols[ix["base_pair_location"]]))
            if key not in needed:
                continue
            ea = cols[ix["effect_allele"]].upper()
            oa = cols[ix["other_allele"]].upper()
            if len(ea) != 1 or len(oa) != 1:
                continue
            try:
                found[key] = (ea, oa,
                              float(cols[ix["beta"]]),
                              float(cols[ix["standard_error"]]),
                              float(cols[ix["effect_allele_frequency"]]))
            except ValueError:
                continue
    return found


def main() -> None:
    """Harmonise the instrument set with the secondary outcome and estimate."""
    PATHS.ensure()

    instruments = pd.read_csv(PATHS.tables / INSTRUMENTS_FI, sep="\t")
    chain = Chain(PATHS.data / CHAIN_HG38_TO_HG19)

    needed = lift_instruments(instruments, chain)
    print(f"instruments: {len(instruments)}, mappable back to GRCh37: {len(needed)}",
          flush=True)

    outcome = read_second_outcome(PATHS.data / GWAS_SECONDARY_FILE, needed)
    print(f"secondary outcome: {len(outcome)} positions matched", flush=True)

    rows: List[Dict[str, object]] = []
    for key19, (chrom, pos, ea, oa, beta, se) in needed.items():
        if key19 not in outcome:
            continue
        e2, o2, b2, se2, eaf2 = outcome[key19]
        # eaf2 is passed as the frequency argument: see harmonise_second.
        aligned = harmonise_second(ea, oa, beta, eaf2, e2, o2)
        if aligned is None:
            continue
        rows.append(dict(chr=chrom, pos=pos, fi_beta=aligned[0], fi_se=se,
                         ec2_beta=b2, ec2_se=se2, align=aligned[1]))

    harmonised = pd.DataFrame(rows)
    harmonised.to_csv(PATHS.tables / HARMONISED_SECOND, sep="\t", index=False)
    print(f"harmonised {len(harmonised)} variants", flush=True)

    if len(harmonised) < MIN_VARIANTS:
        print(f"too few usable variants (< {MIN_VARIANTS}); "
              f"the secondary outcome is not evaluated", flush=True)
        return

    b_exp = harmonised.fi_beta.values
    b_out = harmonised.ec2_beta.values
    se_out = harmonised.ec2_se.values

    r_ivw = mr.ivw(b_exp, b_out, se_out, random_effects=True)
    r_egg = mr.egger(b_exp, b_out, se_out)
    r_wm = mr.weighted_median(b_exp, b_out, se_out)

    print(f"\n=== Secondary outcome {SECOND_OUTCOME_NAME} "
          f"(O'Mara 2018), {len(harmonised)} variants ===", flush=True)
    out_rows: List[Dict[str, object]] = []
    for method, res in (("IVW", r_ivw), ("MR-Egger", r_egg),
                        ("Weighted median", r_wm)):
        # The Egger fit reports beta_p (t test); the others report p.
        p = res["p"] if "p" in res else res["beta_p"]
        out_rows.append(dict(outcome=SECOND_OUTCOME_NAME, method=method,
                             beta=res["beta"], se=res["se"], p=p,
                             OR=float(np.exp(res["beta"])), n=len(harmonised)))
        print(f"  {method:16s} beta={res['beta']:+.4f} se={res['se']:.4f} "
              f"OR={np.exp(res['beta']):.3f} p={p:.3g}", flush=True)

    print(f"  heterogeneity Q={r_ivw['Q']:.1f} (df={r_ivw['Q_df']}) "
          f"p={r_ivw['Q_p']:.3g}  I2={r_ivw['I2']:.2f}", flush=True)

    pd.DataFrame(out_rows)[["outcome", "method", "beta", "se", "p", "OR", "n"]].to_csv(
        PATHS.tables / RESULTS_SECOND, sep="\t", index=False)


if __name__ == "__main__":
    main()
