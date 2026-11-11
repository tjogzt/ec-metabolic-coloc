#!/usr/bin/env python3
"""
03_harmonise_datasets.py
========================

Purpose
-------
Pair every selected instrument with its effect on the outcome and with its
effect on the second exposure (BMI), on a common effect allele, and write the
single harmonised table that every downstream MR model reads.

This is the step where the three data sources are reconciled:

    exposure 1  fasting insulin, GCST90002238, GRCh37, no rsIDs
    exposure 2  BMI,             GCST006900,   GRCh37, with rsIDs
    outcome     endometrial cancer, GWAS_PRIMARY_FILE, GRCh38, with rsIDs

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
config.PATHS.tables / instruments_fi.tsv      from 02_select_instruments.py
config.PATHS.tables / instruments_bmi.tsv     from 02_select_instruments.py
config.PATHS.data / FI_EXPOSURE_FILE          full fasting insulin summary stats
config.PATHS.data / BMI_EXPOSURE_FILE         full BMI summary stats
config.PATHS.data / GWAS_PRIMARY_FILE         outcome summary stats (GRCh38)
config.PATHS.data / CHAIN_HG38_TO_HG19        chain file for the reverse lift
config.PARAMS.mr_palindromic_lo               lower bound of the ambiguous
config.PARAMS.mr_palindromic_hi               frequency band (0.42 / 0.58)

Outputs
-------
config.PATHS.tables / harmonised.tsv, one row per variant with both an outcome
effect and a fasting insulin effect:

    chr, pos                        GRCh38 coordinates
    rsid                            outcome-file identifier
    ec_EA, ec_OA, ec_beta, ec_se    outcome effect, on the outcome's effect allele
    ec_eaf
    fi_beta, fi_se, fi_eaf, fi_N    fasting insulin effect, realigned
    fi_align                        how the exposure was realigned (see below)
    bmi_beta, bmi_se, bmi_eaf       BMI effect, realigned (absent when the
    bmi_align                       variant is not in the BMI instrument set)

Notes
-----
Coordinate strategy. The exposure files are on GRCh37 and the outcome is on
GRCh38, so a naive implementation would lift tens of millions of rows. Instead
only the few hundred positions that are actually needed are converted: each
instrument position is mapped GRCh38 -> GRCh37 once, the two exposure files are
streamed and filtered to that set, and the outcome is queried directly at the
GRCh38 positions. Cost scales with the number of instruments, not with the size
of the summary statistics.

Allele alignment. ``harmonise`` returns one of four outcomes, plus the refusal
case:

    same         identical effect/other alleles          beta unchanged
    flip         alleles swapped                         beta negated, EAF -> 1 - EAF
    strand       complementary alleles, same order       beta unchanged
    strand_flip  complementary alleles, swapped order    beta negated, EAF -> 1 - EAF
    None         incompatible, or palindromic and
                 ambiguous (see below)                   variant dropped

Palindromic variants (A/T and C/G) are the only case where a strand flip cannot
be distinguished from an allele swap by looking at the alleles. Those are kept
only when the exposure allele frequency is outside the ambiguous band
[mr_palindromic_lo, mr_palindromic_hi]; inside it, the direction of the effect
is unknowable from the summary statistics alone and the variant is dropped
rather than guessed. This is the convention fixed before the analysis
(PLAN_A1 section 4.2).

Row order. Rows are emitted in genomic order (chromosome numerically, then
position) so that the table is byte-reproducible across runs. The source script
iterated over a Python set and therefore produced an arbitrary order; the set of
rows, the columns and every value are unchanged.

BMI frequency. The BMI record stores the frequency exactly as it appears in the
BMI file, while the fasting insulin record stores the realigned frequency. Only
beta and se enter any estimate, and BMI's frequency is unused downstream, so
this asymmetry is preserved rather than silently "fixed".
"""
from __future__ import annotations

import gzip
import importlib
import pathlib
import sys
from typing import Container, Dict, List, Optional, Tuple

import pandas as pd

# Repository root, so that ``config`` resolves when the file is run directly.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

from config.config import (  # noqa: E402  (import after sys.path edit)
    PATHS,
    PARAMS,
    SEED,
    CHAIN_HG38_TO_HG19,
    GWAS_PRIMARY_FILE,
)

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# The instrument tables are produced by the previous stage; the exposure file
# names and the chain reader live there, so they are imported rather than
# repeated. The stage guards its driver with ``if __name__ == "__main__"``, so
# importing it has no side effects.
_select = importlib.import_module("02_select_instruments")
Chain = _select.Chain
FI_EXPOSURE_FILE = _select.FI_EXPOSURE_FILE
BMI_EXPOSURE_FILE = _select.BMI_EXPOSURE_FILE
INSTRUMENTS_FI = _select.INSTRUMENTS_FI
INSTRUMENTS_BMI = _select.INSTRUMENTS_BMI

HARMONISED = "harmonised.tsv"

# Palindromic allele pairs: a strand flip is indistinguishable from an allele
# swap, so they can only be aligned when the frequency is informative.
PALINDROMIC = {frozenset("AT"), frozenset("CG")}
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


# --------------------------------------------------------------------------
# Allele alignment
# --------------------------------------------------------------------------
def harmonise(ea1: str, oa1: str, b1: float, eaf1: float,
              ea2: str, oa2: str) -> Optional[Tuple[float, float, str]]:
    """Realign exposure 1 onto the effect allele of exposure 2.

    Parameters
    ----------
    ea1, oa1 : str
        Effect and other allele of the variant in exposure 1.
    b1 : float
        Effect size in exposure 1.
    eaf1 : float
        Effect-allele frequency in exposure 1; NaN is tolerated.
    ea2, oa2 : str
        Effect and other allele of the *reference* trait (the outcome here).

    Returns
    -------
    (beta, eaf, label) for an alignable variant, or None when the alleles are
    incompatible or the variant is palindromic with an ambiguous frequency.

    Notes
    -----
    The frequency is complemented on every swap, because the effect allele
    changed identity: freq(other) = 1 - freq(effect). The guard
    ``eaf1 == eaf1`` is True for every value except NaN, so a missing frequency
    is carried through as NaN instead of being turned into 1 - NaN by arithmetic.
    """
    if ea1 == ea2 and oa1 == oa2:
        return b1, eaf1, "same"

    if ea1 == oa2 and oa1 == ea2:
        # Alleles swapped: negate the effect and complement the frequency.
        return -b1, (1 - eaf1 if eaf1 == eaf1 else eaf1), "flip"

    if COMPLEMENT.get(ea1) == ea2 and COMPLEMENT.get(oa1) == oa2:
        # The two files are on opposite strands: same alleles after
        # complementation, so the effect keeps its sign.
        if frozenset((ea1, oa1)) in PALINDROMIC and \
                PARAMS.mr_palindromic_lo <= eaf1 <= PARAMS.mr_palindromic_hi:
            return None
        return b1, eaf1, "strand"

    if COMPLEMENT.get(ea1) == oa2 and COMPLEMENT.get(oa1) == ea2:
        # Opposite strands and swapped alleles: negate and complement.
        if frozenset((ea1, oa1)) in PALINDROMIC and \
                PARAMS.mr_palindromic_lo <= eaf1 <= PARAMS.mr_palindromic_hi:
            return None
        return -b1, (1 - eaf1 if eaf1 == eaf1 else eaf1), "strand_flip"

    return None     # incompatible allele pair (for example a triallelic site)


# --------------------------------------------------------------------------
# Streamed readers, each restricted to the positions actually needed
# --------------------------------------------------------------------------
def read_fasting_insulin(path: pathlib.Path,
                         needed: Container[Tuple[str, int]]
                         ) -> Dict[Tuple[str, int], Tuple[str, str, float, float, float, float]]:
    """Collect fasting insulin effects at the required GRCh37 positions.

    Parameters
    ----------
    path : GRCh37 summary statistics, gzip TSV.
    needed : collection of GRCh37 (chrom, pos) pairs to retain; the mapping
        built by the driver also carries the GRCh38 partner, which is unused here.

    Returns
    -------
    dict keyed by GRCh37 position, holding (EA, OA, beta, se, eaf, N).

    Notes
    -----
    The whole file is streamed because the summary statistics are not indexed;
    the position test happens before any parsing, so the cost per irrelevant
    row is one integer conversion and a hash lookup.
    """
    found: Dict[Tuple[str, int], Tuple[str, str, float, float, float, float]] = {}
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {name: i for i, name in enumerate(header)}
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            key = (cols[ix["chromosome"]], int(cols[ix["base_pair_location"]]))
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
                              float(cols[ix["effect_allele_frequency"]]),
                              float(cols[ix["sample_size"]]))
            except ValueError:
                continue
    return found


def read_bmi(path: pathlib.Path,
             needed: Container[Tuple[str, int]]
             ) -> Dict[Tuple[str, int], Tuple[str, str, float, float, float]]:
    """Collect BMI effects at the required GRCh37 positions.

    Parameters
    ----------
    path : GRCh37 summary statistics, gzip TSV.
    needed : collection of GRCh37 (chrom, pos) pairs to retain.

    Returns
    -------
    dict keyed by GRCh37 position, holding (EA, OA, beta, se, eaf).

    Notes
    -----
    No sample size is carried: the BMI file has none per variant, and BMI
    enters as a second exposure rather than as an instrument source here.
    """
    found: Dict[Tuple[str, int], Tuple[str, str, float, float, float]] = {}
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {name: i for i, name in enumerate(header)}
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            key = (cols[ix["CHR"]], int(cols[ix["POS"]]))
            if key not in needed:
                continue
            ea = cols[ix["Tested_Allele"]].upper()
            oa = cols[ix["Other_Allele"]].upper()
            if len(ea) != 1 or len(oa) != 1:
                continue
            try:
                found[key] = (ea, oa,
                              float(cols[ix["BETA"]]),
                              float(cols[ix["SE"]]),
                              float(cols[ix["Freq_Tested_Allele_in_HRS"]]))
            except ValueError:
                continue
    return found


def read_outcome(path: pathlib.Path,
                 needed: Container[Tuple[str, int]]
                 ) -> Dict[Tuple[str, int], Tuple[str, str, float, float, float, str]]:
    """Collect outcome effects at the required GRCh38 positions.

    Parameters
    ----------
    path : GRCh38 summary statistics, gzip TSV.
    needed : collection of GRCh38 (chrom, pos) pairs to retain.

    Returns
    -------
    dict keyed by GRCh38 position, holding (EA, OA, beta, se, eaf, rsid).

    Notes
    -----
    The primary outcome file labels chromosomes without the "chr" prefix, which
    is why no stripping is done here; the secondary outcome in
    06_second_outcome.py does carry the prefix and strips it.
    """
    found: Dict[Tuple[str, int], Tuple[str, str, float, float, float, str]] = {}
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {name: i for i, name in enumerate(header)}
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            key = (cols[ix["chromosome"]], int(cols[ix["base_pair_location"]]))
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
                              float(cols[ix["effect_allele_frequency"]]),
                              cols[ix["rsid"]])
            except ValueError:
                continue
    return found


def _chrom_sort_key(chrom: str) -> Tuple[int, object]:
    """Sort key placing autosomes in numeric order and the rest afterwards."""
    return (0, int(chrom)) if chrom.isdigit() else (1, chrom)


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------
def main() -> None:
    """Build the harmonised table and report the yield of each join."""
    PATHS.ensure()

    fi_ins = pd.read_csv(PATHS.tables / INSTRUMENTS_FI, sep="\t")
    bmi_ins = pd.read_csv(PATHS.tables / INSTRUMENTS_BMI, sep="\t")
    print(f"instruments: fasting insulin {len(fi_ins)}; BMI {len(bmi_ins)}", flush=True)

    fi_pos = {(str(row.chr), int(row.pos)) for row in fi_ins.itertuples()}
    bmi_pos = {(str(row.chr), int(row.pos)) for row in bmi_ins.itertuples()}
    all_pos = fi_pos | bmi_pos
    print(f"positions needed for both exposures (union): {len(all_pos)}", flush=True)

    # ---- 1. GRCh38 -> GRCh37 for the two exposure lookups ---------------
    chain = Chain(PATHS.data / CHAIN_HG38_TO_HG19)
    need19: Dict[Tuple[str, int], Tuple[str, int]] = {}   # GRCh37 -> GRCh38
    rev19: Dict[Tuple[str, int], Tuple[str, int]] = {}    # GRCh38 -> GRCh37
    for chrom, pos in all_pos:
        mapped = chain.lift(chrom, pos)
        if mapped is None:
            continue
        key19 = (str(mapped[0]), int(mapped[1]))
        need19[key19] = (chrom, pos)
        rev19[(chrom, pos)] = key19
    print(f"positions that map back to GRCh37: {len(need19)}", flush=True)

    # ---- 2. effects at those positions ----------------------------------
    fi_rows = read_fasting_insulin(PATHS.data / FI_EXPOSURE_FILE, need19)
    print(f"fasting insulin: {len(fi_rows)} positions matched", flush=True)
    bmi_rows = read_bmi(PATHS.data / BMI_EXPOSURE_FILE, need19)
    print(f"BMI: {len(bmi_rows)} positions matched", flush=True)
    ec_rows = read_outcome(PATHS.data / GWAS_PRIMARY_FILE, all_pos)
    print(f"outcome: {len(ec_rows)} positions matched", flush=True)

    # ---- 3. join, realigning each exposure onto the outcome allele ------
    rows: List[Dict[str, object]] = []
    # Iterating in genomic order keeps the output table reproducible; a Python
    # set alone has no stable order across processes.
    for hg38_key in sorted(all_pos, key=lambda kv: (_chrom_sort_key(kv[0]), kv[1])):
        if hg38_key not in ec_rows:
            continue                       # outcome has no effect: nothing to model
        ec_ea, ec_oa, ec_beta, ec_se, ec_eaf, rsid = ec_rows[hg38_key]

        key19 = rev19.get(hg38_key)
        if key19 is None:
            continue                       # instrument region absent from GRCh37

        record: Dict[str, object] = dict(
            chr=hg38_key[0], pos=hg38_key[1], rsid=rsid,
            ec_EA=ec_ea, ec_OA=ec_oa, ec_beta=ec_beta, ec_se=ec_se, ec_eaf=ec_eaf,
        )

        # Fasting insulin is the exposure of record: a variant that cannot be
        # aligned, or is missing from the file, is dropped from the table.
        usable = True
        if key19 in fi_rows:
            f_ea, f_oa, f_beta, f_se, f_eaf, f_n = fi_rows[key19]
            aligned = harmonise(f_ea, f_oa, f_beta, f_eaf, ec_ea, ec_oa)
            if aligned is None:
                usable = False
            else:
                record.update(fi_beta=aligned[0], fi_se=f_se, fi_eaf=aligned[1],
                              fi_N=f_n, fi_align=aligned[2])
        else:
            usable = False

        # BMI is a second exposure, not an instrument source: a variant that
        # cannot be aligned for BMI is still usable for the univariable fasting
        # insulin model, so its absence does not drop the row.
        if key19 in bmi_rows:
            m_ea, m_oa, m_beta, m_se, m_eaf = bmi_rows[key19]
            aligned = harmonise(m_ea, m_oa, m_beta, m_eaf, ec_ea, ec_oa)
            if aligned is not None:
                record.update(bmi_beta=aligned[0], bmi_se=float(m_se),
                              bmi_eaf=float(m_eaf), bmi_align=aligned[2])

        if usable:
            rows.append(record)

    harmonised = pd.DataFrame(rows)
    harmonised.to_csv(PATHS.tables / HARMONISED, sep="\t", index=False)

    print(f"\nharmonised {len(harmonised)} variants", flush=True)
    if len(harmonised):
        print(f"  with a fasting insulin effect: {harmonised.fi_beta.notna().sum()}", flush=True)
        print(f"  with a BMI effect: {harmonised.bmi_beta.notna().sum()}", flush=True)
        print("  both exposures present (for MVMR): "
              f"{(harmonised.fi_beta.notna() & harmonised.bmi_beta.notna()).sum()}", flush=True)
        print(harmonised[["chr", "pos", "rsid", "ec_beta", "fi_beta", "bmi_beta"]]
              .head(5).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
