#!/usr/bin/env python3
"""
02_select_instruments.py
========================

Purpose
-------
Select genetic instruments for the two exposures of the metabolic arm of the
study, fasting insulin and body mass index (BMI), from their genome-wide
association summary statistics.

The output of this stage is one small table per exposure with one row per
independent variant, in GRCh38 coordinates, holding everything the downstream
stages need: effect allele, other allele, per-allele effect, standard error,
p-value and allele frequency.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
config.PATHS.data / FI_EXPOSURE_FILE       GCST90002238, fasting insulin, GRCh37,
                                           no rsIDs, one row per variant
config.PATHS.data / BMI_EXPOSURE_FILE      GCST006900, BMI (Locke et al. 2018,
                                           GIANT + UK Biobank), GRCh37, with rsIDs
config.PATHS.data / CHAIN_HG19_TO_HG38     UCSC chain file for the coordinate
                                           conversion (exposures are GRCh37, the
                                           outcome of this study is GRCh38)
config.PARAMS.mr_p_thresh                  instrument significance threshold
config.PARAMS.mr_min_n                     minimum per-variant sample size
config.PARAMS.mr_clump_bp                  independence window

Selection rules, in order of application:

    1. p < mr_p_thresh (genome-wide significance);
    2. effect and other allele are single ACGT bases, which drops indels,
       multi-allelic records and malformed rows;
    3. beta, standard error and allele frequency parse as numbers;
    4. fasting insulin only: sample size >= mr_min_n;
    5. the variant lifts over to GRCh38;
    6. distance pruning at mr_clump_bp, p-value priority.

Outputs
-------
config.PATHS.tables / instruments_fi.tsv     chr19, pos19, EA, OA, beta, se, p,
                                             eaf, N, chr, pos
config.PATHS.tables / instruments_bmi.tsv    chr19, pos19, SNP, EA, OA, beta, se,
                                             p, eaf, chr, pos

Both tables keep the GRCh37 coordinates next to the lifted GRCh38 ones, so that
03_harmonise_datasets.py can query the GRCh37 exposure files without repeating
the lift-over.

Notes
-----
Why the sample-size floor for fasting insulin. GCST90002238 is a meta-analysis
of many cohorts and contains a large homogeneous block of rows with N = 7,428
(23.28% of its genome-wide significant rows) contributed by a single
underpowered study. Those rows are not independent of the rest, they carry
inflated effect sizes, and they would swamp a p-value-based instrument set. The
N >= mr_min_n filter isolates and removes that block. BMI is exempt because the
BMI file carries no per-variant sample size: its N is the fixed meta-analysis
sample size of the source study, so there is nothing to filter on.

Why distance pruning and not LD clumping. No 1000 Genomes LD panel is shipped
with this package, so independence is enforced by a distance window
(mr_clump_bp, 10 Mb in the reference run) with p-value priority: variants are
visited from most to least significant and a variant is kept only when no
already-kept variant on the same chromosome lies within the window. This is a
weaker criterion than LD clumping and is stated as such in the manuscript; the
10 Mb window is deliberately wider than linkage disequilibrium typically
extends, so the instrument set is conservative.

The coordinate conversion is performed by the stage-local ``Chain`` reader
documented below. The algorithm is identical to
``code/02_colocalisation/01_chain_liftover.py``; it is duplicated here so that
this stage runs on its own, because the original source imported the external
``liftover`` package, which is not installed in this environment.
"""
from __future__ import annotations

import gzip
import pathlib
import sys
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

# Repository root, so that ``config`` resolves when the file is run directly.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

from config.config import (  # noqa: E402  (import after sys.path edit)
    PATHS,
    PARAMS,
    SEED,
    CHAIN_HG19_TO_HG38,
)

# Exposure summary statistics. Declared here because config.py, which owns every
# path and accession, does not yet carry the two exposure accessions; the file
# names are resolved against config.PATHS.data and are never absolute.
FI_EXPOSURE_FILE = "GCST90002238.tsv.gz"
BMI_EXPOSURE_FILE = "GCST006900_bmi.txt.gz"

# Output table names, unchanged from the source analysis so that downstream
# stages and shipped results keep referring to the same files.
INSTRUMENTS_FI = "instruments_fi.tsv"
INSTRUMENTS_BMI = "instruments_bmi.tsv"

# Single-base substitutions only; anything else is dropped before lift-over.
BASES = frozenset("ACGT")


# --------------------------------------------------------------------------
# Minimal UCSC chain reader
# --------------------------------------------------------------------------
class Chain:
    """Forward-block lookup over one UCSC chain file.

    Only whole-block lookups are implemented, which is sufficient for
    single-base substitutions: a variant whose position falls inside an aligned
    block receives the block offset, and a variant outside every block is
    reported as unmappable rather than being placed by interpolation.

    Attributes
    ----------
    blocks : dict
        source chromosome -> sorted list of block tuples
        (t_start, t_end, target_chrom, q_start, q_strand, q_size, block_size),
        all coordinates 0-based half-open.
    """

    def __init__(self, path: "str | pathlib.Path"):
        self.blocks: Dict[str, List[Tuple[int, int, str, int, str, int, int]]] = defaultdict(list)
        self._parse(str(path))
        for chrom in self.blocks:
            self.blocks[chrom].sort()

    def _parse(self, path: str) -> None:
        """Read the chain file into ``self.blocks``.

        A chain begins with a header line (``chain score tName tSize tStrand
        tStart tEnd qName qSize qStrand qStart qEnd id``) followed by block
        lines of the form ``size dt dq``. The running target and query offsets
        advance by the block size plus the corresponding gap; a two-field line
        marks the last block of the chain.
        """
        with gzip.open(path, "rt") as fh:
            in_chain = False
            t_name = q_name = ""
            q_strand = "+"
            q_size = 0
            t_pos = q_pos = 0
            for line in fh:
                if line.startswith("chain"):
                    p = line.split()
                    t_name = p[2].replace("chr", "")      # target chromosome
                    q_name = p[7].replace("chr", "")      # query chromosome
                    q_size = int(p[8])
                    q_strand = p[9]
                    t_pos, q_pos = int(p[5]), int(p[10])
                    in_chain = True
                elif line.startswith("#") or not line.strip():
                    continue
                elif in_chain:
                    p = line.split()
                    size = int(p[0])
                    self.blocks[t_name].append(
                        (t_pos, t_pos + size, q_name, q_pos, q_strand, q_size, size)
                    )
                    if len(p) >= 3:
                        t_pos += size + int(p[1])
                        q_pos += size + int(p[2])
                    else:
                        in_chain = False      # last block of this chain

    def lift(self, chrom: str, pos: int) -> Optional[Tuple[str, int]]:
        """Map a 1-based source position to the target assembly.

        Parameters
        ----------
        chrom : str
            Source chromosome, with or without the "chr" prefix.
        pos : int
            1-based source position.

        Returns
        -------
        (target_chrom, target_pos), 1-based, or None when the position is not
        covered by any block of the chain.

        Notes
        -----
        The block is located by binary search on block starts (the blocks are
        sorted by construction), then the within-block offset is applied. For a
        reverse-strand block the chain stores a 0-based coordinate counted from
        the far end of the query chromosome, so the forward 1-based coordinate
        is q_size - q_start - offset.
        """
        chrom = str(chrom).replace("chr", "")
        pos0 = pos - 1
        blocks = self.blocks.get(chrom)
        if not blocks:
            return None

        # Index of the last block whose start is <= pos0.
        lo, hi = 0, len(blocks) - 1
        while lo <= hi:
            mid = (lo + hi) // 2
            if blocks[mid][0] <= pos0:
                lo = mid + 1
            else:
                hi = mid - 1
        if hi < 0:
            return None

        t0, t1, q_name, q0, q_strand, q_size, _ = blocks[hi]
        if not (t0 <= pos0 < t1):
            return None                 # falls in an unaligned gap
        offset = pos0 - t0
        if q_strand == "+":
            return q_name, q0 + offset + 1
        return q_name, q_size - q0 - offset


# --------------------------------------------------------------------------
# Readers
# --------------------------------------------------------------------------
def load_fasting_insulin(path: pathlib.Path) -> pd.DataFrame:
    """Stream the fasting insulin summary statistics and apply rules 1 to 4.

    Parameters
    ----------
    path : pathlib.Path
        GCST90002238 gzip TSV.

    Returns
    -------
    DataFrame with columns chr19, pos19, EA, OA, beta, se, p, eaf, N.

    Notes
    -----
    The file is streamed line by line rather than read into pandas: the full
    table is several gigabytes, and every filter is applied before a row is
    materialised. Missing or unparsable numeric fields drop the row instead of
    raising, because the source file encodes missing values as empty strings.
    """
    rows: List[Tuple[object, ...]] = []
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {name: i for i, name in enumerate(header)}
        for line in fh:
            cols = line.rstrip("\n").split("\t")

            # Rule 1: genome-wide significance.
            try:
                p = float(cols[ix["p_value"]])
            except ValueError:
                continue
            if p >= PARAMS.mr_p_thresh:
                continue

            # Rule 2: single-base substitution.
            ea = cols[ix["effect_allele"]].upper()
            oa = cols[ix["other_allele"]].upper()
            if len(ea) != 1 or len(oa) != 1 or ea not in BASES or oa not in BASES:
                continue

            # Rule 3: numeric effect, standard error, frequency, sample size.
            try:
                n = float(cols[ix["sample_size"]])
                beta = float(cols[ix["beta"]])
                se = float(cols[ix["standard_error"]])
                eaf = float(cols[ix["effect_allele_frequency"]])
            except ValueError:
                continue

            # Rule 4: sample-size floor, which excludes the N = 7,428 block.
            if n < PARAMS.mr_min_n:
                continue

            rows.append((cols[ix["chromosome"]], int(cols[ix["base_pair_location"]]),
                         ea, oa, beta, se, p, eaf, n))
    return pd.DataFrame(rows, columns=["chr19", "pos19", "EA", "OA",
                                       "beta", "se", "p", "eaf", "N"])


def load_bmi(path: pathlib.Path) -> pd.DataFrame:
    """Stream the BMI summary statistics and apply rules 1 to 3.

    Parameters
    ----------
    path : pathlib.Path
        GCST006900 BMI gzip TSV.

    Returns
    -------
    DataFrame with columns chr19, pos19, SNP, EA, OA, beta, se, p, eaf.

    Notes
    -----
    No sample-size filter: this file has no per-variant N column, so rule 4 does
    not apply (see the module Notes). The rsID is kept because it is the only
    human-readable identifier the exposure file provides.
    """
    rows: List[Tuple[object, ...]] = []
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {name: i for i, name in enumerate(header)}
        for line in fh:
            cols = line.rstrip("\n").split("\t")

            try:
                p = float(cols[ix["P"]])
            except ValueError:
                continue
            if p >= PARAMS.mr_p_thresh:
                continue

            ea = cols[ix["Tested_Allele"]].upper()
            oa = cols[ix["Other_Allele"]].upper()
            if len(ea) != 1 or len(oa) != 1 or ea not in BASES or oa not in BASES:
                continue

            try:
                beta = float(cols[ix["BETA"]])
                se = float(cols[ix["SE"]])
                # Frequencies are taken from the HRS subsample, the only
                # frequency column the release provides.
                eaf = float(cols[ix["Freq_Tested_Allele_in_HRS"]])
            except ValueError:
                continue

            rows.append((cols[ix["CHR"]], int(cols[ix["POS"]]), cols[ix["SNP"]],
                         ea, oa, beta, se, p, eaf))
    return pd.DataFrame(rows, columns=["chr19", "pos19", "SNP", "EA", "OA",
                                       "beta", "se", "p", "eaf"])


# --------------------------------------------------------------------------
# Coordinate conversion and pruning
# --------------------------------------------------------------------------
def lift_frame(frame: pd.DataFrame, chain: Chain,
               chrom_col: str = "chr19", pos_col: str = "pos19") -> pd.DataFrame:
    """Lift every row of ``frame`` from GRCh37 to GRCh38, dropping misses.

    Parameters
    ----------
    frame : DataFrame
        Instrument candidates with GRCh37 coordinates.
    chain : Chain
        Chain reader opened on the hg19ToHg38 file.
    chrom_col, pos_col : str
        Names of the coordinate columns to read.

    Returns
    -------
    DataFrame with all original columns plus ``chr`` and ``pos`` (GRCh38).
    Rows that do not lift are dropped: an unmappable variant cannot be matched
    to the GRCh38 outcome, so keeping it would only defer the failure.

    Notes
    -----
    A dict is built per row rather than assigning into the frame, because the
    surviving row count is unknown until the loop finishes; building a list of
    dicts and constructing the frame once avoids repeated reallocation.
    """
    lifted: List[Dict[str, object]] = []
    for rec in frame.itertuples(index=False):
        mapped = chain.lift(str(getattr(rec, chrom_col)), int(getattr(rec, pos_col)))
        if mapped is None:
            continue
        row = rec._asdict()
        row["chr"] = str(mapped[0])
        row["pos"] = int(mapped[1])
        lifted.append(row)
    return pd.DataFrame(lifted)


def prune_by_distance(frame: pd.DataFrame,
                      window: int = PARAMS.mr_clump_bp) -> pd.DataFrame:
    """Keep one variant per genomic neighbourhood, p-value priority.

    Parameters
    ----------
    frame : DataFrame
        Lifted instrument candidates, needs ``chr``, ``pos`` and ``p``.
    window : int
        Minimum distance in base pairs between two kept variants.

    Returns
    -------
    DataFrame with the retained rows, in order of increasing p.

    Algorithm
    ---------
    Sort by p ascending, then walk the rows: a candidate is kept when no
    already-kept variant on the same chromosome lies closer than ``window``.
    The result is a greedy maximal independent set in which the strongest
    association always wins a contested window.

    Notes
    -----
    This is O(n * kept) and therefore quadratic in the worst case, which is
    acceptable because the candidate set is a few thousand genome-wide
    significant variants, not the whole summary statistic.
    """
    ordered = frame.sort_values("p").reset_index(drop=True)
    kept: List[Tuple[str, int, Any]] = []
    for rec in ordered.itertuples(index=False):
        chrom = str(rec.chr)
        pos = int(rec.pos)
        in_window = False
        for kept_chrom, kept_pos, _ in kept:
            if kept_chrom == chrom and abs(pos - kept_pos) < window:
                in_window = True
                break
        if not in_window:
            kept.append((chrom, pos, rec))
    return pd.DataFrame([row._asdict() for _, _, row in kept])


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------
def main() -> None:
    """Run the selection for both exposures and write the instrument tables."""
    PATHS.ensure()
    chain_path = PATHS.data / CHAIN_HG19_TO_HG38
    if not chain_path.exists():
        raise SystemExit(f"missing chain file: {chain_path}")
    chain = Chain(chain_path)

    fi_raw = load_fasting_insulin(PATHS.data / FI_EXPOSURE_FILE)
    print(f"fasting insulin: {len(fi_raw):,} single-base variants at "
          f"p<{PARAMS.mr_p_thresh:g} with N>={PARAMS.mr_min_n:,}", flush=True)
    fi = prune_by_distance(lift_frame(fi_raw, chain))
    print(f"  {len(fi)} independent instruments after lift-over and "
          f"{PARAMS.mr_clump_bp / 1e6:g} Mb pruning", flush=True)
    fi.to_csv(PATHS.tables / INSTRUMENTS_FI, sep="\t", index=False)

    bmi_raw = load_bmi(PATHS.data / BMI_EXPOSURE_FILE)
    print(f"BMI: {len(bmi_raw):,} single-base variants at "
          f"p<{PARAMS.mr_p_thresh:g}", flush=True)
    bmi = prune_by_distance(lift_frame(bmi_raw, chain))
    print(f"  {len(bmi)} independent instruments after lift-over and "
          f"{PARAMS.mr_clump_bp / 1e6:g} Mb pruning", flush=True)
    bmi.to_csv(PATHS.tables / INSTRUMENTS_BMI, sep="\t", index=False)

    # First five instruments of each set, so a run can be eyeballed in the log.
    for name, table in (("fasting insulin", fi), ("BMI", bmi)):
        print(f"\n{name}: first 5 instruments")
        print(table[["chr", "pos", "EA", "OA", "beta", "se", "p"]]
              .head(5).to_string(index=False))


if __name__ == "__main__":
    main()
