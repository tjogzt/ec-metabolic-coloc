"""
11_palindromic_audit.py
=======================

Purpose
-------
Audit the allele alignment. A palindromic variant (A/T, T/A, C/G, G/C) has the
same alleles on both strands, so a strand flip cannot be detected from the
allele strings alone; if such variants were flipped by mistake their effect
directions would be inverted with no error raised anywhere downstream.

The audit rests on one expectation: if both datasets are on the reference
forward strand, the correlation between the eQTL z and the outcome z inside a
window must be the same for palindromic and non-palindromic variants. A
systematic sign inversion of the palindromic subset drives its correlation
towards zero or below zero while leaving the non-palindromic subset intact.
Restricting further to palindromic variants with a minor allele frequency
between 0.05 and 0.45 removes the ones that are trivially easy to align (very
rare or very common), so the comparison is made where a flip would actually
matter.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
--tissue        GTEx tissue label (default: Liver)
--gene          Ensembl gene id (default: ENSG00000141293, SKAP1)
--chrom, --lead window centre; the gene's eQTL were extracted around this lead
                (default: chr17, 48,127,474)
--pad           window half-width (default: config.PARAMS.region_pad_bp)
--gwas          outcome file in data/ (default: the primary outcome)
--outcome-build build of that outcome, hg38 (default) or hg19. When hg19, every
                eQTL position is carried back with the hg38 -> hg19 chain before
                the join, which is the configuration the first release used with
                the second outcome.
--subset        pre-extracted eQTL TSV to use instead of a tabix query; the
                locus subsets written by 07 are accepted as they are
--out           optional TSV with the per-variant z statistics

Outputs
-------
A table on stdout: group, n, Spearman rho, p, and the fraction of variants with
concordant z signs. Nothing is written unless --out is given.

Notes
-----
The first release of this audit read the liver locus subset, lifted it to
GRCh37 and matched the second outcome. It printed the same four groups, but the
subset had been extracted for the primary outcome, so the window and the
outcome did not belong to the same analysis. The defaults here audit the primary
analysis; --gwas with --outcome-build hg19 reproduces the earlier run.

The frequency band and the palindromic classes are properties of the diagnostic,
not analysis parameters, so they are module constants rather than config entries.
"""
from __future__ import annotations

import argparse
import pathlib
import sys
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# Local package imports; the directory holding this file is on sys.path.
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))       # repository root, for config/

from config.config import (PATHS, PARAMS, TISSUES, GWAS_PRIMARY_FILE,
                           GWAS_SECONDARY_FILE, CHAIN_HG38_TO_HG19)
from importlib import import_module

_gwas = import_module("02_gwas_reader")
_eqtl = import_module("03_eqtl_catalogue_reader")
_lift = import_module("01_chain_liftover")

# Allele pairs that are unchanged by complementing both members.
PALINDROMIC = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}

# Frequency band in which a mistaken flip would change the correlation. Outside
# it one allele is so rare that alignment is unambiguous from frequency alone.
AMBIGUOUS_MAF_LO = 0.05
AMBIGUOUS_MAF_HI = 0.45

# Minimum variants for a group to be reported rather than counted.
MIN_GROUP = 5


def read_window(tissue: str, gene_id: str, chrom: str, lead: int, pad: int,
                subset: Optional[str]) -> pd.DataFrame:
    """Read the eQTL variants of one gene in one window.

    Parameters
    ----------
    tissue : GTEx tissue label, resolved through config.TISSUES.
    gene_id : Ensembl id without version suffix.
    chrom : chromosome of the window.
    lead : window centre, 1-based.
    pad : window half-width in base pairs.
    subset : pre-extracted TSV, or None to query the local eQTL file.

    Returns
    -------
    DataFrame with the columns of 03_eqtl_catalogue_reader.OUT_COLS, filtered to
    finite beta and se.

    Raises
    ------
    SystemExit
        When the tissue or the local eQTL file is unavailable. An empty return
        value would be indistinguishable from "this gene has no eQTL", which is
        precisely the ambiguity the audit must not create.
    """
    if subset:
        path = pathlib.Path(subset)
        if not path.is_absolute():
            path = PATHS.cache / path
        if not path.exists():
            raise SystemExit(f"missing subset: {path}")
        df = pd.read_csv(path, sep="\t")
    else:
        qtd = TISSUES.get(tissue)
        if qtd is None:
            raise SystemExit(f"unknown tissue {tissue}; add it to config.TISSUES")
        local = PATHS.data / f"{qtd}.all.tsv.gz"
        if not local.exists():
            raise SystemExit(f"missing eQTL file: {local.name}")
        df = _eqtl.qtl_region(chrom, max(1, lead - pad), lead + pad,
                              str(local), _eqtl.column_positions(local))

    df = df[df["gene_id"] == gene_id].dropna(subset=["position", "beta", "se"])
    df = df[(df["maf"] >= PARAMS.eqtl_maf_min)
            & (df["maf"] <= PARAMS.eqtl_maf_max) & (df["se"] > 0)]
    if df.empty:
        raise SystemExit(f"no usable eQTL variant for {gene_id} in "
                         f"{chrom}:{lead - pad}-{lead + pad}")
    return df


def to_outcome_build(df: pd.DataFrame, chrom: str,
                     build: str) -> pd.DataFrame:
    """Express the eQTL positions on the outcome build.

    Parameters
    ----------
    df : eQTL rows, GRCh38.
    chrom : chromosome of the rows.
    build : "hg38" to keep the coordinates, "hg19" to lift them back.

    Returns
    -------
    The frame with a ``pos_out`` column; variants that do not lift are dropped.

    Notes
    -----
    A lift that silently returned the input would hide build errors, so a
    missing chain file stops the run instead of falling back to hg38.
    """
    df = df.copy()
    if build == "hg38":
        df["pos_out"] = df["position"].astype(int)
        return df

    chain_path = PATHS.data / CHAIN_HG38_TO_HG19
    if not chain_path.exists():
        raise SystemExit(f"missing chain file: {chain_path.name}")
    down = _lift.Chain(str(chain_path))

    lifted = np.full(len(df), np.nan)
    for j, pos in enumerate(df["position"].to_numpy()):
        mapped = down.lift(chrom, int(pos))
        if mapped is not None:
            lifted[j] = mapped[1]
    df["pos_out"] = lifted
    df = df.dropna(subset=["pos_out"])
    if df.empty:
        raise SystemExit("no eQTL variant lifts to the outcome build; "
                         "check --outcome-build")
    df["pos_out"] = df["pos_out"].astype(int)
    return df


def pair_with_outcome(df: pd.DataFrame, store: Dict[str, pd.DataFrame],
                      chrom: str) -> pd.DataFrame:
    """Join the eQTL rows to the outcome on position and alleles.

    Parameters
    ----------
    df : eQTL rows carrying ``pos_out``.
    store : per-chromosome outcome table from 02_gwas_reader.load_gwas.
    chrom : chromosome of the rows.

    Returns
    -------
    DataFrame with ref, alt, maf, z_e, z_g, pos.

    Notes
    -----
    The z statistics are computed here, not taken from the files: the outcome
    file's p-value and the eQTL file's p-value are on different scales, while
    beta/se is directly comparable once the alleles are oriented. A negative
    beta always refers to the alt allele, in both datasets.
    """
    rows: List[tuple] = []
    for x in df.itertuples(index=False):
        ref = str(x.ref).upper()
        alt = str(x.alt).upper()
        hit = _gwas.lookup(store, chrom, int(x.pos_out), ref, alt)
        if hit is None:
            continue
        beta_g, se_g, _p = hit
        if not (np.isfinite(beta_g) and np.isfinite(se_g)) or se_g <= 0:
            continue
        rows.append((ref, alt, float(x.maf), float(x.beta) / float(x.se),
                     float(beta_g) / float(se_g), int(x.pos_out)))
    return pd.DataFrame(rows, columns=["ref", "alt", "maf", "z_e", "z_g", "pos"])


def summarise(d: pd.DataFrame) -> pd.DataFrame:
    """Correlation of the two z statistics within each allele class.

    Parameters
    ----------
    d : output of :func:`pair_with_outcome`.

    Returns
    -------
    DataFrame with label, n, rho, p, and the fraction of variants with the same
    z sign, one row per class plus the whole window.

    Notes
    -----
    The fourth class is the informative one: if palindromic variants were being
    flipped, the 0.05-0.45 frequency band would show it, because the other two
    classes are dominated by variants whose orientation is unambiguous.
    """
    palindromic = np.array([(r, a) in PALINDROMIC
                            for r, a in zip(d["ref"], d["alt"])])
    ambiguous = palindromic & d["maf"].between(AMBIGUOUS_MAF_LO,
                                               AMBIGUOUS_MAF_HI).to_numpy()

    classes = [
        ("all", np.ones(len(d), dtype=bool)),
        ("non-palindromic", ~palindromic),
        ("palindromic, all", palindromic),
        (f"palindromic, MAF {AMBIGUOUS_MAF_LO}-{AMBIGUOUS_MAF_HI}", ambiguous),
    ]

    rows: List[Dict[str, object]] = []
    for label, mask in classes:
        sub = d[mask]
        if len(sub) < MIN_GROUP:
            rows.append({"label": label, "n": len(sub), "rho": np.nan,
                         "p": np.nan, "same_sign": np.nan})
            continue
        rho, p = spearmanr(sub["z_e"], sub["z_g"])
        same = float((np.sign(sub["z_e"]) == np.sign(sub["z_g"])).mean())
        rows.append({"label": label, "n": len(sub), "rho": float(rho),
                     "p": float(p), "same_sign": same})
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tissue", default="Liver")
    ap.add_argument("--gene", default="ENSG00000141293", help="Ensembl id")
    ap.add_argument("--chrom", default="17")
    ap.add_argument("--lead", type=int, default=48_127_474,
                    help="window centre, 1-based")
    ap.add_argument("--pad", type=int, default=PARAMS.region_pad_bp)
    ap.add_argument("--gwas", default=GWAS_PRIMARY_FILE)
    ap.add_argument("--outcome-build", default="hg38", choices=["hg38", "hg19"])
    ap.add_argument("--subset", default=None,
                    help="pre-extracted eQTL TSV instead of a tabix query")
    ap.add_argument("--out", default=None,
                    help="optional TSV with the per-variant statistics")
    args = ap.parse_args()

    PATHS.ensure()

    gwas_path = pathlib.Path(args.gwas)
    if not gwas_path.is_absolute():
        gwas_path = PATHS.data / gwas_path
    if not gwas_path.exists():
        raise SystemExit(f"missing outcome file: {gwas_path}")

    qtl = read_window(args.tissue, args.gene, args.chrom, args.lead, args.pad,
                      args.subset)
    qtl = to_outcome_build(qtl, args.chrom, args.outcome_build)

    store = _gwas.load_gwas(gwas_path)
    d = pair_with_outcome(qtl, store, args.chrom)
    if len(d) < MIN_GROUP:
        raise SystemExit(f"only {len(d)} variants paired; the window, the gene "
                         f"or the outcome build is wrong")

    print(f"{args.gene} in {args.tissue}: window "
          f"{args.chrom}:{args.lead - args.pad:,}-{args.lead + args.pad:,}, "
          f"outcome build {args.outcome_build}", flush=True)
    print(f"paired variants: {len(d):,}", flush=True)

    table = summarise(d)
    for row in table.itertuples(index=False):
        if np.isnan(row.rho):
            print(f"  {row.label:34s} n={row.n:5d}  too few variants", flush=True)
        else:
            print(f"  {row.label:34s} n={row.n:5d}  rho={row.rho:+.3f} "
                  f"(p={row.p:.1e})  same sign={100 * row.same_sign:.1f}%",
                  flush=True)

    # A collapsed correlation in the ambiguous class, with the non-palindromic
    # class intact, is the signature of a systematic orientation error.
    base = table[table["label"] == "non-palindromic"]["rho"]
    amb = table[table["label"].str.startswith("palindromic, MAF")]["rho"]
    if len(base) and len(amb) and np.isfinite(base.iloc[0]) and np.isfinite(amb.iloc[0]):
        delta = float(base.iloc[0]) - float(amb.iloc[0])
        verdict = "consistent" if abs(delta) < 0.15 else "CHECK ORIENTATION"
        print(f"  non-palindromic minus ambiguous rho: {delta:+.3f} -> {verdict}",
              flush=True)

    if args.out:
        dest = pathlib.Path(args.out)
        if not dest.is_absolute():
            dest = PATHS.tables / dest
        d.to_csv(dest, sep="\t", index=False)
        print(f"wrote {dest.name}", flush=True)


if __name__ == "__main__":
    main()
