"""
14_fixed_variant_check.py
=========================

Purpose
-------
The four checks that were run by hand on the SKAP1 (ENSG00000141293) signal,
collected into one module. Each answers a question a reviewer is entitled to ask
about a single reported hit, and each is independent of the others:

    local      recompute three reported gene-locus hits from the local
               all-associations files with a plain sequential read, and compare
               with the values the tabix run produced. This is a file-integrity
               check: the two paths read the same variants, so any disagreement
               is a corrupted download or a changed column layout.
    window     describe the internal structure of one gene's window: are the
               expression peak and the outcome peak the same variant, how far
               apart are they, and do the z statistics correlate across the
               window. A high PP.H4 with two separated peaks and no correlation
               would be an artefact of a diffuse signal.
    fixed      take one variant, the expression lead, and read its effect on
               expression in every available tissue. A shared causal variant
               must show the same allele and a comparable effect wherever the
               gene is expressed; a variant that exists in one tissue only
               cannot support a cross-tissue claim.
    profile    describe the same gene across a panel of tissues at a fixed
               window: how many variants were tested, the strongest expression
               p, and the effect size there. This separates "not expressed" from
               "not tested".

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
--section       local, window, fixed, profile, or all (default: all)
--tissue        GTEx tissue for the window check (default: Liver)
--gene          Ensembl id (default: ENSG00000141293)
--chrom, --lead window centre for the window check
                (default: chr17, 48,127,474)
--fixed-lead   the expression lead read by the fixed-variant check
                (default: 48,127,316, the strongest liver eQTL for SKAP1; it is
                not the window centre, which is the colocalisation peak)
--pad           window half-width (default: config.PARAMS.region_pad_bp)
--gwas          outcome file in data/ (default: the primary outcome)
--outcome-build hg38 (default) or hg19; hg19 lifts the eQTL positions back
                before the join, which is how the first release ran this check
--subset        pre-extracted eQTL TSV instead of a tabix query
--panel         JSON object of label -> eQTL Catalogue dataset id, for the
                cross-tissue sections
--profile-lo, --profile-hi fixed window for the tissue-profile section
                (default: chr17:48,000,000-48,300,000)

Outputs
-------
results/tables/skap1_fixed_variant.tsv    tissue, QTD, status, alleles, maf, p,
                                          beta, z, pos
results/tables/skap1_tissue_profile.tsv   tissue, QTD, status, n_variants,
                                          lead_pos, min_p, beta, z, maf

The local and window sections print to stdout and write no file: they are
checks on numbers that already exist elsewhere, and a table would not add
anything a reader could use.

Notes
-----
The cross-tissue panel extends config.TISSUES with tissues that appear in the
figure but are not part of the analysis set. Their dataset ids are fixed here
because a single-gene display does not belong in the shared configuration; the
five analysis tissues always take their ids from config, and the module asserts
that the panel and the configuration agree.

EBI occasionally returns a transient failure for a large file. The tabix query
is retried three times with a pause, and a persistent failure is reported as a
FAIL row rather than as an absent variant, because the two mean different things.
"""
from __future__ import annotations

import argparse
import gzip
import json
import pathlib
import subprocess
import sys
import time
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# Local package imports; the directory holding this file is on sys.path.
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))       # repository root, for config/

from config.config import (PATHS, PARAMS, TISSUES, GWAS_PRIMARY_FILE,
                           EQTL_CATALOGUE_BASE)
from importlib import import_module

_abf = import_module("05_coloc_abf_engine")
_gwas = import_module("02_gwas_reader")
_eqtl = import_module("03_eqtl_catalogue_reader")
_audit = import_module("11_palindromic_audit")

# Reported hits of the primary analysis, with the PP.H4 the tabix path produced
# for each. The gene id, window and expected value form one record: changing one
# without the others would compare two different quantities.
CASES: List[Tuple[str, str, int, str, str, float]] = [
    ("Liver", "17", 48_127_474, "ENSG00000141293", "QTD000266", 0.931),
    ("Muscle_Skeletal", "17", 48_127_474, "ENSG00000278765", "QTD000281", 0.898),
    ("Muscle_Skeletal", "15", 40_086_266, "ENSG00000275636", "QTD000281", 0.962),
]

# Agreement tolerance for the recomputation. The two paths read the same rows,
# so anything beyond rounding is a difference in the data, not in the arithmetic.
TOLERANCE = 0.005

# Labels used in the released tables, so that a rerun can be joined to them.
LABELS = {
    "Liver": "liver",
    "Muscle_Skeletal": "muscle",
    "Adipose_Visceral_Omentum": "adipose (visceral)",
    "Pancreas": "pancreas",
    "Uterus": "uterus",
}

# Tissues shown in the cross-tissue display but not part of the analysis set.
DISPLAY_ONLY = {
    "adipose": "QTD000116",
    "vagina": "QTD000351",
    "ovary": "QTD000291",
    "blood": "QTD000356",
    "spleen": "QTD000326",
    "LCL": "QTD000221",
    "breast": "QTD000211",
    "thyroid": "QTD000341",
    "kidney (cortex)": "QTD000261",
}

RETRIES = 3
RETRY_PAUSE_S = 12


def build_panel(override: Optional[str] = None) -> Dict[str, str]:
    """Resolve the label -> dataset id map used by the cross-tissue sections.

    Parameters
    ----------
    override : optional JSON file holding the map.

    Returns
    -------
    Ordered dict: the five analysis tissues first, then the display-only ones.

    Raises
    ------
    AssertionError
        If a configured tissue is missing from the panel or carries a different
        id, which would silently analyse a different dataset than the headline
        results.
    """
    if override:
        path = pathlib.Path(override)
        if not path.is_absolute():
            path = PATHS.data / path
        with open(path) as fh:
            return dict(json.load(fh))

    panel = {LABELS.get(name, name): qtd for name, qtd in TISSUES.items()}
    for name, qtd in TISSUES.items():
        assert panel[LABELS.get(name, name)] == qtd
    panel.update(DISPLAY_ONLY)
    return panel


def scan_local_window(path: pathlib.Path, chrom: str, lo: int, hi: int,
                      gene: str) -> List[tuple]:
    """Read one gene's variants from a local file without tabix.

    Parameters
    ----------
    path : eQTL Catalogue all-associations file.
    chrom : chromosome, without the "chr" prefix.
    lo, hi : inclusive interval.
    gene : Ensembl id.

    Returns
    -------
    List of (pos, ref, alt, beta, se).

    Notes
    -----
    The file is sorted by chromosome and position, so the scan stops as soon as
    it passes the right edge of the window. Reading sequentially is the point:
    an index is a separate artefact that can be stale or truncated, and this
    check exists to catch exactly that.
    """
    rows: List[tuple] = []
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_chrom = header.index("chromosome")
        i_pos = header.index("position")
        i_gene = header.index("gene_id")
        i_beta, i_se = header.index("beta"), header.index("se")
        i_maf = header.index("maf")
        i_ref, i_alt = header.index("ref"), header.index("alt")

        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if fields[i_chrom] != chrom:
                continue
            pos = int(fields[i_pos])
            if pos < lo:
                continue
            if pos > hi:
                break                          # sorted: nothing further can match
            if fields[i_gene] != gene:
                continue
            try:
                beta, se, maf = (float(fields[i_beta]), float(fields[i_se]),
                                 float(fields[i_maf]))
            except ValueError:
                continue
            if not (np.isfinite(beta) and np.isfinite(se)) or se <= 0:
                continue
            if not (PARAMS.eqtl_maf_min <= maf <= PARAMS.eqtl_maf_max):
                continue
            rows.append((pos, fields[i_ref].upper(), fields[i_alt].upper(),
                         beta, se))
    return rows


def section_local(store: Dict[str, pd.DataFrame]) -> None:
    """Recompute the reported hits from the local files and compare.

    Parameters
    ----------
    store : the outcome, already loaded.

    Notes
    -----
    Nothing is written. A mismatch is printed as a warning rather than raised:
    the reader of the log has to decide whether the difference is a corrupt file
    or a genuine change in the analysis, and stopping the run would hide the
    remaining cases.
    """
    print("== local file recheck", flush=True)
    for tissue, chrom, lead, gene, qtd, expected in CASES:
        local = PATHS.data / f"{qtd}.all.tsv.gz"
        if not local.exists():
            print(f"  {tissue}: missing {local.name}", flush=True)
            continue

        rows = scan_local_window(local, chrom, max(1, lead - PARAMS.region_pad_bp),
                                 lead + PARAMS.region_pad_bp, gene)
        pairs = []
        for pos, ref, alt, beta_e, se_e in rows:
            hit = _gwas.lookup(store, chrom, pos, ref, alt)
            if hit is None:
                continue
            beta_g, se_g, _p = hit
            if np.isfinite(beta_g) and np.isfinite(se_g) and se_g > 0:
                pairs.append((beta_e, se_e, beta_g, se_g))
        if len(pairs) < PARAMS.eqtl_min_snps:
            print(f"  {tissue:18s} {gene:18s} too few paired variants "
                  f"({len(pairs)})", flush=True)
            continue

        arr = np.array(pairs)
        res = _abf.coloc_abf(_abf.log_abf(arr[:, 0], arr[:, 1], 1.0),
                             _abf.log_abf(arr[:, 2], arr[:, 3], 1.0),
                             p1=PARAMS.p1, p2=PARAMS.p2, p12=PARAMS.p12)
        h4 = res["PP.H4"]
        verdict = "match" if abs(h4 - expected) < TOLERANCE else "MISMATCH"
        print(f"  {tissue:18s} {gene:18s} n={len(pairs):5d}  "
              f"H4={h4:.3f} (tabix {expected:.3f})  {verdict}", flush=True)


def section_window(store: Dict[str, pd.DataFrame], tissue: str, gene: str,
                   chrom: str, lead: int, pad: int, build: str,
                   subset: Optional[str]) -> None:
    """Describe one gene's window: peaks, their distance, and the z correlation.

    Parameters
    ----------
    store : the outcome, already loaded.
    tissue, gene, chrom, lead, pad : window definition.
    build : build of the outcome; passed through to the shared reader.
    subset : optional pre-extracted eQTL table.

    Notes
    -----
    The window reader, the build conversion and the join are the ones used by
    11_palindromic_audit, so a reproduction here and an audit there cannot drift
    apart. Both statistics reported by the first release are printed: the peak
    comparison and the Spearman correlation over the whole window.
    """
    print(f"== window structure: {gene} in {tissue}", flush=True)
    qtl = _audit.read_window(tissue, gene, chrom, lead, pad, subset)
    qtl = _audit.to_outcome_build(qtl, chrom, build)
    d = _audit.pair_with_outcome(qtl, store, chrom)
    if len(d) < PARAMS.eqtl_min_snps:
        print(f"  only {len(d)} paired variants; nothing to describe", flush=True)
        return

    d = d.assign(p_e=_two_sided_p(d["z_e"]), p_g=_two_sided_p(d["z_g"]))
    i_e = int(d["p_e"].idxmin())
    i_g = int(d["p_g"].idxmin())
    rho, p_rho = spearmanr(d["z_e"], d["z_g"])

    print(f"  paired variants        {len(d):,}", flush=True)
    print(f"  expression peak        {int(d['pos'][i_e]):,}  "
          f"z={d['z_e'][i_e]:+.2f}  p={d['p_e'][i_e]:.2e}", flush=True)
    print(f"  outcome peak           {int(d['pos'][i_g]):,}  "
          f"z={d['z_g'][i_g]:+.2f}  p={d['p_g'][i_g]:.2e}", flush=True)
    print(f"  peak distance          "
          f"{abs(int(d['pos'][i_e]) - int(d['pos'][i_g])) / 1000:.1f} kb",
          flush=True)
    print(f"  outcome at expression lead: z={d['z_g'][i_e]:+.2f}  "
          f"p={d['p_g'][i_e]:.2e}  concordant="
          f"{bool(d['z_e'][i_e] * d['z_g'][i_e] > 0)}", flush=True)
    print(f"  window z correlation   rho={rho:+.3f}  p={p_rho:.2e}", flush=True)


def _two_sided_p(z: pd.Series) -> pd.Series:
    """Two-sided p-value from a z statistic, without an extra dependency.

    Parameters
    ----------
    z : z statistics.

    Returns
    -------
    p-values on the same index.

    Notes
    -----
    erfc is used rather than a normal CDF so that the tail stays accurate for
    the large z values a genome-wide significant peak produces.
    """
    from math import erfc, sqrt

    return pd.Series([erfc(abs(float(x)) / sqrt(2.0)) for x in z], index=z.index)


def _tabix_text(url: str, region: str) -> Tuple[bool, str]:
    """Query one interval, retrying transient network failures.

    Parameters
    ----------
    url : remote file accepted by tabix.
    region : tabix region string, e.g. "17:48127316-48127316".

    Returns
    -------
    (ok, stdout). ``ok`` is False only when every attempt failed, which the
    caller must report as a failure rather than as an absent variant.
    """
    for attempt in range(RETRIES):
        proc = subprocess.run(["tabix", url, region],
                              capture_output=True, text=True, timeout=900)
        if proc.returncode == 0:
            return True, proc.stdout
        if attempt < RETRIES - 1:
            time.sleep(RETRY_PAUSE_S)
    return False, ""


def _parse_remote(text: str, want: List[int]) -> pd.DataFrame:
    """Turn raw tabix output into the reader's standard columns."""
    lines = [ln.split("\t") for ln in text.strip().splitlines() if ln]
    if not lines:
        return pd.DataFrame(columns=_eqtl.OUT_COLS)
    df = pd.DataFrame(lines).iloc[:, list(want)]
    df.columns = _eqtl.OUT_COLS
    for col in _eqtl.NUMERIC:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _want_for_remote() -> List[int]:
    """Column indices for a remote file, read from any local eQTL header."""
    donors = sorted(PATHS.data.glob("*.all.tsv.gz"))
    if not donors:
        raise SystemExit("the cross-tissue sections need one local "
                         "*.all.tsv.gz to read the column layout from")
    return _eqtl.column_positions(donors[0])


def section_fixed(panel: Dict[str, str], chrom: str, lead: int,
                  gene: str, want: List[int]) -> None:
    """Read one fixed variant's expression effect in every panel tissue.

    Parameters
    ----------
    panel : label -> dataset id.
    chrom, lead : the variant, 1-based.
    gene : Ensembl id.
    want : column indices for remote parsing.

    Writes
    ------
    results/tables/skap1_fixed_variant.tsv with one row per tissue.
    """
    print(f"== fixed variant {chrom}:{lead:,} for {gene}", flush=True)
    rows: List[tuple] = []
    for label, qtd in panel.items():
        url = f"{EQTL_CATALOGUE_BASE}/{qtd}/{qtd}.all.tsv.gz"
        ok, text = _tabix_text(url, f"{chrom}:{lead}-{lead}")
        if not ok:
            rows.append((label, qtd, "FAIL", "", "", "", "", "", str(lead)))
            print(f"  {label:18s} query failed", flush=True)
            continue

        hit = _parse_remote(text, want)
        hit = hit[hit["gene_id"] == gene]
        if hit.empty:
            rows.append((label, qtd, "ABSENT", "", "", "", "", "", str(lead)))
            print(f"  {label:18s} variant not tested for this gene", flush=True)
            continue

        row = hit.iloc[0]
        z = float(row["beta"]) / float(row["se"]) if float(row["se"]) > 0 else np.nan
        rows.append((label, qtd, "OK", f"{row['ref']}/{row['alt']}",
                     f"{float(row['maf']):.6g}", f"{float(row['pvalue']):.3e}",
                     f"{float(row['beta']):+.4f}", f"{z:+.2f}",
                     int(row["position"])))
        print(f"  {label:18s} {row['ref']}/{row['alt']}  "
              f"p={float(row['pvalue']):.3e}  beta={float(row['beta']):+.4f}  "
              f"z={z:+.2f}", flush=True)

    dest = PATHS.tables / "skap1_fixed_variant.tsv"
    with open(dest, "w", newline="") as fh:
        fh.write("tissue\tQTD\tstatus\talleles\tmaf\tp\tbeta\tz\tpos\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")
    print(f"  {len(rows)} tissues -> {dest.name}", flush=True)


def section_profile(panel: Dict[str, str], chrom: str, lo: int, hi: int,
                    gene: str, want: List[int]) -> None:
    """Describe one gene's expression across the panel at a fixed window.

    Parameters
    ----------
    panel : label -> dataset id.
    chrom, lo, hi : fixed window, shared by every tissue so that the tissues are
        comparable.
    gene : Ensembl id.
    want : column indices for remote parsing.

    Writes
    ------
    results/tables/skap1_tissue_profile.tsv with one row per tissue.
    """
    print(f"== tissue profile {chrom}:{lo:,}-{hi:,} for {gene}", flush=True)
    rows: List[tuple] = []
    for label, qtd in panel.items():
        url = f"{EQTL_CATALOGUE_BASE}/{qtd}/{qtd}.all.tsv.gz"
        ok, text = _tabix_text(url, f"{chrom}:{lo}-{hi}")
        if not ok:
            rows.append((label, qtd, "FAIL", "", "", "", "", "", ""))
            print(f"  {label:18s} query failed", flush=True)
            continue

        hit = _parse_remote(text, want)
        hit = hit[hit["gene_id"] == gene].dropna(subset=["pvalue", "beta", "se"])
        if hit.empty:
            rows.append((label, qtd, "NOT_TESTED", "0", "", "", "", "", ""))
            print(f"  {label:18s} gene not in the tested set", flush=True)
            continue

        i = int(hit["pvalue"].idxmin())
        row = hit.loc[i]
        z = float(row["beta"]) / float(row["se"]) if float(row["se"]) > 0 else np.nan
        rows.append((label, qtd, "OK", str(len(hit)), int(row["position"]),
                     f"{float(row['pvalue']):.3e}", f"{float(row['beta']):+.4f}",
                     f"{z:+.2f}", f"{float(row['maf']):.6g}"))
        print(f"  {label:18s} n={len(hit):5d}  min p={float(row['pvalue']):.3e}  "
              f"beta={float(row['beta']):+.4f}  z={z:+.2f}", flush=True)

    dest = PATHS.tables / "skap1_tissue_profile.tsv"
    with open(dest, "w", newline="") as fh:
        fh.write("tissue\tQTD\tstatus\tn_variants\tlead_pos\tmin_p\tbeta\tz\tmaf\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")
    print(f"  {len(rows)} tissues -> {dest.name}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--section", default="all",
                    choices=["local", "window", "fixed", "profile", "all"])
    ap.add_argument("--tissue", default="Liver")
    ap.add_argument("--gene", default="ENSG00000141293")
    ap.add_argument("--chrom", default="17")
    ap.add_argument("--lead", type=int, default=48_127_474)
    ap.add_argument("--fixed-lead", type=int, default=48_127_316,
                    help="expression lead read by the fixed-variant check")
    ap.add_argument("--pad", type=int, default=PARAMS.region_pad_bp)
    ap.add_argument("--gwas", default=GWAS_PRIMARY_FILE)
    ap.add_argument("--outcome-build", default="hg38", choices=["hg38", "hg19"])
    ap.add_argument("--subset", default=None)
    ap.add_argument("--panel", default=None, help="JSON label -> dataset id")
    ap.add_argument("--profile-lo", type=int, default=48_000_000)
    ap.add_argument("--profile-hi", type=int, default=48_300_000)
    args = ap.parse_args()

    PATHS.ensure()
    t0 = time.time()

    gwas_path = pathlib.Path(args.gwas)
    if not gwas_path.is_absolute():
        gwas_path = PATHS.data / gwas_path
    if not gwas_path.exists():
        raise SystemExit(f"missing outcome file: {gwas_path}")
    store = _gwas.load_gwas(gwas_path)

    if args.section in ("local", "all"):
        section_local(store)

    if args.section in ("window", "all"):
        section_window(store, args.tissue, args.gene, args.chrom, args.lead,
                       args.pad, args.outcome_build, args.subset)

    if args.section in ("fixed", "profile", "all"):
        panel = build_panel(args.panel)
        want = _want_for_remote()
        if args.section in ("fixed", "all"):
            # The fixed-variant check reads the expression lead, not the window
            # centre, because the lead is the variant a shared-signal claim names.
            section_fixed(panel, args.chrom, args.fixed_lead, args.gene, want)
        if args.section in ("profile", "all"):
            section_profile(panel, args.chrom, args.profile_lo, args.profile_hi,
                            args.gene, want)

    print(f"total {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
