"""
07_run_colocalisation.py
========================

Purpose
-------
Run the ABF colocalisation test between the endometrial cancer risk signal and
cis-eQTLs for every gene that lies inside a risk locus, for one GTEx tissue.

For each locus the script takes every gene with a significant cis-eQTL in that
tissue, restricts the eQTL table to the padded window, joins the two datasets on
position and allele, and tests whether one shared causal variant is more likely
than two distinct ones.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
--tissue      GTEx tissue label; its eQTL dataset id comes from config.TISSUES
              (override with --eqtl-id)
--gwas        outcome file name inside the data directory
              (default: config.GWAS_PRIMARY_FILE)
--remote      query the EBI tabix endpoint instead of a local file; needed for
              datasets that were never mirrored locally
all other thresholds come from config.PARAMS and are not repeated here.

Outputs
-------
results/tables/coloc_<tissue>.tsv with one row per gene-locus pair:

    PP.H0 .. PP.H4   posterior probabilities of the five configurations
    n_snps           variants entering the test
    top_snp_share    posterior mass on the single best shared variant
    tissue, locus    GTEx tissue and "chrom:lead_position"
    gene_id          Ensembl id of the eQTL gene
    best_gwas_p      strongest risk association inside the window

Notes
-----
The risk loci are defined once from the GWAS and are therefore identical across
tissues; only the eQTL side changes. This is what makes the per-tissue counts
comparable in Table 1.

The eQTL file for adipose tissue is not mirrored locally in this package, so
that tissue must be run with --remote. Passing a missing local path without
--remote raises rather than returning zero tests.
"""
from __future__ import annotations

import argparse
import pathlib
import sys
import time
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

# Local package imports; the directory holding this file is on sys.path.
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))       # repository root, for config/

from config.config import PATHS, PARAMS, TISSUES, GWAS_PRIMARY_FILE, EQTL_CATALOGUE_BASE
from importlib import import_module

_abf = import_module("05_coloc_abf_engine")
_gwas = import_module("02_gwas_reader")
_eqtl = import_module("03_eqtl_catalogue_reader")
_loci = import_module("06_define_risk_loci")


def fetch_regions_remote(url: str, regions: Dict[str, List[Tuple[int, int]]],
                         out_path: pathlib.Path, want) -> int:
    """Pull the loci windows from a remote tabix endpoint.

    Parameters
    ----------
    url : remote .tsv.gz accepted by tabix.
    regions : chromosome -> list of (start, end) intervals.
    out_path : destination TSV, written with a header.
    want : column indices resolved from the local header layout.

    Returns
    -------
    Number of rows retrieved.

    Notes
    -----
    Used when the dataset was not mirrored locally. It issues one tabix call
    per interval, so the cost scales with the number of loci rather than with
    the size of the full table.
    """
    import subprocess

    kept = 0
    with open(out_path, "w") as out:
        out.write("\t".join(_eqtl.OUT_COLS) + "\n")
        for chrom, spans in regions.items():
            for lo, hi in spans:
                proc = subprocess.run(["tabix", url, f"{chrom}:{lo}-{hi}"],
                                      capture_output=True, text=True)
                for line in proc.stdout.strip().split("\n"):
                    if not line:
                        continue
                    fields = line.split("\t")
                    out.write("\t".join(fields[j] for j in want) + "\n")
                    kept += 1
    return kept


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tissue", required=True, help="GTEx tissue label")
    ap.add_argument("--gwas", default=GWAS_PRIMARY_FILE, help="outcome file in data/")
    ap.add_argument("--eqtl-id", default=None, help="override the eQTL dataset id")
    ap.add_argument("--remote", action="store_true",
                    help="use the EBI tabix endpoint instead of a local file")
    ap.add_argument("--force-rescan", action="store_true",
                    help="rebuild the locus subset even if a cached copy exists")
    args = ap.parse_args()

    PATHS.ensure()
    t0 = time.time()

    gwas_path = PATHS.data / args.gwas
    if not gwas_path.exists():
        raise SystemExit(f"missing outcome file: {gwas_path}")

    qtl_id = args.eqtl_id or TISSUES.get(args.tissue)
    if qtl_id is None:
        raise SystemExit(f"unknown tissue {args.tissue}; add it to config.TISSUES")

    # ---- 1. risk loci, identical for every tissue -----------------------
    store = _gwas.load_gwas(gwas_path)
    n_variants = sum(len(v) for v in store.values())
    print(f"[1/4] outcome loaded: {n_variants:,} variants "
          f"({time.time() - t0:.0f}s)", flush=True)

    loci = _gwas.peaks(store, PARAMS.gwas_p_thresh, PARAMS.locus_gap_bp)
    # Coordinates are already in the build of the outcome file; no lift-over.
    print(f"      {len(loci)} loci at p<{PARAMS.gwas_p_thresh:g}", flush=True)

    # ---- 2. eGenes inside each locus -----------------------------------
    tar_path = PATHS.data / "GTEx_Analysis_v8_eQTL.tar"
    if not tar_path.exists():
        raise SystemExit(f"missing GTEx archive: {tar_path}")
    genes_by_locus = _loci.genes_near_loci(tar_path, args.tissue, loci,
                                           PARAMS.region_pad_bp)
    print("[2/4] candidate eGenes per locus: "
          + str({i: len(v) for i, v in genes_by_locus.items()}), flush=True)

    # ---- 3. restrict the eQTL table to the loci windows ----------------
    regions: Dict[str, List[Tuple[int, int]]] = {}
    for row in loci.itertuples(index=False):
        regions.setdefault(str(row.chrom), []).append(
            (max(1, int(row.lead) - PARAMS.region_pad_bp),
             int(row.lead) + PARAMS.region_pad_bp))
    for chrom in regions:
        regions[chrom].sort()

    local_eqtl = PATHS.data / f"{qtl_id}.all.tsv.gz"
    subset = PATHS.cache / f"loci_{qtl_id}_{args.tissue}.tsv"
    subset.parent.mkdir(parents=True, exist_ok=True)

    if args.remote:
        url = f"{EQTL_CATALOGUE_BASE}/{qtl_id}/{qtl_id}.all.tsv.gz"
        if args.force_rescan or not subset.exists():
            n = fetch_regions_remote(url, regions, subset, _eqtl.column_positions(local_eqtl)
                                     if local_eqtl.exists() else list(range(len(_eqtl.OUT_COLS))))
            print(f"[3/4] remote fetch: {n:,} rows", flush=True)
    else:
        if not local_eqtl.exists():
            raise SystemExit(f"missing eQTL file: {local_eqtl}\n"
                             f"  run again with --remote, or place the file in data/")
        if args.force_rescan or not subset.exists():
            n = _eqtl.stream_regions(local_eqtl, regions, subset)
            print(f"[3/4] streamed: kept {n:,} rows ({time.time() - t0:.0f}s)", flush=True)
        else:
            print(f"[3/4] reusing cached subset {subset.name}", flush=True)

    qtl = pd.read_csv(subset, sep="\t")
    qtl = qtl.dropna(subset=["beta", "se"])
    qtl = qtl[(qtl["maf"] >= PARAMS.eqtl_maf_min)
              & (qtl["maf"] <= PARAMS.eqtl_maf_max) & (qtl["se"] > 0)]
    print(f"      {len(qtl):,} usable eQTL records", flush=True)

    # ---- 4. colocalise gene by gene ------------------------------------
    records = []
    for i, row in loci.iterrows():
        genes = genes_by_locus.get(i, [])
        if not genes:
            continue
        lead = int(row.lead)
        near = qtl[(qtl["chromosome"].astype(str) == str(row.chrom))
                   & (qtl["position"].between(lead - PARAMS.region_pad_bp,
                                              lead + PARAMS.region_pad_bp))]
        if near.empty:
            continue
        for gene in genes:
            sub = near[near["gene_id"] == gene]
            if len(sub) < PARAMS.eqtl_min_snps:
                continue
            pairs = []
            for x in sub.itertuples(index=False):
                hit = _gwas.lookup(store, str(row.chrom), int(x.position),
                                   str(x.ref).upper(), str(x.alt).upper())
                if hit is None:
                    continue
                beta_g, se_g, p_g = hit
                if np.isfinite(beta_g) and np.isfinite(se_g) and se_g > 0:
                    pairs.append((x.beta, x.se, beta_g, se_g, p_g))
            if len(pairs) < PARAMS.eqtl_min_snps:
                continue
            d = pd.DataFrame(pairs, columns=["b_e", "se_e", "b_g", "se_g", "p_g"])
            res = _abf.coloc_abf(
                _abf.log_abf(d["b_e"].values, d["se_e"].values, 1.0),
                _abf.log_abf(d["b_g"].values, d["se_g"].values, 1.0),
                p1=PARAMS.p1, p2=PARAMS.p2, p12=PARAMS.p12,
            )
            records.append({
                **{k: v for k, v in res.items() if k not in ("n_snps", "top_snp_share")},
                "n_snps": len(d),
                "top_snp_share": res["top_snp_share"],
                "tissue": args.tissue,
                "locus": f"{row.chrom}:{lead}",
                "gene_id": gene,
                "best_gwas_p": float(d["p_g"].min()),
            })

    if not records:
        raise SystemExit("no testable gene-locus pairs; check the tissue label and the window")

    out = pd.DataFrame(records)[
        ["PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "n_snps", "top_snp_share",
         "tissue", "locus", "gene_id", "best_gwas_p"]
    ].sort_values("PP.H4", ascending=False)

    dest = PATHS.tables / f"coloc_{args.tissue}.tsv"
    out.to_csv(dest, sep="\t", index=False)

    print(f"[4/4] {len(out)} tests -> {dest.name}", flush=True)
    print(f"      PP.H4>0.8: {(out['PP.H4'] > PARAMS.pp_h4_high).sum()}"
          f"  PP.H4>0.5: {(out['PP.H4'] > 0.5).sum()}"
          f"  max PP.H4: {out['PP.H4'].max():.3f}"
          f"  PP.H3>0.8: {(out['PP.H3'] > PARAMS.pp_h3_high).sum()}", flush=True)
    print(f"      total {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
