"""
09_dual_outcome.py
==================

Purpose
-------
Rerun the colocalisation layer against the pre-registered second outcome,
GCST006465 (O'Mara et al. 2018, GRCh37), and write the result beside the
primary-outcome table so the two can be compared row by row.

This is a consistency check, not a replication. GCST006465 and the primary
outcome GCST90454186 share the ECAC core samples, so agreement only means "the
positive signal survives a change of outcome dataset". The manuscript must not
describe the result of this script as replication.

Why the joins are built this way
--------------------------------
The two sides live on different builds, so both directions of the UCSC chain are
used inside a single run:

    hg19 -> hg38   locus lead, to place the eQTL window on the expression build
    hg38 -> hg19   every eQTL variant in that window, to match the outcome

Matching is on (position, ref, alt) rather than on rsIDs. An rsID merge silently
drops variants whose identifier is missing or differs between releases, and it
hides orientation errors; the allele-aware reader in 02_gwas_reader flips the
sign when the outcome stores the alleles in the opposite order.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
--gwas        outcome file name inside the data directory
              (default: config.GWAS_SECONDARY_FILE, GRCh37)
--tissue      one GTEx tissue label; omitted means every entry of config.TISSUES
--min-snps    minimum paired variants for a gene to be tested
              (default: config.PARAMS.eqtl_min_snps)
--remote      query the EBI tabix endpoint; an all-associations file that was
              never mirrored locally can only be read this way
--force-rescan rebuild the cached locus subset instead of reusing it

Loci are discovered with the thresholds in config.PARAMS (gwas_p_thresh,
locus_gap_bp) and the window is config.PARAMS.region_pad_bp; none are repeated
on the command line.

Outputs
-------
results/tables/coloc_omara_confirm.tsv, one row per gene-locus pair:

    PP.H0 .. PP.H4   posterior probabilities of the five configurations
    n_snps           variants entering the test
    tissue, locus    GTEx tissue and "lifted chrom:hg19 lead position"
    gene_id          Ensembl id of the eQTL gene
    best_gwas_p      strongest association inside the window
    alleles_matched  variants whose alleles agreed after the reverse lift
    alleles_dropped  variants in the window that could not be paired

Notes
-----
``alleles_dropped`` counts variants that are either untested in the outcome or
stored with incompatible alleles; the two cases are indistinguishable from the
outcome file alone. A large value therefore points at an allele-coding mismatch
between the two releases, not at a shortage of variants.

The ``locus`` label keeps the hg19 coordinate because that is the coordinate
the earlier release of this table used; changing it would break the join with
results/tables/table1_coloc_by_tissue.tsv.

A local eQTL file is read in one streaming pass and the locus subset is cached
under PATHS.cache, rather than queried per locus with tabix: several GTEx v8
files in this release were distributed without an index, and on those tabix
returns nothing while exiting successfully.

``top_snp_share`` is not written here. This path was run without the shared-
variant posterior diagnostic, and writing a placeholder column would invite the
reader to compare it with the diagnostic reported by 07.
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

from config.config import (PATHS, PARAMS, TISSUES, GWAS_SECONDARY_FILE,
                           EQTL_CATALOGUE_BASE, CHAIN_HG19_TO_HG38,
                           CHAIN_HG38_TO_HG19)
from importlib import import_module

_abf = import_module("05_coloc_abf_engine")
_gwas = import_module("02_gwas_reader")
_eqtl = import_module("03_eqtl_catalogue_reader")
_gtex = import_module("13_gtex_v11_reader")
_lift = import_module("01_chain_liftover")

# The significant-pairs archive is a single tar; only the wanted member is read.
GTEX_TAR = "GTEx_Analysis_v8_eQTL.tar"

OUT_COLS = ["PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "n_snps", "tissue",
            "locus", "gene_id", "best_gwas_p", "alleles_matched",
            "alleles_dropped"]


def load_chain(name: str):
    """Open one UCSC chain file from the data directory.

    Parameters
    ----------
    name : file name of the chain, taken from config.

    Returns
    -------
    A ``Chain`` instance (01_chain_liftover).

    Raises
    ------
    SystemExit
        If the chain is absent. A missing chain must stop the run: silently
        skipping the lift would either drop every locus or, worse, match the two
        builds by position coincidence.
    """
    path = PATHS.data / name
    if not path.exists():
        raise SystemExit(f"missing chain file: {path}\n"
                         f"  stage 02/00_download_inputs.sh fetches it")
    return _lift.Chain(str(path))


def lift_loci(store: Dict[str, pd.DataFrame], up) -> pd.DataFrame:
    """Cluster the second outcome and carry each lead forward to GRCh38.

    Parameters
    ----------
    store : per-chromosome outcome table from 02_gwas_reader.load_gwas, in the
        build of the file itself (GRCh37 here).
    up : hg19 -> hg38 chain.

    Returns
    -------
    DataFrame with chrom, start, end, lead (all GRCh38) and lead19, p (GRCh37).
    Loci whose lead cannot be lifted are dropped.

    Notes
    -----
    ``start`` and ``end`` are both the lifted lead: the window is later built by
    adding config.PARAMS.region_pad_bp on either side, which is what the module
    06 helper expects.
    """
    # A store with no variant below the threshold makes 02_gwas_reader.peaks
    # fail on an empty frame, which surfaces as a KeyError far from the cause.
    # Checking here turns it into a statement about the outcome file.
    n_sig = sum(int((d["pval"] < PARAMS.gwas_p_thresh).sum())
                for d in store.values())
    if n_sig == 0:
        raise SystemExit(
            f"no variant reaches p<{PARAMS.gwas_p_thresh:g} in this outcome; "
            f"there is no locus to carry forward")

    loc19 = _gwas.peaks(store, PARAMS.gwas_p_thresh, PARAMS.locus_gap_bp)
    rows: List[Dict[str, object]] = []
    for row in loc19.itertuples(index=False):
        mapped = up.lift(str(row.chrom), int(row.lead))
        if mapped is None:
            continue
        chrom38, lead38 = str(mapped[0]), int(mapped[1])
        rows.append({"chrom": chrom38, "start": lead38, "end": lead38,
                     "lead": lead38, "p": float(row.p), "lead19": int(row.lead)})
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "lead", "p",
                                       "lead19"])


def locus_regions(loc38: pd.DataFrame, pad: int) -> Dict[str, List[Tuple[int, int]]]:
    """Collect the eQTL windows of every lifted locus.

    Parameters
    ----------
    loc38 : output of :func:`lift_loci`.
    pad : half-width in base pairs.

    Returns
    -------
    dict chromosome -> sorted list of (start, end) intervals, ready for
    03_eqtl_catalogue_reader.stream_regions.
    """
    regions: Dict[str, List[Tuple[int, int]]] = {}
    for row in loc38.itertuples(index=False):
        lead = int(row.lead)
        regions.setdefault(str(row.chrom), []).append(
            (max(1, lead - pad), lead + pad))
    for chrom in regions:
        regions[chrom].sort()
    return regions


def want_columns(src: str) -> List[int]:
    """Column indices of the fields retained by the eQTL reader.

    Parameters
    ----------
    src : local path or remote URL.

    Returns
    -------
    List of 0-based indices aligned with 03_eqtl_catalogue_reader.OUT_COLS.

    Notes
    -----
    The layout is read from a header, never assumed. A remote URL cannot be
    opened for its header, so the header of any local eQTL file is used as the
    donor: every GTEx v8 release of the eQTL Catalogue shares one column order.
    Lacking both, the run stops rather than reading the wrong columns.
    """
    if not str(src).startswith(("http://", "https://")):
        return _eqtl.column_positions(src)

    donors = sorted(PATHS.data.glob("*.all.tsv.gz"))
    if not donors:
        raise SystemExit("remote eQTL source: at least one local "
                         "*.all.tsv.gz is needed to read the column layout")
    return _eqtl.column_positions(donors[0])


def load_tissue_variants(loc38: pd.DataFrame, qtd: str, tissue: str,
                         allow_remote: bool, force_rescan: bool) -> pd.DataFrame:
    """Read every eQTL row that falls in a locus window, for one tissue.

    Parameters
    ----------
    loc38 : lifted loci, which define the windows.
    qtd, tissue : dataset id and its label.
    allow_remote : whether the EBI endpoint may be used instead of a local file.
    force_rescan : rebuild the cached subset even if one exists.

    Returns
    -------
    DataFrame with the columns of 03_eqtl_catalogue_reader.OUT_COLS, all loci
    concatenated.

    Notes
    -----
    A local file is streamed in one pass and the subset is cached, rather than
    queried per locus with tabix. Most of the GTEx v8 files in this release were
    distributed without a tabix index; on those, tabix exits successfully with no
    output, which is indistinguishable from a window with no variants. Streaming
    cannot fail that way, and the cache makes a rerun cheap.

    The remote path has to go locus by locus, because there is nothing local to
    stream over.
    """
    regions = locus_regions(loc38, PARAMS.region_pad_bp)
    local = PATHS.data / f"{qtd}.all.tsv.gz"

    if allow_remote and not local.exists():
        url = f"{EQTL_CATALOGUE_BASE}/{qtd}/{qtd}.all.tsv.gz"
        want = want_columns(url)
        frames = []
        for chrom, spans in regions.items():
            for lo, hi in spans:
                part = _eqtl.qtl_region(chrom, lo, hi, url, want)
                if not part.empty:
                    frames.append(part)
        if not frames:
            raise SystemExit(f"remote query returned nothing for {qtd}; "
                             f"the endpoint may be unavailable or the index "
                             f"missing")
        return pd.concat(frames, ignore_index=True)

    if not local.exists():
        raise SystemExit(f"missing eQTL file: {local.name}\n"
                         f"  run again with --remote, or let 15_locus_screen.py "
                         f"download it first")

    subset = PATHS.cache / f"dual_{qtd}_{tissue}.tsv"
    if force_rescan or not subset.exists():
        n = _eqtl.stream_regions(local, regions, subset)
        print(f"      streamed {n:,} rows -> {subset.name}", flush=True)
    return pd.read_csv(subset, sep="\t")


def tests_for_tissue(tissue: str, loc38: pd.DataFrame,
                     store: Dict[str, pd.DataFrame], down,
                     qtl_all: pd.DataFrame) -> List[Dict[str, object]]:
    """Run every gene-locus test for one tissue.

    Parameters
    ----------
    tissue : GTEx label.
    loc38 : output of :func:`lift_loci`.
    store : the second outcome, GRCh37.
    down : hg38 -> hg19 chain, applied to each eQTL variant.
    qtl_all : every eQTL row inside any locus window, from
        :func:`load_tissue_variants`.

    Returns
    -------
    List of records, one per testable gene-locus pair.

    Notes
    -----
    The reverse lift runs once per variant over the whole subset, and genes in
    one locus share the resulting variant set. That is what makes the tests
    within a locus comparable, and it is the same convention 07 applies on the
    forward build.
    """
    tar_path = PATHS.data / GTEX_TAR
    if not tar_path.exists():
        raise SystemExit(f"missing GTEx archive: {tar_path.name}")

    genes_by_locus = _gtex.genes_near_loci(tar_path, tissue, loc38,
                                           PARAMS.region_pad_bp)

    qtl_all = qtl_all.dropna(subset=["position", "beta", "se"])
    qtl_all = qtl_all[(qtl_all["maf"] >= PARAMS.eqtl_maf_min)
                      & (qtl_all["maf"] <= PARAMS.eqtl_maf_max)
                      & (qtl_all["se"] > 0)]
    if qtl_all.empty:
        return []

    # One reverse-lift pass over the cached subset, keyed by
    # (chromosome, position) so that a position is lifted at most once.
    lifted: Dict[Tuple[str, int], int] = {}
    pos19 = np.full(len(qtl_all), np.nan)
    for j, (chrom, pos) in enumerate(zip(qtl_all["chromosome"].astype(str),
                                         qtl_all["position"].to_numpy())):
        key = (chrom, int(pos))
        if key not in lifted:
            mapped = down.lift(chrom, key[1])
            if mapped is not None:
                lifted[key] = int(mapped[1])
        if key in lifted:
            pos19[j] = lifted[key]
    qtl_all = qtl_all.assign(pos19=pos19).dropna(subset=["pos19"])
    if qtl_all.empty:
        return []
    qtl_all["pos19"] = qtl_all["pos19"].astype(int)

    records: List[Dict[str, object]] = []
    for i, row in loc38.iterrows():
        genes = genes_by_locus.get(i, [])
        if not genes:
            continue

        chrom, lead = str(row["chrom"]), int(row["lead"])
        window = qtl_all[
            (qtl_all["chromosome"].astype(str) == chrom)
            & (qtl_all["pos19"].between(lead - PARAMS.region_pad_bp,
                                        lead + PARAMS.region_pad_bp))]
        if window.empty:
            continue

        for gene in genes:
            sub = window[window["gene_id"] == gene]
            if len(sub) < PARAMS.eqtl_min_snps:
                continue

            pairs = []
            for x in sub.itertuples(index=False):
                hit = _gwas.lookup(store, chrom, int(x.pos19),
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
            res = {k: v for k, v in res.items()
                   if k not in ("n_snps", "top_snp_share")}
            records.append({
                **res,
                "n_snps": len(d),
                "tissue": tissue,
                "locus": f"{chrom}:{int(row['lead19'])}",
                "gene_id": gene,
                "best_gwas_p": float(d["p_g"].min()),
                "alleles_matched": len(d),
                "alleles_dropped": int(len(sub) - len(d)),
            })
    return records


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gwas", default=GWAS_SECONDARY_FILE,
                    help="outcome file in data/ (default: the second outcome)")
    ap.add_argument("--tissue", default=None,
                    help="single GTEx tissue label (default: every configured tissue)")
    ap.add_argument("--min-snps", type=int, default=PARAMS.eqtl_min_snps,
                    help="minimum paired variants per test")
    ap.add_argument("--remote", action="store_true",
                    help="use the EBI tabix endpoint instead of a local file")
    ap.add_argument("--force-rescan", action="store_true",
                    help="rebuild the cached locus subset even if one exists")
    args = ap.parse_args()

    PATHS.ensure()
    t0 = time.time()

    if args.tissue is not None and args.tissue not in TISSUES:
        raise SystemExit(f"unknown tissue {args.tissue}; add it to config.TISSUES")
    tissues = ({args.tissue: TISSUES[args.tissue]} if args.tissue
               else dict(TISSUES))

    gwas_path = pathlib.Path(args.gwas)
    if not gwas_path.is_absolute():
        gwas_path = PATHS.data / gwas_path
    if not gwas_path.exists():
        raise SystemExit(f"missing outcome file: {gwas_path}")

    # ---- 1. outcome and loci, in the build of the file ------------------
    store = _gwas.load_gwas(gwas_path)
    n_var = sum(len(v) for v in store.values())
    print(f"[1/4] outcome loaded: {n_var:,} variants "
          f"({time.time() - t0:.0f}s)", flush=True)

    up = load_chain(CHAIN_HG19_TO_HG38)
    loc38 = lift_loci(store, up)
    if loc38.empty:
        raise SystemExit("no locus survived the lift to GRCh38")
    print(f"      {len(loc38)} loci at p<{PARAMS.gwas_p_thresh:g}, "
          f"{len(loc38)} lifted to GRCh38", flush=True)

    # ---- 2. per-tissue tests -------------------------------------------
    down = load_chain(CHAIN_HG38_TO_HG19)
    records: List[Dict[str, object]] = []
    for tissue, qtd in tissues.items():
        t_tissue = time.time()
        qtl_all = load_tissue_variants(loc38, qtd, tissue, args.remote,
                                       args.force_rescan)
        recs = tests_for_tissue(tissue, loc38, store, down, qtl_all)
        records.extend(recs)
        print(f"[2/4] {tissue}: {len(recs)} tests "
              f"({time.time() - t_tissue:.0f}s)", flush=True)

    if not records:
        raise SystemExit("no testable gene-locus pairs\n"
                         "  the outcome window may be empty, the tissue label "
                         "may be wrong, or the eQTL file may lack variants")

    # ---- 3. write ------------------------------------------------------
    out = pd.DataFrame(records)[OUT_COLS]
    dest = PATHS.tables / "coloc_omara_confirm.tsv"
    out.to_csv(dest, sep="\t", index=False)
    print(f"[3/4] {len(out)} tests -> {dest.name}", flush=True)

    # ---- 4. summary ----------------------------------------------------
    print(f"[4/4] PP.H4>{PARAMS.pp_h4_high:g}: "
          f"{(out['PP.H4'] > PARAMS.pp_h4_high).sum()}"
          f"  PP.H4>0.5: {(out['PP.H4'] > 0.5).sum()}"
          f"  max PP.H4: {out['PP.H4'].max():.3f}"
          f"  PP.H3>{PARAMS.pp_h3_high:g}: "
          f"{(out['PP.H3'] > PARAMS.pp_h3_high).sum()}", flush=True)

    hits = out[out["PP.H4"] > PARAMS.pp_h4_high]
    if len(hits):
        print("      shared-signal hits:", flush=True)
        # Column-wise iteration: the PP.* names are not valid attribute names,
        # so itertuples would mangle them.
        for tissue, locus, gene, h4, n in zip(hits["tissue"], hits["locus"],
                                              hits["gene_id"], hits["PP.H4"],
                                              hits["n_snps"]):
            print(f"        {tissue:24s} {locus:16s} {gene:20s} "
                  f"H4={h4:.3f} n={n}", flush=True)
    print(f"      total {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
