"""
15_locus_screen.py
==================

Purpose
-------
The toolkit that sits around the colocalisation test: bring one tissue's eQTL
data in, reduce a whole transcriptome to a candidate list, run the region-level
screen, check the leads against the second outcome, and compare two result
tables. Each job is a subcommand because they are run at different times and for
different reasons; keeping them in one file means the variant matching and the
variant-set construction exist once.

    download   fetch one eQTL Catalogue dataset, verify its size against the
               remote object, export the eGene list from the GTEx archive and
               build the gene window index
    screen     first pass for one tissue: does the lead eQTL variant of each
               gene associate with the outcome, and in which direction. This is
               the no-LD approximation of SMR, not colocalisation; its job is to
               cut thousands of genes down to a candidate list
    coloc      region-level colocalisation on each gene's significant variant
               set, with the permutation null of 12 available as an option
    replicate  take the lead variants into the second outcome, which is on
               GRCh37, and report direction concordance
    compare    print the summary statistics of two result tables side by side

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
--tissue    GTEx tissue label or one of the short aliases (uterus, liver,
            muscle, pancreas, adipose_visceral); the dataset id always comes
            from config.TISSUES
--gwas      outcome file name inside the data directory; the layout is detected
            from its header by 02_gwas_reader
--qtd       dataset id, only needed by `download` for a tissue that is not in
            config.TISSUES
--min-snps  minimum variants for a gene to be tested (default from config.PARAMS)
--min-af    minimum minor allele frequency for a lead variant (screen only)
--null      number of permutation replicates for `coloc` (0 = skip)
--limit     first N genes only, for a timing check

Outputs
-------
results/tables/screen_<tissue>.tsv            one row per gene lead variant
results/tables/coloc_<tissue>.tsv             one row per gene tested
results/tables/screen_<tissue>_gcst006465.tsv one row per gene lead variant

``download`` writes into the data directory and into PATHS.cache: it is the only
subcommand that fetches anything, and it never overwrites a complete file.

Notes
-----
The region-level screen is biased, deliberately and in a known direction. It
uses each gene's *significant* variants, whereas standard colocalisation uses
every variant that was tested; a truncated variant set inflates PP.H4 because
every retained variant is already expression-associated. The absolute values
must therefore not be quoted as colocalisation evidence. What is valid is the
ordering, and only after the permutation calibration in 12 has removed the
gene-specific inflation.

``screen`` no longer writes an allele-frequency column. The shared outcome reader
keeps position, alleles, effect, standard error and p-value only; the screening
decision uses direction and p-value, and carrying a frequency column filled with
placeholders would be worse than omitting it. The released table has the column,
so a re-run produces a narrower file; the join key is ``gene_id``.

Allele matching here is direct or swapped, which is all that one build requires.
Complement matching was part of the earlier implementation because it also
served a cross-build join; that path now lives in `replicate`, where the lift
actually happens.
"""
from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys
import time
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats

# Local package imports; the directory holding this file is on sys.path.
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))       # repository root, for config/

from config.config import (PATHS, PARAMS, TISSUES, GWAS_PRIMARY_FILE,
                           GWAS_SECONDARY_FILE, EQTL_CATALOGUE_BASE,
                           CHAIN_HG38_TO_HG19, SEED)
from importlib import import_module

_abf = import_module("05_coloc_abf_engine")
_gwas = import_module("02_gwas_reader")
_eqtl = import_module("03_eqtl_catalogue_reader")
_gtex = import_module("13_gtex_v11_reader")
_lift = import_module("01_chain_liftover")

# Short aliases accepted on the command line, mapped onto config labels.
ALIASES = {
    "uterus": "Uterus",
    "liver": "Liver",
    "muscle": "Muscle_Skeletal",
    "pancreas": "Pancreas",
    "adipose_visceral": "Adipose_Visceral_Omentum",
    "adipose": "Adipose_Visceral_Omentum",
}

GTEX_TAR = "GTEx_Analysis_v8_eQTL.tar"

# Columns of the region-level table, in the order the released file used.
COLOC_COLS = ["PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "n_snps",
              "top_snp_share", "best_gwas_p", "top_sign_concordant", "gene_id"]

# Columns of the lead-variant screen.
SCREEN_COLS = ["gene_id", "variant_id", "eQTL_slope", "eQTL_p", "in_gwas",
               "gwas_beta", "gwas_se", "gwas_p", "allele_mismatch"]

# Complement map, used only by the cross-build lead check.
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def resolve_tissue(label: str) -> Tuple[str, str]:
    """Map a label or alias onto (config label, dataset id).

    Parameters
    ----------
    label : config label or short alias.

    Returns
    -------
    (tissue, qtd) where tissue is the key used in config.TISSUES.

    Raises
    ------
    SystemExit
        When the tissue is unknown, so that a typo cannot silently produce an
        empty result table.
    """
    tissue = ALIASES.get(label.lower() if isinstance(label, str) else label, label)
    qtd = TISSUES.get(tissue)
    if qtd is None:
        raise SystemExit(f"unknown tissue {label}; known labels: "
                         f"{sorted(TISSUES)} (aliases: {sorted(ALIASES)})")
    return tissue, qtd


def resolve_input(name: str) -> pathlib.Path:
    """Resolve a user-supplied file name against the data directory."""
    path = pathlib.Path(name)
    return path if path.is_absolute() else PATHS.data / path


def load_outcome(name: str) -> Dict[str, pd.DataFrame]:
    """Load an outcome file, stopping when it is absent."""
    path = resolve_input(name)
    if not path.exists():
        raise SystemExit(f"missing outcome file: {path}")
    store = _gwas.load_gwas(path)
    print(f"  outcome: {sum(len(v) for v in store.values()):,} variants",
          flush=True)
    return store


def gtex_tar() -> pathlib.Path:
    """Path to the GTEx significant-pairs archive, checked for existence."""
    path = PATHS.data / GTEX_TAR
    if not path.exists():
        raise SystemExit(f"missing GTEx archive: {path.name}")
    return path


def gene_variants(tissue: str):
    """Group one tissue's significant pairs by gene.

    Parameters
    ----------
    tissue : GTEx tissue label.

    Returns
    -------
    A groupby over gene_id, with pos populated and rows lacking a usable
    position removed.
    """
    pairs = _gtex.split_variant(_gtex.load_tissue(gtex_tar(), tissue))
    pairs = pairs[pairs["pos"].notna()]
    return pairs.groupby("gene_id")


def build_snp_table(gene_df: pd.DataFrame,
                    store: Dict[str, pd.DataFrame]) -> Optional[pd.DataFrame]:
    """Align one gene's significant eQTL variants with the outcome.

    Parameters
    ----------
    gene_df : rows of one gene, after 13_gtex_v11_reader.split_variant.
    store : per-chromosome outcome table.

    Returns
    -------
    DataFrame with variant_id, b_e, se_e, b_g, se_g, p_g, or None when fewer
    than two variants can be paired.

    Notes
    -----
    The earlier version fetched a tabix interval expanded by 10 kb around the
    gene and built a position-allele lookup from it. The expansion was an
    artefact of interval queries: a position-keyed lookup needs no window, and
    the variant set is defined by the eQTL rows, so nothing outside them is
    reachable either way.

    Two variants is the floor here rather than a scientific threshold; the ABF
    test needs at least that to distinguish H3 from H4, and the calling script
    applies its own minimum.
    """
    rows = []
    for x in gene_df.itertuples(index=False):
        chrom = str(x.chrom)
        hit = _gwas.lookup(store, chrom, int(x.pos), str(x.ref).upper(),
                           str(x.alt).upper())
        if hit is None:
            continue
        beta_g, se_g, p_g = hit
        if not (np.isfinite(beta_g) and np.isfinite(se_g)) or se_g <= 0:
            continue
        rows.append((x.variant_id, x.slope, x.slope_se, beta_g, se_g, p_g))
    if len(rows) < 2:
        return None
    return pd.DataFrame(rows, columns=["variant_id", "b_e", "se_e", "b_g",
                                       "se_g", "p_g"])


def coloc_one(df: pd.DataFrame) -> Dict[str, object]:
    """Run the ABF test on one aligned variant table.

    Parameters
    ----------
    df : output of :func:`build_snp_table`.

    Returns
    -------
    dict of PP.H0..H4, n_snps, top_snp_share, best_gwas_p and
    top_sign_concordant.

    Notes
    -----
    ``top_sign_concordant`` records whether the strongest outcome variant has
    the same effect direction on expression. A high PP.H4 with a discordant sign
    at the top variant is the pattern produced by two overlapping but opposite
    signals, and is worth seeing before the number is believed.
    """
    l1 = _abf.log_abf(df["b_e"].values, df["se_e"].values, 1.0)
    l2 = _abf.log_abf(df["b_g"].values, df["se_g"].values, 1.0)
    res = _abf.coloc_abf(l1, l2, p1=PARAMS.p1, p2=PARAMS.p2, p12=PARAMS.p12)
    top = df.loc[df["p_g"].idxmin()]
    return {
        "PP.H0": res["PP.H0"], "PP.H1": res["PP.H1"], "PP.H2": res["PP.H2"],
        "PP.H3": res["PP.H3"], "PP.H4": res["PP.H4"],
        "n_snps": len(df),
        "top_snp_share": res["top_snp_share"],
        "best_gwas_p": float(df["p_g"].min()),
        "top_sign_concordant": bool(top["b_e"] * top["b_g"] > 0),
    }


# --------------------------------------------------------------------------
# download
# --------------------------------------------------------------------------
def remote_size(url: str) -> int:
    """Content-Length of a remote file, or 0 when it cannot be read."""
    proc = subprocess.run(["curl", "-sI", "-m", "60", url],
                          capture_output=True, text=True)
    for line in proc.stdout.splitlines():
        if line.lower().startswith("content-length"):
            try:
                return int(line.split(":", 1)[1].strip())
            except ValueError:
                return 0
    return 0


def cmd_download(args: argparse.Namespace) -> None:
    """Fetch one dataset and prepare the files the other stages read.

    Parameters
    ----------
    args : parsed command line (tissue, qtd, coloc).

    Notes
    -----
    Completeness is judged against the remote Content-Length, never against the
    presence of a file: a partially downloaded bgzip still opens, and tabix then
    returns zero rows with a zero exit status, which is indistinguishable from
    "no signal in this window" and has already caused one empty analysis.

    The earlier version finished by calling a colocalisation script that lived
    in a different directory. That script is superseded by 07_run_colocalisation,
    so the hand-off is now explicit: ``--coloc`` runs it, otherwise the command
    to run is printed.
    """
    tissue, qtd = (resolve_tissue(args.tissue) if args.tissue
                   else (None, args.qtd))
    if qtd is None:
        raise SystemExit("give --tissue or --qtd")
    label = tissue or qtd

    gz = PATHS.data / f"{qtd}.all.tsv.gz"
    url = f"{EQTL_CATALOGUE_BASE}/{qtd}/{qtd}.all.tsv.gz"

    expected = remote_size(url)
    complete = gz.exists() and expected and gz.stat().st_size == expected
    if complete:
        print(f"[1/4] {gz.name} already complete "
              f"({gz.stat().st_size / 1e9:.2f} GB)", flush=True)
    else:
        if gz.exists():
            print(f"[1/4] resuming {gz.name}: "
                  f"{gz.stat().st_size:,} of {expected:,} bytes", flush=True)
        else:
            print(f"[1/4] downloading {qtd}", flush=True)
        subprocess.run(["curl", "-sS", "-C", "-", "--retry", "5",
                        "--retry-delay", "5", "-m", "5400", "-o", str(gz), url],
                       check=True)
        # The tabix index is optional for the streaming reader, so a failure here
        # must not abort a download that is otherwise complete.
        subprocess.run(["curl", "-sS", "-o", f"{gz}.tbi", f"{url}.tbi"],
                       check=False)

    if not gz.exists():
        raise SystemExit(f"download failed: {gz}")
    print(f"[1/4] {gz.name}: {gz.stat().st_size / 1e9:.2f} GB", flush=True)

    # ---- eGene list, from the GTEx archive rather than from the file ----
    tar_path = gtex_tar()
    egenes = PATHS.cache / f"egenes_{label}.txt"
    genes = _gtex.tissue_genes(_gtex.load_tissue(tar_path, label))
    egenes.write_text("\n".join(genes) + "\n")
    print(f"[2/4] {len(genes):,} eGenes from {label} -> {egenes.name}", flush=True)

    # ---- gene window index ---------------------------------------------
    index_path = PATHS.cache / f"gene_index_{qtd}.tsv"
    if not index_path.exists() or index_path.stat().st_size == 0:
        print("[3/4] building the gene window index", flush=True)
        _eqtl.build_gene_index(gz).to_csv(index_path, sep="\t", index=False)
    else:
        print(f"[3/4] reusing {index_path.name}", flush=True)

    # ---- hand-off -------------------------------------------------------
    if args.coloc:
        script = HERE / "07_run_colocalisation.py"
        print(f"[4/4] running {script.name} --tissue {label}", flush=True)
        subprocess.run([sys.executable, str(script), "--tissue", label],
                       check=True)
    else:
        print(f"[4/4] next: {sys.executable} "
              f"{HERE / '07_run_colocalisation.py'} --tissue {label}", flush=True)


# --------------------------------------------------------------------------
# screen
# --------------------------------------------------------------------------
def cmd_screen(args: argparse.Namespace) -> None:
    """Association of each gene's lead eQTL variant with the outcome.

    Parameters
    ----------
    args : parsed command line.

    Notes
    -----
    This is a screen, not a test: no correction is applied and no threshold is
    claimed. Its output exists to be narrowed by ``coloc``, and the binomial
    test on the direction of effect is the only inference it supports.
    """
    tissue, _qtd = resolve_tissue(args.tissue)
    t0 = time.time()

    pairs = _gtex.split_variant(_gtex.load_tissue(gtex_tar(), tissue))
    print(f"[1/4] {len(pairs):,} significant pairs, "
          f"{pairs['gene_id'].nunique():,} genes", flush=True)

    leads = _gtex.lead_per_gene(pairs)
    if "maf" not in leads.columns:
        raise SystemExit("the GTEx archive member has no maf column; the lead "
                         "frequency filter cannot be applied, and testing "
                         "every lead regardless of frequency would change the "
                         "screen")
    leads = leads[leads["maf"].between(args.min_af, 1 - args.min_af)
                  & leads["pos"].notna()]
    if args.limit:
        leads = leads.head(args.limit)
    print(f"[2/4] one lead per gene: {len(leads):,}", flush=True)

    store = load_outcome(args.gwas)
    rows: List[Dict[str, object]] = []
    for k, r in enumerate(leads.itertuples(index=False), 1):
        base = {"gene_id": r.gene_id, "variant_id": r.variant_id,
                "eQTL_slope": r.slope, "eQTL_p": r.pval_nominal}
        hit = _gwas.lookup(store, str(r.chrom), int(r.pos),
                           str(r.ref).upper(), str(r.alt).upper())
        if hit is None:
            rows.append(dict(base, in_gwas=False))
            continue
        beta_g, se_g, p_g = hit
        rows.append(dict(base, in_gwas=True, gwas_beta=beta_g, gwas_se=se_g,
                         gwas_p=p_g))
        if k % 500 == 0:
            print(f"      {k}/{len(leads)} ({time.time() - t0:.0f}s)",
                  flush=True)

    # reindex first: a run with no matched lead would otherwise have no
    # gwas_beta column at all, and the flanking statements would raise.
    df = pd.DataFrame(rows).reindex(columns=SCREEN_COLS)
    df["in_gwas"] = df["in_gwas"].fillna(False).astype(bool)
    df["allele_mismatch"] = df["in_gwas"] & df["gwas_beta"].isna()

    sig = df[df["in_gwas"] & (df["gwas_p"] < 0.05)]
    dest = PATHS.tables / f"screen_{args.tissue}.tsv"
    df.to_csv(dest, sep="\t", index=False)

    print(f"[3/4] {len(df):,} genes -> {dest.name}", flush=True)
    print(f"[4/4] lead variant in outcome: {int(df['in_gwas'].sum()):,}", flush=True)
    print(f"      outcome p<0.05: {len(sig):,}", flush=True)
    if len(sig):
        concordant = int((sig["eQTL_slope"] * sig["gwas_beta"] > 0).sum())
        p_dir = stats.binomtest(concordant, len(sig), 0.5).pvalue
        print(f"      same direction: {concordant}/{len(sig)} "
              f"({100 * concordant / len(sig):.1f}%), binomial p={p_dir:.4g}",
              flush=True)
    print(f"      total {time.time() - t0:.0f}s", flush=True)


# --------------------------------------------------------------------------
# coloc
# --------------------------------------------------------------------------
def cmd_coloc(args: argparse.Namespace) -> None:
    """Region-level colocalisation for every gene of one tissue.

    Parameters
    ----------
    args : parsed command line.

    Notes
    -----
    The null distribution, when requested, permutes the expression effects
    inside each gene's own variant set and asks how often the resulting PP.H4
    exceeds 0.8. It is an operating-characteristic check on the screen, not a
    p-value for any single gene; the per-gene calibration is 12.
    """
    tissue, _qtd = resolve_tissue(args.tissue)
    t0 = time.time()

    groups = gene_variants(tissue)
    keys = list(groups.groups.keys())
    if args.limit:
        keys = keys[:args.limit]
    print(f"[1/3] {tissue}: {len(keys):,} genes", flush=True)

    store = load_outcome(args.gwas)
    records: List[Dict[str, object]] = []
    tables: Dict[str, pd.DataFrame] = {}
    for i, gene in enumerate(keys, 1):
        d = build_snp_table(groups.get_group(gene), store)
        if d is None or len(d) < args.min_snps:
            continue
        try:
            res = coloc_one(d)
        except AssertionError:
            continue
        res["gene_id"] = gene
        records.append(res)
        tables[gene] = d
        if i % 400 == 0:
            print(f"      {i}/{len(keys)} ({time.time() - t0:.0f}s)", flush=True)

    if not records:
        raise SystemExit("no gene passed the variant filter; check --gwas and "
                         "--min-snps")
    res = pd.DataFrame(records).reindex(columns=COLOC_COLS)
    dest = PATHS.tables / f"coloc_{args.tissue}.tsv"
    res.to_csv(dest, sep="\t", index=False)

    print(f"[2/3] {len(res):,} genes -> {dest.name}", flush=True)
    for thr in (0.5, PARAMS.pp_h4_high, 0.9):
        print(f"      PP.H4 > {thr:g}: {int((res['PP.H4'] > thr).sum()):,}",
              flush=True)
    print(f"      PP.H3 > {PARAMS.pp_h3_high:g}: "
          f"{int((res['PP.H3'] > PARAMS.pp_h3_high).sum()):,}", flush=True)
    print(f"      PP.H4 median {res['PP.H4'].median():.3f}; "
          f"n_snps median {res['n_snps'].median():.0f}", flush=True)

    if args.null:
        print(f"[3/3] permutation null: {args.null} replicates", flush=True)
        rng = np.random.default_rng(SEED)
        sample = list(tables.values())[:min(300, len(tables))]
        rates = []
        for _ in range(args.null):
            hits = 0
            for d in sample:
                dd = d.copy()
                dd["b_e"] = rng.permutation(dd["b_e"].values)
                dd["se_e"] = rng.permutation(dd["se_e"].values)
                hits += coloc_one(dd)["PP.H4"] > PARAMS.pp_h4_high
            rates.append(hits / max(1, len(sample)))
        print(f"      fraction with PP.H4>{PARAMS.pp_h4_high:g}: "
              f"mean {np.mean(rates):.4f}, 95th percentile "
              f"{np.quantile(rates, 0.95):.4f}", flush=True)
        pd.Series(rates).to_csv(dest.with_suffix(".null.csv"), index=False)
    else:
        print("[3/3] permutation null skipped (--null 0)", flush=True)

    print(f"      total {time.time() - t0:.0f}s", flush=True)


# --------------------------------------------------------------------------
# replicate
# --------------------------------------------------------------------------
def cmd_replicate(args: argparse.Namespace) -> None:
    """Carry one tissue's lead variants into the GRCh37 outcome.

    Parameters
    ----------
    args : parsed command line.

    Notes
    -----
    Four allele conventions are accepted: identity, swap, complement and reverse
    complement. They are needed here because the lift changes the build of the
    coordinates, not the strand notation of the alleles, and the two releases do
    not spell the alleles the same way. The sign of the effect is flipped for the
    two swapped cases and kept for the two complement cases, so the reported beta
    always refers to the effect allele of the second outcome.
    """
    tissue, _qtd = resolve_tissue(args.tissue)
    t0 = time.time()

    pairs = _gtex.split_variant(_gtex.load_tissue(gtex_tar(), tissue))
    leads = _gtex.lead_per_gene(pairs[pairs["pos"].notna()])
    print(f"[1/4] {len(leads):,} lead variants", flush=True)

    chain_path = PATHS.data / CHAIN_HG38_TO_HG19
    if not chain_path.exists():
        raise SystemExit(f"missing chain file: {chain_path.name}")
    down = _lift.Chain(str(chain_path))

    chrom19, pos19 = [], []
    for r in leads.itertuples(index=False):
        mapped = down.lift(str(r.chrom), int(r.pos))
        chrom19.append(mapped[0] if mapped else None)
        pos19.append(mapped[1] if mapped else None)
    leads = leads.assign(hg19_chrom=chrom19, hg19_pos=pos19)
    n_ok = int(leads["hg19_pos"].notna().sum())
    print(f"      lifted: {n_ok:,}/{len(leads):,} "
          f"({100 * n_ok / len(leads):.1f}%)", flush=True)
    leads = leads[leads["hg19_pos"].notna()].copy()

    # ---- one streaming pass over the second outcome ---------------------
    gwas2 = resolve_input(args.gwas2)
    if not gwas2.exists():
        raise SystemExit(f"missing outcome file: {gwas2}")
    want = set(zip(leads["hg19_chrom"].astype(str),
                   leads["hg19_pos"].astype(int)))
    print(f"[2/4] scanning {gwas2.name} for {len(want):,} positions", flush=True)

    hits: List[pd.DataFrame] = []
    n_rows = 0
    usecols = ["chromosome", "base_pair_location", "effect_allele",
               "other_allele", "beta", "standard_error", "p_value", "variant_id"]
    for chunk in pd.read_csv(gwas2, sep="\t", usecols=usecols,
                             chunksize=1_000_000,
                             dtype={"chromosome": str}):
        n_rows += len(chunk)
        key = list(zip(chunk["chromosome"].str.replace("chr", "", regex=False),
                       chunk["base_pair_location"]))
        mask = [k in want for k in key]
        if any(mask):
            sub = chunk.loc[mask].copy()
            sub["chr19"] = [k[0] for k, m in zip(key, mask) if m]
            sub["pos19"] = [k[1] for k, m in zip(key, mask) if m]
            hits.append(sub)
    g = pd.concat(hits) if hits else pd.DataFrame(columns=usecols)
    print(f"      scanned {n_rows:,} rows, matched {len(g):,}", flush=True)

    # ---- allele-aware join ----------------------------------------------
    rows: List[Dict[str, object]] = []
    for r in leads.itertuples(index=False):
        sel = g[(g["chr19"] == str(r.hg19_chrom))
                & (g["pos19"] == int(r.hg19_pos))]
        base = {"gene_id": r.gene_id, "variant_id": r.variant_id,
                "eQTL_slope": r.slope, "eQTL_p": r.pval_nominal}
        if sel.empty:
            rows.append(dict(base, in_gwas2=False))
            continue

        got = None
        for s in sel.itertuples(index=False):
            ea = str(s.effect_allele).upper()
            oa = str(s.other_allele).upper()
            ref = str(r.ref).upper()
            alt = str(r.alt).upper()
            if (ea, oa) == (alt, ref):
                got = (s.beta, s.standard_error, s.p_value, s.variant_id)
            elif (ea, oa) == (ref, alt):
                got = (-s.beta, s.standard_error, s.p_value, s.variant_id)
            elif (ea, oa) == (alt.translate(COMPLEMENT), ref.translate(COMPLEMENT)):
                got = (s.beta, s.standard_error, s.p_value, s.variant_id)
            elif (ea, oa) == (ref.translate(COMPLEMENT), alt.translate(COMPLEMENT)):
                got = (-s.beta, s.standard_error, s.p_value, s.variant_id)
            if got:
                break

        if got is None:
            rows.append(dict(base, in_gwas2=True, allele_mismatch=True))
        else:
            rows.append(dict(base, in_gwas2=True, rsid=got[3],
                             gwas2_beta=got[0], gwas2_se=got[1],
                             gwas2_p=got[2]))

    df = pd.DataFrame(rows)
    dest = PATHS.tables / f"screen_{args.tissue}_gcst006465.tsv"
    df.to_csv(dest, sep="\t", index=False)
    print(f"[3/4] {len(df):,} genes -> {dest.name}", flush=True)

    # ---- direction concordance ------------------------------------------
    print("[4/4] direction concordance", flush=True)
    if "gwas2_p" in df.columns:
        d = df[df["gwas2_p"].notna()]
        print(f"      aligned: {len(d):,}/{len(df):,}", flush=True)
        for thr in (0.05, 1e-3):
            s = d[d["gwas2_p"] < thr]
            if not len(s):
                continue
            c = int((s["eQTL_slope"] * s["gwas2_beta"] > 0).sum())
            p = stats.binomtest(c, len(s), 0.5).pvalue
            print(f"      outcome p<{thr:g}: n={len(s):4d}  same direction "
                  f"{c}/{len(s)} = {100 * c / len(s):5.1f}%  "
                  f"binomial p={p:.3g}", flush=True)
    print(f"      total {time.time() - t0:.0f}s", flush=True)


# --------------------------------------------------------------------------
# compare
# --------------------------------------------------------------------------
def cmd_compare(args: argparse.Namespace) -> None:
    """Print summary statistics of two region-level result tables.

    Parameters
    ----------
    args : parsed command line (new, old).

    Notes
    -----
    Written to compare the outcome used for the headline analysis with the
    earlier FinnGen endpoint, which is why the default read is asymmetric: the
    new table is read first and its top rows are printed.
    """
    new_path = PATHS.tables / args.new
    old_path = PATHS.tables / args.old
    for path in (new_path, old_path):
        if not path.exists():
            raise SystemExit(f"missing table: {path.name}")

    new = pd.read_csv(new_path, sep="\t")
    old = pd.read_csv(old_path, sep="\t")

    print(f"new {new_path.name}: n={len(new)}  "
          f"H4max={new['PP.H4'].max():.3f}  "
          f"H4>0.5={int((new['PP.H4'] > 0.5).sum())}  "
          f"H3>0.8={int((new['PP.H3'] > PARAMS.pp_h3_high).sum())}  "
          f"best outcome p median={new['best_gwas_p'].median():.2e}  "
          f"n_snps median={new['n_snps'].median():.0f}", flush=True)
    print(f"old {old_path.name}: n={len(old)}  "
          f"H4max={old['PP.H4'].max():.3f}  "
          f"H4>0.5={int((old['PP.H4'] > 0.5).sum())}  "
          f"H3>0.8={int((old['PP.H3'] > PARAMS.pp_h3_high).sum())}  "
          f"best outcome p median={old['best_gwas_p'].median():.2e}", flush=True)

    cols = [c for c in ("locus", "PP.H3", "PP.H4", "n_snps", "best_gwas_p")
            if c in new.columns]
    print("\ntop rows, new table:", flush=True)
    print(new.nlargest(5, "PP.H4")[cols].to_string(index=False), flush=True)


# --------------------------------------------------------------------------
# entry point
# --------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="command", required=True)

    p = sub.add_parser("download", help="fetch one dataset and prepare its inputs")
    p.add_argument("--tissue", default=None, help="GTEx label or alias")
    p.add_argument("--qtd", default=None, help="dataset id for a tissue not in config")
    p.add_argument("--coloc", action="store_true",
                   help="continue into 07_run_colocalisation.py")
    p.set_defaults(func=cmd_download)

    p = sub.add_parser("screen", help="lead-variant association screen")
    p.add_argument("--tissue", default="Uterus")
    p.add_argument("--gwas", default=GWAS_PRIMARY_FILE)
    p.add_argument("--min-af", type=float, default=PARAMS.eqtl_maf_min)
    p.add_argument("--limit", type=int, default=0)
    p.set_defaults(func=cmd_screen)

    p = sub.add_parser("coloc", help="region-level colocalisation screen")
    p.add_argument("--tissue", default="Uterus")
    p.add_argument("--gwas", default=GWAS_PRIMARY_FILE)
    p.add_argument("--min-snps", type=int, default=3)
    p.add_argument("--null", type=int, default=0,
                   help="permutation replicates for the operating characteristic")
    p.add_argument("--limit", type=int, default=0)
    p.set_defaults(func=cmd_coloc)

    p = sub.add_parser("replicate", help="leads in the second outcome")
    p.add_argument("--tissue", default="Uterus")
    p.add_argument("--gwas2", default=GWAS_SECONDARY_FILE)
    p.set_defaults(func=cmd_replicate)

    p = sub.add_parser("compare", help="two result tables side by side")
    p.add_argument("--new", default="coloc_Uterus.tsv")
    p.add_argument("--old", default="coloc_loci_Uterus.tsv")
    p.set_defaults(func=cmd_compare)

    args = ap.parse_args()
    PATHS.ensure()
    args.func(args)


if __name__ == "__main__":
    main()
