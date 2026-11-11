"""
12_permutation_null.py
======================

Purpose
-------
Calibrate the region-level colocalisation screen by permutation, gene by gene.

Running the ABF test on the *significant* variant set of each gene inflates
PP.H4: every variant kept is an expression-associated variant, so an incidental
fluctuation in the outcome is more likely to fall on the same variants. The
inflation is not a property of the data but of the truncated variant set, and it
differs between genes, so a raw PP.H4 cannot be compared across genes.

The calibration permutes the outcome effects *within each gene's own variant
set*, which keeps the marginal distribution of the outcome effects and the
number of variants fixed, and asks where the observed PP.H4 sits in that
gene-specific null. The resulting p-value is calibrated conditional on the
variant set, which is exactly the comparison the screen needs.

Two null modes are offered:

    shift   circular shift of the outcome vector along genomic position. The
            spatial clustering of the outcome signal (its LD blocks) survives;
            only its pairing with the expression effects is broken.
    perm    full permutation. Destroys the clustering as well, which makes the
            null too wide and the test conservative.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
--tissue     GTEx tissue label whose significant pairs define the genes
             (default: Uterus)
--gwas       outcome file in data/ (default: the primary outcome; pass the
             FinnGen endpoint file to reproduce the exploratory run)
--B          permutations per gene (default: 50)
--min-snps   minimum variants for a gene to enter (default: 5)
--null-mode  shift (default) or perm
--limit      test only the first N genes, for a timing check
--seed       RNG seed (default: config.SEED)

Outputs
-------
results/tables/g_a_candidates.tsv, one row per gene:

    gene_id, n_snps, PP.H4_obs, null_mean, null_q95, excess, p_perm,
    best_gwas_p, q_BH

p_perm = (1 + #{null >= observed}) / (B + 1), so it can never be zero; q_BH is
the Benjamini-Hochberg adjustment over the genes actually tested.

Notes
-----
This is a permutation test conditional on a restricted variant set. Its power is
lower than a standard coloc run on all tested variants, and its null hypothesis
is "no pairing between the expression and outcome effects on this variant set",
not "no shared causal variant". The Methods section has to say both, otherwise
the p-values will be read as if they came from a standard colocalisation test.

The output is a candidate ranking, not a set of colocalisation calls: a gene
that survives here still has to be tested on the full variant set by 07.
"""
from __future__ import annotations

import argparse
import pathlib
import sys
import time
from typing import Dict, List

import numpy as np
import pandas as pd

# Local package imports; the directory holding this file is on sys.path.
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))       # repository root, for config/

from config.config import PATHS, PARAMS, TISSUES, GWAS_PRIMARY_FILE, SEED
from importlib import import_module

_abf = import_module("05_coloc_abf_engine")
_gwas = import_module("02_gwas_reader")
_gtex = import_module("13_gtex_v11_reader")
_screen = import_module("15_locus_screen")

OUT_COLS = ["gene_id", "n_snps", "PP.H4_obs", "null_mean", "null_q95", "excess",
            "p_perm", "best_gwas_p", "q_BH"]


def observed_and_null(df: pd.DataFrame, n_perm: int, mode: str,
                      rng: np.random.Generator) -> Dict[str, float]:
    """Observed PP.H4 and its permutation distribution for one gene.

    Parameters
    ----------
    df : aligned variant table from 15_locus_screen.build_snp_table, sorted by
        variant id so that the shift is stable across calls.
    n_perm : permutations.
    mode : "shift" or "perm".
    rng : generator shared across genes, so one seed reproduces the whole run.

    Returns
    -------
    dict with PP.H4_obs, null_mean, null_q95, excess and p_perm.
    """
    b_e = df["b_e"].to_numpy(dtype=float)
    se_e = df["se_e"].to_numpy(dtype=float)
    b_g = df["b_g"].to_numpy(dtype=float)
    se_g = df["se_g"].to_numpy(dtype=float)

    # The expression side never changes: only the pairing is permuted.
    l1 = _abf.log_abf(b_e, se_e, 1.0)
    obs = _abf.coloc_abf(l1, _abf.log_abf(b_g, se_g, 1.0),
                         p1=PARAMS.p1, p2=PARAMS.p2, p12=PARAMS.p12)["PP.H4"]

    null = np.empty(n_perm, dtype=float)
    n = len(b_g)
    for i in range(n_perm):
        if mode == "shift":
            # Circular shift: a rotation, never the identity, so every replicate
            # is a genuine re-pairing while the outcome's run structure survives.
            k = int(rng.integers(1, n))
            bg, sg = np.roll(b_g, k), np.roll(se_g, k)
        else:
            idx = rng.permutation(n)
            bg, sg = b_g[idx], se_g[idx]
        null[i] = _abf.coloc_abf(l1, _abf.log_abf(bg, sg, 1.0),
                                 p1=PARAMS.p1, p2=PARAMS.p2,
                                 p12=PARAMS.p12)["PP.H4"]

    # The +1 in both terms keeps the p-value strictly positive and unbiased for a
    # finite number of permutations.
    p_perm = (1.0 + float((null >= obs - 1e-12).sum())) / (n_perm + 1.0)
    return {"PP.H4_obs": float(obs),
            "null_mean": float(null.mean()),
            "null_q95": float(np.quantile(null, 0.95)),
            "excess": float(obs) - float(null.mean()),
            "p_perm": p_perm}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tissue", default="Uterus")
    ap.add_argument("--gwas", default=GWAS_PRIMARY_FILE)
    ap.add_argument("--B", type=int, default=50, dest="n_perm")
    ap.add_argument("--min-snps", type=int, default=5)
    ap.add_argument("--null-mode", default="shift", choices=["shift", "perm"],
                    help="shift keeps the outcome's LD clustering; perm does not")
    ap.add_argument("--limit", type=int, default=0,
                    help="test only the first N genes")
    ap.add_argument("--seed", type=int, default=SEED)
    args = ap.parse_args()

    PATHS.ensure()
    t0 = time.time()

    gwas_path = pathlib.Path(args.gwas)
    if not gwas_path.is_absolute():
        gwas_path = PATHS.data / gwas_path
    if not gwas_path.exists():
        raise SystemExit(f"missing outcome file: {gwas_path}")

    tar_path = PATHS.data / "GTEx_Analysis_v8_eQTL.tar"
    if not tar_path.exists():
        raise SystemExit(f"missing GTEx archive: {tar_path.name}")

    # ---- 1. genes and their variant sets --------------------------------
    pairs = _gtex.split_variant(_gtex.load_tissue(tar_path, args.tissue))
    pairs = pairs[pairs["pos"].notna()]
    genes = pairs.groupby("gene_id")
    keys = list(genes.groups.keys())
    if args.limit:
        keys = keys[:args.limit]
    print(f"[1/2] {len(keys):,} genes, {args.n_perm} permutations per gene "
          f"({args.null_mode} null)", flush=True)

    store = _gwas.load_gwas(gwas_path)
    print(f"      outcome: {sum(len(v) for v in store.values()):,} variants",
          flush=True)

    rng = np.random.default_rng(args.seed)
    rows: List[Dict[str, object]] = []
    for i, gene in enumerate(keys, 1):
        d = _screen.build_snp_table(genes.get_group(gene), store)
        if d is None or len(d) < args.min_snps:
            continue
        # Sorted by variant id so that a repeated run permutes the same rows.
        d = d.sort_values("variant_id").reset_index(drop=True)
        row = {"gene_id": gene, "n_snps": len(d),
               "best_gwas_p": float(d["p_g"].min())}
        row.update(observed_and_null(d, args.n_perm, args.null_mode, rng))
        rows.append(row)
        if i % 400 == 0:
            print(f"      {i}/{len(keys)} ({time.time() - t0:.0f}s)", flush=True)

    if not rows:
        raise SystemExit("no gene produced a testable variant set; check the "
                         "tissue label and --min-snps")

    # ---- 2. Benjamini-Hochberg over the genes tested --------------------
    res = pd.DataFrame(rows).sort_values("p_perm").reset_index(drop=True)
    m = len(res)
    res["q_BH"] = (res["p_perm"] * m / (np.arange(m) + 1)).cummin().clip(upper=1)

    out = res[OUT_COLS]
    dest = PATHS.tables / "g_a_candidates.tsv"
    out.to_csv(dest, sep="\t", index=False)

    print(f"[2/2] {m:,} genes -> {dest.name}", flush=True)
    print(f"      p<0.05        {int((out['p_perm'] < 0.05).sum()):,}", flush=True)
    print(f"      p<0.05/m      {int((out['p_perm'] < 0.05 / m).sum()):,}",
          flush=True)
    print(f"      q<0.10        {int((out['q_BH'] < 0.10).sum()):,}", flush=True)
    print(f"      q<0.05        {int((out['q_BH'] < 0.05).sum()):,}", flush=True)
    print(f"      total {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
