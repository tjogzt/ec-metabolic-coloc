#!/usr/bin/env python3
"""
04_normal_tissue_null.py
========================

Purpose
-------
Build the same decile-matched null distribution for the solid-tissue normals
that stage 03 builds for the tumours, so that a correlation seen in normal
tissue can be compared against a null of its own.

Normal-tissue groups are small (tens of samples), and a Spearman rho estimated
from tens of samples is wide. Without a null of matching size and expression
content, a positive rho in a normal group could not be attributed to the
programme rather than to the small sample.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
<PATHS.cache>/expr.npz, <PATHS.cache>/meta.tsv   written by 02
results/tables/program_symbols.tsv               written by 01
config.PARAMS.n_random_partitions                null-module count
config.PARAMS.min_samples_normal                minimum group size
config.PARAMS.alpha                              significance level for the tally
config.MODULES                                   fixed proliferation genes

Outputs
-------
results/tables/null_normal.tsv
    One row per normal-tissue group:
    cancer_type, n, obs, null_mean, null_sd, null_p95, z, p_emp.

Usage
-----
    python3 code/06_cross_cancer/04_normal_tissue_null.py

Notes
-----
* Only the proliferation axis is tested here. Adipose-tissue markers are close
  to absent in most normal epithelia, and the adipose contrast is settled on the
  tumour side in stage 03.
* The null is built exactly as in stage 03: modules equal in size to the real
  programme, drawn from the same expression deciles, split into the same UP and
  DOWN arm sizes. The decile background is estimated once from the whole matrix
  and then held fixed, so tumour and normal nulls are directly comparable.
* The decile-matching block is repeated from stage 03 on purpose instead of
  being imported: each stage script stays runnable on its own, and the shared
  logic is short enough that duplication is cheaper than a module whose only
  consumer is these two entry points. Any change to the matching rule must be
  applied to both files.
* Normal samples come from the TCGA matrix itself, not from GTEx; see stage 02.
"""
from __future__ import annotations

import sys
import pathlib
import time
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from config.config import MODULES, PARAMS, PATHS, SEED  # noqa: E402

PROG_TABLE = "program_symbols.tsv"


def zmat(X: np.ndarray, cols: List[int]) -> np.ndarray:
    """Z-transform every gene across the samples of one group.

    Zero-variance genes inside the group get a standard deviation of 1.0, which
    makes them constant and therefore silent in every score.
    """
    sub = X[:, cols]
    mu = sub.mean(axis=1, keepdims=True)
    sd = sub.std(axis=1, keepdims=True)
    sd[sd == 0] = 1.0
    return (sub - mu) / sd


def build_background(X: np.ndarray, genes: np.ndarray, gene_idx: Dict[str, int],
                     prog_set: set) -> tuple:
    """Split the non-programme transcriptome into ten expression deciles.

    Returns the decile pools and the decile occupancy of the real programme,
    which is the histogram the random modules have to reproduce.
    """
    pool = np.array([g for g in genes if g not in prog_set])
    pool_mean = X[[gene_idx[g] for g in pool]].mean(axis=1)
    deciles = pd.qcut(pool_mean, 10, labels=False, duplicates="drop")
    pool_by_dec = {d: pool[deciles == d] for d in np.unique(deciles)}

    prog_mean = X[[gene_idx[g] for g in prog_set]].mean(axis=1)
    edges = np.quantile(pool_mean, np.linspace(0, 1, 11))[1:-1]
    dist = np.bincount(np.digitize(prog_mean, edges), minlength=10)
    return pool_by_dec, dist


def draw_module(pool_by_dec: Dict[int, np.ndarray], want: Dict[int, int],
                n_up: int, n_dn: int, rng: np.random.Generator):
    """Draw one size- and expression-matched module and cut it into two arms.

    ``want[d]`` genes are taken from decile ``d``; the resulting pool is shuffled
    and split into ``n_up`` and ``n_dn`` genes. Matching both the per-decile
    counts and the arm sizes is what makes the null score comparable to the
    observed one.

    Returns
    -------
    (up genes, down genes), or None when the pools fall short of the required
    size.
    """
    draw = []
    for d, k in want.items():
        cand = pool_by_dec.get(d, np.array([]))
        if len(cand) == 0:
            continue
        draw.extend(rng.choice(cand, size=min(k, len(cand)), replace=False).tolist())
    if len(draw) < n_up + n_dn:
        return None
    rng.shuffle(draw)
    return draw[:n_up], draw[n_up:n_up + n_dn]


def main() -> None:
    PATHS.ensure()
    t0 = time.time()
    rng = np.random.default_rng(SEED)

    min_samples = PARAMS.min_samples_normal
    n_null = PARAMS.n_random_partitions
    alpha = PARAMS.alpha

    npz_path = PATHS.cache / "expr.npz"
    meta_path = PATHS.cache / "meta.tsv"
    prog_path = PATHS.tables / PROG_TABLE
    for p in (npz_path, meta_path, prog_path):
        if not p.exists():
            raise SystemExit(f"missing input: {p}\nrun 01 and 02 first")

    z = np.load(npz_path, allow_pickle=True)
    X = z["X"]
    genes = np.array([str(g) for g in z["genes"]])
    samples = np.array([str(s) for s in z["samples"]])
    gene_idx = {g: i for i, g in enumerate(genes)}
    meta = pd.read_csv(meta_path, sep="\t").set_index("sample_id")
    col_of = {s: i for i, s in enumerate(samples)}

    # ---- programme arms and the proliferation axis ----------------------
    pmap = pd.read_csv(prog_path, sep="\t")
    pmap = pmap[pmap.SYMBOL.isin(gene_idx)]
    up = [g for g in pmap.loc[pmap.direction == "UP", "SYMBOL"] if g in gene_idx]
    down = [g for g in pmap.loc[pmap.direction == "DOWN", "SYMBOL"] if g in gene_idx]
    prolif = [g for g in MODULES["proliferation"] if g in gene_idx]

    if not up or not down or not prolif:
        raise SystemExit("an arm or the proliferation module is empty; check the symbol match")

    # ---- matched background, estimated once -----------------------------
    prog_set = set(up) | set(down)
    pool_by_dec, dist = build_background(X, genes, gene_idx, prog_set)
    want = {d: int(v) for d, v in enumerate(dist) if v > 0}
    n_up, n_dn = len(up), len(down)
    print(f"Programme arms: UP {n_up} / DOWN {n_dn}; decile occupancy {dist.tolist()}",
          flush=True)

    # ---- null per normal-tissue group -----------------------------------
    norm = meta[meta.sample_type == "Solid Tissue Normal"]
    if norm.empty:
        raise SystemExit("no solid tissue normal sample in the phenotype table")

    rows = []
    for ct, grp in norm.groupby("cancer_type"):
        cols = [col_of[s] for s in grp.index if s in col_of]
        if len(cols) < min_samples:
            continue

        # Z-scored once per group and reused by all random modules below, so the
        # null loop costs O(n) rather than O(n x n_null).
        Z = zmat(X, cols)
        sp = Z[[gene_idx[g] for g in up]].mean(0) - Z[[gene_idx[g] for g in down]].mean(0)
        pr = Z[[gene_idx[g] for g in prolif]].mean(0)
        obs = stats.spearmanr(sp, pr).statistic

        nulls = []
        for _ in range(n_null):
            drawn = draw_module(pool_by_dec, want, n_up, n_dn, rng)
            if drawn is None:
                continue
            ru, rd = drawn
            # Same arithmetic as the observed score, on the matched random arms.
            sc = (Z[[gene_idx[g] for g in ru]].mean(0)
                  - Z[[gene_idx[g] for g in rd]].mean(0))
            nulls.append(stats.spearmanr(sc, pr).statistic)

        if not nulls:
            continue
        nulls = np.array(nulls)
        rows.append(dict(cancer_type=ct, n=len(cols), obs=obs, null_mean=nulls.mean(),
                         null_sd=nulls.std(ddof=1), null_p95=np.percentile(nulls, 95),
                         z=(obs - nulls.mean()) / nulls.std(ddof=1),
                         p_emp=(np.sum(nulls >= obs) + 1) / (len(nulls) + 1)))
        print(f"  {ct:42s} n={len(cols):4d} obs={obs:+.3f} "
              f"null={nulls.mean():+.3f}+/-{nulls.std(ddof=1):.3f} "
              f"z={rows[-1]['z']:+.2f} p={rows[-1]['p_emp']:.3f}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(PATHS.tables / "null_normal.tsv", sep="\t", index=False)

    if len(df):
        print(f"\nNormal groups {len(df)}; groups with obs above null p95 "
              f"{(df.obs > df.null_p95).sum()}; groups with p<{alpha} "
              f"{(df.p_emp < alpha).sum()}; median z {df.z.median():+.2f}")
    else:
        print("\nNo normal-tissue group reached the minimum size; table written empty.")
    print(f"Done in {time.time() - t0:.0f}s.", flush=True)


if __name__ == "__main__":
    main()
