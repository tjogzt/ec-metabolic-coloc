#!/usr/bin/env python3
"""
03_score_across_cancers.py
==========================

Purpose
-------
Score the insulin-resistance metabolic programme in every TCGA cancer type and
test whether that score tracks two reference axes: proliferation and adipose
differentiation.

The question is whether the programme is a general property of dividing tumour
cells rather than of metabolic tissue. Two comparisons answer it: the score
correlates with a fixed proliferation module, and the same correlation survives
a null built from random modules that are matched to the programme in size and
in expression level.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
<PATHS.cache>/expr.npz, <PATHS.cache>/meta.tsv   written by 02
results/tables/program_symbols.tsv               written by 01
config.PARAMS.n_random_partitions                null-module count
config.PARAMS.min_samples_tumour / _normal      minimum n per block
config.PARAMS.alpha                              significance level used when
                                                 tallying empirical p-values
config.MODULES                                   fixed proliferation / adipose genes

Outputs
-------
results/tables/correlation.tsv
    One row per cancer type and set (tumour block first, then normal block):
    set, cancer_type, n, rho_prolif, lo_p, hi_p, rho_adipo, lo_a, hi_a.
results/tables/null.tsv
    One row per tumour cancer type: observed rho, null mean/sd/p95, z and the
    one-sided empirical p, for both the proliferation and adipose axes.
results/tables/summary.json
    Counts, the median rho, the fraction of cancer types above the coupling
    threshold, and the programme coverage.

Usage
-----
    python3 code/06_cross_cancer/03_score_across_cancers.py

Notes
-----
* Scoring is done inside each cancer type: every gene is z-transformed across
  the samples of that cancer type alone, then
  score = mean(z over UP genes) - mean(z over DOWN genes). No cross-cancer
  normalisation happens anywhere, so a correlation cannot be manufactured by
  comparing cancer types that differ in baseline expression.
* The z-transformation is computed once per cancer type and reused by all
  random modules. Recomputing it inside the null loop would make the stage
  O(n x 100) instead of O(n) and is the reason a previous implementation ran
  for hours.
* Random modules are equal in size to the real programme (separately for the UP
  and DOWN arms) and are drawn from expression deciles that reproduce the
  programme's own decile occupancy. Without the size and expression match the
  null is not comparable to the observed score.
* The normal-tissue control uses the solid-tissue normals shipped inside the
  same TCGA matrix, never GTEx: same platform, same processing, so the contrast
  isolates tissue state rather than resource.
* Gene symbols are the match key throughout. The matrix is symbol-keyed; an
  Entrez-based join returns an empty set, which is how an earlier version of
  this stage silently produced no results.
* ``min_samples_tumour`` (50) gates the tumour block and the null loop;
  ``min_samples_normal`` (20) gates the solid-tissue-normal block. The two
  thresholds differ, as they did in the original analysis.
  alike. The pre-port analysis used 50 tumours and 20 normals; the single
  configured floor replaces both, which is the only intended numerical
  difference from the source script.
* The bootstrap resampling count (1000) is fixed here rather than read from
  config: config.PARAMS.n_bootstrap belongs to the stage-3 burden analysis and
  is a different quantity.
"""
from __future__ import annotations

import json
import pathlib
import sys
import time
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from config.config import MODULES, PARAMS, PATHS, SEED  # noqa: E402

# Resamples used for the confidence interval around each Spearman rho.
N_BOOTSTRAP = 1000

# Fraction of cancer types with rho above this value that counts as "coupling is
# universal"; below it, heterogeneity has to be reported instead.
COUPLING_RHO = 0.4
COUPLING_FRACTION = 0.8

PROG_TABLE = "program_symbols.tsv"


def zmat(X: np.ndarray, cols: List[int]) -> np.ndarray:
    """Z-transform every gene across the samples of one cancer type.

    Parameters
    ----------
    X : genes x samples matrix.
    cols : sample column indices belonging to this cancer type.

    Returns
    -------
    genes x ``len(cols)`` array of per-cancer z-scores.

    Notes
    -----
    A gene with zero variance inside the cancer type would divide by zero; its
    standard deviation is replaced by 1.0, which turns that gene into a constant
    and contributes nothing to any score. This is the only special case.
    """
    sub = X[:, cols]
    mu = sub.mean(axis=1, keepdims=True)
    sd = sub.std(axis=1, keepdims=True)
    sd[sd == 0] = 1.0
    return (sub - mu) / sd


def scores(Z: np.ndarray, gene_idx: Dict[str, int], up: List[str],
           down: List[str], prolif: List[str],
           adipo: List[str]) -> Dict[str, Optional[np.ndarray]]:
    """Programme score and the two reference-axis scores, per sample.

    The programme score is the mean z of the UP arm minus the mean z of the DOWN
    arm, so that the two directions are combined with opposite signs instead of
    being averaged into a cancellation.

    Parameters
    ----------
    Z : z-scored genes x samples matrix for one cancer type.
    gene_idx : gene symbol -> row index in ``Z``.
    up, down : programme arms, already restricted to genes present in the matrix.
    prolif, adipo : reference modules, same restriction.

    Returns
    -------
    dict with keys ``prog``, ``prolif`` and ``adipo``; ``adipo`` is None when no
    adipose gene survived the restriction.
    """
    def arm(genes: List[str]) -> np.ndarray:
        return Z[[gene_idx[g] for g in genes]].mean(axis=0)

    return {
        "prog": arm(up) - arm(down),
        "prolif": arm(prolif),
        "adipo": arm(adipo) if adipo else None,
    }


def boot_ci(a: np.ndarray, b: np.ndarray, rng: np.random.Generator,
            n: int = N_BOOTSTRAP):
    """Spearman rho with a resampling confidence interval.

    Samples are drawn with replacement and the correlation recomputed; draws
    that yield a non-finite rho (possible when a resample is constant) are
    dropped before the percentile.

    Returns
    -------
    (rho, lower 2.5% bound, upper 97.5% bound)
    """
    r = stats.spearmanr(a, b).statistic
    bs = []
    N = len(a)
    for _ in range(n):
        k = rng.integers(0, N, N)
        v = stats.spearmanr(a[k], b[k]).statistic
        if np.isfinite(v):
            bs.append(v)
    lo, hi = np.percentile(bs, [2.5, 97.5]) if bs else (np.nan, np.nan)
    return r, lo, hi


def run_set(label: str, min_n: int, sub_meta: pd.DataFrame, X: np.ndarray,
            gene_idx: Dict[str, int], col_of: Dict[str, int], rng: np.random.Generator,
            up: List[str], down: List[str], prolif: List[str],
            adipo: List[str]) -> pd.DataFrame:
    """Correlate the programme score with both reference axes, one row per cancer type.

    Cancer types with fewer than ``min_n`` samples are skipped: below that floor
    a Spearman estimate is dominated by sampling noise. Because cancer types are
    visited in sorted order and the random generator is shared, the bootstrap
    draws are reproducible only if the call order of this function is preserved.
    """
    rows = []
    for ct, grp in sub_meta.groupby("cancer_type"):
        cols = [col_of[s] for s in grp.index if s in col_of]
        if len(cols) < min_n:
            continue
        # One z-transformation per cancer type, used for every score below.
        Z = zmat(X, cols)
        sc = scores(Z, gene_idx, up, down, prolif, adipo)
        r_p, lo_p, hi_p = boot_ci(sc["prog"], sc["prolif"], rng)
        if sc["adipo"] is not None:
            r_a, lo_a, hi_a = boot_ci(sc["prog"], sc["adipo"], rng)
        else:
            r_a = lo_a = hi_a = np.nan
        rows.append(dict(set=label, cancer_type=ct, n=len(cols),
                         rho_prolif=r_p, lo_p=lo_p, hi_p=hi_p,
                         rho_adipo=r_a, lo_a=lo_a, hi_a=hi_a))
        print(f"  [{label}] {ct:45s} n={len(cols):5d} "
              f"rho_prolif={r_p:+.3f} rho_adipo={r_a:+.3f}", flush=True)
    return pd.DataFrame(rows)


def build_background(X: np.ndarray, genes: np.ndarray, gene_idx: Dict[str, int],
                     prog_set: set) -> tuple:
    """Split the non-programme transcriptome into ten expression deciles.

    The deciles are computed on the pooled sample mean of each gene and define
    the pool from which random modules are drawn. Matching on decile rather than
    on exact expression keeps the pool large enough to sample from while still
    removing the dominant confounder: highly expressed genes carry more
    measurement precision and different variance.

    Returns
    -------
    (pool_by_decile, programme_decile_distribution)
    """
    pool = np.array([g for g in genes if g not in prog_set])
    pool_idx = np.array([gene_idx[g] for g in pool])
    pool_mean = X[pool_idx].mean(axis=1)
    deciles = pd.qcut(pool_mean, 10, labels=False, duplicates="drop")
    pool_by_dec = {d: pool[deciles == d] for d in np.unique(deciles)}

    # Where the real programme sits in that same decile system: this histogram,
    # not a uniform draw, is the target the random modules must reproduce.
    prog_mean = X[[gene_idx[g] for g in prog_set]].mean(axis=1)
    edges = np.quantile(pool_mean, np.linspace(0, 1, 11))[1:-1]
    prog_decile_dist = np.bincount(np.digitize(prog_mean, edges), minlength=10)
    return pool_by_dec, prog_decile_dist


def draw_module(pool_by_dec: Dict[int, np.ndarray], want: Dict[int, int],
                n_up: int, n_dn: int, rng: np.random.Generator):
    """Draw one size- and expression-matched random module and split it.

    The draw takes ``want[d]`` genes from decile ``d`` -- exactly the number the
    real programme contributes to that decile -- then shuffles the pool and cuts
    it into an UP arm of ``n_up`` and a DOWN arm of ``n_dn``. Size and split
    therefore match the real programme by construction, and the only thing that
    differs between the observed and the null score is which genes are used.

    Returns
    -------
    (up genes, down genes) or None when the pool could not supply enough genes.
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

    npz_path = PATHS.cache / "expr.npz"
    meta_path = PATHS.cache / "meta.tsv"
    prog_path = PATHS.tables / PROG_TABLE
    for p in (npz_path, meta_path, prog_path):
        if not p.exists():
            raise SystemExit(f"missing input: {p}\nrun 01 and 02 first")

    # Tumour and normal blocks were filtered at different sizes in the
    # original analysis; keeping them separate preserves the published set.
    min_tumour = PARAMS.min_samples_tumour
    min_normal = PARAMS.min_samples_normal
    n_null = PARAMS.n_random_partitions
    alpha = PARAMS.alpha

    # ---- cached matrix and metadata -------------------------------------
    z = np.load(npz_path, allow_pickle=True)
    X = z["X"]
    genes = np.array([str(g) for g in z["genes"]])
    samples = np.array([str(s) for s in z["samples"]])
    gene_idx = {g: i for i, g in enumerate(genes)}
    meta = pd.read_csv(meta_path, sep="\t").set_index("sample_id")
    col_of = {s: i for i, s in enumerate(samples)}

    # ---- programme arms and reference modules ---------------------------
    pmap = pd.read_csv(prog_path, sep="\t")
    n_prog_total = len(pmap)
    if n_prog_total == 0:
        raise SystemExit(f"programme table is empty: {prog_path}")

    # Symbol matching: the matrix is keyed by HUGO symbol, so this is the only
    # join that can succeed.
    pmap = pmap[pmap.SYMBOL.isin(gene_idx)]
    up = [g for g in pmap.loc[pmap.direction == "UP", "SYMBOL"] if g in gene_idx]
    down = [g for g in pmap.loc[pmap.direction == "DOWN", "SYMBOL"] if g in gene_idx]
    prolif = [g for g in MODULES["proliferation"] if g in gene_idx]
    adipo = [g for g in MODULES["adipose"] if g in gene_idx]

    n_scored = len(up) + len(down)
    print(f"Programme scoreable: {n_scored}/{n_prog_total} "
          f"(UP {len(up)} / DOWN {len(down)}), coverage "
          f"{100 * n_scored / n_prog_total:.1f}%", flush=True)

    missing_prolif = set(MODULES["proliferation"]) - set(prolif)
    print(f"Proliferation module {len(prolif)}/{len(MODULES['proliferation'])}"
          f"{' (missing ' + ','.join(sorted(missing_prolif)) + ')' if missing_prolif else ''}",
          flush=True)
    missing_adipo = sorted(set(MODULES["adipose"]) - set(adipo))
    print(f"Adipose module {len(adipo)}/{len(MODULES['adipose'])}, "
          f"missing {','.join(missing_adipo) or 'none'}", flush=True)

    if not up or not down:
        raise SystemExit("one programme arm has no scoreable gene; check the symbol match")

    # ---- decile-matched background --------------------------------------
    prog_set = set(up) | set(down)
    pool_by_dec, prog_decile_dist = build_background(X, genes, gene_idx, prog_set)
    print("Programme genes by expression decile:", prog_decile_dist.tolist(), flush=True)

    # ---- main scan and normal-tissue control -----------------------------
    print("\n=== Main analysis: primary tumours ===", flush=True)
    res_t = run_set("tumor", min_tumour,
                    meta[meta.sample_type == "Primary Tumor"], X, gene_idx, col_of,
                    rng, up, down, prolif, adipo)

    print("\n=== Control: TCGA solid tissue normals ===", flush=True)
    res_n = run_set("normal", min_normal,
                    meta[meta.sample_type == "Solid Tissue Normal"], X, gene_idx, col_of,
                    rng, up, down, prolif, adipo)

    pd.concat([res_t, res_n]).to_csv(PATHS.tables / "correlation.tsv",
                                     sep="\t", index=False)

    # ---- null distribution ----------------------------------------------
    print(f"\n=== Null distribution ({n_null} size-matched random modules) ===", flush=True)
    want = {d: int(v) for d, v in enumerate(prog_decile_dist) if v > 0}
    n_up, n_dn = len(up), len(down)
    null_rows = []

    for ct, grp in meta[meta.sample_type == "Primary Tumor"].groupby("cancer_type"):
        cols = [col_of[s] for s in grp.index if s in col_of]
        if len(cols) < min_tumour:
            continue
        # Z-scored once and reused by every random module below: this is what
        # keeps the null loop at O(n) instead of O(n x n_null).
        Z = zmat(X, cols)
        sc = scores(Z, gene_idx, up, down, prolif, adipo)
        obs = stats.spearmanr(sc["prog"], sc["prolif"]).statistic

        null_rhos, null_rhos_a = [], []
        for _ in range(n_null):
            drawn = draw_module(pool_by_dec, want, n_up, n_dn, rng)
            if drawn is None:
                continue
            ru, rd = drawn
            # Same arithmetic as the observed score, on the matched random arms.
            sc_r = Z[[gene_idx[g] for g in ru]].mean(0) - Z[[gene_idx[g] for g in rd]].mean(0)
            null_rhos.append(stats.spearmanr(sc_r, sc["prolif"]).statistic)
            if sc["adipo"] is not None:
                null_rhos_a.append(stats.spearmanr(sc_r, sc["adipo"]).statistic)

        if not null_rhos:
            continue
        null_rhos = np.array(null_rhos)
        zz = (obs - null_rhos.mean()) / null_rhos.std(ddof=1)
        # One-sided empirical p with the +1 correction of Davison & Hinkley, so
        # an observed value above every null draw never reports p = 0.
        p_emp = (np.sum(null_rhos >= obs) + 1) / (len(null_rhos) + 1)

        obs_a = (stats.spearmanr(sc["prog"], sc["adipo"]).statistic
                 if sc["adipo"] is not None else np.nan)
        if null_rhos_a:
            null_rhos_a = np.array(null_rhos_a)
            zz_a = (obs_a - null_rhos_a.mean()) / null_rhos_a.std(ddof=1)
            p_a = (np.sum(null_rhos_a >= obs_a) + 1) / (len(null_rhos_a) + 1)
            na, nsa, p95a = (null_rhos_a.mean(), null_rhos_a.std(ddof=1),
                             np.percentile(null_rhos_a, 95))
        else:
            zz_a = p_a = na = nsa = p95a = np.nan

        null_rows.append(dict(cancer_type=ct, n=len(cols), obs=obs,
                              null_mean=null_rhos.mean(), null_sd=null_rhos.std(ddof=1),
                              z=zz, p_emp=p_emp, null_p95=np.percentile(null_rhos, 95),
                              obs_adipo=obs_a, null_mean_adipo=na, null_sd_adipo=nsa,
                              z_adipo=zz_a, p_emp_adipo=p_a, null_p95_adipo=p95a))
        print(f"  {ct:45s} prolif obs={obs:+.3f} z={zz:+.2f} p={p_emp:.3f} | "
              f"adipo obs={obs_a:+.3f} z={zz_a:+.2f} p={p_a:.3f}", flush=True)

    if null_rows:
        pd.DataFrame(null_rows).to_csv(PATHS.tables / "null.tsv", sep="\t", index=False)

    # ---- verdict ---------------------------------------------------------
    if len(res_t):
        r = res_t.rho_prolif.values
        frac = float(np.mean(r > COUPLING_RHO))
        n_null_sig = sum(1 for row in null_rows if row["p_emp"] < alpha)
        print("\n=== Main criterion ===")
        print(f"Tumour cancer types {len(r)}; fraction with rho>{COUPLING_RHO} {frac:.1%}; "
              f"median rho {np.median(r):+.3f}")
        if len(res_n):
            print(f"Normal tissues {len(res_n)}; median rho "
                  f"{np.median(res_n.rho_prolif):+.3f}")
        print(f"Empirical p < {alpha} on the proliferation axis: {n_null_sig}/{len(null_rows)}"
              f" cancer types")
        print(f"Verdict: {'universal coupling' if frac >= COUPLING_FRACTION else 'not universal (heterogeneity must be reported)'}")

        json.dump({"n_tumor": len(res_t), "n_normal": len(res_n),
                   "frac_rho_gt_0.4": frac,
                   "median_rho_tumor": float(np.median(r)),
                   "median_rho_normal": (float(np.median(res_n.rho_prolif))
                                         if len(res_n) else None),
                   "prog_coverage": f"{n_scored}/{n_prog_total}",
                   "alpha": alpha,
                   "n_p_emp_below_alpha": n_null_sig},
                  open(PATHS.tables / "summary.json", "w"), indent=2)

    print(f"Done in {time.time() - t0:.0f}s.", flush=True)


if __name__ == "__main__":
    main()
