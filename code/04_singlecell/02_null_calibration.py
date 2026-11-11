"""
02_null_calibration.py
======================

Purpose
-------
Rebuild the single-cell null distribution so that it is isomorphic with the
statistic it calibrates.

The earlier implementation compared two objects that live on different scales.
The observed statistic was a *difference*,

    score(UP arm) - score(DOWN arm), averaged per sample,
    then malignant mean minus normal mean,

whereas the null was the per-sample *level* of a single random module, with no
up/down split. A contrast was therefore tested against a displacement: the
locations and the widths of the two distributions differ, and the resulting
p_emp carries no interpretation.

This script reconstructs the null from the same differential statistic. Each
size-matched random module is split at random into two arms whose sizes follow
the programme's own up:down ratio (v2: 608:254 before expression filtering),
the same UP - DOWN difference is formed, and the same sample-level contrast
(malignant mean minus normal mean) is taken. The observed value is placed in
that null.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
--programme   programme table indexed by gene with a `direction` column of
              UP/DOWN (default: PATHS.tables / metabolic_program_v2.tsv,
              862 genes)
--modules     size-matched random modules with `module` and `gene`
              (default: PATHS.tables / random_modules_v2.tsv)
--atlas       integrated h5ad, opened in backed mode
              (default: PATHS.data / ec_integrated_v1.h5ad)
--out         destination for the summary table
              (default: PATHS.tables / singlecell_null_matched.tsv)

Outputs
-------
The summary table, one row, columns in this order:
    observed, null_mean, null_sd, null_p95, z, p_upper, p_lower, p_two_sided,
    n_null, n_up, n_down, mwu_p

A companion file with the same stem plus `_values.tsv`, one column
`null_delta`, holding every null value so that the distribution can be
re-plotted or re-summarised without re-running the atlas.

Notes
-----
Observed value. The malignant-minus-normal difference of the programme score is
negative (-0.059): the programme is *lower* in malignant epithelium, opposite to
the programme direction. The sample-level Mann-Whitney p is 0.11 and the
standardised position in the matched null is z = -3.8. This negative result is
the one the main text cites; the one-sided down-tail p is 0.00, i.e. below the
1/100 resolution of a 100-module null, and must be quoted as such rather than
as a precise zero.

Arms. k_up = round(|module| * n_up / (n_up + n_dn)) genes go to the UP arm and
the remainder to the DOWN arm. Module genes that are absent from the atlas are
dropped before splitting, so a module can be smaller than the programme; the
relative arm sizes are preserved instead of discarding the module.

Module matching. ``random_modules_v2.tsv`` holds 100 modules of 862 genes,
size-matched to the 862-gene programme and split in the same up/down ratio. The
v1 module file is matched to the 146-gene v1 programme and is not a substitute:
the observed contrast would then be compared against contrasts over modules four
times smaller.

Guards. A module needs at least 100 detected genes and at least 20 genes per
arm; fewer than 10 usable modules aborts, because a null of that size cannot
support a tail probability. An arm needs at least five genes for its score to
be defined, the same floor as the projection script.

Statistical caveats. p_two_sided doubles the smaller tail and is truncated at
1.0; it is not an exact two-sided test. z standardises with the null SD at
ddof=1. The modules are permutations of the same size-matched gene pool used
for the projection, so the null asks whether the programme's contrast exceeds
what any gene set of that size shows, not whether it exceeds random genomic
draws. The random split uses config.SEED and nothing else, which is the only
literal constant allowed to be fixed in this package.
"""
from __future__ import annotations

import argparse
import pathlib
import sys
from typing import List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import mannwhitneyu

# Repository root on sys.path, so that `config` resolves from any working
# directory rather than only from the repository root.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

from config.config import PARAMS, PATHS, SEED

# Minimum genes for a score, matching the floor used for the real arms.
MIN_ARM_GENES = 5
# A module is only usable for the null if it is large enough to be split.
MIN_MODULE_GENES = 100
# Each arm must keep enough genes to average over something.
MIN_GENES_PER_ARM = 20
# Below this many usable modules, the tail probability is not resolvable.
MIN_USABLE_MODULES = 10


def _parse_args() -> argparse.Namespace:
    """Command-line arguments; every default is resolved from config."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--programme", type=pathlib.Path,
                    default=PATHS.tables / "metabolic_program_v2.tsv",
                    help="programme table indexed by gene, with a direction column")
    ap.add_argument("--modules", type=pathlib.Path,
                    default=PATHS.tables / "random_modules_v2.tsv",
                    help="size-matched random gene modules")
    ap.add_argument("--atlas", type=pathlib.Path,
                    default=PATHS.data / "ec_integrated_v1.h5ad",
                    help="integrated single-cell atlas in h5ad format")
    ap.add_argument("--out", type=pathlib.Path,
                    default=PATHS.tables / "singlecell_null_matched.tsv",
                    help="summary table for the calibration")
    return ap.parse_args()


def arm_score(gl: List[str], X: sp.csr_matrix, mu: np.ndarray, sd: np.ndarray,
              gidx: dict) -> Optional[np.ndarray]:
    """Per-cell standardised mean expression for one arm of a gene list.

    Parameters
    ----------
    gl : genes of the arm.
    X : cells x detected-genes sparse matrix, columns aligned with ``gidx``.
    mu, sd : per-gene mean and floored standard deviation across all cells.
    gidx : gene name -> column index of ``X``.

    Returns
    -------
    ndarray of length n_cells, or None when fewer than five arm genes are
    detected. Genes absent from the atlas are ignored rather than failing.

    Notes
    ------
    The score is mean_j((x_ij - mu_j) / sd_j), rewritten as x.w - c so that one
    sparse product covers every cell. This is the same definition used for the
    real programme, which is what makes the null comparable to the observation.
    """
    gj = [gidx[g] for g in gl if g in gidx]
    if len(gj) < MIN_ARM_GENES:
        return None
    w = 1.0 / (sd[gj] * len(gj))
    c = float(np.sum(mu[gj] / sd[gj]) / len(gj))
    return np.asarray((X[:, gj] @ w).ravel()) - c


def sample_delta(score: np.ndarray, samples: np.ndarray, classes: np.ndarray,
                 epi: np.ndarray) -> Tuple[float, np.ndarray, np.ndarray]:
    """Sample-level malignant-minus-normal difference of a per-cell score.

    Parameters
    ----------
    score : per-cell score over the whole atlas.
    samples, classes : per-cell sample id and sample class, whole atlas.
    epi : positions of the epithelial cells that enter the aggregation.

    Returns
    -------
    delta : float
        Mean of the malignant sample means minus the mean of the normal sample
        means, or NaN if either group is empty.
    mal, nor : the per-sample means behind that difference, returned so that a
        Mann-Whitney test can be run on the same units.

    Notes
    ------
    Cells enter only through their sample mean. Aggregating to samples first is
    what keeps the test from treating thousands of cells from one donor as
    independent replicates.
    """
    ps = (pd.DataFrame({"v": score[epi], "s": samples[epi], "c": classes[epi]})
          .groupby(["s", "c"], observed=True)["v"].mean().dropna().reset_index())
    mal = ps.loc[ps.c == "Malignant", "v"].dropna().values
    nor = ps.loc[ps.c == "Normal", "v"].dropna().values
    if mal.size == 0 or nor.size == 0:
        return float("nan"), mal, nor
    return float(mal.mean() - nor.mean()), mal, nor


def main() -> None:
    args = _parse_args()
    PATHS.ensure()

    prog = pd.read_csv(args.programme, sep="\t", index_col=0)
    rm = pd.read_csv(args.modules, sep="\t")
    up = list(prog.index[prog.direction == "UP"])
    dn = list(prog.index[prog.direction == "DOWN"])
    print(f"Programme: UP {len(up)} / DOWN {len(dn)}; "
          f"random modules {rm.module.nunique()} "
          f"(PARAMS.n_random_modules={PARAMS.n_random_modules})")
    if rm.module.nunique() != PARAMS.n_random_modules:
        print(f"  warning: {rm.module.nunique()} modules supplied, "
              f"{PARAMS.n_random_modules} configured")

    print(f"Reading atlas (backed): {args.atlas}")
    a = ad.read_h5ad(args.atlas, backed="r")
    obs = a.obs[["sample", "sample_class", "coarse_label"]].copy()
    # Epithelium only: the programme is being projected onto the compartment
    # the bulk signature was derived from, and the two labels are compared.
    epi_mask = obs.coarse_label.isin(["MalignantEpithelial", "NormalEpithelial"]).values
    epi = np.where(epi_mask)[0]
    print(f"{a.n_obs:,} cells x {a.n_vars:,} genes; epithelial {len(epi):,}")

    # Sample support required by the stage-4 parameters. Reported, never
    # enforced: the malignant vs normal design is fixed before this script.
    epi_obs = obs.iloc[epi]
    support = (epi_obs.groupby(["sample", "sample_class"], observed=True).size()
               .rename("n_epithelial_cells").reset_index())
    thin = support[support.n_epithelial_cells < PARAMS.sc_min_cells_per_sample]
    print(f"Epithelial cells per sample: min "
          f"{int(support.n_epithelial_cells.min()):,}, max "
          f"{int(support.n_epithelial_cells.max()):,}; {len(thin)} sample(s) "
          f"below PARAMS.sc_min_cells_per_sample={PARAMS.sc_min_cells_per_sample}")
    for _, r in thin.iterrows():
        print(f"  under-supported: {r['sample']} ({r['sample_class']}), "
              f"n={int(r['n_epithelial_cells']):,}")
    for cls, cnt in support.groupby("sample_class", observed=True).size().items():
        flag = "meets" if cnt >= PARAMS.sc_min_samples else "below"
        print(f"  {cls}: {cnt} samples, {flag} "
              f"PARAMS.sc_min_samples={PARAMS.sc_min_samples}")

    # Gene universe: both programme arms plus every module gene, in one pass.
    mods = {k: list(g.gene) for k, g in rm.groupby("module")}
    need = set(up) | set(dn)
    for v in mods.values():
        need |= set(v)
    idx = a.var_names.get_indexer(list(need))
    keep = idx >= 0
    genes_ok = [g for g, k in zip(list(need), keep) if k]
    X = sp.csr_matrix(a[:, genes_ok].to_memory().X)
    mu = np.asarray(X.mean(axis=0)).ravel()
    sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
    sd = np.sqrt(np.maximum(sq - mu ** 2, 1e-8))
    gidx = {g: i for i, g in enumerate(genes_ok)}
    print(f"Detected genes usable for scoring: {len(genes_ok):,}")

    # ---- observed statistic ----------------------------------------------
    up_s, dn_s = arm_score(up, X, mu, sd, gidx), arm_score(dn, X, mu, sd, gidx)
    if up_s is None or dn_s is None:
        raise SystemExit("one programme arm has too few detected genes")
    ir = up_s - dn_s
    n_up = len([g for g in up if g in gidx])
    n_dn = len([g for g in dn if g in gidx])
    print(f"Arms actually scored: UP {n_up} / DOWN {n_dn}")

    samples = obs["sample"].values
    classes = obs["sample_class"].values

    obs_d, mal, nor = sample_delta(ir, samples, classes, epi)
    u = mannwhitneyu(mal, nor, alternative="two-sided")
    print(f"\nObserved (matched statistic): delta={obs_d:+.4f}  "
          f"MWU p={u.pvalue:.4f}")

    # ---- matched null: split every module in the programme's arm ratio ---
    rng = np.random.default_rng(SEED)
    null_vals: List[float] = []
    ratio = n_up / (n_up + n_dn)          # actual UP share of the programme
    for nm, gl in mods.items():
        g = [x for x in gl if x in gidx]
        if len(g) < MIN_MODULE_GENES:
            continue
        # Expression filtering can leave a module smaller than the programme.
        # Splitting it at the programme's ratio keeps the two arms in the same
        # relative proportion, rather than discarding the module entirely.
        k_up = int(round(len(g) * ratio))
        k_dn = len(g) - k_up
        if k_up < MIN_GENES_PER_ARM or k_dn < MIN_GENES_PER_ARM:
            continue
        perm = rng.permutation(len(g))
        a_up = [g[i] for i in perm[:k_up]]
        a_dn = [g[i] for i in perm[k_up:]]
        s_up, s_dn = arm_score(a_up, X, mu, sd, gidx), arm_score(a_dn, X, mu, sd, gidx)
        if s_up is None or s_dn is None:
            continue
        d, _, _ = sample_delta(s_up - s_dn, samples, classes, epi)
        null_vals.append(d)

    null = np.array(null_vals, dtype=float)
    if null.size < MIN_USABLE_MODULES:
        raise SystemExit(f"only {null.size} usable random modules, "
                         f"too few to build a null distribution")
    null_sd = float(null.std(ddof=1))
    p_upper = float((null >= obs_d).mean())
    p_lower = float((null <= obs_d).mean())
    p_two = float(min(1.0, 2 * min(p_upper, p_lower)))
    z = float((obs_d - null.mean()) / null_sd)
    null_p95 = float(np.percentile(null, 95))
    print(f"Matched null ({len(null)} modules): mean {null.mean():+.4f}, "
          f"SD {null_sd:.4f}, 95th percentile {null_p95:+.4f}")
    print(f"Observed z = {z:+.2f}; upper-tail p={p_upper:.3f}, "
          f"lower-tail p={p_lower:.3f}, two-sided p={p_two:.3f}")
    print(f"Sample-level MWU p = {float(u.pvalue):.3f} at alpha = "
          f"{PARAMS.alpha:g} "
          f"({'significant' if float(u.pvalue) < PARAMS.alpha else 'not significant'})")
    if z < 0:
        print("Direction: the programme score is lower in malignant epithelium, "
              "i.e. opposite to the programme direction.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([dict(observed=obs_d, null_mean=null.mean(), null_sd=null_sd,
                       null_p95=null_p95, z=z,
                       p_upper=p_upper, p_lower=p_lower, p_two_sided=p_two,
                       n_null=len(null), n_up=n_up, n_down=n_dn,
                       mwu_p=float(u.pvalue))]).to_csv(
        args.out, sep="\t", index=False)
    values_path = pathlib.Path(str(args.out).replace(".tsv", "_values.tsv"))
    np.savetxt(str(values_path), null, header="null_delta")
    print(f"\nWrote {args.out} and {values_path}")


if __name__ == "__main__":
    main()
