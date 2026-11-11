"""
01_project_programme.py
=======================

Purpose
-------
Project the insulin-resistance (IR) transcriptional programme, the endometrial
cancer (EC) disease signature and a fixed proliferation score onto the
integrated endometrial single-cell atlas, then ask whether the programme score
separates malignant from normal epithelium.

Three questions are answered, in this order:

    1. descriptive   mean score per sample x coarse cell type, all lineages
    2. inferential   sample-level comparison of epithelial scores, malignant
                     versus normal samples
    3. confounding   the same comparison after removing the proliferating half
                     of the epithelial cells, plus the programme-proliferation
                     correlation

The atlas holds 141,316 cells x 36,601 genes from 11 samples: four normal
(NE1-NE4) and seven malignant (HEEC/MEEC/LEEC), assembled from GSE278879
(11 snRNA samples) and GSE251923 (2 scRNA samples).

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
--atlas       integrated h5ad, opened in backed mode
              (default: PATHS.data / ec_integrated_v1.h5ad)
--programme   programme table indexed by gene, with a `direction` column taking
              the values UP or DOWN
              (default: PATHS.tables / metabolic_program_v2.tsv, the 862-gene
              programme: 608 up, 254 down)
--signature   compact disease signature with `gene` and `direction` columns
              (default: PATHS.tables / EC_signature_v2_prolifadj_compact.tsv)
--modules     random modules with `module` and `gene` columns
              (default: PATHS.tables / random_modules_v2.tsv)
--tag         suffix appended to both output file names (default: empty)

Outputs
-------
PATHS.tables/state_program_table<tag>.tsv
    One row per sample x coarse cell type. Columns, in this order:
    sample, sample_class, coarse_label, n, IR, IR_UP, IR_DN, ECsig, Prolif

PATHS.tables/singlecell_IR_test<tag>.tsv
    A single row: delta_IR, null_mean, null_p95, p_emp

Notes
-----
Scoring. For a gene set with detected members j the per-cell score is
mean_j((x_ij - mu_j) / sd_j), where mu_j and sd_j are the mean and standard
deviation of gene j across every cell of the atlas. It is evaluated as one
sparse product, x.w - c, with w_j = 1 / (sd_j * |set|) and
c = mean_j(mu_j / sd_j), so no dense per-cell intermediate is formed. sd is
floored at 1e-8 to keep near-constant genes from dominating the average. A set
with fewer than five detected members is skipped with a message rather than
scored on a handful of genes.

Sign convention. IR_score = IR_UP - IR_DOWN, so a positive value points along
the insulin-resistant direction of the programme as derived in bulk tissue.
The proliferation score is carried as a covariate only: cycling tumour cells
dominate the epithelial compartment, so an apparent programme shift could be a
cell-cycle artefact.

Reported outcome. In the run the manuscript cites, the epithelial programme
score is *lower* in malignant than in normal epithelium, i.e. opposite to the
programme direction: the malignant-minus-normal difference is -0.059, the
sample-level Mann-Whitney p is 0.11 (7 vs 4 samples), and the matched null of
02_null_calibration.py places the observation 3.8 SD below the null mean. That
negative result is what the main text reports; this script must not be read as
a positive projection of the programme onto malignant cells.

Why the empirical p here is not the headline. The random modules scored below
are whole modules, so their null describes the *level* of a random gene-set
score, whereas the observed statistic is a malignant-normal *difference* at the
sample level. The two are not on the same scale, so p_emp is reported only
because the original script reported it. 02_null_calibration.py rebuilds the
null from the differential statistic itself (each random module split into two
arms of the same relative size as the programme) and is the number to quote.

Module matching. ``random_modules_v2.tsv`` holds 100 modules of 862 genes, i.e.
size-matched to the 862-gene programme and split in the same up/down ratio. The
v1 module file is matched to the 146-gene v1 programme instead and must not be
substituted here: a null built from modules four times smaller than the
programme does not calibrate the programme's own score.

Input provenance. The v2 programme is the 862-gene programme shipped with the
package. The earlier analysis tree defaulted this script to the 146-gene v1
programme in w2/results, which gave a positive but non-significant difference
(+0.056, empirical p 0.91); the v2 numbers were produced there through
PROG_PATH / RM_PATH / OUT_TAG environment overrides. Those overrides are
replaced here by config-derived defaults plus the arguments above.
"""
from __future__ import annotations

import argparse
import pathlib
import sys
from typing import Dict, List, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import mannwhitneyu, spearmanr

# Repository root on sys.path, so that `config` resolves from any working
# directory rather than only from the repository root.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

from config.config import PARAMS, PATHS

# Fixed proliferation markers. Used as a covariate for the confounding check,
# never as a readout of the programme itself.
CYCLING: List[str] = [
    "MKI67", "PCNA", "CCNB1", "CCNB2", "CCNA2", "CDK1", "CDC20", "BUB1", "BUB1B",
    "TOP2A", "AURKA", "AURKB", "PLK1", "TPX2", "UBE2C", "RRM2", "TYMS", "MCM2",
    "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "CENPF", "KIF11", "KIF23", "ANLN",
    "PRC1", "ECT2", "NCAPG", "ASPM", "DLGAP5", "CENPE", "SMC4", "NUSAP1",
]


def _parse_args() -> argparse.Namespace:
    """Command-line arguments; every default is resolved from config."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--atlas", type=pathlib.Path,
                    default=PATHS.data / "ec_integrated_v1.h5ad",
                    help="integrated single-cell atlas in h5ad format")
    ap.add_argument("--programme", type=pathlib.Path,
                    default=PATHS.tables / "metabolic_program_v2.tsv",
                    help="programme table indexed by gene, with a direction column")
    ap.add_argument("--signature", type=pathlib.Path,
                    default=PATHS.tables / "EC_signature_v2_prolifadj_compact.tsv",
                    help="compact disease signature table")
    ap.add_argument("--modules", type=pathlib.Path,
                    default=PATHS.tables / "random_modules_v2.tsv",
                    help="size-matched random gene modules")
    ap.add_argument("--tag", default="", help="suffix for the output file names")
    return ap.parse_args()


def build_sets(programme: pathlib.Path, signature: pathlib.Path,
               modules: pathlib.Path) -> Tuple[Dict[str, List[str]], int]:
    """Assemble the gene sets scored on the atlas.

    Parameters
    ----------
    programme : path
        Programme table indexed by gene with a ``direction`` column of UP/DOWN.
    signature : path
        Disease signature with ``gene`` and ``direction`` columns.
    modules : path
        Random modules with ``module`` and ``gene`` columns. They form the
        empirical null and are scored exactly like the real sets.

    Returns
    -------
    sets : dict
        ``{set name: [gene, ...]}`` in a fixed order: IR_UP, IR_DOWN,
        ECsig_UP, ECsig_DOWN, Prolif, then one entry per random module
        (R000, R001, ...).
    n_sets : int
        Number of sets in the dict, including the random modules.

    Notes
    ------
    Gene order inside a set is irrelevant to the score, which is a mean.
    Duplicated genes inside a set would be weighted twice, so the upstream
    tables are expected to be deduplicated.
    """
    prog = pd.read_csv(programme, sep="\t", index_col=0)
    sig = pd.read_csv(signature, sep="\t")
    rm = pd.read_csv(modules, sep="\t")

    sets: Dict[str, List[str]] = {
        "IR_UP": list(prog.index[prog.direction == "UP"]),
        "IR_DOWN": list(prog.index[prog.direction == "DOWN"]),
        "ECsig_UP": list(sig.gene[sig.direction == "UP"]),
        "ECsig_DOWN": list(sig.gene[sig.direction == "DOWN"]),
        "Prolif": CYCLING,
    }
    # Random modules keep their own labels; the leading "R" identifies them
    # later, which is how the null set is separated from the real sets.
    for name, g in rm.groupby("module"):
        sets[name] = list(g.gene)
    return sets, len(sets)


def sample_level_test(frame: pd.DataFrame, col: str,
                      case: str = "Malignant",
                      control: str = "Normal") -> Tuple[float, float, float]:
    """Compare one per-sample quantity between two sample classes.

    Parameters
    ----------
    frame : DataFrame
        One row per sample, already aggregated; must carry ``sample_class``.
    col : column holding the per-sample value, typically a mean score.
    case, control : labels of the compared classes in ``sample_class``.

    Returns
    -------
    delta : float
        Mean difference, case minus control.
    p : float
        Two-sided Mann-Whitney U p value. With four to seven samples per group
        this is the only test appropriate at this level.
    effect : float
        Cohen's d using the pooled within-group standard deviation
        (``ddof=1``). The 1e-12 term guards against a zero denominator.

    Notes
    ------
    The unit of observation is the sample. Testing on cells instead would
    treat thousands of cells from one donor as independent replicates and
    return p values that shrink with sequencing depth rather than with
    biological replication.
    """
    x = frame.loc[frame.sample_class == case, col]
    y = frame.loc[frame.sample_class == control, col]
    delta = x.mean() - y.mean()
    p = mannwhitneyu(x, y).pvalue
    pooled_sd = np.sqrt((x.var(ddof=1) + y.var(ddof=1)) / 2) + 1e-12
    return delta, float(p), float(delta / pooled_sd)


def main() -> None:
    args = _parse_args()
    PATHS.ensure()

    print(f"Reading atlas (backed): {args.atlas}")
    a = ad.read_h5ad(args.atlas, backed="r")
    obs = a.obs[["sample", "sample_class", "coarse_label"]].copy()
    # One class label per sample, so the design can be stated from the data
    # rather than assumed in the code.
    classes = obs.groupby("sample", observed=True)["sample_class"].first()
    print(f"{a.n_obs:,} cells x {a.n_vars:,} genes; {classes.size} samples "
          f"({int((classes == 'Malignant').sum())} malignant, "
          f"{int((classes == 'Normal').sum())} normal)")

    sets, n_sets = build_sets(args.programme, args.signature, args.modules)
    n_rand = len(sets) - 5            # five real sets are added before the modules
    print(f"Gene sets: {n_sets} in total, {n_rand} random modules "
          f"(PARAMS.n_random_modules={PARAMS.n_random_modules})")
    if n_rand != PARAMS.n_random_modules:
        print(f"  warning: {n_rand} modules supplied, "
              f"{PARAMS.n_random_modules} configured")

    # Union of all sets, restricted to what this atlas actually measures.
    genes = sorted({g for v in sets.values() for g in v})
    idx = a.var_names.get_indexer(genes)
    ok = idx >= 0
    genes_ok = [g for g, k in zip(genes, ok) if k]
    print(f"Union of gene sets: {len(genes)} genes, {len(genes_ok)} detected")

    # Pull the union sub-matrix once, in memory; only these columns are needed.
    X = a[:, genes_ok].to_memory().X
    X = sp.csr_matrix(X)
    print(f"Sub-matrix {X.shape}, non-zeros {X.nnz:,}")

    mu = np.asarray(X.mean(axis=0)).ravel()
    sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
    sd = np.sqrt(np.maximum(sq - mu ** 2, 1e-8))
    gidx = {g: i for i, g in enumerate(genes_ok)}

    # score_i(set) = mean_j((x_ij - mu_j) / sd_j)
    #              = sum_j x_ij / (sd_j * |set|) - mean_j(mu_j / sd_j)
    # so each set becomes one column of a sparse weight matrix and all scores
    # are produced by a single product X @ W.
    rows: List[int] = []
    cols: List[int] = []
    vals: List[float] = []
    consts: List[float] = []
    names: List[str] = []
    for nm, gl in sets.items():
        gj = [gidx[g] for g in gl if g in gidx]
        if len(gj) < 5:
            print(f"  skipped {nm}: only {len(gj)} detected member genes")
            continue
        w = 1.0 / (sd[gj] * len(gj))
        rows += gj
        cols += [len(names)] * len(gj)
        vals += list(w)
        consts.append(float(np.sum(mu[gj] / sd[gj]) / len(gj)))
        names.append(nm)
    W = sp.csr_matrix((vals, (rows, cols)), shape=(len(genes_ok), len(names)))
    S = np.asarray((X @ W).toarray()) - np.array(consts)
    scores = pd.DataFrame(S, columns=names, index=np.arange(S.shape[0]))
    print(f"Scored {len(names)} sets")

    # Cell-level frame: metadata plus the composite programme and signature
    # scores. Row order follows the atlas, so the reset index is positional.
    df = obs.reset_index(drop=True)
    df["IR_UP"] = scores.IR_UP.values
    df["IR_DOWN"] = scores.IR_DOWN.values
    df["IR_score"] = df.IR_UP - df.IR_DOWN
    df["ECsig_score"] = scores.ECsig_UP - scores.ECsig_DOWN
    df["Prolif"] = scores.Prolif
    rand_names = [n for n in names if n.startswith("R")]

    # ---- cell-level description: sample x coarse cell type means ----------
    cellmean = (df.groupby(["sample", "sample_class", "coarse_label"], observed=True)
                .agg(n=("IR_score", "size"), IR=("IR_score", "mean"),
                     IR_UP=("IR_UP", "mean"), IR_DN=("IR_DOWN", "mean"),
                     ECsig=("ECsig_score", "mean"), Prolif=("Prolif", "mean"))
                .reset_index())
    out_cell = PATHS.tables / f"state_program_table{args.tag}.tsv"
    cellmean.to_csv(out_cell, sep="\t", index=False)
    print("\n=== Per sample x cell type (epithelial rows) ===")
    # "Epithelial" matches both the normal Epithelial and MalignantEpithelial
    # labels, which is intended: the two classes are compared below.
    ep = cellmean[cellmean.coarse_label.str.contains("Epithelial")]
    print(ep.round(3).to_string(index=False))

    # ---- sample-level test: epithelial score, malignant vs normal ---------
    print("\n=== Sample-level test (epithelial cells, malignant vs normal) ===")
    ep_df = df[df.coarse_label.str.contains("Epithelial")].copy()
    per_sample = ep_df.groupby(["sample", "sample_class"], observed=True).agg(
        n=("IR_score", "size"), IR=("IR_score", "mean"),
        ECsig=("ECsig_score", "mean"), Prolif=("Prolif", "mean")).reset_index()
    print(per_sample.round(3).to_string(index=False))

    # Sample support required by the stage-4 parameters. Reported, never
    # enforced: dropping a sample here would silently change the malignant vs
    # normal design that the sample-level test rests on.
    support = (ep_df.groupby(["sample", "sample_class"], observed=True).size()
               .rename("n_epithelial_cells").reset_index())
    thin = support[support.n_epithelial_cells < PARAMS.sc_min_cells_per_sample]
    print(f"\nEpithelial cells per sample: min {int(support.n_epithelial_cells.min()):,}, "
          f"max {int(support.n_epithelial_cells.max()):,}; "
          f"{len(thin)} sample(s) below "
          f"PARAMS.sc_min_cells_per_sample={PARAMS.sc_min_cells_per_sample}")
    for _, r in thin.iterrows():
        print(f"  under-supported: {r['sample']} ({r['sample_class']}), "
              f"n={int(r['n_epithelial_cells']):,}")
    for cls, cnt in support.groupby("sample_class", observed=True).size().items():
        flag = "meets" if cnt >= PARAMS.sc_min_samples else "below"
        print(f"  {cls}: {cnt} samples, {flag} "
              f"PARAMS.sc_min_samples={PARAMS.sc_min_samples}")

    rows_out = []
    for label, col in [("All epithelium: IR", "IR"),
                       ("All epithelium: disease signature", "ECsig"),
                       ("All epithelium: proliferation", "Prolif")]:
        d, p, es = sample_level_test(per_sample, col)
        rows_out.append({"comparison": label,
                         "delta_malignant_minus_normal": round(d, 3),
                         "mwu_p": round(p, 3),
                         "cohens_d": round(es, 2)})

    # Non-proliferating epithelium: keep the lower half of the proliferation
    # distribution and repeat the test, so that any residual difference cannot
    # be attributed to cycling cells.
    thr = ep_df.Prolif.median()
    quies = ep_df[ep_df.Prolif <= thr]
    qps = quies.groupby(["sample", "sample_class"], observed=True).agg(
        n=("IR_score", "size"), IR=("IR_score", "mean"),
        ECsig=("ECsig_score", "mean"), Prolif=("Prolif", "mean")).reset_index()
    d_q, p_q, es_q = sample_level_test(qps, "IR")
    rows_out.append({"comparison": "Non-proliferating epithelium: IR",
                     "delta_malignant_minus_normal": round(d_q, 3),
                     "mwu_p": round(p_q, 3),
                     "cohens_d": round(es_q, 2)})
    print("\n" + pd.DataFrame(rows_out).to_string(index=False))

    # ---- confounding with proliferation ----------------------------------
    # Cell-level Spearman correlation, descriptive: it quantifies how far the
    # programme score co-varies with the cell cycle in the epithelial cells.
    rho, p_rho = spearmanr(ep_df.IR_score, ep_df.Prolif)
    print(f"\nIR score vs proliferation in epithelial cells: rho={rho:.3f} "
          f"(p={p_rho:.1e}, n={len(ep_df):,} cells)")

    # ---- random-module null at the sample level --------------------------
    # Each random module is reduced to the same malignant-minus-normal
    # difference. Note that this is a *level* statistic, not the programme's
    # differential statistic; see the module docstring.
    null_vals: List[float] = []
    for nm in rand_names:
        ep_df["_s"] = scores[nm].values[ep_df.index]
        ps = (ep_df.groupby(["sample", "sample_class"], observed=True)["_s"]
              .mean().reset_index())
        x = ps.loc[ps.sample_class == "Malignant", "_s"].mean()
        y = ps.loc[ps.sample_class == "Normal", "_s"].mean()
        null_vals.append(x - y)
    null = np.array(null_vals, dtype=float)
    d_ir = sample_level_test(per_sample, "IR")[0]
    p_ir = sample_level_test(per_sample, "IR")[1]
    null_p95 = float(np.percentile(null, 95))
    p_emp = float((null >= d_ir).mean())
    print(f"\nRandom modules ({len(null)} sets), sample-level difference: "
          f"mean {null.mean():+.4f}, 95th percentile {null_p95:+.4f}, "
          f"max {null.max():+.4f}")
    print(f"Observed IR difference {d_ir:+.4f}; empirical p = {p_emp:.3f}")
    print(f"Sample-level MWU p = {p_ir:.3f} at alpha = {PARAMS.alpha:g} "
          f"({'significant' if p_ir < PARAMS.alpha else 'not significant'})")

    out_test = PATHS.tables / f"singlecell_IR_test{args.tag}.tsv"
    pd.DataFrame({"delta_IR": [d_ir], "null_mean": [null.mean()],
                  "null_p95": [null_p95], "p_emp": [p_emp]}).to_csv(
        out_test, sep="\t", index=False)
    print(f"\nWrote {out_cell} and {out_test}")


if __name__ == "__main__":
    main()
