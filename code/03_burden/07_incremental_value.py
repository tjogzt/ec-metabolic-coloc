"""
07_incremental_value.py
=======================

Purpose
-------
Ask whether the insulin-resistance programme adds predictive information about
lymph-node metastasis in endometrioid carcinoma, over and above a base clinical
model, using an endpoint and a cohort that are independent of the score's
derivation.

Four steps, previously four scripts, now one file because they share the matrix
decoder, the scoring rule and the cross-validation scheme -- and because they
must agree on all three for the comparison to mean anything:

    main         base model vs base + programme vs base + disease signature,
                 with expression-matched random modules as the null. The
                 disease signature is the positive control: it must move the
                 AUC if the apparatus works at all.
    adjusted     the same comparison after regressing the score on the
                 per-sample proliferation score, i.e. asking whether anything
                 survives once cell-cycle content is removed.
    bootstrap    stratified bootstrap confidence interval for delta-AUC, and the
                 smallest increment the design could have detected.
    diagnostics  univariate AUCs and correlations with the clinical covariates,
                 used to explain a null result rather than assert one.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
PATHS.data/GSE120490_series_matrix.txt.gz   145 tumours, 45 metastasis events
PATHS.data/GPL570_annot.tsv                 probe-to-symbol map
PATHS.tables/IR_program_v1.tsv              programme (``--programme`` overrides)
PATHS.tables/random_modules_v1.tsv          null modules (``--modules`` overrides)
PATHS.tables/incremental_null<tag>.tsv      written by --step main, read by
                                            --step bootstrap and --step adjusted

config.SEED                    seed for the CV splits and the bootstrap
config.PARAMS.n_bootstrap      bootstrap resamples
config.PARAMS.program_de_padj  not used; kept out of this stage deliberately

Outputs
-------
PATHS.tables/incremental_value<tag>.tsv      point estimates, one row
PATHS.tables/incremental_null<tag>.tsv       one row per random module
PATHS.tables/incremental_adjusted.tsv        proliferation-residualised result
PATHS.tables/incremental_boot<tag>.tsv       bootstrap intervals

Usage
-----
    python code/03_burden/07_incremental_value.py                    # main
    python code/03_burden/07_incremental_value.py --step all
    python code/03_burden/07_incremental_value.py --step bootstrap \
        --programme metabolic_program_v2.tsv --modules random_modules_v2.tsv \
        --tag _v2

Notes
-----
The endpoint decoder is shared by all four steps. The original scripts differed:
the main step normalised the label before comparing it and the other three
compared it exactly. One decoder is used here, so the four steps cannot silently
disagree about which samples are events.

The size guard in ``--step main`` is not a formality. A null module of a
different size than the programme has a different variance, and delta-AUC would
then partly measure set size; the run stops instead.

Seeds: the CV splits use ``config.SEED``; each bootstrap replicate offsets it so
that replicates differ from one another while the whole run stays reproducible.

Unused imports and intermediate arrays from the original four scripts (a
point-biserial helper, one single-score design matrix that was never consumed)
are not carried over.
"""
from __future__ import annotations

import argparse
import gzip
import io
import pathlib
import sys
from importlib import import_module
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import (RepeatedStratifiedKFold, StratifiedKFold,
                                     cross_val_score)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

# Repository root for `config`; this directory for the sibling stages.
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))

from config.config import PARAMS, PATHS, SEED

# The cell-cycle panel is defined once, in the stage that audits the signature.
CYCLING_GENES = import_module("05_build_programme").CYCLING_GENES

# --------------------------------------------------------------------------
# Module constants
# --------------------------------------------------------------------------
# Cohort that carries the endpoint. 145 endometrioid carcinomas, 45 with lymph
# node metastasis; the score comes from a different cohort entirely.
ENDPOINT_SERIES = "GSE120490_series_matrix.txt.gz"

# Clinical covariates of the base model, as named in the series matrix. The
# spelling "matastasis" is the one in the GEO record.
FIELD_EVENT = "matastasis"
FIELD_GRADE = "grade"
FIELD_MYOMETRIAL = "mi50pluspercent"

# Cross-validation scheme.
CV_SPLITS = 5
CV_REPEATS = 20               # point estimate, as pre-registered
CV_REPEATS_BOOT = 10          # bootstrap replicates, for run time only
BOOT_SEED_BASE = 1000         # replicate seeds are SEED + this offset + index

# Model settings. Fixed for every comparison so that delta-AUC reflects the
# feature block, not a change of algorithm.
LOGREG_MAX_ITER = 2000
LOGREG_C = 1.0
CV_DESCRIPTION = f"{CV_SPLITS}-fold x{CV_REPEATS} repeats (bootstrap), " \
                 f"point estimate {CV_SPLITS}x{CV_REPEATS}"

# A gene set needs at least this many measured members to be scored; below it
# the mean is dominated by whichever genes happen to be present.
MIN_GENES_IN_SET = 5

# Percentiles of the bootstrap distribution reported as the interval.
BOOT_PERCENTILES = (2.5, 97.5)

# z-multiplier giving roughly 80% power at a two-sided 5% level; used to convert
# the null spread into a detection floor.
DETECTABLE_Z = 2.8

# Pre-registered minimum increment for the axis to be called positive. This is a
# decision threshold, not a significance level, so it is not PARAMS.alpha.
PREREGISTERED_DELTA_AUC = 0.05


# --------------------------------------------------------------------------
# Matrix decoding and scoring
# --------------------------------------------------------------------------
def read_matrix(path: pathlib.Path) -> Tuple[pd.DataFrame, Dict[str, List[List[str]]]]:
    """Read a GEO series matrix and its sample metadata.

    Parameters
    ----------
    path : gzip-compressed series matrix.

    Returns
    -------
    (probes x samples DataFrame, metadata dict keyed by the ``!Sample_*`` field)

    Notes
    -----
    The table sits between two marker lines; everything before the first marker
    is sample metadata, where one field may occupy several lines (one per value
    group). Each field maps to a list of rows, and each row to the values of one
    sample, so that a field can be recovered by row.
    """
    with gzip.open(path, "rt") as fh:
        lines = fh.readlines()
    i0 = next(i for i, l in enumerate(lines)
              if l.startswith("!series_matrix_table_begin"))
    i1 = next(i for i, l in enumerate(lines)
              if l.startswith("!series_matrix_table_end"))
    m = pd.read_csv(io.StringIO("".join(lines[i0 + 1:i1])),
                    sep="\t", quotechar='"', index_col=0)

    meta: Dict[str, List[List[str]]] = {}
    for line in lines[:i0]:
        if line.startswith("!Sample_"):
            key, _, value = line.partition("\t")
            meta.setdefault(key, []).append(
                [x.strip('"') for x in value.rstrip("\n").split("\t")])
    return m, meta


def field(meta: Dict[str, List[List[str]]], name: str) -> List[str]:
    """Extract one sample field, matching the row label by prefix.

    Parameters
    ----------
    meta : metadata dict from :func:`read_matrix`.
    name : field name, matched case-insensitively against the label prefix.

    Returns
    -------
    list of values, one per sample.

    Raises
    ------
    KeyError
        If no row starts with the name, so that a renamed field fails loudly
        instead of yielding an empty covariate.
    """
    chars = meta["!Sample_characteristics_ch1"]
    for row in chars:
        if row[0].lower().startswith(name.lower()):
            # "label: value" -> value; a value without a label passes through.
            return [x.split(": ", 1)[-1] for x in row]
    raise KeyError(name)


def collapse_to_gene(m: pd.DataFrame) -> pd.DataFrame:
    """Collapse probes to gene symbols, keeping the highest-mean probe.

    Parameters
    ----------
    m : probes x samples matrix.

    Returns
    -------
    genes x samples DataFrame.

    Notes
    -----
    Highest mean expression, not highest variance: the latter would select the
    noisiest probe and favour genes with a large between-sample spread, which is
    the quantity the model is trying to predict.
    """
    annot = pd.read_csv(PATHS.data / "GPL570_annot.tsv", sep="\t")
    symbols = annot.set_index("ID")["Gene Symbol"]

    m = m.loc[m.index.intersection(symbols.index)]
    m = m.astype(float)
    mean = m.mean(axis=1)

    df = pd.DataFrame({"sym": symbols.reindex(m.index).values, "mean": mean.values})
    df["idx"] = np.arange(len(df))          # positional index into the matrix
    df = df[~df.sym.isin(["---", ""])].dropna(subset=["sym"])

    best = df.sort_values("mean", ascending=False).drop_duplicates("sym")
    out = m.iloc[best.idx.values]
    out.index = best.sym.values
    return out


def zscore(m: pd.DataFrame) -> pd.DataFrame:
    """Z-score each gene across samples.

    Parameters
    ----------
    m : genes x samples matrix.

    Returns
    -------
    Matrix with zero mean and unit variance per row; constant rows become NaN,
    which excludes them from every mean rather than contributing a zero.
    """
    return m.sub(m.mean(axis=1), axis=0).div(m.std(axis=1).replace(0, np.nan),
                                             axis=0)


def score_set(z: pd.DataFrame, genes: List[str],
              up: Optional[List[str]] = None) -> Optional[pd.Series]:
    """Score one gene set per sample.

    Parameters
    ----------
    z : gene x sample z-scored matrix.
    genes : gene set to score.
    up : genes of the up arm; when given, the score is mean(up) - mean(rest),
         otherwise it is the plain mean of the set.

    Returns
    -------
    Series of per-sample scores, or None when fewer than ``MIN_GENES_IN_SET``
    members are measured or one arm is empty. None is returned rather than a
    degenerate score so that the caller can drop the module.
    """
    present = [g for g in genes if g in z.index]
    if len(present) < MIN_GENES_IN_SET:
        return None
    if up is None:
        return z.loc[present].mean(axis=0)

    up_present = [g for g in up if g in z.index]
    down_present = [g for g in present if g not in up_present]
    if not up_present or not down_present:
        return None
    return z.loc[up_present].mean(axis=0) - z.loc[down_present].mean(axis=0)


def cv_auc(X: np.ndarray, y: np.ndarray, seed: int = SEED,
           n_repeats: int = CV_REPEATS) -> float:
    """Repeated stratified cross-validated AUC of a logistic model.

    Parameters
    ----------
    X : n_samples x n_features design matrix.
    y : binary outcome.
    seed : seed for the split generator.
    n_repeats : number of repeats of the 5-fold split.

    Returns
    -------
    Mean AUC over all folds and repeats.

    Notes
    -----
    Scaling is inside the pipeline, so it is fitted on the training fold only.
    Scaling before splitting would leak the test-fold distribution into the
    model and inflate every AUC by a small, hard-to-detect amount.
    """
    cv = RepeatedStratifiedKFold(n_splits=CV_SPLITS, n_repeats=n_repeats,
                                 random_state=seed)
    aucs = []
    for train, test in cv.split(X, y):
        model = make_pipeline(
            StandardScaler(),
            LogisticRegression(max_iter=LOGREG_MAX_ITER, C=LOGREG_C))
        model.fit(X[train], y[train])
        aucs.append(roc_auc_score(y[test], model.predict_proba(X[test])[:, 1]))
    return float(np.mean(aucs))


def _load_endpoint() -> Tuple[np.ndarray, np.ndarray, np.ndarray, pd.DataFrame]:
    """Build the endpoint and the two clinical covariates.

    Returns
    -------
    (y, grade, myometrial, z) where z is the z-scored gene x sample matrix; the
    last element is computed here because every step needs it.

    Notes
    -----
    Labels are normalised before comparison so that "Yes", "yes" and " Yes" all
    decode as an event. A silently mis-decoded label would flip the meaning of
    every AUC in the file, which is the one error no downstream check can catch.
    """
    raw, meta = read_matrix(PATHS.data / ENDPOINT_SERIES)
    met = field(meta, FIELD_EVENT)
    grade = field(meta, FIELD_GRADE)
    mi = field(meta, FIELD_MYOMETRIAL)

    y = np.array([1 if v.strip().lower() == "yes" else 0 for v in met])
    grade_i = np.array([1 if v == "Yes" else 0 for v in grade])
    mi_i = np.array([1 if v == "Yes" else 0 for v in mi])

    m = collapse_to_gene(raw)
    z = zscore(m)
    return y, grade_i, mi_i, z


def _base_design(prolif: pd.Series, grade: np.ndarray,
                 myometrial: np.ndarray) -> np.ndarray:
    """Assemble the base clinical design matrix.

    Parameters
    ----------
    prolif : per-sample proliferation score.
    grade : high-grade indicator.
    myometrial : myometrial invasion >= 50% indicator.

    Returns
    -------
    n_samples x 3 matrix: proliferation, grade, myometrial invasion.
    """
    return np.column_stack([prolif.values, grade, myometrial])


def _read_programme(name: str) -> pd.DataFrame:
    """Read a programme table from the tables directory."""
    path = PATHS.tables / name
    if not path.exists():
        raise SystemExit(f"missing programme: {path}")
    return pd.read_csv(path, sep="\t", index_col=0)


def _read_modules(name: str) -> pd.DataFrame:
    """Read a random-module table from the tables directory."""
    path = PATHS.tables / name
    if not path.exists():
        raise SystemExit(f"missing random modules: {path}")
    return pd.read_csv(path, sep="\t")


def _write(df: pd.DataFrame, stem: str, tag: str, index: bool = False) -> pathlib.Path:
    """Write a result table as ``<stem><tag>.tsv`` in the tables directory."""
    dest = PATHS.tables / f"{stem}{tag}.tsv"
    df.to_csv(dest, sep="\t", index=index)
    return dest


# --------------------------------------------------------------------------
# Step: main comparison
# --------------------------------------------------------------------------
def step_main(programme: str, modules: str, tag: str) -> None:
    """Compare the base model with each score block against the null modules.

    Parameters
    ----------
    programme : file name of the programme table inside ``PATHS.tables``.
    modules : file name of the random-module table inside ``PATHS.tables``.
    tag : suffix appended to the output file names.

    Returns
    -------
    None. Two tables are written.
    """
    PATHS.ensure()
    y, grade, mi, z = _load_endpoint()
    print(f"samples {len(y)}; events {y.sum()} ({100 * y.mean():.1f}%); "
          f"high grade {grade.sum()}; myometrial >=50% {mi.sum()}")

    prog = _read_programme(programme)
    sig = pd.read_csv(PATHS.tables / "EC_signature_v2_prolifadj_compact.tsv",
                      sep="\t")
    print(f"gene-level matrix {z.shape[0]:,} genes x {z.shape[1]} samples")

    up_genes = list(prog.index[prog.direction == "UP"])
    ir = score_set(z, list(prog.index), up_genes)
    ecsig = score_set(z, list(sig.gene), list(sig.gene[sig.direction == "UP"]))
    # A score is None when one arm of the set has no measured gene. The original
    # scripts dereferenced the result directly and died with an AttributeError;
    # failing here names the cause instead.
    for label, score in (("programme", ir), ("disease signature", ecsig)):
        if score is None:
            raise SystemExit(
                f"{label} cannot be scored on this cohort: fewer than "
                f"{MIN_GENES_IN_SET} genes measured, or one direction arm empty")
    prolif = z.loc[[g for g in CYCLING_GENES if g in z.index]].mean(axis=0)
    print(f"programme genes measured "
          f"{sum(g in z.index for g in prog.index)}/{len(prog)}; "
          f"signature {sum(g in z.index for g in sig.gene)}/{len(sig)}")

    base = _base_design(prolif, grade, mi)
    auc_base = cv_auc(base, y)
    print(f"\nbase (proliferation + grade + myometrial invasion) CV AUC = "
          f"{auc_base:.3f}")

    d_ir = cv_auc(np.column_stack([base, ir.values]), y) - auc_base
    d_sig = cv_auc(np.column_stack([base, ecsig.values]), y) - auc_base
    print(f"  + programme         dAUC = {d_ir:+.4f}")
    print(f"  + disease signature dAUC = {d_sig:+.4f}  <- positive control")

    # ---- null distribution ----------------------------------------------
    rm = _read_modules(modules)
    rm_sizes = rm.groupby("module").size()
    prog_size = len(prog.index)
    if not (rm_sizes == prog_size).all():
        raise SystemExit(
            f"random modules and programme differ in size: modules "
            f"{sorted(rm_sizes.unique())} vs programme {prog_size} "
            f"({modules}). Check that --modules matches --programme.")

    null: List[float] = []
    used: List[str] = []
    for name, g in rm.groupby("module"):
        s = score_set(z, list(g.gene))
        if s is None:
            continue
        null.append(cv_auc(np.column_stack([base, s.values]), y) - auc_base)
        used.append(name)
    if not null:
        raise SystemExit(
            "no random module could be scored on this cohort; the modules and "
            "the endpoint matrix share too few genes, so there is no null "
            f"distribution to compare against ({modules})")
    null_arr = np.array(null)
    p_ir = float((null_arr >= d_ir).mean())
    p_sig = float((null_arr >= d_sig).mean())
    print(f"\nrandom modules {len(null_arr)}: dAUC mean {null_arr.mean():+.4f}, "
          f"95th percentile {np.percentile(null_arr, 95):+.4f}, "
          f"max {null_arr.max():+.4f}")
    print(f"empirical p (programme) = {p_ir:.3f}")
    print(f"empirical p (signature) = {p_sig:.3f}")

    _write(pd.DataFrame({"module": used, "delta_auc": null_arr}),
           "incremental_null", tag)
    _write(pd.DataFrame([{
        "auc_base": auc_base,
        "delta_auc_IR": d_ir,
        "delta_auc_ECsig": d_sig,
        "null_mean": null_arr.mean(),
        "null_p95": np.percentile(null_arr, 95),
        "p_IR": p_ir,
        "p_ECsig": p_sig,
        "n": int(len(y)),
        "events": int(y.sum()),
    }]), "incremental_value", tag)
    print(f"\nwritten incremental_value{tag}.tsv")


# --------------------------------------------------------------------------
# Step: proliferation-adjusted
# --------------------------------------------------------------------------
def _residualise(s: pd.Series, prolif: pd.Series) -> np.ndarray:
    """Remove the linear component of the proliferation score from a score.

    Parameters
    ----------
    s : per-sample score.
    prolif : per-sample proliferation score.

    Returns
    -------
    Residual array.

    Notes
    -----
    A linear fit is enough: the relationship between a mean z-score over a
    tumour cohort's genes and a mean z-score over cell-cycle genes is close to
    linear, and a non-linear adjustment would trade interpretability for a
    marginal gain that the sample size cannot resolve.
    """
    design = prolif.values.reshape(-1, 1)
    lr = LinearRegression().fit(design, s.values)
    return s.values - lr.predict(design)


def step_adjusted(programme: str, modules: str, tag: str) -> None:
    """Repeat the comparison after removing proliferation from every score.

    Parameters
    ----------
    programme, modules, tag : as in :func:`step_main`.

    Returns
    -------
    None. One table is written.
    """
    PATHS.ensure()
    y, grade, mi, z = _load_endpoint()
    prolif = z.loc[[g for g in CYCLING_GENES if g in z.index]].mean(axis=0)

    prog = _read_programme(programme)
    up_g = list(prog.index[prog.direction == "UP"])
    dn_g = list(prog.index[prog.direction == "DOWN"])
    scores = {
        "programme": score_set(z, list(prog.index), up_g),
        "up arm": score_set(z, up_g),
        "down arm": score_set(z, dn_g),
    }

    base = np.column_stack([prolif.values, grade, mi])
    auc_base = cv_auc(base, y)
    print(f"base CV AUC = {auc_base:.3f}\n")
    if scores["programme"] is None:
        raise SystemExit("programme cannot be scored on this cohort: fewer than "
                         f"{MIN_GENES_IN_SET} genes measured, or one arm empty")
    print(f"{'score':12s} {'raw AUC':>9s} {'residual AUC':>13s} "
          f"{'dAUC(base+resid)':>18s}")
    for name, s in scores.items():
        if s is None:
            print(f"{name:12s} not scorable (too few genes measured)")
            continue
        # Direction-free AUC: a score that predicts the endpoint inversely is
        # still informative, and the residualised version may flip sign.
        a_raw = max(roc_auc_score(y, s.values), 1 - roc_auc_score(y, s.values))
        resid = _residualise(s, prolif)
        a_res = max(roc_auc_score(y, resid), 1 - roc_auc_score(y, resid))
        d = cv_auc(np.column_stack([base, resid]), y) - auc_base
        print(f"{name:12s} {a_raw:9.3f} {a_res:13.3f} {d:+18.4f}")

    # ---- null distribution, raw and residualised -------------------------
    rm = _read_modules(modules)
    null_raw, null_res = [], []
    for _, g in rm.groupby("module"):
        s = score_set(z, list(g.gene))
        if s is None:
            continue
        null_raw.append(cv_auc(np.column_stack([base, s.values]), y) - auc_base)
        null_res.append(cv_auc(np.column_stack([base, _residualise(s, prolif)]), y)
                        - auc_base)
    null_raw_arr = np.array(null_raw)
    null_res_arr = np.array(null_res)
    print(f"\nrandom modules ({len(null_res_arr)}):")
    print(f"  raw        dAUC mean {null_raw_arr.mean():+.4f}, "
          f"95th percentile {np.percentile(null_raw_arr, 95):+.4f}")
    print(f"  residual   dAUC mean {null_res_arr.mean():+.4f}, "
          f"95th percentile {np.percentile(null_res_arr, 95):+.4f}")

    resid_total = _residualise(scores["programme"], prolif)
    d_ir = cv_auc(np.column_stack([base, resid_total]), y) - auc_base
    p = float((null_res_arr >= d_ir).mean())
    print(f"\nprogramme, residualised: dAUC = {d_ir:+.4f}, empirical p = {p:.2f}")
    verdict = ("failed (pre-registered threshold dAUC >= "
               f"{PREREGISTERED_DELTA_AUC:g})"
               if d_ir < PREREGISTERED_DELTA_AUC else "passed")
    print(f"incremental-value verdict: {verdict}")

    _write(pd.DataFrame([{
        "auc_base": auc_base,
        "delta_auc_IR_resid": d_ir,
        "p_emp": p,
        "null_res_mean": null_res_arr.mean(),
        "null_res_p95": np.percentile(null_res_arr, 95),
    }]), "incremental_adjusted", "")
    print(f"written {PATHS.tables / 'incremental_adjusted.tsv'}")


# --------------------------------------------------------------------------
# Step: bootstrap interval
# --------------------------------------------------------------------------
def step_bootstrap(programme: str, modules: str, tag: str) -> None:
    """Bootstrap the delta-AUC and convert the null spread into a detection floor.

    Parameters
    ----------
    programme, modules, tag : as in :func:`step_main`.

    Returns
    -------
    None. One table is written.

    Notes
    -----
    A null result is only interpretable together with the smallest increment the
    design could have detected. Both are therefore reported on the same row: the
    interval bounds what is excluded, and the detection floor says how far below
    the point estimate the design was blind.

    Resampling is stratified, so every replicate keeps the event fraction; a
    plain bootstrap produces replicates with no events at all, where the AUC is
    undefined.
    """
    PATHS.ensure()
    y, grade, mi, z = _load_endpoint()
    prolif = z.loc[[g for g in CYCLING_GENES if g in z.index]].mean(axis=0)

    prog = _read_programme(programme)
    up_genes = list(prog.index[prog.direction == "UP"])
    ir = score_set(z, list(prog.index), up_genes)
    if ir is None:
        raise SystemExit("programme cannot be scored on this cohort: fewer than "
                         f"{MIN_GENES_IN_SET} genes measured, or one arm empty")

    base = _base_design(prolif, grade, mi)
    X = np.column_stack([base, ir.values])

    print(f"samples {len(y)}; events {y.sum()}; programme {len(prog)} genes "
          f"({sum(g in z.index for g in prog.index)} measured)")
    a0 = cv_auc(base, y)
    a1 = cv_auc(X, y)
    print(f"base AUC = {a0:.4f}; +programme AUC = {a1:.4f}; "
          f"dAUC = {a1 - a0:+.4f}")

    # ---- stratified bootstrap -------------------------------------------
    rng = np.random.default_rng(SEED)
    idx_pos = np.where(y == 1)[0]
    idx_neg = np.where(y == 0)[0]
    deltas, auc_b, auc_f = [], [], []
    for b in range(PARAMS.n_bootstrap):
        take = np.concatenate([rng.choice(idx_pos, len(idx_pos), replace=True),
                               rng.choice(idx_neg, len(idx_neg), replace=True)])
        Xb, yb = X[take], y[take]
        Bb = Xb[:, :base.shape[1]]
        try:
            b0 = cv_auc(Bb, yb, seed=BOOT_SEED_BASE + b, n_repeats=CV_REPEATS_BOOT)
            b1 = cv_auc(Xb, yb, seed=BOOT_SEED_BASE + b, n_repeats=CV_REPEATS_BOOT)
        except ValueError:
            # A replicate with a single class cannot be cross-validated; rare
            # with stratification, but it must not abort the loop.
            continue
        auc_b.append(b0)
        auc_f.append(b1)
        deltas.append(b1 - b0)
        if (b + 1) % 25 == 0:
            print(f"  bootstrap {b + 1}/{PARAMS.n_bootstrap}  "
                  f"running mean dAUC {np.mean(deltas):+.4f}", flush=True)

    d = np.array(deltas)
    lo, hi = np.percentile(d, BOOT_PERCENTILES)
    print(f"\nbootstrap ({len(d)} resamples) dAUC median {np.median(d):+.4f}, "
          f"95% CI [{lo:+.4f}, {hi:+.4f}]")
    print(f"  AUC base     95% CI "
          f"[{np.percentile(auc_b, BOOT_PERCENTILES[0]):.4f}, "
          f"{np.percentile(auc_b, BOOT_PERCENTILES[1]):.4f}]")
    print(f"  AUC +programme 95% CI "
          f"[{np.percentile(auc_f, BOOT_PERCENTILES[0]):.4f}, "
          f"{np.percentile(auc_f, BOOT_PERCENTILES[1]):.4f}]")

    # ---- detection floor, from the size-matched null ---------------------
    null_path = PATHS.tables / f"incremental_null{tag}.tsv"
    if not null_path.exists():
        raise SystemExit(f"missing null distribution: {null_path}\n"
                         f"  run --step main first, with the same --tag")
    null = pd.read_csv(null_path, sep="\t").delta_auc.values
    sd = float(np.std(null, ddof=1))
    z_obs = (a1 - a0 - null.mean()) / sd
    print(f"\nnull distribution ({len(null)} size-matched modules): "
          f"mean {null.mean():+.4f}, SD {sd:.4f}")
    print(f"  observed dAUC, z against null = {z_obs:+.2f}")
    print(f"  80% power detection floor ~ {DETECTABLE_Z} x SD = "
          f"{DETECTABLE_Z * sd:.4f}")
    print(f"  i.e. increments above {hi:.4f} are excluded (bootstrap upper bound)")

    _write(pd.DataFrame([{
        "auc_base": a0,
        "auc_full": a1,
        "delta_auc": a1 - a0,
        "boot_n": len(d),
        "boot_lo": lo,
        "boot_hi": hi,
        "boot_median": float(np.median(d)),
        "auc_base_lo": float(np.percentile(auc_b, BOOT_PERCENTILES[0])),
        "auc_base_hi": float(np.percentile(auc_b, BOOT_PERCENTILES[1])),
        "auc_full_lo": float(np.percentile(auc_f, BOOT_PERCENTILES[0])),
        "auc_full_hi": float(np.percentile(auc_f, BOOT_PERCENTILES[1])),
        "null_sd": sd,
        "null_mean": float(null.mean()),
        "z_vs_null": z_obs,
        "detectable_80pct": DETECTABLE_Z * sd,
        "n": int(len(y)),
        "events": int(y.sum()),
        "cv": CV_DESCRIPTION,
    }]), "incremental_boot", tag)
    print(f"\nwritten {PATHS.tables / f'incremental_boot{tag}.tsv'}")


# --------------------------------------------------------------------------
# Step: diagnostics
# --------------------------------------------------------------------------
def step_diagnostics(programme: str, modules: str, tag: str) -> None:
    """Explain a null result: univariate signal and mediation by covariates.

    Parameters
    ----------
    programme, modules, tag : accepted for a uniform interface; the module table
        is not used here.

    Returns
    -------
    None. Everything is printed.

    Notes
    -----
    Three explanations are distinguishable and lead to different conclusions:
    a univariate AUC near 0.5 means the score carries no information; a high
    univariate AUC that vanishes once grade and invasion enter means the score
    is mediated by them; an AUC between them with a wide interval means the
    study is underpowered. The tables below separate the three.
    """
    del modules, tag                     # unused in this step by design
    y, grade, mi, z = _load_endpoint()
    prog = _read_programme(programme)
    sig = pd.read_csv(PATHS.tables / "EC_signature_v2_prolifadj_compact.tsv",
                      sep="\t")
    scores = {
        "programme": score_set(z, list(prog.index),
                               list(prog.index[prog.direction == "UP"])),
        "signature": score_set(z, list(sig.gene),
                               list(sig.gene[sig.direction == "UP"])),
        "proliferation": z.loc[[g for g in CYCLING_GENES if g in z.index]].mean(axis=0),
    }

    print("=== univariate AUC against metastasis ===")
    rows = []
    for name, s in scores.items():
        auc = roc_auc_score(y, s.values)
        auc = max(auc, 1 - auc)          # direction-free
        u = mannwhitneyu(s[y == 1], s[y == 0])
        rows.append({"score": name, "AUC": round(auc, 3),
                     "MWU p": f"{u.pvalue:.2g}"})
    print(pd.DataFrame(rows).to_string(index=False))

    print("\n=== score vs clinical covariates (Spearman) ===")
    rows = []
    for name, s in scores.items():
        for label, v in (("grade", grade), ("myometrial >=50%", mi),
                         ("metastasis", y)):
            r, p = spearmanr(s.values, v)
            rows.append({"score": name, "variable": label, "rho": round(r, 3),
                         "p": f"{p:.2g}"})
    print(pd.DataFrame(rows).to_string(index=False))

    print("\n=== each score alone, without the clinical covariates ===")
    cv = StratifiedKFold(CV_SPLITS, shuffle=True, random_state=SEED)
    for name, s in scores.items():
        a = cross_val_score(
            make_pipeline(StandardScaler(),
                          LogisticRegression(max_iter=LOGREG_MAX_ITER)),
            s.values.reshape(-1, 1), y, cv=cv, scoring="roc_auc").mean()
        print(f"  {name:14s} alone CV AUC = {a:.3f}")

    print("\n=== clinical covariates alone ===")
    print(f"  grade + myometrial invasion CV AUC = "
          f"{cv_auc(np.column_stack([grade, mi]), y):.3f}")


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------
def main() -> None:
    """Parse arguments and run the requested step(s)."""
    ap = argparse.ArgumentParser(
        description="Incremental predictive value of the IR programme.")
    ap.add_argument("--step", default="main",
                    choices=("main", "adjusted", "bootstrap", "diagnostics", "all"),
                    help="which analysis to run (default: main)")
    ap.add_argument("--programme", default="IR_program_v1.tsv",
                    help="programme table inside the tables directory")
    ap.add_argument("--modules", default="random_modules_v1.tsv",
                    help="random-module table inside the tables directory")
    ap.add_argument("--tag", default="",
                    help="suffix for the output file names, e.g. _v2")
    args = ap.parse_args()

    PATHS.ensure()
    steps = (("main", step_main), ("adjusted", step_adjusted),
             ("bootstrap", step_bootstrap), ("diagnostics", step_diagnostics))
    chosen = steps if args.step == "all" else [
        (n, f) for n, f in steps if n == args.step]

    for name, fn in chosen:
        print(f"\n===== {name} =====")
        fn(args.programme, args.modules, args.tag)


if __name__ == "__main__":
    main()
