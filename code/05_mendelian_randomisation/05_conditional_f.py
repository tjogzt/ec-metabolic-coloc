#!/usr/bin/env python3
"""
05_conditional_f.py
===================

Purpose
-------
Decide whether the multivariable MR estimates of 04_run_mr.py may be
interpreted, and test how much any single variant drives the univariable
fasting insulin estimate.

Two checks live here because they answer the same question from opposite sides:
the conditional F statistic asks whether each exposure still carries
independent instrument strength once the other exposure is controlled, and the
leave-one-out analysis asks whether the causal estimate survives removal of any
one variant.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
config.PATHS.tables / harmonised.tsv        from 03_harmonise_datasets.py
config.PATHS.tables / instruments_fi.tsv    fasting insulin instrument set
config.PARAMS.mr_min_f                      weak-instrument floor (10)
config.PARAMS.alpha                         two-sided test level for the
                                            leave-one-out significance count

Outputs
-------
config.PATHS.tables / cond_f.tsv, columns exposure, F_cond, R2, K:
two rows, "fasting insulin | BMI" (fasting insulin conditioned on BMI) and
"BMI | fasting insulin".

config.PATHS.tables / loo_fi.tsv, columns dropped_rsid, dropped_chr,
dropped_pos, beta, se, p, OR: one row per fasting insulin variant, holding the
IVW estimate obtained when that variant is omitted.

Notes
-----
The conditional F statistic is the decision rule for the MVMR rows. An exposure
whose conditional F falls below PARAMS.mr_min_f (10) is not identified by the
instrument set once the other exposure is held fixed: the two exposure profiles
are then close to collinear across variants, the split of the total effect
between them is arbitrary, and any conditional estimate must be described as
uninterpretable rather than merely imprecise. In this analysis both orderings
fall below the floor, which is why the MVMR rows are marked invalid in the
results table and greyed out in the supplementary forest plot.

The leave-one-out analysis is run on the fasting insulin instrument set after
restricting to variants that also carry a BMI effect, because the harmonised
table is shared with the MVMR model. That subset is smaller than the
univariable instrument count (22 against 34), so its full-sample estimate is not
identical to the IVW row of mr_results.tsv; the figure script therefore fixes
its reference line at the estimate of this subset, not at the table row.
"""
from __future__ import annotations

import importlib
import pathlib
import sys
from typing import Dict, List

import numpy as np
import pandas as pd

# Repository root, so that ``config`` resolves when the file is run directly.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

from config.config import PATHS, PARAMS, SEED  # noqa: E402  (after sys.path edit)

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# The estimator library has a numeric file name, which cannot appear in an
# import statement; importlib resolves it from the directory added above.
mr = importlib.import_module("01_mr_estimators")

HARMONISED = "harmonised.tsv"
INSTRUMENTS_FI = "instruments_fi.tsv"
COND_F = "cond_f.tsv"
LOO_FI = "loo_fi.tsv"

# Columns that must be present for a variant to enter either the conditional F
# regression or the multivariable model: both exposures, with their standard
# errors, and the outcome with its standard error.
COMPLETE_COLUMNS = ["fi_beta", "fi_se", "bmi_beta", "bmi_se", "ec_beta", "ec_se"]


def conditional_f(b_target, b_other, se_target, se_other) -> Dict[str, float]:
    """Conditional F statistic of one exposure given the other.

    Parameters
    ----------
    b_target : array_like
        Per-variant effect on the exposure whose strength is being tested.
    b_other : array_like
        Per-variant effect on the exposure being conditioned on.
    se_target : array_like
        Standard error of ``b_target``; supplies the regression weights.
    se_other : array_like
        Standard error of ``b_other``. Unused in the computation, kept in the
        signature because the four arrays are always passed as a set.

    Returns
    -------
    dict with F_cond, R2 and K.

    Formula
    -------
    Weighted least squares of the target exposure effect on the conditional one,

        b_target_k = a + b * b_other_k + e_k,   var(e_k) = se_target_k^2,

    with weights w_k = 1 / se_target_k^2. The weighted coefficient of
    determination is

        R^2 = 1 - sum(w * residual^2) / sum(w * (b_target - weighted mean)^2),

    and the conditional F, which asks whether the instruments explain the
    target exposure *after* the other exposure has taken its share, is

        F = (R^2 / (1 - R^2)) * (K - 2).

    The K - 2 is the residual degrees of freedom of a two-parameter weighted
    regression; it is what makes this a conditional statistic rather than the
    marginal F reported in 04_run_mr.py.

    Notes
    -----
    Edge cases are returned as finite values rather than NaN so that the
    decision rule downstream never silently passes: R^2 = 1 (a perfect fit)
    returns infinity, and R^2 <= 0 or an undefined R^2 returns 0. Both outcomes
    are decisive under a < 10 rule, in opposite directions.
    """
    b_target = np.asarray(b_target, dtype=float)
    b_other = np.asarray(b_other, dtype=float)
    se_target = np.asarray(se_target, dtype=float)

    # Weights, then rows scaled by sqrt(w) so that ordinary least squares on the
    # scaled problem is the weighted fit; this avoids a diagonal weight matrix.
    w = 1.0 / se_target ** 2
    sw = np.sqrt(w)
    design = np.column_stack([np.ones(len(b_other)), b_other])
    X = design * sw[:, None]
    y = b_target * sw

    coef, *_ = np.linalg.lstsq(X, y, rcond=None)

    # Residual and total weighted sums of squares, on the original scale.
    pred = design @ coef
    ss_res = float(np.sum(w * (b_target - pred) ** 2))
    y_bar = np.average(b_target, weights=w)
    ss_tot = float(np.sum(w * (b_target - y_bar) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    k = len(b_target)
    if 0 < r2 < 1:
        f_value = (r2 / (1.0 - r2)) * (k - 2)
    else:
        f_value = np.inf if r2 >= 1 else 0.0
    return dict(F_cond=float(f_value), R2=float(r2), K=k)


def leave_one_out(harmonised: pd.DataFrame,
                  instruments: pd.DataFrame) -> pd.DataFrame:
    """Re-estimate the IVW effect with each fasting insulin variant omitted.

    Parameters
    ----------
    harmonised : DataFrame
        Already restricted to complete cases and keyed by ``key``.
    instruments : DataFrame
        Fasting insulin instrument table, used only for its position key.

    Returns
    -------
    DataFrame with dropped_rsid, dropped_chr, dropped_pos, beta, se, p, OR.

    Notes
    -----
    One IVW refit per variant, so the cost is O(K^2) in the vectorised sums; K is
    in the tens here. The alternative, an analytic influence formula, would
    assume the random-effects standard error is differentiable in the dropped
    variant, which it is not (the Q > df condition can switch).
    """
    keys = instruments.chr.astype(str) + ":" + instruments.pos.astype(str)
    subset = harmonised[harmonised.key.isin(keys)].reset_index(drop=True)
    print(f"instruments entering the leave-one-out analysis: {len(subset)}", flush=True)

    full = mr.ivw(subset.fi_beta.values, subset.ec_beta.values, subset.ec_se.values)

    rows: List[Dict[str, object]] = []
    for i in range(len(subset)):
        keep = np.ones(len(subset), dtype=bool)
        keep[i] = False
        res = mr.ivw(subset.fi_beta.values[keep], subset.ec_beta.values[keep],
                     subset.ec_se.values[keep])
        rows.append(dict(dropped_rsid=subset.rsid.iloc[i],
                         dropped_chr=subset.chr.iloc[i],
                         dropped_pos=subset.pos.iloc[i],
                         beta=res["beta"], se=res["se"], p=res["p"],
                         OR=float(np.exp(res["beta"]))))

    loo = pd.DataFrame(rows)
    print(f"  full estimate beta = {full['beta']:+.4f} "
          f"(OR {np.exp(full['beta']):.3f})", flush=True)
    print(f"  leave-one-out beta range {loo.beta.min():+.4f} to "
          f"{loo.beta.max():+.4f}; all p < {PARAMS.alpha:g}: "
          f"{bool((loo.p < PARAMS.alpha).all())}", flush=True)

    # The variant whose omission weakens the estimate most (largest p).
    worst = loo.loc[loo.p.idxmax()]
    print(f"  most influential variant {worst.dropped_rsid} "
          f"(chr{worst.dropped_chr}:{worst.dropped_pos}): "
          f"beta={worst.beta:+.4f} p={worst.p:.3g}", flush=True)
    return loo


def main() -> None:
    """Run the conditional F check, then the leave-one-out analysis."""
    PATHS.ensure()

    harmonised = pd.read_csv(PATHS.tables / HARMONISED, sep="\t")
    harmonised["key"] = (harmonised.chr.astype(str) + ":"
                         + harmonised.pos.astype(str))

    complete = harmonised.dropna(subset=COMPLETE_COLUMNS)
    print(f"variants with both exposures: {len(complete)}", flush=True)

    # ---- 1. conditional F, in both directions ---------------------------
    rows: List[Dict[str, object]] = []
    for label, target, other in (
        ("fasting insulin | BMI", ("fi_beta", "fi_se"), ("bmi_beta", "bmi_se")),
        ("BMI | fasting insulin", ("bmi_beta", "bmi_se"), ("fi_beta", "fi_se")),
    ):
        res = conditional_f(complete[target[0]].values, complete[other[0]].values,
                            complete[target[1]].values, complete[other[1]].values)
        rows.append(dict(exposure=label, F_cond=res["F_cond"], R2=res["R2"],
                         K=res["K"]))
        verdict = (f"UNRELIABLE (< {PARAMS.mr_min_f:g})"
                   if res["F_cond"] < PARAMS.mr_min_f
                   else f"acceptable (>= {PARAMS.mr_min_f:g})")
        print(f"  {label:26s} conditional F = {res['F_cond']:7.2f}  "
              f"R2 = {res['R2']:.3f}  K = {res['K']}  {verdict}", flush=True)

    pd.DataFrame(rows)[["exposure", "F_cond", "R2", "K"]].to_csv(
        PATHS.tables / COND_F, sep="\t", index=False)

    # ---- 2. leave-one-out ----------------------------------------------
    print("\n=== Leave-one-out (fasting insulin, IVW) ===", flush=True)
    fi_instruments = pd.read_csv(PATHS.tables / INSTRUMENTS_FI, sep="\t")
    loo = leave_one_out(complete, fi_instruments)
    loo.to_csv(PATHS.tables / LOO_FI, sep="\t", index=False)


if __name__ == "__main__":
    main()
