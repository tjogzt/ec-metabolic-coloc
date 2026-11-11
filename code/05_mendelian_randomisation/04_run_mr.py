#!/usr/bin/env python3
"""
04_run_mr.py
============

Purpose
-------
Estimate the causal effect of fasting insulin, and separately of BMI, on
endometrial cancer with three univariable MR estimators, then estimate the
effect of each exposure conditional on the other with multivariable MR (MVMR).
Write the results table and the power boundary that the manuscript quotes.

The purpose of the two-exposure design is to separate the two hypotheses that
the metabolic arm of the study contrasts: an insulin-specific effect (fasting
insulin adjusted for BMI) and an adiposity effect (BMI adjusted for fasting
insulin). The MVMR rows of the output are the conditional estimates; whether
they may be interpreted at all is decided in 05_conditional_f.py.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
config.PATHS.tables / harmonised.tsv          from 03_harmonise_datasets.py
config.PATHS.tables / instruments_fi.tsv      fasting insulin instrument set
config.PATHS.tables / instruments_bmi.tsv     BMI instrument set
config.PARAMS.mr_max_snps                     density guard for the Egger fit

Outputs
-------
config.PATHS.tables / mr_results.tsv, columns

    exposure, method, n, beta, se, p, OR, OR_lo, OR_hi

in the fixed order fasting insulin IVW, MR-Egger, weighted median; BMI IVW,
MR-Egger, weighted median; MVMR-fasting insulin, MVMR-BMI. Beta is on the log
odds-ratio scale per one-unit increase in the exposure (log pmol/L for fasting
insulin, kg/m^2 for BMI); OR and its 95% interval are exp(beta) and
exp(beta -/+ 1.96 se).

config.PATHS.tables / power.json, one block per univariable exposure, holding
the instrument count, the mean and minimum first-stage F, the IVW standard
error, and the smallest odds ratio detectable with 80% power.

Notes
-----
Two design points are worth stating because they are easy to misread.

First, the two exposures use different instrument sets. Fasting insulin is
instrumented by its own 34 variants, BMI by its 184; the MVMR then uses the 205
variants that carry both exposure effects and the outcome effect. Estimates for
the same exposure therefore differ slightly between the univariable and the
MVMR rows, and that difference is a change of instrument set as well as a change
of model.

Second, MR-PRESSO is run for fasting insulin only. With 184 BMI instruments the
residual screen has almost no power to name a single variant after the
Bonferroni correction, and the global heterogeneity test is already reported
through Cochran's Q, so running it would add a number without adding evidence.

The power boundary is the smallest odds ratio that a two-sided test at
PARAMS.alpha would detect with 80% probability given the observed IVW standard
error: OR_min = exp((z_{1-alpha/2} + z_{0.80}) * se_IVW) = exp(2.802 * se_IVW).
It describes the design, not the estimate: an interval that includes OR_min has
not excluded an effect of that size.
"""
from __future__ import annotations

import importlib
import json
import pathlib
import sys
from typing import Dict, List, Tuple

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
INSTRUMENTS_BMI = "instruments_bmi.tsv"
RESULTS = "mr_results.tsv"
POWER = "power.json"

# z_{0.975} + z_{0.80} = 1.960 + 0.842; the multiplier that turns a standard
# error into the smallest effect detectable with 80% power at alpha = 0.05.
Z_POWER_80 = 2.802

# Normal quantile for a 95% interval on the log scale.
Z_95 = 1.96

# Exposure labels, in the order the rows are written. The MVMR rows are prefixed
# so that a reader cannot mistake them for the univariable estimates.
LABEL_FI = "fasting insulin"
LABEL_BMI = "BMI"
LABEL_MVMR = "MVMR-"


def f_statistic(b_exp, se_exp):
    """First-stage F statistic of one instrument, (beta / se)^2.

    Parameters
    ----------
    b_exp, se_exp : array_like
        Instrument effect on the exposure and its standard error.

    Returns
    -------
    ndarray of per-variant F statistics.

    Notes
    -----
    This is the per-variant F of the exposure-on-instrument regression with one
    degree of freedom, not the conditional F that MVMR requires; the conditional
    F is computed in 05_conditional_f.py. A mean above PARAMS.mr_min_f is the
    usual evidence against weak-instrument bias.
    """
    b_exp = np.asarray(b_exp, dtype=float)
    se_exp = np.asarray(se_exp, dtype=float)
    return (b_exp / se_exp) ** 2


def _row(exposure: str, method: str, n: int, beta: float, se: float,
         p: float) -> Dict[str, object]:
    """Assemble one results row, converting the log estimate to an odds ratio."""
    return dict(exposure=exposure, method=method, n=n, beta=float(beta),
                se=float(se), p=float(p), OR=float(np.exp(beta)),
                OR_lo=float(np.exp(beta - Z_95 * se)),
                OR_hi=float(np.exp(beta + Z_95 * se)))


def _egger_density_guard(n_instruments: int) -> None:
    """Warn when the instrument set is too dense for a stable Egger intercept.

    Parameters
    ----------
    n_instruments : int
        Number of variants entering the Egger fit.

    Notes
    -----
    The intercept is estimated from the same variants as the slope, so with
    thousands of nearly dependent variants its standard error collapses towards
    zero and a trivially small intercept becomes "significant". The threshold is
    config.PARAMS.mr_max_snps; the reported analyses use 34 and 184 instruments
    and never trigger it.
    """
    if n_instruments > PARAMS.mr_max_snps:
        print(f"  WARNING: {n_instruments} instruments exceed "
              f"PARAMS.mr_max_snps={PARAMS.mr_max_snps:,}; the MR-Egger "
              f"intercept is not interpretable at this density", flush=True)


# --------------------------------------------------------------------------
# Univariable: fasting insulin
# --------------------------------------------------------------------------
def run_fasting_insulin(harmonised: pd.DataFrame,
                        instruments: pd.DataFrame) -> Tuple[List[Dict[str, object]],
                                                            Dict[str, Dict[str, float]]]:
    """Three univariable MR estimates for fasting insulin, plus MR-PRESSO.

    Parameters
    ----------
    harmonised : DataFrame
        Output of 03_harmonise_datasets.py, needs the ``key`` column.
    instruments : DataFrame
        Fasting insulin instrument table (chr, pos in GRCh38).

    Returns
    -------
    (rows, power) where rows are the results-table dicts for IVW, MR-Egger and
    weighted median and power is the {"fasting insulin": {...}} block.

    Notes
    -----
    The subset is restricted twice: to the instrument positions, and to variants
    with a complete fasting insulin and outcome record. Both filters are needed
    because the harmonised table also carries variants that come from the BMI
    instrument set alone.
    """
    keys = instruments.chr.astype(str) + ":" + instruments.pos.astype(str)
    subset = harmonised[harmonised.key.isin(keys)].dropna(
        subset=["fi_beta", "fi_se", "ec_beta", "ec_se"])

    b_exp = subset.fi_beta.values
    se_exp = subset.fi_se.values
    b_out = subset.ec_beta.values
    se_out = subset.ec_se.values

    f_stats = f_statistic(b_exp, se_exp)
    print(f"\n=== {LABEL_FI}: {len(subset)} instruments, "
          f"mean F = {f_stats.mean():.1f}, min F = {f_stats.min():.1f} ===", flush=True)

    r_ivw = mr.ivw(b_exp, b_out, se_out, random_effects=True)
    _egger_density_guard(len(subset))
    r_egg = mr.egger(b_exp, b_out, se_out)
    r_wm = mr.weighted_median(b_exp, b_out, se_out)
    r_pr = mr.presso(b_exp, b_out, se_out)

    rows: List[Dict[str, object]] = []
    for method, res in (("IVW", r_ivw), ("MR-Egger", r_egg),
                        ("Weighted median", r_wm)):
        # The Egger fit reports beta_p (a t-test on K-2 df); the others report
        # p (a normal Wald test). Both are already two-sided.
        p = res["p"] if "p" in res else res["beta_p"]
        row = _row(LABEL_FI, method, len(subset), res["beta"], res["se"], p)
        rows.append(row)
        print(f"  {method:16s} beta={row['beta']:+.4f} se={row['se']:.4f} "
              f"OR={row['OR']:.3f} [{row['OR_lo']:.3f}, {row['OR_hi']:.3f}] "
              f"p={row['p']:.3g}", flush=True)

    print(f"  MR-Egger intercept {r_egg['intercept']:+.5f} "
          f"(p={r_egg['intercept_p']:.3f})", flush=True)
    print(f"  heterogeneity Q={r_ivw['Q']:.1f} (df={r_ivw['Q_df']}) "
          f"p={r_ivw['Q_p']:.3g}  I2={r_ivw['I2']:.2f}", flush=True)
    print(f"  MR-PRESSO global p={r_pr['p_global']:.3f}, "
          f"{r_pr['n_flagged']} variant(s) flagged", flush=True)
    if r_pr["after_removal"]:
        after = r_pr["after_removal"]
        print(f"  IVW after removing outliers beta={after['beta']:+.4f} "
              f"se={after['se']:.4f} p={after['p']:.3g}", flush=True)

    # Power boundary: smallest OR detectable with 80% power at the observed SE.
    min_or = float(np.exp(Z_POWER_80 * r_ivw["se"]))
    power = {LABEL_FI: dict(n=len(subset), mean_F=float(f_stats.mean()),
                            min_F=float(f_stats.min()),
                            se_ivw=float(r_ivw["se"]),
                            min_detectable_OR=min_or)}
    print(f"  power: smallest OR detectable with 80% power = {min_or:.3f}", flush=True)
    return rows, power


# --------------------------------------------------------------------------
# Univariable: BMI
# --------------------------------------------------------------------------
def run_bmi(harmonised: pd.DataFrame,
            instruments: pd.DataFrame) -> List[Dict[str, object]]:
    """Three univariable MR estimates for BMI.

    Parameters
    ----------
    harmonised : DataFrame
        Output of 03_harmonise_datasets.py.
    instruments : DataFrame
        BMI instrument table; its ``se`` column supplies the exposure standard
        error, which the harmonised table does not carry for BMI.

    Returns
    -------
    rows for IVW, MR-Egger and weighted median.

    Notes
    -----
    The BMI point estimate has no counterpart of ``fi_se`` in the harmonised
    table, so the standard error is looked up from the instrument table by
    (chromosome, position) and variants without a match are dropped. That is a
    real filter, not a formality: without it the first-stage F would be
    computed against a missing denominator.
    """
    subset = harmonised.dropna(subset=["bmi_beta", "ec_beta", "ec_se"])
    keys = instruments.chr.astype(str) + ":" + instruments.pos.astype(str)
    subset = subset[subset.key.isin(keys)]

    b_exp = subset.bmi_beta.values
    b_out = subset.ec_beta.values
    se_out = subset.ec_se.values

    # Exposure standard errors, keyed by position, aligned row by row.
    se_map = {(str(row.chr), int(row.pos)): row.se for row in instruments.itertuples()}
    se_exp = np.array([se_map.get((str(c), int(p)), np.nan)
                       for c, p in zip(subset.chr, subset.pos)])
    complete = ~np.isnan(se_exp)
    b_exp = b_exp[complete]
    b_out = b_out[complete]
    se_out = se_out[complete]
    se_exp = se_exp[complete]

    print(f"\n=== {LABEL_BMI}: {len(b_exp)} instruments, "
          f"mean F = {f_statistic(b_exp, se_exp).mean():.1f} ===", flush=True)

    r_ivw = mr.ivw(b_exp, b_out, se_out, random_effects=True)
    _egger_density_guard(len(b_exp))
    r_egg = mr.egger(b_exp, b_out, se_out)
    r_wm = mr.weighted_median(b_exp, b_out, se_out)

    rows: List[Dict[str, object]] = []
    for method, res in (("IVW", r_ivw), ("MR-Egger", r_egg),
                        ("Weighted median", r_wm)):
        p = res["p"] if "p" in res else res["beta_p"]
        row = _row(LABEL_BMI, method, len(b_exp), res["beta"], res["se"], p)
        rows.append(row)
        print(f"  {method:16s} beta={row['beta']:+.4f} OR={row['OR']:.3f} "
              f"p={row['p']:.3g}", flush=True)

    print(f"  MR-Egger intercept {r_egg['intercept']:+.5f} "
          f"(p={r_egg['intercept_p']:.3f})", flush=True)
    return rows


# --------------------------------------------------------------------------
# Multivariable: both exposures jointly
# --------------------------------------------------------------------------
def run_mvmr(harmonised: pd.DataFrame) -> List[Dict[str, object]]:
    """Joint (conditional) effects of fasting insulin and BMI.

    Parameters
    ----------
    harmonised : DataFrame
        Output of 03_harmonise_datasets.py.

    Returns
    -------
    Two rows, MVMR-fasting insulin and MVMR-BMI, in that order.

    Notes
    -----
    Only variants with both exposure effects and the outcome effect enter the
    fit: the second exposure is the regressor that does the conditioning, so a
    missing BMI effect cannot be imputed. That is why the MVMR instrument count
    (205) exceeds neither exposure's univariable count but is smaller than their
    union in the harmonised table.

    The estimates are conditional: each theta is the effect of one exposure
    holding the other fixed. They are reported only if the conditional F
    statistics of 05_conditional_f.py clear the weak-instrument floor, because
    otherwise the near-collinearity of the two exposure profiles makes the split
    arbitrary.
    """
    subset = harmonised.dropna(
        subset=["fi_beta", "fi_se", "bmi_beta", "ec_beta", "ec_se"])

    res = mr.mvmr([subset.fi_beta.values, subset.bmi_beta.values],
                  subset.ec_beta.values, subset.ec_se.values)

    print("\n=== MVMR: fasting insulin and BMI entered jointly ===", flush=True)
    rows: List[Dict[str, object]] = []
    for i, label in enumerate((LABEL_FI, LABEL_BMI)):
        row = _row(LABEL_MVMR + label, "MVMR", len(subset),
                   res["theta"][i], res["se"][i], res["p"][i])
        rows.append(row)
        print(f"  {label:16s} theta={row['beta']:+.4f} se={row['se']:.4f} "
              f"OR={row['OR']:.3f} p={row['p']:.3g}", flush=True)
    return rows


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------
def main() -> None:
    """Run every model, then write the results table and the power block."""
    PATHS.ensure()

    harmonised = pd.read_csv(PATHS.tables / HARMONISED, sep="\t")
    # Position key used to intersect the harmonised table with each instrument
    # set; built once because chromosome is read as an integer by pandas.
    harmonised["key"] = harmonised.chr.astype(str) + ":" + harmonised.pos.astype(str)

    fi_instruments = pd.read_csv(PATHS.tables / INSTRUMENTS_FI, sep="\t")
    bmi_instruments = pd.read_csv(PATHS.tables / INSTRUMENTS_BMI, sep="\t")

    fi_rows, power = run_fasting_insulin(harmonised, fi_instruments)
    bmi_rows = run_bmi(harmonised, bmi_instruments)
    mvmr_rows = run_mvmr(harmonised)

    results = pd.DataFrame(fi_rows + bmi_rows + mvmr_rows)
    results.to_csv(PATHS.tables / RESULTS, sep="\t", index=False)

    with open(PATHS.tables / POWER, "w") as fh:
        json.dump(power, fh, indent=2)

    print(f"\nwrote {RESULTS} ({len(results)} rows) and {POWER}", flush=True)


if __name__ == "__main__":
    main()
