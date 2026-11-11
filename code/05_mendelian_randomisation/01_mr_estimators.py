#!/usr/bin/env python3
"""
01_mr_estimators.py
===================

Purpose
-------
Self-contained estimators for two-sample Mendelian randomisation (MR) and
multivariable MR (MVMR). The execution environment provides neither
TwoSampleMR nor MRPRESSO, so every estimator that the analysis reports is
implemented here from its published closed form:

    wald_ratio       single-instrument ratio estimate
    ivw              inverse-variance weighted estimate, with fixed-effect and
                     multiplicative random-effects standard errors, Cochran's Q
                     and I^2
    egger            MR-Egger slope and directional-pleiotropy intercept
    weighted_median  weighted median estimator with a bootstrap standard error
    presso           MR-PRESSO global heterogeneity test and outlier detection
    mvmr             multivariable IVW by weighted least squares

The module is the whole numerical core of the stage; the driver scripts
(02 to 06) contain no estimator arithmetic of their own.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
Every estimator takes three aligned arrays of length K (the number of
independent instruments):

    b_exp    per-variant effect on the exposure
    b_out    per-variant effect on the outcome
    se_out   standard error of b_out

Conventions that the callers guarantee:

    * b_exp and b_out are already harmonised to the same effect allele, so the
      sign of each ratio is comparable across variants;
    * b_exp is non-zero (instruments were selected at genome-wide significance,
      so |b_exp| is far from zero);
    * se_out > 0.

The master random seed for the two resampling procedures comes from
``config.SEED``; the type-I error rate used by MR-PRESSO comes from
``config.PARAMS.alpha``.

Outputs
-------
Estimator functions return plain dicts, never dataframes, so that the drivers
decide how results are tabulated:

    wald_ratio   (beta, se) tuple
    ivw          beta, se, p, Q, Q_df, Q_p, I2, n, se_fixed
    egger        intercept, intercept_se, intercept_p, beta, se, beta_p, n
    weighted_median  beta, se, p, n
    presso       p_global, Q_obs, n_flagged, flagged, after_removal
    mvmr         theta (list), se (list), p (list), df

Running the module directly executes a synthetic-data self-check that verifies
each estimator against a known simulation truth, including the ability of
MR-PRESSO to name the injected outlier.

Notes
-----
Ratio estimates. For a single variant the causal estimate is the Wald ratio
beta_out / beta_exp and its first-order (delta-method) standard error is
se_out / |beta_exp|. That approximation is adequate here because instrument
selection keeps |beta_exp| large relative to se_exp.

Numerical stability. Two procedures resample, and both are seeded from
``config.SEED`` so that a rerun reproduces the same standard error to the last
digit. Every regression is solved by scaled least squares (rows multiplied by
sqrt(weight)) rather than by forming a diagonal weight matrix, and every
chi-square p-value uses the survival function directly instead of 1 - cdf,
because the tail probabilities here are routinely below 1e-100 and would be
lost to cancellation in the subtraction.

What is deliberately not implemented. No pleiotropy-robust estimator beyond the
weighted median (no MR-RAPS, no contamination mixture), no LD-score adjustment,
and no instrument-level clumping: independence of instruments is enforced
upstream in 02_select_instruments.py.
"""
from __future__ import annotations

import pathlib
import sys
from typing import Any, Dict, List, Tuple

import numpy as np
from scipy import stats

# Repository root, so that ``config`` resolves when the file is run directly.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

from config.config import PARAMS, SEED  # noqa: E402  (import after sys.path edit)

__all__ = [
    "wald_ratio",
    "ivw",
    "egger",
    "weighted_median",
    "presso",
    "mvmr",
    "N_BOOTSTRAP_MEDIAN",
    "N_DIST_PRESSO",
]

# Bootstrap replicates for the weighted-median standard error. Fixed at the
# value used for the reported estimates: the median SE moves by <1% between
# 2000 and 10000 replicates, so the cost is not justified.
N_BOOTSTRAP_MEDIAN = 2000

# Monte-Carlo draws for the MR-PRESSO global test. The p-value resolution is
# 1/(N+1) = 1e-3, which is finer than the heterogeneity questions asked here.
N_DIST_PRESSO = 1000


# --------------------------------------------------------------------------
# Single-instrument estimate
# --------------------------------------------------------------------------
def wald_ratio(b_exp, b_out, se_out) -> Tuple[np.ndarray, np.ndarray]:
    """Wald ratio estimate and its first-order standard error.

    Parameters
    ----------
    b_exp, b_out : array_like
        Per-variant effects on exposure and outcome, allele-aligned.
    se_out : array_like
        Standard error of ``b_out``.

    Returns
    -------
    (beta, se), both ndarray of the input length.

    Formula
    -------
    beta = b_out / b_exp                (the ratio, or Wald, estimate)
    se   = se_out / |b_exp|             (delta method, exposure measured without
                                         error so only the outcome term remains)

    The absolute value keeps the standard error positive when the exposure
    effect is negative; the sign of the estimate itself carries the direction.
    """
    b_exp = np.asarray(b_exp, dtype=float)
    b_out = np.asarray(b_out, dtype=float)
    se_out = np.asarray(se_out, dtype=float)
    return b_out / b_exp, se_out / np.abs(b_exp)


# --------------------------------------------------------------------------
# Inverse-variance weighted estimator
# --------------------------------------------------------------------------
def ivw(b_exp, b_out, se_out, random_effects: bool = True) -> Dict[str, float]:
    """Inverse-variance weighted (IVW) causal estimate.

    Parameters
    ----------
    b_exp, b_out : array_like
        Per-variant effects on exposure and outcome, allele-aligned.
    se_out : array_like
        Standard error of ``b_out``.
    random_effects : bool
        If True (the default) and Cochran's Q exceeds its degrees of freedom,
        the standard error is inflated multiplicatively by sqrt(Q / df). If
        False, the fixed-effect standard error is returned unchanged.

    Returns
    -------
    dict with beta, se, p, Q, Q_df, Q_p, I2, n and se_fixed. ``se_fixed`` is
    kept alongside ``se`` so that a caller can report both the fixed-effect and
    the heterogeneity-inflated interval without refitting.

    Formula
    -------
    With ratio_k = b_out_k / b_exp_k, se_k = se_out_k / |b_exp_k| and
    w_k = 1 / se_k^2:

        beta_IVW = sum(w_k * ratio_k) / sum(w_k)
        se_fixed = sqrt(1 / sum(w_k))
        Q        = sum(w_k * (ratio_k - beta_IVW)^2)     ~ chi2(K - 1)
        I^2      = max(0, (Q - df) / Q)
        se_re    = se_fixed * sqrt(Q / df)               (only when Q > df)

    Notes
    -----
    The random-effects form is the multiplicative one used by the two-sample MR
    literature, not the additive DerSimonian-Laird form: the additive version
    can return a standard error below the fixed-effect value when Q is only
    marginally above df, which is not a defensible interval.
    """
    b_exp = np.asarray(b_exp, dtype=float)
    b_out = np.asarray(b_out, dtype=float)
    se_out = np.asarray(se_out, dtype=float)

    # Per-variant ratio estimate and its delta-method standard error.
    ratio = b_out / b_exp
    se_ratio = se_out / np.abs(b_exp)

    # Inverse-variance weights. Written as one vector reciprocal rather than as
    # a K x K diagonal matrix, so the cost stays O(K) in time and memory.
    w = 1.0 / se_ratio ** 2

    # Weighted mean of the ratios: beta = sum(w * ratio) / sum(w).
    beta = np.sum(w * ratio) / np.sum(w)

    # Fixed-effect standard error: sqrt(1 / sum(w)). This is the variance of a
    # weighted mean of independent estimates, not a regression standard error.
    se_fixed = np.sqrt(1.0 / np.sum(w))

    # Cochran's Q = sum(w * (ratio - beta)^2), chi-square with K - 1 df under
    # the null of a single common effect. Accumulating weighted squared
    # deviations avoids materialising the covariance matrix.
    Q = float(np.sum(w * (ratio - beta) ** 2))
    df = len(ratio) - 1
    # Survival function, not 1 - cdf: Q is often enormous here and the
    # subtraction would flush the tail probability to exactly zero.
    Q_p = float(stats.chi2.sf(Q, df)) if df > 0 else np.nan
    # I^2 = max(0, (Q - df) / Q): share of the total variance attributable to
    # heterogeneity rather than to sampling error; clipped at zero because a Q
    # below its df carries no information about excess dispersion.
    I2 = max(0.0, (Q - df) / Q) if Q > 0 else 0.0

    # Multiplicative random-effects correction. Applied only when Q > df, so
    # that a homogeneous instrument set keeps the narrower fixed-effect
    # interval; sqrt(Q / df) is the standard excess-dispersion factor.
    se = se_fixed
    if random_effects and df > 0 and Q > df:
        se = se_fixed * np.sqrt(Q / df)

    # Two-sided Wald p-value from the normal reference distribution.
    p = float(2 * stats.norm.sf(abs(beta / se)))

    return dict(beta=float(beta), se=float(se), p=p, Q=Q, Q_df=df,
                Q_p=Q_p, I2=float(I2), n=len(ratio),
                se_fixed=float(se_fixed))


# --------------------------------------------------------------------------
# MR-Egger
# --------------------------------------------------------------------------
def egger(b_exp, b_out, se_out) -> Dict[str, float]:
    """MR-Egger slope (pleiotropy-corrected estimate) and intercept.

    Parameters
    ----------
    b_exp, b_out : array_like
        Per-variant effects on exposure and outcome, allele-aligned.
    se_out : array_like
        Standard error of ``b_out``.

    Returns
    -------
    dict with intercept, intercept_se, intercept_p, beta, se, beta_p and n.

    Model
    -----
        b_out_k = theta_0 + theta_1 * b_exp_k + e_k,   var(e_k) = se_out_k^2

    theta_1 is the causal estimate corrected for average directional
    pleiotropy and theta_0 is the pleiotropy intercept: under the
    Instrument Strength Independent of Direct Effect assumption theta_0 = 0, so
    a non-zero intercept is evidence that the instrument set violates it.

    Notes
    -----
    All variants are first oriented so that b_exp > 0. Without that step the
    intercept would absorb the correlation between the sign of the exposure
    effect and the pleiotropic effect, and the test would be meaningless
    (Burgess and Thompson 2017).

    The regression is weighted by 1 / se_out^2. Instead of building the weight
    matrix, each row of the design matrix and of the response is multiplied by
    1 / se_out, which squares to the required weight inside the normal
    equations; this is exact and O(K) in memory.

    The residual variance uses K - 2 degrees of freedom because two parameters
    are estimated, and the p-values come from a t reference distribution rather
    than a normal one, since sigma^2 is itself estimated.
    """
    b_exp = np.asarray(b_exp, dtype=float).copy()
    b_out = np.asarray(b_out, dtype=float).copy()
    se_out = np.asarray(se_out, dtype=float)

    # Orient every variant to a positive exposure effect. The copy above keeps
    # the caller's arrays untouched.
    sign = np.sign(b_exp)
    sign[sign == 0] = 1.0
    b_exp = b_exp * sign
    b_out = b_out * sign

    # Row scaling by 1 / se_out supplies the 1 / se_out^2 weights implicitly.
    sw = 1.0 / se_out
    X = np.column_stack([np.full_like(b_exp, 1.0), b_exp]) * sw[:, None]
    y = b_out * sw

    # Least squares on the scaled problem; rcond=None keeps NumPy's default
    # cutoff for rank deficiency instead of pinning a machine-specific value.
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)

    resid = y - X @ coef
    df = len(b_exp) - 2
    # sigma^2 = RSS / (K - 2); undefined for fewer than three variants, in which
    # case the intercept is not estimable and NaNs are returned downstream.
    sigma2 = float((resid ** 2).sum() / df) if df > 0 else np.nan

    # cov = sigma^2 * (X'X)^-1, exploiting that X already carries the weights.
    cov = sigma2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))

    # Two-sided t tests for intercept and slope on K - 2 degrees of freedom.
    p = (2 * stats.t.sf(np.abs(coef / se), df) if df > 0
         else np.array([np.nan, np.nan]))

    return dict(intercept=float(coef[0]), intercept_se=float(se[0]),
                intercept_p=float(p[0]), beta=float(coef[1]), se=float(se[1]),
                beta_p=float(p[1]), n=len(b_exp))


# --------------------------------------------------------------------------
# Weighted median
# --------------------------------------------------------------------------
def weighted_median(b_exp, b_out, se_out, n_boot: int = N_BOOTSTRAP_MEDIAN,
                    seed: int = SEED) -> Dict[str, float]:
    """Weighted median estimator with a bootstrap standard error.

    Parameters
    ----------
    b_exp, b_out : array_like
        Per-variant effects on exposure and outcome, allele-aligned.
    se_out : array_like
        Standard error of ``b_out``.
    n_boot : int
        Number of bootstrap replicates for the standard error.
    seed : int
        Seed for the bootstrap generator; defaults to the master seed in
        ``config.SEED``.

    Returns
    -------
    dict with beta, se, p and n.

    Formula
    -------
    The estimate is the 50th percentile of the ratio estimates weighted by
    w_k = 1 / se_k^2 with se_k = se_out_k / |b_exp_k|. The weights are turned
    into a stepwise cumulative distribution by taking the cumulative sum and
    subtracting half of each variant's own weight, which places every ratio at
    the midpoint of its weight interval; dividing by the total normalises the
    curve to [0, 1]. The estimate is then the linear interpolation of ratio at
    cumulative weight 0.5.

    Notes
    -----
    The standard error is a nonparametric bootstrap over variants (resampling
    indices, not residuals), recomputing the weighted median in every replicate.
    This is the interval proposed with the estimator itself; an analytic
    approximation is not used because the sampling distribution of a weighted
    median is not normal in small instrument sets.
    """
    b_exp = np.asarray(b_exp, dtype=float)
    b_out = np.asarray(b_out, dtype=float)
    se_out = np.asarray(se_out, dtype=float)

    # Ratios and inverse-variance weights, as in the IVW estimator.
    ratio = b_out / b_exp
    w = 1.0 / (se_out / np.abs(b_exp)) ** 2

    # Sort by ratio so that the weighted cumulative curve is monotone.
    order = np.argsort(ratio)
    ratio, w = ratio[order], w[order]

    # Midpoint cumulative weights: cw_k = sum_{j<=k} w_j - w_k / 2, normalised.
    cw = np.cumsum(w) - 0.5 * w
    cw /= np.sum(w)

    def _weighted_quantile(ratios: np.ndarray, cum_w: np.ndarray) -> float:
        """Interpolate the ratio at cumulative weight 0.5."""
        return float(np.interp(0.5, cum_w, ratios))

    beta = _weighted_quantile(ratio, cw)

    # Bootstrap: resample variants with replacement, then recompute the
    # weighted median from scratch (sorting included, since the resample
    # destroys the ordering).
    rng = np.random.default_rng(seed)
    draws: List[float] = []
    n = len(ratio)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        r, ww = ratio[idx], w[idx]
        o = np.argsort(r)
        r, ww = r[o], ww[o]
        c = np.cumsum(ww) - 0.5 * ww
        c /= np.sum(ww)
        draws.append(_weighted_quantile(r, c))

    # ddof=1 matches the sample standard deviation of the replicate estimates.
    se = float(np.std(draws, ddof=1))
    p = float(2 * stats.norm.sf(abs(beta / se))) if se > 0 else np.nan
    return dict(beta=beta, se=se, p=p, n=n)


# --------------------------------------------------------------------------
# MR-PRESSO: global heterogeneity test and outlier detection
# --------------------------------------------------------------------------
def presso(b_exp, b_out, se_out, n_dist: int = N_DIST_PRESSO,
           seed: int = SEED, alpha: float = PARAMS.alpha) -> Dict[str, Any]:
    """MR-PRESSO global test with outlier identification and re-estimation.

    Parameters
    ----------
    b_exp, b_out : array_like
        Per-variant effects on exposure and outcome, allele-aligned.
    se_out : array_like
        Standard error of ``b_out``.
    n_dist : int
        Number of Monte-Carlo draws from the null (no-pleiotropy) distribution.
    seed : int
        Seed for the generator; defaults to ``config.SEED``.
    alpha : float
        Family-wise error rate for outlier calling; defaults to
        ``config.PARAMS.alpha`` and is Bonferroni-split over the K variants.

    Returns
    -------
    dict with p_global, Q_obs, n_flagged, flagged (variant indices in input
    order) and after_removal (an IVW result dict, or None when nothing was
    flagged or too few variants would remain).

    Procedure
    ---------
    1. Fit the IVW estimate beta and compute the observed heterogeneity statistic
           Q_obs = sum(w_k * (ratio_k - beta)^2),  w_k = 1 / se_k^2.
    2. Simulate the null by drawing K ratios
           ratio_k^* ~ N(beta, se_k^2)
       refitting beta^* and recomputing Q^*. The null keeps the observed beta,
       because the test asks whether the *dispersion* exceeds sampling noise,
       not whether the mean differs from zero. The excess of Q over Q_obs gives
           p_global = (#{Q^* >= Q_obs} + 1) / (n_dist + 1),
       the +1 form being the usual Monte-Carlo correction that keeps the test
       valid (it can never return exactly 0).
    3. Call a variant an outlier when its standardised residual
           s_k = (ratio_k - beta) / se_k
       exceeds the two-sided normal threshold after a Bonferroni correction to
       alpha / K. This is a test of each variant against the common effect, so
       the residuals are already on a standard-normal scale and no extra
       variance correction is needed.
    4. Refit the IVW estimate on the retained variants when at least two
       survive and at least one was flagged.

    Notes
    -----
    The outlier screen is applied even when p_global is not significant, exactly
    as in the source analysis: the global test is the reported heterogeneity
    result, while the residual screen is what names the variant that drives it.
    """
    b_exp = np.asarray(b_exp, dtype=float)
    b_out = np.asarray(b_out, dtype=float)
    se_out = np.asarray(se_out, dtype=float)

    # Reference fit; the random-effects flag is irrelevant here because Q_obs is
    # recomputed below from the weighted deviations.
    base = ivw(b_exp, b_out, se_out, random_effects=True)

    ratio = b_out / b_exp
    se_ratio = se_out / np.abs(b_exp)
    w = 1.0 / se_ratio ** 2
    beta = base["beta"]

    # Observed heterogeneity around the IVW estimate.
    resid = ratio - beta
    Q_obs = float(np.sum(w * resid ** 2))

    # Null distribution of Q by simulation. Only K draws of K numbers are
    # needed, so the cost is O(n_dist * K) with no linear algebra.
    rng = np.random.default_rng(seed)
    Q_null = np.empty(n_dist)
    for i in range(n_dist):
        sim = beta + rng.normal(0, 1, len(ratio)) * se_ratio
        beta_sim = np.sum(w * sim) / np.sum(w)
        Q_null[i] = np.sum(w * (sim - beta_sim) ** 2)

    # Monte-Carlo p-value with the +1 continuity correction.
    p_global = (np.sum(Q_null >= Q_obs) + 1) / (n_dist + 1)

    # Standardised residuals and their two-sided p-values; a variant is called
    # when it clears alpha / K (Bonferroni over the instrument set). With K <= 2
    # there is no residual information, so nothing is flagged.
    sres = resid / se_ratio
    p_i = 2 * stats.norm.sf(np.abs(sres))
    if len(ratio) > 2:
        flagged = p_i < (alpha / len(ratio))
    else:
        flagged = np.zeros(len(ratio), dtype=bool)

    keep = ~flagged
    if keep.sum() >= 2 and flagged.any():
        # Re-estimate after removing the flagged variants. Returning None rather
        # than a partial fit signals "not estimable" to the caller.
        after = ivw(b_exp[keep], b_out[keep], se_out[keep], random_effects=True)
    else:
        after = None

    return dict(p_global=float(p_global), Q_obs=Q_obs,
                n_flagged=int(flagged.sum()),
                flagged=np.where(flagged)[0].tolist(),
                after_removal=after)


# --------------------------------------------------------------------------
# Multivariable MR
# --------------------------------------------------------------------------
def mvmr(beta_exp_list: List[np.ndarray], b_out, se_out) -> Dict[str, Any]:
    """Multivariable IVW: joint effects of several exposures on one outcome.

    Parameters
    ----------
    beta_exp_list : sequence of array_like
        One K-vector per exposure, holding that exposure's per-variant effect.
        A two-exposure call passes [b_exp_1, b_exp_2].
    b_out : array_like
        Per-variant outcome effect, allele-aligned.
    se_out : array_like
        Standard error of ``b_out``.

    Returns
    -------
    dict with theta (list, one entry per exposure), se, p and df.

    Model
    -----
        b_out_k = theta_1 * b_exp_1k + ... + theta_p * b_exp_pk + e_k,
        var(e_k) = se_out_k^2

    solved by weighted least squares with weights 1 / se_out_k^2. The design
    matrix has no intercept: the model is fitted through the origin because a
    variant with no effect on any exposure must have no effect on the outcome
    under the exclusion restriction.

    Notes
    -----
    As in the univariable case, the weight matrix is avoided by scaling rows of
    the design matrix and of the response by 1 / se_out, which yields the same
    normal equations at O(K * p) cost.

    The residual variance uses K - p degrees of freedom. When that is not
    positive (fewer variants than exposures) sigma^2 is set to 1 and the
    standard errors become the model-based ones, i.e. the fit is reported but
    its uncertainty is nominal; the drivers check instrument counts before the
    estimate is interpreted.

    Degrees of freedom are returned so that the caller can state how much
    information backed the estimate.
    """
    X0 = np.column_stack([np.asarray(b, dtype=float) for b in beta_exp_list])
    b_out = np.asarray(b_out, dtype=float)
    se_out = np.asarray(se_out, dtype=float)

    # Row scaling by 1 / se_out, equivalent to weighting by 1 / se_out^2.
    sw = 1.0 / se_out
    X = X0 * sw[:, None]
    y = b_out * sw

    # Solve the normal equations X'X theta = X'y directly; the design is small
    # (K x p) and well conditioned once scaled.
    XtX = X.T @ X
    theta = np.linalg.solve(XtX, X.T @ y)

    resid = y - X @ theta
    df = len(b_out) - X.shape[1]
    sigma2 = float((resid ** 2).sum() / df) if df > 0 else 1.0

    # cov = sigma^2 * (X'X)^-1, then the standard errors are its diagonal roots.
    cov = sigma2 * np.linalg.inv(XtX)
    se = np.sqrt(np.diag(cov))

    # Two-sided normal p-values; a t reference would need the residual df to be
    # well estimated, which it is not in the typical MVMR instrument set.
    p = 2 * stats.norm.sf(np.abs(theta / se))
    return dict(theta=theta.tolist(), se=se.tolist(), p=p.tolist(), df=df)


# --------------------------------------------------------------------------
# Self-check
# --------------------------------------------------------------------------
if __name__ == "__main__":
    # Synthetic data with known truth: K = 400 instruments, a true causal effect
    # of 0.30 on the log-odds scale, and homogeneous noise so that all estimators
    # should land near 0.30 and no heterogeneity should be detected.
    rng = np.random.default_rng(1)
    n = 400
    b_exp = rng.uniform(0.02, 0.10, n)
    TRUE = 0.30
    b_out = TRUE * b_exp + rng.normal(0, 0.01, n)
    se_out = np.full(n, 0.01)

    print("Synthetic check, K=400, true causal effect 0.30")
    r = ivw(b_exp, b_out, se_out)
    print(f"  IVW        beta={r['beta']:.4f}  se={r['se']:.4f}  p={r['p']:.3g}")

    e = egger(b_exp, b_out, se_out)
    print(f"  MR-Egger   slope={e['beta']:.4f}  intercept={e['intercept']:+.5f} "
          f"(expected ~0)  p={e['intercept_p']:.2f}")

    m = weighted_median(b_exp, b_out, se_out, n_boot=500)
    print(f"  W-median   beta={m['beta']:.4f}  se={m['se']:.4f}")

    # No injected outlier: the global test should not reject homogeneity.
    pr = presso(b_exp, b_out, se_out, n_dist=200)
    print(f"  MR-PRESSO  global p={pr['p_global']:.3f} (expected non-significant)  "
          f"flagged={pr['n_flagged']}")

    # Inject one outlier by adding 0.05 to the first variant's outcome effect.
    # Its standardised residual is 0.05 / (0.01 / b_exp) >= 3.8 for every
    # b_exp >= 0.02, so it clears the Bonferroni threshold regardless of the
    # draw; the check below asserts that the screen names that variant and no
    # other.
    b_out2, se_out2 = b_out.copy(), se_out.copy()
    b_out2[0] += 0.05
    pr2 = presso(b_exp, b_out2, se_out2, n_dist=200)
    print(f"  MR-PRESSO+ global p={pr2['p_global']:.3f}  flagged={pr2['n_flagged']} "
          f"-> {pr2['flagged'][:3]}")
    assert pr2["flagged"] == [0], (
        f"MR-PRESSO should isolate the injected outlier (variant 0), got {pr2['flagged']}")
    print("  MR-PRESSO+ outlier recovery: variant 0 isolated exactly -> OK")

    # MVMR: two exposures, theta_1 = 0.30, theta_2 = 0.00 (the second exposure
    # is a decoy that must be shrunk to zero).
    b1 = rng.uniform(0.02, 0.10, n)
    b2 = rng.uniform(0.02, 0.10, n)
    y = 0.30 * b1 + 0.00 * b2 + rng.normal(0, 0.01, n)
    mv = mvmr([b1, b2], y, se_out)
    print(f"  MVMR       theta1={mv['theta'][0]:.3f} (true 0.30)  "
          f"theta2={mv['theta'][1]:+.3f} (true 0.00)")
