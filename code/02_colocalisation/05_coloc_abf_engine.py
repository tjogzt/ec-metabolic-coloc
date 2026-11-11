"""
05_coloc_abf_engine.py
======================

Purpose
-------
Pure-NumPy implementation of the approximate Bayes factor (ABF) colocalisation
test of Giambartolomei et al. (2014), using the Wakefield (2009) ABF.

The module computes the five posterior probabilities PP.H0 through PP.H4 for a
pair of traits measured on the same set of variants:

    H0  no association with either trait
    H1  association with trait 1 only
    H2  association with trait 2 only
    H3  association with both traits, two distinct causal variants
    H4  association with both traits, one shared causal variant

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
log_abf()            beta, se and the trait standard deviation for one variant set
coloc_abf()          two aligned arrays of log ABF (same variants, same order)
coloc_from_stats()   beta/se arrays directly, for convenience

Prior probabilities are arguments, defaulting to the coloc package defaults
(p1 = p2 = 1e-4, p12 = 1e-5). Callers pass the values held in config.PARAMS.

Outputs
-------
A dict with keys PP.H0 ... PP.H4, n_snps and top_snp_share (the posterior mass
of the single most likely shared variant, a diagnostic for diffuse signals).

Notes
-----
The implementation mirrors coloc::coloc.abf: the effect-size prior standard
deviation is 0.15 * sd(Y), every trait is assumed to carry at most one causal
variant, and H3 is computed over ordered pairs with the diagonal excluded so
that it stays mutually exclusive with H4.

Reference
---------
Giambartolomei C, et al. PLoS Genet 2014;10:e1004383.
Wakefield J. Genet Epidemiol 2009;33:79-86.
"""
from __future__ import annotations

import math
from typing import Dict

import numpy as np

# Default priors, identical to coloc::coloc.abf.
P1_DEFAULT = 1e-4
P2_DEFAULT = 1e-4
P12_DEFAULT = 1e-5

# Prior standard deviation of the per-variant effect, as a fraction of sd(Y).
PRIOR_SD_FRACTION = 0.15


def log_abf(beta: np.ndarray, se: np.ndarray, sd_y: float,
            prior_sd: float = PRIOR_SD_FRACTION) -> np.ndarray:
    """Log approximate Bayes factor for one variant set (Wakefield 2009).

    Parameters
    ----------
    beta : array_like
        Effect size of each variant.
    se : array_like
        Standard error of each effect size; must be strictly positive.
    sd_y : float
        Standard deviation of the trait, used to scale the prior.
    prior_sd : float, optional
        Prior SD of the effect as a fraction of ``sd_y`` (default 0.15).

    Returns
    -------
    ndarray of log ABF, same shape as ``beta``.

    Derivation
    ----------
    With W = (prior_sd * sd_y)^2 and V = se^2, the Wakefield ABF is

        ABF = sqrt(V / (V + W)) * exp(z^2 * W / (2 * (V + W)))

    Writing r = W / (V + W), this is sqrt(1 - r) * exp(r * z^2 / 2), whose log is
    0.5 * log(1 - r) + 0.5 * r * z^2. The log form avoids overflow when z^2 is
    large, which happens routinely at genome-wide significant variants.
    """
    beta = np.asarray(beta, dtype=float)
    se = np.asarray(se, dtype=float)

    W = (prior_sd * sd_y) ** 2       # prior variance of the effect
    V = se ** 2                      # sampling variance
    r = W / (W + V)                  # shrinkage factor
    z2 = (beta / se) ** 2            # squared z statistic
    return 0.5 * np.log(1.0 - r) + 0.5 * r * z2


def _logsum(x: np.ndarray) -> float:
    """Numerically stable log(sum(exp(x))).

    Shifting by the maximum before exponentiation prevents underflow when the
    log Bayes factors are very negative, which is the normal case for variants
    that show no association.
    """
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return -np.inf
    m = np.max(x)
    if not np.isfinite(m):
        return -np.inf
    return float(m + math.log(np.sum(np.exp(x - m))))


def coloc_abf(l1: np.ndarray, l2: np.ndarray,
              p1: float = P1_DEFAULT, p2: float = P2_DEFAULT,
              p12: float = P12_DEFAULT) -> Dict[str, float]:
    """Posterior probabilities PP.H0 to PP.H4 from two aligned log-ABF arrays.

    Parameters
    ----------
    l1, l2 : array_like
        Log ABF of the same variants in trait 1 and trait 2. The caller is
        responsible for having harmonised the variants and their effect alleles
        beforehand; length and order must match.
    p1, p2 : float
        Prior probability that trait 1 (respectively trait 2) has a signal.
    p12 : float
        Prior probability that the two traits share one causal variant.

    Returns
    -------
    dict with PP.H0 ... PP.H4, n_snps and top_snp_share.

    Raises
    ------
    AssertionError
        If the two arrays differ in shape or are empty.
    """
    l1 = np.asarray(l1, dtype=float)
    l2 = np.asarray(l2, dtype=float)
    assert l1.shape == l2.shape and l1.size > 0, "log-ABF arrays must align"

    # H1 and H2: one trait has a signal, the other does not.
    lh1 = _logsum(l1) + math.log(p1) + math.log(1.0 - p2)
    lh2 = _logsum(l2) + math.log(p2) + math.log(1.0 - p1)

    # H3: both traits have a signal at two distinct variants. Sum over ordered
    # pairs, then remove the diagonal, which is the H4 configuration.
    off = np.add.outer(l1, l2)
    off = np.where(np.eye(len(l1), dtype=bool), -np.inf, off)
    lh3 = math.log(p1) + math.log(p2) + _logsum(off.ravel())

    # H4: one shared causal variant, so the log ABF adds across traits.
    lh4 = _logsum(l1 + l2) + math.log(p12)

    # H0: neither trait has a signal.
    lh0 = math.log(1.0 - p1 - p2 - p12)

    all_l = np.array([lh0, lh1, lh2, lh3, lh4], dtype=float)
    pp = np.exp(all_l - _logsum(all_l))

    # Posterior mass concentrated on the single best shared variant. A low
    # value means the H4 support is spread over many variants, which is worth
    # reporting alongside PP.H4 rather than hiding.
    shared = l1 + l2
    top_snp_share = float(np.max(np.exp(shared - _logsum(shared))))

    return {
        "PP.H0": float(pp[0]),
        "PP.H1": float(pp[1]),
        "PP.H2": float(pp[2]),
        "PP.H3": float(pp[3]),
        "PP.H4": float(pp[4]),
        "n_snps": int(len(l1)),
        "top_snp_share": top_snp_share,
    }


def coloc_from_stats(beta1, se1, beta2, se2,
                     sd_y1: float, sd_y2: float, **kwargs) -> Dict[str, float]:
    """Convenience wrapper: compute both log-ABF arrays, then colocalise.

    Parameters
    ----------
    beta1, se1 : array_like
        Effects and standard errors for trait 1.
    beta2, se2 : array_like
        Effects and standard errors for trait 2, aligned to the same variants.
    sd_y1, sd_y2 : float
        Trait standard deviations used to scale each prior.
    **kwargs
        Forwarded to :func:`coloc_abf` (p1, p2, p12).

    Returns
    -------
    The dict returned by :func:`coloc_abf`.
    """
    l1 = log_abf(beta1, se1, sd_y1)
    l2 = log_abf(beta2, se2, sd_y2)
    return coloc_abf(l1, l2, **kwargs)
