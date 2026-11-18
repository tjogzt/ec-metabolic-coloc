"""Numerical primitives shared by colocalisation stages."""

from __future__ import annotations

import math

import numpy as np


def logsum(x: np.ndarray) -> float:
    """Numerically stable log(sum(exp(x))).

    Shifting by the maximum before exponentiation prevents underflow when the
    log Bayes factors are very negative, which is the normal case for variants
    that show no association. An empty or non-finite input returns ``-inf``.
    """
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return -np.inf
    m = float(np.max(x))
    if not np.isfinite(m):
        return -np.inf
    return float(m + math.log(np.sum(np.exp(x - m))))
