"""
coloc_abf.py — Colocalisation ABF (approximate Bayes factor) 的纯 Python 实现

参考：Giambartolomei et al. 2014 PLoS Genet（coloc）；Wakefield 2009 Genet Epidemiol（ABF）。
与 R 包 coloc::coloc.abf 的默认设置一致：
  - 效应量先验标准差 = 0.15 * sd(Y)
  - 先验概率 p1 = p2 = 1e-4, p12 = 1e-5
  - 假定每个性状至多一个因果变异

用法：
    res = coloc_abf(df, sdY1=1.0, sdY2=1.0)
其中 df 需包含列：beta1, se1, beta2, se2（同一批 SNP，已对齐效应等位基因）。
"""
from __future__ import annotations
import math
from typing import Iterable

import numpy as np

P1 = 1e-4      # 先验：性状 1 有信号
P2 = 1e-4      # 先验：性状 2 有信号
P12 = 1e-5     # 先验：两者共享同一因果变异


def log_abf(beta: np.ndarray, se: np.ndarray, sd_y: float, prior_sd: float = 0.15) -> np.ndarray:
    """单个 SNP 的对数近似贝叶斯因子（Wakefield ABF）。"""
    beta = np.asarray(beta, dtype=float)
    se = np.asarray(se, dtype=float)
    W = (prior_sd * sd_y) ** 2       # 先验方差
    V = se ** 2                      # 抽样方差
    r = W / (W + V)
    z2 = (beta / se) ** 2
    # Wakefield ABF: ABF = sqrt(V/(V+W)) * exp(z^2 * W / (2(V+W)))
    #               = sqrt(1-r) * exp(r*z^2/2)
    return 0.5 * np.log(1 - r) + 0.5 * r * z2


def _logsum(x: np.ndarray) -> float:
    """数值稳定的 log(sum(exp(x)))。"""
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return -np.inf
    m = np.max(x)
    if not np.isfinite(m):
        return -np.inf
    return m + math.log(np.sum(np.exp(x - m)))


def coloc_abf(l1: np.ndarray, l2: np.ndarray, p1: float = P1, p2: float = P2,
              p12: float = P12) -> dict:
    """
    由两套对数 ABF 计算 PP.H0--H4。

    l1, l2: 同一批 SNP 在两个性状上的 log ABF。
    """
    l1 = np.asarray(l1, dtype=float)
    l2 = np.asarray(l2, dtype=float)
    assert l1.shape == l2.shape and l1.size > 0

    # H1: 仅性状 1 关联（对 SNP 求和后乘以 "仅一个性状有关联" 的先验）
    lh1 = _logsum(l1) + math.log(p1) + math.log(1 - p2)
    lh2 = _logsum(l2) + math.log(p2) + math.log(1 - p1)
    # H3: 两者都关联但因果变异不同（两个独立集合，互斥）
    off = np.add.outer(l1, l2)
    mask = np.eye(len(l1), dtype=bool)
    off = np.where(mask, -np.inf, off)          # 排除 i == j（那是 H4 的定义）
    lh3 = math.log(p1) + math.log(p2) + _logsum(off.ravel())
    # H4: 共享同一因果变异
    lh4 = _logsum(l1 + l2) + math.log(p12)
    # H0
    lh0 = math.log(1 - p1 - p2 - p12)

    all_l = np.array([lh0, lh1, lh2, lh3, lh4])
    denom = _logsum(all_l)
    pp = np.exp(all_l - denom)
    return {"PP.H0": float(pp[0]), "PP.H1": float(pp[1]), "PP.H2": float(pp[2]),
            "PP.H3": float(pp[3]), "PP.H4": float(pp[4]),
            "n_snps": int(len(l1)),
            "top_snp_share": float(np.max(np.exp(l1 + l2 - _logsum(l1 + l2))))}


def coloc_from_stats(beta1, se1, beta2, se2, sd_y1: float, sd_y2: float, **kw) -> dict:
    """便捷入口：直接由 beta/se 计算。"""
    l1 = log_abf(beta1, se1, sd_y1)
    l2 = log_abf(beta2, se2, sd_y2)
    return coloc_abf(l1, l2, **kw)
