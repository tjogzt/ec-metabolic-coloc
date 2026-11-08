"""
p0_concordance_v2.py — P0-2 修复：方向一致率改用与结论相称的统计量。

问题：原口径把 862 个基因全部当作独立伯努利试验。但验证队列里绝大多数基因
没有可测信号，其符号近似掷硬币，会把一致率人为拉向 50%。

改为三个统计量并报：
  (1) 全基因符号一致率（原口径，保留以便对照）
  (2) 响应基因子集（验证队列名义 p<0.05）的一致率 —— 主统计量
  (3) 全部基因 logFC 的秩相关（连续量，不受符号二分影响）
"""
from __future__ import annotations

import os as _os
import sys as _sys
import pathlib as _pathlib

_root = _pathlib.Path(__file__).resolve()
while not (_root / "config.py").exists():
    if _root.parent == _root:
        raise RuntimeError("config.py not found; set ECIR_ROOT")
    _root = _root.parent
if str(_root) not in _sys.path:
    _sys.path.insert(0, str(_root))
import config as CFG
import pathlib
import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr

W4 = CFG.ROOT
prog = pd.read_csv(W4 / "results/IR_state_program_v2.tsv", sep="\t", index_col=0)

COHORTS = [
    ("GSE13070 muscle, clamp-defined IR vs IS", "DE_GSE13070_muscle_baseline_IR.tsv"),
    ("GSE18732 muscle, type 2 diabetes vs NGT", "DE_GSE18732_muscle_DM_state.tsv"),
    ("GSE15773 adipose, IR vs IS", "DE_GSE15773_adipose_IR_state.tsv"),
]


def wilson(k, n, z=1.96):
    p = k / n
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return c - h, c + h


rows = []
print(f"{'队列':46s} {'全基因一致率':>14s} {'响应基因(p<.05)':>18s} {'秩相关':>18s}")
for label, fn in COHORTS:
    de = pd.read_csv(W4 / fn, sep="\t", index_col=0)
    m = prog[["logFC_adipose_obese"]].join(
        de[["logFC", "P.Value", "adj.P.Val"]], how="inner").dropna()
    m["same"] = np.sign(m.logFC_adipose_obese) == np.sign(m.logFC)

    k, n = int(m.same.sum()), len(m)
    lo, hi = wilson(k, n)
    p_all = binomtest(k, n, 0.5, alternative="greater").pvalue

    resp = m[m["P.Value"] < 0.05]
    k2, n2 = int(resp.same.sum()), len(resp)
    if n2:
        lo2, hi2 = wilson(k2, n2)
        p_resp = binomtest(k2, n2, 0.5, alternative="greater").pvalue
    else:
        lo2 = hi2 = p_resp = np.nan

    rho, p_rho = spearmanr(m.logFC_adipose_obese, m.logFC)

    print(f"{label:46s} {100*k/n:5.1f}% ({k}/{n}) {100*k2/n2:7.1f}% ({k2}/{n2}) "
          f"rho={rho:+.3f} (p={p_rho:.1e})")
    rows.append(dict(cohort=label, n_overlap=n,
                     conc_all=100*k/n, ci_all_lo=100*lo, ci_all_hi=100*hi, p_all=p_all,
                     n_responsive=n2, conc_resp=100*k2/n2 if n2 else np.nan,
                     ci_resp_lo=100*lo2 if n2 else np.nan,
                     ci_resp_hi=100*hi2 if n2 else np.nan,
                     p_resp=p_resp, rho=rho, p_rho=p_rho,
                     fdr_sig_in_program=int((m["adj.P.Val"] < 0.05).sum()),
                     fdr_sig_genomewide=int((de["adj.P.Val"] < 0.05).sum())))

out = pd.DataFrame(rows)
f = W4 / "results/concordance_v2.tsv"
out.to_csv(f, sep="\t", index=False)
print(f"\n已写 {f}")
print(out[["cohort", "conc_all", "conc_resp", "rho", "fdr_sig_in_program",
           "fdr_sig_genomewide"]].to_string(index=False))
