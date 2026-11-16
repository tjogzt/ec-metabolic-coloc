#!/usr/bin/env python3
"""R1-2/R3-2 补充分析：代理层的两个缺失统计对象。
(1) 基因块 bootstrap 95% CI（对程序基因重抽样，重算 Spearman ρ）；
(2) 100 套表达匹配随机模块的零分布（与观测同估计量），报零分布统计量与经验 p。
输出：w4/results/proxy_uncertainty.tsv
"""
import pathlib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

W4 = pathlib.Path(__file__).resolve().parent
RNG = np.random.default_rng(20261116)
N_BOOT = 5000

prog = pd.read_csv(W4 / "results/IR_state_program_v2.tsv", sep="\t", index_col=0)
disco_fc = prog["logFC_adipose_obese"].dropna()
disco_all = pd.read_csv(W4 / "DE_GSE25401_adipose_obese_state.tsv", sep="\t", index_col=0)["logFC"]
n_prog = len(disco_fc)
print(f"programme genes with discovery logFC: {n_prog}")

COHORTS = [
    ("GSE13070 muscle, clamp-defined IR vs IS", "DE_GSE13070_muscle_baseline_IR.tsv"),
    ("GSE18732 muscle, type 2 diabetes vs NGT", "DE_GSE18732_muscle_DM_state.tsv"),
    ("GSE15773 adipose, IR vs IS", "DE_GSE15773_adipose_IR_state.tsv"),
]

rows = []
for label, fn in COHORTS:
    de = pd.read_csv(W4 / fn, sep="\t", index_col=0)
    m = pd.DataFrame({"fc_d": disco_fc}).join(de[["logFC"]], how="inner").dropna()
    print(f"{label}: overlap {len(m)}")
    fc_d = m["fc_d"].values
    fc_c = m["logFC"].values

    obs_rho, obs_p = spearmanr(fc_d, fc_c)

    # (1) 基因块 bootstrap：对程序基因有放回重抽样
    idx = np.arange(len(m))
    boots = np.empty(N_BOOT)
    for b in range(N_BOOT):
        ib = RNG.integers(0, len(m), size=len(m))
        r, _ = spearmanr(fc_d[ib], fc_c[ib])
        boots[b] = r
    lo, hi = np.percentile(boots, [2.5, 97.5])

    # (2) 匹配随机模块零分布：100 套等大模块
    mods = pd.read_csv(
        W4.parent.parent / "ec-metabolic-coloc/results/tables/random_modules_v2.tsv",
        sep="\t")
    null_rhos = []
    for module, g in mods.groupby("module"):
        genes = g["gene"].values
        sub = pd.DataFrame({"fc_d": disco_all.reindex(genes)}).join(
            de[["logFC"]], how="inner").dropna()
        if len(sub) < 50:
            continue
        r, _ = spearmanr(sub["fc_d"], sub["logFC"])
        null_rhos.append(r)
    null_rhos = np.array(null_rhos)
    p_emp = float((null_rhos >= obs_rho).mean())  # 单侧上尾，与观测方向一致

    rows.append(dict(
        cohort=label, n=len(m),
        rho=obs_rho, p_rho=obs_p,
        boot_lo=lo, boot_hi=hi,
        null_mean=null_rhos.mean(), null_sd=null_rhos.std(),
        null_p95=np.percentile(null_rhos, 95),
        null_n=len(null_rhos), p_emp=p_emp))
    print(f"  rho={obs_rho:+.4f} (p={obs_p:.3g}) boot CI [{lo:+.3f},{hi:+.3f}] "
          f"null mean {null_rhos.mean():+.4f} SD {null_rhos.std():.4f} "
          f"p95 {np.percentile(null_rhos,95):+.4f} p_emp={p_emp:.3g}")

out = pd.DataFrame(rows)
out.to_csv(W4 / "results/proxy_uncertainty.tsv", sep="\t", index=False)
print(f"\nwritten {W4 / 'results/proxy_uncertainty.tsv'}")
