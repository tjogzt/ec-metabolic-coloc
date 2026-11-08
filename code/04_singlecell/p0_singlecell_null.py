"""
p0_singlecell_null.py — P0-3 修复：零分布必须与被检验统计量同构。

原实现：观测统计量 = score(UP臂) − score(DOWN臂)（差分），
        零分布      = 单个随机模块的均值（不分方向）。
两者位置与尺度都不同，故 p_emp 不可解释——把"对比"拿去比"整体位移"。

本脚本按同一统计量重建零分布：每个随机模块随机拆成与程序等大的 UP/DOWN 两臂
（v2 为 608/254），算同一个差分，再在样本层取（恶性均值 − 正常均值）。

用法：
  PROG_PATH=...w4/results/IR_state_program_v2.tsv \
  RM_PATH=...w2/results/random_modules_v2.tsv \
  OUT=...w3/results/singlecell_null_matched_v2.tsv \
  /opt/anaconda3/bin/python3 w3/p0_singlecell_null.py
"""
from __future__ import annotations
import os, pathlib
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from scipy.stats import mannwhitneyu

HERE = pathlib.Path(__file__).resolve().parent
W2 = HERE.parent / "w2"
H5 = os.environ.get(
    "H5",
    "/Volumes/tjogzt4T/data/scrnaseq_reference/EC_REF_v1/ec_integrated_v1.h5ad")

CYCLING = ["MKI67", "PCNA", "CCNB1", "CCNB2", "CCNA2", "CDK1", "CDC20", "BUB1", "BUB1B",
           "TOP2A", "AURKA", "AURKB", "PLK1", "TPX2", "UBE2C", "RRM2", "TYMS", "MCM2",
           "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "CENPF", "KIF11", "KIF23", "ANLN",
           "PRC1", "ECT2", "NCAPG", "ASPM", "DLGAP5", "CENPE", "SMC4", "NUSAP1"]

PROG = os.environ.get("PROG_PATH", str(HERE.parent / "w4/results/IR_state_program_v2.tsv"))
RM = os.environ.get("RM_PATH", str(W2 / "results/random_modules_v2.tsv"))
OUT = os.environ.get("OUT", str(HERE / "results/singlecell_null_matched.tsv"))
SEED = 49


def main() -> None:
    prog = pd.read_csv(PROG, sep="\t", index_col=0)
    rm = pd.read_csv(RM, sep="\t")
    up = list(prog.index[prog.direction == "UP"])
    dn = list(prog.index[prog.direction == "DOWN"])
    print(f"程序：UP {len(up)} / DOWN {len(dn)}；随机模块 {rm.module.nunique()} 套")

    a = ad.read_h5ad(H5, backed="r")
    obs = a.obs[["sample", "sample_class", "coarse_label"]].copy()
    epi_mask = obs.coarse_label.isin(["MalignantEpithelial", "NormalEpithelial"]).values
    epi = np.where(epi_mask)[0]
    print(f"{a.n_obs:,} 细胞 × {a.n_vars:,} 基因；上皮 {len(epi):,}")

    mods = {k: list(g.gene) for k, g in rm.groupby("module")}
    need = set(up) | set(dn)
    for v in mods.values():
        need |= set(v)
    idx = a.var_names.get_indexer(list(need))
    keep = idx >= 0
    genes_ok = [g for g, k in zip(list(need), keep) if k]
    X = sp.csr_matrix(a[:, genes_ok].to_memory().X)
    mu = np.asarray(X.mean(axis=0)).ravel()
    sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
    sd = np.sqrt(np.maximum(sq - mu**2, 1e-8))
    gidx = {g: i for i, g in enumerate(genes_ok)}
    print(f"可用基因 {len(genes_ok):,}")

    def arm_score(gl):
        gj = [gidx[g] for g in gl if g in gidx]
        if len(gj) < 5:
            return None
        w = 1.0 / (sd[gj] * len(gj))
        c = float(np.sum(mu[gj] / sd[gj]) / len(gj))
        return np.asarray((X[:, gj] @ w).ravel()) - c

    up_s, dn_s = arm_score(up), arm_score(dn)
    if up_s is None or dn_s is None:
        raise SystemExit("程序某一臂基因不足")
    ir = up_s - dn_s
    n_up, n_dn = len([g for g in up if g in gidx]), len([g for g in dn if g in gidx])
    print(f"实际进入打分的臂：UP {n_up} / DOWN {n_dn}")

    def sample_delta(score):
        ps = (pd.DataFrame({"v": score[epi], "s": obs.iloc[epi]["sample"].values,
                            "c": obs.iloc[epi]["sample_class"].values})
              .groupby(["s", "c"], observed=True)["v"].mean().dropna().reset_index())
        mal = ps.loc[ps.c == "Malignant", "v"].dropna().values
        nor = ps.loc[ps.c == "Normal", "v"].dropna().values
        if mal.size == 0 or nor.size == 0:
            return np.nan, mal, nor
        return mal.mean() - nor.mean(), mal, nor

    obs_d, mal, nor = sample_delta(ir)
    u = mannwhitneyu(mal, nor, alternative="two-sided")
    print(f"\n观测（v2，匹配统计量）：Δ={obs_d:+.4f}  MWU p={u.pvalue:.4f}")

    rng = np.random.default_rng(SEED)
    null = []
    ratio = n_up / (n_up + n_dn)          # 程序两臂的实际比例
    for nm, gl in mods.items():
        g = [x for x in gl if x in gidx]
        if len(g) < 100:
            continue
        # 模块可用基因数可能少于程序（表达过滤后），按同一比例拆分，
        # 保证两臂相对大小与程序一致，而不是放弃该模块。
        k_up = int(round(len(g) * ratio))
        k_dn = len(g) - k_up
        if k_up < 20 or k_dn < 20:
            continue
        perm = rng.permutation(len(g))
        a_up = [g[i] for i in perm[:k_up]]
        a_dn = [g[i] for i in perm[k_up:]]
        su, sdw = arm_score(a_up), arm_score(a_dn)
        if su is None or sdw is None:
            continue
        d, _, _ = sample_delta(su - sdw)
        null.append(d)

    null = np.array(null)
    if null.size < 10:
        raise SystemExit(f"有效随机模块只有 {null.size} 套，不足以构建零分布")
    p_upper = float((null >= obs_d).mean())
    p_lower = float((null <= obs_d).mean())
    p_two = float(min(1.0, 2 * min(p_upper, p_lower)))
    z = (obs_d - null.mean()) / null.std(ddof=1)
    print(f"匹配零分布（{len(null)} 套）：均值 {null.mean():+.4f}，SD {null.std(ddof=1):.4f}，"
          f"95% 分位 {np.percentile(null, 95):+.4f}")
    print(f"观测 z = {z:+.2f}；上尾 p={p_upper:.3f}，下尾 p={p_lower:.3f}，双侧 p={p_two:.3f}")

    pd.DataFrame([dict(observed=obs_d, null_mean=null.mean(), null_sd=null.std(ddof=1),
                       null_p95=float(np.percentile(null, 95)), z=z,
                       p_upper=p_upper, p_lower=p_lower, p_two_sided=p_two,
                       n_null=len(null), n_up=n_up, n_down=n_dn,
                       mwu_p=float(u.pvalue))]).to_csv(OUT, sep="\t", index=False)
    np.savetxt(str(OUT).replace(".tsv", "_values.tsv"), null, header="null_delta")
    print(f"\n已写 {OUT}")


if __name__ == "__main__":
    main()
