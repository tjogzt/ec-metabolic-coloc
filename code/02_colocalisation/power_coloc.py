"""
power_coloc.py — 共定位的统计功效曲线：eQTL 强度 → 能否达到 PP.H4>0.8

背景：GTEx 子宫 eGene 的最强 eQTL z 中位数只有 5.28（90% 分位 8.29，
仅 11.6% 超过 8）。若真实共享因果变异存在，这样的 eQTL 强度还能不能把
PP.H4 推到 0.8？本脚本用真实变异集做蒙特卡洛回答这个问题。

做法：取真实区域的全部变异（真实位置／等位基因／se），
      在**同一个** SNP 上同时放置 eQTL 与 GWAS 的真实效应（共享因果），
      只改变 eQTL 的 z 强度与 GWAS 的 z 强度，重复 200 次，统计：
        - PP.H4 的中位数与四分位
        - 达到 PP.H4>0.8 的比例（= 检出率）
这是"该组织能否回答该问题"的直接证据，而不是"有没有共享因果"。
"""
from __future__ import annotations
import pathlib, sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from coloc_eqtlcat import qtl_region
from coloc_abf import coloc_abf, log_abf

DATA = pathlib.Path(__file__).resolve().parents[1] / "data"
EQTL = DATA / "QTD000346.all.tsv.gz"
CHROM, LO, HI = "8", 38_600_000, 39_700_000
GENE = "ENSG00000197140"
RNG = np.random.default_rng(7)


def load_region() -> pd.DataFrame:
    d = qtl_region(CHROM, LO, HI, EQTL)
    d = d[d.gene_id == GENE].dropna(subset=["position", "beta", "se", "pvalue"])
    d = d[(d.maf >= 0.01) & (d.maf <= 0.99) & (d.se > 0)]
    d = d.rename(columns={"position": "pos"})
    return d.sort_values("pos").reset_index(drop=True)


def one_trial(d: pd.DataFrame, z_eqtl: float, z_gwas: float, causal_idx: int) -> dict:
    """在同一 SNP 上放两个真实效应（共享因果），其余为噪声。"""
    n = len(d)
    se_e = d["se"].values.copy()
    be = RNG.normal(0, se_e)
    be[causal_idx] = z_eqtl * se_e[causal_idx]

    se_g = se_e * RNG.uniform(0.7, 1.3, n)
    bg = RNG.normal(0, se_g)
    bg[causal_idx] = z_gwas * se_g[causal_idx]

    r = coloc_abf(log_abf(be, se_e, 1.0), log_abf(bg, se_g, 1.0))
    r["n"] = n
    return r


def main():
    d = load_region()
    causal_idx = int(np.argmin(d["pvalue"].values))   # 真实的 eQTL lead
    print(f"区域 {GENE}：{len(d):,} 个变异，共享因果放在 {int(d['pos'].iloc[causal_idx]):,}\n")

    z_grid = [4, 5.3, 6.5, 8, 10, 13]
    rows = []
    for z_e in z_grid:
        for z_g in (10.0, 15.0, 25.0):
            h4 = [one_trial(d, z_e, z_g, causal_idx)["PP.H4"] for _ in range(60)]
            h3 = [one_trial(d, z_e, z_g, causal_idx)["PP.H3"] for _ in range(60)]
            rows.append({"eQTL_z": z_e, "GWAS_z": z_g,
                         "H4_中位": round(float(np.median(h4)), 3),
                         "H4_P25": round(float(np.percentile(h4, 25)), 3),
                         "H4_P75": round(float(np.percentile(h4, 75)), 3),
                         "检出率_H4>0.8": f"{100*np.mean(np.array(h4) > 0.8):.0f}%",
                         "H3_中位": round(float(np.median(h3)), 3)})
    out = pd.DataFrame(rows)
    print(out.to_string(index=False))
    p = pathlib.Path("results/coloc_power.tsv")
    out.to_csv(p, sep="\t", index=False)
    print(f"\n写表 {p}")


if __name__ == "__main__":
    main()
