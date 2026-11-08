"""
coloc_perm.py — 置换校准的区域共定位

背景：直接用 GTEx **显著变异集**做 coloc 会系统性高估 PP.H4（实测：观测 24.0%，
置换零分布 20.8%，95% 分位 23.2%）。原因是 SNP 集被截断——留下的都是 eQTL 强关联位点，
任何一个偶然的 GWAS 波动都更容易落在同一批 SNP 上。

做法：对每个基因，在其**自身 SNP 集内**打乱 GWAS 效应（保持边际分布），
得到该基因条件零分布下的 PP.H4；用观测值在零分布中的位置给出经验 p 值。
这样得到的排序在“给定该 SNP 集”的条件下是校准的，可以跨基因比较。

⚠ 仍需在 Methods 声明：这是条件于受限 SNP 集的置换检验，检验力低于标准 coloc，
   零假设是“在该 SNP 集上，eQTL 与 GWAS 效应无配对关系”。

用法：
  python coloc_perm.py --tissue Uterus --B 50 --out results/g_a_candidates.tsv
"""
from __future__ import annotations
import argparse, pathlib, sys, time

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from io_gtex import load_tissue, split_variant
from io_finngen import FinnGen
from coloc_abf import coloc_abf, log_abf
from coloc_region import build_snp_table

DATA = pathlib.Path(__file__).resolve().parents[1] / "data"
GWAS = DATA / "finngen_R12_C3_CORPUS_UTERI_EXALLC.gz"
TAR = DATA / "GTEx_Analysis_v8_eQTL.tar"


def h4(df, l1_cache=None, rng=None, permute=False, mode="shift"):
    b, se = df["b_e"].values, df["se_e"].values
    bg, sg = df["b_g"].values, df["se_g"].values
    if permute:
        if mode == "shift":
            # 环形平移：保留 GWAS 信号的空间聚集结构（LD 簇），只打断与 eQTL 的配对
            k = int(rng.integers(1, len(bg)))
            bg, sg = np.roll(bg, k), np.roll(sg, k)
        else:
            idx = rng.permutation(len(bg))
            bg, sg = bg[idx], sg[idx]
    l1 = log_abf(b, se, 1.0)
    l2 = log_abf(bg, sg, 1.0)
    return coloc_abf(l1, l2)["PP.H4"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", default="Uterus")
    ap.add_argument("--B", type=int, default=50)
    ap.add_argument("--min_snps", type=int, default=5)
    ap.add_argument("--out", default="results/g_a_candidates.tsv")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--null_mode", default="shift", choices=["shift", "perm"],
                    help="shift=按位置环形平移（保留 LD 簇结构，推荐）；perm=完全打乱")
    a = ap.parse_args()

    t0 = time.time()
    rng = np.random.default_rng(49)
    pairs = split_variant(load_tissue(TAR, a.tissue))
    pairs = pairs[pairs["pos"].notna()]
    genes = pairs.groupby("gene_id")
    keys = list(genes.groups.keys())
    if a.limit:
        keys = keys[: a.limit]
    print(f"[1/2] {len(keys):,} 个基因，B={a.B} 次零分布（模式 {a.null_mode}）")

    fg = FinnGen(GWAS)
    rows = []
    for i, gid in enumerate(keys, 1):
        d = build_snp_table(genes.get_group(gid), fg)
        if d is None or len(d) < a.min_snps:
            continue
        d = d.sort_values("variant_id").reset_index(drop=True)
        obs = h4(d)
        null = np.array([h4(d, rng=rng, permute=True, mode=a.null_mode)
                         for _ in range(a.B)])
        p = (1 + int((null >= obs - 1e-12).sum())) / (a.B + 1)
        rows.append({"gene_id": gid, "n_snps": len(d), "PP.H4_obs": obs,
                     "null_mean": float(null.mean()),
                     "null_q95": float(np.quantile(null, 0.95)),
                     "excess": obs - float(null.mean()),
                     "p_perm": p,
                     "best_gwas_p": float(d["p_g"].min())})
        if i % 400 == 0:
            print(f"      … {i}/{len(keys)}  ({time.time()-t0:.0f}s)")

    res = pd.DataFrame(rows)
    # BH-FDR
    res = res.sort_values("p_perm").reset_index(drop=True)
    m = len(res)
    res["q_BH"] = (res["p_perm"] * m / (np.arange(m) + 1)).cummin().clip(upper=1)

    out = pathlib.Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out, sep="\t", index=False)
    print(f"\n[2/2] 写表 {out}：{m:,} 个基因")
    print(f"      p<0.05   {int((res.p_perm < 0.05).sum()):,}")
    print(f"      p<0.05/基因数 {int((res.p_perm < 0.05/m).sum()):,}")
    print(f"      q<0.10   {int((res.q_BH < 0.10).sum()):,}")
    print(f"      q<0.05   {int((res.q_BH < 0.05).sum()):,}")
    print(f"      用时 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
