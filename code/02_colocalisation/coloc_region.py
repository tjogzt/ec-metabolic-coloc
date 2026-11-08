"""
coloc_region.py — 按基因做区域共定位（受限 SNP 集版本）

输入：GTEx v8 某组织的全部显著 variant-gene pairs（每基因数十至数百个变异）、
      给定 GWAS（FinnGen tabix）。
对每个基因：
  1. 取其显著变异的最小--最大坐标，一次 tabix 拉整个区域；
  2. 按 (pos, ref, alt) 与 eQTL 变异对齐，统一效应方向；
  3. 以 eQTL 显著变异集合为 SNP 集，计算 ABF 并给出 PP.H0--H4。

⚠ 口径说明（必须写进 Methods）：
   标准 coloc 用区域内**全部被检验过的变异**；这里只有“显著对”，SNP 集被截断，
   PP.H3 与 PP.H4 的绝对值会有偏。因此本步骤用于**排序与候选择优**，
   并以“置换零分布”校准假阳性率（见 --null 参数），不单独作为因果证据。

用法：
  python coloc_region.py --tissue Uterus --out results/coloc_uterus.tsv
  python coloc_region.py --tissue Uterus --null 200          # 置换零分布
"""
from __future__ import annotations
import argparse, pathlib, sys, time

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from io_gtex import load_tissue, split_variant
from io_finngen import FinnGen, harmonise
from coloc_abf import coloc_abf, log_abf

DATA = pathlib.Path(__file__).resolve().parents[1] / "data"
GWAS = DATA / "finngen_R12_C3_CORPUS_UTERI_EXALLC.gz"
TAR = DATA / "GTEx_Analysis_v8_eQTL.tar"


def build_snp_table(gene_df: pd.DataFrame, fg: FinnGen):
    """把某基因的 eQTL 显著变异与 GWAS 区域结果合并，返回对齐后的 beta/se。"""
    lo, hi = int(gene_df["pos"].min()), int(gene_df["pos"].max())
    chrom = str(gene_df["chrom"].iloc[0])
    # 扩边 10kb，保证能取到同位置的所有 GWAS 记录
    g = fg.region(chrom, max(1, lo - 10000), hi + 10000)
    if g.empty:
        return None
    g = g.copy()
    g["key"] = g["pos"].astype(str) + "_" + g["ref"] + "_" + g["alt"]
    g["key_swap"] = g["pos"].astype(str) + "_" + g["alt"] + "_" + g["ref"]
    lut = {}
    for r in g.itertuples(index=False):
        lut[r.key] = (r.beta, r.sebeta, r.pval)
        lut[r.key_swap] = (-r.beta, r.sebeta, r.pval)

    rows = []
    for r in gene_df.itertuples(index=False):
        k = f"{int(r.pos)}_{r.ref}_{r.alt}"
        if k not in lut:
            continue
        b, se, p = lut[k]
        if not np.isfinite(b) or not np.isfinite(se) or se <= 0:
            continue
        rows.append((r.variant_id, r.slope, r.slope_se, b, se, p))
    if len(rows) < 2:
        return None
    return pd.DataFrame(rows, columns=["variant_id", "b_e", "se_e", "b_g", "se_g", "p_g"])


def coloc_one(df: pd.DataFrame, sd_e=1.0, sd_g=1.0):
    l1 = log_abf(df["b_e"].values, df["se_e"].values, sd_e)
    l2 = log_abf(df["b_g"].values, df["se_g"].values, sd_g)
    res = coloc_abf(l1, l2)
    res["n_snps"] = len(df)
    res["best_gwas_p"] = float(df["p_g"].min())
    top = df.loc[df["p_g"].idxmin()]
    res["top_sign_concordant"] = bool(top["b_e"] * top["b_g"] > 0)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", default="Uterus")
    ap.add_argument("--out", default="results/coloc_uterus.tsv")
    ap.add_argument("--min_snps", type=int, default=3)
    ap.add_argument("--null", type=int, default=0, help="置换次数（0=不跑零分布）")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    t0 = time.time()
    pairs = split_variant(load_tissue(TAR, a.tissue))
    pairs = pairs[pairs["pos"].notna()]
    genes = pairs.groupby("gene_id")
    keys = list(genes.groups.keys())
    if a.limit:
        keys = keys[: a.limit]
    print(f"[1/3] {a.tissue}：{len(pairs):,} 显著对 / {pairs.gene_id.nunique():,} 基因"
          f"；本次处理 {len(keys):,}")

    fg = FinnGen(GWAS)
    recs = []
    for i, gid in enumerate(keys, 1):
        d = build_snp_table(genes.get_group(gid), fg)
        if d is None or len(d) < a.min_snps:
            continue
        try:
            r = coloc_one(d)
        except Exception:
            continue
        r["gene_id"] = gid
        recs.append(r)
        if i % 400 == 0:
            print(f"      … {i}/{len(keys)}  ({time.time()-t0:.0f}s)")

    res = pd.DataFrame(recs)
    out = pathlib.Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out, sep="\t", index=False)
    print(f"\n[2/3] 写表 {out}：{len(res):,} 个基因完成 coloc")

    for thr in (0.5, 0.8, 0.9):
        print(f"      PP.H4 > {thr}: {int((res['PP.H4'] > thr).sum()):,}")
    print(f"      PP.H3 > 0.8: {int((res['PP.H3'] > 0.8).sum()):,}")
    print(f"      PP.H4 中位 {res['PP.H4'].median():.3f}；"
          f"n_snps 中位 {res['n_snps'].median():.0f}")

    if a.null:
        print(f"\n[3/3] 置换零分布（{a.null} 次，打乱 eQTL 基因与 SNP 的配对）")
        rng = np.random.default_rng(49)
        null_h4 = []
        sub = res.head(min(300, len(res)))
        pre = []
        for gid in sub["gene_id"]:
            d = build_snp_table(genes.get_group(gid), fg)
            if d is not None and len(d) >= a.min_snps:
                pre.append(d)
        for k in range(a.null):
            hits = 0
            for d in pre:
                dd = d.copy()
                dd["b_e"] = rng.permutation(dd["b_e"].values)
                dd["se_e"] = rng.permutation(dd["se_e"].values)
                r = coloc_one(dd)
                hits += r["PP.H4"] > 0.8
            null_h4.append(hits / max(1, len(pre)))
        print(f"      零分布 PP.H4>0.8 比例：均值 {np.mean(null_h4):.4f}，"
              f"95% 分位 {np.quantile(null_h4, 0.95):.4f}")
        pd.Series(null_h4).to_csv(out.with_suffix(".null.csv"), index=False)

    print(f"\n用时 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
