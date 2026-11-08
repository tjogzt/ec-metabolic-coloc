"""
screen.py — 一级筛选：GTEx 子宫 eQTL 的 lead 变异是否在 EC GWAS 中呈关联

这一步不是 colocalisation，只是“同一个变异是否同时影响表达与疾病”的初筛
（相当于 SMR 的无 LD 近似）。它的作用是：
  1) 把 3,577 个子宫 eGene 压到几十到几百个候选；
  2) 为后续真正的区域 coloc 提供候选区域清单。

用法：
  python screen.py --tissue Uterus --out results/screen_uterus.tsv [--limit 100]
"""
from __future__ import annotations
import argparse, json, math, pathlib, sys, time

import pandas as pd
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from io_gtex import load_tissue, split_variant, lead_per_gene
from io_finngen import FinnGen, harmonise

DATA = pathlib.Path(__file__).resolve().parents[1] / "data"
GWAS = DATA / "finngen_R12_C3_CORPUS_UTERI_EXALLC.gz"
TAR = DATA / "GTEx_Analysis_v8_eQTL.tar"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", default="Uterus")
    ap.add_argument("--out", default="results/screen_uterus.tsv")
    ap.add_argument("--limit", type=int, default=0, help="仅测前 N 个基因（调试用）")
    ap.add_argument("--min_af", type=float, default=0.01, help="跳过极低频变异")
    a = ap.parse_args()

    t0 = time.time()
    print(f"[1/4] 读取 {a.tissue} 显著 eQTL 对 …")
    pairs = split_variant(load_tissue(TAR, a.tissue))
    print(f"      显著对 {len(pairs):,} 条；基因 {pairs['gene_id'].nunique():,} 个")

    leads = lead_per_gene(pairs)
    leads = leads[(leads["maf"].between(a.min_af, 1 - a.min_af)) & leads["pos"].notna()]
    if a.limit:
        leads = leads.head(a.limit)
    print(f"[2/4] 每个基因取 lead 变异：{len(leads):,} 个")

    fg = FinnGen(GWAS)
    rows = []
    for k, r in enumerate(leads.itertuples(index=False), 1):
        hit = fg.lookup(r.chrom, int(r.pos))
        if hit.empty:
            rows.append({"gene_id": r.gene_id, "variant_id": r.variant_id,
                         "eQTL_slope": r.slope, "eQTL_p": r.pval_nominal,
                         "in_gwas": False})
            continue
        h = harmonise(hit, r.ref, r.alt)
        if h.empty:
            rows.append({"gene_id": r.gene_id, "variant_id": r.variant_id,
                         "eQTL_slope": r.slope, "eQTL_p": r.pval_nominal,
                         "in_gwas": True, "allele_mismatch": True})
            continue
        b = h.iloc[0]
        rows.append({"gene_id": r.gene_id, "variant_id": r.variant_id,
                     "eQTL_slope": r.slope, "eQTL_p": r.pval_nominal,
                     "in_gwas": True, "gwas_beta": b.beta_h, "gwas_se": b.se_h,
                     "gwas_p": b.pval, "gwas_af": b.af_alt})
        if k % 500 == 0:
            print(f"      … {k}/{len(leads)}  ({time.time()-t0:.0f}s)")

    df = pd.DataFrame(rows)
    df["in_gwas"] = df["in_gwas"].fillna(False).astype(bool)
    sig = df[(df["in_gwas"]) & (df.get("gwas_p", pd.Series(dtype=float)) < 0.05)]
    conc = ((sig["eQTL_slope"] * sig["gwas_beta"]) > 0).sum()

    out = pathlib.Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)

    print(f"\n[3/4] 写表：{out}  ({len(df):,} 行)")
    print(f"[4/4] 汇总")
    print(f"  基因数                 {len(df):,}")
    print(f"  lead 变异在 GWAS 中     {int(df['in_gwas'].sum()):,}")
    print(f"  GWAS p < 0.05          {len(sig):,}")
    if len(sig):
        print(f"  方向一致               {int(conc):,} / {len(sig)}  "
              f"({conc/len(sig)*100:.1f}%)")
        p_dir = stats.binomtest(int(conc), len(sig), 0.5).pvalue
        print(f"  方向一致性二项检验 p   {p_dir:.4g}")
    print(f"  用时 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
