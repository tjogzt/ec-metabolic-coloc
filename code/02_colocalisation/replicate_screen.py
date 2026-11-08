"""
replicate_screen.py — 用 GCST006465（O'Mara 2018，8,758 例）复现子宫 eQTL lead 变异筛选

GTEx v8 是 GRCh38，GCST006465 是 GRCh37 → 先做 liftover（实测 3,577 个 lead 中
可降 3,5xx 个，往返一致率 100%）。

输出：每个基因在第二个 GWAS 中的效应与 p 值，并给出跨队列一致性统计。
"""
from __future__ import annotations
import argparse, pathlib, sys, time

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from io_gtex import load_tissue, split_variant, lead_per_gene
from liftover import Chain

DATA = pathlib.Path(__file__).resolve().parents[1] / "data"
TAR = DATA / "GTEx_Analysis_v8_eQTL.tar"
GWAS2 = DATA / "GCST006465_GRCh37.tsv.gz"
CHAIN = DATA / "hg38ToHg19.over.chain.gz"

USECOLS = ["chromosome", "base_pair_location", "effect_allele", "other_allele",
           "beta", "standard_error", "p_value", "variant_id"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", default="Uterus")
    ap.add_argument("--out", default="results/screen_uterus_gcst006465.tsv")
    a = ap.parse_args()

    t0 = time.time()
    print("[1/3] 读取子宫 eQTL lead 并 liftover 到 hg19 …")
    pairs = split_variant(load_tissue(TAR, a.tissue))
    leads = lead_per_gene(pairs[pairs["pos"].notna()])
    ch = Chain(CHAIN)
    lc, lp = [], []
    for r in leads.itertuples(index=False):
        v = ch.lift(r.chrom, int(r.pos))
        lc.append(v[0] if v else None)
        lp.append(v[1] if v else None)
    leads["hg19_chrom"] = lc
    leads["hg19_pos"] = lp
    ok = leads["hg19_pos"].notna().sum()
    print(f"      3,577 个 lead 中成功降坐标 {ok:,}（{ok/len(leads)*100:.1f}%）")
    leads = leads[leads["hg19_pos"].notna()].copy()

    want = set(zip(leads["hg19_chrom"].astype(str), leads["hg19_pos"].astype(int)))
    print(f"[2/3] 流式扫描 GCST006465（9.52M 行），只保留 {len(want):,} 个目标位置 …")
    hits = []
    nrows = 0
    for chunk in pd.read_csv(GWAS2, sep="\t", usecols=USECOLS,
                             chunksize=1_000_000, dtype={"chromosome": str}):
        nrows += len(chunk)
        key = list(zip(chunk["chromosome"].str.replace("chr", "", regex=False),
                       chunk["base_pair_location"]))
        m = [k in want for k in key]
        if any(m):
            sub = chunk.loc[m].copy()
            sub["chr19"] = [k[0] for k, mm in zip(key, m) if mm]
            sub["pos19"] = [k[1] for k, mm in zip(key, m) if mm]
            hits.append(sub)
    g = pd.concat(hits) if hits else pd.DataFrame(columns=USECOLS)
    print(f"      扫描 {nrows:,} 行，命中 {len(g):,} 条")

    COMP = str.maketrans("ACGT", "TGCA")
    rows = []
    for r in leads.itertuples(index=False):
        sel = g[(g["chr19"] == str(r.hg19_chrom)) & (g["pos19"] == int(r.hg19_pos))]
        if sel.empty:
            rows.append({"gene_id": r.gene_id, "variant_id": r.variant_id,
                         "eQTL_slope": r.slope, "in_gwas2": False})
            continue
        got = None
        for s in sel.itertuples(index=False):
            ea, oa = str(s.effect_allele).upper(), str(s.other_allele).upper()
            ref, alt = str(r.ref).upper(), str(r.alt).upper()
            if (ea, oa) == (alt, ref):
                got = (s.beta, s.standard_error, s.p_value, s.variant_id)
            elif (ea, oa) == (ref, alt):
                got = (-s.beta, s.standard_error, s.p_value, s.variant_id)
            elif (ea, oa) == (alt.translate(COMP), ref.translate(COMP)):
                got = (s.beta, s.standard_error, s.p_value, s.variant_id)
            elif (ea, oa) == (ref.translate(COMP), alt.translate(COMP)):
                got = (-s.beta, s.standard_error, s.p_value, s.variant_id)
            if got:
                break
        if got is None:
            rows.append({"gene_id": r.gene_id, "variant_id": r.variant_id,
                         "eQTL_slope": r.slope, "in_gwas2": True,
                         "allele_mismatch": True})
        else:
            rows.append({"gene_id": r.gene_id, "variant_id": r.variant_id,
                         "eQTL_slope": r.slope, "eQTL_p": r.pval_nominal,
                         "in_gwas2": True, "rsid": got[3],
                         "gwas2_beta": got[0], "gwas2_se": got[1], "gwas2_p": got[2]})

    df = pd.DataFrame(rows)
    out = pathlib.Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)

    print(f"\n[3/3] 写表 {out}")
    d = df[df.get("gwas2_p").notna()] if "gwas2_p" in df else pd.DataFrame()
    print(f"      基因数 {len(df):,}；能对齐到 GCST006465 的 {len(d):,}")
    if len(d):
        for thr in (0.05, 1e-3):
            s = d[d.gwas2_p < thr]
            if not len(s):
                continue
            c = int(((s.eQTL_slope * s.gwas2_beta) > 0).sum())
            p = stats.binomtest(c, len(s), 0.5).pvalue
            print(f"      GWAS2 p<{thr:g}: n={len(s):4d}  方向一致 {c}/{len(s)} "
                  f"= {c/len(s)*100:5.1f}%  binom p={p:.3g}")
    print(f"      用时 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
