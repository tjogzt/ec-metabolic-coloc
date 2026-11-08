"""
p0_coloc_omara.py — P0-1 修复：把遗传层的结局换成 Methods 声称的 GCST006465（O'Mara 2018）。

原分析误用 FinnGen R12（最小 p=2.3e-7，全基因组显著位点 0），却按 GWAS z=10–25 做功效模拟。
本脚本：
  1. 从 GCST006465（GRCh37，313 个 p<5e-8）重新聚类位点；
  2. 位点 lead 抬到 GRCh38 以定位 eQTL 区域；
  3. 区域内 eQTL 变异（GRCh38）逐条降到 GRCh37，与 GWAS 按位置+等位基因配对上；
  4. 跑与原文完全相同的 ABF coloc。
输出列与 coloc_loci_*.tsv 保持一致，便于直接替换。
"""
from __future__ import annotations
import argparse, pathlib, sys, time
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
ROOT = HERE.parent
DATA = ROOT / "data"
RES = HERE / "results"

from liftover import Chain
from coloc_abf import coloc_abf, log_abf
from coloc_eqtlcat import qtl_region
from coloc_loci import tissue_egenes_in_loci, TISSUE_QTD

GWAS = DATA / "GCST006465_GRCh37.tsv.gz"
CH_UP = DATA / "hg19ToHg38.over.chain.gz"      # 位点 lead：hg19 -> hg38
CH_DN = DATA / "hg38ToHg19.over.chain.gz"      # eQTL 变异：hg38 -> hg19


def load_gwas(path=GWAS) -> dict:
    """把 GWAS 按染色体装进内存，供按位置二分查询。"""
    store = {}
    for ch in pd.read_csv(path, sep="\t", comment="#", low_memory=False,
                          usecols=["variant_id", "base_pair_location", "chromosome",
                                   "effect_allele", "other_allele", "beta",
                                   "standard_error", "p_value"],
                          chunksize=2_000_000):
        ch = ch.rename(columns={"base_pair_location": "pos", "chromosome": "chrom",
                                "effect_allele": "alt", "other_allele": "ref",
                                "standard_error": "se", "p_value": "pval"})
        ch["chrom"] = ch["chrom"].astype(str).str.replace("chr", "", regex=False)
        ch = ch[ch["ref"].str.len().eq(1) & ch["alt"].str.len().eq(1)]
        ch = ch[np.isfinite(ch["beta"]) & np.isfinite(ch["se"]) & (ch["se"] > 0)]
        for c, sub in ch.groupby("chrom"):
            sub = sub.sort_values("pos")
            rec = store.get(c)
            if rec is None:
                store[c] = sub[["pos", "ref", "alt", "beta", "se", "pval"]].reset_index(drop=True)
            else:
                store[c] = (pd.concat([rec, sub[["pos", "ref", "alt", "beta", "se", "pval"]]])
                            .sort_values("pos").reset_index(drop=True))
    for c in store:
        store[c] = store[c].drop_duplicates(subset=["pos", "ref", "alt"], keep="first")
    return store


def peaks(store: dict, p_thresh: float, gap: int) -> pd.DataFrame:
    rows = []
    for chrom, d in store.items():
        d = d[d["pval"] < p_thresh].sort_values("pos")
        if d.empty:
            continue
        cur = None
        for x in d.itertuples(index=False):
            if cur is None or x.pos - cur["end"] > gap:
                if cur:
                    rows.append(cur)
                cur = {"chrom": chrom, "start": x.pos, "end": x.pos, "lead": x.pos, "p": x.pval}
            else:
                cur["end"] = x.pos
                if x.pval < cur["p"]:
                    cur["p"], cur["lead"] = x.pval, x.pos
        if cur:
            rows.append(cur)
    return pd.DataFrame(rows).sort_values("p").reset_index(drop=True)


def lookup(store, chrom, pos, ref, alt):
    """按 hg19 位置 + 等位基因取 GWAS 行，自动处理链翻转。"""
    d = store.get(str(chrom))
    if d is None:
        return None
    i = int(np.searchsorted(d["pos"].values, pos, side="left"))
    best = None
    while i < len(d) and d["pos"].iat[i] == pos:
        R, A = str(d["ref"].iat[i]).upper(), str(d["alt"].iat[i]).upper()
        if (R, A) == (ref, alt):
            best = (d["beta"].iat[i], d["se"].iat[i], d["pval"].iat[i])
        elif (R, A) == (alt, ref):
            best = (-d["beta"].iat[i], d["se"].iat[i], d["pval"].iat[i])
        if best is not None:
            break
        i += 1
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", required=True)
    ap.add_argument("--pad", type=int, default=600_000)
    ap.add_argument("--min_snps", type=int, default=10)
    ap.add_argument("--p_thresh", type=float, default=5e-8)
    ap.add_argument("--gap", type=int, default=1_000_000)
    a = ap.parse_args()

    t0 = time.time()
    print("载入 GWAS…", flush=True)
    store = load_gwas()
    print(f"  染色体数 {len(store)}，变异合计 {sum(len(v) for v in store.values()):,}"
          f"（{time.time()-t0:.0f}s）", flush=True)

    loc19 = peaks(store, a.p_thresh, a.gap)
    print(f"O'Mara 位点（p<{a.p_thresh:g}）{len(loc19)} 个", flush=True)

    up = Chain(CH_UP)
    loc = []
    for r in loc19.itertuples(index=False):
        m = up.lift(r.chrom, int(r.lead))
        if m is None:
            continue
        loc.append({"chrom": str(m[0]), "start": int(m[1]), "end": int(m[1]),
                    "lead": int(m[1]), "p": r.p, "lead19": int(r.lead)})
    loc = pd.DataFrame(loc)
    print(f"可抬升到 GRCh38 的位点 {len(loc)} 个", flush=True)

    eg = tissue_egenes_in_loci(a.tissue, loc, pad=a.pad)
    print("各位点候选 eGene 数：" + str({i: len(v) for i, v in eg.items()}), flush=True)

    qtd = TISSUE_QTD[a.tissue]
    src = DATA / f"{qtd}.all.tsv.gz"
    if not src.exists():
        raise SystemExit(f"缺少本地 eQTL 文件 {src}")
    dn = Chain(CH_DN)

    recs = []
    for i, r in loc.iterrows():
        genes = eg.get(i, [])
        if not genes:
            continue
        lo, hi = max(1, r.lead - a.pad), r.lead + a.pad
        q = qtl_region(r.chrom, lo, hi, src)
        if q.empty:
            continue
        # 区域内 eQTL 变异一次性降坐标
        q = q.dropna(subset=["beta", "se"])
        q = q[(q.maf >= 0.01) & (q.maf <= 0.99) & (q.se > 0)]
        if q.empty:
            continue
        pos19 = [dn.lift(r.chrom, int(p)) for p in q.position.values]
        q = q.assign(pos19=[x[1] if x else None for x in pos19])
        q = q.dropna(subset=["pos19"])
        q["pos19"] = q["pos19"].astype(int)
        for g in genes:
            sub = q[q.gene_id == g]
            if len(sub) < a.min_snps:
                continue
            rows = []
            for x in sub.itertuples(index=False):
                hit = lookup(store, r.chrom, int(x.pos19),
                             str(x.ref).upper(), str(x.alt).upper())
                if hit is None:
                    continue
                b, se, p = hit
                if np.isfinite(b) and np.isfinite(se) and se > 0:
                    rows.append((x.beta, x.se, b, se, p))
            if len(rows) < a.min_snps:
                continue
            d = pd.DataFrame(rows, columns=["b_e", "se_e", "b_g", "se_g", "p_g"])
            res = coloc_abf(log_abf(d.b_e.values, d.se_e.values, 1.0),
                            log_abf(d.b_g.values, d.se_g.values, 1.0))
            res = {k: v for k, v in res.items() if k not in ("n_snps", "top_snp_share")}
            recs.append(dict(**res, n_snps=len(d),
                             top_snp_share=np.nan, tissue=a.tissue,
                             locus=f"{r.chrom}:{r.lead19}", gene_id=g,
                             best_gwas_p=float(d.p_g.min())))
    out = pd.DataFrame(recs)
    RES.mkdir(exist_ok=True)
    if out.empty:
        print("无结果", flush=True)
        return
    out = out[["PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "n_snps",
               "top_snp_share", "tissue", "locus", "gene_id", "best_gwas_p"]]
    f = RES / f"coloc_omara_{a.tissue}.tsv"
    out.to_csv(f, sep="\t", index=False)
    print(f"\n{len(out)} 次检验 → {f}", flush=True)
    print(f"  PP.H4>0.8: {(out['PP.H4']>0.8).sum()} | >0.5: {(out['PP.H4']>0.5).sum()}"
          f" | 最大 {out['PP.H4'].max():.3f}", flush=True)
    print(f"  PP.H3>0.8: {(out['PP.H3']>0.8).sum()}", flush=True)
    print(f"总耗时 {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
