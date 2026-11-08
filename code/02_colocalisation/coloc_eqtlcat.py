"""
coloc_eqtlcat.py — 用 eQTL Catalogue 的区域全变异汇总统计做标准 coloc

数据：QTS000015 / QTD000346 = GTEx v8 uterus（n=129，quant_method=ge），
      文件 QTD000346.all.tsv.gz（2.9 GB，bgzip+tabix），含每个基因 cis 窗口内**全部**
      被检验过的变异（不是只有显著对）——这正是标准 coloc.abf 需要的输入。

列：molecular_trait_id chromosome position ref alt variant ma_samples maf pvalue
    beta se type ac an r2 molecular_trait_object_id gene_id median_tpm rsid

与 FinnGen（GRCh38，本地 tabix）同坐标，不需要 liftOver。

用法：
  python coloc_eqtlcat.py --min_snps 10 --out results/coloc_full_uterus.tsv
"""
from __future__ import annotations
import argparse, pathlib, subprocess, sys, time

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from io_finngen import FinnGen
from coloc_abf import coloc_abf, log_abf

DATA = pathlib.Path(__file__).resolve().parents[1] / "data"
EQTL = DATA / "QTD000346.all.tsv.gz"
REMOTE = ("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
          "QTS000015/QTD000346/QTD000346.all.tsv.gz")
GWAS = DATA / "finngen_R12_C3_CORPUS_UTERI_EXALLC.gz"

ECOLS = ["gene_id", "chromosome", "position", "ref", "alt", "variant", "maf",
         "pvalue", "beta", "se", "median_tpm", "rsid"]
KEEP = [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 17, 18]   # 列序号（0-based）


def qtl_region(chrom: str, start: int, end: int, src=None) -> pd.DataFrame:
    out = subprocess.run(
        ["tabix", str(src or EQTL), f"{chrom}:{start}-{end}"],
        capture_output=True, text=True)
    if out.returncode != 0 or not out.stdout.strip():
        return pd.DataFrame(columns=ECOLS)
    rows = [l.split("\t") for l in out.stdout.strip().split("\n")]
    df = pd.DataFrame(rows).iloc[:, KEEP]
    df.columns = ECOLS
    for c in ("position", "maf", "pvalue", "beta", "se", "median_tpm"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def gene_bounds(gene_id: str) -> tuple[str, int, int] | None:
    """用该基因在 eQTL 文件中的行范围反推其 cis 窗口（一次查询拿到）。"""
    out = subprocess.run(["tabix", str(EQTL), "--list-chromosomes"],
                         capture_output=True, text=True)
    return None  # 见 build_gene_index


def build_gene_index(limit: int = 0) -> pd.DataFrame:
    """流式扫一遍文件，只为建立 gene_id -> (chrom, min_pos, max_pos)。"""
    gene = pd.read_csv(EQTL, sep="\t", usecols=["molecular_trait_id", "chromosome",
                                                "position"], dtype={"chromosome": str},
                       chunksize=2_000_000)
    acc = {}
    for ch in gene:
        g = ch.groupby("molecular_trait_id")
        for gid, sub in g:
            c = sub["chromosome"].iloc[0]
            lo, hi = sub["position"].min(), sub["position"].max()
            if gid in acc:
                c0, lo0, hi0 = acc[gid]
                acc[gid] = (c, min(lo, lo0), max(hi, hi0))
            else:
                acc[gid] = (c, lo, hi)
    df = pd.DataFrame([(g, v[0], v[1], v[2]) for g, v in acc.items()],
                      columns=["gene_id", "chrom", "lo", "hi"])
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min_snps", type=int, default=10)
    ap.add_argument("--out", default="results/coloc_full_uterus.tsv")
    ap.add_argument("--genes", default="", help="限定基因列表文件（每行一个 ENSG）")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--remote", action="store_true", help="直接读 EBI 远端（慢，用于试跑）")
    ap.add_argument("--eqtl", default="", help="指定组织全变异文件（默认子宫）")
    ap.add_argument("--index", default="", help="指定基因索引文件（默认子宫索引）")
    ap.add_argument("--gwas", default="", help="指定结局 GWAS（默认 EC；用于阳性对照）")
    a = ap.parse_args()

    t0 = time.time()
    print("[1/3] 建立基因窗口索引（流式扫描 2.9 GB 文件）…")
    gi_path = pathlib.Path(a.index or "results/gene_index.tsv")
    if gi_path.exists():
        gi = pd.read_csv(gi_path, sep="\t")
        print(f"      复用已有索引：{len(gi):,} 个基因")
    else:
        gi = build_gene_index()
        gi_path.parent.mkdir(parents=True, exist_ok=True)
        gi.to_csv(gi_path, sep="\t", index=False)
        print(f"      建成索引：{len(gi):,} 个基因（{time.time()-t0:.0f}s）")

    if a.genes:
        keep = {l.strip() for l in open(a.genes) if l.strip()}
        gi = gi[gi["gene_id"].isin(keep)]
    if a.limit:
        gi = gi.head(a.limit)

    src = REMOTE if a.remote else pathlib.Path(a.eqtl or EQTL)
    print(f"[2/3] 待分析基因 {len(gi):,}（数据源 {pathlib.Path(src).name}；"
          f"结局 {pathlib.Path(a.gwas or GWAS).name}）")
    fg = FinnGen(a.gwas or GWAS)
    recs = []
    for i, r in enumerate(gi.itertuples(index=False), 1):
        e = qtl_region(r.chrom, int(r.lo), int(r.hi), src)
        # 关键：窗口查询会带回窗口内所有基因的行，必须只保留目标基因
        e = e[e["gene_id"] == r.gene_id]
        e = e.dropna(subset=["position", "beta", "se", "pvalue"])
        e = e[(e["maf"] >= 0.01) & (e["maf"] <= 0.99) & (e["se"] > 0)]
        if len(e) < a.min_snps:
            continue
        g = fg.region(r.chrom, int(e["position"].min()) - 1, int(e["position"].max()) + 1)
        if g.empty:
            continue
        g = g.copy()
        g["key"] = g["pos"].astype(str) + "_" + g["ref"] + "_" + g["alt"]
        g["key_sw"] = g["pos"].astype(str) + "_" + g["alt"] + "_" + g["ref"]
        lut = {}
        for x in g.itertuples(index=False):
            lut[x.key] = (x.beta, x.sebeta, x.pval)
            lut[x.key_sw] = (-x.beta, x.sebeta, x.pval)

        rows = []
        for x in e.itertuples(index=False):
            k = f"{int(x.position)}_{x.ref}_{x.alt}"
            if k in lut:
                b, se, p = lut[k]
                if np.isfinite(b) and np.isfinite(se) and se > 0:
                    rows.append((x.beta, x.se, b, se, p))
        if len(rows) < a.min_snps:
            continue
        d = pd.DataFrame(rows, columns=["b_e", "se_e", "b_g", "se_g", "p_g"])
        l1 = log_abf(d["b_e"].values, d["se_e"].values, 1.0)
        l2 = log_abf(d["b_g"].values, d["se_g"].values, 1.0)
        res = coloc_abf(l1, l2)
        res.update({"gene_id": r.gene_id, "n_snps": len(d),
                    "best_gwas_p": float(d["p_g"].min()),
                    "top_sign_concordant": bool(d.iloc[int(np.argmin(d["p_g"]))]["b_e"] *
                                                d.iloc[int(np.argmin(d["p_g"]))]["b_g"] > 0)})
        recs.append(res)
        if i % 200 == 0:
            print(f"      … {i}/{len(gi)}  （{time.time()-t0:.0f}s）")

    res = pd.DataFrame(recs).sort_values("PP.H4", ascending=False).reset_index(drop=True)
    out = pathlib.Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out, sep="\t", index=False)
    print(f"\n[3/3] 写表 {out}：{len(res):,} 个基因")
    for thr in (0.5, 0.8, 0.9):
        print(f"      PP.H4 > {thr}: {int((res['PP.H4'] > thr).sum()):,}")
    print(f"      PP.H3 > 0.8: {int((res['PP.H3'] > 0.8).sum()):,}")
    print(f"      n_snps 中位 {res['n_snps'].median():.0f}")
    print(f"      用时 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
