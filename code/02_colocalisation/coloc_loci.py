"""
coloc_loci.py — 位点导向的组织面板 coloc（不下载全库，远端按位点取数）

思路：EC 的强信号位点很少（p<1e-5 聚成若干位点），每个位点周围真正需要检验的
只有该组织在该位点有 eQTL 的基因。于是每个位点只发一次远端 tabix 请求（约 1.2 Mb），
一次拿回该位点内所有基因的全变异统计，再逐基因做 coloc。
每组织约 20 次请求、8–15 分钟，而下载全库要 3–7 小时。

用法：python coloc_loci.py --tissue Adipose_Visceral_Omentum
"""
from __future__ import annotations
import argparse, gzip, json, pathlib, subprocess, sys, tarfile, time

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from io_finngen import FinnGen
from coloc_abf import coloc_abf, log_abf
from coloc_eqtlcat import qtl_region

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RES = ROOT / "pipeline" / "results"
GWAS = DATA / "finngen_R12_C3_CORPUS_UTERI_EXALLC.gz"
TAR = DATA / "GTEx_Analysis_v8_eQTL.tar"

# 组织 → eQTL Catalogue 数据集
TISSUE_QTD = {
    "Uterus": "QTD000346",
    "Adipose_Visceral_Omentum": "QTD000121",
    "Liver": "QTD000266",
    "Muscle_Skeletal": "QTD000281",
    "Pancreas": "QTD000296",
    "Ovary": "QTD000291",
    "Whole_Blood": "QTD000356",
}
BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015"


def loci(p_thresh=1e-5, gap=1_000_000) -> pd.DataFrame:
    """FinnGen EC 中 p<p_thresh 的变异按距离聚成位点。"""
    fg = FinnGen(GWAS)
    rows = []
    for chrom in [str(c) for c in range(1, 23)]:
        d = fg.region(chrom, 1, 250_000_000)
        if d.empty:
            continue
        d = d[d["pval"] < p_thresh]
        if d.empty:
            continue
        d = d.sort_values("pos")
        cur = None
        for x in d.itertuples(index=False):
            if cur is None or x.pos - cur["end"] > gap:
                if cur: rows.append(cur)
                cur = {"chrom": chrom, "start": x.pos, "end": x.pos,
                       "lead": x.pos, "p": x.pval}
            else:
                cur["end"] = x.pos
                if x.pval < cur["p"]:
                    cur["p"], cur["lead"] = x.pval, x.pos
        if cur: rows.append(cur)
    return pd.DataFrame(rows)


def tissue_egenes_in_loci(tissue: str, loc: pd.DataFrame, pad=1_000_000) -> dict[str, list[str]]:
    """从 GTEx 显著对里找出落在各位点附近且有 eQTL 的基因。"""
    out = {i: [] for i in range(len(loc))}
    with tarfile.open(TAR, "r") as tf:
        path = next((m for m in tf.getmembers()
                     if "signif_variant_gene_pairs" in m.name
                     and pathlib.Path(m.name).name.startswith(tissue + ".")), None)
        if path is None:
            raise SystemExit(f"tar 中找不到 {tissue}")
        f = gzip.GzipFile(fileobj=tf.extractfile(path))
        for line in f:
            p = line.decode().split("\t")
            if p[0] == "variant_id":
                continue
            pos = int(p[0].split("_")[1])
            for i, r in loc.iterrows():
                if r["chrom"] == p[0].split("_")[0].replace("chr", "") and \
                   r["start"] - pad <= pos <= r["end"] + pad:
                    out[i].append(p[1].split(".")[0])
    return {k: sorted(set(v)) for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", required=True)
    ap.add_argument("--pad", type=int, default=600_000)
    ap.add_argument("--min_snps", type=int, default=10)
    ap.add_argument("--p_thresh", type=float, default=1e-6)
    a = ap.parse_args()

    qtd = TISSUE_QTD[a.tissue]
    src = DATA / f"{qtd}.all.tsv.gz"
    url = f"{BASE}/{qtd}/{qtd}.all.tsv.gz"
    # 必须与远端 Content-Length 比对：残缺的本地文件会让 tabix 静默返回空
    head = subprocess.run(f"curl -sI -m 60 {url} | grep -i '^content-length'",
                          shell=True, capture_output=True, text=True).stdout
    exp = int(head.split(":")[-1].strip()) if head.strip() else 0
    complete = src.exists() and exp and src.stat().st_size == exp
    if complete:
        remote = False
        print(f"本地文件完整：{src.stat().st_size/1e9:.2f} GB")
    else:
        src = url
        remote = True
        if src and (DATA / f"{qtd}.all.tsv.gz").exists():
            print(f"本地文件不完整（{(DATA / f'{qtd}.all.tsv.gz').stat().st_size/1e9:.2f} / {exp/1e9:.2f} GB），改用远端")
    print(f"数据源：{'远端 ' if remote else '本地 '}{src}")

    loc = loci(p_thresh=a.p_thresh)
    print(f"EC 位点 {len(loc)} 个：" +
          ", ".join(f"{r.chrom}:{r.lead}" for r in loc.itertuples(index=False)))
    eg = tissue_egenes_in_loci(a.tissue, loc)
    print("各位点候选 eGene 数：" + str({i: len(v) for i, v in eg.items()}))

    fg = FinnGen(GWAS)
    recs = []
    for i, r in loc.iterrows():
        genes = eg[i]
        if not genes:
            continue
        t0 = time.time()
        lo, hi = max(1, r["lead"] - a.pad), r["lead"] + a.pad
        q = qtl_region(r["chrom"], lo, hi, src)
        gw = fg.region(r["chrom"], lo, hi)
        if q.empty or gw.empty:
            continue
        k2 = gw.pos.astype(str) + "_" + gw.ref + "_" + gw.alt
        k3 = gw.pos.astype(str) + "_" + gw.alt + "_" + gw.ref
        lut = dict(zip(k2, zip(gw.beta, gw.sebeta, gw.pval)))
        lut.update(dict(zip(k3, zip(-gw.beta, gw.sebeta, gw.pval))))
        for g in genes:
            sub = q[q.gene_id == g].dropna(subset=["beta", "se"])
            sub = sub[(sub.maf >= 0.01) & (sub.maf <= 0.99) & (sub.se > 0)]
            if len(sub) < a.min_snps:
                continue
            rows = []
            for x in sub.itertuples(index=False):
                k = f"{int(x.position)}_{x.ref}_{x.alt}"
                if k in lut:
                    b, se, p = lut[k]
                    if np.isfinite(b) and np.isfinite(se) and se > 0:
                        rows.append((x.beta, x.se, b, se, p))
            if len(rows) < a.min_snps:
                continue
            d = pd.DataFrame(rows, columns=["b_e", "se_e", "b_g", "se_g", "p_g"])
            res = coloc_abf(log_abf(d.b_e.values, d.se_e.values, 1.0),
                            log_abf(d.b_g.values, d.se_g.values, 1.0))
            res.update({"tissue": a.tissue, "locus": f"{r.chrom}:{r.lead}",
                        "gene_id": g, "n_snps": len(d),
                        "best_gwas_p": float(d.p_g.min())})
            recs.append(res)
        print(f"  位点 {r.chrom}:{r.lead}  {len(genes)} 基因  用时 {time.time()-t0:.0f}s")

    out = RES / f"coloc_loci_{a.tissue}.tsv"
    if not recs:
        print(f"!! {a.tissue} 没有任何可用基因-位点配对（远端查询可能受限），未写表")
        return
    df = pd.DataFrame(recs).sort_values("PP.H4", ascending=False)
    df.to_csv(out, sep="\t", index=False)
    print(f"\n写出 {out}：{len(df)} 条")
    if len(df):
        for thr in (0.5, 0.8, 0.9):
            print(f"  PP.H4 > {thr}: {int((df['PP.H4'] > thr).sum())}")
        print(f"  PP.H3 > 0.8: {int((df['PP.H3'] > 0.8).sum())}")
        print(df.head(8)[["gene_id", "locus", "n_snps", "PP.H3", "PP.H4"]].to_string(index=False))


if __name__ == "__main__":
    main()
