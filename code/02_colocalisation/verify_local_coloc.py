#!/usr/bin/env python3
"""独立复核：用本地完整 eQTL 文件重算三个关键命中，与远端 tabix 结果对照。"""
import sys, pathlib, gzip
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from p0_coloc_omara import load_gwas, lookup
from coloc_abf import coloc_abf, log_abf

DATA = pathlib.Path(__file__).resolve().parent.parent / "data"
GWAS = DATA / "GCST90454186.h.tsv.gz"
PAD = 600_000

CASES = [
    ("QTD000266", "Liver",  "17", 48_127_474, "ENSG00000141293", 0.931),
    ("QTD000281", "Muscle", "17", 48_127_474, "ENSG00000278765", 0.898),
    ("QTD000281", "Muscle", "15", 40_086_266, "ENSG00000275636", 0.962),
]

store = load_gwas(GWAS)
print("GWAS 载入完成", flush=True)

for qtd, tissue, chrom, lead, gene, expect in CASES:
    lo, hi = max(1, lead - PAD), lead + PAD
    rows = []
    with gzip.open(DATA / f"{qtd}.all.tsv.gz", "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ic, ip, ig = hdr.index("chromosome"), hdr.index("position"), hdr.index("gene_id")
        ib, ise = hdr.index("beta"), hdr.index("se")
        imaf = hdr.index("maf")
        iref, ialt = hdr.index("ref"), hdr.index("alt")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[ic] != chrom:
                continue
            pos = int(p[ip])
            if pos < lo:
                continue
            if pos > hi:
                break            # 文件按位置排序
            if p[ig] != gene:
                continue
            try:
                b, se, maf = float(p[ib]), float(p[ise]), float(p[imaf])
            except ValueError:
                continue
            if not (np.isfinite(b) and np.isfinite(se)) or se <= 0:
                continue
            if not (0.01 <= maf <= 0.99):
                continue
            rows.append((pos, p[iref].upper(), p[ialt].upper(), b, se))

    d = []
    for pos, r, a, b_e, se_e in rows:
        hit = lookup(store, chrom, pos, r, a)
        if hit is None:
            continue
        b_g, se_g, _ = hit
        if np.isfinite(b_g) and np.isfinite(se_g) and se_g > 0:
            d.append((b_e, se_e, b_g, se_g))

    if len(d) < 10:
        print(f"{tissue} {gene}: 匹配不足 ({len(d)})", flush=True)
        continue
    arr = np.array(d)
    res = coloc_abf(log_abf(arr[:, 0], arr[:, 1], 1.0), log_abf(arr[:, 2], arr[:, 3], 1.0))
    h4 = res["PP.H4"]
    ok = "一致" if abs(h4 - expect) < 0.005 else "**不一致**"
    print(f"{tissue:7s} {gene:18s} n={len(d):4d}  H4={h4:.3f} (远端 {expect:.3f})  {ok}", flush=True)
