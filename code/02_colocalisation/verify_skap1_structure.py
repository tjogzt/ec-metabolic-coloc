#!/usr/bin/env python3
"""用新 GWAS + 完整肝 eQTL 重算 SKAP1 内部结构统计量。"""
import sys, pathlib, gzip
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from p0_coloc_omara import load_gwas, lookup
from scipy.stats import spearmanr

DATA = pathlib.Path(__file__).resolve().parent.parent / "data"
GWAS = DATA / "GCST90454186.h.tsv.gz"
CHROM, LEAD, GENE = "17", 48_127_474, "ENSG00000141293"
PAD = 600_000

store = load_gwas(GWAS)
print("GWAS 载入", flush=True)

lo, hi = LEAD - PAD, LEAD + PAD
recs = []
with gzip.open(DATA / "QTD000266.all.tsv.gz", "rt") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    ic, ip, ig = hdr.index("chromosome"), hdr.index("position"), hdr.index("gene_id")
    ib, ise, ipv = hdr.index("beta"), hdr.index("se"), hdr.index("pvalue")
    imaf, iref, ialt = hdr.index("maf"), hdr.index("ref"), hdr.index("alt")
    for line in f:
        p = line.rstrip("\n").split("\t")
        if p[ic] != CHROM:
            continue
        pos = int(p[ip])
        if pos < lo:
            continue
        if pos > hi:
            break
        if p[ig] != GENE:
            continue
        try:
            b, se = float(p[ib]), float(p[ise])
        except ValueError:
            continue
        if not (np.isfinite(b) and np.isfinite(se)) or se <= 0:
            continue
        recs.append((pos, p[iref].upper(), p[ialt].upper(), b, se, float(p[ipv])))

rows = []
for pos, r, a, b_e, se_e, p_e in recs:
    hit = lookup(store, CHROM, pos, r, a)
    if hit is None:
        continue
    b_g, se_g, p_g = hit
    rows.append((pos, b_e, se_e, p_e, b_g, se_g, p_g))

d = pd.DataFrame(rows, columns=["pos", "b_e", "se_e", "p_e", "b_g", "se_g", "p_g"])
print(f"匹配变异 {len(d):,}", flush=True)

z_e = d.b_e / d.se_e
z_g = d.b_g / d.se_g
rho, prho = spearmanr(z_e, z_g)
print(f"eQTL↔risk z 相关 ρ={rho:.3f}  p={prho:.2e}", flush=True)

i_e = int(d.p_e.idxmin()); i_g = int(d.p_g.idxmin())
print(f"eQTL 峰 {d.pos[i_e]} (z={z_e[i_e]:.2f}, p={d.p_e[i_e]:.2e})", flush=True)
print(f"risk 峰 {d.pos[i_g]} (p={d.p_g[i_g]:.2e})", flush=True)
print(f"峰距 {abs(int(d.pos[i_e])-int(d.pos[i_g]))/1000:.1f} kb", flush=True)
print(f"lead eQTL: beta={d.b_e[i_e]:.4f} se={d.se_e[i_e]:.4f} z={z_e[i_e]:.2f} p={d.p_e[i_e]:.2e}", flush=True)
print(f"该变异 risk: beta={d.b_g[i_e]:.4f} p={d.p_g[i_e]:.2e}", flush=True)
