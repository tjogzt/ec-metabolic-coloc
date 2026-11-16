#!/usr/bin/env python3
"""R1-8：在当前（Ramachandran/GCST90454186）数据上重算 17q21.32 链审计三分组。
避免重读 2.9 GB eQTL 全文件：用本运行已缓存的位点文件 _tmp_ramach_QTD000266_loci.tsv。
输出：pipeline/skap1_strand_current.log
"""
import sys, pathlib, gzip
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from p0_coloc_omara import load_gwas

DATA = pathlib.Path(__file__).resolve().parent.parent / "data"
GWAS = DATA / "GCST90454186.h.tsv.gz"
GENE = "ENSG00000141293"          # SKAP1
LEAD = 48_127_474
PAD = 600_000
LO, HI = LEAD - PAD, LEAD + PAD

print("load GWAS ...", flush=True)
store = load_gwas(GWAS)
print("done", flush=True)

# 缓存文件列: gene_id chromosome position ref alt variant maf pvalue beta se median_tpm rsid
df = pd.read_csv(DATA / "_tmp_ramach_QTD000266_loci.tsv", sep="\t")
sub = df[(df.gene_id == GENE) & (df.position >= LO) & (df.position <= HI)].copy()
print(f"eQTL rows in window: {len(sub)}", flush=True)

g17 = store["17"]
g17 = g17[(g17.pos >= LO) & (g17.pos <= HI)]
idx = {}
for r in g17.itertuples(index=False):
    idx.setdefault(r.pos, {})[frozenset((r.ref.upper(), r.alt.upper()))] = (r.beta, r.se, r.pval)

rows = []
for r in sub.itertuples():
    ref, alt = r.ref.upper(), r.alt.upper()
    hit = idx.get(r.position, {}).get(frozenset((ref, alt)))
    if hit is None:
        continue
    b_g, se_g, p_g = hit
    if not (np.isfinite(b_g) and np.isfinite(se_g)) or se_g <= 0:
        continue
    rows.append((r.position, ref, alt, r.beta, r.se, b_g, se_g, r.maf))
d = pd.DataFrame(rows, columns=["pos", "ref", "alt", "b_e", "se_e", "b_g", "se_g", "maf"])
d["z_e"] = d.b_e / d.se_e
d["z_g"] = d.b_g / d.se_g

pal = d[["ref", "alt"]].apply(
    lambda s: sorted((str(s.ref).upper(), str(s.alt).upper())) in (["A", "T"], ["C", "G"]), axis=1)

def rho(x):
    r_, p_ = spearmanr(x.z_e, x.z_g)
    return r_, p_

n_all = len(d)
r_all, p_all = rho(d)
amb = d[pal]; unamb = d[~pal]
r_amb, p_amb = rho(amb)
r_una, p_una = rho(unamb)
amb_int = amb[(amb.maf >= 0.05) & (amb.maf <= 0.45)]
r_int, p_int = rho(amb_int)

with open(pathlib.Path(__file__).resolve().parent / "skap1_strand_current.log", "w") as f:
    L = [
        f"SKAP1 17q21.32 strand audit on CURRENT dataset (GCST90454186)",
        f"matched variants: {n_all}",
        f"all: rho={r_all:.3f} p={p_all:.2e}",
        f"strand-ambiguous (A/T, C/G): {len(amb)} ({100*len(amb)/n_all:.1f}%) rho={r_amb:.3f} p={p_amb:.2e}",
        f"unambiguous: {len(unamb)} rho={r_una:.3f} p={p_una:.2e}",
        f"ambiguous, intermediate MAF 0.05-0.45: {len(amb_int)} rho={r_int:.3f} p={p_int:.2e}",
    ]
    for s in L:
        print(s, flush=True)
        f.write(s + "\n")
print("written skap1_strand_current.log", flush=True)
