#!/usr/bin/env python3
"""B1: GTEx v11 半独立复制检查 — 17q21.32（第2版，修正列名与取数路径）。"""
import pathlib, tarfile
import pandas as pd
import numpy as np

ROOT = pathlib.Path("/Users/taozhu/my researches/gongbing/06_EC_IR")
TAR = ROOT / "data" / "GTEx_Analysis_v11_eQTL_SuSiE.tar"
TMP = pathlib.Path("/tmp/gtex_v11")
TMP.mkdir(exist_ok=True)

TISSUES = ["Liver", "Muscle_Skeletal", "Adipose_Subcutaneous",
           "Adipose_Visceral_Omentum", "Uterus", "Pancreas", "Spleen",
           "Breast_Mammary_Tissue", "Ovary", "Vagina", "Whole_Blood",
           "Thyroid", "Kidney_Cortex"]

def get_parquet(tissue):
    p = TMP / f"GTEx_Analysis_v11_eQTL/{tissue}.v11.eQTLs.SuSiE_summary.parquet"
    if not p.exists():
        with tarfile.open(TAR) as tf:
            member = f"GTEx_Analysis_v11_eQTL/{tissue}.v11.eQTLs.SuSiE_summary.parquet"
            tf.extract(member, TMP)
    return pd.read_parquet(p)

# 1) 主分析基因：17q21.32 各行
ra = pd.read_csv(ROOT / "pipeline/results/coloc_ramach_all.tsv", sep="\t")
l17 = ra[(ra.locus == "17:48127474")]
print("17q21.32 主分析行 (PP.H4):")
print(l17[["tissue", "gene_id", "PP.H3", "PP.H4", "n_snps"]].to_string())

# 2) v8 SKAP1 liver lead 变异（从窗口表取 p 最小行）
win = pd.read_csv(ROOT / "data/_tmp_ramach_QTD000266_loci.tsv", sep="\t")
print("\n_tmp_ramach_QTD000266_loci 列:", list(win.columns))
sk = win[win.gene_id == "ENSG00000141293"].copy()
if len(sk):
    lead = sk.loc[sk.pvalue.idxmin()]
    print("v8 SKAP1 liver lead:", lead[["position", "ref", "alt", "pvalue"]].to_dict())

# 3) 风险 lead 变异（GCST90454186 17q21.32 位点；文件已归档至外置盘）
GWAS_FILE = pathlib.Path("/Volumes/tjogzt4T/archive_gongbing_20261117/06_EC_IR/data/GCST90454186.h.tsv.gz")
import gzip
target = 48127474
hit = None
with gzip.open(GWAS_FILE, "rt") as f:
    hdr = f.readline().strip().split("\t")
    # 列: chromosome base_pair_location effect_allele other_allele beta se eaf p_value rsid ...
    for line in f:
        c = line.strip().split("\t")
        try:
            chrom = c[0]
            pos = int(c[1])
        except (ValueError, IndexError):
            continue
        if chrom == "17" and abs(pos - target) < 2000:
            hit = (hdr, c)
            break
if hit:
    hdr, c = hit
    d = dict(zip(hdr, c))
    print("风险 lead 行:", {k: d.get(k) for k in ["#CHROM", "POS", "ID", "REF", "ALT", "P"] if k in d})

# 4) v11 各组织 SKAP1 / RP5 / CBX1 状态
risk_pos = target
out_rows = []
for t in TISSUES:
    try:
        df = get_parquet(t)
    except Exception as e:
        print(t, "ERROR", e)
        continue
    for name in ["SKAP1"]:
        s = df[df.gene_name == name]
        if len(s) == 0:
            out_rows.append((t, name, 0, None, None))
            continue
        css = s.cs_id.dropna().unique()
        # 找位置距风险 lead 最近的 CS 变异
        s2 = s.copy()
        s2["pos"] = s2.variant_id.str.split("_").str[1].astype(int)
        near = s2.assign(d=lambda x: (x.pos - risk_pos).abs()).nsmallest(3, "d")
        top = s.nlargest(2, "pip")
        out_rows.append((t, name, len(css), "|".join(map(str, css)),
                         "; ".join(f"{r.variant_id.split('_')[1]}:pip={r.pip:.3f},afc={r.afc:.3f}" for r in top.itertuples())))
        print(f"{t} SKAP1: CSs={sorted(css)} n={len(s)}")
        print("  近风险lead的CS变异:", [f"{r.variant_id.split('_')[1]}:pip={r.pip:.3f}" for r in near.itertuples()])
res = pd.DataFrame(out_rows, columns=["tissue", "gene", "n_cs", "cs_ids", "top"])
res.to_csv(ROOT / "pipeline/results/v11_suSie_check.tsv", sep="\t", index=False)
print("\n落盘 pipeline/results/v11_suSie_check.tsv")
