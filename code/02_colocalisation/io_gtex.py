"""
io_gtex.py — 从 GTEx v8 eQTL tar 中抽取指定组织的显著 variant-gene pairs

GTEx_Analysis_v8_eQTL.tar 内含每个组织的 `*.v8.signif_variant_gene_pairs.txt.gz`，
列为：gene_id variant_id tss_distance af ma_samples ma_count pval_nominal slope slope_se

注意：这只是“显著对”（每个基因通常只有 1--数个变异），用于筛选与定位；
真正的 coloc 需要该区域**全部**被检验过的变异，需另取 all-associations 数据。
"""
from __future__ import annotations
import tarfile
from pathlib import Path

import pandas as pd

COLS = ["gene_id", "variant_id", "tss_distance", "ma_samples",
        "ma_count", "maf", "pval_nominal", "slope", "slope_se"]


def list_tissues(tar_path: str | Path) -> list[str]:
    with tarfile.open(tar_path, "r") as tf:
        return sorted({Path(m.name).name.split(".")[0]
                       for m in tf.getmembers() if "signif_variant_gene_pairs" in m.name})


def load_tissue(tar_path: str | Path, tissue: str) -> pd.DataFrame:
    """返回该组织的显著 eQTL 表。tissue 形如 'Uterus'。"""
    with tarfile.open(tar_path, "r") as tf:
        member = None
        for m in tf.getmembers():
            if "signif_variant_gene_pairs" in m.name and Path(m.name).name.startswith(tissue + "."):
                member = m
                break
        if member is None:
            raise FileNotFoundError(f"{tissue} 不在 {tar_path} 中")
        f = tf.extractfile(member)
        df = pd.read_csv(f, sep="\t", compression="gzip", dtype=str)
    df.columns = [c.lower() for c in df.columns]
    for c in ("pval_nominal", "slope", "slope_se", "maf", "tss_distance"):
        if c in df:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def split_variant(df: pd.DataFrame) -> pd.DataFrame:
    """把 chr1_12345_A_G_b38 拆成 chrom/pos/ref/alt。"""
    parts = df["variant_id"].str.split("_", expand=True)
    df = df.copy()
    df["chrom"] = parts[0].str.replace("chr", "", regex=False)
    df["pos"] = pd.to_numeric(parts[1], errors="coerce").astype("Int64")
    df["ref"] = parts[2]
    df["alt"] = parts[3]
    return df


def lead_per_gene(df: pd.DataFrame) -> pd.DataFrame:
    """每个基因取 p 最小的变异作为 lead。"""
    idx = df.groupby("gene_id")["pval_nominal"].idxmin()
    return df.loc[idx].reset_index(drop=True)
