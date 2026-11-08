"""
w3_load.py — 读取 W3 核心数据集，构建统一 AnnData 与细胞级打分

流程：发现文件 → 建 AnnData → QC → 归一化 → 谱系注释 → 打分
打分对象：IR 程序（UP/DOWN 分列）、疾病签名（紧凑版）、增殖分数、100 套随机模块
"""
from __future__ import annotations

import glob
import os
import pathlib
import re
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.io import mmread

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE / "data"
RES = HERE / "results"
RES.mkdir(exist_ok=True)
sc.settings.verbosity = 1

MARKERS = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],
    "Fibroblast": ["DCN", "LUM", "COL1A1", "COL1A2", "PDGFRA"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "PLVAP"],
    "Pericyte/SMC": ["RGS5", "ACTA2", "DES", "MYH11", "TAGLN"],
    "T/NK": ["PTPRC", "CD3D", "CD3E", "CD2", "NKG7", "GNLY"],
    "Myeloid": ["CD68", "CD14", "LYZ", "C1QA", "ITGAX"],
    "B/Plasma": ["MS4A1", "CD79A", "MZB1", "JCHAIN"],
    "Mast": ["TPSAB1", "CPA3", "MS4A2"],
}


def find_10x(sample_dir: pathlib.Path) -> dict | None:
    """在一个样本目录里找 10x 三件套或 h5。

    只看本层（iterdir），不用 rglob：递归会让父目录把子目录的文件也搜刮进来，
    再取「字母序第一个 matrix + 同名 barcodes」就会跨样本配错，
    随后按 matrix 路径去重时真实样本被当成重复丢掉。
    """
    files = [p for p in sample_dir.iterdir() if p.is_file()]
    if not files:
        return None
    h5 = [p for p in files if p.suffix == ".h5"]
    if h5:
        return {"h5": str(sorted(h5)[0])}
    mtx = [p for p in files if re.search(r"matrix\.mtx(\.gz)?$", p.name)]
    if not mtx:
        return None
    mtx = sorted(mtx)[0]
    pref = re.sub(r"matrix\.mtx(\.gz)?$", "", mtx.name)
    d = mtx.parent
    bc = [p for p in files if p.name in (f"{pref}barcodes.tsv.gz", f"{pref}barcodes.tsv")]
    ft = [p for p in files if p.name in (f"{pref}features.tsv.gz", f"{pref}genes.tsv.gz",
                                         f"{pref}features.tsv", f"{pref}genes.tsv")]
    if not bc or not ft:
        return None
    return {"mtx": str(mtx), "barcodes": str(bc[0]), "features": str(ft[0])}


def read_10x(spec: dict) -> ad.AnnData:
    if "h5" in spec:
        return sc.read_10x_h5(spec["h5"], gex_only=True)
    feat = pd.read_csv(spec["features"], sep="\t", header=None)
    sym = feat[1].astype(str)
    X = mmread(spec["mtx"]).T.tocsr().astype(np.float32)
    bc = pd.read_csv(spec["barcodes"], sep="\t", header=None)[0].astype(str)
    a = ad.AnnData(X=X, obs=pd.DataFrame(index=bc), var=pd.DataFrame(index=sym.values))
    a.var_names_make_unique()
    return a


def load_dataset(gse: str) -> list[ad.AnnData]:
    root = DATA / gse
    if not root.exists():
        return []
    # 每个子目录/顶层算一个候选样本
    cands = []
    for d in [root] + [p for p in root.rglob("*") if p.is_dir()]:
        spec = find_10x(d)
        if spec:
            cands.append((d, spec))
    # 去掉父目录重复（父目录会重复找到子目录文件）
    seen = set()
    uniq = []
    for d, spec in cands:
        key = (spec.get("mtx") or spec.get("h5"))
        if key in seen:
            continue
        seen.add(key)
        uniq.append((d, spec))
    out = []
    for d, spec in uniq:
        try:
            a = read_10x(spec)
        except Exception as e:  # noqa: BLE001
            print(f"  [skip] {gse}/{d.name}: {e}", file=sys.stderr)
            continue
        a.obs["dataset"] = gse
        a.obs["sample"] = d.name if d != root else gse
        a.obs_names = [f"{gse}_{d.name}_{b}" for b in a.obs_names]
        out.append(a)
        print(f"  [{gse}] {d.name}: {a.n_obs} 细胞 × {a.n_vars} 基因")
    return out


def main() -> None:
    genes = [g for g in ["GSE278879", "GSE268903", "GSE212217", "GSE251923"]
             if (DATA / g).exists()]
    print("加载：", genes)
    parts = []
    for g in genes:
        parts += load_dataset(g)
    if not parts:
        print("没有可用的表达矩阵", file=sys.stderr)
        sys.exit(1)

    print(f"\n合并 {len(parts)} 个样本…")
    a = ad.concat(parts, join="outer", label="batch_id", index_unique=None, fill_value=0)
    print(f"合并后 {a.n_obs:,} 细胞 × {a.n_vars:,} 基因")
    a.var_names = a.var_names.str.upper()

    # ---- QC ----
    a.var["mt"] = a.var_names.str.startswith(("MT-",))
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None)
    n0 = a.n_obs
    sc.pp.filter_cells(a, min_genes=200)
    sc.pp.filter_genes(a, min_cells=3)
    a = a[a.obs.pct_counts_mt < 20].copy()
    print(f"QC：{n0:,} → {a.n_obs:,} 细胞（保留 {100*a.n_obs/n0:.1f}%），"
          f"中位基因数 {a.obs.n_genes_by_counts.median():.0f}")
    a.write(RES / "atlas_core_qc.h5ad")

    # ---- 归一化 / 降维 / 聚类 ----
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=2000, batch_key="dataset")
    b = a[:, a.var.highly_variable].copy()
    sc.pp.scale(b, max_value=10)
    sc.tl.pca(b, n_comps=30, svd_solver="arpack")
    try:
        import harmonypy  # noqa: F401
        sc.external.pp.harmony_integrate(b, "dataset")
        rep = "X_pca_harmony"
        print("已用 harmony 按数据集整合")
    except Exception as e:  # noqa: BLE001
        rep = "X_pca"
        print(f"harmony 不可用（{e}），改用原始 PCA")
    sc.pp.neighbors(b, n_neighbors=15, use_rep=rep)
    sc.tl.leiden(b, resolution=1.0, key_added="leiden")
    a.obs["leiden"] = b.obs["leiden"].values
    a.obsm["X_pca"] = b.obsm["X_pca"]
    if rep in b.obsm:
        a.obsm[rep] = b.obsm[rep]
    sc.tl.umap(b)
    a.obsm["X_umap"] = b.obsm["X_umap"]
    print(f"聚类：{a.obs.leiden.nunique()} 个簇")

    # ---- 谱系注释 ----
    sc.pp.scale(a, max_value=10) if False else None
    score = pd.DataFrame(index=a.obs_names)
    for lin, mk in MARKERS.items():
        present = [g for g in mk if g in a.var_names]
        if not present:
            continue
        _safe = re.sub(r"[^A-Za-z0-9_.-]", "_", lin)
        sc.tl.score_genes(a, present, score_name=f"lin_{_safe}", use_raw=False if a.raw is None else True)
        score[lin] = a.obs[f"lin_{_safe}"].values
    # anndata 不允许 obs 的分类名或列名含 "/"（如 Pericyte/SMC），写盘前统一替换
    a.obs["lineage"] = [re.sub(r"[^A-Za-z0-9_.-]", "_", str(x)) for x in score.idxmax(axis=1).values]
    a.obs["lineage_conf"] = (score.max(axis=1) - score.apply(
        lambda r: r.nlargest(2).iloc[-1], axis=1)).values
    print("\n谱系分布：")
    print(a.obs.lineage.value_counts().to_string())

    # ---- 三个程序打分 ----
    prog = pd.read_csv(HERE.parent / "w2/results/IR_program_v1.tsv", sep="\t", index_col=0)
    sig = pd.read_csv(HERE.parent / "w2/results/EC_signature_v2_prolifadj_compact.tsv", sep="\t")
    CYCLING = ["MKI67", "PCNA", "CCNB1", "CCNB2", "CCNA2", "CDK1", "CDC20", "BUB1", "BUB1B",
               "TOP2A", "AURKA", "AURKB", "PLK1", "TPX2", "UBE2C", "RRM2", "TYMS", "MCM2",
               "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "CENPF", "KIF11", "KIF23", "ANLN",
               "PRC1", "ECT2", "NCAPG", "ASPM", "DLGAP5", "CENPE", "SMC4", "NUSAP1"]

    up = [g for g in prog.index[prog.direction == "UP"] if g in a.var_names]
    dn = [g for g in prog.index[prog.direction == "DOWN"] if g in a.var_names]
    print(f"\nIR 程序覆盖率：UP {len(up)}/{sum(prog.direction=='UP')}，"
          f"DOWN {len(dn)}/{sum(prog.direction=='DOWN')}")
    sc.tl.score_genes(a, up, score_name="IR_UP", ctrl_size=len(up))
    sc.tl.score_genes(a, dn, score_name="IR_DOWN", ctrl_size=len(dn))
    a.obs["IR_score"] = a.obs.IR_UP - a.obs.IR_DOWN

    sg_up = [g for g in sig.gene[sig.direction == "UP"] if g in a.var_names]
    sg_dn = [g for g in sig.gene[sig.direction == "DOWN"] if g in a.var_names]
    print(f"疾病签名覆盖率：UP {len(sg_up)}/{sum(sig.direction=='UP')}，"
          f"DOWN {len(sg_dn)}/{sum(sig.direction=='DOWN')}")
    sc.tl.score_genes(a, sg_up, score_name="ECsig_UP", ctrl_size=len(sg_up))
    sc.tl.score_genes(a, sg_dn, score_name="ECsig_DOWN", ctrl_size=len(sg_dn))
    a.obs["ECsig_score"] = a.obs.ECsig_UP - a.obs.ECsig_DOWN

    cyc = [g for g in CYCLING if g in a.var_names]
    sc.tl.score_genes(a, cyc, score_name="Prolif", ctrl_size=len(cyc))
    print(f"增殖分数用 {len(cyc)} 个基因")

    a.write(RES / "atlas_core_scored.h5ad")
    print(f"\n已写 {RES/'atlas_core_scored.h5ad'}")

    # 汇总
    s = (a.obs.groupby(["dataset", "lineage"], observed=True)
         .agg(n=("IR_score", "size"),
              IR=("IR_score", "mean"),
              ECsig=("ECsig_score", "mean"),
              Prolif=("Prolif", "mean"))
         .round(3).reset_index())
    s.to_csv(RES / "state_program_table.tsv", sep="\t", index=False)
    print("\n状态—程序对应表：")
    print(s.to_string(index=False))


if __name__ == "__main__":
    main()
