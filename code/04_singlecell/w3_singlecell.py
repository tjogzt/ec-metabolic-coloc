"""
w3_singlecell.py — W3 单细胞投影：把 IR 程序与疾病签名打到细胞状态上

数据：本地 EC_REF_v1/ec_integrated_v1.h5ad
      141,316 细胞 × 36,601 基因；来源 GSE278879(11 例 snRNA) + GSE251923(2 例 scRNA)
      正常 4 例（NE1-4）、恶性 7 例（HEEC/MEEC/LEEC）

统计纪律：
  · 细胞数是伪重复——所有组间检验在**样本层**（7 vs 4）做，细胞层只作描述
  · 增殖混杂单独处理：先看 IR 分数与增殖的相关，再在非增殖细胞中复算
  · 随机模块（100 套）给出零分布，判定 IR 效应是否超过同规模随机基因集
"""
from __future__ import annotations

import os
import pathlib

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import mannwhitneyu, spearmanr

HERE = pathlib.Path(__file__).resolve().parent
W2 = HERE.parent / "w2"
RES = HERE / "results"
RES.mkdir(exist_ok=True)
H5 = "/Volumes/tjogzt4T/data/scrnaseq_reference/EC_REF_v1/ec_integrated_v1.h5ad"

CYCLING = ["MKI67", "PCNA", "CCNB1", "CCNB2", "CCNA2", "CDK1", "CDC20", "BUB1", "BUB1B",
           "TOP2A", "AURKA", "AURKB", "PLK1", "TPX2", "UBE2C", "RRM2", "TYMS", "MCM2",
           "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "CENPF", "KIF11", "KIF23", "ANLN",
           "PRC1", "ECT2", "NCAPG", "ASPM", "DLGAP5", "CENPE", "SMC4", "NUSAP1"]


def build_sets() -> tuple[dict[str, list[str]], int]:
    prog_path = os.environ.get("PROG_PATH", str(W2 / "results/IR_program_v1.tsv"))
    prog = pd.read_csv(prog_path, sep="\t", index_col=0)
    sig_path = os.environ.get("SIG_PATH", str(W2 / "results/EC_signature_v2_prolifadj_compact.tsv"))
    sig = pd.read_csv(sig_path, sep="\t")
    rm = pd.read_csv(os.environ.get("RM_PATH", str(W2 / "results/random_modules_v1.tsv")), sep="\t")

    sets: dict[str, list[str]] = {
        "IR_UP": list(prog.index[prog.direction == "UP"]),
        "IR_DOWN": list(prog.index[prog.direction == "DOWN"]),
        "ECsig_UP": list(sig.gene[sig.direction == "UP"]),
        "ECsig_DOWN": list(sig.gene[sig.direction == "DOWN"]),
        "Prolif": CYCLING,
    }
    for name, g in rm.groupby("module"):
        sets[name] = list(g.gene)
    return sets, len(sets)


def main() -> None:
    print("读取 h5ad（backed）…")
    a = ad.read_h5ad(H5, backed="r")
    obs = a.obs[["sample", "sample_class", "coarse_label"]].copy()
    print(f"{a.n_obs:,} 细胞 × {a.n_vars:,} 基因")

    sets, n_sets = build_sets()
    genes = sorted({g for v in sets.values() for g in v})
    idx = a.var_names.get_indexer(genes)
    ok = idx >= 0
    genes_ok = [g for g, k in zip(genes, ok) if k]
    idx_ok = idx[ok]
    print(f"基因集并集 {len(genes)}，其中有表达检测的 {len(genes_ok)}")

    # 一次性取出所需列
    X = a[:, genes_ok].to_memory().X
    X = sp.csr_matrix(X)
    print(f"子矩阵 {X.shape}，非零 {X.nnz:,}")

    mu = np.asarray(X.mean(axis=0)).ravel()
    sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
    sd = np.sqrt(np.maximum(sq - mu ** 2, 1e-8))
    gidx = {g: i for i, g in enumerate(genes_ok)}

    # 每个集合：score_i = mean_j((x_ij - mu_j)/sd_j) = x_i·w - c
    cols, rows, vals, consts = [], [], [], []
    names = []
    for k, (nm, gl) in enumerate(sets.items()):
        gj = [gidx[g] for g in gl if g in gidx]
        if len(gj) < 5:
            continue
        w = 1.0 / (sd[gj] * len(gj))
        rows += gj
        cols += [len(names)] * len(gj)
        vals += list(w)
        consts.append(float(np.sum(mu[gj] / sd[gj]) / len(gj)))
        names.append(nm)
    W = sp.csr_matrix((vals, (rows, cols)), shape=(len(genes_ok), len(names)))
    S = np.asarray((X @ W).toarray()) - np.array(consts)
    scores = pd.DataFrame(S, columns=names, index=np.arange(S.shape[0]))
    print(f"打分完成：{len(names)} 个集合")

    df = obs.reset_index(drop=True)
    df["IR_UP"] = scores.IR_UP.values
    df["IR_DOWN"] = scores.IR_DOWN.values
    df["IR_score"] = df.IR_UP - df.IR_DOWN
    df["ECsig_score"] = scores.ECsig_UP - scores.ECsig_DOWN
    df["Prolif"] = scores.Prolif
    rand_names = [n for n in names if n.startswith("R")]

    # ---- 细胞层描述（样本×细胞类型均值） ----
    cellmean = (df.groupby(["sample", "sample_class", "coarse_label"], observed=True)
                .agg(n=("IR_score", "size"), IR=("IR_score", "mean"),
                     IR_UP=("IR_UP", "mean"), IR_DN=("IR_DOWN", "mean"),
                     ECsig=("ECsig_score", "mean"), Prolif=("Prolif", "mean"))
                .reset_index())
    tag = os.environ.get("OUT_TAG", "")
    cellmean.to_csv(RES / f"state_program_table{tag}.tsv", sep="\t", index=False)
    print("\n=== 样本×细胞类型（上皮） ===")
    ep = cellmean[cellmean.coarse_label.str.contains("Epithelial")]
    print(ep.round(3).to_string(index=False))

    # ---- 样本层检验：上皮细胞的 IR 分数，恶性 vs 正常 ----
    print("\n=== 样本层检验（上皮细胞，恶性 7 例 vs 正常 4 例）===")
    ep_df = df[df.coarse_label.str.contains("Epithelial")]
    per_sample = ep_df.groupby(["sample", "sample_class"], observed=True).agg(
        n=("IR_score", "size"), IR=("IR_score", "mean"),
        ECsig=("ECsig_score", "mean"), Prolif=("Prolif", "mean")).reset_index()
    print(per_sample.round(3).to_string(index=False))

    def test(col: str, frame=per_sample) -> tuple[float, float, float]:
        x = frame.loc[frame.sample_class == "Malignant", col]
        y = frame.loc[frame.sample_class == "Normal", col]
        return x.mean() - y.mean(), mannwhitneyu(x, y).pvalue, (
            x.mean() - y.mean()) / (np.sqrt((x.var(ddof=1) + y.var(ddof=1)) / 2) + 1e-12)

    rows = []
    for label, col, frame in [("全部上皮：IR", "IR", per_sample),
                              ("全部上皮：疾病签名", "ECsig", per_sample),
                              ("全部上皮：增殖", "Prolif", per_sample)]:
        d, p, es = test(col, frame)
        rows.append({"比较": label, "差值(恶性−正常)": round(d, 3),
                     "MWU p": round(p, 3), "Cohen d": round(es, 2)})

    # 非增殖上皮（增殖分数下半）
    thr = ep_df.Prolif.median()
    quies = ep_df[ep_df.Prolif <= thr]
    qps = quies.groupby(["sample", "sample_class"], observed=True).agg(
        n=("IR_score", "size"), IR=("IR_score", "mean"),
        ECsig=("ECsig_score", "mean"), Prolif=("Prolif", "mean")).reset_index()
    d, p, es = test("IR", qps)
    rows.append({"比较": "非增殖上皮：IR", "差值(恶性−正常)": round(d, 3),
                 "MWU p": round(p, 3), "Cohen d": round(es, 2)})
    print("\n" + pd.DataFrame(rows).to_string(index=False))

    # ---- 增殖共线 ----
    r, p = spearmanr(ep_df.IR_score, ep_df.Prolif)
    print(f"\n上皮细胞中 IR 分数与增殖的相关：ρ={r:.3f}（p={p:.1e}，n={len(ep_df):,} 细胞）")

    # ---- 随机模块零分布（样本层） ----
    null = []
    for nm in rand_names:
        ep_df["_s"] = scores[nm].values[ep_df.index]
        ps = ep_df.groupby(["sample", "sample_class"], observed=True)["_s"].mean().reset_index()
        x = ps.loc[ps.sample_class == "Malignant", "_s"].mean()
        y = ps.loc[ps.sample_class == "Normal", "_s"].mean()
        null.append(x - y)
    null = np.array(null)
    d_ir = test("IR")[0]
    print(f"\n随机模块（{len(null)} 套）样本层差值：均值 {null.mean():+.4f}，"
          f"95% 分位 {np.percentile(null, 95):+.4f}，最大 {null.max():+.4f}")
    print(f"IR 实测差值 {d_ir:+.4f}，经验 p = {(null >= d_ir).mean():.3f}")

    pd.DataFrame({"delta_IR": [d_ir], "null_mean": [null.mean()],
                  "null_p95": [np.percentile(null, 95)],
                  "p_emp": [(null >= d_ir).mean()]}).to_csv(
        RES / f"singlecell_IR_test{tag}.tsv", sep="\t", index=False)
    print(f"\n已写 {RES/'state_program_table.tsv'} 与 {RES/'singlecell_IR_test.tsv'}")


if __name__ == "__main__":
    main()
