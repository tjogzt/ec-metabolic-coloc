"""
w2_program.py — 构建并锁定胰岛素响应程序（IR program v1）

输入：
  DE_GSE24215_insulin_vs_basal.tsv   骨骼肌，10 受试者配对，胰岛素 vs 基础
  DE_GSE24422_insulin_adipocyte.tsv  脂肪细胞，±胰岛素
  data/GSE17025_series_matrix.txt.gz 内膜癌组织（用于"是否在内膜表达"过滤）
输出：
  IR_program_v1.tsv + .sha256 + random_modules_v1.tsv
"""
from __future__ import annotations
import hashlib
import pathlib
import re
import subprocess

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE / "data"
RNG = np.random.default_rng(20261102)

# ---------- 1. 读两套 DE ----------
A = pd.read_csv(HERE / "DE_GSE24215_insulin_vs_basal.tsv", sep="\t", index_col=0)
B = pd.read_csv(HERE / "DE_GSE24422_insulin_adipocyte.tsv", sep="\t", index_col=0)
print(f"GSE24215（骨骼肌）基因 {len(A):,}；GSE24422（脂肪细胞）基因 {len(B):,}")

common = A.index.intersection(B.index)
df = pd.DataFrame({
    "logFC_muscle": A.loc[common, "logFC"],
    "FDR_muscle": A.loc[common, "adj.P.Val"],
    "logFC_adip": B.loc[common, "logFC"],
    "FDR_adip": B.loc[common, "adj.P.Val"],
}, index=common)

# ---------- 2. 内膜表达过滤（GSE17025，103 例内膜癌组织，GPL570） ----------
def ec_expression() -> pd.Series:
    """GSE17025（103 例内膜癌组织，GPL570）基因级平均表达，log2 转换。"""
    import gzip
    import io
    f = DATA / "GSE17025_series_matrix.txt.gz"
    with gzip.open(f, "rt") as fh:
        lines = fh.readlines()
    i0 = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_begin"))
    i1 = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_end"))
    m = pd.read_csv(io.StringIO("".join(lines[i0 + 1:i1])),
                    sep="\t", quotechar='"', index_col=0)
    m = m.apply(pd.to_numeric, errors="coerce")
    m = np.log2(m.clip(lower=0) + 1)
    annot = pd.read_csv(DATA / "GPL570_annot.tsv", sep="\t")
    a = annot.set_index("ID")["Gene Symbol"]
    m = m.loc[m.index.intersection(a.index)]
    m["sym"] = a.reindex(m.index).values
    m = m.dropna(subset=["sym"])
    m = m[~m.sym.isin(["---", ""])]
    return m.groupby("sym").mean(numeric_only=True).mean(axis=1)


ec_expr = ec_expression()
thr = ec_expr.quantile(0.25)
print(f"内膜癌组织可测基因 {len(ec_expr):,}；第 25 百分位阈值 {thr:.2f}")

df["expr_in_EC"] = df.index.map(ec_expr)
df["expressed_in_EC"] = df["expr_in_EC"].fillna(-np.inf) > thr

# ---------- 3. 程序分级 ----------
sigA = (df.FDR_muscle < 0.05) & (df.logFC_muscle.abs() >= 0.3)
same_dir = np.sign(df.logFC_muscle) == np.sign(df.logFC_adip)
df["same_direction"] = same_dir

core = sigA & same_dir & (df.FDR_adip < 0.05) & df.expressed_in_EC
ext = sigA & same_dir & df.expressed_in_EC & ~core

df["tier"] = np.where(core, "core", np.where(ext, "extended", "excluded"))
df["direction"] = np.where(df.logFC_muscle > 0, "UP", "DOWN")

prog = df[df.tier.isin(["core", "extended"])].copy()
prog = prog.sort_values(["tier", "FDR_muscle"])
print(f"\n核心程序（两数据集同向且均显著、内膜可测）：{int(core.sum())} 个基因")
print(f"扩展程序（骨骼肌显著同向、内膜可测）      ：{int(ext.sum())} 个基因")
print(f"  UP {int((prog.direction=='UP').sum())} / DOWN {int((prog.direction=='DOWN').sum())}")

out = HERE / "results"
out.mkdir(exist_ok=True)
cols = ["direction", "tier", "logFC_muscle", "FDR_muscle", "logFC_adip", "FDR_adip",
        "same_direction", "expr_in_EC", "expressed_in_EC"]
prog[cols].to_csv(out / "IR_program_v1.tsv", sep="\t", float_format="%.6g")
h = hashlib.sha256((out / "IR_program_v1.tsv").read_bytes()).hexdigest()
(out / "IR_program_v1.sha256").write_text(f"{h}  IR_program_v1.tsv\n")
print(f"\n哈希 {h}")

# ---------- 4. 表达匹配的随机模块 ----------
expr_all = ec_expr.dropna()
deciles = pd.qcut(expr_all.rank(method="first"), 10, labels=False)
pool = {}   # decile -> gene list
prog_expr = expr_all.reindex(prog.index).dropna()
for g, e in prog_expr.items():
    d = int(deciles.get(g, -1))
    if d >= 0:
        pool.setdefault(d, []).append(g)
allpool = {}
for g, e in expr_all.items():
    d = int(deciles.get(g, -1))
    if d >= 0 and g not in set(prog.index):
        allpool.setdefault(d, []).append(g)

need = prog_expr.groupby(deciles.reindex(prog_expr.index)).size().to_dict()
recs = []
for rep in range(100):
    picks = []
    for d, k in need.items():
        cand = allpool.get(int(d), [])
        if len(cand) < k:
            continue
        picks += list(RNG.choice(cand, size=k, replace=False))
    for g in picks:
        recs.append({"module": f"R{rep:03d}", "gene": g})
rm = pd.DataFrame(recs).drop_duplicates()
rm.to_csv(out / "random_modules_v1.tsv", sep="\t", index=False)
print(f"随机模块：{rm.module.nunique()} 套，中位基因数 "
      f"{int(rm.groupby('module').size().median())}（程序 {len(prog)} 个基因，按表达十分位匹配）")
