"""复核 bulk 层 IR(v2) 评分与增殖评分的相关，以及控制增殖后的增量预测。
沿用 w3_incremental 的读取/折叠/标准化函数，保证与锁定管线同口径。
"""
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import w3_incremental as W

raw, meta = W.read_matrix(W.DATA / "GSE120490_series_matrix.txt.gz")
m = W.collapse(raw); z = W.zscore(m)

prog2 = pd.read_csv(W.HERE.parent / "w4/results/IR_state_program_v2.tsv", sep="\t", index_col=0)
print("v2 程序列名:", list(prog2.columns)[:5], "行数", len(prog2))
up = list(prog2.index[prog2.direction == "UP"]) if "direction" in prog2.columns else []

ir = W.score_set(z, list(prog2.index), up)
prolif = z.loc[[g for g in W.CYCLING if g in z.index]].mean(axis=0)

common = ir.notna() & prolif.notna()
r, p = spearmanr(ir[common], prolif[common])
print(f"\nIR(v2) vs 增殖  Spearman rho = {r:+.3f}  p = {p:.3e}  n = {int(common.sum())}")

# v1 同样算一遍以对照
prog1 = pd.read_csv(W.W2 / "results/IR_program_v1.tsv", sep="\t", index_col=0)
ir1 = W.score_set(z, list(prog1.index), list(prog1.index[prog1.direction == "UP"]))
c2 = ir1.notna() & prolif.notna()
r1, p1 = spearmanr(ir1[c2], prolif[c2])
print(f"IR(v1) vs 增殖  Spearman rho = {r1:+.3f}  p = {p1:.3e}  n = {int(c2.sum())}")
