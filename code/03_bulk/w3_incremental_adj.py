"""
w3_incremental_adj.py — 轴 C 收尾：增殖残差化后的 IR 评分是否仍有独立信息

做法：
  1. 把 IR 评分对增殖分数做线性回归，取残差（IR_resid）
  2. 检验 IR_resid 对转移的单变量 AUC
  3. 检验基座 + IR_resid 的 CV ΔAUC，并与 100 套随机模块的残差化版本比较
  4. 分臂重复（上调臂 / 下调臂）
"""
from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import roc_auc_score

import w3_incremental as W

HERE = pathlib.Path(__file__).resolve().parent
RES = HERE / "results"

raw, meta = W.read_matrix(W.DATA / "GSE120490_series_matrix.txt.gz")
chars = meta["!Sample_characteristics_ch1"]
def field(n):
    for r in chars:
        if r[0].lower().startswith(n.lower()):
            return [x.split(": ", 1)[-1] for x in r]
    raise KeyError(n)

met = np.array([1 if v == "Yes" else 0 for v in field("matastasis")])
grade = np.array([1 if v == "Yes" else 0 for v in field("grade")])
mi = np.array([1 if v == "Yes" else 0 for v in field("mi50pluspercent")])
m = W.collapse(raw); z = W.zscore(m)
prolif = z.loc[[g for g in W.CYCLING if g in z.index]].mean(axis=0)

prog = pd.read_csv(W.W2 / "results/IR_program_v1.tsv", sep="\t", index_col=0)
up_g = list(prog.index[prog.direction == "UP"])
dn_g = list(prog.index[prog.direction == "DOWN"])
scores = {
    "IR总分": W.score_set(z, list(prog.index), up_g),
    "IR上调臂": W.score_set(z, up_g),
    "IR下调臂": W.score_set(z, dn_g),
}


def resid(s: pd.Series) -> np.ndarray:
    lr = LinearRegression().fit(prolif.values.reshape(-1, 1), s.values)
    return s.values - lr.predict(prolif.values.reshape(-1, 1))


base = np.column_stack([prolif.values, grade, mi])
auc_base = W.cv_auc(base, met)
print(f"基座 CV AUC = {auc_base:.3f}\n")
print(f"{'评分':10s} {'原AUC':>8s} {'残差AUC':>9s} {'ΔAUC(基座+残差)':>18s}")
for k, s in scores.items():
    a_raw = max(roc_auc_score(met, s.values), 1 - roc_auc_score(met, s.values))
    r = resid(s)
    a_res = max(roc_auc_score(met, r), 1 - roc_auc_score(met, r))
    d = W.cv_auc(np.column_stack([base, r]), met) - auc_base
    print(f"{k:10s} {a_raw:10.3f} {a_res:9.3f} {d:+16.4f}")

# 随机模块的残差化零分布
rm = pd.read_csv(W.W2 / "results/random_modules_v1.tsv", sep="\t")
null_raw, null_res = [], []
for name, g in rm.groupby("module"):
    s = W.score_set(z, list(g.gene))
    if s is None:
        continue
    null_raw.append(W.cv_auc(np.column_stack([base, s.values]), met) - auc_base)
    null_res.append(W.cv_auc(np.column_stack([base, resid(s)]), met) - auc_base)
null_res = np.array(null_res); null_raw = np.array(null_raw)
print(f"\n随机模块（{len(null_res)} 套）：")
print(f"  原始    ΔAUC 均值 {null_raw.mean():+.4f}，95% 分位 {np.percentile(null_raw,95):+.4f}")
print(f"  残差化  ΔAUC 均值 {null_res.mean():+.4f}，95% 分位 {np.percentile(null_res,95):+.4f}")

d_ir = W.cv_auc(np.column_stack([base, resid(scores['IR总分'])]), met) - auc_base
p = float((null_res >= d_ir).mean())
print(f"\nIR 总分残差化 ΔAUC = {d_ir:+.4f}，经验 p = {p:.2f}")
verdict = "不通过（预注册阈值 ΔAUC≥0.05）" if d_ir < 0.05 else "通过"
print(f"轴 C 判定：{verdict}")

pd.DataFrame([{
    "auc_base": auc_base, "delta_auc_IR_resid": d_ir, "p_emp": p,
    "null_res_mean": null_res.mean(), "null_res_p95": np.percentile(null_res, 95),
}]).to_csv(RES / "incremental_adjusted.tsv", sep="\t", index=False)
print(f"已写 {RES/'incremental_adjusted.tsv'}")
