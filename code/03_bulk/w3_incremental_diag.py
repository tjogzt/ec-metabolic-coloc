"""
w3_incremental_diag.py — 诊断：为什么阳性对照（疾病签名）也不增加预测

三种可能：
  (a) 评分与终点本身无关            → 单变量 AUC ≈ 0.5
  (b) 评分有信号但被分级/肌层浸润完全中介 → 单变量 AUC > 0.6，加入基座后消失，且与分级强相关
  (c) 评分方向正确但样本量不足      → 单变量 AUC 0.55–0.60，ΔAUC 置信区间跨越 0
"""
from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, pointbiserialr, spearmanr
from sklearn.metrics import roc_auc_score

import w3_incremental as W

HERE = pathlib.Path(__file__).resolve().parent
RES = HERE / "results"

raw, meta = W.read_matrix(W.DATA / "GSE120490_series_matrix.txt.gz")
chars = meta["!Sample_characteristics_ch1"]
def field(name):
    for row in chars:
        if row[0].lower().startswith(name.lower()):
            return [x.split(": ", 1)[-1] for x in row]
    raise KeyError(name)

met = np.array([1 if v == "Yes" else 0 for v in field("matastasis")])
grade = np.array([1 if v == "Yes" else 0 for v in field("grade")])
mi = np.array([1 if v == "Yes" else 0 for v in field("mi50pluspercent")])

m = W.collapse(raw); z = W.zscore(m)
prog = pd.read_csv(W.W2 / "results/IR_program_v1.tsv", sep="\t", index_col=0)
sig = pd.read_csv(W.W2 / "results/EC_signature_v2_prolifadj_compact.tsv", sep="\t")
sc = {
    "IR": W.score_set(z, list(prog.index), list(prog.index[prog.direction == "UP"])),
    "ECsig": W.score_set(z, list(sig.gene), list(sig.gene[sig.direction == "UP"])),
    "Prolif": z.loc[[g for g in W.CYCLING if g in z.index]].mean(axis=0),
}

print("=== 单变量 AUC（对转移） ===")
rows = []
for k, s in sc.items():
    auc = roc_auc_score(met, s.values)
    auc = max(auc, 1 - auc)
    u = mannwhitneyu(s[met == 1], s[met == 0])
    rows.append({"评分": k, "AUC": round(auc, 3), "MWU p": f"{u.pvalue:.2g}"})
print(pd.DataFrame(rows).to_string(index=False))

print("\n=== 评分与临床变量的关系（Spearman） ===")
rows = []
for k, s in sc.items():
    for nm, v in [("分级", grade), ("肌层浸润≥50%", mi), ("转移", met)]:
        r, p = spearmanr(s.values, v)
        rows.append({"评分": k, "变量": nm, "rho": round(r, 3), "p": f"{p:.2g}"})
print(pd.DataFrame(rows).to_string(index=False))

print("\n=== 去掉分级与肌层浸润后，评分能否单独预测转移 ===")
X = np.column_stack([sc["ECsig"].values])
from sklearn.model_selection import cross_val_score, StratifiedKFold
cv = StratifiedKFold(5, shuffle=True, random_state=49)
for k, s in sc.items():
    a = cross_val_score(W.make_pipeline(W.StandardScaler(), W.LogisticRegression(max_iter=2000)),
                        s.values.reshape(-1, 1), met, cv=cv, scoring="roc_auc").mean()
    print(f"  {k:8s} 单独 CV AUC = {a:.3f}")

print("\n=== 反向：临床变量单独 ===")
base = np.column_stack([grade, mi])
print(f"  分级+肌层浸润 CV AUC = {W.cv_auc(base, met):.3f}")
