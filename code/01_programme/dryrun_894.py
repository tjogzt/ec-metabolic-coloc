#!/usr/bin/env python3

import os as _os
import sys as _sys
import pathlib as _pathlib

_root = _pathlib.Path(__file__).resolve()
while not (_root / "config.py").exists():
    if _root.parent == _root:
        raise RuntimeError("config.py not found; set ECIR_ROOT")
    _root = _root.parent
if str(_root) not in _sys.path:
    _sys.path.insert(0, str(_root))
import config as CFG
"""894 程序的干跑：在不触碰已锁定 862 文件的前提下，按同一冻结规则重建程序，
并计算三队列一致率，用于判断 48 个基因的变动是否改变结论。

输出写到 w4/results/dryrun_894/ 下，不覆盖任何既有文件。
"""
import gzip
import io
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr, binomtest

ROOT = CFG.ROOT
W2, W4 = ROOT / "w2", ROOT / "w4"
OUT = W4 / "results" / "dryrun_894"
OUT.mkdir(parents=True, exist_ok=True)

# ---------- 1. 按冻结规则重建发现集 ----------
D = pd.read_csv(W4 / "DE_GSE25401_adipose_obese_state.tsv", sep="\t", index_col=0)
sig = D[(D["adj.P.Val"] < 0.05) & (D.logFC.abs() >= 0.3)]

with gzip.open(W2 / "data/GSE17025_series_matrix.txt.gz", "rt") as fh:
    lines = fh.readlines()
i0 = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_begin"))
i1 = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_end"))
m = pd.read_csv(io.StringIO("".join(lines[i0 + 1:i1])), sep="\t", quotechar='"', index_col=0)
m = np.log2(m.apply(pd.to_numeric, errors="coerce").clip(lower=0) + 1)

ann = pd.read_csv(W2 / "data/GPL570_annot.tsv", sep="\t").set_index("ID")["Gene Symbol"]
m = m.loc[m.index.intersection(ann.index)]
m["sym"] = ann.reindex(m.index).values
m = m.dropna(subset=["sym"])
m = m[~m.sym.isin(["---", ""])]
ec = m.groupby("sym").mean(numeric_only=True).mean(axis=1)
thr = ec.quantile(0.25)

new = sig[sig.index.map(ec).fillna(-np.inf) > thr].copy()
new["direction"] = np.where(new.logFC > 0, "UP", "DOWN")
new.index.name = "gene"
new.to_csv(OUT / "IR_state_program_894_dryrun.tsv", sep="\t")
print(f"894 干跑程序：{len(new)} 基因（UP {(new.direction=='UP').sum()} / DOWN {(new.direction=='DOWN').sum()}）")

old = pd.read_csv(W4 / "results/IR_state_program_v2.tsv", sep="\t", index_col=0)
# 锁定文件的列名是 logFC_adipose_obese，干跑生成的是 logFC —— 统一成 old_fc
_fc = [c for c in old.columns if c.lower().startswith("logfc")][0]
old = old.rename(columns={_fc: "logFC"})

# ---------- 2. 三队列一致率（与 p0_concordance_v2.py 同口径）----------
COH = [
    ("GSE15773", "DE_GSE15773_adipose_IR_state.tsv", "脂肪 IR vs IS"),
    ("GSE18732", "DE_GSE18732_muscle_DM_state.tsv", "肌肉 T2D vs NGT"),
    ("GSE13070", "DE_GSE13070_muscle_baseline_IR.tsv", "肌肉 钳夹 IR vs IS"),
]
rows = []
for tag, fn, label in COH:
    d = pd.read_csv(W4 / fn, sep="\t", index_col=0)
    for name, prog in (("862", old), ("894", new)):
        common = prog.index.intersection(d.index)
        rec = d.loc[common]
        pdir = np.sign(prog.loc[common, "logFC"].values)
        vdir = np.sign(rec.logFC.values)
        agree = int((pdir == vdir).sum())
        n = len(common)
        conc_all = 100 * agree / n
        resp = rec[rec["P.Value"] < 0.05]
        rdir = np.sign(resp.logFC.values)
        prd = np.sign(prog.loc[resp.index, "logFC"].values)
        ra = int((rdir == prd).sum())
        rn = len(resp)
        rho = spearmanr(prog.loc[common, "logFC"], rec.logFC).statistic
        rows.append(dict(program=name, cohort=tag, label=label, n=n,
                         conc_all=round(conc_all, 2), agree=agree,
                         resp_n=rn, resp_agree=ra,
                         conc_resp=round(100 * ra / rn, 2) if rn else np.nan,
                         rho=round(rho, 4)))

res = pd.DataFrame(rows)
res.to_csv(OUT / "concordance_862_vs_894.tsv", sep="\t", index=False)
pd.set_option("display.width", 160)
print("\n=== 862 vs 894 一致率对照 ===")
print(res[["cohort", "program", "n", "conc_all", "resp_n", "conc_resp", "rho"]].to_string(index=False))

# ---------- 3. 直接对比：结论是否改变 ----------
print("\n=== 差异 ===")
for tag, _, _ in COH:
    a = res[(res.cohort == tag) & (res.program == "862")].iloc[0]
    b = res[(res.cohort == tag) & (res.program == "894")].iloc[0]
    print(f"{tag}: 全基因 {a.conc_all}→{b.conc_all} ({b.conc_all-a.conc_all:+.2f}pp) | "
          f"响应 {a.conc_resp}({a.resp_agree}/{a.resp_n})→{b.conc_resp}({b.resp_agree}/{b.resp_n}) | "
          f"ρ {a.rho}→{b.rho} ({b.rho-a.rho:+.4f})")

r862 = [res[(res.cohort == t) & (res.program == "862")].iloc[0].rho for t, _, _ in COH]
r894 = [res[(res.cohort == t) & (res.program == "894")].iloc[0].rho for t, _, _ in COH]
print(f"\n梯度 862: {r862[0]} > {r862[1]} > {r862[2]}  → {all(r862[i] > r862[i+1] for i in range(2))}")
print(f"梯度 894: {r894[0]} > {r894[1]} > {r894[2]}  → {all(r894[i] > r894[i+1] for i in range(2))}")
