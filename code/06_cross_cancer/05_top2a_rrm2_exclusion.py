#!/usr/bin/env python3
"""A4 步骤 5：TOP2A/RRM2 剔除敏感性（R18 起作为投稿数字的落盘出处）。

用途：补充材料 S3 "Results in normal tissue" 与主稿对应句引用
"excluding both leaves the normal-tissue median correlation at +0.477 versus +0.483"。
此前 S3 曾引 "+0.466 versus +0.473 (Entrez-keyed reconstruction)"，该组数字
无法在本项目数据中复现（R18 主席核验），已替换为本脚本输出的可复现值。

规则与 03_analyse.py 完全一致：癌型内逐基因 z-score，score=mean(z_UP)-mean(z_DOWN)，
增殖模块 10 基因预写死；肿瘤 n>=50（30 型），正常 n>=20（13 型）。
剔除 = 从程序基因中移除 TOP2A、RRM2（与增殖模块共享的 2 基因）。

Author  : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created : 2026-11-20 (R18)
"""
import pathlib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

BASE = pathlib.Path("/Users/taozhu/my researches/gongbing/06_EC_IR")
CACHE = BASE / "a4/cache"
OUT = BASE / "a4/results"

z = np.load(CACHE / "expr.npz", allow_pickle=True)
X = z["X"]
genes = np.array([str(g) for g in z["genes"]])
samples = np.array([str(s) for s in z["samples"]])
gidx = {g: i for i, g in enumerate(genes)}
meta = pd.read_csv(CACHE / "meta.tsv", sep="\t").set_index("sample_id")
col_of = {s: i for i, s in enumerate(samples)}

pmap = pd.read_csv(OUT / "prog_entrez.tsv", sep="\t")
pmap = pmap[pmap.SYMBOL.isin(gidx)]
UP = pmap[pmap.direction == "UP"].SYMBOL.tolist()
DN = pmap[pmap.direction == "DOWN"].SYMBOL.tolist()
PROLIF = ["MKI67", "PCNA", "TOP2A", "CCNB1", "CDK1", "BUB1", "AURKA", "PLK1", "CENPE", "RRM2"]
prolif = [g for g in PROLIF if g in gidx]


def zmat(cols):
    sub = X[:, cols]
    mu = sub.mean(axis=1, keepdims=True)
    sd = sub.std(axis=1, keepdims=True)
    sd[sd == 0] = 1.0
    return (sub - mu) / sd


def median_rho(up_l, dn_l, sample_type, min_n):
    rhos = []
    for _ct, grp in meta[meta.sample_type == sample_type].groupby("cancer_type"):
        cols = [col_of[s] for s in grp.index if s in col_of]
        if len(cols) < min_n:
            continue
        Z = zmat(cols)
        su = Z[[gidx[g] for g in up_l if g in gidx]].mean(axis=0)
        sd_ = Z[[gidx[g] for g in dn_l if g in gidx]].mean(axis=0)
        sp = Z[[gidx[g] for g in prolif]].mean(axis=0)
        rhos.append(spearmanr(su - sd_, sp).statistic)
    return float(np.median(rhos)), len(rhos)


mt, nt = median_rho(UP, DN, "Primary Tumor", 50)
mn, nn = median_rho(UP, DN, "Solid Tissue Normal", 20)
up2 = [g for g in UP if g not in {"TOP2A", "RRM2"}]
dn2 = [g for g in DN if g not in {"TOP2A", "RRM2"}]
mt2, _ = median_rho(up2, dn2, "Primary Tumor", 50)
mn2, _ = median_rho(up2, dn2, "Solid Tissue Normal", 20)

print(f"full:            tumour {mt:+.4f} ({nt} types) | normal {mn:+.4f} ({nn} types)")
print(f"excl TOP2A/RRM2: tumour {mt2:+.4f}            | normal {mn2:+.4f}")
print("S3 引用值：+0.477 versus +0.483（normal）。")
