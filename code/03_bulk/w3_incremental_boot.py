"""
w3_incremental_boot.py — 增量价值模型的 ΔAUC 置信区间与可检测下限（回应 MoA 评审 D4）

原脚本只报点估计 + 随机模块经验 p，没有区间、也没说清"阴性"能不能支撑"无增量价值"。
本脚本补三件事：
  1. 对样本做分层 bootstrap，给 ΔAUC 的 95% 百分位区间
  2. 用 100 套等大随机模块的 ΔAUC 离散度给出可检测下限（80% 功效 ≈ 2.8 SD）
  3. 明确观测值落在零分布的哪个位置（z 与经验 p）

用法（复用 v2 版程序与随机模块）：
  PROG_PATH=.../IR_state_program_v2.tsv RM_PATH=.../random_modules_v2.tsv \
  OUT_TAG=_v2 BOOT=200 /opt/anaconda3/bin/python3 w3_incremental_boot.py
"""
from __future__ import annotations

import os
import pathlib
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from w3_incremental import (read_matrix, collapse, zscore, score_set, DATA, W2, RES, CYCLING)

BOOT = int(os.environ.get("BOOT", "200"))
TAG = os.environ.get("OUT_TAG", "_v2")


def cv_auc(X, y, seed=49, n_repeats=10):
    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=n_repeats, random_state=seed)
    aucs = []
    for tr, te in cv.split(X, y):
        m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1.0))
        m.fit(X[tr], y[tr])
        aucs.append(roc_auc_score(y[te], m.predict_proba(X[te])[:, 1]))
    return float(np.mean(aucs))


def main() -> None:
    raw, meta = read_matrix(DATA / "GSE120490_series_matrix.txt.gz")
    chars = meta["!Sample_characteristics_ch1"]

    def field(name):
        for row in chars:
            if row[0].lower().startswith(name.lower()):
                return [x.split(": ", 1)[-1] for x in row]
        raise KeyError(name)

    met = field("matastasis")
    grade = field("grade")
    mi = field("mi50pluspercent")
    y = np.array([1 if v.strip().lower() == "yes" else 0 for v in met])

    m = collapse(raw)
    z = zscore(m)

    prog_path = os.environ.get("PROG_PATH", str(W2.parent / "w4/results/IR_state_program_v2.tsv"))
    prog = pd.read_csv(prog_path, sep="\t", index_col=0)
    up_genes = list(prog.index[prog.direction == "UP"])
    ir = score_set(z, list(prog.index), up_genes)

    base = np.column_stack([
        z.loc[[g for g in CYCLING if g in z.index]].mean(axis=0).values,
        np.array([1 if v == "Yes" else 0 for v in grade]),
        np.array([1 if v == "Yes" else 0 for v in mi]),
    ])
    X = np.column_stack([base, ir.values])

    print(f"样本 {len(y)}；事件 {y.sum()}；程序 {len(prog)} 基因（可用 {sum(g in z.index for g in prog.index)}）")
    a0 = cv_auc(base, y)
    a1 = cv_auc(X, y)
    print(f"基座 AUC = {a0:.4f}；+IR AUC = {a1:.4f}；ΔAUC = {a1-a0:+.4f}")

    # ── 分层 bootstrap ───────────────────────────────────────────────
    rng = np.random.default_rng(49)
    idx_pos = np.where(y == 1)[0]
    idx_neg = np.where(y == 0)[0]
    deltas, auc_b, auc_f = [], [], []
    for b in range(BOOT):
        take = np.concatenate([
            rng.choice(idx_pos, len(idx_pos), replace=True),
            rng.choice(idx_neg, len(idx_neg), replace=True),
        ])
        Xb, yb = X[take], y[take]
        Bb = Xb[:, :base.shape[1]]
        try:
            b0 = cv_auc(Bb, yb, seed=1000 + b)
            b1 = cv_auc(Xb, yb, seed=1000 + b)
        except ValueError:
            continue
        auc_b.append(b0); auc_f.append(b1); deltas.append(b1 - b0)
        if (b + 1) % 25 == 0:
            print(f"  bootstrap {b+1}/{BOOT}  当前 ΔAUC 均值 {np.mean(deltas):+.4f}", flush=True)

    d = np.array(deltas)
    lo, hi = np.percentile(d, [2.5, 97.5])
    print(f"\nbootstrap（{len(d)} 次）ΔAUC 中位 {np.median(d):+.4f}，"
          f"95% CI [{lo:+.4f}, {hi:+.4f}]")
    print(f"  AUC 基座 95% CI [{np.percentile(auc_b,2.5):.4f}, {np.percentile(auc_b,97.5):.4f}]")
    print(f"  AUC +IR  95% CI [{np.percentile(auc_f,2.5):.4f}, {np.percentile(auc_f,97.5):.4f}]")

    # ── 可检测下限（基于等大随机模块的零分布离散度）───────────────────
    rm_path = os.environ.get("RM_PATH", str(W2 / "results/random_modules_v2.tsv"))
    rm = pd.read_csv(rm_path, sep="\t")
    null = pd.read_csv(RES / f"incremental_null{TAG}.tsv", sep="\t").delta_auc.values
    sd = float(np.std(null, ddof=1))
    z_obs = (a1 - a0 - null.mean()) / sd
    print(f"\n零分布（{len(null)} 套等大随机模块）：均值 {null.mean():+.4f}，SD {sd:.4f}")
    print(f"  观测 ΔAUC 相对零分布 z = {z_obs:+.2f}")
    print(f"  80% 功效可检测下限 ≈ 2.8 × SD = {2.8*sd:.4f}")
    print(f"  即：本次可排除 > {hi:.4f} 量级的增量（按 bootstrap 上界）")

    out = pd.DataFrame([{
        "auc_base": a0, "auc_full": a1, "delta_auc": a1 - a0,
        "boot_n": len(d), "boot_lo": lo, "boot_hi": hi,
        "boot_median": float(np.median(d)),
        "auc_base_lo": float(np.percentile(auc_b, 2.5)), "auc_base_hi": float(np.percentile(auc_b, 97.5)),
        "auc_full_lo": float(np.percentile(auc_f, 2.5)), "auc_full_hi": float(np.percentile(auc_f, 97.5)),
        "null_sd": sd, "null_mean": float(null.mean()), "z_vs_null": z_obs,
        "detectable_80pct": 2.8 * sd,
        "n": int(len(y)), "events": int(y.sum()),
        "cv": "5-fold x10 repeats (bootstrap), point estimate 5x20",
    }])
    f = RES / f"incremental_boot{TAG}.tsv"
    out.to_csv(f, sep="\t", index=False)
    print(f"\n已写 {f}")


if __name__ == "__main__":
    main()
