#!/usr/bin/env python3
"""B2-a: 为 894 基因版程序构建 100 套表达十分位匹配随机模块（与 862 版同法）。"""
import pathlib, sys
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent   # w4
sys.path.insert(0, str(HERE))
from w4_program_v2 import ec_expression            # GSE17025 内膜表达（同 862 版匹配池）

PROG = HERE / "results/dryrun_894/IR_state_program_894_dryrun.tsv"
OUTD = HERE / "results/dryrun_894"
OUT = OUTD / "random_modules_894.tsv"
RNG = np.random.default_rng(20261118)

prog = pd.read_csv(PROG, sep="\t", index_col=0)
df = prog[["direction"]].copy()
print(f"894 程序：UP {int((df.direction=='UP').sum())} / DOWN {int((df.direction=='DOWN').sum())}，共 {len(df)}")

expr = ec_expression().dropna()
dec = pd.qcut(expr.rank(method="first"), 10, labels=False)
prog_expr = expr.reindex(df.index).dropna()
need = prog_expr.groupby(dec.reindex(prog_expr.index)).size().to_dict()
prog_set = set(df.index)
pool = {}
for g in expr.index:
    if g in prog_set:
        continue
    pool.setdefault(int(dec.get(g, -1)), []).append(g)
recs = []
for rep in range(100):
    picks = []
    for d, k in need.items():
        cand = pool.get(int(d), [])
        if len(cand) < k:
            continue
        picks += list(RNG.choice(cand, size=k, replace=False))
    for g in picks:
        recs.append({"module": f"R{rep:03d}", "gene": g})
rm = pd.DataFrame(recs).drop_duplicates()
sizes = rm.groupby("module").size()
bad = sizes[sizes != len(df)]
if len(bad):
    raise SystemExit(f"随机模块尺寸不符：{len(bad)} 套不是 {len(df)} 基因")
rm.to_csv(OUT, sep="\t", index=False)
print(f"随机模块 {rm.module.nunique()} 套 × {int(sizes.median())} 基因 → {OUT}")
