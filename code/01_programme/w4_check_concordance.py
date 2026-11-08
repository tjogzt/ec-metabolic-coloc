
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
"""复核 v2 程序在三个外部队列中的方向一致率（附录 J 用数）。"""
import pandas as pd, numpy as np
from scipy.stats import binomtest

W4 = str(CFG.ROOT)
prog = pd.read_csv(f"{W4}/results/IR_state_program_v2.tsv", sep="\t", index_col=0)
print("程序:", len(prog), "列:", list(prog.columns))

def conc(file, label, fc_col=None):
    d = pd.read_csv(file, sep="\t", index_col=0)
    cand = fc_col or [c for c in d.columns if "logFC" in c][0]
    d = d[[cand]].dropna()
    common = prog.index.intersection(d.index)
    same = int((np.sign(d.loc[common, cand]) == np.sign(prog.loc[common, "logFC_adipose_obese"])).sum())
    n = len(common)
    p = binomtest(same, n, 0.5, alternative="greater").pvalue
    print(f"{label:10s} 重叠 {n:4d}  同向 {same:4d} ({100*same/n:.1f}%)  二项 p = {p:.2e}")

conc(f"{W4}/DE_GSE18732_muscle_DM_state.tsv", "GSE18732")
conc(f"{W4}/DE_GSE15773_adipose_IR_state.tsv", "GSE15773")
conc(f"{W4}/DE_GSE13070_muscle_baseline_IR.tsv", "GSE13070")
