
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
"""主席复核：51.9% 的方向一致率是否被无信号基因稀释。
限定到钳夹队列中有响应的基因，方向一致率是否显著高于 50%。
口径与 w4_check_concordance.py 完全一致（v2 程序 vs GSE13070 DE）。
"""
import pandas as pd, numpy as np
from scipy.stats import binomtest, spearmanr

W4 = str(CFG.ROOT)
prog = pd.read_csv(f"{W4}/results/IR_state_program_v2.tsv", sep="\t", index_col=0)
de = pd.read_csv(f"{W4}/DE_GSE13070_muscle_baseline_IR.tsv", sep="\t", index_col=0)
fc = [c for c in de.columns if "logFC" in c][0]
fdr = [c for c in de.columns if "adj" in c.lower() or "FDR" in c][0]
pcol = [c for c in de.columns if c.lower() == "p.value" or c.lower() == "pvalue"][0]

m = prog[["logFC_adipose_obese"]].join(de[[fc, fdr, pcol]], how="inner").dropna()
m["same"] = np.sign(m["logFC_adipose_obese"]) == np.sign(m[fc])
print(f"重叠基因 {len(m)}；全口径同向 {m.same.sum()} ({100*m.same.mean():.2f}%)；"
      f"单侧 p={binomtest(int(m.same.sum()), len(m), 0.5, alternative='greater').pvalue:.4f}")

for label, sub in [
    ("FDR<0.05", m[m[fdr] < 0.05]),
    ("名义 p<0.001", m[m[pcol] < 1e-3]),
    ("|logFC|>=0.3", m[m[fc].abs() >= 0.3]),
    ("FDR<0.05 且 |logFC|>=0.3", m[(m[fdr] < 0.05) & (m[fc].abs() >= 0.3)]),
]:
    k, n = int(sub.same.sum()), len(sub)
    if n == 0:
        print(f"{label:26s} n=0"); continue
    p1 = binomtest(k, n, 0.5, alternative="greater").pvalue
    print(f"{label:26s} {k}/{n} = {100*k/n:5.1f}%   单侧 p={p1:.4g}")

rho, pr = spearmanr(m["logFC_adipose_obese"], m[fc])
print(f"\n全 862 基因 logFC 秩相关 Spearman rho={rho:+.3f}, p={pr:.4g}")
