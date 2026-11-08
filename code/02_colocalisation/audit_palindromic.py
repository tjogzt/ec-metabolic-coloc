
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
"""审查用：检验等位基因对齐是否把回文位点的效应方向弄反。

逻辑：若两个数据集都在参考正链上，回文位点无需翻转，其 eQTL z 与 GWAS z 的相关性
应与非回文位点相同。若回文位点被错误翻转，其相关性会显著变号或归零。
"""
import sys, pathlib
import numpy as np, pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(CFG.CODE / "02_colocalisation"))
from liftover import Chain
from p0_coloc_omara import load_gwas, lookup

DATA = CFG.DATA
tmp = DATA / "_tmp_QTD000266_loci.tsv"
dn = Chain(DATA / "hg38ToHg19.over.chain.gz")
store = load_gwas()

q = pd.read_csv(tmp, sep="\t")
s = q[q.gene_id == "ENSG00000141293"].dropna(subset=["beta", "se"]).copy()
s = s[(s.maf >= 0.01) & (s.maf <= 0.99) & (s.se > 0)]
s["pos19"] = [dn.lift("17", int(p))[1] if dn.lift("17", int(p)) else np.nan for p in s.position]
s = s.dropna(subset=["pos19"]); s["pos19"] = s.pos19.astype(int)

rows = []
for x in s.itertuples(index=False):
    hit = lookup(store, "17", int(x.pos19), str(x.ref).upper(), str(x.alt).upper())
    if hit is None:
        continue
    b, se, p = hit
    rows.append((str(x.ref).upper(), str(x.alt).upper(), x.maf,
                 x.beta / x.se, b / se, int(x.pos19)))
d = pd.DataFrame(rows, columns=["ref", "alt", "maf", "z_e", "z_g", "pos"])

pal = ((d.ref == "A") & (d.alt == "T")) | ((d.ref == "T") & (d.alt == "A")) | \
      ((d.ref == "C") & (d.alt == "G")) | ((d.ref == "G") & (d.alt == "C"))
amb = pal & d.maf.between(0.05, 0.45)

print(f"SKAP1 窗口：配对 {len(d)} 个变异")
for lab, m in [("全部", np.ones(len(d), bool)), ("非回文", ~pal.values),
               ("回文(全部)", pal.values), ("回文且MAF 0.05-0.45", amb.values)]:
    sub = d[m]
    if len(sub) < 5:
        print(f"  {lab:22s} n={len(sub):5d} 样本不足"); continue
    rho, p = spearmanr(sub.z_e, sub.z_g)
    agree = float(((np.sign(sub.z_e) == np.sign(sub.z_g))).mean())
    print(f"  {lab:22s} n={len(sub):5d}  rho={rho:+.3f} (p={p:.1e})  同号比例={100*agree:.1f}%")

# 全窗口基线
rho, p = spearmanr(d.z_e, d.z_g)
print(f"\n参考：全窗口 rho={rho:+.3f}")
