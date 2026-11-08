
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
"""验证 SKAP1(liver) 的 PP.H4=0.922 是否来自两个峰真的重合。"""
import sys, pathlib
import numpy as np, pandas as pd

sys.path.insert(0, str(CFG.CODE / "02_colocalisation"))
from liftover import Chain
from p0_coloc_omara import load_gwas, lookup

DATA = CFG.DATA
CH_DN = DATA / "hg38ToHg19.over.chain.gz"

q = pd.read_csv(DATA / "_tmp_QTD000266_loci.tsv", sep="\t")
s = q[q.gene_id == "ENSG00000141293"].copy()
print("SKAP1 记录", len(s))

dn = Chain(CH_DN)
s["pos19"] = [dn.lift("17", int(p))[1] if dn.lift("17", int(p)) else np.nan
              for p in s.position]
s = s.dropna(subset=["pos19"])
s["pos19"] = s.pos19.astype(int)

store = load_gwas()
rows = []
for x in s.itertuples(index=False):
    hit = lookup(store, "17", int(x.pos19), str(x.ref).upper(), str(x.alt).upper())
    if hit is None:
        continue
    b, se, p = hit
    rows.append((int(x.pos19), x.pvalue, x.beta / x.se, p, b / se))
d = pd.DataFrame(rows, columns=["pos19", "eqtl_p", "eqtl_z", "gwas_p", "gwas_z"])
print(f"可配对 {len(d)} 个变异")

e_lead = d.loc[d.eqtl_p.idxmin()]
g_lead = d.loc[d.gwas_p.idxmin()]
print(f"\neQTL lead : pos19={e_lead.pos19:,}  z={e_lead.eqtl_z:+.2f}  (p={e_lead.eqtl_p:.2e})")
print(f"GWAS lead : pos19={g_lead.pos19:,}  z={g_lead.gwas_z:+.2f}  (p={g_lead.gwas_p:.2e})")
print(f"两峰间距  : {abs(int(e_lead.pos19)-int(g_lead.pos19))/1000:.1f} kb")
print(f"eQTL lead 处的 GWAS: z={e_lead.gwas_z:+.2f} (p={e_lead.gwas_p:.2e})")
from scipy.stats import spearmanr
r, p = spearmanr(d.eqtl_z, d.gwas_z)
print(f"全窗口 z 相关: rho={r:+.3f} (p={p:.1e})")
