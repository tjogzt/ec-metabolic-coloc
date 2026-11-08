
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
import pandas as pd
RES = str(CFG.RESULTS)
new = pd.read_csv(f"{RES}/coloc_omara_Uterus.tsv", sep="\t")
old = pd.read_csv(f"{RES}/coloc_loci_Uterus.tsv", sep="\t")

print("new O'Mara : n=%d  H4max=%.3f  H4>0.5=%d  H3>0.8=%d  gwas_p_med=%.2e  n_snps_med=%d"
      % (len(new), new["PP.H4"].max(), (new["PP.H4"] > 0.5).sum(),
         (new["PP.H3"] > 0.8).sum(), new.best_gwas_p.median(), new.n_snps.median()))
print("old FinnGen: n=%d  H4max=%.3f  H4>0.5=%d  H3>0.8=%d  gwas_p_med=%.2e"
      % (len(old), old["PP.H4"].max(), (old["PP.H4"] > 0.5).sum(),
         (old["PP.H3"] > 0.8).sum(), old.best_gwas_p.median()))
print()
print("top H4 (new):")
print(new.nlargest(5, "PP.H4")[["locus", "PP.H3", "PP.H4", "n_snps", "best_gwas_p"]].to_string(index=False))
