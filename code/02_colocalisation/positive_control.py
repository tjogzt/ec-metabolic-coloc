"""
positive_control.py — 阳性对照：用完全相同的 coloc 代码，检验已知有共定位的
基因-性状组合，证明管线不是“什么都测不出来”。

核心对照（论文主图）：
  同一基因 SLC30A8、同一组织 胰腺、同一段区域
    T2D  → 期望 PP.H4 高（SLC30A8 是 β 细胞锌转运体，经典 T2D 基因）
    EC   → 实测 PP.H3 = 0.958（子宫）/ 低 H4
若 T2D 也测不出 H4，说明管线坏了，EC 的阴性不可采信。

用法：python positive_control.py --tissue Pancreas --gwas ../data/finngen_R12_T2D.gz
"""
from __future__ import annotations
import argparse, pathlib, sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from io_finngen import FinnGen
from coloc_abf import coloc_abf, log_abf
from coloc_eqtlcat import qtl_region
from coloc_loci import TISSUE_QTD, BASE, DATA

CANONICAL = {           # 经典 T2D 基因（ENSG，无版本号）
    "SLC30A8":   "ENSG00000164756",
    "KCNJ11":    "ENSG00000187486",
    "ABCC8":     "ENSG00000006071",
    "GCK":       "ENSG00000106633",
    "HNF1A":     "ENSG00000135100",
    "HNF1B":     "ENSG00000275410",
    "TCF7L2":    "ENSG00000148737",
    "PPARG":     "ENSG00000132170",
    "KCNQ1":     "ENSG00000053918",
    "IGF2BP2":   "ENSG00000073756",
    "MTNR1B":    "ENSG00000134640",
    "CDKN2B":    "ENSG00000147883",
}


BMI_GENES = {           # 肥胖/BMI 经典位点（多为调控型变异，适合做阳性对照）
    "FTO":    "ENSG00000140718",
    "MC4R":   "ENSG00000166603",
    "SH2B1":  "ENSG00000178188",
    "TMEM18": "ENSG00000151353",
    "BDNF":   "ENSG00000176697",
    "PCSK1":  "ENSG00000175426",
    "GCKR":   "ENSG00000084734",
    "PPARG":  "ENSG00000132170",
    "LEPR":   "ENSG00000116678",
    "IRX3":   "ENSG00000166444",
    "IRX5":   "ENSG00000176842",
}


def run(gwas_path: str, tissue: str, genes: dict[str, str], min_snps=10) -> pd.DataFrame:
    qtd = TISSUE_QTD[tissue]
    local = DATA / f"{qtd}.all.tsv.gz"
    url = f"{BASE}/{qtd}/{qtd}.all.tsv.gz"
    head = __import__("subprocess").run(
        f"curl -sI -m 60 {url} | grep -i '^content-length'",
        shell=True, capture_output=True, text=True).stdout
    exp = int(head.split(":")[-1].strip()) if head.strip() else 0
    src = local if (local.exists() and exp and local.stat().st_size == exp) else url
    print(f"数据源：{src if isinstance(src, str) else src.name}")

    idx = pd.read_csv(DATA.parent / "pipeline" / "results" / "gene_index.tsv", sep="\t")
    fg = FinnGen(gwas_path)
    out = []
    for sym, g in genes.items():
        row = idx[idx.gene_id == g]
        if row.empty:
            print(f"  {sym:8s} 不在基因索引中，跳过")
            continue
        r = row.iloc[0]
        q = qtl_region(r.chrom, int(r.lo), int(r.hi), src)
        q = q[q.gene_id == g].dropna(subset=["beta", "se", "pvalue"])
        q = q[(q.maf >= 0.01) & (q.maf <= 0.99) & (q.se > 0)]
        if len(q) < min_snps:
            print(f"  {sym:8s} eQTL 变异不足（{len(q)}）")
            continue
        gw = fg.region(r.chrom, int(q.position.min()), int(q.position.max()))
        if gw.empty:
            continue
        lut = {}
        for x in gw.itertuples(index=False):
            lut[f"{x.pos}_{x.ref}_{x.alt}"] = (x.beta, x.sebeta, x.pval)
            lut[f"{x.pos}_{x.alt}_{x.ref}"] = (-x.beta, x.sebeta, x.pval)
        rows = []
        for x in q.itertuples(index=False):
            k = f"{int(x.position)}_{x.ref}_{x.alt}"
            if k in lut:
                b, se, p = lut[k]
                if np.isfinite(b) and np.isfinite(se) and se > 0:
                    rows.append((x.beta, x.se, b, se, p))
        if len(rows) < min_snps:
            print(f"  {sym:8s} 可配对变异不足（{len(rows)}）")
            continue
        d = pd.DataFrame(rows, columns=["b_e", "se_e", "b_g", "se_g", "p_g"])
        res = coloc_abf(log_abf(d.b_e.values, d.se_e.values, 1.0),
                        log_abf(d.b_g.values, d.se_g.values, 1.0))
        res.update({"gene": sym, "gene_id": g, "tissue": tissue, "n_snps": len(d),
                    "best_gwas_p": float(d.p_g.min())})
        out.append(res)
        print(f"  {sym:8s} n={len(d):5d}  H4={res['PP.H4']:.3f}  H3={res['PP.H3']:.3f}  "
              f"最强 GWAS p={d.p_g.min():.1e}")
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", default="Pancreas")
    ap.add_argument("--gwas", default=str(DATA / "finngen_R12_T2D.gz"))
    ap.add_argument("--out", default="results/positive_control_t2d.tsv")
    ap.add_argument("--set", default="t2d", choices=["t2d", "bmi"])
    a = ap.parse_args()
    df = run(a.gwas, a.tissue, BMI_GENES if a.set == "bmi" else CANONICAL)
    if len(df):
        p = pathlib.Path(a.out); p.parent.mkdir(parents=True, exist_ok=True)
        df.sort_values("PP.H4", ascending=False).to_csv(p, sep="\t", index=False)
        print(f"\n写表 {p}：{len(df)} 个基因；PP.H4>0.8 有 {int((df['PP.H4']>0.8).sum())} 个")


if __name__ == "__main__":
    main()
