"""
sensitivity_r1m3.py — R1-M3 敏感性：肝 SKAP1 17q21.32 单因果变异假设检验

背景（nature-reviewer R1-M3）：主分析 PP.H4=0.93 基于"每性状至多一个因果变异"模型；
GTEx v11 SuSiE 显示肝 SKAP1 有两个可信集（CS1 53 变异、CS2 60 变异）。
本脚本做无 LD 的条件化近似敏感性：
  1) 基线复现（全窗口，p12=1e-5，须精确等于 0.930889 以验证数据通路）
  2) CS1-only / CS2-only / CS1∪CS2-only
  3) 窗口去 CS1 / 窗口去 CS2
  4) 剔除 top 变异（48,127,474）/ 剔除 v8 lead（48,127,316）
  5) 仅 top 变异单点

数据源（远端 tabix，均在 2026-11-20 提取）：
  eQTL: https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000241/QTD000241.all.tsv.gz
        （GTEx v8 liver，beta 相对 ALT）
  GWAS: https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90454001-GCST90455000/
        GCST90454186/harmonised/GCST90454186.h.tsv.gz（beta 相对 effect_allele）
匹配语义与 coloc_robustness.py 完全一致：位置+ref/alt（大写），翻转取 -beta。
"""
from __future__ import annotations
import io, math, subprocess, tarfile, sys, time
import numpy as np
import pandas as pd

sys.path.insert(0, "/Users/taozhu/my researches/gongbing/06_EC_IR/pipeline")
from coloc_abf import log_abf
from coloc_robustness import posteriors  # noqa: E402  (同口径五假设后验)

EQTL_URL = ("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/"
            "QTD000266/QTD000266.all.tsv.gz")
EQTL_CACHE = ("/Users/taozhu/my researches/gongbing/06_EC_IR/data/"
              "_tmp_ramach_QTD000266_loci.tsv")   # 原始分析用缓存，若在盘直接复用
GWAS_URL = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
            "GCST90454001-GCST90455000/GCST90454186/harmonised/GCST90454186.h.tsv.gz")
GENE = "ENSG00000141293"           # SKAP1
LEAD19 = 48_127_474
PAD = 600_000
LO, HI = LEAD19 - PAD, LEAD19 + PAD
P12 = 1e-5
TAR = ("/Users/taozhu/my researches/gongbing/06_EC_IR/data/"
       "GTEx_Analysis_v11_eQTL_SuSiE.tar")
OUT = ("/Users/taozhu/my researches/gongbing/06_EC_IR/pipeline/results/"
       "r1m3_liver_skap1_sensitivity.tsv")


def tabix(url: str, region: str) -> str:
    out = subprocess.run(["tabix", url, region], capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(f"tabix failed {url} {region}: {out.stderr[:200]}")
    return out.stdout


def eqtl_window() -> pd.DataFrame:
    import os
    if os.path.exists(EQTL_CACHE):
        d = pd.read_csv(EQTL_CACHE, sep="\t")
        d = d[["gene_id", "chromosome", "position", "ref", "alt", "maf", "beta", "se"]]
        d.columns = ["gene_id", "chrom", "pos", "ref", "alt", "maf", "beta", "se"]
        d = d[(d.gene_id == GENE) & (d.maf.between(0.01, 0.99)) & (d.se > 0)
              & d.beta.notna() & d.se.notna()]
        d["ref"] = d.ref.str.upper()
        d["alt"] = d.alt.str.upper()
        return d.reset_index(drop=True)
    txt = tabix(EQTL_URL, f"17:{LO}-{HI}")
    rows = [l.split("\t") for l in txt.strip().split("\n") if l]
    d = pd.DataFrame(rows).iloc[:, [0, 1, 2, 3, 4, 7, 9, 10]]
    d.columns = ["gene_id", "chrom", "pos", "ref", "alt", "maf", "beta", "se"]
    for c in ("pos", "maf", "beta", "se"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d[(d.gene_id == GENE) & (d.maf.between(0.01, 0.99)) & (d.se > 0)
          & d.beta.notna() & d.se.notna()]
    d["ref"] = d.ref.str.upper()
    d["alt"] = d.alt.str.upper()
    return d.reset_index(drop=True)


def gwas_window() -> pd.DataFrame:
    txt = tabix(GWAS_URL, f"17:{LO}-{HI}")
    rows = [l.split("\t") for l in txt.strip().split("\n") if l]
    # 列(0-based): 0 chrom,1 pos,2 effect_allele,3 other_allele,4 beta,5 se,7 p_value
    d = pd.DataFrame(rows).iloc[:, [0, 1, 2, 3, 4, 5, 7]]
    d.columns = ["chrom", "pos", "effect", "other", "beta", "se", "pval"]
    for c in ("pos", "beta", "se", "pval"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d[(d.effect.str.len() == 1) & (d.other.str.len() == 1)
          & d.beta.notna() & d.se.notna() & (d.se > 0)]
    d["effect"] = d.effect.str.upper()
    d["other"] = d.other.str.upper()
    d = d.drop_duplicates(subset=["pos", "effect", "other"], keep="first")
    return d.reset_index(drop=True)


def v11_cs() -> dict[int, set[str]]:
    """GTEx v11 肝 SKAP1 可信集成员：cs_id -> {chr17_POS_REF_ALT_b38}"""
    tar = tarfile.open(TAR)
    liv = [n for n in tar.getnames() if "/Liver." in n][0]
    f = tar.extractfile(liv)
    p = pd.read_parquet(io.BytesIO(f.read()))
    s = p[p.gene_name == "SKAP1"]
    return {int(c): set(g.variant_id) for c, g in s.groupby("cs_id")}


def match_gwas(dg: pd.DataFrame, pos: int, ref: str, alt: str) -> tuple | None:
    """与 coloc_robustness.lookup 同语义：位置 + 等位基因，先直向再翻转。"""
    sub = dg[dg.pos.eq(pos)]
    m = sub[(sub.effect == alt) & (sub.other == ref)]
    if len(m):
        r = m.iloc[0]
        return float(r.beta), float(r.se)
    m = sub[(sub.effect == ref) & (sub.other == alt)]
    if len(m):
        r = m.iloc[0]
        return -float(r.beta), float(r.se)
    return None


def run(name: str, e: pd.DataFrame, dg: pd.DataFrame) -> dict:
    l1, l2, s1, s2, poss = [], [], [], [], []
    for x in e.itertuples(index=False):
        hit = match_gwas(dg, int(x.pos), str(x.ref), str(x.alt))
        if hit is None:
            continue
        b, s = hit
        l1.append(float(x.beta))
        s1.append(float(x.se))
        l2.append(b)
        s2.append(s)
        poss.append(int(x.pos))
    a1 = log_abf(np.array(l1), np.array(s1), 1.0)
    a2 = log_abf(np.array(l2), np.array(s2), 1.0)
    po = posteriors(a1, a2, P12)
    return {"run": name, "n": int(len(l1)),
            "H4": float(po["H4"]), "H3": float(po["H3"]),
            "cond_H4": float(po["H4"] / (po["H3"] + po["H4"])),
            "npos": poss}


def main() -> None:
    t0 = time.time()
    print("fetching eQTL window ...", flush=True)
    e = eqtl_window()
    print(f"eQTL SKAP1 rows in window: {len(e)}", flush=True)
    print("fetching GWAS window ...", flush=True)
    dg = gwas_window()
    print(f"GWAS rows in window: {len(dg)}", flush=True)
    cs = v11_cs()
    print(f"v11 CS sizes: { {k: len(v) for k, v in cs.items()} }", flush=True)

    def cs_mask(e: pd.DataFrame, cid: int) -> pd.DataFrame:
        ids = cs[cid]
        return e[[(f"chr17_{x.pos}_{x.ref}_{x.alt}_b38" in ids)
                  for x in e.itertuples()]].reset_index(drop=True)

    e1, e2 = cs_mask(e, 1), cs_mask(e, 2)
    rows = [run("baseline_full_window", e, dg),
            run("CS1_only", e1, dg),
            run("CS2_only", e2, dg),
            run("CS1_union_CS2", pd.concat([e1, e2]).drop_duplicates("pos"), dg),
            run("window_minus_CS1", e.drop(e1.index), dg),
            run("window_minus_CS2", e.drop(e2.index), dg),
            run("drop_top_48127474", e[e.pos != 48_127_474], dg),
            run("drop_v8lead_48127316", e[e.pos != 48_127_316], dg)]
    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep="\t", index=False)
    print(out.to_string(index=False))
    print(f"\n基线对照：robustness 表 H4=0.930889 cond=0.935 n=2805（p12=1e-5）")
    print(f"耗时 {time.time()-t0:.0f}s → {OUT}")


if __name__ == "__main__":
    main()
