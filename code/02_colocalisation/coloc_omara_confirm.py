"""
coloc_omara_confirm.py — 用预注册指定的结局 GCST006465（O'Mara 2018）重跑共定位层。

定位：这不是独立复制。GCST006465 与新用的 GCST90454186 共享 ECAC 核心样本，
两套结果的一致性是"换结局数据集后阳性是否还在"的**一致性检验**，措辞不得写成 replicated。

与 p0_coloc_omara.py 的差别：
  - 不用 tabix（QTD000266 等文件历来无 .tbi，tabix 会静默返回空）
  - 不做逐变异坐标抬升，改用 rsID 配对 + 等位基因定向
  - 复用已缓存的 hg38 区域文件 _tmp_ramach_<QTD>_loci.tsv
"""
from __future__ import annotations
import gzip, pathlib, sys, time
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
ROOT = HERE.parent
DATA = ROOT / "data"
RES = HERE / "results"

from liftover import Chain
from coloc_abf import coloc_abf, log_abf
from coloc_loci import tissue_egenes_in_loci, TISSUE_QTD
from p0_coloc_omara import peaks

PAD = 600_000
MIN_SNPS = 10


def load_omara_by_rsid(path: pathlib.Path) -> pd.DataFrame:
    """GCST006465 → 以 rsID 为键的轻量表（effect/other 均大写）。"""
    rows = []
    with gzip.open(path, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(hdr)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            try:
                rs = p[ix["variant_id"]]
                b = float(p[ix["beta"]]); se = float(p[ix["standard_error"]])
                pv = float(p[ix["p_value"]])
            except (ValueError, IndexError, KeyError):
                continue
            if not rs.startswith("rs") or not np.isfinite(b) or not np.isfinite(se) or se <= 0:
                continue
            rows.append((rs, p[ix["other_allele"]].upper(), p[ix["effect_allele"]].upper(),
                         b, se, pv))
    df = pd.DataFrame(rows, columns=["rsid", "other", "effect", "g_beta", "g_se", "g_p"])
    df = df.drop_duplicates("rsid", keep="first")
    return df


def main():
    t0 = time.time()
    print("载入 GCST006465…", flush=True)
    gw = load_omara_by_rsid(DATA / "GCST006465_GRCh37.tsv.gz")
    print(f"  {len(gw):,} 个 rsID", flush=True)

    # 位点：按 p<5e-8 + 1Mb 间隔聚类（GRCh37），lead 抬到 hg38 以定位 eQTL 窗口
    recs = []
    with gzip.open(DATA / "GCST006465_GRCh37.tsv.gz", "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(hdr)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            try:
                c = str(p[ix["chromosome"]]).replace("chr", "")
                pos = int(p[ix["base_pair_location"]])
                pv = float(p[ix["p_value"]])
            except (ValueError, IndexError):
                continue
            recs.append((c, pos, pv))
    g37 = pd.DataFrame(recs, columns=["chrom", "pos", "pval"])
    # 简化版聚类：p<5e-8，同染色体间隔 1Mb 内归为同一座
    g37 = g37[g37.pval < 5e-8].sort_values(["chrom", "pos"])
    loci = []
    for c, d in g37.groupby("chrom"):
        d = d.sort_values("pos")
        cur = None
        for x in d.itertuples(index=False):
            if cur is None or x.pos - cur["end"] > 1_000_000:
                if cur:
                    loci.append(cur)
                cur = {"chrom": c, "start": x.pos, "end": x.pos,
                       "lead": x.pos, "p": x.pval}
            else:
                cur["end"] = x.pos
                if x.pval < cur["p"]:
                    cur["p"] = x.pval; cur["lead"] = x.pos
        if cur:
            loci.append(cur)
    loc = pd.DataFrame(loci)
    print(f"  O'Mara 位点 {len(loc)} 个", flush=True)

    up = Chain(str(DATA / "hg19ToHg38.over.chain.gz"))
    lifted = []
    for r in loc.itertuples(index=False):
        m = up.lift(str(r.chrom), int(r.lead))
        if m is None:
            continue
        lifted.append({"chrom": str(m[0]), "start": int(m[1]), "end": int(m[1]),
                       "lead": int(m[1]), "p": r.p, "lead19": int(r.lead)})
    loc = pd.DataFrame(lifted)
    print(f"  可抬升至 hg38：{len(loc)} 个", flush=True)

    all_out = []
    for tissue, qtd in TISSUE_QTD.items():
        tmp = DATA / f"_tmp_ramach_{qtd}_loci.tsv"
        if not tmp.exists():
            print(f"  跳过 {tissue}（无缓存）", flush=True)
            continue
        tt = time.time()
        eg = tissue_egenes_in_loci(tissue, loc, pad=PAD)
        q = pd.read_csv(tmp, sep="\t")
        q = q.dropna(subset=["beta", "se"])
        q = q[(q.maf >= 0.01) & (q.maf <= 0.99) & (q.se > 0)]
        q = q.assign(pos19=q["position"].astype(int))
        n_here = 0
        for i, r in loc.iterrows():
            genes = eg.get(i, [])
            if not genes:
                continue
            near = q[(q.chromosome.astype(str) == str(r.chrom)) &
                     (q.pos19.between(r.lead - PAD, r.lead + PAD))]
            if near.empty:
                continue
            for g in genes:
                sub = near[near.gene_id == g]
                if len(sub) < MIN_SNPS:
                    continue
                m = sub.merge(gw, on="rsid", how="inner")
                if len(m) < MIN_SNPS:
                    continue
                ref = m.ref.astype(str).str.upper()
                alt = m.alt.astype(str).str.upper()
                oth = m.other.astype(str).str.upper()
                eff = m.effect.astype(str).str.upper()
                fwd = (ref == oth) & (alt == eff)
                rev = (ref == eff) & (alt == oth)
                keep = fwd | rev
                if keep.sum() < MIN_SNPS:
                    continue
                b_g = np.where(fwd[keep].values, m.g_beta[keep].values, -m.g_beta[keep].values)
                res = coloc_abf(log_abf(m.beta[keep].values, m.se[keep].values, 1.0),
                                log_abf(b_g, m.g_se[keep].values, 1.0))
                res = {k: v for k, v in res.items() if k not in ("n_snps", "top_snp_share")}
                all_out.append(dict(**res, n_snps=int(keep.sum()), top_snp_share=np.nan,
                                    tissue=tissue, locus=f"{r.chrom}:{r.lead19}",
                                    gene_id=g, best_gwas_p=float(m.g_p[keep].min()),
                                    alleles_matched=int(keep.sum()),
                                    alleles_dropped=int(len(m) - keep.sum())))
                n_here += 1
        print(f"  {tissue}: {n_here} 次检验（{time.time()-tt:.0f}s）", flush=True)

    out = pd.DataFrame(all_out)
    if out.empty:
        print("无结果", flush=True); return
    cols = ["PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "n_snps", "tissue", "locus",
            "gene_id", "best_gwas_p", "alleles_matched", "alleles_dropped"]
    out = out[cols]
    f = RES / "coloc_omara_confirm.tsv"
    out.to_csv(f, sep="\t", index=False)
    print(f"\n{len(out)} 次检验 → {f.name}", flush=True)
    print(f"  PP.H4>0.8: {(out['PP.H4']>0.8).sum()}  >0.5: {(out['PP.H4']>0.5).sum()}"
          f"  最大 {out['PP.H4'].max():.3f}  PP.H3>0.8: {(out['PP.H3']>0.8).sum()}", flush=True)
    print("\n=== 超过 0.8 的命中 ===", flush=True)
    for _, r in out[out["PP.H4"] > 0.8].iterrows():
        print(f"  {r.tissue:24s} {r.locus:16s} {r.gene_id:20s} "
              f"H4={r['PP.H4']:.3f} n={r.n_snps}", flush=True)
    print(f"\n总耗时 {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
