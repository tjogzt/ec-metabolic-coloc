"""
p0_coloc_omara2.py — P0-1（改进版）：不依赖 tabix 索引。

问题：QTD000121/266/281/296 四个 eQTL 文件没有 .tbi 索引，tabix 静默返回空，
导致 coloc 结果为 0 条却不报错。

做法：每个组织只扫描整个文件一次，按 11 个位点区间（±pad）过滤出小表，
其余步骤与原 coloc_loci 完全相同（ABF + 等位基因对齐）。
"""
from __future__ import annotations
import argparse, gzip, pathlib, sys, time
from bisect import bisect_right
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
ROOT = HERE.parent
DATA = ROOT / "data"
RES = HERE / "results"

from liftover import Chain
from coloc_abf import coloc_abf, log_abf
from coloc_loci import tissue_egenes_in_loci, TISSUE_QTD, BASE
from p0_coloc_omara import load_gwas, peaks, lookup

GWAS = DATA / "GCST90454186.h.tsv.gz"
CH_UP = DATA / "hg19ToHg38.over.chain.gz"
CH_DN = DATA / "hg38ToHg19.over.chain.gz"

COLS = ["gene_id", "chromosome", "position", "ref", "alt", "variant", "maf",
        "pvalue", "beta", "se", "median_tpm", "rsid"]
KEEP = [16, 1, 2, 3, 4, 5, 7, 8, 9, 10, 17, 18]   # 0-based


def fetch_regions_remote(url: str, regions: dict, out_path) -> int:
    """用远端 tabix 按区间取数（本地文件截断时的替代路径）。"""
    import subprocess
    kept = 0
    with open(out_path, "w") as o:
        o.write("\t".join(COLS) + "\n")
        for c, rngs in regions.items():
            for lo, hi in rngs:
                out = subprocess.run(["tabix", url, f"{c}:{lo}-{hi}"],
                                     capture_output=True, text=True)
                for line in out.stdout.strip().split("\n"):
                    if not line:
                        continue
                    parts = line.split("\t")
                    o.write("\t".join(parts[j] for j in KEEP) + "\n")
                    kept += 1
    return kept


def stream_regions(path, regions: dict, out_path) -> int:
    """单次扫描，保留落在任一位点区间的行。返回保留行数。"""
    starts = {c: [r[0] for r in v] for c, v in regions.items()}
    kept = 0
    with gzip.open(path, "rt") as f, open(out_path, "w") as o:
        o.write("\t".join(COLS) + "\n")
        f.readline()  # header
        for line in f:
            p = line.split("\t", 4)
            c = p[1]
            rng = regions.get(c)
            if rng is None:
                continue
            try:
                pos = int(p[2])
            except ValueError:
                continue
            i = bisect_right(starts[c], pos) - 1
            if i >= 0 and pos <= rng[i][1]:
                parts = line.rstrip("\n").split("\t")
                o.write("\t".join(parts[j] for j in KEEP) + "\n")
                kept += 1
    return kept


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", required=True)
    ap.add_argument("--pad", type=int, default=600_000)
    ap.add_argument("--min_snps", type=int, default=10)
    ap.add_argument("--p_thresh", type=float, default=5e-8)
    ap.add_argument("--gap", type=int, default=1_000_000)
    ap.add_argument("--remote", action="store_true", help="用 EBI 远端 tabix 取区间")
    a = ap.parse_args()
    t0 = time.time()

    store = load_gwas(GWAS)
    print(f"GWAS 载入完成：{sum(len(v) for v in store.values()):,} 变异", flush=True)

    loc = peaks(store, a.p_thresh, a.gap)
    # 全程 hg38：新 GWAS 与 eQTL Catalogue 同坐标，无需 liftOver
    loc = loc.rename(columns={"lead": "lead38"})
    loc["lead"] = loc["lead38"]
    loc["lead19"] = loc["lead38"]     # 下游字段名沿用，实际为 hg38
    print(f"位点 {len(loc)} 个（p<{a.p_thresh:g}）", flush=True)

    eg = tissue_egenes_in_loci(a.tissue, loc, pad=a.pad)
    print("候选 eGene 数：" + str({i: len(v) for i, v in eg.items()}), flush=True)

    regions: dict[str, list] = {}
    for r in loc.itertuples(index=False):
        regions.setdefault(str(r.chrom), []).append(
            (max(1, r.lead - a.pad), r.lead + a.pad))
    for c in regions:
        regions[c].sort()

    qtd = TISSUE_QTD[a.tissue]
    src = DATA / f"{qtd}.all.tsv.gz"
    tmp = DATA / f"_tmp_ramach_{qtd}_loci.tsv"   # 独立命名，避免复用旧位点窗口
    if not tmp.exists():
        if a.remote:
            url = f"{BASE}/{qtd}/{qtd}.all.tsv.gz"
            n = fetch_regions_remote(url, regions, tmp)
            print(f"远端取回 {n:,} 行 → {tmp.name}（{time.time()-t0:.0f}s）", flush=True)
        else:
            n = stream_regions(src, regions, tmp)
            print(f"扫描完成：保留 {n:,} 行 → {tmp.name}（{time.time()-t0:.0f}s）", flush=True)
    q = pd.read_csv(tmp, sep="\t")
    print(f"区域内 eQTL 记录 {len(q):,}", flush=True)

    q = q.dropna(subset=["beta", "se"])
    q = q[(q.maf >= 0.01) & (q.maf <= 0.99) & (q.se > 0)]
    q = q.assign(pos19=q["position"].astype(int))   # 全程 hg38
    print(f"可降坐标 {len(q):,}", flush=True)

    recs = []
    for i, r in loc.iterrows():
        genes = eg.get(i, [])
        if not genes:
            continue
        near = q[(q.chromosome.astype(str) == str(r.chrom)) &
                 (q.pos19.between(r.lead19 - a.pad, r.lead19 + a.pad))]
        if near.empty:
            continue
        for g in genes:
            sub = near[near.gene_id == g]
            if len(sub) < a.min_snps:
                continue
            rows = []
            for x in sub.itertuples(index=False):
                hit = lookup(store, str(r.chrom), int(x.pos19),
                             str(x.ref).upper(), str(x.alt).upper())
                if hit is None:
                    continue
                b, se, p = hit
                if np.isfinite(b) and np.isfinite(se) and se > 0:
                    rows.append((x.beta, x.se, b, se, p))
            if len(rows) < a.min_snps:
                continue
            d = pd.DataFrame(rows, columns=["b_e", "se_e", "b_g", "se_g", "p_g"])
            res = coloc_abf(log_abf(d.b_e.values, d.se_e.values, 1.0),
                            log_abf(d.b_g.values, d.se_g.values, 1.0))
            res = {k: v for k, v in res.items() if k not in ("n_snps", "top_snp_share")}
            recs.append(dict(**res, n_snps=len(d), top_snp_share=np.nan,
                             tissue=a.tissue, locus=f"{r.chrom}:{r.lead19}",
                             gene_id=g, best_gwas_p=float(d.p_g.min())))

    if not recs:
        print("无结果（检查区间与基因匹配）", flush=True)
        return
    out = pd.DataFrame(recs)[["PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "n_snps",
                              "top_snp_share", "tissue", "locus", "gene_id", "best_gwas_p"]]
    f = RES / f"coloc_ramach_{a.tissue}.tsv"
    out.to_csv(f, sep="\t", index=False)
    print(f"\n{len(out)} 次检验 → {f.name}", flush=True)
    print(f"  H4>0.8={int((out['PP.H4']>0.8).sum())}  H4>0.5={int((out['PP.H4']>0.5).sum())}"
          f"  maxH4={out['PP.H4'].max():.3f}  H3>0.8={int((out['PP.H3']>0.8).sum())}", flush=True)
    print(f"总耗时 {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
