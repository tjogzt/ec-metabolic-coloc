"""
coloc_robustness.py — 共定位层稳健性套件（回应 MoA 评审 D2）

复现 p0_coloc_ramach.py 的每一次检验，并额外输出：
  1. p12 先验敏感性：p12 ∈ {1e-6, 5e-6, 1e-5, 5e-5, 1e-4}（p1=p2=1e-4 固定）
  2. 条件后验 PP.H4/(PP.H3+PP.H4)：区分"共享"与"仅双向关联"
  3. H4 隐含的 95% 可信集：SNP 后验 ∝ ABF1×ABF2，按降序累加到 0.95，报告大小与是否含 lead
  4. 贝叶斯 FDR：阈值 PP.H4>0.8 下的期望假发现数 Σ(1-PP.H4)

注意：可信集由 coloc 的单因果假设导出，不等同于 SuSiE 可信集，输出列名与措辞须区分。
"""
from __future__ import annotations
import gzip, math, pathlib, sys, time
from bisect import bisect_right
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
ROOT = HERE.parent
DATA = ROOT / "data"
RES = HERE / "results"

from coloc_abf import log_abf
from coloc_loci import tissue_egenes_in_loci, TISSUE_QTD
from p0_coloc_omara import load_gwas, peaks, lookup

GWAS = DATA / "GCST90454186.h.tsv.gz"
PAD = 600_000
MIN_SNPS = 10
P_THRESH = 5e-8
GAP = 1_000_000

P1 = P2 = 1e-4
P12_GRID = [1e-6, 5e-6, 1e-5, 5e-5, 1e-4]


def _logsum(x: np.ndarray) -> float:
    m = np.max(x)
    if not np.isfinite(m):
        return -np.inf
    return float(m + np.log(np.sum(np.exp(x - m))))


def posteriors(l1: np.ndarray, l2: np.ndarray, p12: float,
               p1: float = P1, p2: float = P2) -> dict:
    """ABF 后验，与 coloc_abf 同口径但允许自由 p12。

    H3 的 off-diagonal 求和用恒等式避免 n×n 矩阵：
        Σ_{i≠j} exp(l1_i + l2_j) = [Σ_i exp(l1_i)][Σ_j exp(l2_j)] − Σ_i exp(l1_i + l2_i)
    """
    lh0 = math.log(max(1e-300, 1 - p1 - p2 - p12))
    L1, L2, L4 = _logsum(l1), _logsum(l2), _logsum(l1 + l2)
    lh1 = L1 + math.log(p1) + math.log(1 - p2)
    lh2 = L2 + math.log(p2) + math.log(1 - p1)
    # Δ = log(Σe^{l1+l2} / (Σe^{l1}·Σe^{l2})) ≤ 0，用 log1p 保持数值稳定
    delta = min(0.0, L4 - L1 - L2)
    lh3 = (math.log(p1) + math.log(p2) + L1 + L2 + math.log1p(-math.exp(delta))
           if delta < 0 else -math.inf)
    lh4 = L4 + math.log(p12)
    v = np.array([lh0, lh1, lh2, lh3, lh4])
    v = np.exp(v - _logsum(v))
    return {"H0": v[0], "H1": v[1], "H2": v[2], "H3": v[3], "H4": v[4]}


def credible_set(l1: np.ndarray, l2: np.ndarray, idx, thr: float = 0.95):
    """H4 隐含可信集：SNP 后验 ∝ exp(l1+l2)。"""
    w = l1 + l2
    w = np.exp(w - np.max(w))
    w = w / w.sum()
    order = np.argsort(-w)
    cum, members = 0.0, []
    for j in order:
        members.append(int(idx[j]))
        cum += float(w[j])
        if cum >= thr:
            break
    return members, float(w[order[0]]), int(order[0])


def main():
    t0 = time.time()
    store = load_gwas(GWAS)
    print(f"GWAS 载入：{sum(len(v) for v in store.values()):,} 变异", flush=True)
    loc = peaks(store, P_THRESH, GAP)
    loc["lead19"] = loc["lead"]
    print(f"位点 {len(loc)} 个", flush=True)

    rows = []
    for tissue, qtd in TISSUE_QTD.items():
        tmp = DATA / f"_tmp_ramach_{qtd}_loci.tsv"
        if not tmp.exists():
            print(f"跳过 {tissue}（无缓存 {tmp.name}）", flush=True)
            continue
        eg = tissue_egenes_in_loci(tissue, loc, pad=PAD)
        q = pd.read_csv(tmp, sep="\t")
        q = q.dropna(subset=["beta", "se"])
        q = q[(q.maf >= 0.01) & (q.maf <= 0.99) & (q.se > 0)]
        q = q.assign(pos19=q["position"].astype(int))
        n_done = 0
        for i, r in loc.iterrows():
            genes = eg.get(i, [])
            if not genes:
                continue
            near = q[(q.chromosome.astype(str) == str(r.chrom)) &
                     (q.pos19.between(r.lead19 - PAD, r.lead19 + PAD))]
            if near.empty:
                continue
            for g in genes:
                sub = near[near.gene_id == g]
                if len(sub) < MIN_SNPS:
                    continue
                be, se_e, bg, se_g, pos, rs = [], [], [], [], [], []
                for x in sub.itertuples(index=False):
                    hit = lookup(store, str(r.chrom), int(x.pos19),
                                 str(x.ref).upper(), str(x.alt).upper())
                    if hit is None:
                        continue
                    b, se, _p = hit
                    if np.isfinite(b) and np.isfinite(se) and se > 0:
                        be.append(x.beta); se_e.append(x.se)
                        bg.append(b); se_g.append(se)
                        pos.append(int(x.pos19)); rs.append(str(getattr(x, "rsid", "")))
                if len(be) < MIN_SNPS:
                    continue
                l1 = log_abf(np.array(be), np.array(se_e), 1.0)
                l2 = log_abf(np.array(bg), np.array(se_g), 1.0)
                mem, top_w, top_j = credible_set(l1, l2, range(len(be)))
                base = dict(tissue=tissue, locus=f"{r.chrom}:{r.lead19}", gene_id=g,
                            n_snps=len(be),
                            cred_set_size=len(mem),
                            cred_set_top_pos=pos[top_j],
                            cred_set_top_rsid=rs[top_j],
                            cred_set_top_post=round(top_w, 4),
                            cred_set_positions=";".join(str(pos[j]) for j in mem))
                for p12 in P12_GRID:
                    po = posteriors(l1, l2, p12)
                    rows.append(dict(base, p12=p12, PP_H3=po["H3"], PP_H4=po["H4"],
                                     cond_H4=(po["H4"] / (po["H3"] + po["H4"])
                                              if (po["H3"] + po["H4"]) > 0 else np.nan)))
                n_done += 1
        print(f"  {tissue}: {n_done} 次检验（{time.time()-t0:.0f}s）", flush=True)

    out = pd.DataFrame(rows)
    f = RES / "coloc_robustness.tsv"
    out.to_csv(f, sep="\t", index=False)
    print(f"\n{len(out)} 行 → {f.name}", flush=True)

    base = out[out.p12 == 1e-5]
    print(f"\n基线 p12=1e-5：检验 {len(base)} 次，H4>0.8 = {(base.PP_H4>0.8).sum()}, "
          f"H3>0.8 = {(base.PP_H3>0.8).sum()}", flush=True)

    print("\n=== p12 敏感性：各阈值下命中数 ===", flush=True)
    for p12 in P12_GRID:
        s = out[out.p12 == p12]
        print(f"  p12={p12:<8g} H4>0.8: {(s.PP_H4>0.8).sum():3d}   "
              f"H4>0.5: {(s.PP_H4>0.5).sum():3d}   "
              f"贝叶斯FDR@0.8 期望假发现: {s.loc[s.PP_H4>0.8,'PP_H4'].rsub(1).sum():.3f}", flush=True)

    hits = base[base.PP_H4 > 0.8]
    print("\n=== 基线三个命中在各 p12 下的 PP.H4 ===", flush=True)
    for _, h in hits.iterrows():
        s = out[(out.tissue == h.tissue) & (out.locus == h.locus) & (out.gene_id == h.gene_id)]
        vals = " ".join(f"{r.p12:g}:{r.PP_H4:.3f}" for r in s.itertuples())
        print(f"  {h.tissue:26s} {h.locus:16s} H4={h.PP_H4:.3f} "
              f"cond={h.cond_H4:.3f} 可信集={h.cred_set_size} "
              f"top={h.cred_set_top_pos}({h.cred_set_top_post})", flush=True)
        print(f"      p12 扫描 → {vals}", flush=True)
    print(f"\n总耗时 {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
