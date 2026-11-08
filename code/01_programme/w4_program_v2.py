"""
w4_program_v2.py — 按预注册（主文件 + A1 + A2）构造 IR 状态程序 v2（慢性口径）

发现集（A2）：GSE25401 白色脂肪，肥胖(30) vs 非肥胖(26)
方向检验：GSE18732 肌肉 T2D vs NGT、GSE15773 脂肪 IR vs IS（均为功效不足的辅助集，只报一致率）

阈值全部来自预注册，不得在此调整：FDR<0.05 且 |log2FC|>=0.3。
同时按程序规模生成 100 套表达十分位匹配的随机模块（零分布须与程序等大，否则比较不公平）。
"""
from __future__ import annotations

import gzip
import hashlib
import io
import pathlib

import numpy as np
import pandas as pd
from scipy.stats import binomtest

HERE = pathlib.Path(__file__).resolve().parent
W2 = HERE.parent / "w2"
RES = HERE / "results"
RES.mkdir(exist_ok=True)
RNG = np.random.default_rng(20261103)

FDR, LFC = 0.05, 0.3          # 预注册阈值，写死
CYCLING = {"MKI67", "PCNA", "CCNB1", "CCNB2", "CCNA2", "CDK1", "CDC20", "BUB1", "BUB1B",
           "TOP2A", "AURKA", "AURKB", "PLK1", "TPX2", "UBE2C", "RRM2", "TYMS", "MCM2",
           "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "CENPF", "KIF11", "KIF23", "ANLN",
           "PRC1", "ECT2", "NCAPG", "ASPM", "DLGAP5", "CENPE", "SMC4", "NUSAP1"}


def ec_expression() -> pd.Series:
    f = W2 / "data/GSE17025_series_matrix.txt.gz"
    with gzip.open(f, "rt") as fh:
        lines = fh.readlines()
    i0 = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_begin"))
    i1 = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_end"))
    m = pd.read_csv(io.StringIO("".join(lines[i0 + 1:i1])), sep="\t", quotechar='"', index_col=0)
    m = np.log2(m.apply(pd.to_numeric, errors="coerce").clip(lower=0) + 1)
    annot = pd.read_csv(W2 / "data/GPL570_annot.tsv", sep="\t").set_index("ID")["Gene Symbol"]
    m = m.loc[m.index.intersection(annot.index)]
    m["sym"] = annot.reindex(m.index).values
    m = m.dropna(subset=["sym"])
    m = m[~m.sym.isin(["---", ""])]
    return m.groupby("sym").mean(numeric_only=True).mean(axis=1)


def main() -> None:
    disp = HERE / "DE_GSE25401_adipose_obese_state.tsv"
    if not disp.exists():
        raise SystemExit("缺少 GSE25401 DE 结果，先跑 w4_de.R")
    D = pd.read_csv(disp, sep="\t", index_col=0)
    print(f"发现集 GSE25401 基因 {len(D):,}")
    sig = D[(D["adj.P.Val"] < FDR) & (D.logFC.abs() >= LFC)]
    print(f"显著（FDR<{FDR}，|logFC|>={LFC}）：{len(sig)}"
          f"（UP {int((sig.logFC > 0).sum())} / DOWN {int((sig.logFC < 0).sum())}）")

    # ---- 内膜可测过滤（规则同 v1） ----
    ec = ec_expression()
    thr = ec.quantile(0.25)
    df = pd.DataFrame({"logFC_adipose_obese": sig.logFC, "FDR": sig["adj.P.Val"]}, index=sig.index)
    df["expr_in_EC"] = df.index.map(ec)
    n0 = len(df)
    df = df[df.expr_in_EC.fillna(-np.inf) > thr]
    df["direction"] = np.where(df.logFC_adipose_obese > 0, "UP", "DOWN")
    print(f"内膜表达过滤：{n0} → {len(df)}（UP {int((df.direction=='UP').sum())} / "
          f"DOWN {int((df.direction=='DOWN').sum())}）")

    # ---- 增殖纪律（预注册第 6 条） ----
    cyc = sorted(set(df.index) & CYCLING)
    frac = len(cyc) / max(len(df), 1)
    print(f"细胞周期基因 {len(cyc)}/{len(df)} = {100*frac:.1f}%"
          f"{'  ⚠ 超过 10%，正文须披露' if frac > 0.10 else '  ✓ 未超 10%'}")
    if cyc:
        print("  命中：", ", ".join(cyc))

    # ---- 方向一致性（辅助集功效不足，只报一致率，不作门控） ----
    for tag, f in [("肌肉 GSE18732", HERE / "DE_GSE18732_muscle_DM_state.tsv"),
                   ("脂肪 GSE15773", HERE / "DE_GSE15773_adipose_IR_state.tsv")]:
        if not f.exists():
            continue
        O = pd.read_csv(f, sep="\t", index_col=0)
        ov = df.index.intersection(O.index)
        same = np.sign(df.loc[ov, "logFC_adipose_obese"]) == np.sign(O.loc[ov, "logFC"])
        p = binomtest(int(same.sum()), len(ov), 0.5).pvalue
        print(f"  与 {tag}：重叠 {len(ov)}，同向 {int(same.sum())} "
              f"（{100*same.mean():.1f}%），p = {p:.2g}")

    df.to_csv(RES / "IR_state_program_v2.tsv", sep="\t", float_format="%.6g")
    h = hashlib.sha256((RES / "IR_state_program_v2.tsv").read_bytes()).hexdigest()
    (RES / "IR_state_program_v2.sha256").write_text(f"{h}  IR_state_program_v2.tsv\n")
    print(f"\n锁定哈希 {h}")

    # ---- 与程序等大的随机模块（零分布） ----
    expr = ec.dropna()
    dec = pd.qcut(expr.rank(method="first"), 10, labels=False)
    prog_expr = expr.reindex(df.index).dropna()
    need = prog_expr.groupby(dec.reindex(prog_expr.index)).size().to_dict()
    prog_set = set(df.index)
    pool: dict[int, list[str]] = {}
    for g in expr.index:
        if g in prog_set:
            continue
        pool.setdefault(int(dec.get(g, -1)), []).append(g)
    recs = []
    for rep in range(100):
        picks: list[str] = []
        for d, k in need.items():
            cand = pool.get(int(d), [])
            if len(cand) < k:
                continue
            picks += list(RNG.choice(cand, size=k, replace=False))
        for g in picks:
            recs.append({"module": f"R{rep:03d}", "gene": g})
    rm = pd.DataFrame(recs).drop_duplicates()
    # 零分布必须与程序等大：某个表达十分位候选池不足时上面会 continue，
    # 模块会静默变小，而下游"与程序等大"的比较前提随即失效。
    sizes = rm.groupby("module").size()
    bad = sizes[sizes != len(df)]
    if len(bad):
        raise SystemExit(
            f"随机模块尺寸不符：{len(bad)} 套不是 {len(df)} 基因"
            f"（最小 {sizes.min()}，最大 {sizes.max()}）；有十分位候选池不足")
    rm.to_csv(W2 / "results/random_modules_v2.tsv", sep="\t", index=False)
    print(f"随机模块 {rm.module.nunique()} 套，每套 "
          f"{int(rm.groupby('module').size().median())} 基因（与程序 {len(df)} 等大，表达十分位匹配）")
    print(f"→ {W2/'results/random_modules_v2.tsv'}")


if __name__ == "__main__":
    main()
