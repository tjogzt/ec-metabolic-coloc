"""
w3_incremental.py — 轴 C：组织内 IR 评分对独立结局的增量预测

终点：GSE120490（145 例内膜样癌）的转移（matastasis: Yes/No），
      与评分来源（GSE17025 差异分析）完全独立。

做法：重复分层 5 折交叉验证的 AUC，比较
  基座 = 增殖分数 + 分级 + 肌层浸润≥50%
  基座 + IR 评分            ← 待检
  基座 + 疾病签名           ← 阳性对照（应显著增加）
  基座 + 随机模块（100 套）  ← 零分布，判定 ΔAUC 的经验 p
"""
from __future__ import annotations

import io
import gzip
import os
import pathlib

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

HERE = pathlib.Path(__file__).resolve().parent
W2 = HERE.parent / "w2"
DATA = W2 / "data"
RES = HERE / "results"
RES.mkdir(exist_ok=True)

CYCLING = ["MKI67", "PCNA", "CCNB1", "CCNB2", "CCNA2", "CDK1", "CDC20", "BUB1", "BUB1B",
           "TOP2A", "AURKA", "AURKB", "PLK1", "TPX2", "UBE2C", "RRM2", "TYMS", "MCM2",
           "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "CENPF", "KIF11", "KIF23", "ANLN",
           "PRC1", "ECT2", "NCAPG", "ASPM", "DLGAP5", "CENPE", "SMC4", "NUSAP1"]


def read_matrix(f: pathlib.Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    with gzip.open(f, "rt") as fh:
        lines = fh.readlines()
    i0 = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_begin"))
    i1 = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_end"))
    m = pd.read_csv(io.StringIO("".join(lines[i0 + 1:i1])), sep="\t", quotechar='"', index_col=0)
    # 样本注释
    meta = {}
    for l in lines[:i0]:
        if l.startswith("!Sample_"):
            k, _, v = l.partition("\t")
            meta.setdefault(k, []).append([x.strip('"') for x in v.rstrip("\n").split("\t")])
    return m, meta


def collapse(m: pd.DataFrame) -> pd.DataFrame:
    annot = pd.read_csv(DATA / "GPL570_annot.tsv", sep="\t")
    a = annot.set_index("ID")["Gene Symbol"]
    m = m.loc[m.index.intersection(a.index)]
    m = m.astype(float)
    mean = m.mean(axis=1)
    df = pd.DataFrame({"sym": a.reindex(m.index).values, "mean": mean.values})
    df["idx"] = np.arange(len(df))
    df = df[~df.sym.isin(["---", ""])].dropna(subset=["sym"])
    best = df.sort_values("mean", ascending=False).drop_duplicates("sym")
    out = m.iloc[best.idx.values]
    out.index = best.sym.values
    return out


def zscore(m: pd.DataFrame) -> pd.DataFrame:
    return m.sub(m.mean(axis=1), axis=0).div(m.std(axis=1).replace(0, np.nan), axis=0)


def score_set(z: pd.DataFrame, genes: list[str], up: list[str] | None = None) -> pd.Series | None:
    g = [x for x in genes if x in z.index]
    if len(g) < 5:
        return None
    if up is None:
        return z.loc[g].mean(axis=0)
    u = [x for x in up if x in z.index]
    d = [x for x in g if x not in u]
    if not u or not d:
        return None
    return z.loc[u].mean(axis=0) - z.loc[d].mean(axis=0)


def cv_auc(X: np.ndarray, y: np.ndarray, seed: int = 49) -> float:
    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=20, random_state=seed)
    aucs = []
    for tr, te in cv.split(X, y):
        m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1.0))
        m.fit(X[tr], y[tr])
        aucs.append(roc_auc_score(y[te], m.predict_proba(X[te])[:, 1]))
    return float(np.mean(aucs))


def main() -> None:
    raw, meta = read_matrix(DATA / "GSE120490_series_matrix.txt.gz")
    chars = meta["!Sample_characteristics_ch1"]
    def field(name: str) -> list[str]:
        for row in chars:
            if row[0].lower().startswith(name.lower()):
                return [x.split(": ", 1)[-1] for x in row]
        raise KeyError(name)

    met = field("matastasis")
    grade = field("grade")
    mi = field("mi50pluspercent")
    y = np.array([1 if v.strip().lower() == "yes" else 0 for v in met])
    print(f"样本 {len(y)}；转移 {y.sum()}（{100*y.mean():.1f}%）；"
          f"高级别 {sum(v=='Yes' for v in grade)}；肌层浸润≥50% {sum(v=='Yes' for v in mi)}")

    m = collapse(raw)
    print(f"基因级矩阵 {m.shape[0]:,} 基因 × {m.shape[1]} 样本")
    z = zscore(m)

    prog_path = os.environ.get("PROG_PATH", str(W2 / "results/IR_program_v1.tsv"))
    prog = pd.read_csv(prog_path, sep="\t", index_col=0)
    up_genes = list(prog.index[prog.direction == "UP"])
    ir = score_set(z, list(prog.index), up_genes)
    sig = pd.read_csv(W2 / "results/EC_signature_v2_prolifadj_compact.tsv", sep="\t")
    ecsig = score_set(z, list(sig.gene), list(sig.gene[sig.direction == "UP"]))
    prolif = z.loc[[g for g in CYCLING if g in z.index]].mean(axis=0)
    print(f"IR 评分可用基因 {sum(g in z.index for g in prog.index)}/{len(prog)}；"
          f"疾病签名 {sum(g in z.index for g in sig.gene)}/{len(sig)}")

    base = np.column_stack([
        prolif.values,
        np.array([1 if v == "Yes" else 0 for v in grade]),
        np.array([1 if v == "Yes" else 0 for v in mi]),
    ])
    auc_base = cv_auc(base, y)
    print(f"\n基座（增殖+分级+肌层浸润）CV AUC = {auc_base:.3f}")

    d_ir = cv_auc(np.column_stack([base, ir.values]), y) - auc_base
    d_sig = cv_auc(np.column_stack([base, ecsig.values]), y) - auc_base
    print(f"  + IR 评分     ΔAUC = {d_ir:+.4f}")
    print(f"  + 疾病签名    ΔAUC = {d_sig:+.4f}  ← 阳性对照")

    rm_path = os.environ.get("RM_PATH", str(W2 / "results/random_modules_v1.tsv"))
    rm = pd.read_csv(rm_path, sep="\t")
    # 零分布模块必须与当前程序等大，否则"随机模块"不是同一尺度的对照。
    # RM_PATH 缺省指向 v1 模块；若程序换成了 v2 而漏设 RM_PATH，这里会静默错配。
    rm_sizes = rm.groupby("module").size()
    prog_size = len(prog.index)
    if not (rm_sizes == prog_size).all():
        raise SystemExit(
            f"随机模块与程序不等大：模块 {sorted(rm_sizes.unique())} vs 程序 {prog_size}"
            f"（{rm_path}）。请确认 RM_PATH 与 PROG_PATH 配套。")
    null = []
    for name, g in rm.groupby("module"):
        s = score_set(z, list(g.gene))
        if s is None:
            continue
        null.append(cv_auc(np.column_stack([base, s.values]), y) - auc_base)
    null = np.array(null)
    p_ir = float((null >= d_ir).mean())
    p_sig = float((null >= d_sig).mean())
    print(f"\n随机模块 {len(null)} 套：ΔAUC 均值 {null.mean():+.4f}，"
          f"95% 分位 {np.percentile(null,95):+.4f}，最大 {null.max():+.4f}")
    print(f"经验 p（IR）  = {p_ir:.3f}")
    print(f"经验 p（签名）= {p_sig:.3f}")

    pd.DataFrame({"module": rm.module.unique()[:len(null)], "delta_auc": null}) \
        .to_csv(RES / f"incremental_null{os.environ.get('OUT_TAG','')}.tsv", sep="\t", index=False)
    pd.DataFrame([{
        "auc_base": auc_base, "delta_auc_IR": d_ir, "delta_auc_ECsig": d_sig,
        "null_mean": null.mean(), "null_p95": np.percentile(null, 95),
        "p_IR": p_ir, "p_ECsig": p_sig, "n": int(len(y)), "events": int(y.sum()),
    }]).to_csv(RES / f"incremental_result{os.environ.get('OUT_TAG','')}.tsv", sep="\t", index=False)
    print(f"\n已写 {RES}/incremental_result{os.environ.get('OUT_TAG','')}.tsv")


if __name__ == "__main__":
    main()
