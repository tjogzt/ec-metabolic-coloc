"""
w2_finalize.py — 锁定 W2 全部产物 + 独立队列复现检查
"""
from __future__ import annotations
import hashlib
import pathlib

import numpy as np
import pandas as pd
from scipy.stats import binomtest

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

prog = pd.read_csv(OUT / "IR_program_v1.tsv", sep="\t", index_col=0)
sig_full = pd.read_csv(HERE / "EC_signature_v2_prolifadj_full.tsv", sep="\t")
sig_comp = pd.read_csv(HERE / "EC_signature_v2_prolifadj_compact.tsv", sep="\t")

sig_full.to_csv(OUT / "EC_signature_v2_prolifadj_full.tsv", sep="\t", index=False)
sig_comp.to_csv(OUT / "EC_signature_v2_prolifadj_compact.tsv", sep="\t", index=False)

# ---- 独立队列复现（GSE196033，57 癌 vs 30 内膜） ----
rep = pd.read_csv(HERE / "DE_GSE196033_tumor_vs_normal.tsv", sep="\t", index_col=0)
for name, s in [("全量", sig_full), ("紧凑版", sig_comp)]:
    ov = s[s.gene.isin(rep.index)]
    a = np.sign(ov.logFC.values)
    b = np.sign(rep.loc[ov.gene, "logFC"].values)
    agree = int((a == b).sum())
    p = binomtest(agree, len(ov), 0.5).pvalue
    same_sig = int(((a == b) & (rep.loc[ov.gene, "adj.P.Val"].values < 0.05)).sum())
    print(f"{name}：重叠 {len(ov)}，方向一致 {agree}（{100*agree/len(ov):.1f}%，二项 p={p:.1e}），"
          f"同向且复现队列 FDR<0.05：{same_sig}")

# ---- 哈希 ----
manifest = []
for f in ["IR_program_v1.tsv", "EC_signature_v2_prolifadj_full.tsv",
          "EC_signature_v2_prolifadj_compact.tsv", "random_modules_v1.tsv"]:
    h = hashlib.sha256((OUT / f).read_bytes()).hexdigest()
    manifest.append(f"{h}  {f}")
(OUT / "W2_lock_manifest.txt").write_text("\n".join(manifest) + "\n")
print("\n" + "\n".join(manifest))

# ---- 程序与签名的重叠（用于后续 ΔS 的解释） ----
ir_genes = set(prog.index)
for name, s in [("全量", sig_full), ("紧凑版", sig_comp)]:
    ov = ir_genes & set(s.gene)
    print(f"\nIR 程序与疾病签名（{name}）重叠 {len(ov)} 个：" / None if False else
          f"\nIR 程序与疾病签名（{name}）重叠 {len(ov)} 个：{', '.join(sorted(ov)[:20])}")
