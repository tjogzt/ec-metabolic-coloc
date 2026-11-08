"""
w2_lock.py — 锁定 W2 两个签名对象，并做必要的质量检查

产物（全部带哈希）：
  IR_program_v1.tsv            胰岛素响应程序（核心 78 + 扩展 68）
  EC_signature_v1_full.tsv     疾病签名（全量 3,588）
  EC_signature_v1_compact.tsv  疾病签名（紧凑版：上下调各 200，主分析用）
  W2_lock_manifest.txt         哈希清单
"""
from __future__ import annotations
import hashlib
import pathlib

import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

sig = pd.read_csv(HERE / "EC_signature_v1.tsv", sep="\t")
sig = sig.drop_duplicates("gene")
sig["abs"] = sig.logFC.abs()
full = sig.sort_values("abs", ascending=False)
full.to_csv(OUT / "EC_signature_v1_full.tsv", sep="\t", index=False, float_format="%.6g")

compact = pd.concat([
    full[full.direction == "UP"].nlargest(200, "abs"),
    full[full.direction == "DOWN"].nlargest(200, "abs"),
]).drop(columns=["abs"]).sort_values(["direction", "logFC"], ascending=[True, False])
compact.to_csv(OUT / "EC_signature_v1_compact.tsv", sep="\t", index=False, float_format="%.6g")

print(f"疾病签名：全量 {len(full)}（UP {int((full.direction=='UP').sum())} / "
      f"DOWN {int((full.direction=='DOWN').sum())}）")
print(f"          紧凑版 {len(compact)}（上下调各 200，主分析用）")

# ---- 增殖混杂检查：签名与经典细胞周期基因的重叠 ----
CYCLING = {"MKI67", "PCNA", "CCNB1", "CCNB2", "CCNA2", "CDK1", "CDC20", "BUB1", "BUB1B",
           "TOP2A", "AURKA", "AURKB", "PLK1", "TPX2", "UBE2C", "RRM2", "TYMS", "MCM2",
           "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "CENPF", "KIF11", "KIF23", "ANLN",
           "PRC1", "ECT2", "NCAPG", "ASPM", "DLGAP5", "CENPE", "SMC4", "NUSAP1"}
hit_c = sorted(CYCLING & set(compact.gene))
hit_f = sorted(CYCLING & set(full.gene))
print(f"\n增殖基因命中：紧凑版 {len(hit_c)}/{len(CYCLING)}，全量 {len(hit_f)}/{len(CYCLING)}")
print("  紧凑版命中：", ", ".join(hit_c) if hit_c else "无")
print("  全量命中  ：", ", ".join(hit_f) if hit_f else "无")
print("  判读：" + ("紧凑版未受增殖基因主导" if len(hit_c) <= 3 else
                  "⚠ 紧凑版仍含较多增殖基因，需在 ΔS 中做增殖校正"))

# ---- 哈希清单 ----
manifest = []
for f in ["IR_program_v1.tsv", "EC_signature_v1_full.tsv", "EC_signature_v1_compact.tsv"]:
    h = hashlib.sha256((OUT / f).read_bytes()).hexdigest()
    manifest.append(f"{h}  {f}")
(OUT / "W2_lock_manifest.txt").write_text("\n".join(manifest) + "\n")
print("\n哈希清单：")
print("\n".join(manifest))
