"""
w4_recon.py — 核实"慢性胰岛素抵抗状态"候选数据集
标准：① 人体在体组织（非细胞系）② 肥胖 vs 瘦 或 T2D vs 正常
      ③ 有 series matrix ④ 每侧 ≥15 例
"""
from __future__ import annotations

import re
import time
import urllib.request

UA = {"User-Agent": "Mozilla/5.0 (research script)"}
CAND = {
    "GSE13070": "骨骼肌 胰岛素敏感 vs 抵抗",
    "GSE25401": "骨骼肌 T2D vs 正常",
    "GSE18732": "骨骼肌 T2D vs 正常（含胰岛素输注）",
    "GSE16415": "脂肪组织 肥胖 vs 瘦",
    "GSE15773": "脂肪组织 肥胖 vs 瘦",
    "GSE20966": "胰岛 T2D vs 正常",
    "GSE15653": "肝 肥胖 vs 瘦",
    "GSE55650": "脂肪 胰岛素抵抗",
    "GSE38642": "脂肪/肌肉 胰岛素敏感",
    "GSE62164": "骨骼肌 胰岛素抵抗",
}


def get(url: str, tries: int = 2) -> str:
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=45) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:  # noqa: BLE001
            time.sleep(2)
    return ""


for gse, desc in CAND.items():
    txt = get(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ=self&form=text&view=brief")
    if not txt:
        print(f"{gse:10s} ✗ 取不到")
        continue
    title = re.search(r"!Series_title = (.*)", txt)
    summ = re.search(r"!Series_summary = (.*)", txt)
    n = len(re.findall(r"!Series_sample_id = ", txt))
    plats = sorted(set(re.findall(r"!Series_platform_id = (\S+)", txt)))
    mats = re.findall(r"!Series_supplementary_file = (\S+)", txt)
    gsm = get(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ=gsm&form=text&view=brief")
    titles = re.findall(r"!Sample_title = (.*)", gsm)[:6]
    print(f"\n{gse}  n={n:3d}  平台={','.join(plats)}  补充文件={len(mats)}")
    print(f"  用途: {desc}")
    print(f"  标题: {title.group(1)[:95] if title else '—'}")
    print(f"  摘要: {(summ.group(1)[:150] if summ else '—')}")
    print(f"  样本名样例: {' | '.join(t.strip()[:28] for t in titles)}")
    time.sleep(0.3)
