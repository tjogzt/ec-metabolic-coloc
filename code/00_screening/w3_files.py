"""
w3_files.py — 核对核心候选的补充文件格式（能否直接取到表达矩阵）
"""
from __future__ import annotations

import time
import urllib.request
import re

UA = {"User-Agent": "Mozilla/5.0 (research script)"}
CORE = {
    "GSE278879": "内膜样癌 snRNA-seq（最对靶）",
    "GSE268903": "MMRd 内膜癌 新辅助免疫治疗 I 期",
    "GSE251923": "内膜癌 抗 PD-1 TME + 空间",
    "GSE225691": "Netrin-1 内膜癌",
    "GSE212217": "MMR 缺陷（含内膜癌）",
    "GSE299623": "内膜癌肉瘤",
    "GSE111976": "正常内膜图谱（最大）",
    "GSE250130": "黄体期内膜时序",
    "GSE247962": "黄体期内膜",
    "GSE205481": "TRIM28 内膜",
}


def get(url: str, tries: int = 3) -> str:
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:  # noqa: BLE001
            time.sleep(2 + 2 * i)
    return ""


for gse, desc in CORE.items():
    txt = get(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ=self&form=text&view=brief")
    files = re.findall(r"!Series_supplementary_file = (\S+)", txt)
    sizes = re.findall(r"!Series_supplementary_file = (\S+)\s*\n!Series_supplementary_file_size = (\S+)", txt)
    gpl = re.findall(r"!Series_platform_id = (\S+)", txt)
    print(f"\n=== {gse}  {desc}")
    print(f"    平台 {','.join(sorted(set(gpl)))}  补充文件 {len(files)} 个")
    for f in files[:6]:
        name = f.rsplit("/", 1)[-1]
        kind = ("10x-mtx" if re.search(r"(matrix|barcodes|features)", name, re.I)
                else "h5/h5ad" if re.search(r"\.(h5|h5ad)$", name, re.I)
                else "tar/gz" if re.search(r"\.(tar|tgz)(\.gz)?$", name, re.I)
                else "other")
        print(f"      - {name}  [{kind}]")
    if len(files) > 6:
        print(f"      … 另 {len(files) - 6} 个")
    time.sleep(0.35)
