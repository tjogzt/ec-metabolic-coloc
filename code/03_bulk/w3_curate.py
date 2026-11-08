"""
w3_curate.py — 对候选 Series 逐条取摘要并分类（v2：修多行字段 + 分类优先级）

分类优先级（先特异后泛化）：
  EC_tumor → endometriosis → organoid → endometrium → decidua → off_topic
"""
from __future__ import annotations

import re
import time
import urllib.request

import pandas as pd

UA = {"User-Agent": "Mozilla/5.0 (research script)"}


def get(url: str, tries: int = 3) -> str:
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:  # noqa: BLE001
            time.sleep(2 + 2 * i)
    return ""


def series_fields(gse: str) -> dict:
    """把 acc.cgi 文本解析成字段→值，多行续行拼接。"""
    txt = get(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ=self&form=text&view=brief")
    out: dict[str, str] = {}
    cur = None
    for line in txt.splitlines():
        m = re.match(r"!(\w+) = (.*)", line)
        if m:
            cur = m.group(1)
            val = m.group(2).strip()
            out[cur] = (out[cur] + ";" + val) if cur in out else val
        elif line.startswith("!") or not line.strip():
            cur = None
        elif cur:
            out[cur] += " " + line.strip()
    return out


RULES = [
    ("EC_tumor", r"endometrial (cancer|carcinoma|tumou?r|neoplas)|endometrioid|uterine (cancer|carcinoma|serous)|carcinosarcoma"),
    ("endometriosis", r"endometrio(sis|tic)|menstrual (blood|effluent|fluid)|retrograde menstruation"),
    ("organoid", r"organoid|assembloid|blastoid|cell line|cultured |E6/E7|TERT-1"),
    ("endometrium", r"endometri(um|al)|uterus|uterine|Asherman|recurrent implantation failure|endometrial polyp|PCOS|polycystic"),
    ("decidua", r"decidual|trophoblast|placenta|implantation|blastocyst|embryo"),
]
OFF = r"microglia|APOE4|hippocamp|pan-cancer|natural killer cells|neutrophil|leukemia|triple-negative|colorectal|glioma"


def classify(title: str, summary: str, design: str) -> str:
    blob = f"{title} {summary} {design}".lower()
    for name, pat in RULES:
        if re.search(pat, blob):
            return name
    if re.search(OFF, blob):
        return "off_topic"
    return "off_topic"


cand = pd.read_csv("w3_candidates.tsv", sep="\t")
rows = []
for _, r in cand.iterrows():
    f = series_fields(r["gse"])
    nid = f.get("Series_sample_id", "")
    n = len([x for x in nid.split(";") if x.strip()])
    cls = classify(r["title"], f.get("Series_summary", ""), f.get("Series_overall_design", ""))
    rows.append({"gse": r["gse"], "n": n, "cls": cls,
                 "type": f.get("Series_type", ""), "title": r["title"][:95]})
    print(f"{r['gse']}  n={n:>5}  {cls:<14} {r['title'][:58]}")
    time.sleep(0.3)

d = pd.DataFrame(rows).sort_values(["cls", "n"], ascending=[True, False])
d.to_csv("w3_curated.tsv", sep="\t", index=False)
print("\n=== 分类汇总 ===")
print(d.cls.value_counts().to_string())
for c in ["EC_tumor", "endometrium", "endometriosis", "organoid"]:
    sub = d[d.cls == c]
    print(f"\n--- {c}（{len(sub)} 个，样本合计 {int(sub.n.sum())}）---")
    print(sub[["gse", "n", "title"]].to_string(index=False, max_colwidth=76))
