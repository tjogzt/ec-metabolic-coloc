"""
w3_recon.py — 核实"40-60 个内膜单细胞数据集"这一说法

策略：GEO 多组检索式并集 → 逐条取标题/样本数/平台 →
过滤出真正的单细胞转录组（10x/ Smart-seq / Drop-seq），再人工判读。
"""
from __future__ import annotations

import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = {"User-Agent": "Mozilla/5.0 (research script)"}


def get(url: str, tries: int = 3) -> str:
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 + 2 * i)
    print(f"[warn] {last}", file=sys.stderr)
    return ""


QUERIES = [
    '(endometrial cancer[Title] OR endometrial carcinoma[Title] OR endometrium[Title] OR endometrial[Title]) AND (single cell[Title] OR single-cell[Title] OR scRNA[Title] OR single-cell RNA[Title] OR "single cell RNA"[Title])',
    'endometrial[All Fields] AND (10x[All Fields] OR "10X Genomics"[All Fields] OR "single-cell RNA sequencing"[All Fields]) AND gse[Entry Type]',
    '(endometriosis[Title] OR endometrium[Title]) AND (single-cell[Title] OR "single cell"[Title])',
]

ids: set[str] = set()
for q in QUERIES:
    url = (f"{EUTILS}/esearch.fcgi?db=gds&retmax=500&retmode=json&term="
           + urllib.parse.quote(q + ' AND "Homo sapiens"[Organism] AND gse[Entry Type]'))
    txt = get(url)
    if '"idlist"' in txt:
        import json
        got = json.loads(txt)["esearchresult"]["idlist"]
        print(f"检索式命中 {len(got)}：{q[:70]}…")
        ids |= set(got)
    time.sleep(0.4)

print(f"\n并集去重后 {len(ids)} 个 Series，开始取摘要…\n")
rows = []
ids_l = sorted(ids)
for i in range(0, len(ids_l), 100):
    chunk = ",".join(ids_l[i:i + 100])
    txt = get(f"{EUTILS}/esummary.fcgi?db=gds&retmode=json&id={chunk}")
    if not txt:
        continue
    import json
    try:
        d = json.loads(txt)["result"]
    except Exception:  # noqa: BLE001
        continue
    for k, v in d.items():
        if k == "uids":
            continue
        rows.append({
            "gse": v.get("accession", k),
            "n": v.get("n_samples", 0),
            "type": v.get("gdstype", ""),
            "title": v.get("title", ""),
            "pdat": v.get("pdat", ""),
            "taxon": v.get("taxon", ""),
        })
    time.sleep(0.4)

import pandas as pd
df = pd.DataFrame(rows)
df = df[df.taxon.str.contains("Homo sapiens", na=False)]
sc = df[df.type.str.contains("Expression profiling by high throughput sequencing", na=False)]
print(f"人类 Series {len(df)} 个；其中高通量测序表达谱 {len(sc)} 个\n")
sc = sc.sort_values("n", ascending=False)
pd.set_option("display.width", 200)
pd.set_option("display.max_colwidth", 85)
print(sc[["gse", "n", "pdat", "title"]].to_string(index=False))
sc.to_csv("w3_candidates.tsv", sep="\t", index=False)
print(f"\n已写 w3_candidates.tsv（{len(sc)} 行）")
