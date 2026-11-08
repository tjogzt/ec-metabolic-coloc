"""
w2_recon.py — W2 前置核实
1) 四个候选 EC bulk 队列的实际元数据（平台、样本、supplementary 文件大小）
2) 搜索人类胰岛素/IGF-1 刺激转录组数据集（用于定义"胰岛素响应程序"）
3) 搜索胰岛、脂肪单细胞图谱（计划书 W2 明写要用）
"""
import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = {"User-Agent": "Mozilla/5.0 (research; contact tjogzt@gmail.com)"}


def get(url: str, tries: int = 3) -> str:
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            if i == tries - 1:
                raise
            time.sleep(2)
    return ""


def esearch(term: str, n: int = 40) -> list[str]:
    url = f"{EUTILS}/esearch.fcgi?db=gds&retmax={n}&term={urllib.parse.quote(term)}"
    return re.findall(r"<Id>(\d+)</Id>", get(url))


def esummary(ids: list[str]) -> list[dict]:
    if not ids:
        return []
    out = []
    for i in range(0, len(ids), 80):
        chunk = ids[i:i + 80]
        url = f"{EUTILS}/esummary.fcgi?db=gds&id={','.join(chunk)}"
        xml = get(url)
        root = ET.fromstring(xml)
        for doc in root.iter("DocSum"):
            d = {}
            for item in doc.iter("Item"):
                d[item.get("Name")] = (item.text or "").strip()
            out.append(d)
    return out


def show(rows, fields=("accession", "title", "n_samples", "gdstype", "taxon", "suppfile")):
    for r in rows:
        print({k: r.get(k, "") for k in fields if k in r})


print("=" * 70)
print("1) 四个候选 EC bulk 队列")
print("=" * 70)
for gse in ["GSE14860", "GSE293591", "GSE120490", "GSE17025"]:
    rows = esummary(esearch(f"{gse}[Accession]", 3))
    for r in rows:
        if r.get("accession") == gse:
            print(f"{gse}: n={r.get('n_samples')} | {r.get('gdstype')} | {r.get('taxon')}")
            print(f"   标题: {r.get('title')[:110]}")
            print(f"   平台: {r.get('gpl')} | supp: {r.get('suppfile','—')[:60]}")
    time.sleep(0.4)

print()
print("=" * 70)
print("2) 胰岛素／IGF-1 刺激转录组（人类）")
print("=" * 70)
terms = [
    'insulin[Title] AND (adipocyte[Title] OR adipose[Title]) AND "Homo sapiens"[Organism] AND gse[Entry Type]',
    'insulin[Title] AND (endometri*[Title]) AND "Homo sapiens"[Organism] AND gse[Entry Type]',
    'IGF-1[Title] AND (adipocyte[Title] OR endometri*[Title]) AND "Homo sapiens"[Organism] AND gse[Entry Type]',
    'insulin resistance[Title] AND (adipose[Title] OR endometrium[Title]) AND "Homo sapiens"[Organism] AND gse[Entry Type]',
]
seen = set()
for t in terms:
    rows = esummary(esearch(t, 25))
    print(f"\n--- {t[:70]} ---")
    for r in rows:
        acc = r.get("accession", "")
        if acc in seen:
            continue
        seen.add(acc)
        print(f"  {acc} n={r.get('n_samples')} | {r.get('title')[:95]}")
    time.sleep(0.4)

print()
print("=" * 70)
print("3) 胰岛／脂肪单细胞图谱")
print("=" * 70)
for t in ['(pancreatic islet[Title] OR islet[Title]) AND single cell[Title] AND "Homo sapiens"[Organism] AND gse[Entry Type]',
          '(adipose[Title]) AND single cell[Title] AND "Homo sapiens"[Organism] AND gse[Entry Type]']:
    rows = esummary(esearch(t, 20))
    print(f"\n--- {t[:70]} ---")
    for r in rows[:10]:
        print(f"  {r.get('accession')} n={r.get('n_samples')} | {r.get('title')[:95]}")
    time.sleep(0.4)
