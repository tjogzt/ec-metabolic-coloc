"""
w2_cohorts.py — 重新清点可用的内膜癌 bulk 表达队列（逐个核实样本构成）
只认两件事：① gdstype 含 Expression profiling；② 样本标题真的是内膜癌。
"""
import re
import subprocess
import time

E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def curl(url: str, t: int = 90) -> str:
    return subprocess.run(["curl", "-sS", "-m", str(t), url],
                          capture_output=True, text=True).stdout


def esearch(term: str, n: int = 60) -> list[str]:
    import urllib.parse
    url = f"{E}/esearch.fcgi?db=gds&retmax={n}&term={urllib.parse.quote(term)}"
    return re.findall(r"<Id>(\d+)</Id>", curl(url))


def esummary(ids: list[str]) -> list[dict]:
    import xml.etree.ElementTree as ET
    out = []
    for i in range(0, len(ids), 100):
        url = f"{E}/esummary.fcgi?db=gds&id={','.join(ids[i:i+100])}"
        root = ET.fromstring(curl(url))
        for doc in root.iter("DocSum"):
            d = {}
            for item in doc.iter("Item"):
                d[item.get("Name")] = (item.text or "").strip()
            out.append(d)
    return out


term = ('(endometrial cancer[Title] OR endometrial carcinoma[Title] OR endometrioid[Title]) '
        'AND "Homo sapiens"[Organism] AND gse[Entry Type]')
ids = esearch(term, 120)
rows = esummary(ids)
time.sleep(0.5)

cand = []
for r in rows:
    acc = r.get("accession", "")
    n = r.get("n_samples", "0")
    typ = r.get("gdstype", "")
    try:
        n = int(n)
    except ValueError:
        n = 0
    if "Expression profiling" in typ and n >= 40:
        cand.append((acc, n, r.get("title", "")[:88], typ.split(";")[0]))

cand.sort(key=lambda x: -x[1])
print(f"命中 {len(rows)} 条，其中表达谱且 n≥40 的 {len(cand)} 条：\n")
print(f"{'GSE':<14}{'n':>5}  标题")
for acc, n, title, typ in cand:
    print(f"{acc:<14}{n:>5}  {title}")

# 对前 12 个做样本标题抽查，识别混杂队列
print("\n" + "=" * 74)
print("样本标题抽查（每队列取 4 条，判断是否真是内膜癌）")
print("=" * 74)
for acc, n, title, _ in cand[:12]:
    txt = curl(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}"
               f"&targ=gsm&form=text&view=brief", t=180)
    titles = re.findall(r"!Sample_title = (.+)", txt)
    org = set(re.findall(r"!Sample_organism_ch1 = (.+)", txt))
    print(f"\n{acc} (n={len(titles)}, 物种={','.join(sorted(org))[:24]})")
    print("   ", " | ".join(t.strip() for t in titles[:4]))
    time.sleep(0.6)
