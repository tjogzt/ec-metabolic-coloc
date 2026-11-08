"""
w2_insulin.py — 找可定义"胰岛素响应程序"的人类扰动数据
只认：人类、有胰岛素或 IGF-1 处理组 + 对照、表达谱
"""
import re
import subprocess
import urllib.parse


def curl(url: str, t: int = 60) -> str:
    return subprocess.run(["curl", "-sS", "-m", str(t), url],
                          capture_output=True, text=True).stdout


def esearch(term: str, n: int = 60):
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
           f"?db=gds&retmax={n}&term={urllib.parse.quote(term)}")
    return re.findall(r"<Id>(\d+)</Id>", curl(url))


def brief(acc: str) -> dict:
    url = (f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
           f"?acc={acc}&targ=self&form=text&view=brief")
    txt = curl(url, 90)
    d = {"acc": acc, "n": 0, "type": [], "title": "", "summary": ""}
    for line in txt.splitlines():
        if not line.startswith("!"):
            continue
        k, _, v = line[1:].partition(" = ")
        k, v = k.strip(), v.strip()
        if k == "Series_sample_id":
            d["n"] += 1
        elif k == "Series_type":
            d["type"].append(v)
        elif k == "Series_title":
            d["title"] = v
        elif k == "Series_summary" and not d["summary"]:
            d["summary"] = v
    return d


QUERIES = [
    ('insulin treatment response[Title] AND "Homo sapiens"[Organism] AND gse[Entry Type]', 40),
    ('(insulin[Title] OR IGF-1[Title] OR IGF1[Title]) AND (adipocyte[Title] OR myotube[Title] '
     'OR skeletal muscle[Title]) AND "Homo sapiens"[Organism] AND gse[Entry Type]', 40),
    ('(insulin[Title] OR IGF1[Title] OR IGF-1[Title]) AND endometri*[Title] '
     'AND "Homo sapiens"[Organism] AND gse[Entry Type]', 40),
    ('hyperinsulinemia[Title] AND "Homo sapiens"[Organism] AND gse[Entry Type]', 20),
]

seen, out = set(), []
for term, n in QUERIES:
    ids = esearch(term, n)
    print(f"\n--- {term[:76]} → {len(ids)} 条 ---")
    for raw in ids:
        acc = f"GSE{raw[3:]}"
        if acc in seen:
            continue
        seen.add(acc)
        d = brief(acc)
        if d["n"] < 4 or not any("Expression profiling" in t for t in d["type"]):
            continue
        out.append(d)
        print(f"  {acc:<12} n={d['n']:<4} {d['title'][:88]}")

print(f"\n\n可用候选 {len(out)} 个，按设计相关性排序：")
kw = ("insulin", "igf", "homa", "resistan")
out.sort(key=lambda d: (-sum(k in (d["title"] + d["summary"]).lower() for k in kw), -d["n"]))
for d in out[:15]:
    print(f"\n{d['acc']} n={d['n']}  {d['title'][:100]}")
    print(f"   {d['summary'][:200]}")
