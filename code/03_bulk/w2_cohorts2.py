"""
w2_cohorts2.py — 逐个拉 GEO 文本统计真实样本数（不用 esummary，它把嵌套 Sample 节点混进来）
"""
import re
import subprocess
import urllib.parse


def curl(url: str, t: int = 60) -> str:
    return subprocess.run(["curl", "-sS", "-m", str(t), url],
                          capture_output=True, text=True).stdout


def esearch(term: str, n: int):
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
           f"?db=gds&retmax={n}&term={urllib.parse.quote(term)}")
    return re.findall(r"<Id>(\d+)</Id>", curl(url))


def series_brief(acc: str) -> dict:
    url = (f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
           f"?acc={acc}&targ=self&form=text&view=brief")
    txt = curl(url)
    d = {"acc": acc, "n": 0, "type": [], "title": "", "pubmed": ""}
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
        elif k == "Series_pubmed_id":
            d["pubmed"] = v
    return d


def main():
    term = ('(endometrial cancer[Title] OR endometrial carcinoma[Title] OR endometrioid[Title] '
            'OR endometrial neoplasms[Title]) AND "Homo sapiens"[Organism] AND gse[Entry Type]')
    ids = esearch(term, 200)
    accs = [f"GSE{i[3:]}" for i in ids]     # 200012345 -> GSE12345
    print(f"检索命中 {len(accs)} 个 Series，逐个核实…\n")
    rows = []
    for i, acc in enumerate(accs, 1):
        d = series_brief(acc)
        if any("Expression profiling" in t for t in d["type"]) and d["n"] >= 40:
            rows.append(d)
        if i % 40 == 0:
            print(f"  … {i}/{len(accs)}，已入选 {len(rows)}")
    rows.sort(key=lambda x: -x["n"])
    print(f"\n表达谱且 n≥40 的队列 {len(rows)} 个：\n")
    print(f"{'GSE':<13}{'n':>5}  {'PMID':<10}标题")
    for d in rows:
        print(f"{d['acc']:<13}{d['n']:>5}  {d['pubmed']:<10}{d['title'][:80]}")
    with open("cohort_candidates.tsv", "w") as f:
        f.write("gse\tn\tpmid\ttitle\ttypes\n")
        for d in rows:
            f.write(f"{d['acc']}\t{d['n']}\t{d['pubmed']}\t{d['title']}\t{';'.join(d['type'])}\n")
    print("\n写入 cohort_candidates.tsv")


main()
