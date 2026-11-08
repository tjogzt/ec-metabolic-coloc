"""
w2_recon2.py — 用 GEO acc.cgi 文本接口核实四个队列的真实样本构成
（eutils 偶发 500，且 esummary 的 n_samples 会混入非表达样本）
"""
import re
import subprocess
import sys
import time

SERIES = ["GSE14860", "GSE293591", "GSE120490", "GSE17025"]


def fetch(acc: str, view: str = "brief") -> str:
    url = (f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
           f"?acc={acc}&targ=self&form=text&view={view}")
    r = subprocess.run(["curl", "-sS", "-m", "60", url], capture_output=True, text=True)
    return r.stdout


for acc in SERIES:
    txt = fetch(acc)
    if not txt.strip():
        print(f"{acc}: 取回为空")
        continue
    fields = {}
    for line in txt.splitlines():
        if line.startswith("!"):
            k, _, v = line[1:].partition(" = ")
            fields.setdefault(k.strip(), []).append(v.strip())
    samples = fields.get("Series_sample_id", [])
    types = fields.get("Series_type", [])
    title = (fields.get("Series_title") or [""])[0]
    plat = fields.get("Series_platform_id", [])
    print(f"\n{acc}  |  样本 {len(samples)}  |  平台 {','.join(plat)}")
    print(f"  标题: {title[:100]}")
    print(f"  类型: {'; '.join(types)}")
    summ = (fields.get("Series_summary") or [""])[0]
    print(f"  摘要: {summ[:150]}")
    time.sleep(0.5)

# GSE14860 单独看表达样本数：用 targ=gsm 全量抓，统计 title 里的特征
print("\n" + "=" * 60)
print("GSE14860 样本级核查（表达 vs SNP）")
url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14860&targ=gsm&form=text&view=brief"
txt = subprocess.run(["curl", "-sS", "-m", "120", url], capture_output=True, text=True).stdout
n_sup = len(re.findall(r"\^SAMPLE", txt))
titles = re.findall(r"!Sample_title = (.+)", txt)
lib = re.findall(r"!Sample_library_strategy = (.+)", txt)
print("样本数:", n_sup)
print("前 5 个标题:", [t.strip() for t in titles[:5]])
print("文库策略分布:", {x: lib.count(x) for x in set(lib)})
