"""
w2_annot.py — 候选 EC 队列的可用性面：supplementary 文件 + 样本注释字段
决定它们能否承担"独立复现"（需要：表达矩阵 + 组织学分型/分期/生存）
"""
import re
import subprocess
import time

CAND = ["GSE120490", "GSE17025", "GSE106191", "GSE196033", "GSE56026",
        "GSE216872", "GSE14860", "GSE178671", "GSE32507", "GSE21882"]


def curl(url: str, t: int = 120) -> str:
    return subprocess.run(["curl", "-sS", "-m", str(t), url],
                          capture_output=True, text=True).stdout


for acc in CAND:
    self_txt = curl(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
                    f"?acc={acc}&targ=self&form=text&view=brief")
    titles = [l.split(" = ", 1)[-1].strip() for l in self_txt.splitlines()
              if l.startswith("!Series_title")]
    types = [l.split(" = ", 1)[-1].strip() for l in self_txt.splitlines()
             if l.startswith("!Series_type")]
    n = sum(1 for l in self_txt.splitlines() if l.startswith("!Series_sample_id"))

    # supplementary（quick= 取文件列表）
    sup = curl(f"https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/")
    files = re.findall(r'href="([^"]+\.(?:gz|tar|txt|csv|xlsx|tsv))"[^<]*</a>\s*([\d\.]+[KMG]?)',
                       sup)

    # 注释字段
    gsm = curl(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
               f"?acc={acc}&targ=gsm&form=text&view=brief", 180)
    keys = {}
    for line in gsm.splitlines():
        if line.startswith("!Sample_characteristics_ch1"):
            v = line.split(" = ", 1)[-1]
            k = v.split(":")[0].strip() if ":" in v else v.strip()
            keys[k] = keys.get(k, 0) + 1
    nsample = len(re.findall(r"\^SAMPLE", gsm))

    print(f"\n{'='*70}\n{acc}  n={n}  (GSM {nsample})")
    print(f"  标题: {(titles or [''])[0][:95]}")
    print(f"  类型: {'; '.join(types)}")
    print(f"  补充文件: {[(f, s) for f, s in files][:6] if files else '（无或需展开）'}")
    top = sorted(keys.items(), key=lambda x: -x[1])[:9]
    print(f"  注释字段: {dict((k[:22], v) for k, v in top)}")
    time.sleep(0.5)
