"""
w2_recon3.py — 逐样本核实：哪些 GSM 是"表达谱"，哪些是 SNP/其他
结论直接决定签名队列能不能写"三队列合并 ≥300 例"
"""
import collections
import re
import subprocess


def fetch_gsm(acc: str) -> str:
    url = (f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
           f"?acc={acc}&targ=gsm&form=text&view=brief")
    return subprocess.run(["curl", "-sS", "-m", "180", url],
                          capture_output=True, text=True).stdout


def parse(txt: str) -> list[dict]:
    blocks = txt.split("^SAMPLE = ")[1:]
    out = []
    for b in blocks:
        lines = b.splitlines()
        d = {"gsm": lines[0].strip()}
        for line in lines:
            if line.startswith("!Sample_title"):
                d["title"] = line.split(" = ", 1)[-1].strip()
            elif line.startswith("!Sample_platform_id"):
                d["gpl"] = line.split(" = ", 1)[-1].strip()
            elif line.startswith("!Sample_organism_ch1"):
                d["org"] = line.split(" = ", 1)[-1].strip()
            elif line.startswith("!Sample_source_name_ch1"):
                d["src"] = line.split(" = ", 1)[-1].strip()
        out.append(d)
    return out


for acc in ["GSE14860", "GSE293591", "GSE120490", "GSE17025"]:
    txt = fetch_gsm(acc)
    samples = parse(txt)
    print(f"\n{'='*66}\n{acc}: 共 {len(samples)} 个 GSM")
    by_plat = collections.Counter(s.get("gpl", "?") for s in samples)
    print("  平台分布:", dict(by_plat))
    for gpl, n in by_plat.most_common():
        ex = [s["title"] for s in samples if s.get("gpl") == gpl][:3]
        print(f"   {gpl} (n={n}) 例: {ex}")
