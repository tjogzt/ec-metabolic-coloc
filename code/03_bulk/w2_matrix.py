"""
w2_matrix.py — 关键可用性检查：有没有 series matrix（已标准化的表达矩阵）
没有的话只剩 CEL/RAW.tar，而本机 R 不可用，处理成本会大幅上升。
"""
import re
import subprocess

CAND = ["GSE120490", "GSE17025", "GSE106191", "GSE196033", "GSE56026",
        "GSE216872", "GSE14860", "GSE178671", "GSE32507", "GSE21882"]

for acc in CAND:
    stub = acc[:-3] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{acc}/matrix/"
    out = subprocess.run(["curl", "-sS", "-m", "60", url],
                         capture_output=True, text=True).stdout
    files = re.findall(r'href="([^"]*series_matrix[^"]*)"[^<]*</a>\s*([\d\.]+[KMG]?)', out)
    print(f"{acc:<12} {'; '.join(f'{f} ({s})' for f, s in files) if files else '✗ 无 series matrix'}")
