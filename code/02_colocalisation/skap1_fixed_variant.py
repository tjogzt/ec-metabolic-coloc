#!/usr/bin/env python3
"""固定肝 lead eQTL 变异，逐组织取 SKAP1 的效应（同一变异，可跨组织比较）。"""
import subprocess, json, pathlib, csv, time

D = CFG.RESULTS
ROOT = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015"
GENE = "ENSG00000141293"
CHROM = "17"
LEAD = 48_127_316          # liver 窗口内 SKAP1 最强 eQTL 位置
I_GENE, I_POS, I_REF, I_ALT, I_P, I_B, I_SE, I_MAF = 0, 2, 3, 4, 8, 9, 10, 7

TISS = {"liver":"QTD000266","adipose":"QTD000116","adipose (visceral)":"QTD000121",
        "muscle":"QTD000281","pancreas":"QTD000296","uterus":"QTD000346",
        "vagina":"QTD000351","ovary":"QTD000291","blood":"QTD000356",
        "spleen":"QTD000326","LCL":"QTD000221","breast":"QTD000211",
        "thyroid":"QTD000341","kidney (cortex)":"QTD000261"}
tissues = TISS
rows = []

for tis, qtd in tissues.items():
    url = f"{ROOT}/{qtd}/{qtd}.all.tsv.gz"
    out = None
    for attempt in range(3):
        try:
            p = subprocess.run(["tabix", url, f"{CHROM}:{LEAD}-{LEAD}"],
                               capture_output=True, text=True, timeout=240)
            if p.returncode == 0:
                out = p.stdout
                break
            raise RuntimeError(p.stderr[:120])
        except Exception as e:
            if attempt == 2:
                print(f"{tis:26s} 失败", flush=True)
            else:
                time.sleep(12)
    if out is None:
        rows.append((tis, qtd, "FAIL", "", "", "", "", "", ""))
        continue

    hit = None
    for ln in out.splitlines():
        f = ln.split("\t")
        if len(f) < 11 or f[I_GENE] != GENE or f[I_POS] != str(LEAD):
            continue
        hit = f
        break

    if hit is None:
        rows.append((tis, qtd, "ABSENT", "", "", "", "", "", ""))
        print(f"{tis:26s} 该变异未检验 SKAP1", flush=True)
    else:
        pv, b, se = float(hit[I_P]), float(hit[I_B]), float(hit[I_SE])
        z = b / se if se > 0 else float("nan")
        rows.append((tis, qtd, "OK", f"{hit[I_REF]}/{hit[I_ALT]}", hit[I_MAF],
                     f"{pv:.3e}", f"{b:+.4f}", f"{z:+.2f}", hit[I_POS]))
        print(f"{tis:26s} {hit[I_REF]}/{hit[I_ALT]}  p={pv:.3e}  beta={b:+.4f}  z={z:+.2f}", flush=True)

outf = D / "skap1_fixed_variant.tsv"
with open(outf, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["tissue", "QTD", "status", "alleles", "maf", "p", "beta", "z", "pos"])
    w.writerows(rows)
print("\n写入", outf)
