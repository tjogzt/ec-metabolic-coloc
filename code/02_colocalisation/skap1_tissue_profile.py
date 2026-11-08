#!/usr/bin/env python3
"""SKAP1 (ENSG00000141293) 的 eQTL 跨组织分布。
统一窗口 chr17:48,000,000-48,300,000，远端 tabix。
列（eQTL Catalogue r7, 19 列）：0 molecular_trait_id(=gene) 2 position 7 maf 8 pvalue 9 beta 10 se
"""
import subprocess, json, pathlib, csv, time

D = CFG.RESULTS
ROOT = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015"
GENE = "ENSG00000141293"
CHROM, LO, HI = "17", 48_000_000, 48_300_000
I_GENE, I_POS, I_MAF, I_P, I_B, I_SE = 0, 2, 7, 8, 9, 10

tissues = json.load(open(D / "gtex_tissue.json"))
rows = []

for tis, qtd in tissues.items():
    url = f"{ROOT}/{qtd}/{qtd}.all.tsv.gz"
    out = None
    for attempt in range(3):
        try:
            p = subprocess.run(["tabix", url, f"{CHROM}:{LO}-{HI}"],
                               capture_output=True, text=True, timeout=900)
            if p.returncode == 0:
                out = p.stdout
                break
            raise RuntimeError(p.stderr[:120])
        except Exception as e:
            if attempt == 2:
                print(f"{tis:24s} 失败 {str(e)[:70]}", flush=True)
            else:
                time.sleep(12)
    if out is None:
        rows.append((tis, qtd, "FAIL", "", "", "", "", ""))
        continue

    n, best = 0, None
    for ln in out.splitlines():
        f = ln.split("\t")
        if len(f) < 11 or f[I_GENE] != GENE:
            continue
        try:
            pv, b, se = float(f[I_P]), float(f[I_B]), float(f[I_SE])
        except ValueError:
            continue
        n += 1
        if best is None or pv < best[0]:
            best = (pv, b, se, f[I_POS], f[I_MAF])

    if best is None:
        rows.append((tis, qtd, "NOT_TESTED", "0", "", "", "", ""))
        print(f"{tis:24s} SKAP1 不在检验集", flush=True)
    else:
        pv, b, se, pos, maf = best
        z = b / se if se > 0 else float("nan")
        rows.append((tis, qtd, "OK", str(n), pos, f"{pv:.3e}", f"{b:+.4f}", f"{z:+.2f}", maf))
        print(f"{tis:24s} n={n:5d}  min_p={pv:.3e}  beta={b:+.4f}  z={z:+.2f}", flush=True)

outf = D / "skap1_tissue_profile.tsv"
with open(outf, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["tissue", "QTD", "status", "n_variants", "lead_pos", "min_p", "beta", "z", "maf"])
    w.writerows(rows)
print("\n写入", outf, "| 行数", len(rows))
