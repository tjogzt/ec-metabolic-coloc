"""
run_tissue.py — 单组织端到端：下载 eQTL Catalogue 全变异 → 建索引 → 标准 coloc

用法：python run_tissue.py QTD000121 adipose_visceral
（study 固定 QTS000015 = GTEx v8；eGene 列表从本地 GTEx tar 抽取，保证与组织一致）
"""
import argparse, pathlib, subprocess, sys, tarfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RES = ROOT / "pipeline" / "results"
TAR = DATA / "GTEx_Analysis_v8_eQTL.tar"
BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015"


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, check=True, **kw)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("qtd")
    ap.add_argument("tissue")
    a = ap.parse_args()

    gz = DATA / f"{a.qtd}.all.tsv.gz"
    # 用远端 Content-Length 判断完整性，避免把并发下载中的半成品当完整文件
    head = subprocess.run(f"curl -sI -m 60 {BASE}/{a.qtd}/{a.qtd}.all.tsv.gz | grep -i content-length",
                          shell=True, capture_output=True, text=True).stdout
    expect = int(head.split(":")[-1].strip()) if head.strip() else 0
    complete = gz.exists() and expect and gz.stat().st_size == expect
    if not complete:
        print(f"[1/4] 下载 {a.qtd} …", flush=True)
        sh(f'curl -sS -C - --retry 5 --retry-delay 5 -m 5400 -o "{gz}" {BASE}/{a.qtd}/{a.qtd}.all.tsv.gz')
        sh(f'curl -sS -o "{gz}.tbi" {BASE}/{a.qtd}/{a.qtd}.all.tsv.gz.tbi')
    print(f"[1/4] 文件 {gz.stat().st_size/1e9:.2f} GB", flush=True)

    # eGene 列表（与 GTEx 官方显著对一致）
    tissue_key = {"adipose_visceral": "Adipose_Visceral_Omentum",
                  "adipose_subcutaneous": "Adipose_Subcutaneous",
                  "liver": "Liver", "muscle": "Muscle_Skeletal",
                  "pancreas": "Pancreas", "blood": "Whole_Blood",
                  "ovary": "Ovary", "uterus": "Uterus"}[a.tissue]
    eg = RES / f"egenes_{a.tissue}.txt"
    print(f"[2/4] 抽取 {tissue_key} 的 eGene …", flush=True)
    genes = set()
    with tarfile.open(TAR, "r") as tf:
        for m in tf.getmembers():
            if "signif_variant_gene_pairs" in m.name and pathlib.Path(m.name).name.startswith(tissue_key + "."):
                f = tf.extractfile(m)
                next(f, None)
                for line in f:
                    g = line.decode().split("\t")[1].split(".")[0]
                    genes.add(g)
                break
    eg.write_text("\n".join(sorted(genes)) + "\n")
    print(f"      {len(genes):,} 个 eGene", flush=True)

    # 索引
    idx = RES / f"index_{a.tissue}.tsv"
    if not idx.exists() or idx.stat().st_size == 0:
        print("[3/4] 建基因窗口索引 …", flush=True)
        sh(f"""gunzip -c "{gz}" | awk -F'\\t' 'NR>1{{g=$1;c=$2;p=$3; if(!(g in lo)){{ch[g]=c;lo[g]=p;hi[g]=p}} else {{if(p<lo[g])lo[g]=p; if(p>hi[g])hi[g]=p}}}} END{{print "gene_id\\tchrom\\tlo\\thi"; for(g in lo) print g"\\t"ch[g]"\\t"lo[g]"\\t"hi[g]}}' > "{idx}" """)

    print("[4/4] coloc …", flush=True)
    sh(f'cd {ROOT}/pipeline && /opt/anaconda3/bin/python3 coloc_eqtlcat.py '
       f'--eqtl "{gz}" --index "{idx}" --genes "{eg}" --min_snps 10 '
       f'--out "results/coloc_{a.tissue}.tsv"')
    print("DONE", a.tissue, flush=True)


if __name__ == "__main__":
    main()
