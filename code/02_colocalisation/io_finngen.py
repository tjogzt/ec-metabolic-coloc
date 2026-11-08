"""
io_finngen.py — FinnGen 汇总统计的读取与区域查询（基于 tabix）

FinnGen R12 格式：bgzip + tabix，列名见 HEADER。
坐标与 GTEx v8 一致，均为 GRCh38，无需 liftOver。
"""
from __future__ import annotations
import subprocess
import pathlib
from pathlib import Path

import pandas as pd

# FinnGen R12：列数因性状类型而异——
#   二元性状 13 列（多 af_alt_cases / af_alt_controls），数量性状 11 列（只有 af_alt）。
# 因此不能按位置硬编码，必须由表头决定。
COLUMNS = ["chrom", "pos", "ref", "alt", "rsids", "nearest_genes",
           "pval", "mlogp", "beta", "sebeta",
           "af_alt", "af_alt_cases", "af_alt_controls"]


class FinnGen:
    def __init__(self, path: str | Path):
        self.path = str(path)
        if not pathlib.Path(self.path).exists():
            raise FileNotFoundError(f"GWAS 文件不存在：{self.path}")
        # 先确认索引可用、染色体列表非空。被截断的 bgzip 会让 tabix 每次都静默返回 0 行，
        # 与"该区域确实没有关联"无法区分——本项目已因此得到过一次空结果。
        self.chroms = self._list_chroms()
        self.columns = self._read_header()

    def _list_chroms(self) -> list[str]:
        r = subprocess.run(["tabix", "-l", self.path],
                           capture_output=True, text=True, timeout=300)
        if r.returncode != 0:
            raise RuntimeError(f"tabix -l 失败（文件或索引损坏）：{self.path}\n{r.stderr.strip()}")
        chroms = [c.strip() for c in r.stdout.split() if c.strip()]
        if not chroms:
            raise RuntimeError(f"tabix -l 返回空染色体列表（文件可能被截断）：{self.path}")
        return chroms

    def _resolve_chrom(self, chrom: str) -> str:
        """兼容 chr1 / 1 两种命名；都不在文件里就报错，而不是静默返回空。"""
        if chrom in self.chroms:
            return chrom
        alt = f"chr{chrom}"
        if alt in self.chroms:
            return alt
        raise ValueError(
            f"染色体 {chrom} 不在 {pathlib.Path(self.path).name} 中"
            f"（实际前几个：{self.chroms[:5]}）")

    def _read_header(self) -> list[str]:
        """从 gz 表头取列名。用 head 截断，避免把整个文件解压进内存。"""
        try:
            out = subprocess.run(f'gunzip -c "{self.path}" | head -1',
                                 shell=True, capture_output=True, text=True,
                                 check=True, timeout=300).stdout.strip()
        except Exception:
            return COLUMNS[:11]
        names = [c.strip().lstrip("#") for c in out.split("\t")]
        return names if names and names[0] == "chrom" else COLUMNS[:11]

    def region(self, chrom: str, start: int, end: int) -> pd.DataFrame:
        chrom = self._resolve_chrom(str(chrom).replace("chr", ""))
        r = subprocess.run(["tabix", self.path, f"{chrom}:{start}-{end}"],
                           capture_output=True, text=True, timeout=300)
        if r.returncode != 0:
            raise RuntimeError(f"tabix 查询失败 {chrom}:{start}-{end}\n{r.stderr.strip()}")
        out = r.stdout
        if not out.strip():
            return pd.DataFrame(columns=self.columns)
        rows = [l.split("\t") for l in out.strip().split("\n")]
        df = pd.DataFrame(rows, columns=self.columns)
        for c in ("pos", "pval", "mlogp", "beta", "sebeta", "af_alt"):
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors="coerce")
        return df

    def lookup(self, chrom: str, pos: int, window: int = 0) -> pd.DataFrame:
        return self.region(chrom, max(1, pos - window), pos + window)


COMP = str.maketrans("ACGT", "TGCA")


def harmonise(df: pd.DataFrame, ref: str, alt: str) -> pd.DataFrame:
    """
    把 df 的效应统一到 ref/alt 参照序列上。
    返回带 beta_h / se_h 的副本；无法对齐的行被丢弃（记录在 .attrs）。
    """
    ref, alt = ref.upper(), alt.upper()
    keep, betas, ses = [], [], []
    for i, r in enumerate(df.itertuples(index=False)):
        r_ref, r_alt = str(r.ref).upper(), str(r.alt).upper()
        if (r_ref, r_alt) == (ref, alt):
            keep.append(i); betas.append(r.beta); ses.append(r.sebeta)
        elif (r_ref, r_alt) == (alt, ref):
            keep.append(i); betas.append(-r.beta); ses.append(r.sebeta)
        elif (r_ref, r_alt) == (ref.translate(COMP), alt.translate(COMP)):
            keep.append(i); betas.append(r.beta); ses.append(r.sebeta)
        elif (r_ref, r_alt) == (alt.translate(COMP), ref.translate(COMP)):
            keep.append(i); betas.append(-r.beta); ses.append(r.sebeta)
    out = df.iloc[keep].copy().reset_index(drop=True)
    out["beta_h"] = betas
    out["se_h"] = ses
    out.attrs["harmonised"] = len(keep)
    out.attrs["dropped"] = len(df) - len(keep)
    return out
