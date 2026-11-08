"""
liftover.py — 最小的 UCSC chain 转换实现（只做坐标映射，不处理 indel 边界）

GTEx v8 是 GRCh38，GCST006465 是 GRCh37，需要把 eQTL 变异坐标降到 hg19。
这里只实现正向块查找，足以应付 SNP（1bp 变异）。
"""
from __future__ import annotations
import gzip, pathlib
from collections import defaultdict


class Chain:
    def __init__(self, path: str | pathlib.Path):
        self.blocks = defaultdict(list)   # tName -> [(tStart, tEnd, qName, qStart, qStrand)]
        with gzip.open(path, "rt") as f:
            cur = None
            for line in f:
                if line.startswith("chain"):
                    p = line.split()
                    cur = {"tName": p[2].replace("chr", ""),
                           "tStart": int(p[5]), "tEnd": int(p[6]),
                           "qName": p[7].replace("chr", ""),
                           "qStart": int(p[10]), "qSize": int(p[8]),
                           "qStrand": p[9]}
                    t, q = cur["tStart"], cur["qStart"]
                elif line.startswith("#") or not line.strip():
                    continue
                elif cur is not None:
                    p = line.split()
                    size = int(p[0])
                    self.blocks[cur["tName"]].append(
                        (t, t + size, cur["qName"], q, cur["qStrand"],
                         cur["qSize"], size))
                    if len(p) >= 3:
                        t += size + int(p[1])
                        q += size + int(p[2])
                    else:
                        cur = None
        for k in self.blocks:
            self.blocks[k].sort()
        self._last = {}

    def lift(self, chrom: str, pos: int):
        """pos 为 1-based。返回 (qchrom, qpos) 或 None。"""
        chrom = str(chrom).replace("chr", "")
        pos0 = pos - 1
        bl = self.blocks.get(chrom)
        if not bl:
            return None
        # 二分：找到起始 <= pos0 的最后一个块
        lo, hi = 0, len(bl) - 1
        while lo <= hi:
            mid = (lo + hi) // 2
            if bl[mid][0] <= pos0:
                lo = mid + 1
            else:
                hi = mid - 1
        if hi < 0:
            return None
        t0, t1, qname, q0, qstrand, qsize, blk = bl[hi]
        if not (t0 <= pos0 < t1):
            return None
        off = pos0 - t0
        if qstrand == "+":
            return qname, q0 + off + 1
        # 反向块：chain 里的 q0 是负链坐标（0-based），
        # 负链坐标 r 对应正链 0-based 坐标 qSize-1-r，1-based 即 qSize-r
        return qname, qsize - q0 - off
