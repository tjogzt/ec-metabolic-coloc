"""
01_chain_liftover.py
====================

Purpose
-------
Minimal UCSC chain-file reader that maps single-base coordinates between
genome builds. Used because the eQTL Catalogue release is on GRCh38 while the
secondary risk dataset is on GRCh37, and the two must be joined on coordinates
rather than on rsIDs.

Only whole-block lookups are implemented, which is sufficient for SNVs. Indels
that straddle a block boundary are not handled; callers must have filtered to
single-base substitutions before lifting.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
Chain gzip file (for example hg19ToHg38.over.chain.gz). Paths are resolved by
the caller from config.PATHS.data; nothing is embedded here.

Outputs
-------
Chain.lift(chrom, pos) returns (target_chrom, target_pos) on the other build,
or None when the position falls outside every block of the chain.

Notes
-----
Chromosome names have the "chr" prefix stripped so that the reader matches both
UCSC-style ("chr17") and Ensembl-style ("17") identifiers.
"""
from __future__ import annotations

import gzip
import pathlib
from collections import defaultdict
from typing import Dict, List, Optional, Tuple


class Chain:
    """Forward-block lookup over one UCSC chain file.

    Attributes
    ----------
    blocks : dict
        source chromosome -> sorted list of block tuples
        (t_start, t_end, target_chrom, q_start, q_strand, q_size, block_size).
        All coordinates are 0-based half-open.
    """

    def __init__(self, path: "str | pathlib.Path"):
        # Each block is stored with the metadata needed to invert the strand.
        self.blocks: Dict[str, List[Tuple[int, int, str, int, str, int, int]]] = defaultdict(list)
        self._parse(str(path))
        for chrom in self.blocks:
            self.blocks[chrom].sort()

    def _parse(self, path: str) -> None:
        """Read the chain file and fill ``self.blocks``.

        A chain header line starts with "chain" and describes the alignment.
        Each following non-comment line is one block: size, gap on the target,
        gap on the query. The running target and query offsets advance by
        block size plus the corresponding gap.
        """
        with gzip.open(path, "rt") as fh:
            current = None
            t_pos = q_pos = 0          # running offsets; defined before the header is parsed
            for line in fh:
                if line.startswith("chain"):
                    p = line.split()
                    current = {
                        "t_name": p[2].replace("chr", ""),      # target chromosome
                        "t_start": int(p[5]),
                        "t_end": int(p[6]),
                        "q_name": p[7].replace("chr", ""),      # query chromosome
                        "q_start": int(p[10]),
                        "q_size": int(p[8]),
                        "q_strand": p[9],
                    }
                    t_pos, q_pos = current["t_start"], current["q_start"]
                elif line.startswith("#") or not line.strip():
                    continue
                elif current is not None:
                    p = line.split()
                    size = int(p[0])
                    self.blocks[current["t_name"]].append(
                        (t_pos, t_pos + size, current["q_name"], q_pos,
                         current["q_strand"], current["q_size"], size)
                    )
                    if len(p) >= 3:
                        # More blocks follow: advance past this block and its gaps.
                        t_pos += size + int(p[1])
                        q_pos += size + int(p[2])
                    else:
                        current = None                          # last block of the chain

    def lift(self, chrom: str, pos: int) -> Optional[Tuple[str, int]]:
        """Map a 1-based position on the source build to the target build.

        Parameters
        ----------
        chrom : str
            Source chromosome, with or without the "chr" prefix.
        pos : int
            1-based source position.

        Returns
        -------
        (target_chrom, target_pos) with target_pos 1-based, or None if the
        position is not covered by any block of the chain.

        Notes
        -----
        The block is found by binary search over block starts, then the offset
        within the block is applied. For reverse-strand blocks the chain stores
        a 0-based coordinate counted from the end of the query chromosome, so
        the 1-based forward coordinate is q_size - q_start - offset.
        """
        chrom = str(chrom).replace("chr", "")
        pos0 = pos - 1
        blocks = self.blocks.get(chrom)
        if not blocks:
            return None

        # Last block whose start is <= pos0.
        lo, hi = 0, len(blocks) - 1
        while lo <= hi:
            mid = (lo + hi) // 2
            if blocks[mid][0] <= pos0:
                lo = mid + 1
            else:
                hi = mid - 1
        if hi < 0:
            return None

        t0, t1, q_name, q0, q_strand, q_size, _ = blocks[hi]
        if not (t0 <= pos0 < t1):
            return None
        offset = pos0 - t0
        if q_strand == "+":
            return q_name, q0 + offset + 1
        return q_name, q_size - q0 - offset
