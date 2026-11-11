"""
08_rank_stability.py
====================

Purpose
-------
Measure how much of the programme score is explained by proliferation, for both
programme versions, by rank correlation between the two per-sample scores on the
endometrioid cohort.

This is the quantity that decides whether the incremental-value result can be
read at face value. A programme that tracks cell-cycle content will predict any
proliferation-driven endpoint, so its delta-AUC over a base model that already
contains a proliferation score is constrained before any model is fitted. If the
correlation is high, the adjusted analysis in stage 07 is the primary result; if
it is low, the two are complementary.

The functions, the collapse rule and the z-scoring are imported from stage 07
rather than re-implemented, so this check cannot drift onto a different scale
from the model it is meant to constrain.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
PATHS.data/GSE120490_series_matrix.txt.gz      cohort and expression matrix
PATHS.data/GPL570_annot.tsv                    probe-to-symbol map
PATHS.tables/metabolic_program_v2.tsv          programme v2 (shipped)
PATHS.tables/IR_program_v1.tsv                 programme v1 (stage 05)

Outputs
-------
None. Two correlations are printed; the values are quoted in the text, and
writing a table here would create a second source for a number that belongs to
stage 07's adjusted result.

Usage
-----
    python code/03_burden/08_rank_stability.py

Notes
-----
Spearman, not Pearson: both scores are means of z-scored genes, and their joint
distribution has heavy tails from samples with a few extreme genes. The rank
correlation is what the manuscript reports, so it is what is computed here.
"""
from __future__ import annotations

import pathlib
import sys
from importlib import import_module
from typing import List, Optional

import pandas as pd
from scipy.stats import spearmanr

# Repository root for `config`; this directory for the sibling stages.
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))

from config.config import PATHS

# Stage 07 owns the matrix decoder and the scoring rule; importing them keeps
# this check on exactly the same scale as the incremental-value model.
_incremental = import_module("07_incremental_value")


def _up_arm(prog: pd.DataFrame) -> Optional[List[str]]:
    """Return the up-regulated arm of a programme table, if it records direction.

    Parameters
    ----------
    prog : programme table indexed by gene symbol.

    Returns
    -------
    list of up genes, or None when the table carries no direction column, in
    which case the caller falls back to a plain mean score.
    """
    if "direction" not in prog.columns:
        return None
    return list(prog.index[prog.direction == "UP"])


def _correlate(z: pd.DataFrame, prog: pd.DataFrame, prolif: pd.Series,
               label: str) -> None:
    """Print the rank correlation between one programme score and proliferation.

    Parameters
    ----------
    z : gene x sample z-scored matrix.
    prog : programme table.
    prolif : per-sample proliferation score.
    label : name used in the printed line.

    Returns
    -------
    None.
    """
    score = _incremental.score_set(z, list(prog.index), _up_arm(prog))
    if score is None:
        print(f"{label}: not scorable (too few genes measurable in this cohort)")
        return
    common = score.notna() & prolif.notna()
    rho, p = spearmanr(score[common], prolif[common])
    print(f"{label}  Spearman rho = {rho:+.3f}  p = {p:.3e}  "
          f"n = {int(common.sum())}")


def main() -> None:
    """Score both programme versions on one cohort and report the correlations."""
    raw, _ = _incremental.read_matrix(PATHS.data / _incremental.ENDPOINT_SERIES)
    m = _incremental.collapse_to_gene(raw)
    z = _incremental.zscore(m)
    prolif = z.loc[[g for g in _incremental.CYCLING_GENES if g in z.index]].mean(axis=0)

    # Version 2: the 862-gene programme quoted in the manuscript.
    v2 = pd.read_csv(PATHS.tables / "metabolic_program_v2.tsv", sep="\t",
                     index_col=0)
    print(f"programme v2: {len(v2)} genes, columns {list(v2.columns)[:6]}")
    _correlate(z, v2, prolif, "programme v2 vs proliferation")

    # Version 1: the earlier, insulin-only programme, kept for comparison so
    # that a change in correlation can be attributed to the derivation change.
    v1 = pd.read_csv(PATHS.tables / "IR_program_v1.tsv", sep="\t", index_col=0)
    _correlate(z, v1, prolif, "programme v1 vs proliferation")


if __name__ == "__main__":
    main()
