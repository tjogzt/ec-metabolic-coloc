#!/usr/bin/env python3
"""
05b_build_programme_v2.py
=========================

Purpose
-------
Rebuild the locked 862-gene obesity-associated programme and its 100 matched
random modules from the GSE25401 differential-expression table.

The programme is defined without any endometrial-cancer outcome label. Genes
are selected on FDR < 0.05 and |log2FC| >= 0.3 in the obese-versus-non-obese
contrast, then restricted to genes expressed above the lower quartile in
endometrial tissue (GSE17025). Random modules are drawn from the same expression
deciles and are exactly the programme size; a smaller module would not be a
valid size-matched null.

Inputs
------
<DATA>/DE_GSE25401_adipose_obese_state.tsv   GSE25401 obese vs non-obese DE
<DATA>/GSE17025_series_matrix.txt.gz         endometrial expression reference
<DATA>/GPL570_annot.tsv                      probe-to-symbol map
<TABLES>/DE_GSE18732_muscle_DM_state.tsv     optional direction check
<TABLES>/DE_GSE15773_adipose_IR_state.tsv    optional direction check

Outputs
-------
<TABLES>/metabolic_program_v2.tsv            862 genes, direction and logFC
<TABLES>/metabolic_program_v2.sha256         SHA-256 of the locked table
<TABLES>/random_modules_v2.tsv               100 modules of equal size

Notes
-----
This is the builder for the definition reported in the manuscript. The earlier
insulin-response definition is built by 05_build_programme.py and is retained
for provenance. The random-module seed is fixed in config.PROGRAMME_V2_SEED so
that the null distribution stays identical to the locked version.

The locked table is never overwritten silently. If the current DE table yields
a different gene set, the rebuild is written with a ``_rebuild`` suffix, the
added and removed genes are printed, and the process exits 2.
"""
from __future__ import annotations

import gzip
import hashlib
import io
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy.stats import binomtest

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from config.config import PATHS, PARAMS, PROGRAMME_V2_SEED  # noqa: E402

# Preregistered thresholds. They match 05_build_programme.INSULIN_LOGFC_CUT;
# they are analysis parameters, not result values.
FDR_CUT = PARAMS.program_de_padj
LOGFC_CUT = 0.3

# Thirty-five canonical cell-cycle genes. Reported for disclosure, not used to
# remove genes from the programme.
CYCLING = {
    "MKI67", "PCNA", "CCNB1", "CCNB2", "CCNA2", "CDK1", "CDC20", "BUB1",
    "BUB1B", "TOP2A", "AURKA", "AURKB", "PLK1", "TPX2", "UBE2C", "RRM2",
    "TYMS", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "CENPF",
    "KIF11", "KIF23", "ANLN", "PRC1", "ECT2", "NCAPG", "ASPM", "DLGAP5",
    "CENPE", "SMC4", "NUSAP1",
}


def ec_expression() -> pd.Series:
    """Mean endometrial expression per gene symbol, from GSE17025."""
    path = PATHS.data / "GSE17025_series_matrix.txt.gz"
    if not path.exists():
        raise SystemExit(f"missing endometrial reference: {path}")
    with gzip.open(path, "rt") as handle:
        lines = handle.readlines()
    start = next(i for i, line in enumerate(lines)
                 if line.startswith("!series_matrix_table_begin"))
    end = next(i for i, line in enumerate(lines)
               if line.startswith("!series_matrix_table_end"))
    matrix = pd.read_csv(io.StringIO("".join(lines[start + 1:end])),
                         sep="\t", quotechar='"', index_col=0)
    matrix = np.log2(matrix.apply(pd.to_numeric, errors="coerce").clip(lower=0) + 1)

    annot_path = PATHS.data / "GPL570_annot.tsv"
    if not annot_path.exists():
        raise SystemExit(f"missing platform annotation: {annot_path}")
    annot = pd.read_csv(annot_path, sep="\t").set_index("ID")["Gene Symbol"]
    matrix = matrix.loc[matrix.index.intersection(annot.index)]
    matrix["sym"] = annot.reindex(matrix.index).values
    matrix = matrix.dropna(subset=["sym"])
    matrix = matrix[~matrix.sym.isin(["---", ""])]
    return matrix.groupby("sym").mean(numeric_only=True).mean(axis=1)


def main() -> int:
    de_path = PATHS.tables / "DE_GSE25401_adipose_obese_state.tsv"
    if not de_path.exists():
        raise SystemExit(
            f"missing GSE25401 DE table: {de_path}\n"
            f"  run code/03_burden/02b_state_programme_de.R first")
    de = pd.read_csv(de_path, sep="\t", index_col=0)
    print(f"GSE25401 genes: {len(de):,}")

    significant = de[(de["adj.P.Val"] < FDR_CUT)
                     & (de["logFC"].abs() >= LOGFC_CUT)]
    print(f"significant (FDR<{FDR_CUT}, |logFC|>={LOGFC_CUT}): "
          f"{len(significant)} (UP {int((significant.logFC > 0).sum())} / "
          f"DOWN {int((significant.logFC < 0).sum())})")

    expression = ec_expression()
    threshold = expression.quantile(0.25)
    programme = pd.DataFrame(
        {"logFC_adipose_obese": significant.logFC,
         "FDR": significant["adj.P.Val"]},
        index=significant.index,
    )
    programme["expr_in_EC"] = programme.index.map(expression)
    before = len(programme)
    programme = programme[programme.expr_in_EC.fillna(-np.inf) > threshold]
    programme["direction"] = np.where(
        programme.logFC_adipose_obese > 0, "UP", "DOWN")
    print(f"endometrial expression filter: {before} -> {len(programme)} "
          f"(UP {int((programme.direction == 'UP').sum())} / "
          f"DOWN {int((programme.direction == 'DOWN').sum())})")

    cycle_hits = sorted(set(programme.index) & CYCLING)
    fraction = len(cycle_hits) / max(len(programme), 1)
    print(f"cell-cycle genes {len(cycle_hits)}/{len(programme)} = "
          f"{100 * fraction:.1f}%")
    if cycle_hits:
        print("  hits:", ", ".join(cycle_hits))

    for label, path in (
        ("GSE18732 muscle", PATHS.tables / "DE_GSE18732_muscle_DM_state.tsv"),
        ("GSE15773 adipose", PATHS.tables / "DE_GSE15773_adipose_IR_state.tsv"),
    ):
        if not path.exists():
            continue
        other = pd.read_csv(path, sep="\t", index_col=0)
        overlap = programme.index.intersection(other.index)
        same = (np.sign(programme.loc[overlap, "logFC_adipose_obese"])
                == np.sign(other.loc[overlap, "logFC"]))
        p = binomtest(int(same.sum()), len(overlap), 0.5).pvalue
        print(f"  vs {label}: overlap {len(overlap)}, same direction "
              f"{int(same.sum())} ({100 * same.mean():.1f}%), p={p:.2g}")

    locked_path = PATHS.tables / "metabolic_program_v2.tsv"
    same_as_locked = False
    if locked_path.exists():
        locked = pd.read_csv(locked_path, sep="\t", index_col=0)
        same_as_locked = (
            set(locked.index) == set(programme.index)
            and np.allclose(locked["logFC_adipose_obese"],
                            programme["logFC_adipose_obese"],
                            rtol=0, atol=1e-6, equal_nan=True)
            and np.allclose(locked["FDR"], programme["FDR"],
                            rtol=0, atol=1e-12, equal_nan=True)
            and (locked["direction"].reindex(programme.index)
                 == programme["direction"]).all()
        )

    if same_as_locked or not locked_path.exists():
        program_path = locked_path
        module_path = PATHS.tables / "random_modules_v2.tsv"
        status = 0
    else:
        program_path = PATHS.tables / "metabolic_program_v2_rebuild.tsv"
        module_path = PATHS.tables / "random_modules_v2_rebuild.tsv"
        status = 2
        added = sorted(set(programme.index) - set(locked.index))
        removed = sorted(set(locked.index) - set(programme.index))
        print(f"\nWARNING: the current DE table no longer reproduces the locked "
              f"programme ({len(programme)} rebuilt vs {len(locked)} locked; "
              f"{len(added)} added, {len(removed)} removed).")
        print("  Writing rebuild outputs with a _rebuild suffix; the locked "
              "table and its modules are left untouched.")
        print(f"  added: {', '.join(added[:12])}{' ...' if len(added) > 12 else ''}")
        print(f"  removed: {', '.join(removed[:12])}{' ...' if len(removed) > 12 else ''}")

    programme.to_csv(program_path, sep="\t", float_format="%.6g")
    digest = hashlib.sha256(program_path.read_bytes()).hexdigest()
    program_path.with_suffix(".sha256").write_text(
        f"{digest}  {program_path.name}\n")
    print(f"programme hash {digest} -> {program_path.name}")

    rng = np.random.default_rng(PROGRAMME_V2_SEED)
    expr = expression.dropna()
    decile = pd.qcut(expr.rank(method="first"), 10, labels=False)
    programme_expr = expr.reindex(programme.index).dropna()
    need = programme_expr.groupby(decile.reindex(programme_expr.index)).size().to_dict()
    programme_set = set(programme.index)
    pool: dict[int, list[str]] = {}
    for gene in expr.index:
        if gene in programme_set:
            continue
        pool.setdefault(int(decile.get(gene, -1)), []).append(gene)

    records = []
    for replicate in range(100):
        picks: list[str] = []
        for d, k in need.items():
            candidates = pool.get(int(d), [])
            if len(candidates) < k:
                continue
            picks += list(rng.choice(candidates, size=k, replace=False))
        for gene in picks:
            records.append({"module": f"R{replicate:03d}", "gene": gene})

    modules = pd.DataFrame(records).drop_duplicates()
    sizes = modules.groupby("module").size()
    bad = sizes[sizes != len(programme)]
    if len(bad):
        raise SystemExit(
            f"random modules of wrong size: {len(bad)} not {len(programme)} "
            f"genes (min {sizes.min()}, max {sizes.max()}); a decile pool is "
            f"too small")
    modules.to_csv(module_path, sep="\t", index=False)
    print(f"random modules: {modules.module.nunique()} sets, "
          f"{int(sizes.median())} genes each (programme {len(programme)}); "
          f"written {module_path.name}")
    return status


if __name__ == "__main__":
    sys.exit(main())
