"""
05_build_programme.py
=====================

Purpose
-------
Build, audit and lock the two signature objects the rest of the package scores:
the insulin-response programme and the endometrial-cancer disease signature.

Three steps, in the order the objects are frozen:

    build    combine the two insulin contrasts into a tiered programme, filter
             to genes that are measurable in endometrial tissue, and draw
             expression-matched random modules of the same size.
    lock     turn the disease signature produced by stage 02 into a full and a
             compact form, audit how much of the compact form is cell cycle
             driven, and record SHA-256 hashes.
    finalise rewrite the proliferation-adjusted signatures into the tables
             directory, check them against the independent cohort, and write the
             manifest that the manuscript quotes.

The programme is defined without any outcome label: the two insulin datasets
carry no tumour information, and the endometrial filter uses only the tissue
expression distribution. That is what allows the later incremental-value test to
be called independent of the endpoint.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
PATHS.tables/DE_GSE24215_insulin_vs_basal.tsv      stage 02, muscle
PATHS.tables/DE_GSE24422_insulin_adipocyte.tsv     stage 02, adipocytes
PATHS.tables/EC_signature_v1.tsv                   stage 02, discovery signature
PATHS.tables/EC_signature_v2_prolifadj_full.tsv    stage 02
PATHS.tables/EC_signature_v2_prolifadj_compact.tsv stage 02
PATHS.tables/DE_GSE196033_tumor_vs_normal.tsv      stage 02, reproduction cohort
PATHS.data/GSE17025_series_matrix.txt.gz           tissue-expression filter
PATHS.data/GPL570_annot.tsv                        probe-to-symbol map

config.SEED                          the only random seed in the package
config.PARAMS.program_de_padj        DE significance for both insulin arms
config.PARAMS.n_random_modules       number of expression-matched null modules

Outputs
-------
PATHS.tables/IR_program_v1.tsv           direction, tier, four effect columns,
                                         two QC columns; one row per gene
PATHS.tables/IR_program_v1.sha256        hash sidecar in `sha256sum` format
PATHS.tables/random_modules_v1.tsv       module, gene
PATHS.tables/EC_signature_v1_full.tsv    gene, direction, logFC, FDR
PATHS.tables/EC_signature_v1_compact.tsv same, top 200 per direction
PATHS.tables/EC_signature_v2_prolifadj_*.tsv    normalised in place
PATHS.tables/W2_lock_manifest.txt        hash manifest of the locked objects

Usage
-----
    python code/03_burden/05_build_programme.py

Notes
-----
``CYCLING_GENES`` is defined here because this is the first module that needs
the cell-cycle panel (for the confound audit); stages 06-08 import it from this
file so that a single definition decides what "proliferation" means.

The random modules are drawn from the same expression decile as each programme
gene, so a null module has the same size and a similar expression profile as the
real programme. A module that is shorter than the programme cannot serve as a
size-matched control, and the caller (07) refuses to score with one; the subset
is not padded to hide this.

Seed note: the original scripts carried a literal seed (20261102) at this exact
site. Every draw now comes from config.SEED, so the realised modules differ from
the ones shipped in the first pass while the procedure is unchanged.
"""
from __future__ import annotations

import gzip
import hashlib
import io
import pathlib
import sys
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

# Repository root, so that `config` resolves when the file is run as a script.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

from config.config import PARAMS, PATHS, SEED

# --------------------------------------------------------------------------
# Module constants
# --------------------------------------------------------------------------
# Effect-size cut on the insulin contrast. Smaller than the tumour cut because a
# paired muscle biopsy and a cultured adipocyte experiment resolve smaller
# shifts than a 91-vs-11 tumour comparison.
INSULIN_LOGFC_CUT = 0.3

# Expression filter: a gene is "measurable in endometrium" if its mean level in
# GSE17025 exceeds the 25th percentile of all measured genes. The quartile is
# deliberately permissive -- this removes dark features, not lowly expressed but
# real transcripts.
EXPRESSION_QUANTILE = 0.25

# Deciles used to match random modules to the programme's expression profile.
N_EXPRESSION_DECILES = 10

# Tier labels, in the order they are stored in the programme table.
TIERS = ("core", "extended")

# Genes per direction in the compact signature.
COMPACT_PER_DIRECTION = 200

# Audit threshold: a compact signature built from a tumour contrast will always
# contain some cell-cycle genes; more than three means the contrast is being
# driven by proliferation and the adjusted analysis is mandatory.
MAX_PROLIF_IN_COMPACT = 3

# Thirty-five canonical cell-cycle genes. Used to score proliferation and to
# audit the signature; defined once, imported by stages 06-08.
CYCLING_GENES: Tuple[str, ...] = (
    "MKI67", "PCNA", "CCNB1", "CCNB2", "CCNA2", "CDK1", "CDC20", "BUB1", "BUB1B",
    "TOP2A", "AURKA", "AURKB", "PLK1", "TPX2", "UBE2C", "RRM2", "TYMS", "MCM2",
    "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "CENPF", "KIF11", "KIF23", "ANLN",
    "PRC1", "ECT2", "NCAPG", "ASPM", "DLGAP5", "CENPE", "SMC4", "NUSAP1",
)

# Column order of the programme table. Frozen: downstream stages index by name,
# but the shipped table is read by eye and its order is part of the record.
PROGRAMME_COLUMNS = ("direction", "tier", "logFC_muscle", "FDR_muscle",
                     "logFC_adip", "FDR_adip", "same_direction", "expr_in_EC",
                     "expressed_in_EC")


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def sha256_of(path: pathlib.Path) -> str:
    """Return the hex SHA-256 of a file.

    Parameters
    ----------
    path : file to hash.

    Returns
    -------
    Lowercase hex digest.
    """
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_tsv(name: str, index_col: int | None = 0) -> pd.DataFrame:
    """Read a table from the tables directory, failing with a hint if absent.

    Parameters
    ----------
    name : file name inside ``PATHS.tables``.
    index_col : column to use as the row index, or None to keep a RangeIndex.
    """
    path = PATHS.tables / name
    if not path.exists():
        raise SystemExit(f"missing input table: {path}\n"
                         f"  run code/03_burden/02_differential_expression.R first")
    return pd.read_csv(path, sep="\t", index_col=index_col)


def ec_expression() -> pd.Series:
    """Mean expression of every gene in the endometrial tumour cohort.

    Returns
    -------
    Series indexed by gene symbol, one value per gene (mean over samples of the
    gene-level log2 expression).

    Notes
    -----
    GSE17025 is used as the expression reference because it is the largest
    endometrioid cohort in the package and its platform (GPL570) annotates
    almost every transcript. The matrix is parsed by hand rather than through
    GEOquery so that this filter can run in the Python half of the pipeline
    without an R call.

    Values are clipped at zero before the log because the series matrix carries
    a handful of slightly negative intensities after normalisation, and log2 of
    a negative number would silently become NaN and drop the gene.
    """
    path = PATHS.data / "GSE17025_series_matrix.txt.gz"
    if not path.exists():
        raise SystemExit(f"missing expression matrix: {path}")
    with gzip.open(path, "rt") as fh:
        lines = fh.readlines()
    i0 = next(i for i, l in enumerate(lines)
              if l.startswith("!series_matrix_table_begin"))
    i1 = next(i for i, l in enumerate(lines)
              if l.startswith("!series_matrix_table_end"))
    m = pd.read_csv(io.StringIO("".join(lines[i0 + 1:i1])),
                    sep="\t", quotechar='"', index_col=0)
    m = m.apply(pd.to_numeric, errors="coerce")
    m = np.log2(m.clip(lower=0) + 1)

    # Probe -> symbol, then mean over probes per symbol.
    annot = pd.read_csv(PATHS.data / "GPL570_annot.tsv", sep="\t")
    symbols = annot.set_index("ID")["Gene Symbol"]
    m = m.loc[m.index.intersection(symbols.index)]
    m["sym"] = symbols.reindex(m.index).values
    m = m.dropna(subset=["sym"])
    m = m[~m.sym.isin(["---", ""])]
    return m.groupby("sym").mean(numeric_only=True).mean(axis=1)


# --------------------------------------------------------------------------
# Step 1: build the programme
# --------------------------------------------------------------------------
def build_programme() -> pd.DataFrame:
    """Derive the tiered insulin-response programme and draw null modules.

    Returns
    -------
    The programme table as written, indexed by gene symbol.

    Notes
    -----
    Two tiers. ``core`` requires significance and the same direction in both
    insulin datasets; ``extended`` requires significance and direction agreement
    in muscle only, provided the gene is measurable in endometrium. Both tiers
    are reported, and the compact disease signature is what the incremental
    model scores, so the tier boundary changes the gene count but not the
    endpoint.
    """
    PATHS.ensure()
    A = read_tsv("DE_GSE24215_insulin_vs_basal.tsv")
    B = read_tsv("DE_GSE24422_insulin_adipocyte.tsv")
    print(f"GSE24215 (skeletal muscle) genes {len(A):,}; "
          f"GSE24422 (adipocytes) genes {len(B):,}")

    common = A.index.intersection(B.index)
    df = pd.DataFrame({
        "logFC_muscle": A.loc[common, "logFC"],
        "FDR_muscle": A.loc[common, "adj.P.Val"],
        "logFC_adip": B.loc[common, "logFC"],
        "FDR_adip": B.loc[common, "adj.P.Val"],
    }, index=common)

    # ---- tissue filter: is the gene measurable in endometrium? -----------
    ec_expr = ec_expression()
    threshold = ec_expr.quantile(EXPRESSION_QUANTILE)
    print(f"genes measurable in endometrial tissue {len(ec_expr):,}; "
          f"q{int(EXPRESSION_QUANTILE * 100)} threshold {threshold:.2f}")

    df["expr_in_EC"] = df.index.map(ec_expr)
    # A gene absent from the reference matrix cannot be shown to be expressed,
    # so it is mapped to -inf rather than to the mean.
    df["expressed_in_EC"] = df["expr_in_EC"].fillna(-np.inf) > threshold

    # ---- tier assignment ------------------------------------------------
    sig_muscle = (df.FDR_muscle < PARAMS.program_de_padj) & \
                 (df.logFC_muscle.abs() >= INSULIN_LOGFC_CUT)
    same_dir = np.sign(df.logFC_muscle) == np.sign(df.logFC_adip)
    df["same_direction"] = same_dir

    core = (sig_muscle & same_dir
            & (df.FDR_adip < PARAMS.program_de_padj) & df.expressed_in_EC)
    ext = sig_muscle & same_dir & df.expressed_in_EC & ~core

    df["tier"] = np.where(core, "core", np.where(ext, "extended", "excluded"))
    df["direction"] = np.where(df.logFC_muscle > 0, "UP", "DOWN")

    prog = df[df.tier.isin(TIERS)].copy()
    prog = prog.sort_values(["tier", "FDR_muscle"])
    print(f"\ncore (significant and concordant in both, measurable in EC): "
          f"{int(core.sum())} genes")
    print(f"extended (muscle only, measurable in EC)                  : "
          f"{int(ext.sum())} genes")
    print(f"  UP {int((prog.direction == 'UP').sum())} / "
          f"DOWN {int((prog.direction == 'DOWN').sum())}")

    dest = PATHS.tables / "IR_program_v1.tsv"
    prog[list(PROGRAMME_COLUMNS)].to_csv(dest, sep="\t", float_format="%.6g")
    digest = sha256_of(dest)
    (PATHS.tables / "IR_program_v1.sha256").write_text(f"{digest}  IR_program_v1.tsv\n")
    print(f"\nsha256 {digest}")

    _draw_random_modules(prog, ec_expr)
    return prog


def _draw_random_modules(prog: pd.DataFrame, ec_expr: pd.Series) -> pd.DataFrame:
    """Draw expression-matched random gene modules of the programme's size.

    Parameters
    ----------
    prog : the programme table (genes to exclude from the candidate pool).
    ec_expr : per-gene mean expression used to define the deciles.

    Returns
    -------
    The module table written to disk: columns ``module`` and ``gene``.

    Notes
    -----
    For each expression decile, as many random genes are drawn as the programme
    contributes genes to that decile. Matching on expression is what makes the
    null modules comparable: an unmatched random set would be systematically
    less expressed, and score differences would then reflect expression level
    rather than programme content.

    A decile with too few candidates is skipped, which yields a module shorter
    than the programme; the caller detects and rejects that rather than padding
    it with genes from another decile.
    """
    expr_all = ec_expr.dropna()
    deciles = pd.qcut(expr_all.rank(method="first"), N_EXPRESSION_DECILES,
                      labels=False)

    # Candidate pool: every measured gene outside the programme, by decile.
    programme_genes = set(prog.index)
    allpool: Dict[int, List[str]] = {}
    for gene in expr_all.index:
        d = int(deciles.get(gene, -1))
        if d >= 0 and gene not in programme_genes:
            allpool.setdefault(d, []).append(gene)

    # How many genes each decile has to contribute.
    prog_expr = expr_all.reindex(prog.index).dropna()
    need = prog_expr.groupby(deciles.reindex(prog_expr.index)).size().to_dict()

    rng = np.random.default_rng(SEED)
    recs: List[dict] = []
    for rep in range(PARAMS.n_random_modules):
        picks: List[str] = []
        for d, k in need.items():
            cand = allpool.get(int(d), [])
            if len(cand) < k:
                continue
            picks += list(rng.choice(cand, size=k, replace=False))
        for gene in picks:
            recs.append({"module": f"R{rep:03d}", "gene": gene})

    # A decile that cannot supply its quota leaves every module short, and an
    # empty draw leaves no table at all. Both mean the expression reference is
    # too thin to match on, so they are reported instead of written out.
    if not recs:
        raise SystemExit(
            "no random module could be drawn: for every expression decile the "
            "candidate pool outside the programme is smaller than the programme's "
            f"own contribution ({dict(sorted((int(k), int(v)) for k, v in need.items()))})")

    rm = pd.DataFrame(recs, columns=["module", "gene"]).drop_duplicates()
    rm.to_csv(PATHS.tables / "random_modules_v1.tsv", sep="\t", index=False)
    print(f"random modules: {rm.module.nunique()} sets, median size "
          f"{int(rm.groupby('module').size().median())} "
          f"(programme {len(prog)} genes, matched on expression decile)")
    return rm


# --------------------------------------------------------------------------
# Step 2: lock the signatures
# --------------------------------------------------------------------------
def lock_signatures() -> pd.DataFrame:
    """Write the full and compact disease signatures and audit the compact one.

    Returns
    -------
    The compact signature table.

    Raises
    ------
    SystemExit
        If the compact signature is dominated by cell-cycle genes. Failing here
        is intentional: the incremental model would otherwise attribute a
        proliferation effect to the disease signature.
    """
    PATHS.ensure()
    sig = read_tsv("EC_signature_v1.tsv", index_col=None)
    sig = sig.drop_duplicates("gene")
    sig["abs"] = sig.logFC.abs()
    full = sig.sort_values("abs", ascending=False)
    full.to_csv(PATHS.tables / "EC_signature_v1_full.tsv", sep="\t",
                index=False, float_format="%.6g")

    # Compact form: the strongest genes on each side, then re-sorted by
    # direction and descending effect so that the table reads as two blocks.
    compact = pd.concat([
        full[full.direction == "UP"].nlargest(COMPACT_PER_DIRECTION, "abs"),
        full[full.direction == "DOWN"].nlargest(COMPACT_PER_DIRECTION, "abs"),
    ]).drop(columns=["abs"]).sort_values(["direction", "logFC"],
                                         ascending=[True, False])
    compact.to_csv(PATHS.tables / "EC_signature_v1_compact.tsv", sep="\t",
                   index=False, float_format="%.6g")

    print(f"disease signature: full {len(full)} "
          f"(UP {int((full.direction == 'UP').sum())} / "
          f"DOWN {int((full.direction == 'DOWN').sum())})")
    print(f"                   compact {len(compact)} "
          f"({COMPACT_PER_DIRECTION} per direction, used by the models)")

    # ---- proliferation confound audit -----------------------------------
    hit_c = sorted(set(CYCLING_GENES) & set(compact.gene))
    hit_f = sorted(set(CYCLING_GENES) & set(full.gene))
    print(f"\ncell-cycle genes present: compact {len(hit_c)}/{len(CYCLING_GENES)}, "
          f"full {len(hit_f)}/{len(CYCLING_GENES)}")
    print("  compact:", ", ".join(hit_c) if hit_c else "none")
    print("  full   :", ", ".join(hit_f) if hit_f else "none")
    verdict = ("compact signature is not proliferation-dominated"
               if len(hit_c) <= MAX_PROLIF_IN_COMPACT else
               "compact signature still carries many cell-cycle genes; the "
               "proliferation-adjusted delta-S is required")
    print("  verdict: " + verdict)
    if len(hit_c) > MAX_PROLIF_IN_COMPACT:
        raise SystemExit("proliferation confound: " + verdict)

    _write_manifest(["IR_program_v1.tsv", "EC_signature_v1_full.tsv",
                     "EC_signature_v1_compact.tsv"])
    return compact


def _write_manifest(files: Sequence[str]) -> List[str]:
    """Hash a list of tables and write the manifest.

    Parameters
    ----------
    files : file names inside ``PATHS.tables``.

    Returns
    -------
    The manifest lines, in ``sha256sum`` format.
    """
    lines = []
    for name in files:
        path = PATHS.tables / name
        if not path.exists():
            raise SystemExit(f"cannot hash missing table: {path}")
        lines.append(f"{sha256_of(path)}  {name}")
    (PATHS.tables / "W2_lock_manifest.txt").write_text("\n".join(lines) + "\n")
    print("\nmanifest:")
    print("\n".join(lines))
    return lines


# --------------------------------------------------------------------------
# Step 3: finalise and check against the independent cohort
# --------------------------------------------------------------------------
def finalise() -> None:
    """Normalise the adjusted signatures, verify reproduction, rewrite manifest.

    Returns
    -------
    None.

    Notes
    -----
    In the original layout the adjusted signatures were produced in a working
    directory and copied into the results tree here. Stage 02 now writes them
    straight into the tables directory, so source and destination coincide; the
    step is kept because it re-serialises the floats at a fixed precision, which
    is what makes the manifest stable across machines.
    """
    PATHS.ensure()
    prog = read_tsv("IR_program_v1.tsv")
    sig_full = read_tsv("EC_signature_v2_prolifadj_full.tsv", index_col=None)
    sig_comp = read_tsv("EC_signature_v2_prolifadj_compact.tsv", index_col=None)

    sig_full.to_csv(PATHS.tables / "EC_signature_v2_prolifadj_full.tsv",
                    sep="\t", index=False, float_format="%.6g")
    sig_comp.to_csv(PATHS.tables / "EC_signature_v2_prolifadj_compact.tsv",
                    sep="\t", index=False, float_format="%.6g")

    # ---- independent reproduction cohort --------------------------------
    rep = read_tsv("DE_GSE196033_tumor_vs_normal.tsv")
    for name, s in (("full", sig_full), ("compact", sig_comp)):
        ov = s[s.gene.isin(rep.index)]
        a = np.sign(ov.logFC.values)
        b = np.sign(rep.loc[ov.gene, "logFC"].values)
        agree = int((a == b).sum())
        p = binomtest_two_sided(agree, len(ov))
        same_sig = int(((a == b)
                        & (rep.loc[ov.gene, "adj.P.Val"].values
                           < PARAMS.program_de_padj)).sum())
        print(f"{name}: overlap {len(ov)}, same direction {agree} "
              f"({100 * agree / len(ov):.1f}%, binomial p={p:.1e}), "
              f"same direction and FDR<{PARAMS.program_de_padj} here: {same_sig}")

    _write_manifest(["IR_program_v1.tsv",
                     "EC_signature_v2_prolifadj_full.tsv",
                     "EC_signature_v2_prolifadj_compact.tsv",
                     "random_modules_v1.tsv"])

    # ---- overlap between the programme and the disease signature --------
    # Reported for the interpretation of the incremental value: an overlap of
    # zero would mean the two scores cannot be measuring the same biology.
    ir_genes = set(prog.index)
    for name, s in (("full", sig_full), ("compact", sig_comp)):
        ov = ir_genes & set(s.gene)
        print(f"\nprogramme and disease signature ({name}) share {len(ov)} genes: "
              f"{', '.join(sorted(ov)[:20])}")


def binomtest_two_sided(k: int, n: int) -> float:
    """Two-sided exact binomial p-value against p = 0.5.

    Parameters
    ----------
    k : number of successes.
    n : number of trials.

    Returns
    -------
    p-value.

    Notes
    -----
    Two-sided, unlike the concordance table in stage 01, because this check asks
    only whether the reproduction differs from chance at all -- the direction of
    the deviation is already known from the discovery cohort.
    """
    from scipy.stats import binomtest
    if n == 0:
        return float("nan")
    return binomtest(k, n, 0.5).pvalue


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------
def main() -> None:
    """Run the three steps in the order in which the objects are frozen."""
    print(f"tables directory: {PATHS.tables}\n")
    prog = build_programme()
    print(f"\nprogramme: {len(prog)} genes\n")
    lock_signatures()
    print()
    finalise()
    print("\ndone")


if __name__ == "__main__":
    main()
