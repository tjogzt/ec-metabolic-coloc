"""
01_concordance_pipeline.py
==========================

Purpose
-------
Measure how consistently an insulin-resistance transcriptional programme is
reproduced, in direction, across the three GEO cohorts, and report the three
statistics the manuscript quotes side by side:

    (1) sign concordance over every gene shared with the cohort -- the original
        statistic, kept for comparison;
    (2) sign concordance restricted to genes that actually respond in the
        validation cohort (nominal p below the significance level) -- the
        primary statistic;
    (3) Spearman correlation of the two log fold-change vectors, which is not
        affected by dichotomising the sign.

Statistic (1) is dominated by genes that carry no measurable signal: their sign
is close to a coin flip, so the estimate is pulled towards 50% by construction.
Statistic (2) asks the question the claim needs -- among the genes the cohort
does detect, is the direction right more often than chance?

Two audits that were separate scripts in the first pass are merged here because
they recompute the same quantity and must never drift apart:

    --check      recomputes the all-gene concordance with independently detected
                 column names (the appendix audit);
    --dilution   restricts the clamp cohort to strata of increasingly strict
                 response and asks whether the all-gene estimate is diluted.

Both are read-only: they print and write nothing.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
PATHS.tables/metabolic_program_v2.tsv
    the locked programme; ``logFC_adipose_obese`` is the reference effect every
    cohort is compared against.
PATHS.tables/DE_<cohort>.tsv
    limma output for the three concordance cohorts: one row per gene, columns
    logFC, P.Value, adj.P.Val. These are NOT produced by stage 02, which fits
    the insulin-response and tumour contrasts; they come from the chronic-IR DE
    run (code/01_programme/w4_de.R for GSE15773 and GSE18732, w4_gse13070.R for
    the clamp cohort) and must be placed in the tables directory before this
    script runs.
config.COHORTS
    the cohort set. Membership is asserted against COHORT_ORDER below, so a
    cohort added to config fails loudly here instead of being silently skipped.
config.PARAMS.alpha
    nominal detection level used to define the responsive subset.
config.PARAMS.program_fdr
    FDR level used for the two count columns.

Outputs
-------
PATHS.tables/concordance_by_cohort.tsv, one row per cohort, columns in this
order:

    cohort, n_overlap, conc_all, ci_all_lo, ci_all_hi, p_all, n_responsive,
    conc_resp, ci_resp_lo, ci_resp_hi, p_resp, rho, p_rho,
    fdr_sig_in_program, fdr_sig_genomewide

Proportions are stored in percent (0-100), matching the table the figures
consume. Rows follow the manuscript order, clamp cohort first.

Usage
-----
    python code/03_burden/01_concordance_pipeline.py
    python code/03_burden/01_concordance_pipeline.py --check
    python code/03_burden/01_concordance_pipeline.py --dilution

Notes
-----
Confidence intervals are Wilson score intervals rather than normal
approximations: the observed proportions sit near 50% with a few hundred genes,
which is precisely where the Wald interval is worst.

The one-sided binomial test is used throughout -- the hypothesis is directional
("more concordant than chance"), not merely "different from chance".
"""
from __future__ import annotations

import argparse
import pathlib
import sys
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr

# Repository root, so that `config` resolves when the file is run as a script.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

from config.config import COHORTS, PARAMS, PATHS

# --------------------------------------------------------------------------
# Module constants
# --------------------------------------------------------------------------
# The programme table is a shipped, hash-locked artefact; its file name is part
# of the analysis definition, not a tunable.
PROGRAMME_FILE = "metabolic_program_v2.tsv"

# Column of the programme table used as the reference effect. It is the adipose
# obesity contrast, i.e. the tissue in which the programme was derived.
REFERENCE_COLUMN = "logFC_adipose_obese"

# Row order of the output table, clamp cohort first, as in the manuscript.
# Membership is governed by config.COHORTS and checked in _cohort_plan().
COHORT_ORDER: Tuple[str, ...] = ("GSE13070", "GSE18732", "GSE15773")

# Per-cohort strings: the manuscript label and the DE file produced by stage 02.
COHORT_LABEL: Dict[str, str] = {
    "GSE13070": "GSE13070 muscle, clamp-defined IR vs IS",
    "GSE18732": "GSE18732 muscle, type 2 diabetes vs NGT",
    "GSE15773": "GSE15773 adipose, IR vs IS",
}
COHORT_DE_FILE: Dict[str, str] = {
    "GSE13070": "DE_GSE13070_muscle_baseline_IR.tsv",
    "GSE18732": "DE_GSE18732_muscle_DM_state.tsv",
    "GSE15773": "DE_GSE15773_adipose_IR_state.tsv",
}

# Two-sided 95% normal quantile used by the Wilson interval.
Z_95 = 1.96

# Pre-specified small-effect boundary on the log2 fold change (0.3 log2 units,
# roughly a 23% change). Used by the dilution audit on the clamp contrast.
SMALL_EFFECT_CUT = 0.3

# Stricter nominal level for the "clearly responsive" stratum of the dilution
# audit; deliberately tighter than PARAMS.alpha so that the stratum is small
# and uncontaminated by borderline calls.
DILUTION_NOMINAL = 1e-3


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------
def _cohort_plan() -> List[Tuple[str, str, str]]:
    """Return (accession, label, DE file name) in manuscript order.

    Returns
    -------
    list of triples, one per cohort, ordered as in ``COHORT_ORDER``.

    Raises
    ------
    SystemExit
        If the configured cohort set and the wired-up cohort set disagree, in
        either direction. Silent skipping is the failure mode this prevents: a
        cohort that config knows about but this table does not would drop out of
        the figure without any warning.
    """
    configured = {c.gse for c in COHORTS}
    wired = set(COHORT_ORDER)
    if configured != wired:
        missing = sorted(configured - wired)
        extra = sorted(wired - configured)
        raise SystemExit(
            f"cohort wiring mismatch: in config.COHORTS but not here {missing}; "
            f"here but not in config.COHORTS {extra}")
    return [(gse, COHORT_LABEL[gse], COHORT_DE_FILE[gse]) for gse in COHORT_ORDER]


def _load_programme() -> pd.DataFrame:
    """Read the locked programme table.

    Returns
    -------
    DataFrame indexed by gene symbol with at least ``REFERENCE_COLUMN``.

    Raises
    ------
    SystemExit
        If the file or the reference column is missing.
    """
    path = PATHS.tables / PROGRAMME_FILE
    if not path.exists():
        raise SystemExit(f"missing programme table: {path}\n"
                         f"  run code/03_burden/05_build_programme.py first")
    prog = pd.read_csv(path, sep="\t", index_col=0)
    if REFERENCE_COLUMN not in prog.columns:
        raise SystemExit(f"{path.name}: no '{REFERENCE_COLUMN}' column; "
                         f"found {list(prog.columns)}")
    return prog


def _read_de(name: str) -> Tuple[pd.DataFrame, pathlib.Path]:
    """Read one limma DE table written by stage 02.

    Parameters
    ----------
    name : file name inside PATHS.tables.

    Returns
    -------
    (DataFrame indexed by gene symbol, resolved path).
    """
    path = PATHS.tables / name
    if not path.exists():
        raise SystemExit(
            f"missing DE table: {path}\n"
            f"  the concordance cohorts are not fitted by stage 02; copy the "
            f"limma output (columns logFC, P.Value, adj.P.Val) produced by "
            f"code/01_programme/w4_de.R and w4_gse13070.R into the tables directory")
    return pd.read_csv(path, sep="\t", index_col=0), path


def _detect_column(columns: Sequence[str], predicate, what: str,
                   path: pathlib.Path) -> str:
    """Pick the first column honouring ``predicate``, or fail with the layout.

    Parameters
    ----------
    columns : column names of the table being audited.
    predicate : callable applied to each name.
    what : human-readable description of the expected column.
    path : source file, quoted in the error message.

    Notes
    -----
    The audits deliberately re-derive column names instead of reusing the names
    from the main pass. That is the point of the audit: if stage 02 changes a
    column heading, the main pass and the audit stop agreeing.
    """
    hits = [c for c in columns if predicate(c)]
    if not hits:
        raise SystemExit(f"{path.name}: no {what} column among {list(columns)}")
    return hits[0]


def wilson(k: int, n: int, z: float = Z_95) -> Tuple[float, float]:
    """Wilson score interval for a binomial proportion.

    Parameters
    ----------
    k : number of concordant genes.
    n : number of genes compared.
    z : normal quantile, 1.96 by default (two-sided 95%).

    Returns
    -------
    (lower, upper) bounds as proportions, not percentages.

    Notes
    -----
    With n in the hundreds and p near 0.5 the Wald interval is anti-conservative,
    and it degenerates to zero width when k is 0 or n; the score interval has
    neither problem.
    """
    p = k / n
    d = 1 + z ** 2 / n                                  # variance inflation
    centre = (p + z ** 2 / (2 * n)) / d                 # shrunken centre
    half = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / d
    return centre - half, centre + half


# --------------------------------------------------------------------------
# Main table
# --------------------------------------------------------------------------
def build_concordance_table() -> pd.DataFrame:
    """Compute the per-cohort concordance table and write it to the tables dir.

    Returns
    -------
    The table that was written, in memory.

    Notes
    -----
    Every row needs three quantities computed on one joined matrix: the
    all-gene sign concordance, the concordance among genes that respond in the
    cohort, and the rank correlation of the effects. Joining first keeps the
    three on exactly the same gene set; computing them on differently filtered
    subsets is the mistake this layout prevents.
    """
    PATHS.ensure()
    prog = _load_programme()
    plan = _cohort_plan()

    rows: List[dict] = []
    print(f"{'cohort':46s} {'all genes':>16s} {'responsive (p<alpha)':>22s} "
          f"{'rho':>18s}")

    for gse, label, de_name in plan:
        de, _ = _read_de(de_name)

        # Inner join: genes present in both the programme and this cohort.
        m = prog[[REFERENCE_COLUMN]].join(
            de[["logFC", "P.Value", "adj.P.Val"]], how="inner").dropna()
        m["same"] = np.sign(m[REFERENCE_COLUMN]) == np.sign(m.logFC)

        # (1) all shared genes
        k, n = int(m.same.sum()), len(m)
        if n == 0:
            raise SystemExit(f"{gse}: no gene shared with the programme")
        lo, hi = wilson(k, n)
        p_all = binomtest(k, n, 0.5, alternative="greater").pvalue

        # (2) genes that respond in the validation cohort
        resp = m[m["P.Value"] < PARAMS.alpha]
        k2, n2 = int(resp.same.sum()), len(resp)
        if n2:
            lo2, hi2 = wilson(k2, n2)
            p_resp = binomtest(k2, n2, 0.5, alternative="greater").pvalue
            conc_resp = 100 * k2 / n2
            resp_text = f"{conc_resp:5.1f}% ({k2}/{n2})"
        else:
            # No responsive gene: report the row as undefined rather than
            # dividing by zero, and keep the cohort in the table.
            lo2 = hi2 = p_resp = conc_resp = np.nan
            resp_text = f"n/a ({k2}/{n2})"

        # (3) rank correlation of the two effect vectors
        rho, p_rho = spearmanr(m[REFERENCE_COLUMN], m.logFC)

        print(f"{label:46s} {100 * k / n:5.1f}% ({k}/{n}) {resp_text:>14s} "
              f"rho={rho:+.3f} (p={p_rho:.1e})")

        rows.append(dict(
            cohort=label,
            n_overlap=n,
            conc_all=100 * k / n,
            ci_all_lo=100 * lo,
            ci_all_hi=100 * hi,
            p_all=p_all,
            n_responsive=n2,
            conc_resp=conc_resp,
            ci_resp_lo=100 * lo2 if n2 else np.nan,
            ci_resp_hi=100 * hi2 if n2 else np.nan,
            p_resp=p_resp,
            rho=rho,
            p_rho=p_rho,
            # Significant genes inside the overlap, and across the whole cohort
            fdr_sig_in_program=int((m["adj.P.Val"] < PARAMS.program_fdr).sum()),
            fdr_sig_genomewide=int((de["adj.P.Val"] < PARAMS.program_fdr).sum()),
        ))

    out = pd.DataFrame(rows)
    dest = PATHS.tables / "concordance_by_cohort.tsv"
    out.to_csv(dest, sep="\t", index=False)
    print(f"\nwritten {dest}")

    # Compact recap: the four columns a reader reads first.
    print(out[["cohort", "conc_all", "conc_resp", "rho",
               "fdr_sig_in_program", "fdr_sig_genomewide"]].to_string(index=False))
    return out


# --------------------------------------------------------------------------
# Audit 1: independent re-derivation of the all-gene concordance
# --------------------------------------------------------------------------
def verify_concordance() -> None:
    """Recompute the all-gene sign concordance with detected column names.

    Returns
    -------
    None. Everything is printed.

    Notes
    -----
    The fold-change column is located by substring, not by name, so the check
    still runs against a DE table written by a different limma call. The counts
    must equal ``n_overlap`` and ``conc_all`` in the main table; the shipped
    package reports both so that the two can be compared by eye.
    """
    prog = _load_programme()
    print(f"programme: {len(prog)} genes; columns {list(prog.columns)}")

    for gse, label, de_name in _cohort_plan():
        de, path = _read_de(de_name)
        fc = _detect_column(de.columns, lambda c: "logFC" in c,
                            "log fold-change", path)
        d = de[[fc]].dropna()

        common = prog.index.intersection(d.index)
        same = int((np.sign(d.loc[common, fc])
                    == np.sign(prog.loc[common, REFERENCE_COLUMN])).sum())
        n = len(common)
        if n == 0:
            print(f"{gse:10s} overlap 0; skipped")
            continue
        p = binomtest(same, n, 0.5, alternative="greater").pvalue
        print(f"{label:46s} overlap {n:5d}  same sign {same:5d} "
              f"({100 * same / n:5.1f}%)  one-sided p = {p:.2e}")


# --------------------------------------------------------------------------
# Audit 2: is the all-gene estimate diluted by non-responding genes?
# --------------------------------------------------------------------------
def check_dilution() -> None:
    """Stratify the clamp cohort by response strength and re-measure concordance.

    Returns
    -------
    None. Everything is printed.

    Notes
    -----
    If the all-gene estimate is diluted by genes with no measurable signal, the
    concordance must rise monotonically as the stratum tightens. Restricting to
    genes with both a small FDR and a non-trivial effect is the strictest
    stratum available without conditioning on the direction itself, which would
    make the test circular.

    Only the clamp cohort is audited: it is the reference standard for insulin
    resistance, so a dilution argument has to be settled on it rather than on a
    proxy cohort.
    """
    gse = "GSE13070"
    prog = _load_programme()
    de, path = _read_de(COHORT_DE_FILE[gse])

    fc = _detect_column(de.columns, lambda c: "logFC" in c, "log fold-change", path)
    fdr = _detect_column(de.columns,
                         lambda c: "adj" in c.lower() or "FDR" in c,
                         "adjusted p-value", path)
    pcol = _detect_column(de.columns,
                          lambda c: c.lower() in ("p.value", "pvalue"),
                          "nominal p-value", path)

    m = prog[[REFERENCE_COLUMN]].join(de[[fc, fdr, pcol]], how="inner").dropna()
    m["same"] = np.sign(m[REFERENCE_COLUMN]) == np.sign(m[fc])

    overall = binomtest(int(m.same.sum()), len(m), 0.5,
                        alternative="greater").pvalue
    print(f"{gse}: overlap {len(m)}; all genes concordant "
          f"{int(m.same.sum())} ({100 * m.same.mean():.2f}%); "
          f"one-sided p = {overall:.4f}")

    strata = (
        (f"FDR < {PARAMS.program_fdr:g}", m[m[fdr] < PARAMS.program_fdr]),
        (f"nominal p < {DILUTION_NOMINAL:g}", m[m[pcol] < DILUTION_NOMINAL]),
        (f"|logFC| >= {SMALL_EFFECT_CUT:g}", m[m[fc].abs() >= SMALL_EFFECT_CUT]),
        (f"FDR < {PARAMS.program_fdr:g} and |logFC| >= {SMALL_EFFECT_CUT:g}",
         m[(m[fdr] < PARAMS.program_fdr) & (m[fc].abs() >= SMALL_EFFECT_CUT)]),
    )
    for label, sub in strata:
        k, n = int(sub.same.sum()), len(sub)
        if n == 0:
            print(f"{label:34s} n=0")
            continue
        p = binomtest(k, n, 0.5, alternative="greater").pvalue
        print(f"{label:34s} {k}/{n} = {100 * k / n:5.1f}%   one-sided p = {p:.4g}")

    rho, pr = spearmanr(m[REFERENCE_COLUMN], m[fc])
    print(f"\nall {len(m)} genes, log fold-change rank correlation "
          f"rho = {rho:+.3f}, p = {pr:.4g}")


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------
def main() -> None:
    """Parse arguments and dispatch to the table builder or one of the audits."""
    ap = argparse.ArgumentParser(
        description="Cohort concordance for the insulin-resistance programme.")
    ap.add_argument("--check", action="store_true",
                    help="re-derive the all-gene concordance (no output written)")
    ap.add_argument("--dilution", action="store_true",
                    help="stratify the clamp cohort by response strength "
                         "(no output written)")
    args = ap.parse_args()

    if args.check:
        verify_concordance()
        check_dilution()          # both audits always run together
        return
    if args.dilution:
        check_dilution()
        return
    build_concordance_table()


if __name__ == "__main__":
    main()
