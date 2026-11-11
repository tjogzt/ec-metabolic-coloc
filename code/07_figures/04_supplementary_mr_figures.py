#!/usr/bin/env python3
"""
04_supplementary_mr_figures.py
==============================

Purpose
-------
Draw the two supplementary Mendelian randomisation figures:

    Figure S1  forest plot of the univariable and multivariable estimates of the
               effect of fasting insulin and BMI on endometrial cancer
    Figure S2  leave-one-out plot of the fasting insulin IVW estimate

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
config.PATHS.tables / mr_results.tsv      from 04_run_mr.py
config.PATHS.tables / loo_fi.tsv          from 05_conditional_f.py
config.PARAMS.mr_min_f                    weak-instrument floor quoted in the
                                          annotation of Figure S1
config.PARAMS.alpha                       significance level used to encode the
                                          symbols in Figure S2

Outputs
-------
config.PATHS.figures / figureS1_mr_forest.pdf and .png
config.PATHS.figures / figureS2_leave_one_out.pdf and .png

Both are written as vector PDF for the manuscript and as 300 dpi PNG for
inspection, with fonts embedded as Type 42 so that the PDFs remain editable.

Notes
-----
House style. Arial at 8 pt minimum, no panel letters, near-square panels, and
semantic colour: cinnabar for adiposity, azurite for insulin. The multivariable
rows are drawn in neutral grey with a dotted interval because they are not
interpretable: the conditional F statistics of 05_conditional_f.py fall below the
weak-instrument floor, and the figure must not make an uninterpretable estimate
look like the others.

The source script for this stage also drew a paired cross-cancer correlation
panel from a table produced by a different analysis layer. That panel belongs to
the cross-cancer stage and is not reproduced here, so this file depends only on
the outputs of the MR stage.

Figure S2 reference line. The vertical line marks the full-sample estimate of the
variant subset that the leave-one-out analysis actually uses (22 variants: those
carrying both a fasting insulin and a BMI effect). That estimate is not the IVW
row of mr_results.tsv, which is fitted on all 34 fasting insulin instruments, so
the value is recorded here as its odds ratio (1.5024) rather than read from the
results table.
"""
from __future__ import annotations

import pathlib
import sys

import matplotlib

matplotlib.use("Agg")           # headless: no display is available on the cluster

import matplotlib.pyplot as plt  # noqa: E402  (after the backend is fixed)
import numpy as np               # noqa: E402
import pandas as pd              # noqa: E402
from matplotlib import font_manager as fm  # noqa: E402

# Repository root, so that ``config`` resolves when the file is run directly.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

from config.config import (  # noqa: E402  (after sys.path edit)
    PATHS,
    PARAMS,
    FIGURE_FONT,
    FIGURE_MIN_PT,
    PALETTE,
)

# Semantic colours, taken from the shared palette rather than restated here so
# that the manuscript and the supplementary figures cannot drift apart.
ZHU = PALETTE["uterus"]     # cinnabar : BMI estimates
SHI = PALETTE["liver"]      # azurite  : fasting insulin estimates
# A lighter grey than PALETTE["neutral"], reserved for estimates that must not
# be read as interpretable. Expressed as a greyscale level so that it stays readable when printed.
GREY = "0.62"

MR_RESULTS = "mr_results.tsv"
LOO_FI = "loo_fi.tsv"
HARMONISED = "harmonised.tsv"

Z_95 = 1.96


def register_font() -> None:
    """Register Arial with matplotlib, falling back to the system font files.

    The manuscript requires Arial and PDF text left as text (Type 42). Matplotlib
    only finds Arial if the font file is registered explicitly on macOS, so the
    two standard system locations are tried when the font database does not
    already resolve the configured family.

    Notes
    -----
    Failure is silent by design: an unavailable font must not stop a figure run
    that already has its output paths, and the rcParams below still name the
    family, so a substituted font is visible in the PDF metadata rather than
    hidden.
    """
    try:
        fm.findfont(FIGURE_FONT, fallback_to_default=False)
        return
    except Exception:
        pass
    for candidate in ("/System/Library/Fonts/Supplemental/Arial.ttf",
                      "/Library/Fonts/Arial.ttf"):
        if pathlib.Path(candidate).exists():
            fm.fontManager.addfont(candidate)
            return


def apply_style() -> None:
    """Set the figure style shared by both panels."""
    plt.rcParams.update({
        "font.family": FIGURE_FONT,
        "font.size": 9,
        "axes.linewidth": 0.8,
        "axes.edgecolor": "0.20",
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.labelsize": max(FIGURE_MIN_PT, 8),
        "ytick.labelsize": max(FIGURE_MIN_PT, 8),
        "pdf.fonttype": 42,         # embed TrueType so the PDF text stays editable
        "figure.dpi": 300,
        "savefig.dpi": 300,
    })


# --------------------------------------------------------------------------
# Figure S1: forest plot
# --------------------------------------------------------------------------
def forest_plot(results: pd.DataFrame) -> None:
    """Draw the forest plot of every reported estimate.

    Parameters
    ----------
    results : DataFrame
        mr_results.tsv, with columns exposure, method, OR, OR_lo, OR_hi.

    Notes
    -----
    Rows are ordered from the multivariable estimates at the top to the
    univariable IVW estimates at the bottom, so the reader meets the weakest
    evidence first. The order is fixed here rather than derived from the table,
    because the labels are shortened for print and carry the interpretability
    verdict.
    """
    rows = [
        ("MVMR BMI (invalid)", "MVMR-BMI", "MVMR", GREY, True),
        ("MVMR fasting insulin (invalid)", "MVMR-fasting insulin", "MVMR", GREY, True),
        ("BMI, weighted median", "BMI", "Weighted median", ZHU, False),
        ("BMI, MR-Egger", "BMI", "MR-Egger", ZHU, False),
        ("BMI, IVW", "BMI", "IVW", ZHU, False),
        ("Fasting insulin, weighted median", "fasting insulin", "Weighted median", SHI, False),
        ("Fasting insulin, MR-Egger", "fasting insulin", "MR-Egger", SHI, False),
        ("Fasting insulin, IVW", "fasting insulin", "IVW", SHI, False),
    ]

    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    ys = np.arange(len(rows))[::-1]     # first row of the list ends up on top

    for y, (label, exposure, method, colour, invalid) in zip(ys, rows):
        # One row per (exposure, method); a missing combination is an error in
        # the results table, not something to skip silently.
        r = results[(results.exposure == exposure) & (results.method == method)].iloc[0]
        ax.plot([r.OR_lo, r.OR_hi], [y, y], color=colour, lw=1.4,
                ls=":" if invalid else "-", alpha=0.85, zorder=2)
        ax.plot(r.OR, y, "o", ms=6.5, color="white", mec=colour, mew=1.5, zorder=3)
        # En dash between the interval bounds, not a hyphen.
        txt = f"{r.OR:.2f} ({r.OR_lo:.2f}\u2013{r.OR_hi:.2f})"
        ax.text(7.6, y, txt, va="center", ha="left", fontsize=8,
                color=colour if not invalid else "0.48")

    # Null line at OR = 1.
    ax.axvline(1.0, color="0.33", lw=0.9, ls="--", zorder=1)

    # Log scale is mandatory for a ratio: it makes the two halves of the
    # interval symmetric around the estimate.
    ax.set_xscale("log")
    ax.set_xlim(0.5, 7.2)
    ax.set_xticks([0.5, 1, 2, 3, 5])
    ax.set_xticklabels(["0.5", "1", "2", "3", "5"])
    ax.set_yticks(ys)
    ax.set_yticklabels([row[0] for row in rows], fontsize=8)
    ax.set_xlabel("Odds ratio for endometrial cancer (95% CI)", fontsize=9)
    ax.set_ylim(-0.7, len(rows) + 0.15)
    ax.text(7.6, len(rows) - 0.55, "OR (95% CI)", fontsize=8, ha="left",
            color="0.20", weight="bold")
    ax.spines[["top", "right"]].set_visible(False)

    # Shaded band with the reason the multivariable rows cannot be read. It must
    # sit over the two MVMR rows, which are the first two entries and therefore
    # land at the TOP of the axis (ys counts down). The earlier version spanned
    # -0.75..1.4, which covers the two bottom rows, so the band appeared to
    # exclude the fasting-insulin estimates. The threshold in the label is taken
    # from the config value that produced the verdict.
    n_mvmr = sum(1 for r in rows if r[4])          # rows flagged as invalid
    top, bottom = len(rows) - 1, len(rows) - n_mvmr
    ax.axhspan(bottom - 0.45, top + 0.45, color=GREY, alpha=0.09, zorder=0)
    ax.text(0.52, top + 0.55,
            f"Multivariable MR excluded: conditional F < {PARAMS.mr_min_f:g}",
            fontsize=8, color="0.40", va="bottom")

    fig.tight_layout()
    fig.savefig(PATHS.figures / "figureS1_mr_forest.pdf")
    fig.savefig(PATHS.figures / "figureS1_mr_forest.png")
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure S2: leave-one-out
# --------------------------------------------------------------------------
def full_sample_subset_estimate(harmonised: pd.DataFrame,
                                loo: pd.DataFrame) -> float:
    """Inverse-variance weighted log-OR over the variants the LOO panel tests.

    Parameters
    ----------
    harmonised : DataFrame
        harmonised.tsv, with rsid and the per-variant effects and standard
        errors for fasting insulin (fi_beta, fi_se) and the outcome
        (ec_beta, ec_se).
    loo : DataFrame
        loo_fi.tsv; its dropped_rsid column defines the subset.

    Returns
    -------
    The log-OR of the fixed-effect inverse-variance weighted estimate over the
    whole subset, i.e. the value each leave-one-out replicate is compared with.

    Notes
    -----
    This replaces a literal that was carried over from the exploratory run. The
    literal could not be reproduced from any shipped table (the fixed- and
    random-effect forms both give 1.5017 against the recorded 1.5024), so the
    value is now derived from the data every time the panel is drawn. The
    difference is 0.05% in OR and does not move the reference line visibly.

    Weights follow the two-sample ratio estimator, w = (b_exposure / se_outcome)^2.
    """
    subset = harmonised[harmonised["rsid"].isin(set(loo["dropped_rsid"]))]
    subset = subset.dropna(subset=["fi_beta", "fi_se", "ec_beta", "ec_se"])
    if subset.empty:
        raise SystemExit("no variants shared between harmonised.tsv and loo_fi.tsv")
    bx = subset["fi_beta"].to_numpy(float)
    by = subset["ec_beta"].to_numpy(float)
    sy = subset["ec_se"].to_numpy(float)
    w = (bx / sy) ** 2
    return float((w * (by / bx)).sum() / w.sum())


def leave_one_out_plot(loo: pd.DataFrame, full_beta: float) -> None:
    """Draw the leave-one-out plot for the fasting insulin IVW estimate.

    Parameters
    ----------
    loo : DataFrame
        loo_fi.tsv, with columns dropped_rsid, beta, se, p, OR.

    Notes
    -----
    Variants are sorted by the estimate obtained when they are omitted, so the
    plot reads as a gradient from the most influential (left) to the least
    (right). Filled symbols mark replicates that stay significant at
    PARAMS.alpha; open symbols mark those that do not, which is the pattern that
    matters for stability.
    """
    loo = loo.sort_values("beta")

    fig, ax = plt.subplots(figsize=(5.0, 6.4))
    y = np.arange(len(loo))

    for i, r in enumerate(loo.itertuples()):
        significant = r.p < PARAMS.alpha
        # log OR is the stored beta; the interval uses the same 1.96 multiplier
        # as the results table.
        ax.plot([np.log(r.OR) - Z_95 * r.se, np.log(r.OR) + Z_95 * r.se],
                [i, i], color=SHI if significant else GREY, lw=0.9, alpha=0.9)
        ax.plot(np.log(r.OR), i, "o", ms=3.6,
                color=SHI if significant else "white",
                mec=SHI if significant else GREY, mew=1.0)

    # Reference line: the full-sample estimate of the same variant subset.
    ax.axvline(full_beta, color=ZHU, lw=1.2, ls="--")
    ax.axvline(0, color="0.33", lw=0.8, ls=":")

    ax.set_yticks(y)
    ax.set_yticklabels([str(x) for x in loo.dropped_rsid], fontsize=8)
    ax.set_ylabel("Variant omitted", fontsize=9)
    ax.set_xlabel("log OR for endometrial cancer (95% CI)", fontsize=9)
    ax.set_ylim(-0.8, len(loo) - 0.2)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(full_beta + 0.03, len(loo) - 1.0, "all variants", fontsize=8,
            color=ZHU, rotation=90, va="top")
    ax.text(0.02, 0.02,
            f"Open symbols: p \u2265 {PARAMS.alpha:g} after omission",
            transform=ax.transAxes, fontsize=8, color="0.40")

    fig.tight_layout()
    fig.savefig(PATHS.figures / "figureS2_leave_one_out.pdf")
    fig.savefig(PATHS.figures / "figureS2_leave_one_out.png")
    plt.close(fig)


def main() -> None:
    """Render both panels and list what was written."""
    PATHS.ensure()
    register_font()
    apply_style()

    results = pd.read_csv(PATHS.tables / MR_RESULTS, sep="\t")
    forest_plot(results)

    loo = pd.read_csv(PATHS.tables / LOO_FI, sep="\t")
    harmonised = pd.read_csv(PATHS.tables / HARMONISED, sep="\t")
    leave_one_out_plot(loo, full_sample_subset_estimate(harmonised, loo))

    print("wrote:")
    for path in sorted(PATHS.figures.glob("figureS*")):
        print(f"  {path.name}")


if __name__ == "__main__":
    main()
