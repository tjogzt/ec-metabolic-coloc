#!/usr/bin/env python3
"""
02_figure1_design.py
====================

Purpose
-------
Draw Figure 1, the study-design schematic: four evidence bands (programme
definition, colocalisation, bulk proxy validity, single-cell projection), each
with a headline number, a criterion line and a source tag, plus a terminal row
carrying the finding and the boundary of what may be claimed.

Created : 2026-11-11
Version : 2.1.0 (cohort correlations are read from the result table)

Inputs
------
config/config.py                          PALETTE, FIGURE_FONT, FIGURE_MIN_PT,
                                          PATHS.figures
results/tables/concordance_by_cohort.tsv  rank correlations for the three cohorts

The three Layer 3 correlations are read from the concordance table at run time
and checked against the values the manuscript reports. The remaining schematic
labels (programme size, colocalisation counts, single-cell z) restate values
from other result tables and are not derived here.

Outputs
-------
results/figures/fig1_design.pdf
results/figures/fig1_design.png

Notes
-----
* Colours are read from PALETTE by semantic role and are never spelled out as
  hex values here: the red band is the target organ / gold-standard phenotype,
  the blue band is the metabolic tissue and the green band is muscle. Band
  washes are the same palette colours blended towards white, so a palette
  change carries through the figure instead of leaving stale tints.
* Every size passed to matplotlib is clamped to FIGURE_MIN_PT, so the schematic
  still clears 8 pt at the 6.2 in manuscript placement.
* The figure is drawn in a 100 x 80 data coordinate system with the axis off;
  all placement below is in those units, which is what keeps the text inside
  its band at every font size.
* Band text is wrapped by hand: no line may overflow the band rectangle, and
  the right-hand badge is right-aligned to the canvas edge.
"""

from __future__ import annotations

import os
import pathlib
import sys

# Repository root is two levels above this file (code/07_figures/ -> repo root).
# ECIR_ROOT is pinned from the file location rather than inherited from the
# shell, so a stale variable cannot redirect the figure into another tree;
# ECIR_DATA and ECIR_RESULTS remain free to move data and outputs.
_ROOT = pathlib.Path(__file__).resolve().parents[2]
os.environ["ECIR_ROOT"] = str(_ROOT)
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from config import config as CFG  # noqa: E402

import pandas as pd  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402

# ------------------------------------------------------------------ style ----
plt.rcParams.update({
    "font.family": CFG.FIGURE_FONT,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def pt(size: float) -> float:
    """Clamp a requested type size to the configured minimum."""
    return max(float(size), float(CFG.FIGURE_MIN_PT))


def wash(colour: str, strength: float = 0.92) -> tuple:
    """Blend a palette colour towards white.

    ``strength`` is the white fraction, so 0.92 gives the pale band fill that
    separates adjacent bands without introducing a new hue.
    """
    r, g, b = to_rgb(colour)
    return (r + (1 - r) * strength, g + (1 - g) * strength, b + (1 - b) * strength)


# Semantic roles, resolved once from the config palette.
SHA = CFG.PALETTE["uterus"]     # red  : target organ / gold-standard phenotype
SHI = CFG.PALETTE["liver"]      # blue : metabolic tissue
DIAN = CFG.PALETTE["muscle"]    # green: skeletal muscle
GREY = CFG.PALETTE["neutral"]   # neutral : reference / boundary

# Structural greys belong to typography, not to the data palette: they are given
# as grey levels so no second colour system is introduced.
INK = "0.10"
BODY = "0.20"
MUTED = "0.33"
RULE = "0.78"
PANEL = "0.97"

TINT = {"sha": wash(SHA), "dian": wash(DIAN)}

# Concordance-table rows in figure order. Each accession maps to the label shown
# under the stem, the x position of the stem in data units, and the stem colour.
COHORT_KEYS: dict[str, tuple[str, float, str]] = {
    "GSE15773": ("adipose", 74.0, SHA),
    "GSE18732": ("diabetes", 83.0, SHA),
    "GSE13070": ("clamp", 92.0, GREY),
}

# Values the schematic must reproduce. The run stops if the table drifts, so a
# changed estimate cannot be drawn under an unchanged label.
EXPECTED_RHO = {"adipose": 0.313, "diabetes": 0.202, "clamp": 0.098}
EXPECTED_LABEL = {"adipose": "+0.31", "diabetes": "+0.20", "clamp": "+0.10"}


def read_rank_correlations() -> dict[str, float]:
    """Read the three cohort rank correlations from the concordance table."""
    path = CFG.PATHS.tables / "concordance_by_cohort.tsv"
    if not path.exists():
        raise SystemExit(f"missing concordance table: {path}")
    frame = pd.read_csv(path, sep="\t")
    values: dict[str, float] = {}
    for accession, (name, _x, _colour) in COHORT_KEYS.items():
        hit = frame[frame["cohort"].astype(str).str.startswith(accession)]
        if len(hit) != 1:
            raise SystemExit(
                f"expected one row for {accession} in {path}, found {len(hit)}"
            )
        values[name] = float(hit["rho"].iloc[0])
    return values

def main() -> None:
    """Render and save the Figure 1 schematic."""
    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 80)
    ax.axis("off")

    TAG_W = 15.0
    SPINE_X, SPINE_W = 16.6, 1.6
    TXT_X = 19.6

    rho = read_rank_correlations()

    # y origin, height, headline, band name, spine colour, wash, criteria, badge, source
    BANDS = [
        (60, 13, "862", "Layer 1  Programme definition", SHA, TINT["sha"],
         ["FDR<0.05, |logFC|≥0.3, endometrial expression gate; 608 up / 254 down",
          "preregistered, hashed, and carried unchanged into Layers 3 and 4"],
         "genes in programme", "Discovery\ncontrast"),
        (45, 13, "3 / 503", "Layer 2  Colocalisation", SHA, TINT["sha"],
         ["Gene-locus tests across 19 loci and five tissues; PP.H4>0.8 threshold",
          "73 of 95 locus-tissue cells testable; one signal survives a change of dataset"],
         "two at 17q21.32", "Risk GWAS\n+ cis-eQTL"),
        (30, 13, f"+{rho['adipose']:.2f}", "Layer 3  Bulk proxy validity", DIAN, TINT["dian"],
         ["Rank correlation with three independent cohorts",
          "adipose strong, diabetes intermediate, clamp near zero"],
         None, "(from Layer 1)"),
        (15, 13, "−3.8", "Layer 4  Single-cell projection", DIAN, TINT["dian"],
         ["11 samples, 141 316 cells; tumour vs normal epithelium, sample-level test",
          "against 100 matched random modules"],
         None, "(from Layer 1)"),
    ]

    for y0, h, head, label, col, tint, crit, res, src in BANDS:
        ax.add_patch(Rectangle((TAG_W + 0.8, y0), 100 - TAG_W - 1.3, h,
                               facecolor=tint, edgecolor="none", zorder=0))
        ax.add_patch(Rectangle((SPINE_X, y0), SPINE_W, h, facecolor=col,
                               edgecolor="none", zorder=2))
        ax.add_patch(FancyBboxPatch((0.8, y0 + h / 2 - 2.5), TAG_W - 2.0, 5.0,
                                    boxstyle="round,pad=0.22,rounding_size=0.9",
                                    linewidth=0.8, edgecolor=GREY,
                                    facecolor="white", zorder=2))
        ax.text(0.8 + (TAG_W - 2.0) / 2, y0 + h / 2, src, ha="center", va="center",
                fontsize=pt(8), color=MUTED, linespacing=1.35, zorder=3)

        ax.text(TXT_X, y0 + h - 3.6, head, ha="left", va="center",
                fontsize=pt(12), fontweight="bold", color=col, zorder=3)
        ax.text(TXT_X + 13.0, y0 + h - 3.6, label, ha="left", va="center",
                fontsize=pt(9.5), fontweight="bold", color=INK, zorder=3)
        ax.text(TXT_X, y0 + h - 8.4, crit[0], ha="left", va="center",
                fontsize=pt(8), color=BODY, zorder=3)
        ax.text(TXT_X, y0 + h - 11.6, crit[1], ha="left", va="center",
                fontsize=pt(8), color=BODY, zorder=3)
        if res:
            ax.text(99, y0 + h - 3.6, res, ha="right", va="center",
                    fontsize=pt(8), color=BODY, zorder=3)

    # ------------------------------------------------- Layer 3 micro-panel ------
    # Three rank correlations share one zero baseline, so their signs are read by
    # eye; the band itself carries no axis furniture.
    y3 = 30
    base = y3 + 3.0            # zero baseline, lower part of the band
    ax.plot([66, 99], [base, base], color=RULE, lw=0.7, zorder=3)
    ax.text(66, y3 + 12.1, "Spearman ρ with each cohort", ha="left", va="center",
            fontsize=pt(8), color=MUTED, zorder=3)
    PTS = []
    for accession, (name, x, colour) in COHORT_KEYS.items():
        coordinate = round(rho[name], 3)
        label = f"{rho[name]:+.2f}"
        if coordinate != EXPECTED_RHO[name] or label != EXPECTED_LABEL[name]:
            raise SystemExit(
                f"{name}: rho={rho[name]:.6f} renders as {coordinate:.3f}/"
                f"{label}, expected {EXPECTED_RHO[name]:.3f}/"
                f"{EXPECTED_LABEL[name]}"
            )
        PTS.append((x, coordinate, label, colour, name))
    for x, v, lab, c, nm in PTS:
        yy = base + v / 0.35 * 5.0          # keeps every stem between y=33 and y=39
        ax.plot([x, x], [base, yy], color=c, lw=0.8, zorder=4)
        ax.plot([x], [yy], "o", ms=4.0, color=c, zorder=5)
        ax.text(x, yy + 1.2, lab, ha="center", va="center", fontsize=pt(8),
                color=BODY, zorder=5)
        ax.text(x, base - 1.6, nm, ha="center", va="center", fontsize=pt(8),
                color=MUTED, zorder=5)

    # ------------------------------------------------------------ terminal ------
    ax.add_patch(Rectangle((0.8, 1.5), 61.5, 11.5, facecolor=TINT["sha"],
                           edgecolor="none", zorder=0))
    ax.add_patch(Rectangle((0.8, 1.5), 1.6, 11.5, facecolor=SHA,
                           edgecolor="none", zorder=2))
    ax.text(4.5, 12.1, "Finding", ha="left", va="center", fontsize=pt(9.5),
            fontweight="bold", color=INK, zorder=3)
    ax.text(4.5, 5.9, "One shared causal signal at 17q21.32,\n"
                      "in liver (SKAP1, PP.H4 = 0.93) and muscle (0.90);\n"
                      "the 57-variant credible set identifies a signal,\n"
                      "not a variant; a third signal at 15q15.1 did not survive",
            ha="left", va="center", fontsize=pt(8), color=BODY, linespacing=1.4, zorder=3)

    ax.add_patch(Rectangle((64.5, 1.5), 34.7, 11.5, facecolor=PANEL,
                           edgecolor="none", zorder=0))
    ax.add_patch(Rectangle((64.5, 1.5), 1.6, 11.5, facecolor=GREY,
                           edgecolor="none", zorder=2))
    ax.text(68.2, 12.1, "Boundary", ha="left", va="center", fontsize=pt(9.5),
            fontweight="bold", color=INK, zorder=3)
    ax.text(68.2, 5.9, "Uterine eQTL is myometrium, not\n"
                       "epithelium; clamp cohort 51 vs 18\n"
                       "limits precision",
            ha="left", va="center", fontsize=pt(8), color=BODY, linespacing=1.4, zorder=3)

    OUT = CFG.PATHS.figures
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "fig1_design.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(OUT / "fig1_design.png", bbox_inches="tight", dpi=600)
    print("wrote", OUT / "fig1_design.pdf")


if __name__ == "__main__":
    main()
