#!/usr/bin/env python3

import os as _os
import sys as _sys
import pathlib as _pathlib

_root = _pathlib.Path(__file__).resolve()
while not (_root / "config.py").exists():
    if _root.parent == _root:
        raise RuntimeError("config.py not found; set ECIR_ROOT")
    _root = _root.parent
if str(_root) not in _sys.path:
    _sys.path.insert(0, str(_root))
import config as CFG
"""Fig 1 v3：证据横带版式，改为通栏文字布局以消除列宽冲突。

每带三行：
  第一行：头条数字（12pt 粗体，语义色）+ 层名（9.5pt 粗体）+ 右侧结果徽章
  第二、三行：判据（8pt），通栏左对齐
Layer 3 右侧留出微图槽（三队列秩相关，共用零基线）。
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})

SHA, SHI, DIAN = "#C23531", "#3D6BA8", "#177CB0"
GREY, INK = "#8A8A8A", "#1A1A1A"
TINT = {"sha": "#FBEDEB", "dian": "#EAF2F8"}

fig, ax = plt.subplots(figsize=(7.2, 5.8))
ax.set_xlim(0, 100); ax.set_ylim(0, 80); ax.axis("off")

TAG_W = 15.0
SPINE_X, SPINE_W = 16.6, 1.6
TXT_X = 19.6

BANDS = [
    (60, 13, "862", "Layer 1  Programme definition", SHA, TINT["sha"],
     ["FDR<0.05, |logFC|$\\geq$0.3, endometrial expression gate; 608 up / 254 down",
      "preregistered, hashed, and carried unchanged into Layers 3 and 4"],
     "genes in programme", "Discovery\ncontrast"),
    (45, 13, "2 / 259", "Layer 2  Colocalisation", SHA, TINT["sha"],
     ["Gene-locus tests across 10 loci and five tissues; PP.H4$>$0.8 threshold",
      "35 of 50 locus-tissue cells testable, adipose contributing one locus"],
     "both at 17q21.32", "Risk GWAS\n+ cis-eQTL"),
    (30, 13, "$+$0.31", "Layer 3  Bulk proxy validity", DIAN, TINT["dian"],
     ["Rank correlation with three independent cohorts",
      "adipose strong, diabetes intermediate, clamp near zero"],
     None, "(from Layer 1)"),
    (15, 13, "$-$3.8", "Layer 4  Single-cell projection", DIAN, TINT["dian"],
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
            fontsize=8, color="#555555", linespacing=1.35, zorder=3)

    ax.text(TXT_X, y0 + h - 3.6, head, ha="left", va="center",
            fontsize=12, fontweight="bold", color=col, zorder=3)
    ax.text(TXT_X + 13.0, y0 + h - 3.6, label, ha="left", va="center",
            fontsize=9.5, fontweight="bold", color=INK, zorder=3)
    ax.text(TXT_X, y0 + h - 8.4, crit[0], ha="left", va="center",
            fontsize=8, color="#333333", zorder=3)
    ax.text(TXT_X, y0 + h - 11.6, crit[1], ha="left", va="center",
            fontsize=8, color="#333333", zorder=3)
    if res:
        ax.text(99, y0 + h - 3.6, res, ha="right", va="center",
                fontsize=8, color="#333333", zorder=3)

# ---------- Layer 3 微图（右侧槽位，x 66..99） ----------
y3 = 30
base = y3 + 3.0            # 零基线，位于带内下部
ax.plot([66, 99], [base, base], color="#C8C8C8", lw=0.7, zorder=3)
ax.text(66, y3 + 12.1, "Spearman $\\rho$ with each cohort", ha="left", va="center",
        fontsize=8, color="#555555", zorder=3)
PTS = [(74.0, 0.313, "+0.31", SHA, "adipose"),
       (83.0, 0.202, "+0.20", SHA, "diabetes"),
       (92.0, 0.098, "+0.10", GREY, "clamp")]
for x, v, lab, c, nm in PTS:
    yy = base + v / 0.35 * 5.0          # 全部落在 y=33~39 之间
    ax.plot([x, x], [base, yy], color=c, lw=0.8, zorder=4)
    ax.plot([x], [yy], "o", ms=4.0, color=c, zorder=5)
    ax.text(x, yy + 1.2, lab, ha="center", va="center", fontsize=8,
            color="#333333", zorder=5)
    ax.text(x, base - 1.6, nm, ha="center", va="center", fontsize=8,
            color="#555555", zorder=5)

# ---------- 底部终端 ----------
ax.add_patch(Rectangle((0.8, 1.5), 61.5, 11.5, facecolor=TINT["sha"], edgecolor="none", zorder=0))
ax.add_patch(Rectangle((0.8, 1.5), 1.6, 11.5, facecolor=SHA, edgecolor="none", zorder=2))
ax.text(4.5, 11.5, "Finding", ha="left", va="center", fontsize=9.5,
        fontweight="bold", color=INK, zorder=3)
ax.text(4.5, 5.6, "One shared causal variant at 17q21.32, colocalised in liver (SKAP1,\n"
                  "PP.H4 = 0.92) and muscle (PP.H4 = 0.87); the other nine loci are\n"
                  "compatible with distinct causal variants",
        ha="left", va="center", fontsize=8, color="#333333", linespacing=1.55, zorder=3)

ax.add_patch(Rectangle((64.5, 1.5), 34.7, 11.5, facecolor="#F7F7F7", edgecolor="none", zorder=0))
ax.add_patch(Rectangle((64.5, 1.5), 1.6, 11.5, facecolor=GREY, edgecolor="none", zorder=2))
ax.text(68.2, 11.5, "Boundary", ha="left", va="center", fontsize=9.5,
        fontweight="bold", color=INK, zorder=3)
ax.text(68.2, 5.6, "Uterine eQTL is myometrium, not epithelium;\n"
                   "clamp cohort 51 vs 18 limits precision",
        ha="left", va="center", fontsize=8, color="#333333", linespacing=1.55, zorder=3)

OUT = CFG.FIGURES
fig.savefig(OUT / "fig1_design.pdf", bbox_inches="tight", dpi=300)
fig.savefig(OUT / "fig1_design.png", bbox_inches="tight", dpi=600)
print("已写", OUT / "fig1_design.pdf")
