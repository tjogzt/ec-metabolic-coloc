#!/usr/bin/env Rscript
# 01_main_figures_1_2_3.R
# =======================
#
# Purpose
# -------
# Draw the main manuscript figures from the shipped result tables:
#   fig2_coloc.pdf          posterior landscape at 19 risk loci (panels A-C)
#   fig3_incremental.pdf    incremental value of the programme for metastasis
#   fig4_concordance.pdf    sign concordance with three independent cohorts
# A cairo PDF plus a 600 dpi PNG preview is written for each figure.
#
# Author  : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
# Created : 2026-11-11
# Version : 2.0.0 (replication package restructure)
#
# Inputs
# ------
# results/tables/coloc_all_tests.tsv         503 gene-locus tests, 19 loci, 5 tissues
# results/tables/incremental_value.tsv       observed incremental AUC, reported definition
# results/tables/incremental_null.tsv        100 expression-matched random modules
# results/tables/concordance_by_cohort.tsv   sign concordance, three cohorts
# results/script_inputs/incremental_result.tsv and incremental_null.tsv
#                                            superseded programme definition (see Notes)
# config/config.R                            PALETTE, FIGURE_FONT, FIGURE_MIN_PT, paths
#
# Outputs
# -------
# results/figures/fig2_coloc.{pdf,png}
# results/figures/fig3_incremental.{pdf,png}
# results/figures/fig4_concordance.{pdf,png}
#
# Notes
# -----
# * Figure style is frozen: Arial, no type below FIGURE_MIN_PT, near-square
#   canvases, panel letters in uppercase type: figure 2 carries A/B/C in the
#   panel titles, figure 4 carries A/B in the facet strips, consistent with
#   the journal requirement for composite figures.
# * Colours are semantic, never ordinal: they are read from PALETTE in
#   config/config.R and re-labelled to tissues. Nothing here redefines a hex
#   value, so a change of palette propagates to every figure at once.
# * Point size in panel A encodes the number of matched variants, so a dense
#   locus is not read as a stronger signal.
# * The superseded 608-gene definition of the programme is the second facet of
#   fig3_incremental. It is not part of the v2 table set: its two inputs still
#   live in results/script_inputs, the directory the pre-restructure script
#   read. They are located relative to TABLES (never as an absolute path) and
#   the facet is dropped with a warning if they are absent, so the reported
#   definition always reproduces.
# * Nothing in this script writes outside results/figures.

# ------------------------------------------------------------ bootstrap ----
# Resolve the repository root from the script location so the file can be run
# from any working directory, then let config.R own every path below it.
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(script_path) == 0) script_path <- "code/07_figures/01_main_figures_1_2_3.R"
root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE)
Sys.setenv(ECIR_ROOT = root)
source(file.path(root, "config", "config.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(scales)
})

# ggplot2 expresses geom_text() sizes in millimetres; one point is 72.27/25.4 mm.
.pt <- 72.27 / 25.4

# -------------------------------------------------------------- palette ----
# Semantic keys, not series positions. TIS_COL keeps the tissue spellings used
# in the tables (GTEx-style labels would map to NA and silently drop data).
GREY <- PALETTE[["neutral"]]

TIS_COL <- c(
  "Uterus"                   = PALETTE[["uterus"]],
  "Liver"                    = PALETTE[["liver"]],
  "Muscle_Skeletal"          = PALETTE[["muscle"]],
  "Adipose_Visceral_Omentum" = PALETTE[["adipose"]],
  "Pancreas"                 = PALETTE[["pancreas"]]
)
TIS_LAB <- c(
  "Uterus"                   = "Uterus",
  "Liver"                    = "Liver",
  "Muscle_Skeletal"          = "Muscle",
  "Adipose_Visceral_Omentum" = "Adipose",
  "Pancreas"                 = "Pancreas"
)
TIS_LV <- names(TIS_COL)

theme_ec <- function(base = 8.5) {
  stopifnot(base >= FIGURE_MIN_PT)
  theme_minimal(base_size = base, base_family = FIGURE_FONT) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      axis.text          = element_text(colour = "grey20", size = base),
      axis.title         = element_text(colour = "grey10", size = base + 0.5),
      plot.title         = element_text(face = "bold", size = base + 1, hjust = 0),
      plot.subtitle      = element_text(colour = "grey35", size = base - 0.3, hjust = 0),
      legend.title       = element_blank(),
      legend.text        = element_text(size = max(base - 0.5, FIGURE_MIN_PT)),
      legend.key.size    = unit(0.9, "lines"),
      legend.position    = "right",
      plot.margin        = margin(5, 8, 5, 5)
    )
}

# Emit the vector original and the raster preview side by side.
save_pair <- function(plot, stem, width, height) {
  ggsave(file.path(FIGURES, paste0(stem, ".pdf")), plot,
         width = width, height = height, device = cairo_pdf)
  ggsave(file.path(FIGURES, paste0(stem, ".png")), plot,
         width = width, height = height, dpi = 600, type = "cairo")
  cat(sprintf("wrote %s.{pdf,png}  (%.1f x %.1f in)\n", stem, width, height))
}

# ===========================================================================
# Figure 2 --- colocalisation across 19 risk loci and five tissues
# ===========================================================================
d <- read.delim(file.path(TABLES, "coloc_all_tests.tsv"))
stopifnot(nrow(d) > 0, all(c("PP.H3", "PP.H4", "locus", "tissue") %in% names(d)))

# Locus labels are chromosome:megabase, ordered by genomic position.
parts <- do.call(rbind, strsplit(d$locus, ":"))
d$chr <- as.integer(parts[, 1])
d$pos <- as.integer(parts[, 2])
d$loc_lab <- sprintf("%s:%.0f", d$chr, d$pos / 1e6)
lv <- unique(d[order(d$chr, d$pos), c("loc_lab", "chr", "pos")])
d$loc_lab <- factor(d$loc_lab, levels = lv$loc_lab)

# Stable gene symbols for the four genes that carry a result.
SYM <- c("ENSG00000141293" = "SKAP1", "ENSG00000278765" = "RP5-890E16.5",
         "ENSG00000275636" = "RP11-521C20.5", "ENSG00000108468" = "CBX1")
d$sym <- ifelse(d$gene_id %in% names(SYM), SYM[d$gene_id], d$gene_id)

stopifnot(all(d$tissue %in% names(TIS_COL)), all(d$tissue %in% names(TIS_LAB)))
d$tis <- factor(d$tissue, levels = TIS_LV)

# ---------------------------------------------------------------------------
# Locus-tissue cells that were never tested. Drawing them as crosses is the
# difference between "no evidence" and "no test", which the panel must show.
# ---------------------------------------------------------------------------
grid <- expand.grid(loc_lab = levels(d$loc_lab), tis = TIS_LV, stringsAsFactors = FALSE)
tested <- unique(d[, c("loc_lab", "tis")])
tested$loc_lab <- as.character(tested$loc_lab)
untested <- grid[!paste(grid$loc_lab, grid$tis) %in% paste(tested$loc_lab, tested$tis), ]
untested$loc_lab <- factor(untested$loc_lab, levels = levels(d$loc_lab))
n_cells <- nrow(grid)
n_tested <- nrow(tested)
cat(sprintf("locus-tissue cells: %d/%d tested, %d not tested\n",
            n_tested, n_cells, nrow(untested)))

# ---------------- panel A: every test, one dot per gene-locus pair ----------
set.seed(SEED)
pA <- ggplot(d, aes(loc_lab, PP.H4)) +
  geom_hline(yintercept = 0.8, linetype = "22", colour = GREY, linewidth = 0.4) +
  geom_hline(yintercept = 0.5, linetype = "11", colour = "grey80", linewidth = 0.3) +
  geom_point(aes(colour = tissue, size = n_snps),
             position = position_jitter(width = 0.22, height = 0, seed = SEED),
             alpha = 0.85) +
  geom_point(data = untested, aes(loc_lab, y = 0), shape = 4, size = 1.5,
             colour = "grey60", inherit.aes = FALSE) +
  annotate("text", x = 0.65, y = 0.83, label = "0.8 threshold", size = 3.4,
           colour = GREY, hjust = 0, vjust = 0) +
  scale_colour_manual(values = TIS_COL, labels = TIS_LAB, name = NULL) +
  scale_size_continuous(range = c(0.7, 3.0), guide = "none") +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.2)) +
  labs(title = sprintf("A   Shared causal variation across %d risk loci",
                       length(levels(d$loc_lab))),
       subtitle = sprintf("%d tests; point size = matched variants\n\u00d7 = not tested (%d of %d locus-tissue cells)",
                          nrow(d), n_cells - n_tested, n_cells),
       x = "Locus (chromosome:position, Mb)",
       y = "PP.H4 (shared causal variant)") +
  theme_ec(9.6) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9.6),
        plot.subtitle = element_text(lineheight = 1.35))

# ---------------- panel B: the tests worth naming --------------------------
top <- subset(d, PP.H4 > 0.5)
top <- top[order(top$PP.H4), ]
top$lab <- factor(sprintf("%s \u00b7 %s", top$sym, TIS_LAB[top$tissue]),
                  levels = rev(sprintf("%s \u00b7 %s", top$sym, TIS_LAB[top$tissue])))
top$hl <- top$PP.H4 > 0.8

pB <- ggplot(top, aes(PP.H4, lab)) +
  geom_vline(xintercept = 0.8, linetype = "22", colour = GREY, linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = PP.H4, yend = lab, colour = tissue), linewidth = 0.7) +
  geom_point(aes(colour = tissue, shape = hl), size = 2.4) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  scale_colour_manual(values = TIS_COL, guide = "none") +
  scale_x_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.2)) +
  labs(title = "B   Tests above PP.H4 0.5", x = "PP.H4", y = NULL) +
  theme_ec(9.6) +
  theme(axis.text.y = element_text(size = 9.6))

# ---------------- panel C: distinct-causal-variant posterior ---------------
agg <- aggregate(PP.H3 ~ locus + loc_lab + tis, data = d, FUN = max)
pC <- ggplot(agg, aes(loc_lab, PP.H3)) +
  geom_hline(yintercept = 0.8, linetype = "22", colour = GREY, linewidth = 0.4) +
  geom_point(aes(colour = tis), size = 2.0, alpha = 0.9,
             position = position_jitter(width = 0.18, height = 0, seed = SEED)) +
  scale_colour_manual(values = TIS_COL, labels = TIS_LAB, name = NULL) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.2)) +
  labs(title = "C   Distinct causal variants",
       subtitle = sprintf("%d locus-by-tissue combinations; %d above 0.8",
                          nrow(agg), sum(agg$PP.H3 > 0.8)),
       x = "Locus (chromosome:position, Mb)",
       y = "PP.H3 (two distinct causal variants)") +
  theme_ec(9.6) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9.6),
        plot.subtitle = element_text(lineheight = 1.35))

# One shared tissue legend under the three panels.
leg <- get_legend(pA + theme(legend.position = "bottom",
                             legend.direction = "horizontal"))

fig1 <- plot_grid(
  pA + theme(legend.position = "none"),
  plot_grid(pB, pC + theme(legend.position = "none"), nrow = 1,
            rel_widths = c(0.95, 1.05)),
  leg,
  ncol = 1, rel_heights = c(1.00, 0.92, 0.10))

save_pair(fig1, "fig2_coloc", width = 7.2, height = 5.7)

# ===========================================================================
# Figure 3 --- incremental value of the programme over the base model
# ===========================================================================
# The reported definition ships as a table. The superseded 608-gene definition
# survives only in results/script_inputs; because both definitions use a file
# called incremental_null.tsv, each series lists its own explicit pair.
legacy_dir <- file.path(dirname(TABLES), "script_inputs")

series <- list(
  list(label    = "Obesity-associated programme",
       observed = file.path(TABLES, "incremental_value.tsv"),
       null     = file.path(TABLES, "incremental_null.tsv")),
  list(label    = "Earlier definition",
       observed = file.path(legacy_dir, "incremental_result.tsv"),
       null     = file.path(legacy_dir, "incremental_null.tsv"))
)
series <- Filter(function(s) file.exists(s$observed) && file.exists(s$null), series)
if (length(series) < 2L) {
  warning("superseded programme definition not found under ",
          legacy_dir, "; fig3_incremental shows the reported definition only",
          call. = FALSE)
}

nd <- do.call(rbind, lapply(series, function(s) {
  x <- read.delim(s$null)
  data.frame(delta = x$delta_auc, defn = s$label)
}))
obs <- do.call(rbind, lapply(series, function(s) {
  x <- read.delim(s$observed)
  data.frame(delta = x$delta_auc_IR, p = x$p_IR, n = x$n,
             events = x$events, defn = s$label)
}))

obs$defn <- factor(obs$defn, levels = vapply(series, `[[`, character(1), "label"))
nd$defn  <- factor(nd$defn,  levels = levels(obs$defn))
obs$lab  <- sprintf("observed \u0394AUC = %+.4f\np = %.2f", obs$delta, obs$p)

# The subtitle states the cohort and the null size once; both facets are built
# from the same 145 tumours and the same number of random modules.
n_null <- sum(nd$defn == levels(nd$defn)[1])
n_case <- obs$n[1]
n_event <- obs$events[1]

# Journal requirement: denote individual panels of a composite figure with
# uppercase letters. These two facets are two panels, so each strip carries
# its letter rather than relying on a separate legend box.
defn_lab <- stats::setNames(
  paste0(LETTERS[seq_along(levels(obs$defn))], "   ", levels(obs$defn)),
  levels(obs$defn))

p2 <- ggplot(nd, aes(delta)) +
  geom_histogram(bins = 22, fill = "grey80", colour = "white", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_vline(data = obs, aes(xintercept = delta), colour = PALETTE[["uterus"]],
             linewidth = 0.7) +
  geom_text(data = obs, aes(x = delta, y = Inf, label = lab),
            hjust = -0.06, vjust = 1.5, size = FIGURE_MIN_PT / .pt,
            colour = PALETTE[["uterus"]], lineheight = 0.95) +
  facet_wrap(~ defn, nrow = 1, scales = "free_y",
             labeller = labeller(defn = defn_lab)) +
  scale_x_continuous(labels = function(x) sprintf("%.3f", x)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.30))) +
  labs(title = "Incremental value of the programme for predicting metastasis",
       subtitle = sprintf("%d endometrioid carcinomas, %d events; grey = %d expression-matched random modules",
                          n_case, n_event, n_null),
       x = "\u0394AUC versus base model", y = "Random modules") +
  theme_ec() +
  theme(legend.position = "none",
        strip.text = element_text(size = 8.2, face = "bold", hjust = 0),
        strip.background = element_blank(),
        panel.spacing = unit(1.2, "lines"))

save_pair(p2, "fig3_incremental", width = 5.6, height = 2.9)

# ===========================================================================
# Figure 4 --- concordance of programme scores with three cohorts
# ===========================================================================
cc <- read.delim(file.path(TABLES, "concordance_by_cohort.tsv"))

lvl  <- c("GSE13070 muscle, clamp-defined IR vs IS",
          "GSE18732 muscle, type 2 diabetes vs NGT",
          "GSE15773 adipose, IR vs IS")
nice <- c("Clamp-defined IR\n(muscle)", "Type 2 diabetes\n(muscle)", "IR vs IS\n(adipose)")

stopifnot(all(lvl %in% cc$cohort))
cc$cohort <- factor(cc$cohort, levels = lvl)
cc$grp <- factor(nice[match(cc$cohort, lvl)], levels = nice)

long <- rbind(
  data.frame(grp = cc$grp, set = "All programme genes",
             val = cc$conc_all, lo = cc$ci_all_lo, hi = cc$ci_all_hi),
  data.frame(grp = cc$grp, set = "Responsive genes (p<0.05)",
             val = cc$conc_resp, lo = cc$ci_resp_lo, hi = cc$ci_resp_hi))
long$set <- factor(long$set, levels = c("All programme genes",
                                        "Responsive genes (p<0.05)"))

rho_lab <- data.frame(grp = cc$grp, lab = sprintf("\u03c1 = %+.3f", cc$rho))

# Series colours follow the phenotype hierarchy fixed in the shipped figure:
# the clamp cohort carries the gold-standard cinnabar, the diagnosed cohort the
# azurite of the metabolic tissue, the adipose contrast the indigo of muscle.
series_col <- c(PALETTE[["uterus"]], PALETTE[["liver"]], PALETTE[["muscle"]])
names(series_col) <- nice

p3 <- ggplot(long, aes(grp, val, colour = grp, shape = set, group = set)) +
  geom_hline(yintercept = 50, linetype = "22", colour = GREY, linewidth = 0.4) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, linewidth = 0.4,
                position = position_dodge(width = 0.42)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.42)) +
  geom_text(data = rho_lab, aes(x = grp, y = 41, label = lab), inherit.aes = FALSE,
            size = FIGURE_MIN_PT / .pt, colour = "grey25") +
  annotate("text", x = 3.45, y = 51.6, label = "chance (50%)",
           size = FIGURE_MIN_PT / .pt, colour = GREY, hjust = 1) +
  scale_colour_manual(values = series_col, guide = "none") +
  scale_shape_manual(values = c(16, 1), name = NULL) +
  scale_y_continuous(limits = c(38, 90), breaks = seq(40, 90, 10),
                     labels = function(x) paste0(x, "%")) +
  labs(title = "Concordance of programme scores with three independent cohorts",
       subtitle = "Sign agreement (95% CI); filled = all genes, open = responsive subset",
       x = NULL, y = "Direction agreement") +
  theme_ec() +
  theme(legend.position = "top",
        legend.margin = margin(0, 0, -4, 0),
        axis.text.x = element_text(size = 8.0))

save_pair(p3, "fig4_concordance", width = 5.5, height = 4.2)

# ---------------------------------------------------------------- report ----
cat("main figures written to", FIGURES, ":\n")
cat("  fig2 panel A:", sum(d$PP.H4 > 0.8), "tests >0.8 /", nrow(d), "total\n")
cat("  fig2 panel B:", nrow(top), "tests >0.5\n")
cat("  fig2 panel C:", nrow(agg), "combos,", sum(agg$PP.H3 > 0.8), "with H3>0.8\n")
cat("  fig3 observed:", sprintf("%+.4f / %+.4f",
                                obs$delta[1], obs$delta[nrow(obs)]), "\n")
cat("  fig4 cohorts:", nrow(cc), "\n")
