#!/usr/bin/env Rscript
# EC_IR v3 图 1–3（按 P0 修复后的新结论重画）
# 数据：
#   Fig1  pipeline/results/coloc_ramach_all.tsv   (503 次检验 / 19 位点 / 5 组织)
#   Fig2  w3/results/incremental_{result,null}_v2.tsv  + incremental_{result,null}.tsv
#   Fig3  w4/results/concordance_v2.tsv
# 规范：Arial >= 8pt、中国风色系按语义、近正方形、文字不截断

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(scales)
})

ROOT <- Sys.getenv("ECIR_ROOT", unset = normalizePath("."))
OUT  <- file.path(ROOT, "results/figures")

CINNABAR <- "#C23531"   # 朱砂：目标器官（子宫）/ 金标准表型
STONE    <- "#3D6BA8"   # 石青：代谢组织（肝）
INDIGO   <- "#177CB0"   # 靛青：骨骼肌
STONELT  <- "#7BA7C7"   # 石青浅：内脏脂肪
ROUGE    <- "#9D2933"   # 胭脂：胰腺
GREY     <- "#6E6E6E"

# 文件中的实际写法（勿改成 GTEx 风格，否则映射为 NA 会静默丢数据）
TIS_COL <- c("Uterus" = CINNABAR, "Liver" = STONE, "Muscle_Skeletal" = INDIGO,
             "Adipose_Visceral_Omentum" = STONELT, "Pancreas" = ROUGE)
TIS_LAB <- c("Uterus" = "Uterus", "Liver" = "Liver", "Muscle_Skeletal" = "Muscle",
             "Adipose_Visceral_Omentum" = "Adipose", "Pancreas" = "Pancreas")

theme_ec <- function(base = 8.5) {
  theme_minimal(base_size = base, base_family = "Arial") +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      axis.text          = element_text(colour = "grey20", size = base),
      axis.title         = element_text(colour = "grey10", size = base + 0.5),
      plot.title         = element_text(face = "bold", size = base + 1, hjust = 0),
      plot.subtitle      = element_text(colour = "grey35", size = base - 0.3, hjust = 0),
      legend.title       = element_blank(),
      legend.text        = element_text(size = base - 0.5),
      legend.key.size    = unit(0.9, "lines"),
      legend.position    = "right",
      plot.margin        = margin(5, 8, 5, 5)
    )
}

## ---------------- 数据 ----------------
d <- read.delim(file.path(ROOT, "results/coloc_ramach_all.tsv"))
parts <- do.call(rbind, strsplit(d$locus, ":"))
d$chr <- as.integer(parts[, 1]); d$pos <- as.integer(parts[, 2])
d$loc_lab <- sprintf("%s:%.0f", d$chr, d$pos / 1e6)
lv <- unique(d[order(d$chr, d$pos), c("loc_lab", "chr", "pos")])
d$loc_lab <- factor(d$loc_lab, levels = lv$loc_lab)

SYM <- c("ENSG00000141293" = "SKAP1", "ENSG00000278765" = "RP5-890E16.5",
         "ENSG00000275636" = "RP11-521C20.5", "ENSG00000108468" = "CBX1")
d$sym <- ifelse(d$gene_id %in% names(SYM), SYM[d$gene_id], d$gene_id)
stopifnot(all(d$tissue %in% names(TIS_COL)), all(d$tissue %in% names(TIS_LAB)))
TIS_LV <- c("Uterus", "Liver", "Muscle_Skeletal", "Adipose_Visceral_Omentum", "Pancreas")
d$tis <- factor(d$tissue, levels = TIS_LV)

## ---------------- 未检验的位点-组织格子 ----------------
grid <- expand.grid(loc_lab = levels(d$loc_lab), tis = TIS_LV, stringsAsFactors = FALSE)
tested <- unique(d[, c("loc_lab", "tis")]); tested$loc_lab <- as.character(tested$loc_lab)
untested <- grid[!paste(grid$loc_lab, grid$tis) %in% paste(tested$loc_lab, tested$tis), ]
untested$loc_lab <- factor(untested$loc_lab, levels = levels(d$loc_lab))
n_cells <- nrow(grid); n_tested <- nrow(tested)
cat(sprintf("位点-组织格子：%d/%d 有检验，%d 格未检验\n", n_tested, n_cells, nrow(untested)))

## ---------------- Fig 1A ----------------
set.seed(49)
pA <- ggplot(d, aes(loc_lab, PP.H4)) +
  geom_hline(yintercept = 0.8, linetype = "22", colour = GREY, linewidth = 0.4) +
  geom_hline(yintercept = 0.5, linetype = "11", colour = "grey80", linewidth = 0.3) +
  geom_point(aes(colour = tissue, size = n_snps),
             position = position_jitter(width = 0.22, height = 0, seed = 49),
             alpha = 0.85) +
  geom_point(data = untested, aes(loc_lab, y = 0), shape = 4, size = 1.5,
             colour = "grey60", inherit.aes = FALSE) +
  annotate("text", x = 0.65, y = 0.83, label = "0.8 threshold", size = 2.7,
           colour = GREY, hjust = 0, vjust = 0) +
  scale_colour_manual(values = TIS_COL, labels = TIS_LAB, name = NULL) +
  scale_size_continuous(range = c(0.7, 3.0), guide = "none") +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.2)) +
  labs(title = sprintf("A   Shared causal variation across %d risk loci", length(levels(d$loc_lab))),
       subtitle = sprintf("%d tests; point size = matched variants; \u00d7 = not tested (%d of %d locus-tissue cells)",
                          nrow(d), n_cells - n_tested, n_cells),
       x = "Locus (chromosome:position, Mb)", y = "PP.H4 (shared causal variant)") +
  theme_ec() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.8))

## ---------------- Fig 1B ----------------
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
  labs(title = "B   Tests above PP.H4 0.5",
       x = "PP.H4", y = NULL) +
  theme_ec() +
  theme(axis.text.y = element_text(size = 7.4))

## ---------------- Fig 1C ----------------
agg <- aggregate(PP.H3 ~ locus + loc_lab + tis, data = d, FUN = max)
pC <- ggplot(agg, aes(loc_lab, PP.H3)) +
  geom_hline(yintercept = 0.8, linetype = "22", colour = GREY, linewidth = 0.4) +
  geom_point(aes(colour = tis), size = 2.0, alpha = 0.9,
             position = position_jitter(width = 0.18, height = 0, seed = 49)) +
  scale_colour_manual(values = TIS_COL, labels = TIS_LAB, name = NULL) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.2)) +
  labs(title = "C   Distinct causal variants",
       subtitle = sprintf("%d locus-by-tissue combinations; %d above 0.8",
                          nrow(agg), sum(agg$PP.H3 > 0.8)),
       x = "Locus (chromosome:position, Mb)",
       y = "PP.H3 (two distinct causal variants)") +
  theme_ec() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.8))

leg <- get_legend(pA + theme(legend.position = "bottom",
                             legend.direction = "horizontal"))

fig1 <- plot_grid(
  pA + theme(legend.position = "none"),
  plot_grid(pB, pC + theme(legend.position = "none"), nrow = 1, rel_widths = c(0.95, 1.05)),
  leg,
  ncol = 1, rel_heights = c(1.00, 0.92, 0.10))

ggsave(file.path(OUT, "fig2_coloc.pdf"), fig1, width = 7.2, height = 4.9, device = cairo_pdf)
ggsave(file.path(OUT, "fig2_coloc.png"), fig1, width = 7.2, height = 4.9, dpi = 600, type = "cairo")

## ---------------- Fig 2 ----------------
n2 <- read.delim(file.path(ROOT, "results/incremental_null_v2.tsv"))
n1 <- read.delim(file.path(ROOT, "results/incremental_null.tsv"))
r2 <- read.delim(file.path(ROOT, "results/incremental_result_v2.tsv"))
r1 <- read.delim(file.path(ROOT, "results/incremental_result.tsv"))

nd <- rbind(
  data.frame(delta = n2$delta_auc, defn = "Obesity-associated programme"),
  data.frame(delta = n1$delta_auc, defn = "Earlier definition"))
obs <- data.frame(
  delta = c(r2$delta_auc_IR, r1$delta_auc_IR),
  p     = c(r2$p_IR, r1$p_IR),
  defn  = c("Obesity-associated programme", "Earlier definition"))
obs$defn <- factor(obs$defn, levels = c("Obesity-associated programme", "Earlier definition"))
nd$defn  <- factor(nd$defn,  levels = levels(obs$defn))
obs$lab <- sprintf("observed \u0394AUC = %+.4f\np = %.2f", obs$delta, obs$p)

p2 <- ggplot(nd, aes(delta)) +
  geom_histogram(bins = 22, fill = "grey80", colour = "white", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_vline(data = obs, aes(xintercept = delta), colour = CINNABAR, linewidth = 0.7) +
  geom_text(data = obs, aes(x = delta, y = Inf, label = lab),
            hjust = -0.06, vjust = 1.5, size = 2.6, colour = CINNABAR, lineheight = 0.95) +
  facet_wrap(~ defn, nrow = 1, scales = "free_y") +
  scale_x_continuous(labels = function(x) sprintf("%.3f", x)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.30))) +
  labs(title = "Incremental value of the programme for predicting metastasis",
       subtitle = "145 endometrioid carcinomas, 45 events; grey = 100 expression-matched random modules",
       x = "\u0394AUC versus base model", y = "Random modules") +
  theme_ec() +
  theme(legend.position = "none",
        strip.text = element_text(size = 8.2, face = "bold", hjust = 0),
        strip.background = element_blank(),
        panel.spacing = unit(1.2, "lines"))

ggsave(file.path(OUT, "fig3_incremental.pdf"), p2, width = 5.6, height = 2.9, device = cairo_pdf)
ggsave(file.path(OUT, "fig3_incremental.png"), p2, width = 5.6, height = 2.9, dpi = 600, type = "cairo")

## ---------------- Fig 3 ----------------
cc <- read.delim(file.path(ROOT, "results/concordance_v2.tsv"))
lvl <- c("GSE13070 muscle, clamp-defined IR vs IS",
         "GSE18732 muscle, type 2 diabetes vs NGT",
         "GSE15773 adipose, IR vs IS")
nice <- c("Clamp-defined IR\n(muscle)", "Type 2 diabetes\n(muscle)", "IR vs IS\n(adipose)")
cc$cohort <- factor(cc$cohort, levels = lvl)
cc$grp <- factor(nice[match(cc$cohort, lvl)], levels = nice)

long <- rbind(
  data.frame(grp = cc$grp, set = "All programme genes",
             val = cc$conc_all, lo = cc$ci_all_lo, hi = cc$ci_all_hi),
  data.frame(grp = cc$grp, set = "Responsive genes (p<0.05)",
             val = cc$conc_resp, lo = cc$ci_resp_lo, hi = cc$ci_resp_hi))
long$set <- factor(long$set, levels = c("All programme genes", "Responsive genes (p<0.05)"))

rho_lab <- data.frame(grp = cc$grp, lab = sprintf("\u03c1 = %+.3f", cc$rho))

p3 <- ggplot(long, aes(grp, val, colour = grp, shape = set, group = set)) +
  geom_hline(yintercept = 50, linetype = "22", colour = GREY, linewidth = 0.4) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, linewidth = 0.4,
                position = position_dodge(width = 0.42)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.42)) +
  geom_text(data = rho_lab, aes(x = grp, y = 41, label = lab), inherit.aes = FALSE,
            size = 2.7, colour = "grey25") +
  annotate("text", x = 3.45, y = 51.6, label = "chance (50%)", size = 2.5,
           colour = GREY, hjust = 1) +
  scale_colour_manual(values = c(CINNABAR, STONE, INDIGO), guide = "none") +
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

ggsave(file.path(OUT, "fig4_concordance.pdf"), p3, width = 5.5, height = 4.2, device = cairo_pdf)
ggsave(file.path(OUT, "fig4_concordance.png"), p3, width = 5.5, height = 4.2, dpi = 600, type = "cairo")

cat("图 1–3 已生成：\n")
cat("  Fig1", sum(d$PP.H4 > 0.8), "tests >0.8 /", nrow(d), "total\n")
cat("  Fig1B", nrow(top), "tests >0.5\n")
cat("  Fig1C", nrow(agg), "combos,", sum(agg$PP.H3 > 0.8), "with H3>0.8\n")
cat("  Fig2 observed", sprintf("%+.4f / %+.4f", r2$delta_auc_IR, r1$delta_auc_IR), "\n")
cat("  Fig3 cohorts", nrow(cc), "\n")
