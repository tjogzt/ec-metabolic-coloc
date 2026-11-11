# 02_differential_expression.R
# ============================
#
# Purpose
# -------
# Fit every differential-expression model the bulk layer needs, in the order in
# which the programme and the disease signature were derived:
#
#   A. GSE24215  skeletal muscle, insulin vs baseline, paired within subject
#                -> the insulin-response arm of the transcriptional programme
#   B. GSE24422  cultured adipocytes, insulin vs control
#                -> the second arm, used only to require a shared direction
#   C. direction agreement between A and B (printed, plus an RDS cache)
#   D. GSE17025  endometrioid carcinoma vs normal endometrium (discovery)
#                -> EC signature v1
#   E. GSE196033 tumour vs normal (independent reproduction cohort)
#   F. GSE17025 again, with a per-sample proliferation score as a covariate
#                -> EC signature v2 (proliferation-adjusted), full and compact
#
# The insulin contrasts and the tumour contrasts use different effect-size
# cuts, so the cuts live in separate named constants below rather than being
# collapsed into one: a 0.3 log2 change is a meaningful insulin response in a
# paired muscle biopsy, and a meaningless one in a 91-vs-11 tumour comparison.
#
# Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
# Created  : 2026-11-11
# Version  : 2.0.0
#
# Inputs
# ------
# Environment variable ECIR_ROOT (optional) pointing at the repository root;
# otherwise the root is derived from this script's location.
#
#   <DATA>/GSE24215_series_matrix.txt.gz   GPL6480
#   <DATA>/GSE24422_series_matrix.txt.gz   GPL570
#   <DATA>/GSE17025_series_matrix.txt.gz   GPL570
#   <DATA>/GSE196033_series_matrix.txt.gz  (RNA-seq, already log-scaled counts)
#   <DATA>/GPL570_annot.tsv                platform annotation, see 03_annot*.R
#   <DATA>/GPL6480_annot.tsv
#
# All thresholds come from config/config.R (PARAMS); see the constants block for
# the two cuts that are specific to a single contrast.
#
# Outputs
# -------
# <TABLES>/DE_GSE24215_insulin_vs_basal.tsv
# <TABLES>/DE_GSE24422_insulin_adipocyte.tsv
# <TABLES>/DE_GSE17025_tumor_vs_normal.tsv
# <TABLES>/DE_GSE196033_tumor_vs_normal.tsv
# <TABLES>/DE_GSE17025_tumor_vs_normal_prolifadj.tsv
# <TABLES>/EC_signature_v1.tsv                     gene, direction, logFC, FDR
# <TABLES>/EC_signature_v2_prolifadj_full.tsv      same four columns
# <TABLES>/EC_signature_v2_prolifadj_compact.tsv   top 200 per direction
# <RESULTS>/cache/ir_de_cache.rds                  ttA, ttB, exA, exB
#
# Usage
# -----
#   Rscript code/03_burden/02_differential_expression.R
#
# Notes
# -----
# Probe-to-gene collapse keeps, for each symbol, the probe with the highest mean
# expression. That rule is fixed before any contrast is fitted and is applied
# identically to every cohort; choosing the probe per contrast would leak the
# outcome into the feature definition.
#
# Two variables from the original discovery script are not carried over: an
# unused per-sample histology vector parsed from characteristics_ch1 and an
# unused column index derived from it. Neither entered any computation, and the
# histology field is read once, where it is actually used. Two unused arm flags
# in the muscle section (basal / insulin indicators that were superseded by the
# condition vector) are dropped for the same reason.

# ---------------------------------------------------------------- setup ----
suppressMessages({
  library(GEOquery)
  library(limma)
})

# Network timeout for the GEO downloads; getGEO retries inside this budget.
HTTP_TIMEOUT <- 900
options(timeout = HTTP_TIMEOUT)

# Locate the repository, then hand every path decision to config/config.R.
.resolve_root <- function() {
  root <- Sys.getenv("ECIR_ROOT", unset = "")
  if (nzchar(root) && file.exists(file.path(root, "config", "config.R"))) {
    return(normalizePath(root))
  }
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  here <- if (length(hit)) {
    sub("^--file=", "", hit[1])
  } else {
    tryCatch(sys.frame(1)$ofile, error = function(e) "")
  }
  if (!nzchar(here)) stop("cannot locate the repository; set ECIR_ROOT")
  normalizePath(file.path(dirname(normalizePath(here)), "..", ".."),
                mustWork = TRUE)
}

source(file.path(.resolve_root(), "config", "config.R"))

CACHE <- file.path(RESULTS, "cache")
if (!dir.exists(CACHE)) dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------ constants ----
# Effect-size cut for the insulin contrasts. Paired muscle biopsies and
# adipocyte cultures give small but reproducible shifts; 0.3 log2 units is the
# boundary the programme was defined against.
INSULIN_LOGFC_CUT <- 0.3

# Per-direction cap of the compact signature. 200 is what the incremental-value
# model can carry without the feature count overwhelming 145 samples.
COMPACT_PER_DIRECTION <- 200

# Thirty-five canonical cell-cycle genes, used to build the per-sample
# proliferation score and to audit how much of the signature is proliferation.
CYCLING_GENES <- c(
  "MKI67", "PCNA", "CCNB1", "CCNB2", "CCNA2", "CDK1", "CDC20", "BUB1", "BUB1B",
  "TOP2A", "AURKA", "AURKB", "PLK1", "TPX2", "UBE2C", "RRM2", "TYMS", "MCM2",
  "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "CENPF", "KIF11", "KIF23", "ANLN",
  "PRC1", "ECT2", "NCAPG", "ASPM", "DLGAP5", "CENPE", "SMC4", "NUSAP1"
)

# -------------------------------------------------------------- helpers ----
read_series <- function(f) {
  # One series matrix; column names become sample titles so that the contrast
  # labels below can be read off the matrix itself.
  es <- getGEO(filename = f, getGPL = FALSE)
  ex <- exprs(es)
  pd <- pData(es)
  colnames(ex) <- pd[colnames(ex), "title"]
  list(ex = ex, pd = pd)
}

collapse_to_gene <- function(ex, annot, idcol = "ID", symcol = NULL) {
  # Map probes to symbols and keep one probe per symbol.
  #
  # The symbol column is named differently across platforms ("Gene Symbol",
  # "GENE_SYMBOL", "GENE"), so it is detected rather than assumed; an explicit
  # symcol argument always wins.
  if (is.null(symcol) || !(symcol %in% colnames(annot))) {
    candidates <- c("Gene Symbol", "GENE_SYMBOL", "GENE")
    symcol <- candidates[candidates %in% colnames(annot)][1]
  }
  a <- annot[match(rownames(ex), annot[[idcol]]), c(idcol, symcol)]

  # Drop probes with no usable symbol; "---" is the GPL placeholder.
  keep <- !is.na(a[[symcol]]) & a[[symcol]] != "" & a[[symcol]] != "---"
  ex <- ex[keep, , drop = FALSE]
  a <- a[keep, , drop = FALSE]

  # Highest-mean probe per symbol: order by symbol, then descending mean, and
  # take the first row of each symbol block.
  m <- rowMeans(ex, na.rm = TRUE)
  ord <- order(a[[symcol]], -m)
  a <- a[ord, ]
  ex <- ex[ord, , drop = FALSE]
  first <- !duplicated(a[[symcol]])
  ex <- ex[first, , drop = FALSE]
  rownames(ex) <- a[[symcol]][first]
  ex
}

write_de <- function(tt, filename) {
  # limma topTable output, unchanged in structure, written to the tables dir.
  path <- file.path(TABLES, filename)
  write.table(tt, path, sep = "\t", quote = FALSE)
  invisible(path)
}

# ==========================================================================
# A. GSE24215 -- skeletal muscle, insulin vs baseline, paired within subject
# ==========================================================================
insulin_response_muscle <- function() {
  A <- read_series(file.path(DATA, "GSE24215_series_matrix.txt.gz"))
  annot <- read.delim(file.path(DATA, "GPL6480_annot.tsv"), check.names = FALSE)
  exA <- collapse_to_gene(A$ex, annot)
  cat("GSE24215 genes:", nrow(exA), " samples:", ncol(exA), "\n")
  cat("expression range:", round(range(exA, na.rm = TRUE), 2), "\n")

  # Keep the insulin clamp (t0) samples and drop technical replicates; the
  # replicate would otherwise enter the paired design as an extra subject.
  titles <- colnames(exA)
  is_t0 <- grepl("_t0_", titles)
  is_tech <- grepl("techrep", titles)
  exA0 <- exA[, is_t0 & !is_tech, drop = FALSE]

  cond <- ifelse(grepl("basal", colnames(exA0)), "basal", "insulin")
  subj <- paste0("s", sub(".*rep([0-9]+).*", "\\1", colnames(exA0)))
  cat("paired comparison:", sum(cond == "basal"), "basal vs",
      sum(cond == "insulin"), "insulin\n")
  print(table(subj, cond))

  # Subject term absorbs between-person variation; cond is the fixed effect of
  # interest, so the model is a paired analysis with subject as a block.
  design <- model.matrix(~ subj + cond)
  fit <- eBayes(lmFit(exA0, design), trend = TRUE, robust = TRUE)
  ttA <- topTable(fit, coef = "condinsulin", number = Inf, sort.by = "P")
  write_de(ttA, "DE_GSE24215_insulin_vs_basal.tsv")

  cat("GSE24215 DE: FDR<", PARAMS$program_de_padj, " =",
      sum(ttA$adj.P.Val < PARAMS$program_de_padj),
      "; of those |logFC|>=", INSULIN_LOGFC_CUT, " =",
      sum(ttA$adj.P.Val < PARAMS$program_de_padj &
            abs(ttA$logFC) >= INSULIN_LOGFC_CUT), "\n")
  ttA
}

# ==========================================================================
# B. GSE24422 -- cultured adipocytes, insulin vs control (direction check)
# ==========================================================================
insulin_response_adipocyte <- function() {
  B <- read_series(file.path(DATA, "GSE24422_series_matrix.txt.gz"))
  annot <- read.delim(file.path(DATA, "GPL570_annot.tsv"), check.names = FALSE)
  exB <- collapse_to_gene(B$ex, annot)
  cat("\nGSE24422 genes:", nrow(exB), " samples:", ncol(exB), "\n")

  # Sample titles encode the arm: Adip_Ins* is insulin, Adip_rep* is control.
  tb <- colnames(exB)
  grp <- ifelse(grepl("^Adip_Ins", tb), "ins",
                ifelse(grepl("^Adip_rep", tb), "ctl", NA))
  sel <- !is.na(grp)
  exB2 <- exB[, sel, drop = FALSE]
  gB <- factor(grp[sel], levels = c("ctl", "ins"))
  cat("adipocyte arms: control", sum(gB == "ctl"), "insulin", sum(gB == "ins"),
      "\n")

  d2 <- model.matrix(~ gB)
  fitB <- eBayes(lmFit(exB2, d2), trend = TRUE, robust = TRUE)
  ttB <- topTable(fitB, coef = "gBins", number = Inf, sort.by = "P")
  write_de(ttB, "DE_GSE24422_insulin_adipocyte.tsv")

  cat("GSE24422: FDR<", PARAMS$program_de_padj, " =",
      sum(ttB$adj.P.Val < PARAMS$program_de_padj), "\n")
  ttB
}

# ==========================================================================
# C. Direction agreement between the two insulin datasets
# ==========================================================================
insulin_direction_agreement <- function(ttA, ttB) {
  # Genes that respond in muscle are asked to point the same way in adipocytes.
  # No significance requirement is imposed on the adipocyte side here -- the
  # requirement appears later, in the programme tiers.
  common <- intersect(rownames(ttA), rownames(ttB))
  cat("\ngenes shared by both insulin datasets:", length(common), "\n")
  ag <- ttA[common, ]
  bg <- ttB[common, ]

  sigA <- ag$adj.P.Val < PARAMS$program_de_padj & abs(ag$logFC) >= INSULIN_LOGFC_CUT
  cat("significant and |logFC|>=", INSULIN_LOGFC_CUT, " in muscle:",
      sum(sigA), "\n")

  if (sum(sigA) > 0) {
    s <- sign(ag$logFC[sigA])
    t2 <- sign(bg$logFC[sigA])
    agree <- mean(s == t2, na.rm = TRUE)
    cat("same direction in adipocytes:",
        round(100 * agree, 1), "% (binomial p =",
        signif(binom.test(sum(s == t2, na.rm = TRUE),
                          sum(!is.na(t2)), 0.5)$p.value, 3), ")\n")
    cat("same direction and FDR<", PARAMS$program_de_padj, "in adipocytes:",
        sum(s == t2 & bg$adj.P.Val[sigA] < PARAMS$program_de_padj,
            na.rm = TRUE), "\n")
  }
  invisible(NULL)
}

# ==========================================================================
# D. GSE17025 -- endometrioid carcinoma vs normal endometrium (discovery)
# ==========================================================================
ec_discovery <- function() {
  S <- read_series(file.path(DATA, "GSE17025_series_matrix.txt.gz"))
  annot <- read.delim(file.path(DATA, "GPL570_annot.tsv"), check.names = FALSE)
  ex <- collapse_to_gene(S$ex, annot)
  cat("\nGSE17025 genes", nrow(ex), " samples", ncol(ex), "\n")

  # The histology field distinguishes tumour subtypes from normal endometrium;
  # both endometrioid and papillary serous are tumour here.
  histv <- as.character(S$pd[[grep("^histology", colnames(S$pd))[1]]])
  histv <- sub("histology: ", "", histv)
  cat("histology distribution:\n")
  print(table(histv))

  grp <- ifelse(grepl("endometrioid|papillary serous", histv), "tumor",
                ifelse(grepl("endometrium", histv), "normal", NA))
  cat("groups:\n")
  print(table(grp))

  keep <- !is.na(grp)
  ex <- ex[, keep, drop = FALSE]
  grp <- factor(grp[keep], levels = c("normal", "tumor"))

  # Quantile normalisation across arrays before a between-group contrast.
  ex <- normalizeBetweenArrays(ex, method = "quantile")

  design <- model.matrix(~ grp)
  fit <- eBayes(lmFit(ex, design), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = "grptumor", number = Inf, sort.by = "P")
  write_de(tt, "DE_GSE17025_tumor_vs_normal.tsv")

  cat("\nDE genes: FDR<", PARAMS$program_de_padj, "=",
      sum(tt$adj.P.Val < PARAMS$program_de_padj),
      "; FDR<", PARAMS$program_de_padj, "and |logFC|>=", PARAMS$program_lfc_cut,
      "=", sum(tt$adj.P.Val < PARAMS$program_de_padj &
                 abs(tt$logFC) >= PARAMS$program_lfc_cut), "\n")
  cat("top 10 up:", paste(head(rownames(tt)[tt$logFC > 0], 10), collapse = ", "), "\n")
  cat("top 10 down:", paste(head(rownames(tt)[tt$logFC < 0], 10), collapse = ", "), "\n")

  # Signature v1: every gene that passes both the FDR and the effect-size cut,
  # with the direction of the change recorded explicitly.
  sig <- rownames(tt)[tt$adj.P.Val < PARAMS$program_de_padj &
                        abs(tt$logFC) >= PARAMS$program_lfc_cut]
  write.table(
    data.frame(gene = sig,
               direction = ifelse(tt[sig, "logFC"] > 0, "UP", "DOWN"),
               logFC = tt[sig, "logFC"],
               FDR = tt[sig, "adj.P.Val"]),
    file.path(TABLES, "EC_signature_v1.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE)
  cat("\nwritten EC_signature_v1.tsv:", length(sig), "genes\n")

  list(tt = tt, ex = ex, grp = grp, sig = sig)
}

# ==========================================================================
# E. GSE196033 -- independent reproduction cohort
# ==========================================================================
ec_reproduction <- function(tt, sig) {
  T <- read_series(file.path(DATA, "GSE196033_series_matrix.txt.gz"))
  exT <- T$ex

  # The GEO field is spelled "tisue" in this series matrix; matched as
  # published rather than corrected, so the field is found by prefix.
  tv <- as.character(T$pd[[grep("^tisue", colnames(T$pd), ignore.case = TRUE)[1]]])
  tv <- sub("^[^:]*: ", "", tv)
  cat("\nGSE196033 groups:\n")
  print(table(tv))

  gT <- factor(ifelse(grepl("cancer", tv), "tumor", "normal"),
               levels = c("normal", "tumor"))

  # The panel arrives as intensities on a log-like scale; log2(x+1) is applied
  # once, then arrays are quantile normalised.
  exT <- log2(exT + 1)
  exT <- normalizeBetweenArrays(exT, method = "quantile")
  fitT <- eBayes(lmFit(exT, model.matrix(~ gT)), trend = TRUE, robust = TRUE)
  ttT <- topTable(fitT, coef = "gTtumor", number = Inf, sort.by = "P")
  write_de(ttT, "DE_GSE196033_tumor_vs_normal.tsv")

  # Concordance between the discovery signature and this cohort. Only the
  # genes measured on both platforms can be compared.
  ov <- intersect(sig, rownames(ttT))
  cat("\ndiscovery signature:", length(sig), "genes; overlapping this panel:",
      length(ov), "\n")
  if (length(ov) > 10) {
    ag <- sign(tt[ov, "logFC"])
    bg <- sign(ttT[ov, "logFC"])
    cat("same direction:", round(100 * mean(ag == bg), 1), "% (binomial p =",
        signif(binom.test(sum(ag == bg), length(ov), 0.5)$p.value, 3), ")\n")
    cat("same direction and FDR<", PARAMS$program_de_padj, "here:",
        sum(ag == bg & ttT[ov, "adj.P.Val"] < PARAMS$program_de_padj), "\n")
  }
  invisible(ttT)
}

# ==========================================================================
# F. GSE17025 with proliferation as a covariate (signature v2)
# ==========================================================================
ec_proliferation_adjusted <- function(sig_unadj) {
  S <- read_series(file.path(DATA, "GSE17025_series_matrix.txt.gz"))
  annot <- read.delim(file.path(DATA, "GPL570_annot.tsv"), check.names = FALSE)
  ex <- collapse_to_gene(S$ex, annot)

  histv <- sub("histology: ",
               "", as.character(S$pd[[grep("^histology", colnames(S$pd))[1]]]))
  grp <- ifelse(grepl("endometrioid|papillary serous", histv), "tumor",
                ifelse(grepl("endometrium", histv), "normal", NA))
  keep <- !is.na(grp)
  ex <- normalizeBetweenArrays(ex[, keep, drop = FALSE], method = "quantile")
  grp <- factor(grp[keep], levels = c("normal", "tumor"))

  # ---- per-sample proliferation score: mean z of the cell-cycle panel ------
  cc <- intersect(CYCLING_GENES, rownames(ex))
  z <- t(scale(t(ex[cc, , drop = FALSE])))     # z-score each gene across samples
  prolif <- colMeans(z, na.rm = TRUE)
  cat("\nproliferation score: normal mean", round(mean(prolif[grp == "normal"]), 2),
      "; tumour mean", round(mean(prolif[grp == "tumor"]), 2), "\n")
  cat("  group t-test p =", signif(t.test(prolif ~ grp)$p.value, 3), "\n")

  # ---- adjusted model: group + proliferation ------------------------------
  # Any tumour effect that survives conditioning on the proliferation score is
  # not attributable to the cell-cycle expansion that accompanies malignancy.
  design <- model.matrix(~ grp + prolif)
  fit <- eBayes(lmFit(ex, design), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = "grptumor", number = Inf, sort.by = "P")
  write_de(tt, "DE_GSE17025_tumor_vs_normal_prolifadj.tsv")

  hit <- intersect(cc, rownames(tt))
  insig <- hit[tt[hit, "adj.P.Val"] > PARAMS$program_de_padj]
  cat("\ncell-cycle genes still at FDR<", PARAMS$program_de_padj, "after adjustment:",
      length(hit) - length(insig), "/", length(hit), "\n")

  sig <- rownames(tt)[tt$adj.P.Val < PARAMS$program_de_padj &
                        abs(tt$logFC) >= PARAMS$program_lfc_cut]
  cat("adjusted DE genes (FDR<", PARAMS$program_de_padj, "and |logFC|>=",
      PARAMS$program_lfc_cut, "):", length(sig),
      "(UP", sum(tt[sig, "logFC"] > 0), "/ DOWN", sum(tt[sig, "logFC"] < 0), ")\n")
  cat("of which cell-cycle:", sum(sig %in% CYCLING_GENES), "\n")

  # ---- compact form: the strongest COMPACT_PER_DIRECTION genes per side ----
  s <- tt[sig, , drop = FALSE]
  s$gene <- rownames(s)
  s$direction <- ifelse(s$logFC > 0, "UP", "DOWN")
  pick <- function(d) {
    head(s[s$direction == d, ][order(-abs(s$logFC[s$direction == d])), ],
         COMPACT_PER_DIRECTION)
  }
  comp <- rbind(pick("UP"), pick("DOWN"))
  write.table(
    data.frame(gene = comp$gene, direction = comp$direction,
               logFC = comp$logFC, FDR = comp$adj.P.Val),
    file.path(TABLES, "EC_signature_v2_prolifadj_compact.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE)

  # ---- full form: every adjusted DE gene --------------------------------
  write.table(
    data.frame(gene = s$gene, direction = s$direction,
               logFC = s$logFC, FDR = s$adj.P.Val),
    file.path(TABLES, "EC_signature_v2_prolifadj_full.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE)

  cat("\nwritten EC_signature_v2_prolifadj_{compact,full}.tsv\n")
  cat("cell-cycle genes in the compact signature:", sum(comp$gene %in% CYCLING_GENES),
      "\n")

  # Reported for the manuscript: how much the adjustment changed the signature.
  cat("signature size before adjustment:", length(sig_unadj),
      "after adjustment:", length(sig), "\n")
  invisible(tt)
}

# ---------------------------------------------------------------- driver ----
main <- function() {
  cat("repository:", ROOT, "\ndata:", DATA, "\ntables:", TABLES, "\n\n")

  ttA <- insulin_response_muscle()
  ttB <- insulin_response_adipocyte()
  insulin_direction_agreement(ttA, ttB)

  # Cache both insulin fits so that a rerun of the programme builder does not
  # require re-parsing the series matrices.
  saveRDS(list(ttA = ttA, ttB = ttB),
          file.path(CACHE, "ir_de_cache.rds"))
  cat("\ncached both insulin fits\n")

  disc <- ec_discovery()
  ec_reproduction(disc$tt, disc$sig)
  ec_proliferation_adjusted(disc$sig)

  cat("\ndone\n")
}

main()
