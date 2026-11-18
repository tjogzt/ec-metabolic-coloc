#!/usr/bin/env Rscript
# 02b_state_programme_de.R
# =======================
#
# Purpose
# -------
# Rebuild the four differential-expression tables that the locked 862-gene
# programme and the proxy-validity ladder consume:
#
#   GSE15773  adipose, insulin resistant vs sensitive (paired by subject)
#   GSE18732  skeletal muscle, type 2 diabetes vs normal glucose tolerance
#   GSE13070  skeletal muscle, baseline clamp-classified IR vs sensitive
#   GSE25401  adipose, obese vs non-obese (discovery set for the 862 programme)
#
# Inputs
# ------
# <DATA>/GSE15773_series_matrix.txt.gz
# <DATA>/GSE18732_series_matrix.txt.gz
# <DATA>/GSE13070_series_matrix.txt.gz
# <DATA>/GSE25401_series_matrix.txt.gz
# <DATA>/GPL570_annot.tsv
# <DATA>/GPL6244.annot.gz
# <DATA>/ENST2SYMBOL.tsv
# config/config.R   supplies DATA, TABLES, SEED and PARAMS
#
# Outputs
# -------
# <TABLES>/DE_GSE15773_adipose_IR_state.tsv
# <TABLES>/DE_GSE18732_muscle_DM_state.tsv
# <TABLES>/DE_GSE13070_muscle_baseline_IR.tsv
# <TABLES>/DE_GSE25401_adipose_obese_state.tsv
#
# Notes
# -----
# The models mirror the preregistered analysis. GSE15773 is paired: each
# subject contributes an omental and a subcutaneous sample, so the subject is
# the blocking unit and depot is a covariate. GSE18732 has three groups and the
# preregistration specifies DM vs NGT only, so IGT samples are excluded.
# GSE13070 uses the pre-TZD, pre-clamp baseline muscle samples; the +/- suffix
# in those titles denotes treatment response, not insulin sensitivity, so the
# contrast is IRd versus NRd. GSE25401 is the obesity discovery set.

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
})

# ---------------------------------------------------------------- setup ----
script_file <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit)) return(normalizePath(sub("^--file=", "", hit[1]), mustWork = FALSE))
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile)) return(normalizePath(ofile, mustWork = FALSE))
  NULL
}

here <- script_file()
env_root <- Sys.getenv("ECIR_ROOT", unset = "")
if (nzchar(env_root)) {
  root <- normalizePath(env_root, mustWork = TRUE)
} else if (!is.null(here)) {
  root <- normalizePath(file.path(dirname(here), "..", ".."), mustWork = TRUE)
} else {
  stop("cannot locate the repository root: set ECIR_ROOT")
}
source(file.path(root, "config", "config.R"))

options(timeout = 900)

# ------------------------------------------------------------- helpers ----
read_series <- function(path) {
  es <- getGEO(filename = path, getGPL = FALSE)
  ex <- exprs(es)
  pd <- pData(es)
  colnames(ex) <- pd$title
  list(ex = ex, pd = pd)
}

symbol_column <- function(annot) {
  candidates <- c("Gene Symbol", "Gene symbol", "Symbol", "GENE_SYMBOL", "GENE")
  hit <- candidates[candidates %in% colnames(annot)]
  if (!length(hit)) stop("annotation table has no gene-symbol column")
  hit[1]
}

collapse <- function(ex, annot) {
  symcol <- symbol_column(annot)
  a <- annot[match(rownames(ex), annot$ID), c("ID", symcol)]
  ok <- !is.na(a[[symcol]]) & a[[symcol]] != "" & a[[symcol]] != "---"
  ex <- ex[ok, , drop = FALSE]
  a <- a[ok, , drop = FALSE]

  # Multi-symbol entries ("A /// B") must be split; otherwise the composite
  # string stays in the index and genes annotated only in composite form are
  # missing entirely.
  parts <- strsplit(trimws(as.character(a[[symcol]])), "\\s*///\\s*")
  nrep <- lengths(parts)
  if (any(nrep > 1)) {
    ex <- ex[rep(seq_len(nrow(ex)), nrep), , drop = FALSE]
    a <- a[rep(seq_len(nrow(a)), nrep), , drop = FALSE]
    a[[symcol]] <- unlist(parts, use.names = FALSE)
    cat(sprintf("  [collapse] multi-symbol rows %d -> %d\n",
                sum(nrep > 1), sum(nrep)))
  }

  # Raw fluorescence intensities need a log transform before limma.
  if (max(ex, na.rm = TRUE) > 50) {
    cat(sprintf("  [transform] raw intensity (max %.0f) -> log2(x+1)\n",
                max(ex, na.rm = TRUE)))
    ex <- log2(ex + 1)
  }

  m <- rowMeans(ex, na.rm = TRUE)
  o <- order(a[[symcol]], -m)
  a <- a[o, ]
  ex <- ex[o, , drop = FALSE]
  first <- !duplicated(a[[symcol]])
  ex <- ex[first, , drop = FALSE]
  rownames(ex) <- a[[symcol]][first]
  ex
}

read_annot_geo <- function(path) {
  ln <- readLines(gzfile(path), warn = FALSE)
  i <- grep("^ID\t", ln)[1]
  if (is.na(i)) stop("annotation file has no ID<TAB> header")
  read.delim(text = ln[i:length(ln)], check.names = FALSE, quote = "")
}

report <- function(tt, label, suffix = "_state") {
  cat(sprintf("  FDR<0.05 and |logFC|>=0.3: %d\n",
              sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 0.3)))
  dest <- file.path(TABLES, paste0("DE_", label, suffix, ".tsv"))
  write.table(tt, dest, sep = "\t", quote = FALSE)
  cat(sprintf("  written %s\n", dest))
}

# ------------------------------------- GSE15773 adipose, paired design ----
S <- read_series(file.path(DATA, "GSE15773_series_matrix.txt.gz"))
annot570 <- read.delim(file.path(DATA, "GPL570_annot.tsv"), check.names = FALSE)
ex <- collapse(S$ex, annot570)
ttl <- colnames(ex)
grp <- ifelse(grepl("Insulin resistant", ttl, ignore.case = TRUE), "IR",
       ifelse(grepl("Insulin sensitive", ttl, ignore.case = TRUE), "IS", NA))
depot <- ifelse(grepl("omental", ttl, ignore.case = TRUE), "omental",
         ifelse(grepl("subcutaneous", ttl, ignore.case = TRUE), "subcut", NA))
keep <- !is.na(grp) & !is.na(depot)
ex <- normalizeBetweenArrays(ex[, keep, drop = FALSE], method = "quantile")
g <- factor(grp[keep], levels = c("IS", "IR"))
dp <- factor(depot[keep])
cat(sprintf("\n=== GSE15773 adipose ===\ngenes %d, samples %d (IR %d / IS %d; omental %d / subcut %d)\n",
            nrow(ex), sum(keep), sum(g == "IR"), sum(g == "IS"),
            sum(dp == "omental"), sum(dp == "subcut")))

# Each subject contributes two chips; the subject is the block.
subj <- paste(g, sub(".*patient\\s*([0-9]+).*", "\\1", ttl[keep], ignore.case = TRUE), sep = "_")
design <- model.matrix(~ g + dp)
dup <- duplicateCorrelation(ex, design, block = subj)
cons <- dup$consensus
if (!is.finite(cons) || cons < 0) {
  cons <- 0
  cat("  [warning] non-positive within-subject correlation set to 0\n")
}
cat(sprintf("  within-subject consensus = %.3f\n", cons))
fit <- eBayes(lmFit(ex, design, block = subj, correlation = cons),
              trend = TRUE, robust = TRUE)
tt_a <- topTable(fit, coef = "gIR", number = Inf, sort.by = "P")
report(tt_a, "GSE15773_adipose_IR")

# ------------------------------- GSE18732 muscle, DM vs NGT (no IGT) ----
S2 <- read_series(file.path(DATA, "GSE18732_series_matrix.txt.gz"))
enst <- read.delim(file.path(DATA, "ENST2SYMBOL.tsv"), check.names = FALSE)
miss <- !rownames(S2$ex) %in% enst$ID
cat(sprintf("  ENST unmapped %d/%d (%.1f%%)\n",
            sum(miss), length(miss), 100 * mean(miss)))
ex2 <- collapse(S2$ex, enst)
ex2 <- ex2[complete.cases(ex2) & is.finite(rowSums(ex2)), , drop = FALSE]
ttl2 <- colnames(ex2)
g2 <- ifelse(grepl("diabetic", ttl2, ignore.case = TRUE), "DM",
      ifelse(grepl("normal", ttl2, ignore.case = TRUE), "NGT", NA))
keep2 <- !is.na(g2)
ex2 <- normalizeBetweenArrays(ex2[, keep2, drop = FALSE], method = "quantile")
g2 <- factor(g2[keep2], levels = c("NGT", "DM"))
cat(sprintf("\n=== GSE18732 muscle ===\ngenes %d, samples %d (DM %d / NGT %d; IGT excluded)\n",
            nrow(ex2), sum(keep2), sum(g2 == "DM"), sum(g2 == "NGT")))
fit2 <- eBayes(lmFit(ex2, model.matrix(~ g2)), trend = TRUE, robust = TRUE)
tt_b <- topTable(fit2, coef = "g2DM", number = Inf, sort.by = "P")
report(tt_b, "GSE18732_muscle_DM")

# --------------------------- GSE13070 muscle, baseline IRd vs NRd ----
read_plain_series <- function(path) {
  ln <- readLines(gzfile(path), warn = FALSE)
  i0 <- grep("^!series_matrix_table_begin", ln)
  i1 <- grep("^!series_matrix_table_end", ln)
  pd0 <- grep("^!Sample_title", ln, value = TRUE)
  ttl <- strsplit(pd0[1], "\t")[[1]][-1]
  ttl <- gsub('"', "", ttl)
  tb <- read.delim(text = ln[(i0 + 1):(i1 - 1)], check.names = FALSE, quote = "\"")
  ex <- as.matrix(tb[, -1])
  rownames(ex) <- tb[[1]]
  colnames(ex) <- ttl
  list(ex = ex, title = ttl)
}

S3 <- read_plain_series(file.path(DATA, "GSE13070_series_matrix.txt.gz"))
ex3 <- collapse(S3$ex, annot570)
ttl3 <- colnames(ex3)
keep3 <- grepl("^muscle_preTZD_preClamp", ttl3)
ex3 <- ex3[, keep3, drop = FALSE]
ttl3 <- ttl3[keep3]
grp3 <- ifelse(grepl("_IRd", ttl3), "IR",
        ifelse(grepl("_NRd", ttl3), "NR", NA))
stopifnot(!any(is.na(grp3)))
cat(sprintf("\n=== GSE13070 muscle baseline ===\ngenes %d, samples %d (IR %d / NR %d)\n",
            nrow(ex3), length(grp3), sum(grp3 == "IR"), sum(grp3 == "NR")))
ex3 <- ex3[complete.cases(ex3) & is.finite(rowSums(ex3)), , drop = FALSE]
ex3 <- normalizeBetweenArrays(ex3, method = "quantile")
g3 <- factor(grp3, levels = c("NR", "IR"))
fit3 <- eBayes(lmFit(ex3, model.matrix(~ g3)), trend = TRUE, robust = TRUE)
tt_c <- topTable(fit3, coef = "g3IR", number = Inf, sort.by = "P")
report(tt_c, "GSE13070_muscle_baseline_IR", suffix = "")

# ---------------------------- GSE25401 adipose, obese vs non-obese ----
S4 <- read_series(file.path(DATA, "GSE25401_series_matrix.txt.gz"))
annot6244 <- read_annot_geo(file.path(DATA, "GPL6244.annot.gz"))
ex4 <- collapse(S4$ex, annot6244)
ttl4 <- colnames(ex4)
g4 <- ifelse(grepl("Obese", ttl4, ignore.case = TRUE) &
             !grepl("Nonobese", ttl4, ignore.case = TRUE), "Obese",
      ifelse(grepl("Nonobese", ttl4, ignore.case = TRUE), "Nonobese", NA))
keep4 <- !is.na(g4)
ex4 <- normalizeBetweenArrays(ex4[, keep4, drop = FALSE], method = "quantile")
ex4 <- ex4[complete.cases(ex4) & is.finite(rowSums(ex4)), , drop = FALSE]
g4 <- factor(g4[keep4], levels = c("Nonobese", "Obese"))
cat(sprintf("\n=== GSE25401 adipose ===\ngenes %d, samples %d (Obese %d / Nonobese %d)\n",
            nrow(ex4), sum(keep4), sum(g4 == "Obese"), sum(g4 == "Nonobese")))
fit4 <- eBayes(lmFit(ex4, model.matrix(~ g4)), trend = TRUE, robust = TRUE)
tt_d <- topTable(fit4, coef = "g4Obese", number = Inf, sort.by = "P")
report(tt_d, "GSE25401_adipose_obese")

cat("\ndone\n")
