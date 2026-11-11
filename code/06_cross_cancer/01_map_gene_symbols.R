#!/usr/bin/env Rscript
# 01_map_gene_symbols.R
# =====================
#
# Purpose
# -------
# Resolve the insulin-resistance metabolic programme and the two fixed reference
# modules (proliferation, adipose) onto HUGO gene symbols, and write the symbol
# table that stages 03 and 04 consume.
#
# Gene symbols are the only identifier the TCGA pan-cancer matrix carries: its
# rows are labelled with HUGO symbols, not with Entrez ids. An earlier version of
# this stage mapped symbols to Entrez through org.Hs.eg.db, which silently
# produced an empty intersection downstream. The mapping is therefore done
# symbol-to-symbol here, and the Entrez dependency is gone.
#
# Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
# Created  : 2026-11-11
# Version  : 2.0.0
#
# Inputs
# ------
# results/tables/metabolic_program_v2.tsv
#     The frozen 862-gene programme (608 up, 254 down), one row per gene, gene
#     symbols as row names, with at least the columns ``direction`` (UP/DOWN)
#     and ``logFC_adipose_obese``.
# config/config.R
#     Supplies TABLES, DATA, MODULES and SEED; no path or module list is written
#     into this script.
# data/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz (optional)
#     Only its header and first data row are read, to confirm that the matrix is
#     keyed by gene symbol. The script runs without it and says so.
#
# Outputs
# -------
# results/tables/program_symbols.tsv
#     Columns: SYMBOL, direction, logFC. One row per programme gene, in the row
#     order of the source table.
# results/tables/ref_modules.tsv
#     Columns: module, gene, position. Long format, module membership taken
#     verbatim from config.MODULES; nothing is selected from the data.
#
# Usage
# -----
#   Rscript code/06_cross_cancer/01_map_gene_symbols.R
#   ECIR_ROOT=/path/to/repo Rscript code/06_cross_cancer/01_map_gene_symbols.R
#
# Notes
# -----
# * The programme is already symbol-keyed, so this stage is a validation and
#   normalisation step, not an identifier conversion. It enforces three
#   invariants that later stages rely on: unique symbols, a binary direction
#   label, and non-numeric (non-Entrez) identifiers.
# * Duplicate symbols are resolved by keeping the first occurrence, which is what
#   the previous Entrez-based stage did after its AnnotationDbi lookup.
# * The proliferation and adipose modules are frozen in config before any data
#   were inspected and are never derived from expression in this script.

# ---- locate the repository root -------------------------------------------
# Precedence: ECIR_ROOT, then the location of this script. Nothing is hardcoded.
.script_file <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit)) return(normalizePath(sub("^--file=", "", hit[1]), mustWork = FALSE))
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile)) return(normalizePath(ofile, mustWork = FALSE))
  NULL
}

.here <- .script_file()
.env_root <- Sys.getenv("ECIR_ROOT", unset = "")
if (nzchar(.env_root)) {
  .root <- normalizePath(.env_root, mustWork = TRUE)
} else if (!is.null(.here)) {
  # <root>/code/06_cross_cancer/01_map_gene_symbols.R -> up two levels
  .root <- normalizePath(file.path(dirname(.here), "..", ".."), mustWork = TRUE)
} else {
  stop("cannot locate the repository root: set ECIR_ROOT or run this script as a file")
}

source(file.path(.root, "config", "config.R"))

# ---- load the frozen programme --------------------------------------------
prog_file <- file.path(TABLES, "metabolic_program_v2.tsv")
if (!file.exists(prog_file)) {
  stop(sprintf("programme table not found: %s", prog_file))
}

prog <- read.delim(prog_file, row.names = 1, check.names = FALSE,
                   stringsAsFactors = FALSE)

required <- c("direction", "logFC_adipose_obese")
missing_cols <- setdiff(required, colnames(prog))
if (length(missing_cols)) {
  stop(sprintf("programme table lacks column(s): %s",
               paste(missing_cols, collapse = ", ")))
}

genes <- rownames(prog)
if (anyDuplicated(genes)) {
  # Keep the first occurrence so that the row order of the source table survives.
  keep <- !duplicated(genes)
  message(sprintf("dropping %d duplicated symbol row(s)", sum(!keep)))
  prog <- prog[keep, , drop = FALSE]
  genes <- rownames(prog)
}

# Direction must be a clean binary label; anything else would be silently
# averaged into the wrong arm of the score.
bad_dir <- setdiff(unique(prog$direction), c("UP", "DOWN"))
if (length(bad_dir)) {
  stop(sprintf("unexpected direction label(s): %s",
               paste(bad_dir, collapse = ", ")))
}

# Guard the symbol convention: a purely numeric identifier means an Entrez id
# slipped back in, which would match nothing in the TCGA matrix.
numeric_like <- genes[grepl("^[0-9]+$", genes)]
if (length(numeric_like)) {
  stop(sprintf("%d identifier(s) look like Entrez ids, not gene symbols: %s",
               length(numeric_like),
               paste(utils::head(numeric_like, 5), collapse = ", ")))
}
blank <- genes[!nzchar(trimws(genes))]
if (length(blank)) {
  stop("programme table contains empty gene symbols")
}

n_up   <- sum(prog$direction == "UP")
n_down <- sum(prog$direction == "DOWN")
cat(sprintf("Programme: %d genes (UP %d / DOWN %d)\n", length(genes), n_up, n_down))
cat(sprintf("Expression range of logFC_adipose_obese: [%+.2f, %+.2f]\n",
            min(prog$logFC_adipose_obese, na.rm = TRUE),
            max(prog$logFC_adipose_obese, na.rm = TRUE)))

# ---- fixed reference modules from config ----------------------------------
# Membership is written down in config.MODULES before any data were inspected;
# this stage never chooses module genes from expression values.
if (!all(c("proliferation", "adipose") %in% names(MODULES))) {
  stop("config.MODULES must define 'proliferation' and 'adipose'")
}

module_rows <- do.call(rbind, lapply(names(MODULES), function(m) {
  members <- MODULES[[m]]
  data.frame(module = rep(m, length(members)),
             gene = members,
             position = seq_along(members),
             stringsAsFactors = FALSE)
}))

# A module gene that also sits in the programme would be scored twice, once in
# each direction; report it rather than quietly allowing the overlap.
overlap <- intersect(unique(module_rows$gene), genes)
if (length(overlap)) {
  message(sprintf("module genes also present in the programme: %s",
                  paste(overlap, collapse = ", ")))
}

# ---- matrix identifier check ----------------------------------------------
# Cheap probe: read the header and one data row of the expression matrix and
# confirm the row identifier is a symbol. The full matrix is left to stage 02.
expr_file <- file.path(DATA, "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz")
if (file.exists(expr_file)) {
  con <- gzfile(expr_file, "rt")
  header <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1]]
  first  <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1]]
  close(con)
  cat(sprintf("Matrix: %d columns in the header; first row identifier '%s'\n",
              length(header) - 1L, first[1]))
  if (grepl("^[0-9]+$", first[1])) {
    stop(sprintf("matrix is keyed by numeric ids ('%s'), but this package matches on symbols",
                 first[1]))
  }
} else {
  message("expression matrix not present locally; skipping the identifier check")
  message(sprintf("  expected at %s", expr_file))
}

# ---- write outputs ---------------------------------------------------------
out_prog <- data.frame(SYMBOL = genes,
                       direction = prog$direction,
                       logFC = prog$logFC_adipose_obese,
                       stringsAsFactors = FALSE)
write.table(out_prog, file.path(TABLES, "program_symbols.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(module_rows, file.path(TABLES, "ref_modules.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Modules: proliferation %d genes | adipose %d genes (fixed in config)\n",
            length(MODULES$proliferation), length(MODULES$adipose)))
cat(sprintf("Wrote %s (%d rows) and %s (%d rows)\n",
            file.path(TABLES, "program_symbols.tsv"), nrow(out_prog),
            file.path(TABLES, "ref_modules.tsv"), nrow(module_rows)))
