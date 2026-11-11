# 03_annotation_lookup.R
# ======================
#
# Purpose
# -------
# Download the platform annotation tables that every probe-to-gene collapse in
# the bulk layer depends on, and write them as plain TSV next to the expression
# data.
#
# Only the identifier columns are kept. The full GPL record carries tens of
# thousands of probe sequences and is neither reproducible across releases nor
# useful here; what the collapse needs is the probe id plus whichever column
# holds the gene symbol, and that column is named differently on every platform
# ("Gene Symbol" on GPL570, "GENE_SYMBOL" on GPL6480), so several candidates are
# retained and the choice is made downstream by 02_differential_expression.R.
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
# Optional positional arguments: platform accessions to fetch. The default is
# the two platforms the bulk cohorts are annotated on.
#
# Outputs
# -------
# <DATA>/<GPL>_annot.tsv   one row per probe, the ID column first followed by
#                          every identifier column the platform provides.
#
# Usage
# -----
#   Rscript code/03_burden/03_annotation_lookup.R
#   Rscript code/03_burden/03_annotation_lookup.R GPL570
#
# Notes
# -----
# A failed download is reported and skipped, not raised: the caller only needs
# the platforms that its cohorts actually use, and a single retired platform
# should not take the whole stage down. The script exits non-zero if every
# platform failed, so a silent empty run is impossible.

# ---------------------------------------------------------------- setup ----
suppressMessages(library(GEOquery))

# Network timeout for the platform tables, which run to tens of megabytes.
HTTP_TIMEOUT <- 900
options(timeout = HTTP_TIMEOUT)

# Platforms used by the bulk cohorts: GPL570 for GSE17025 / GSE196033 /
# GSE24422, GPL6480 for GSE24215.
DEFAULT_PLATFORMS <- c("GPL570", "GPL6480")

# Identifier columns worth keeping when a platform provides them.
KEEP_COLUMNS <- c("ID", "GB_ACC", "GENE_SYMBOL", "Gene Symbol", "Symbol",
                  "GENE", "ORF", "SEQUENCE")

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

# ---------------------------------------------------------------- helpers ----
fetch_platform <- function(gpl) {
  # Fetch one platform table and write the identifier subset.
  #
  # Returns TRUE on success and FALSE when the download failed, so the driver
  # can decide whether the run as a whole was usable.
  message("fetching ", gpl, " ...")
  g <- tryCatch(getGEO(gpl, destdir = DATA),
                error = function(e) {
                  message("  failed: ", conditionMessage(e))
                  NULL
                })
  if (is.null(g)) return(FALSE)

  tb <- Table(g)
  if (is.null(tb) || nrow(tb) == 0) {
    message("  ", gpl, " carries no annotation table")
    return(FALSE)
  }
  cat(gpl, "columns:", paste(head(colnames(tb), 8), collapse = " | "), "\n")
  cat(gpl, "rows:", nrow(tb), "\n")

  keep <- intersect(KEEP_COLUMNS, colnames(tb))
  cat(gpl, "identifier columns:", paste(keep, collapse = " | "), "\n\n")

  # "ID" is prepended in the original order used by the collapse functions;
  # unique() keeps it from being written twice when it also appears in keep.
  dest <- file.path(DATA, paste0(gpl, "_annot.tsv"))
  write.table(tb[, unique(c("ID", keep)), drop = FALSE],
              file = dest, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("written", dest, "\n")
  TRUE
}

# ---------------------------------------------------------------- driver ----
main <- function(argv) {
  platforms <- if (length(argv)) argv else DEFAULT_PLATFORMS

  # The annotation files must sit beside the series matrices, which the config
  # guarantees only after ensure() has run; create the directory if needed.
  if (!dir.exists(DATA)) dir.create(DATA, recursive = TRUE, showWarnings = FALSE)

  cat("repository:", ROOT, "\ndata:", DATA, "\nplatforms:", paste(platforms, collapse = ", "),
      "\n\n")
  ok <- vapply(platforms, fetch_platform, logical(1))

  if (!any(ok)) stop("no platform annotation could be retrieved")
  if (any(!ok)) {
    message("skipped: ", paste(platforms[!ok], collapse = ", "))
  }
  invisible(NULL)
}

main(commandArgs(trailingOnly = TRUE))
