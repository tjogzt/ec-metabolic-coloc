# config.R
# ========
#
# Purpose
# -------
# R counterpart of config/config.py. Every R script sources this file and must
# not embed paths, accessions or thresholds of its own.
#
# Author  : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
# Created : 2026-11-11
# Version : 2.0.0 (replication package restructure)
#
# Inputs
# ------
# Environment variables (all optional):
#   ECIR_ROOT, ECIR_DATA, ECIR_RESULTS, ECIR_SEED
#
# Outputs
# -------
# ROOT, DATA, RESULTS, TABLES, FIGURES, SEED, PARAMS (named list), COHORTS,
# TISSUES, MODULES, PALETTE
#
# Usage
# -----
#   source(file.path("config", "config.R"))
#   x <- read.delim(file.path(TABLES, "coloc_all_tests.tsv"))

# ---------------------------------------------------------------- paths ----
.this_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)),
                      error = function(e) getwd())
.default_root <- normalizePath(file.path(.this_dir, ".."), mustWork = FALSE)

ROOT    <- normalizePath(Sys.getenv("ECIR_ROOT", unset = .default_root), mustWork = FALSE)
DATA    <- normalizePath(Sys.getenv("ECIR_DATA", unset = file.path(ROOT, "data")), mustWork = FALSE)
RESULTS <- normalizePath(Sys.getenv("ECIR_RESULTS", unset = file.path(ROOT, "results")), mustWork = FALSE)
TABLES  <- file.path(RESULTS, "tables")
FIGURES <- file.path(RESULTS, "figures")

# ------------------------------------------------------- master seed ----
# The only constant fixed in source; everything else is configurable.
SEED <- as.integer(Sys.getenv("ECIR_SEED", unset = "49"))
set.seed(SEED)

# ------------------------------------------------------------ parameters ----
PARAMS <- list(
  # stage 1: state programme
  program_fdr          = as.numeric(Sys.getenv("ECIR_PROGRAM_FDR", "0.05")),
  program_lfc_cut      = as.numeric(Sys.getenv("ECIR_PROGRAM_LFC", "1.0")),
  program_de_padj      = as.numeric(Sys.getenv("ECIR_PROGRAM_DE_PADJ", "0.05")),
  program_min_cohorts  = as.integer(Sys.getenv("ECIR_PROGRAM_MIN_COHORTS", "2")),
  n_random_modules     = as.integer(Sys.getenv("ECIR_N_RANDOM_MODULES", "100")),

  # stage 2: colocalisation
  gwas_p_thresh        = as.numeric(Sys.getenv("ECIR_GWAS_P", "5e-8")),
  locus_gap_bp         = as.integer(Sys.getenv("ECIR_LOCUS_GAP", "1000000")),
  region_pad_bp        = as.integer(Sys.getenv("ECIR_REGION_PAD", "600000")),
  eqtl_min_snps        = as.integer(Sys.getenv("ECIR_EQTL_MIN_SNPS", "10")),
  eqtl_maf_min         = as.numeric(Sys.getenv("ECIR_EQTL_MAF_MIN", "0.01")),
  eqtl_maf_max         = as.numeric(Sys.getenv("ECIR_EQTL_MAF_MAX", "0.99")),
  p12                  = as.numeric(Sys.getenv("ECIR_P12", "1e-5")),
  p12_sensitivity      = c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4),
  m12                  = as.numeric(Sys.getenv("ECIR_M12", "1e-5")),
  pp_h4_high           = as.numeric(Sys.getenv("ECIR_PP_H4_HIGH", "0.8")),
  pp_h3_high           = as.numeric(Sys.getenv("ECIR_PP_H3_HIGH", "0.8")),

  # stage 3: burden
  n_bootstrap          = as.integer(Sys.getenv("ECIR_N_BOOTSTRAP", "200")),
  n_null_draws         = as.integer(Sys.getenv("ECIR_N_NULL_DRAWS", "1000")),
  alpha                = as.numeric(Sys.getenv("ECIR_ALPHA", "0.05")),

  # stage 4: single cell
  sc_min_cells_per_sample = as.integer(Sys.getenv("ECIR_SC_MIN_CELLS", "50")),
  sc_min_samples          = as.integer(Sys.getenv("ECIR_SC_MIN_SAMPLES", "3")),

  # stage 6: cross-cancer
  n_random_partitions  = as.integer(Sys.getenv("ECIR_N_RANDOM_PARTITIONS", "100")),
  min_samples_tumour   = as.integer(Sys.getenv("ECIR_MIN_SAMPLES_TUMOUR", "50")),
  min_samples_normal   = as.integer(Sys.getenv("ECIR_MIN_SAMPLES_NORMAL", "20"))
)

# ---------------------------------------------------- external resources ----
GWAS_PRIMARY_FILE   <- Sys.getenv("ECIR_GWAS_PRIMARY",   "GCST90454186.h.tsv.gz")
GWAS_SECONDARY_FILE <- Sys.getenv("ECIR_GWAS_SECONDARY", "GCST006465_GRCh37.tsv.gz")
GWAS_FINNGEN_FILE   <- Sys.getenv("ECIR_GWAS_FINNGEN",
                                  "finngen_R12_C3_CORPUS_UTERI_EXALLC.gz")
EQTL_CATALOGUE_BASE <- Sys.getenv(
  "ECIR_EQTL_BASE",
  "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015")

# -------------------------------------------- eQTL Catalogue by tissue ----
# Verified against the eQTL Catalogue dataset metadata table
# (data_tables/dataset_metadata_r7.tsv); every entry is GTEx v8, quant_method = ge.
TISSUES <- c(
  Uterus                   = "QTD000346",
  Liver                    = "QTD000266",
  Muscle_Skeletal          = "QTD000281",
  Adipose_Visceral_Omentum = "QTD000121",
  Pancreas                 = "QTD000296"
)

# ------------------------------------------------------------- cohorts ----
# phenotype: how insulin resistance was ascertained in each series.
# gold = TRUE marks the hyperinsulinaemic-euglycaemic clamp cohort, which is
# the reference standard against which the other two are scaled.
COHORTS <- data.frame(
  gse       = c("GSE15773", "GSE18732", "GSE13070"),
  tissue    = c("Adipose_Visceral_Omentum", "Muscle_Skeletal", "Muscle_Skeletal"),
  contrast  = c("IR_state_vs_control", "T2D_vs_control", "clamp_IR_vs_sensitive"),
  phenotype = c("clinical grouping", "diagnosed type 2 diabetes",
                "hyperinsulinaemic-euglycaemic clamp"),
  gold      = c(FALSE, FALSE, TRUE),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------- modules ----
# Frozen before any data were inspected.
MODULES <- list(
  proliferation = c("MKI67", "PCNA", "TOP2A", "CCNB1", "CDK1",
                    "BUB1", "AURKA", "PLK1", "CENPE", "RRM2"),
  adipose       = c("ADIPOQ", "LEP", "PLIN1", "FABP4",
                    "LIPE", "CIDEC", "PPARG", "CEBPA")
)

# ------------------------------------------------------------- palette ----
# Colours carry semantic role, not sequence position.
PALETTE <- c(
  uterus   = "#C23531",   # cinnabar     : target organ / gold-standard phenotype
  liver    = "#3D6BA8",   # azurite      : metabolic tissue
  muscle   = "#177CB0",   # indigo       : skeletal muscle
  adipose  = "#7BA7C7",   # pale azurite : visceral adipose tissue
  pancreas = "#9D2933",   # rouge        : pancreas
  neutral  = "#6E6E6E"    # neutral grey : non-significant / reference
)

FIGURE_FONT   <- Sys.getenv("ECIR_FIG_FONT", "Arial")
FIGURE_MIN_PT <- 8.0

# Ensure output directories exist so that no script has to create them.
for (d in c(RESULTS, TABLES, FIGURES)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
