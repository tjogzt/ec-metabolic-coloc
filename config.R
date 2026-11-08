# Path configuration for the replication package.
# Run scripts from the repository root, or set ECIR_ROOT.
ROOT    <- Sys.getenv("ECIR_ROOT", unset = normalizePath("."))
DATA    <- file.path(ROOT, "data")
RESULTS <- file.path(ROOT, "results")
TABLES  <- file.path(RESULTS, "tables")
FIGURES <- file.path(RESULTS, "figures")
