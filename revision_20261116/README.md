# Revision analyses (2026-11-16)

Added in response to external peer review. All scripts are read-only with
respect to the frozen analysis tables; outputs are listed below.

| Script | Output | What it adds |
|---|---|---|
| `proxy_uncertainty.py` | `proxy_uncertainty.tsv` | Gene-block bootstrap 95% intervals for the three proxy-layer rank correlations (5,000 resamples, seed 20261116) and the matched-module null (100 expression-matched random modules, same estimator as the observed statistic) |
| `gse13070_adipose_arm.R` | stdout (values in Supplementary S7.8) | Within-study tissue-by-phenotype comparison: baseline adipose (IR 35 vs NR 6) and muscle (51 vs 18) per-sample programme scores, exact Mann-Whitney tests, sign concordance, Spearman rho |
| `skap1_strand_current.py` | `skap1_strand_current.log` | Strand-ambiguity audit of the 17q21.32 window recomputed on the primary analysis (GCST90454186): 2,805 matched variants, 394 ambiguous (14.0%), rho -0.628/-0.668/-0.683 |

Errata: the `detectable_80pct` column was removed from
`results/tables/incremental_boot.tsv`; it was computed from the null-module
standard deviation and overstates resolvability (the bootstrap interval is the
correct basis; see Supplementary S7.9).
