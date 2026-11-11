# Analysis flow

This document describes what each stage consumes and produces, in the order the
stages must run. Every path below is resolved from `config/config.py` (Python)
or `config/config.R` (R); no script contains a hard-coded location.

Run every command from the repository root.

```
                       ┌──────────────────────────┐
                       │ config/config.{py,R}     │
                       │ paths, priors, cohorts   │
                       └────────────┬─────────────┘
                                    │
        ┌───────────────────────────┼────────────────────────────┐
        ▼                           ▼                            ▼
┌─────────────────────────┐          ┌─────────────────┐
│ 03_burden               │          │ 02_colocalisation│
│ cohorts → 862-gene      │          │ five tissues,    │
│ programme → incremental │          │ PP.H0-H4         │
└───────┬─────────────────┘          └────────┬─────────┘
        │                           │                            │
        └───────────┬───────────────┴────────────┬───────────────┘
                    ▼                            ▼
            ┌───────────────┐            ┌──────────────────┐
            │ 04_singlecell │            │ 05_mendelian_    │
            │ projection    │            │ randomisation    │
            └───────┬───────┘            └────────┬─────────┘
                    │                             │
                    ▼                             ▼
            ┌───────────────┐            ┌──────────────────┐
            │ 06_cross_     │            │ 07_figures       │
            │ cancer        │            │ Figures 1-4      │
            └───────────────┘            └──────────────────┘
```

## Stage 00 — inputs

Nothing is downloaded automatically. `data/README.md` lists each accession, the
exact file required and the name the scripts expect. Place the files under
`$ECIR_DATA` (default `<repo>/data`).

Verify the tree before running anything:

```bash
python code/00_setup/01_check_inputs.py
```

## Stage 01 — the 862-gene state programme

Implemented in `code/03_burden/`, because the programme is derived from the same
cohort data that later carries the burden analyses.

Inputs: GEO series GSE15773 (adipose), GSE18732 (muscle, diagnosed type 2
diabetes), GSE13070 (muscle, clamp-defined insulin resistance).

The programme is derived from differential expression between insulin-resistant
and insulin-sensitive samples, taken in the direction that is reproducible across
cohorts. Outcome labels are never used, so the programme cannot be tuned towards
the endometriotic cancer endpoint.

Output: `results/tables/metabolic_program_v2.tsv` (608 up, 254 down).

```bash
Rscript code/03_burden/02_differential_expression.R
python  code/03_burden/05_build_programme.py
```

`01_concordance_pipeline.py` computes rank correlation between the programme
score and phenotype in each cohort, which is the proxy-validity ladder reported
in the paper: +0.31 (adipose, clinical grouping), +0.20 (muscle, diagnosed T2D),
+0.10 (muscle, clamp). The ladder is monotone in the strictness of the phenotype,
which is the finding.

## Stage 02 — colocalisation

Inputs: risk GWAS `GCST90454186` (hg38) and eQTL Catalogue GTEx v8 datasets for
uterus, liver, skeletal muscle, visceral adipose and pancreas.

```bash
python code/02_colocalisation/07_run_colocalisation.py --tissue Uterus
python code/02_colocalisation/07_run_colocalisation.py --tissue Liver
python code/02_colocalisation/07_run_colocalisation.py --tissue Muscle_Skeletal
python code/02_colocalisation/07_run_colocalisation.py --tissue Pancreas
python code/02_colocalisation/07_run_colocalisation.py --tissue Adipose_Visceral_Omentum --remote
```

The adipose dataset is not mirrored in this package; `--remote` queries the EBI
tabix endpoint for just the locus windows.

Two follow-up stages constrain how the result may be stated:

```bash
python code/02_colocalisation/09_dual_outcome.py      # repeat against GCST006465
python code/02_colocalisation/08_prior_sensitivity.py # sweep the p12 prior
```

`08_prior_sensitivity.py` requires the locus subsets cached by stage 07, so it
must run after all five tissues.

## Stage 03 — burden and incremental value

Inputs: GSE120490 (145 endometrioid tumours, 45 events) plus the clinical base
model.

Output: `results/tables/incremental_value.tsv` and `incremental_boot.tsv`.

```bash
python code/03_burden/07_incremental_value.py
```

The bootstrap resamples patients, never cells or probes; the statistic is a
difference in cross-validated AUC, so a positive point estimate with a confidence
interval spanning zero is reported as "no evidence of incremental value", never
as equivalence.

## Stage 04 — single-cell projection

Input: the integrated endometrial single-cell atlas (11 samples, 141,316 cells).

```bash
python code/04_singlecell/01_project_programme.py
python code/04_singlecell/02_null_calibration.py
```

The statistical unit is the sample. Scoring individual cells and testing the cell
level difference inflates significance through pseudo-replication, so the
comparison is made on per-sample means.

## Stage 05 — Mendelian randomisation

Inputs: fasting insulin (GCST90002238) and BMI (GCST006900) exposures, risk
outcomes GCST90454186 and GCST006465.

```bash
python code/05_mendelian_randomisation/02_select_instruments.py
python code/05_mendelian_randomisation/03_harmonise_datasets.py
python code/05_mendelian_randomisation/04_run_mr.py
python code/05_mendelian_randomisation/05_conditional_f.py
python code/05_mendelian_randomisation/06_second_outcome.py
```

Two guards matter here. Instruments are restricted to variants with N >= 100,000,
which isolates a homogeneous block of small-sample variants that would otherwise
enter as weak instruments. Conditional F statistics are reported alongside every
multivariable estimate; a conditional F below 10 means the decomposition is not
interpretable, and the package prints that verdict rather than the coefficients
alone.

## Stage 06 — cross-cancer context

Inputs: TCGA PanCan expression matrix and the 862-gene programme.

```bash
Rscript code/06_cross_cancer/01_map_gene_symbols.R
python  code/06_cross_cancer/02_load_tcga_expression.py
python  code/06_cross_cancer/03_score_across_cancers.py
python  code/06_cross_cancer/04_normal_tissue_null.py
```

Random modules are matched to the real programme on size, on expression decile
and on the up/down split. A null distribution built from unmatched modules is
not a valid comparator and would overstate significance.

## Stage 07 — figures

```bash
Rscript code/07_figures/01_main_figures_1_2_3.R
python  code/07_figures/02_figure1_design.py
python  code/07_figures/04_supplementary_mr_figures.py   # supplementary MR panels
```

Figures are vector PDF with PNG previews. Colours come from `PARAMS`/`PALETTE` in
the config and carry a fixed semantic role (target organ, metabolic tissue,
muscle, adipose, pancreas, neutral), not a series index.

## Reproducing the published tables

`results/tables/` ships the exact tables reported in the paper. Re-running the
stages above regenerates them; the pipeline is deterministic given the inputs
apart from the bootstrap and permutation stages, which are seeded from
`config.SEED` (a single value, the only constant fixed in source).

## Environment

| Component | Version used |
|---|---|
| Python | 3.12, `numpy`, `pandas`, `scipy`, `scanpy`, `anndata` |
| R | 4.5, `limma`, `ggplot2`, `cowplot`, `data.table`, `org.Hs.eg.db` |
| External binaries | `tabix`, `bgzip`, `pdftotext` (figure checks only) |

`TwoSampleMR` and `MRPRESSO` are deliberately not required: the estimators are
implemented in `code/05_mendelian_randomisation/01_mr_estimators.py` so that the
exact formulas used are visible and testable.
