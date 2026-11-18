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

## Stage 01 — the programme definitions

Implemented in `code/03_burden/`. Two definitions are kept.

**Reported definition (obesity-associated, 862 genes).** The discovery contrast is
GSE25401 adipose obese versus non-obese. The four DE tables used downstream are
built by `02b_state_programme_de.R`:

```bash
Rscript code/03_burden/02b_state_programme_de.R
python  code/03_burden/05b_build_programme_v2.py
```

`05b_build_programme_v2.py` writes the locked 862-gene table
(`metabolic_program_v2.tsv`, 608 up / 254 down) and its 100 matched modules. The
current DE table yields a 894-gene rebuild, because the DE table was regenerated
after the locked programme was produced. The builder therefore writes
`metabolic_program_v2_rebuild.tsv` and `random_modules_v2_rebuild.tsv` and exits
2 when the gene set differs; it never overwrites the locked table. The 48-gene
difference is documented in `06_EC_IR/w4/dryrun_894.py`.

**Earlier definition (insulin-response, 146 genes).** Built from GSE24215 and
GSE24422:

```bash
Rscript code/03_burden/02_differential_expression.R
python  code/03_burden/05_build_programme.py
```

This writes `IR_program_v1.tsv` and `random_modules_v1.tsv` and is retained for
provenance; it is not the definition reported in the paper.

`01_concordance_pipeline.py` computes rank correlation between the locked
programme score and phenotype in each cohort, which is the proxy-validity ladder
reported in the paper: +0.31 (adipose, clinical grouping), +0.20 (muscle,
diagnosed T2D), +0.10 (muscle, clamp). The ladder is monotone in the strictness
of the phenotype, which is the finding. Its inputs are the three DE tables
written by `02b_state_programme_de.R`.

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

The defaults are the definitions reported in the paper: the locked 862-gene
programme (`metabolic_program_v2.tsv`) and its matched random modules
(`random_modules_v2.tsv`). The earlier insulin-response definition is still
reproducible by passing `--programme IR_program_v1.tsv --modules
random_modules_v1.tsv --tag _v1`.

The bootstrap resamples patients, never cells or probes; the statistic is a
difference in cross-validated AUC, so a positive point estimate with a confidence
interval spanning zero is reported as "no evidence of incremental value", never
as equivalence. The `delta_auc` column of the bootstrap table is one
cross-validation realisation used to centre the interval; the reported point
estimate is `delta_auc_IR` in `incremental_value.tsv`, and the two need not have
the same sign.

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

The retrieved Xena matrix is keyed by Entrez id. `01_map_gene_symbols.R` writes
`data/ENTREZ2SYMBOL.tsv` from `org.Hs.eg.db`, and stage 02 applies it while
streaming, collapsing several Entrez rows onto one symbol by their mean. This is
a reconstruction of the symbol-keyed matrix used for the paper, not the original
object; the scoreable programme fraction and the correlation summaries can shift
slightly.

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
