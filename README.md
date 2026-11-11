# Shared variation at 17q21.32 in endometrial cancer and the limits of obesity-derived transcriptional signatures as proxies for insulin resistance

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.24055114-003399.svg)](https://doi.org/10.5281/zenodo.24055114)

Replication package for a preregistered, public-data-only analysis with three
layers: genetic colocalisation, bulk transcriptional validation, and single-cell
projection.

The analysis asks whether obesity-derived transcriptional signatures can be used
to infer that endometrial cancer shares causal metabolic architecture with
obesity and type 2 diabetes, and whether scoring such a signature in tissue
measures insulin resistance. The answer differs by layer.

## What the analysis reports

| Layer | Result |
|---|---|
| Colocalisation | 503 gene–locus tests in five tissues; 3 exceed PP.H4 > 0.8 at the primary prior, but only the **17q21.32** signal survives a change of risk dataset |
| Programme validity | An obesity-derived programme (862 genes) tracks clamp-defined insulin resistance weakly: rank correlation **+0.10** (p = 0.004), 61.6% sign concordance among responsive genes |
| Incremental value | **ΔAUC = +0.0004** over a base clinical model (95% CI −0.014 to +0.031; p = 0.15, 145 tumours, 45 events); no evidence of value, but increments up to ~0.03 not excluded |
| Single cell | 141,316 cells, 29,592 epithelial; programme score does not separate malignant from normal epithelium (z = −3.8, below the prespecified bar) |

The shared signal at 17q21.32 is the only positive genetic finding: at the lead
eQTL variant the allele that lowers *SKAP1* expression raises endometrial cancer
risk (liver, PP.H4 = 0.93; 2,805 variants in window; ρ = −0.66 between expression
and risk z-scores). The same allele shifts *SKAP1* expression by −0.25 in pancreas
and −0.24 in spleen, but not in adipose, uterus, ovary or blood.

Two robustness checks constrain how that finding may be stated, and both are
reported in the paper. First, the credible set at this locus spans **57 variants**,
with the most probable single variant carrying 4% posterior weight, so the
analysis identifies a shared *signal* rather than a variant. Second, posteriors
depend on the prior: across p12 from 1e-6 to 1e-4 the number of tests above 0.8
ranges from 0 to 7, and at the stricter prior none is called. A signal at 15q15.2
that exceeded the threshold in the primary analysis fell to PP.H4 = 0.15 with
PP.H3 = 0.85 against the preregistered risk dataset, and is not reported as a
finding. The two risk datasets share the Endometrial Cancer Association
Consortium, so agreement between them is consistency across a change of dataset,
not independent replication.

## Repository layout

```
config/
  config.py             single source of truth: paths, priors, cohorts, modules, palette
  config.R              the same values for the R stages
code/
  00_setup/             input verification
  02_colocalisation/    ABF engine, GWAS/eQTL readers, loci, tests, sensitivity
  03_burden/            cohort assembly, differential expression, 862-gene programme,
                        incremental value
  04_singlecell/        single-cell projection and size-matched null
  05_mendelian_randomisation/   instruments, harmonisation, estimators, conditional F
  06_cross_cancer/      programme scoring across TCGA tumour types
  07_figures/           Figures 1–4 and the supplementary MR panels
docs/
  FLOW.md               stage-by-stage description of the analysis
results/
  script_inputs/        exact inputs the figure scripts consume
  tables/               derived tables reported in the paper
  figures/              Figures 1–4 as vector PDF, PNG previews alongside
data/
  README.md             every input, its accession, and how to obtain it
```

Files are named `<index>_<function>.py` (or `.R`) so that the directory listing
is also the run order within a stage.

## Configuration

Every path and parameter is resolved in `config/config.py` (and `config/config.R`
for the R stages). No analysis script contains a hard-coded path, accession or
threshold. The only constant fixed in source is the random seed.

Override anything through the environment:

```bash
export ECIR_ROOT=/path/to/ec-metabolic-coloc   # default: parent of config/
export ECIR_DATA=/path/to/inputs               # default: $ECIR_ROOT/data
export ECIR_RESULTS=/path/to/outputs           # default: $ECIR_ROOT/results
export ECIR_SEED=49                            # master seed
export ECIR_P12=1e-5                           # shared-variant prior
export ECIR_GWAS_P=5e-8                        # locus discovery threshold
```

Run `python code/00_setup/01_check_inputs.py` to see which inputs are in place;
it prints the resolved data directory and the size of every file it finds.

## Data

No raw data are redistributed. Every input is public and re-downloadable;
`data/README.md` lists each accession, the file needed and the local name the
scripts expect.

The one input that is not mirrored locally is the visceral adipose eQTL dataset
(`QTD000121`); stage 02 fetches only the locus windows from the EBI tabix
endpoint when run with `--remote`.

## Reproducing the analysis

See `docs/FLOW.md` for the full stage-by-stage description. Short version:

```bash
python  code/00_setup/01_check_inputs.py

Rscript code/03_burden/02_differential_expression.R
python  code/03_burden/05_build_programme.py

for t in Uterus Liver Muscle_Skeletal Pancreas; do
  python code/02_colocalisation/07_run_colocalisation.py --tissue "$t"
done
python code/02_colocalisation/07_run_colocalisation.py \
       --tissue Adipose_Visceral_Omentum --remote

python code/02_colocalisation/08_prior_sensitivity.py
python code/02_colocalisation/09_dual_outcome.py

python code/03_burden/07_incremental_value.py
python code/04_singlecell/01_project_programme.py
python code/05_mendelian_randomisation/04_run_mr.py

Rscript code/07_figures/01_main_figures_1_2_3.R
```

Scripts are run from the repository root.

Environment: Python 3.12 (`numpy`, `pandas`, `scipy`, `scanpy`, `anndata`),
R 4.5 (`limma`, `ggplot2`, `cowplot`, `data.table`, `org.Hs.eg.db`), and the
`tabix`/`bgzip` binaries. `TwoSampleMR` and `MRPRESSO` are deliberately not
required — the estimators are implemented in
`code/05_mendelian_randomisation/01_mr_estimators.py` so the formulas used are
visible and testable.

## Notes on the colocalisation layer

- Both the eQTL resource and the primary risk dataset are on GRCh38, so variants
  are matched directly by position and alleles; no liftover is applied. The
  secondary risk dataset is on GRCh37 and is lifted, variant by variant.
- Allele orientation is resolved at every join: when the stored alleles are
  swapped relative to the query, the effect is negated. Skipping this step
  inverts the sign of every swapped variant without raising an error.
- Four of the five eQTL files were distributed without a tabix index. On those,
  `tabix` returns an empty result and exits successfully, which would produce
  zero tests rather than a failure; the package therefore streams the file once
  and filters to the locus windows instead of relying on the index.
- The eQTL Catalogue contains only one liver cohort (GTEx), in both release 7 and
  the release 8 pre-release, so an independent liver replication cannot be
  assembled from public processed summary statistics. The substitute is to hold
  the lead variant fixed and read its effect in every tissue
  (`code/02_colocalisation/14_fixed_variant_check.py`).
- The single-causal-variant assumption of the ABF model is checked against
  fine-mapping in GTEx v11 uterus, and the credible sets reported here are those
  implied by the ABF model. They are not SuSiE credible sets and are not
  described as such.

## Licence

MIT. See `LICENSE`.

## Citation

If you use this code or the derived tables, please cite the paper (reference to
be added on publication) and the archived package:

> Hu T, Mo Q, Chen P, Xu C, Wang Y, Sun Q, Zhu T. Replication package: shared risk
> variant at 17q21.32 in endometrial cancer and the limits of obesity-derived
> transcriptional signatures as proxies for insulin resistance. Zenodo.
> https://doi.org/10.5281/zenodo.24055114

The DOI above is the concept DOI and always resolves to the most recent version.
