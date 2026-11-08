# Shared variation at 17q21.32 in endometrial cancer and the limits of obesity-derived transcriptional signatures as proxies for insulin resistance

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.24055114-003399.svg)](https://doi.org/10.5281/zenodo.24055114)

Replication package for a preregistered, public-data-only analysis with three layers:
genetic colocalisation, bulk transcriptional validation, and single-cell projection.

The analysis asks whether obesity-derived transcriptional signatures can be used to
infer that endometrial cancer shares causal metabolic architecture with obesity and
type 2 diabetes, and whether scoring such a signature in tissue measures insulin
resistance. The answer differs by layer.

## What the analysis reports

| Layer | Result |
|---|---|
| Colocalisation | 503 gene–locus tests in five tissues; **3** exceed PP.H4 > 0.8, two at 17q21.32 and one at 15q15.2 |
| Programme validity | An obesity-derived programme (862 genes) tracks clamp-defined insulin resistance weakly: rank correlation **+0.10** (p = 0.004), 61.6% sign concordance among responsive genes |
| Incremental value | **ΔAUC = +0.0004** over a base clinical model (p = 0.15, 145 tumours, 45 events) |
| Single cell | 141,316 cells, 29,592 epithelial; programme score does not separate malignant from normal epithelium (z = −3.8, below the prespecified bar) |

The shared signal at 17q21.32 is the only positive genetic finding: at the lead eQTL
variant the allele that lowers *SKAP1* expression raises endometrial cancer risk
(liver, PP.H4 = 0.93; 2,805 variants in window; ρ = −0.66 between expression and risk
z-scores). The same allele shifts *SKAP1* expression by −0.25 in pancreas and −0.24 in
spleen, but not in adipose, uterus, ovary or blood.

## Repository layout

```
code/
  00_screening/        dataset availability checks (GEO, GTEx, eQTL Catalogue)
  01_programme/        metabolic programme definition and cohort concordance
  02_colocalisation/   ABF colocalisation, variant matching, single-locus checks
  03_bulk/             cohort assembly and incremental-value models
  04_singlecell/       single-cell scoring and size-matched null modules
  05_figures/          figure generation (Figures 1–4)
config.py, config.R    path configuration
results/
  tables/              seven derived tables reported in the paper
  figures/             Figures 1–4 as vector PDF, with PNG previews in figures/png/
data/
  README.md            every input, its accession, and how to obtain it
```

All paths resolve from `config.py` / `config.R` at the repository root. Every script
also honours the `ECIR_ROOT` environment variable, so the tree can be relocated:

```bash
export ECIR_ROOT=/path/to/ec-metabolic-coloc
python code/02_colocalisation/skap1_fixed_variant.py
Rscript code/05_figures/fig1_2_3_v3.R
```

## Reproducing the analysis

The analysis does not ship raw data; every input is public and re-downloadable.
`data/README.md` lists each accession, the file needed, and where it comes from.

R packages: `limma`, `ggplot2`, `patchwork`, `data.table`, `cairo_pdf` device.
Python packages: `numpy`, `pandas`, `scipy`, `scanpy`, `anndata`, `matplotlib`.

Scripts run from the repository root. Long-running steps (colocalisation across five
tissues, single-cell scoring) were executed as background jobs and write to
`results/`; the shipped tables are the outputs of those runs.

## Notes on the colocalisation layer

- Both the eQTL resource and the risk dataset are on GRCh38, so variants are matched
  directly by position and alleles; no liftover is applied.
- Expression reference files supplied to an earlier version of this analysis were
  found to be truncated on disk. `verify_local_coloc.py` re-establishes the three
  above-threshold results from complete local files and reproduces the remote-tabix
  numbers exactly, including the number of matched variants per test.
- The eQTL Catalogue contains only one liver cohort (GTEx), in both release 7 and the
  release 8 pre-release, so an independent liver replication cannot be assembled from
  public processed summary statistics. `skap1_fixed_variant.py` is the substitute:
  it holds the variant fixed and reads its effect in every tissue.

## Licence

MIT. See `LICENSE`.

## Citation

If you use this code or the derived tables, please cite the paper (reference to be added on
publication) and the archived package:

> Hu T, Mo Q, Chen P, Xu C, Wang Y, Sun Q, Zhu T. Replication package: shared risk variant at
> 17q21.32 in endometrial cancer and the limits of obesity-derived transcriptional signatures as
> proxies for insulin resistance. Zenodo. https://doi.org/10.5281/zenodo.24055114

The DOI above is the concept DOI and always resolves to the most recent version.
