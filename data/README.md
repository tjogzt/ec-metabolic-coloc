# Data sources

No raw data are redistributed here. Every input is public and re-downloadable; this
file records what each layer needs, where it comes from, and the local filename the
scripts expect under `$ECIR_ROOT/data/`.

## Genetic layer

| Input | Accession / route | Local name |
|---|---|---|
| Endometrial cancer risk, analysed set | GWAS Catalog **GCST90454186**, harmonised hg38 build: `https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90454001-GCST90455000/GCST90454186/harmonised/` | `GCST90454186.h.tsv.gz` (373 MB) |
| Endometrial cancer risk, preregistered set | GWAS Catalog **GCST006465**, GRCh37 | `GCST006465_GRCh37.tsv.gz` |
| cis-eQTL, five tissues | eQTL Catalogue **QTS000015** (GTEx v8), gene-level quantification: liver `QTD000266`, skeletal muscle `QTD000281`, pancreas `QTD000296`, uterus `QTD000346`, visceral adipose `QTD000121`. Remote tabix: `https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/<QTD>/<QTD>.all.tsv.gz` | `<QTD>.all.tsv.gz` |
| Cross-tissue fixed-variant check | Further GTEx v8 gene-level datasets from QTS000015: adipose `QTD000116`, spleen `QTD000326`, blood `QTD000356`, LCL `QTD000221`, ovary `QTD000291`, vagina `QTD000351`, breast `QTD000211`, thyroid `QTD000341`, kidney cortex `QTD000261` | queried remotely, not stored |
| Fine-mapping (single-causal-set check) | GTEx v11 SuSiE: `https://storage.googleapis.com/adult-gtex/bulk-qtl/v11/susie/GTEx_Analysis_v11_eQTL_SuSiE.tar` | `data/GTEx_Analysis_v11_eQTL_SuSiE/` (Uterus parquet) |
| Instrument checks (not used for the main estimate) | FinnGen R12 endpoints `BMI_IRN`, `T2D`, `C3_CORPUS_UTERI_EXALLC`, `N14_ENDOMETRIOSIS` | `finngen_R12_*.gz` |
| Coordinate conversion | UCSC chain files `hg19ToHg38.over.chain.gz`, `hg38ToHg19.over.chain.gz` | same |

The dataset-to-tissue mapping for eQTL Catalogue is published at
`https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/data_tables/dataset_metadata_r7.tsv`.
Note that the API (`/eqtl/api/v2/datasets`) has been retired; the metadata table above
is the current route.

## Bulk layer

| Input | Accession | Local name |
|---|---|---|
| Endometrioid carcinomas used for the incremental-value model (145 tumours, 45 events) | GEO **GSE120490** | `GSE120490_series_matrix.txt.gz` |
| Adipose, insulin resistant vs insulin sensitive | GEO **GSE15773** | `GSE15773_series_matrix.txt.gz` |
| Skeletal muscle, type 2 diabetes vs normal glucose tolerance | GEO **GSE18732** | `GSE18732_series_matrix.txt.gz` |
| Skeletal muscle, clamp-defined insulin resistant vs sensitive | GEO **GSE13070** | `GSE13070_series_matrix.txt.gz` |
| Platform annotation | GEO **GPL570**, **GPL6244** | `GPL570_annot.tsv`, `GPL6244.annot.gz` |

GSE18732 is annotated on a retired transcript set: 31% of its probes cannot be mapped,
so its gene universe is roughly half that of the other two cohorts. The concordance
estimate for that cohort is therefore less comparable than the others, and the paper
says so.

## Single-cell layer

| Input | Accession | Local name |
|---|---|---|
| Endometrial single-cell transcriptomes, 11 samples, 141,316 cells | GEO **GSE278879** | `GSE278879/GSM*/` |

## Programme definition

The 862-gene programme is shipped as `results/tables/metabolic_program_v2.tsv`
(608 up, 254 down) with a SHA-256 recorded alongside the manuscript. It was derived
without outcome labels; the derivation is in `code/01_programme/`.

## Access notes

- `GCST90454186` harmonised files are GRCh38; no liftover is required for the
  colocalisation layer.
- eQTL Catalogue remote tabix works over HTTPS but individual dataset connections
  occasionally stall. Batch queries in `code/02_colocalisation/` use per-dataset
  timeouts and three retries.
- Two of the three genes above the colocalisation threshold carry GENCODE v26
  identifiers (`RP11-521C20.5`, `RP5-890E16.5`) that current Ensembl and GENCODE
  releases no longer list. This is a property of the annotation, not an error in the
  mapping.
