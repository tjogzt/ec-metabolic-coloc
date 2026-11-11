"""
config.py
=========

Purpose
-------
Single source of truth for every path, accession and analysis parameter used
by the endometrial-cancer / insulin-resistance (EC-IR) replication package.

No analysis script may embed a path, an accession or a threshold. All of them
are resolved here, so that a single edit (or a single environment variable)
re-points the whole pipeline.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0 (replication package restructure)

Inputs
------
Environment variables (all optional, all override the documented defaults):

    ECIR_ROOT                 repository root                (default: parent of this file)
    ECIR_DATA                 directory holding input data    (default: <root>/data)
    ECIR_RESULTS              directory holding outputs       (default: <root>/results)
    ECIR_SEED                 master random seed              (default: 49)
    ECIR_N_THREADS            worker threads for tabix/expr   (default: 4)

Outputs
-------
Importable module exposing:

    PATHS        frozen dataclass of resolved directories
    PARAMS       frozen dataclass of analysis parameters
    SEED         master random seed
    COHORTS      GEO cohorts and their phenotype labels
    TISSUES      GTEx tissues and their eQTL Catalogue dataset ids
    MODULES      fixed gene modules (proliferation, adipose)
    PALETTE      figure colours, assigned by semantic role

Usage
-----
    from config.config import PATHS, PARAMS, SEED, TISSUES

Notes
-----
The only literal constant that is intentionally fixed in source is ``SEED``;
every other value is either resolved from the environment or read from a
primary data source at run time.
"""
from __future__ import annotations

import os
import pathlib
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

# --------------------------------------------------------------------------
# Repository root and derived directories
# --------------------------------------------------------------------------
_DEFAULT_ROOT = pathlib.Path(__file__).resolve().parent.parent


@dataclass(frozen=True)
class Paths:
    """Resolved filesystem locations. Created lazily by ``ensure()``."""

    root: pathlib.Path
    data: pathlib.Path
    results: pathlib.Path
    tables: pathlib.Path
    figures: pathlib.Path
    cache: pathlib.Path

    def ensure(self) -> "Paths":
        """Create the output directories if they do not exist yet."""
        for p in (self.results, self.tables, self.figures, self.cache):
            p.mkdir(parents=True, exist_ok=True)
        return self


def _resolve_paths() -> Paths:
    root = pathlib.Path(os.environ.get("ECIR_ROOT", _DEFAULT_ROOT)).resolve()
    data = pathlib.Path(os.environ.get("ECIR_DATA", root / "data")).resolve()
    results = pathlib.Path(os.environ.get("ECIR_RESULTS", root / "results")).resolve()
    return Paths(
        root=root,
        data=data,
        results=results,
        tables=results / "tables",
        figures=results / "figures",
        cache=root / ".cache",
    )


PATHS = _resolve_paths()

# --------------------------------------------------------------------------
# Master random seed: the single fixed constant in this file
# --------------------------------------------------------------------------
SEED: int = int(os.environ.get("ECIR_SEED", 49))

# --------------------------------------------------------------------------
# Analysis parameters
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Params:
    """Numerical and categorical parameters shared across the pipeline.

    Grouped by analysis stage. Every value can be overridden through the
    environment variables documented in the corresponding comment.
    """

    # ---- stage 1: state programme -------------------------------------
    program_fdr: float = 0.05            # ECIR_PROGRAM_FDR, Benjamini-Hochberg
    program_lfc_cut: float = 1.0         # ECIR_PROGRAM_LFC, |log2FC| cut
    program_min_cohorts: int = 2         # ECIR_PROGRAM_MIN_COHORTS, replica rule
    program_de_padj: float = 0.05        # ECIR_PROGRAM_DE_PADJ
    n_random_modules: int = 100          # ECIR_N_RANDOM_MODULES

    # ---- stage 2: colocalisation --------------------------------------
    gwas_p_thresh: float = 5e-8          # ECIR_GWAS_P, locus discovery
    locus_gap_bp: int = 1_000_000        # ECIR_LOCUS_GAP, lead variant gap
    region_pad_bp: int = 600_000         # ECIR_REGION_PAD, window around lead
    eqtl_min_snps: int = 10              # ECIR_EQTL_MIN_SNPS
    eqtl_maf_min: float = 0.01           # ECIR_EQTL_MAF_MIN
    eqtl_maf_max: float = 0.99           # ECIR_EQTL_MAF_MAX
    p12: float = 1e-5                    # ECIR_P12, prior on shared causal variant
    p1: float = 1e-4                     # ECIR_P1, prior trait 1 has a signal
    p2: float = 1e-4                     # ECIR_P2, prior trait 2 has a signal
    p12_sensitivity: Tuple[float, ...] = (1e-6, 5e-6, 1e-5, 5e-5, 1e-4)
    m12: float = 1e-5                    # ECIR_M12, prior on two distinct variants
    pp_h4_high: float = 0.8              # ECIR_PP_H4_HIGH, "shared" threshold
    pp_h3_high: float = 0.8              # ECIR_PP_H3_HIGH, "distinct" threshold

    # ---- stage 3: burden and incremental value ------------------------
    n_bootstrap: int = 200               # ECIR_N_BOOTSTRAP
    n_null_draws: int = 1000             # ECIR_N_NULL_DRAWS
    alpha: float = 0.05                  # ECIR_ALPHA, two-sided

    # ---- stage 4: single cell -----------------------------------------
    sc_min_cells_per_sample: int = 50    # ECIR_SC_MIN_CELLS
    sc_min_samples: int = 3              # ECIR_SC_MIN_SAMPLES

    # ---- stage 5: Mendelian randomisation -----------------------------
    mr_p_thresh: float = 5e-8            # ECIR_MR_P, instrument selection
    mr_min_f: float = 10.0               # ECIR_MR_MIN_F, weak-instrument floor
    mr_min_n: int = 100_000              # ECIR_MR_MIN_N, isolates the small-N block
    mr_clump_bp: int = 10_000_000        # ECIR_MR_CLUMP_BP, independence window
    mr_palindromic_lo: float = 0.42      # ECIR_MR_PAL_LO
    mr_palindromic_hi: float = 0.58      # ECIR_MR_PAL_HI
    mr_max_snps: int = 10_000            # ECIR_MR_MAX_SNPS, Egger guard

    # ---- stage 6: cross-cancer ----------------------------------------
    n_random_partitions: int = 100       # ECIR_N_RANDOM_PARTITIONS
    min_samples_per_cancer: int = 25     # ECIR_MIN_SAMPLES_CANCER


def _params() -> Params:
    """Build Params from the environment, falling back to the documented defaults."""
    d = Params()
    return Params(
        program_fdr=float(os.environ.get("ECIR_PROGRAM_FDR", d.program_fdr)),
        program_lfc_cut=float(os.environ.get("ECIR_PROGRAM_LFC", d.program_lfc_cut)),
        program_min_cohorts=int(os.environ.get("ECIR_PROGRAM_MIN_COHORTS", d.program_min_cohorts)),
        program_de_padj=float(os.environ.get("ECIR_PROGRAM_DE_PADJ", d.program_de_padj)),
        n_random_modules=int(os.environ.get("ECIR_N_RANDOM_MODULES", d.n_random_modules)),
        gwas_p_thresh=float(os.environ.get("ECIR_GWAS_P", d.gwas_p_thresh)),
        locus_gap_bp=int(os.environ.get("ECIR_LOCUS_GAP", d.locus_gap_bp)),
        region_pad_bp=int(os.environ.get("ECIR_REGION_PAD", d.region_pad_bp)),
        eqtl_min_snps=int(os.environ.get("ECIR_EQTL_MIN_SNPS", d.eqtl_min_snps)),
        eqtl_maf_min=float(os.environ.get("ECIR_EQTL_MAF_MIN", d.eqtl_maf_min)),
        eqtl_maf_max=float(os.environ.get("ECIR_EQTL_MAF_MAX", d.eqtl_maf_max)),
        p12=float(os.environ.get("ECIR_P12", d.p12)),
        p12_sensitivity=d.p12_sensitivity,
        m12=float(os.environ.get("ECIR_M12", d.m12)),
        pp_h4_high=float(os.environ.get("ECIR_PP_H4_HIGH", d.pp_h4_high)),
        pp_h3_high=float(os.environ.get("ECIR_PP_H3_HIGH", d.pp_h3_high)),
        n_bootstrap=int(os.environ.get("ECIR_N_BOOTSTRAP", d.n_bootstrap)),
        n_null_draws=int(os.environ.get("ECIR_N_NULL_DRAWS", d.n_null_draws)),
        alpha=float(os.environ.get("ECIR_ALPHA", d.alpha)),
        sc_min_cells_per_sample=int(os.environ.get("ECIR_SC_MIN_CELLS", d.sc_min_cells_per_sample)),
        sc_min_samples=int(os.environ.get("ECIR_SC_MIN_SAMPLES", d.sc_min_samples)),
        mr_p_thresh=float(os.environ.get("ECIR_MR_P", d.mr_p_thresh)),
        mr_min_f=float(os.environ.get("ECIR_MR_MIN_F", d.mr_min_f)),
        mr_min_n=int(os.environ.get("ECIR_MR_MIN_N", d.mr_min_n)),
        mr_clump_bp=int(os.environ.get("ECIR_MR_CLUMP_BP", d.mr_clump_bp)),
        mr_palindromic_lo=float(os.environ.get("ECIR_MR_PAL_LO", d.mr_palindromic_lo)),
        mr_palindromic_hi=float(os.environ.get("ECIR_MR_PAL_HI", d.mr_palindromic_hi)),
        mr_max_snps=int(os.environ.get("ECIR_MR_MAX_SNPS", d.mr_max_snps)),
        n_random_partitions=int(os.environ.get("ECIR_N_RANDOM_PARTITIONS", d.n_random_partitions)),
        min_samples_per_cancer=int(os.environ.get("ECIR_MIN_SAMPLES_CANCER", d.min_samples_per_cancer)),
    )


PARAMS = _params()
N_THREADS: int = int(os.environ.get("ECIR_N_THREADS", 4))

# --------------------------------------------------------------------------
# External resources
# --------------------------------------------------------------------------
# Primary outcome: meta-analysis of endometrial cancer including UK Biobank,
# Estonian Biobank, FinnGen R6 and a Japanese biobank (hg38).
GWAS_PRIMARY_FILE = os.environ.get("ECIR_GWAS_PRIMARY", "GCST90454186.h.tsv.gz")
# Secondary outcome used for the dual-dataset consistency check (hg19).
GWAS_SECONDARY_FILE = os.environ.get("ECIR_GWAS_SECONDARY", "GCST006465_GRCh37.tsv.gz")
# Outcome deliberately excluded: no genome-wide significant variant in FinnGen R12.
GWAS_FINNGEN_FILE = os.environ.get("ECIR_GWAS_FINNGEN",
                                   "finngen_R12_C3_CORPUS_UTERI_EXALLC.gz")
# eQTL Catalogue summary statistics are fetched over HTTP, never bundled.
EQTL_CATALOGUE_BASE = os.environ.get(
    "ECIR_EQTL_BASE",
    "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015")

# Lift-over chains are downloaded by stage 02/00_download_inputs.sh.
CHAIN_HG19_TO_HG38 = "hg19ToHg38.over.chain.gz"
CHAIN_HG38_TO_HG19 = "hg38ToHg19.over.chain.gz"

# --------------------------------------------------------------------------
# eQTL Catalogue datasets, by GTEx tissue label used in the manuscript
# --------------------------------------------------------------------------
# Verified against the eQTL Catalogue dataset metadata table
# (data_tables/dataset_metadata_r7.tsv); every entry is GTEx v8, quant_method = ge.
TISSUES: Dict[str, str] = {
    "Uterus": "QTD000346",                    # GTEx uterus
    "Liver": "QTD000266",                     # GTEx liver
    "Muscle_Skeletal": "QTD000281",           # GTEx muscle
    "Adipose_Visceral_Omentum": "QTD000121",  # GTEx adipose_visceral
    "Pancreas": "QTD000296",                  # GTEx pancreas
}

# --------------------------------------------------------------------------
# GEO cohorts for the state-programme derivation and the proxy-validity ladder
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Cohort:
    """One GEO series used to derive or validate the transcriptional programme.

    Attributes
    ----------
    gse        : GEO accession.
    tissue     : GTEx-style tissue label.
    contrast   : short identifier of the comparison performed.
    phenotype  : how insulin resistance was ascertained in this cohort.
    gold       : True for the hyperinsulinaemic-euglycaemic clamp cohort.
    """

    gse: str
    tissue: str
    contrast: str
    phenotype: str
    gold: bool = False


COHORTS: Tuple[Cohort, ...] = (
    Cohort("GSE15773", "Adipose_Visceral_Omentum", "IR_state_vs_control",
           "clinical grouping", gold=False),
    Cohort("GSE18732", "Muscle_Skeletal", "T2D_vs_control",
           "diagnosed type 2 diabetes", gold=False),
    Cohort("GSE13070", "Muscle_Skeletal", "clamp_IR_vs_sensitive",
           "hyperinsulinaemic-euglycaemic clamp", gold=True),
)

# Cohorts used only to derive the programme (not part of the proxy ladder).
PROGRAMME_COHORTS: Tuple[Cohort, ...] = COHORTS

# --------------------------------------------------------------------------
# Fixed gene modules, frozen before any data were inspected
# --------------------------------------------------------------------------
MODULES: Dict[str, List[str]] = {
    # Ten canonical proliferation markers, fixed a priori (stage 06).
    "proliferation": [
        "MKI67", "PCNA", "TOP2A", "CCNB1", "CDK1",
        "BUB1", "AURKA", "PLK1", "CENPE", "RRM2",
    ],
    # Eight adipocyte / lipid-handling markers, fixed a priori (stage 06).
    "adipose": [
        "ADIPOQ", "PLIN1", "LEP", "FABP4",
        "LIPE", "CIDEC", "PPARG", "UCP1",
    ],
}

# --------------------------------------------------------------------------
# Figure style: colours carry semantics, not decoration
# --------------------------------------------------------------------------
PALETTE: Dict[str, str] = {
    "uterus": "#C23531",     # cinnabar   : target organ / gold-standard phenotype
    "liver": "#3D6BA8",      # azurite    : metabolic tissue
    "muscle": "#177CB0",     # indigo     : skeletal muscle
    "adipose": "#7BA7C7",    # pale azurite: visceral adipose tissue
    "pancreas": "#9D2933",   # rouge      : pancreas
    "neutral": "#6E6E6E",    # neutral grey: non-significant / reference
}

FIGURE_FONT = os.environ.get("ECIR_FIG_FONT", "Arial")
FIGURE_MIN_PT = 8.0


def data_file(name: str) -> pathlib.Path:
    """Resolve a named input file inside the data directory.

    Parameters
    ----------
    name : file name relative to ``PATHS.data``.

    Returns
    -------
    Absolute path. Existence is *not* checked here; stage scripts assert it.
    """
    return PATHS.data / name
