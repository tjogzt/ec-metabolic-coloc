"""
06_load_programme.py
====================

Purpose
-------
Load the endometrial single-cell atlas, quality-control it, annotate lineages,
and project all three locked scores onto every cell: the insulin-response
programme (up and down arms separately), the disease signature, and the
cell-cycle panel.

This is the step that turns abstract gene lists into per-cell scores. It is
deliberately downstream of 05: the programme is read from the tables directory
as a frozen artefact and never recomputed here, so a change in the scoring code
cannot silently redefine what is being scored.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
PATHS.data/<GSE>/...           10x matrices or h5 files, one directory per series
PATHS.tables/IR_program_v1.tsv                  the programme, with direction
PATHS.tables/EC_signature_v2_prolifadj_compact.tsv   the disease signature

Series that are absent from the data directory are skipped, so the atlas can be
rebuilt from a partial download and the log records exactly which series
entered.

Outputs
-------
PATHS.cache/atlas_core_qc.h5ad       after filtering, before normalisation
PATHS.cache/atlas_core_scored.h5ad   final object: leiden, lineage, three scores
PATHS.tables/state_program_table.tsv dataset x lineage means of the three scores

Usage
-----
    python code/03_burden/06_load_programme.py

Notes
-----
Scoring uses scanpy's ``score_genes``, which subtracts an expression-matched
control set. The control size is set to the size of the gene set itself so that
up and down arms -- which differ in length -- are scored on comparable
control sets; leaving it at the default 50 would make a 608-gene arm and a
254-gene arm incomparable.

The lineage annotation is a maximum-over-markers rule on per-cell marker scores,
with the runner-up score recorded alongside so that ambiguity can be reported
rather than hidden.
"""
from __future__ import annotations

import pathlib
import re
import sys
from importlib import import_module
from typing import Dict, List, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread

# Repository root for `config`; this directory for the sibling stages.
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))

from config.config import PATHS

# The cell-cycle panel is defined once, in the stage that audits the signature
# against it; importing it keeps "proliferation" a single definition.
CYCLING_GENES = import_module("05_build_programme").CYCLING_GENES

# --------------------------------------------------------------------------
# Module constants
# --------------------------------------------------------------------------
# Series that make up the atlas. Kept as a tuple, not read from config: the
# single-cell layer has its own cohort list, and config.COHORTS describes the
# bulk cohorts that define the programme.
ATLAS_SERIES = ("GSE278879", "GSE268903", "GSE212217", "GSE251923")

# Lineage markers, fixed before annotation was inspected.
MARKERS: Dict[str, List[str]] = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],
    "Fibroblast": ["DCN", "LUM", "COL1A1", "COL1A2", "PDGFRA"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "PLVAP"],
    "Pericyte/SMC": ["RGS5", "ACTA2", "DES", "MYH11", "TAGLN"],
    "T/NK": ["PTPRC", "CD3D", "CD3E", "CD2", "NKG7", "GNLY"],
    "Myeloid": ["CD68", "CD14", "LYZ", "C1QA", "ITGAX"],
    "B/Plasma": ["MS4A1", "CD79A", "MZB1", "JCHAIN"],
    "Mast": ["TPSAB1", "CPA3", "MS4A2"],
}

# QC thresholds.
MIN_GENES_PER_CELL = 200     # below this a library is empty, not sparse
MIN_CELLS_PER_GENE = 3       # below this a gene is not measured, just absent
MAX_PCT_MITO = 20.0          # percent of counts in MT- genes

# Normalisation and embedding settings.
TARGET_SUM = 1e4             # counts per cell after rescaling
N_TOP_GENES = 2000           # highly variable genes used for the embedding
N_PCS = 30
N_NEIGHBORS = 15
LEIDEN_RESOLUTION = 1.0
SCALE_MAX = 10               # clip z-scores so that a few cells cannot dominate

# Verbosity of scanpy's own logging; 1 keeps warnings and drops progress bars.
SC_VERBOSITY = 1

sc.settings.verbosity = SC_VERBOSITY


# --------------------------------------------------------------------------
# 10x discovery
# --------------------------------------------------------------------------
def find_10x(sample_dir: pathlib.Path) -> Optional[dict]:
    """Locate a 10x matrix inside one sample directory.

    Parameters
    ----------
    sample_dir : directory to inspect.

    Returns
    -------
    dict describing either an h5 file (``{"h5": path}``) or the matrix triple
    (``{"mtx", "barcodes", "features"}``), or None if the directory holds
    neither.

    Notes
    -----
    Only this level is inspected, never recursively. A recursive search lets a
    parent directory collect its children's files, and the "first matrix in
    alphabetical order plus same-prefix barcodes" rule then pairs a matrix with
    another sample's barcodes; the resulting per-sample de-duplication on the
    matrix path subsequently discards a real sample as a duplicate.
    """
    files = [p for p in sample_dir.iterdir() if p.is_file()]
    if not files:
        return None

    h5 = [p for p in files if p.suffix == ".h5"]
    if h5:
        return {"h5": str(sorted(h5)[0])}

    mtx = [p for p in files if re.search(r"matrix\.mtx(\.gz)?$", p.name)]
    if not mtx:
        return None
    mtx = sorted(mtx)[0]
    prefix = re.sub(r"matrix\.mtx(\.gz)?$", "", mtx.name)

    # Barcodes and features must share the matrix's prefix, otherwise a
    # directory holding two samples would be assembled from both.
    names = {p.name for p in files}
    bc = [n for n in (f"{prefix}barcodes.tsv.gz", f"{prefix}barcodes.tsv") if n in names]
    ft = [n for n in (f"{prefix}features.tsv.gz", f"{prefix}genes.tsv.gz",
                      f"{prefix}features.tsv", f"{prefix}genes.tsv") if n in names]
    if not bc or not ft:
        return None
    return {"mtx": str(mtx),
            "barcodes": str(mtx.parent / bc[0]),
            "features": str(mtx.parent / ft[0])}


def read_10x(spec: dict) -> ad.AnnData:
    """Read one 10x sample into an AnnData object.

    Parameters
    ----------
    spec : output of :func:`find_10x`.

    Returns
    -------
    AnnData with cells as rows and gene symbols as variables.

    Notes
    -----
    The Matrix Market file is genes x cells, hence the transpose. Counts are
    kept as float32 so that the concatenated atlas stays within memory.
    """
    if "h5" in spec:
        return sc.read_10x_h5(spec["h5"], gex_only=True)

    feat = pd.read_csv(spec["features"], sep="\t", header=None)
    symbols = feat[1].astype(str)
    X = mmread(spec["mtx"]).T.tocsr().astype(np.float32)
    barcodes = pd.read_csv(spec["barcodes"], sep="\t", header=None)[0].astype(str)
    a = ad.AnnData(X=X, obs=pd.DataFrame(index=barcodes),
                   var=pd.DataFrame(index=symbols.values))
    a.var_names_make_unique()
    return a


def load_dataset(gse: str) -> List[ad.AnnData]:
    """Load every sample directory of one series.

    Parameters
    ----------
    gse : series accession, i.e. a directory name under ``PATHS.data``.

    Returns
    -------
    list of AnnData objects, one per sample; empty if the series is absent.

    Notes
    -----
    Sample-level de-duplication is on the matrix path. Without it a parent
    directory and its subdirectory both resolve to the same matrix and the
    sample is loaded twice, which doubles its weight in every downstream mean.
    """
    root = PATHS.data / gse
    if not root.exists():
        return []

    # Every directory at or below the series root is a candidate sample.
    cands = []
    for d in [root] + [p for p in root.rglob("*") if p.is_dir()]:
        spec = find_10x(d)
        if spec:
            cands.append((d, spec))

    seen, uniq = set(), []
    for d, spec in cands:
        key = spec.get("mtx") or spec.get("h5")
        if key in seen:
            continue
        seen.add(key)
        uniq.append((d, spec))

    out: List[ad.AnnData] = []
    for d, spec in uniq:
        try:
            a = read_10x(spec)
        except Exception as e:                     # noqa: BLE001
            # One unreadable sample must not lose the rest of the series.
            print(f"  [skip] {gse}/{d.name}: {e}", file=sys.stderr)
            continue
        a.obs["dataset"] = gse
        a.obs["sample"] = d.name if d != root else gse
        a.obs_names = [f"{gse}_{d.name}_{b}" for b in a.obs_names]
        out.append(a)
        print(f"  [{gse}] {d.name}: {a.n_obs} cells x {a.n_vars} genes")
    return out


# --------------------------------------------------------------------------
# Pipeline
# --------------------------------------------------------------------------
def main() -> None:
    """Build the atlas, score it, and write the object and the summary table."""
    PATHS.ensure()

    series = [g for g in ATLAS_SERIES if (PATHS.data / g).exists()]
    print("loading:", series)
    parts: List[ad.AnnData] = []
    for g in series:
        parts += load_dataset(g)
    if not parts:
        print("no expression matrix found under the data directory", file=sys.stderr)
        sys.exit(1)

    print(f"\nconcatenating {len(parts)} samples ...")
    a = ad.concat(parts, join="outer", label="batch_id", index_unique=None,
                  fill_value=0)
    print(f"merged: {a.n_obs:,} cells x {a.n_vars:,} genes")
    a.var_names = a.var_names.str.upper()

    # ---- QC ------------------------------------------------------------
    a.var["mt"] = a.var_names.str.startswith(("MT-",))
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None)
    n0 = a.n_obs
    sc.pp.filter_cells(a, min_genes=MIN_GENES_PER_CELL)
    sc.pp.filter_genes(a, min_cells=MIN_CELLS_PER_GENE)
    a = a[a.obs.pct_counts_mt < MAX_PCT_MITO].copy()
    print(f"QC: {n0:,} -> {a.n_obs:,} cells ({100 * a.n_obs / n0:.1f}% kept), "
          f"median genes per cell {a.obs.n_genes_by_counts.median():.0f}")
    a.write(PATHS.cache / "atlas_core_qc.h5ad")

    # ---- normalisation, embedding, clustering ---------------------------
    sc.pp.normalize_total(a, target_sum=TARGET_SUM)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=N_TOP_GENES, batch_key="dataset")
    b = a[:, a.var.highly_variable].copy()
    sc.pp.scale(b, max_value=SCALE_MAX)
    sc.tl.pca(b, n_comps=N_PCS, svd_solver="arpack")

    # Harmony is optional: it needs harmonypy, which is not a hard dependency of
    # the package. Falling back to uncorrected PCA changes the batch structure
    # of the embedding, not the scores, which are computed on the full matrix.
    # On anndata >= 0.12 scanpy's harmony wrapper raises while writing its result
    # back; the fallback below is the expected outcome there, not an anomaly.
    try:
        import harmonypy                                    # noqa: F401
        sc.external.pp.harmony_integrate(b, "dataset")
        rep = "X_pca_harmony"
        print("integrated across datasets with harmony")
    except Exception as e:                                  # noqa: BLE001
        rep = "X_pca"
        print(f"harmony unavailable ({e}); using uncorrected PCA")

    sc.pp.neighbors(b, n_neighbors=N_NEIGHBORS, use_rep=rep)
    sc.tl.leiden(b, resolution=LEIDEN_RESOLUTION, key_added="leiden")

    # Copy the embedding back onto the full object; the clustered object only
    # holds highly variable genes and cannot be scored.
    a.obs["leiden"] = b.obs["leiden"].values
    a.obsm["X_pca"] = b.obsm["X_pca"]
    if rep in b.obsm:
        a.obsm[rep] = b.obsm[rep]
    sc.tl.umap(b)
    a.obsm["X_umap"] = b.obsm["X_umap"]
    print(f"clustering: {a.obs.leiden.nunique()} clusters")

    # ---- lineage annotation --------------------------------------------
    score = pd.DataFrame(index=a.obs_names)
    for lineage, markers in MARKERS.items():
        present = [g for g in markers if g in a.var_names]
        if not present:
            continue
        # anndata rejects "/" in obs names, so marker-set column names are
        # sanitised; Pericyte/SMC becomes Pericyte_SMC everywhere.
        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", lineage)
        sc.tl.score_genes(a, present, score_name=f"lin_{safe}",
                          use_raw=False if a.raw is None else True)
        score[lineage] = a.obs[f"lin_{safe}"].values

    # No marker set at all means the variable names are not gene symbols, which
    # would also invalidate every other score below; say so rather than letting
    # the row-wise maximum fail on an empty frame.
    if score.shape[1] == 0:
        raise SystemExit("no lineage marker found among the variable names; "
                         "check that the matrices are annotated with gene symbols")

    a.obs["lineage"] = [re.sub(r"[^A-Za-z0-9_.-]", "_", str(x))
                        for x in score.idxmax(axis=1).values]
    # Confidence: margin between the best and the second-best lineage score.
    a.obs["lineage_conf"] = (score.max(axis=1) - score.apply(
        lambda r: r.nlargest(2).iloc[-1], axis=1)).values
    print("\nlineage distribution:")
    print(a.obs.lineage.value_counts().to_string())

    # ---- programme, signature and proliferation scores ------------------
    prog = pd.read_csv(PATHS.tables / "IR_program_v1.tsv", sep="\t", index_col=0)
    sig = pd.read_csv(PATHS.tables / "EC_signature_v2_prolifadj_compact.tsv",
                      sep="\t")

    up = [g for g in prog.index[prog.direction == "UP"] if g in a.var_names]
    dn = [g for g in prog.index[prog.direction == "DOWN"] if g in a.var_names]
    print(f"\nprogramme coverage: UP {len(up)}/{int((prog.direction == 'UP').sum())}, "
          f"DOWN {len(dn)}/{int((prog.direction == 'DOWN').sum())}")
    sc.tl.score_genes(a, up, score_name="IR_UP", ctrl_size=len(up))
    sc.tl.score_genes(a, dn, score_name="IR_DOWN", ctrl_size=len(dn))
    # Net score: the arms are meant to move in opposite directions, so a single
    # difference captures the programme rather than general transcriptional
    # activity.
    a.obs["IR_score"] = a.obs.IR_UP - a.obs.IR_DOWN

    sg_up = [g for g in sig.gene[sig.direction == "UP"] if g in a.var_names]
    sg_dn = [g for g in sig.gene[sig.direction == "DOWN"] if g in a.var_names]
    print(f"signature coverage: UP {len(sg_up)}/{int((sig.direction == 'UP').sum())}, "
          f"DOWN {len(sg_dn)}/{int((sig.direction == 'DOWN').sum())}")
    sc.tl.score_genes(a, sg_up, score_name="ECsig_UP", ctrl_size=len(sg_up))
    sc.tl.score_genes(a, sg_dn, score_name="ECsig_DOWN", ctrl_size=len(sg_dn))
    a.obs["ECsig_score"] = a.obs.ECsig_UP - a.obs.ECsig_DOWN

    cyc = [g for g in CYCLING_GENES if g in a.var_names]
    sc.tl.score_genes(a, cyc, score_name="Prolif", ctrl_size=len(cyc))
    print(f"proliferation score from {len(cyc)} genes")

    a.write(PATHS.cache / "atlas_core_scored.h5ad")
    print(f"\nwritten {PATHS.cache / 'atlas_core_scored.h5ad'}")

    # ---- dataset x lineage summary --------------------------------------
    s = (a.obs.groupby(["dataset", "lineage"], observed=True)
         .agg(n=("IR_score", "size"),
              IR=("IR_score", "mean"),
              ECsig=("ECsig_score", "mean"),
              Prolif=("Prolif", "mean"))
         .round(3).reset_index())
    s.to_csv(PATHS.tables / "state_program_table.tsv", sep="\t", index=False)
    print("\nstate-programme table:")
    print(s.to_string(index=False))


if __name__ == "__main__":
    main()
