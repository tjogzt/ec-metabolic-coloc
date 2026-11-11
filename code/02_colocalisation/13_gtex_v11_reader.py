"""
13_gtex_v11_reader.py
=====================

Purpose
-------
Read the two GTEx files this package needs, and nothing else.

    GTEx v8 significant variant-gene pairs
        one member per tissue inside GTEx_Analysis_v8_eQTL.tar. This is the
        table that defines "the genes with a cis-eQTL in this tissue", which is
        the gene set every stage of the colocalisation analysis is built on. It
        holds significant pairs only, so it can be used to select genes and to
        position them, never to run the ABF test itself.

    GTEx v11 SuSiE credible sets
        a parquet of fine-mapped credible sets. A gene whose expression signal
        resolves to a single credible set is compatible with the single-causal-
        variant assumption that the ABF model makes; a gene with several sets is
        not, and a high PP.H4 there would rest on a mis-specified model. The
        check exists to say which genes those are.

Both readers share the module because they answer one question: what does GTEx
say about expression in this window. Keeping them together also keeps the
variant-id parsing in one place, since both files identify variants as
``chr1_12345_A_G_b38``.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
tar_path      GTEx_Analysis_v8_eQTL.tar, or any archive with the same layout
tissue        GTEx tissue label, e.g. "Uterus"
susie_path    GTEx v11 SuSiE parquet
gene          Ensembl gene id, used to restrict the credible-set table

Outputs
-------
list_tissues()          sorted tissue labels present in the archive
load_tissue()           significant pairs of one tissue
split_variant()         chrom / pos / ref / alt columns derived from variant_id
lead_per_gene()         the most significant variant of each gene
tissue_genes()          sorted Ensembl ids with any significant pair
load_susie()            credible sets of one or more genes
credible_set_summary()  one row per credible set: size, top variant, top pip
is_single_causal_set()  whether a gene is compatible with a single causal variant

Notes
-----
The archive member is streamed rather than extracted: the tar holds one
multi-hundred-megabyte member per tissue and unpacking it would duplicate data
that is already on disk.

The SuSiE column names differ between GTEx releases, so they are resolved through
an alias table (``SUSIE_COLUMNS``) instead of being assumed. An unrecognised
file fails with its own header in the message; silently reading the wrong column
as ``pip`` would make the causal-set check meaningless.
"""
from __future__ import annotations

import gzip
import pathlib
import tarfile
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

# Columns of the GTEx v8 significant-pair member, in file order.
COLS = ["gene_id", "variant_id", "tss_distance", "ma_samples",
        "ma_count", "maf", "pval_nominal", "slope", "slope_se"]

NUMERIC_COLS = ("pval_nominal", "slope", "slope_se", "maf", "tss_distance",
                "ma_samples", "ma_count")

# Accepted spellings of the four SuSiE columns that matter, most specific first.
SUSIE_COLUMNS: Dict[str, Sequence[str]] = {
    "gene": ("molecular_trait_id", "gene_id", "phenotype_id", "trait_id"),
    "variant": ("variant_id", "snp", "variant", "rsid"),
    "pip": ("pip", "posterior_inclusion_probability", "susie_pip"),
    "cs": ("cs_id", "cs", "credible_set", "cs_index"),
}

# The consequence the causal-set check reports.
SINGLE_SET = "single"
MULTI_SET = "multiple"
NO_SET = "none"


def list_tissues(tar_path: "str | pathlib.Path") -> List[str]:
    """Tissue labels present in a GTEx significant-pairs archive.

    Parameters
    ----------
    tar_path : the archive.

    Returns
    -------
    Sorted labels, derived from the member file names.
    """
    with tarfile.open(tar_path, "r") as tf:
        return sorted({pathlib.Path(m.name).name.split(".")[0]
                       for m in tf.getmembers()
                       if "signif_variant_gene_pairs" in m.name})


def _find_member(tar: tarfile.TarFile, tissue: str):
    """Locate the significant-pair member for one tissue, or None."""
    for member in tar.getmembers():
        if "signif_variant_gene_pairs" not in member.name:
            continue
        if pathlib.Path(member.name).name.startswith(tissue + "."):
            return member
    return None


def load_tissue(tar_path: "str | pathlib.Path", tissue: str) -> pd.DataFrame:
    """Read the significant variant-gene pairs of one tissue.

    Parameters
    ----------
    tar_path : GTEx significant-pairs archive.
    tissue : GTEx tissue label, e.g. "Uterus".

    Returns
    -------
    DataFrame with lower-case column names; numeric fields are coerced, so a
    malformed line becomes NaN rather than aborting the read.

    Raises
    ------
    FileNotFoundError
        When the tissue is not in the archive. Returning an empty frame would
        look like "no eQTL in this tissue", which is a different statement.
    """
    with tarfile.open(tar_path, "r") as tf:
        member = _find_member(tf, tissue)
        if member is None:
            raise FileNotFoundError(f"tissue {tissue} not found in {tar_path}")
        stream = tf.extractfile(member)
        df = pd.read_csv(stream, sep="\t", compression="gzip", dtype=str)

    df.columns = [c.lower() for c in df.columns]
    for col in NUMERIC_COLS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def split_variant(df: pd.DataFrame) -> pd.DataFrame:
    """Split ``chr1_12345_A_G_b38`` into usable columns.

    Parameters
    ----------
    df : frame with a ``variant_id`` column.

    Returns
    -------
    Copy of the frame with chrom, pos, ref, alt added. ``pos`` is nullable
    because a malformed identifier must not become position 0.

    Notes
    -----
    The build suffix is dropped: only the coordinates are kept, and the caller
    decides which build they belong to. Every GTEx v8 member is GRCh38, which is
    also the build of the eQTL Catalogue summary statistics.
    """
    parts = df["variant_id"].str.split("_", expand=True, regex=False)
    out = df.copy()
    out["chrom"] = parts[0].str.replace("chr", "", regex=False)
    out["pos"] = pd.to_numeric(parts[1], errors="coerce").astype("Int64")
    out["ref"] = parts[2]
    out["alt"] = parts[3]
    # Drop the Ensembl version suffix, so that a gene is one key across tables.
    # Every downstream comparison joins on gene id, and a versioned id from the
    # archive would never equal the unversioned id in an eQTL Catalogue file.
    if "gene_id" in out.columns:
        out["gene_id"] = out["gene_id"].astype(str).str.split(".", regex=False).str[0]
    return out


def lead_per_gene(df: pd.DataFrame) -> pd.DataFrame:
    """Most significant variant of each gene.

    Parameters
    ----------
    df : significant pairs, after :func:`split_variant`.

    Returns
    -------
    One row per gene, in the input column order.
    """
    idx = df.groupby("gene_id")["pval_nominal"].idxmin()
    return df.loc[idx].reset_index(drop=True)


def tissue_genes(df: pd.DataFrame) -> List[str]:
    """Sorted Ensembl ids with at least one significant pair.

    Parameters
    ----------
    df : output of :func:`load_tissue`.

    Returns
    -------
    Gene ids with the Ensembl version suffix removed, de-duplicated.
    """
    return sorted({str(g).split(".")[0] for g in df["gene_id"].unique()})


def variant_chrom_pos(variant_id: str) -> Tuple[str, int]:
    """Chromosome and 1-based position from a GTEx variant identifier.

    Parameters
    ----------
    variant_id : identifier of the form ``chr17_48127316_G_A_b38``.

    Returns
    -------
    (chromosome without the "chr" prefix, position).

    Raises
    ------
    ValueError
        When the identifier has fewer than two fields, which is the only way an
        identifier can fail to describe a position.

    Notes
    -----
    The number of trailing fields is not fixed: GTEx writes a build suffix
    (``_b38``) while other releases stop after the two alleles. Splitting and
    taking the first two fields accepts both; unpacking a fixed number of fields
    would fail on the form GTEx actually ships.
    """
    parts = str(variant_id).split("_")
    if len(parts) < 2:
        raise ValueError(f"not a variant identifier: {variant_id!r}")
    return parts[0].replace("chr", ""), int(parts[1])


def genes_near_loci(tar_path, tissue: str, loci: pd.DataFrame,
                   pad: int) -> Dict[int, List[str]]:
    """Genes with a significant pair inside a padded window around each locus.

    Parameters
    ----------
    tar_path : GTEx significant-pairs archive.
    tissue : GTEx tissue label.
    loci : DataFrame with columns chrom, start, end.
    pad : window half-width in base pairs.

    Returns
    -------
    dict mapping the locus row index to a sorted list of unique Ensembl ids.

    Raises
    ------
    SystemExit
        When the tissue is absent from the archive. An empty mapping would read
        as "no gene near any locus", which is a different claim.

    Notes
    -----
    Columns are located by name from the member's own header, not by position,
    so a release that reorders variant_id and gene_id is still read correctly.
    """
    windows = [(i, str(row["chrom"]).replace("chr", ""),
                int(row["start"]) - pad, int(row["end"]) + pad)
               for i, row in loci.iterrows()]
    out: Dict[int, List[str]] = {i: [] for i in range(len(loci))}

    with tarfile.open(tar_path, "r") as tar:
        member = _find_member(tar, tissue)
        if member is None:
            raise SystemExit(f"tissue {tissue} not found in "
                             f"{pathlib.Path(tar_path).name}")
        stream = gzip.GzipFile(fileobj=tar.extractfile(member))

        # Positional defaults; the header of this member decides in practice.
        i_var, i_gene = 0, 1
        first = True
        for line in stream:
            fields = line.decode().split("\t")
            if first:                             # first line is the header
                header = [f.lower() for f in fields]
                i_var, i_gene = header.index("variant_id"), header.index("gene_id")
                first = False
                continue
            try:
                chrom, pos = variant_chrom_pos(fields[i_var])
            except (ValueError, IndexError):
                continue
            for i, win_chrom, lo, hi in windows:
                if chrom == win_chrom and lo <= pos <= hi:
                    out[i].append(fields[i_gene].split(".")[0])

    return {k: sorted(set(v)) for k, v in out.items()}


def _resolve_columns(df: pd.DataFrame, path) -> Dict[str, str]:
    """Map the internal column names onto the file's own headers.

    Parameters
    ----------
    df : the loaded SuSiE frame, used for its header.
    path : the source, quoted in the error message.

    Returns
    -------
    dict internal name -> actual column name.

    Raises
    ------
    KeyError
        If a required field has none of its known spellings. The message lists
        the header so that the alias table can be extended in one place.
    """
    resolved: Dict[str, str] = {}
    for internal, aliases in SUSIE_COLUMNS.items():
        for alias in aliases:
            if alias in df.columns:
                resolved[internal] = alias
                break
        else:
            raise KeyError(f"{pathlib.Path(path).name}: no column for '{internal}' "
                           f"(tried {list(aliases)}); header is {list(df.columns)}")
    return resolved


def load_susie(susie_path, gene: Optional[str] = None,
               tissue: Optional[str] = None) -> pd.DataFrame:
    """Read a GTEx SuSiE credible-set table.

    Parameters
    ----------
    susie_path : parquet written by the GTEx fine-mapping pipeline.
    gene : optional Ensembl id; when given, only that gene is returned.
    tissue : optional tissue label; used only when the file carries such a column.

    Returns
    -------
    DataFrame with the standard names gene, variant, pip, cs, and the columns the
    file already had.

    Notes
    -----
    parquet is read through pandas so that the optional pyarrow dependency stays
    optional: importing this module must not require it.
    """
    df = pd.read_parquet(susie_path)
    cols = _resolve_columns(df, susie_path)
    df = df.rename(columns={v: k for k, v in cols.items()})

    df["pip"] = pd.to_numeric(df["pip"], errors="coerce")
    df["cs"] = df["cs"].astype(str)

    if tissue is not None:
        for name in ("tissue", "tissue_id", "dataset_id"):
            if name in df.columns:
                df = df[df[name].astype(str) == str(tissue)]
                break

    if gene is not None:
        key = df["gene"].astype(str).str.split(".", regex=False).str[0]
        df = df[key == str(gene).split(".")[0]]

    return df.reset_index(drop=True)


def credible_set_summary(df: pd.DataFrame, gene: str,
                         min_pip: float = 0.5) -> pd.DataFrame:
    """Summarise the credible sets of one gene.

    Parameters
    ----------
    df : output of :func:`load_susie`, already restricted to the gene.
    gene : Ensembl id, used only for the error message.
    min_pip : a variant counts as a set member at or above this posterior
        inclusion probability.

    Returns
    -------
    DataFrame with cs, n_variants, n_members, top_variant, top_pip, one row per
    credible set, largest first.

    Raises
    ------
    SystemExit
        When the gene has no credible set at all, which means it was not fine
        mapped rather than that it has no signal.
    """
    sub = df[df["gene"].astype(str).str.split(".", regex=False).str[0]
             == str(gene).split(".")[0]]
    if sub.empty:
        raise SystemExit(f"{gene} has no credible set in the SuSiE table; "
                         f"it was not fine mapped, so the check cannot run")

    rows = []
    for cs, group in sub.groupby("cs"):
        top = group.loc[group["pip"].idxmax()]
        rows.append({
            "cs": cs,
            "n_variants": len(group),
            "n_members": int((group["pip"] >= min_pip).sum()),
            "top_variant": str(top["variant"]),
            "top_pip": float(top["pip"]),
        })
    return (pd.DataFrame(rows)
            .sort_values("n_members", ascending=False)
            .reset_index(drop=True))


def is_single_causal_set(df: pd.DataFrame, gene: str,
                         min_pip: float = 0.5) -> str:
    """Classify a gene by how many credible sets it carries.

    Parameters
    ----------
    df : output of :func:`load_susie`, restricted or not.
    gene : Ensembl id.
    min_pip : membership threshold passed to :func:`credible_set_summary`.

    Returns
    -------
    "single", "multiple", or "none" (see the module constants).

    Notes
    -----
    Only the first case licenses a single-causal-variant reading of PP.H4. This
    is a statement about the expression side; the outcome side is not fine
    mapped here, so a "single" verdict is necessary, not sufficient.
    """
    summary = credible_set_summary(df, gene, min_pip=min_pip)
    n = len(summary[summary["n_members"] > 0])
    if n == 0:
        return NO_SET
    return SINGLE_SET if n == 1 else MULTI_SET
