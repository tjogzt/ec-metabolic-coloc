"""
10_positive_control.py
======================

Purpose
-------
Two checks that answer the same question from opposite sides: can this
colocalisation code return a positive at all?

    A. positive control   run the identical ABF test on gene-trait pairs that
                          are known to colocalise, and show that it recovers
                          them. If a canonical type 2 diabetes gene in pancreas
                          also fails to reach a high PP.H4, then the pipeline is
                          broken and a negative result for endometrial cancer
                          cannot be believed.
    B. power curve        map eQTL strength onto the probability of reaching a
                          high PP.H4, holding a shared causal variant fixed by
                          construction. This is evidence about whether the
                          available tissue can answer the question, which is
                          different from evidence about whether a shared causal
                          variant exists.

Both use real variant sets (real positions, alleles and standard errors) and
only simulate effect sizes, so neither depends on the shape of a synthetic LD
structure.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
--section     control, power, or both (default: both)
--tissue      GTEx tissue for the control (default: Pancreas)
--gwas        outcome file in the data directory for the control
              (default: finngen_R12_T2D.gz, the FinnGen R12 type 2 diabetes
              endpoint); the file layout is detected by 02_gwas_reader
--set         canonical (default) or bmi: which gene panel to run
--index       pre-built gene window table; built from the local eQTL file and
              cached under PATHS.cache when omitted
--trials      Monte-Carlo repetitions per (eQTL z, outcome z) cell (default: 60)
--seed        RNG seed for the simulation (default: config.SEED)

The control needs the all-associations eQTL file locally, because the gene
window table is derived from it. This is the same file 07 and 15 read.

Outputs
-------
results/tables/positive_control_t2d.tsv   one row per gene, columns:
    PP.H0 .. PP.H4, n_snps, top_snp_share, gene, gene_id, tissue, best_gwas_p

results/tables/coloc_power.tsv            one row per (eQTL z, outcome z) cell:
    eQTL_z, GWAS_z, H4_median, H4_p25, H4_p75, detect_rate_H4_gt_0.8, H3_median

Notes
-----
The power table used Chinese column labels in its first release. They are
translated here because the replication package is ASCII-only; the estimands
are unchanged.

The earlier version drew two independent batches of 60 trials, one for PP.H4 and
one for PP.H3. H3 and H4 are now read from the same draw: the estimand is
unchanged and the cost halves. The seed is taken from config.SEED rather than
being written into the script, so the whole package shares one seed; the first
release used 7, which only shifts the reported quantiles by Monte-Carlo error.

The shared causal variant is placed at the true eQTL lead of the window, not at
a random variant, because a shared causal variant that is not the expression
lead would not be a shared causal variant in any useful sense. The outcome
standard error is drawn as a multiple of the eQTL standard error in [0.7, 1.3],
which reproduces the range seen between the two real datasets.
"""
from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys
import time
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

# Local package imports; the directory holding this file is on sys.path.
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))       # repository root, for config/

from config.config import PATHS, PARAMS, TISSUES, SEED, EQTL_CATALOGUE_BASE
from importlib import import_module

_abf = import_module("05_coloc_abf_engine")
_gwas = import_module("02_gwas_reader")
_eqtl = import_module("03_eqtl_catalogue_reader")

# Canonical type 2 diabetes genes, fixed a priori from the literature.
CANONICAL: Dict[str, str] = {
    "SLC30A8": "ENSG00000164756",
    "KCNJ11": "ENSG00000187486",
    "ABCC8": "ENSG00000006071",
    "GCK": "ENSG00000106633",
    "HNF1A": "ENSG00000135100",
    "HNF1B": "ENSG00000275410",
    "TCF7L2": "ENSG00000148737",
    "PPARG": "ENSG00000132170",
    "KCNQ1": "ENSG00000053918",
    "IGF2BP2": "ENSG00000073756",
    "MTNR1B": "ENSG00000134640",
    "CDKN2B": "ENSG00000147883",
}

# Body-mass-index loci. Most act through regulatory variants, which makes them a
# fair test of an expression-mediated pipeline.
BMI_GENES: Dict[str, str] = {
    "FTO": "ENSG00000140718",
    "MC4R": "ENSG00000166603",
    "SH2B1": "ENSG00000178188",
    "TMEM18": "ENSG00000151353",
    "BDNF": "ENSG00000176697",
    "PCSK1": "ENSG00000175426",
    "GCKR": "ENSG00000084734",
    "PPARG": "ENSG00000132170",
    "LEPR": "ENSG00000116678",
    "IRX3": "ENSG00000166444",
    "IRX5": "ENSG00000176842",
}

# Power simulation region: the SLC30A8 window in GTEx uterus. The coordinates are
# the cis window the eQTL Catalogue tested, not a hand-picked peak.
POWER_CHROM = "8"
POWER_LO = 38_600_000
POWER_HI = 39_700_000
POWER_GENE = "ENSG00000197140"

POWER_COLS = ["eQTL_z", "GWAS_z", "H4_median", "H4_p25", "H4_p75",
              "detect_rate_H4_gt_0.8", "H3_median"]


def resolve_input(name: str) -> pathlib.Path:
    """Resolve a user-supplied file name against the data directory."""
    path = pathlib.Path(name)
    return path if path.is_absolute() else PATHS.data / path


def remote_size(url: str) -> int:
    """Content-Length of a remote file, or 0 when it cannot be read.

    Parameters
    ----------
    url : HTTP(S) location of the file.

    Returns
    -------
    Size in bytes, or 0 when the request failed. Zero is not treated as an
    error: the run may legitimately be offline, in which case the local file
    cannot be validated but can still be used.
    """
    proc = subprocess.run(["curl", "-sI", "-m", "60", url],
                          capture_output=True, text=True)
    for line in proc.stdout.splitlines():
        if line.lower().startswith("content-length"):
            try:
                return int(line.split(":", 1)[1].strip())
            except ValueError:
                return 0
    return 0


def local_eqtl(qtd: str) -> pathlib.Path:
    """Locate and validate the all-associations file for one dataset.

    Parameters
    ----------
    qtd : eQTL Catalogue dataset id.

    Returns
    -------
    Path to a local file whose size matches the remote object.

    Raises
    ------
    SystemExit
        When the file is absent, or present but truncated. A truncated bgzip
        makes tabix return zero rows and exit successfully, which would look
        like "no eQTL in this tissue" rather than "broken download".
    """
    local = PATHS.data / f"{qtd}.all.tsv.gz"
    if not local.exists():
        raise SystemExit(f"missing eQTL file: {local.name}\n"
                         f"  this check needs the file locally; "
                         f"15_locus_screen.py download {qtd} fetches it")
    expected = remote_size(f"{EQTL_CATALOGUE_BASE}/{qtd}/{qtd}.all.tsv.gz")
    if expected and local.stat().st_size != expected:
        raise SystemExit(
            f"local eQTL file is incomplete: {local.name} has "
            f"{local.stat().st_size:,} of {expected:,} bytes\n"
            f"  re-download it before running; tabix would otherwise return "
            f"an empty result without an error")
    return local


def gene_index(qtd: str, local: pathlib.Path,
               override: Optional[str] = None) -> pd.DataFrame:
    """Return each gene's tested cis window for one dataset.

    Parameters
    ----------
    qtd : eQTL Catalogue dataset id.
    local : the all-associations file; used to build the table when needed.
    override : explicit gene window table (gene_id, chrom, lo, hi).

    Returns
    -------
    DataFrame with columns gene_id, chrom, lo, hi.

    Notes
    -----
    Building the table costs one pass over a multi-gigabyte file, so the result
    is cached under PATHS.cache and reused. The cache is keyed by dataset id, not
    by file size: a re-download of the same dataset keeps the same cis windows.
    """
    if override:
        path = resolve_input(override)
        if not path.exists():
            raise SystemExit(f"missing gene index: {path}")
        return pd.read_csv(path, sep="\t")

    cached = PATHS.cache / f"gene_index_{qtd}.tsv"
    if cached.exists():
        return pd.read_csv(cached, sep="\t")

    print(f"      building gene index for {qtd} (one pass over the file)",
          flush=True)
    idx = _eqtl.build_gene_index(local)
    idx.to_csv(cached, sep="\t", index=False)
    return idx


def control_table(gwas_path: pathlib.Path, tissue: str,
                  genes: Dict[str, str]) -> pd.DataFrame:
    """Run the ABF test on a fixed gene panel.

    Parameters
    ----------
    gwas_path : outcome file; the layout is detected from its header.
    tissue : GTEx tissue label, resolved through config.TISSUES.
    genes : gene symbol -> Ensembl id.

    Returns
    -------
    DataFrame with one row per gene that could be tested, strongest PP.H4 first.

    Notes
    -----
    Every gene is tested with the same code path as the primary analysis, which
    is the whole point of a positive control: a different implementation here
    would not license a conclusion about the primary analysis.
    """
    qtd = TISSUES.get(tissue)
    if qtd is None:
        raise SystemExit(f"unknown tissue {tissue}; add it to config.TISSUES")

    local = local_eqtl(qtd)
    print(f"      eQTL source: {local.name}", flush=True)
    idx = gene_index(qtd, local)
    want = _eqtl.column_positions(local)

    store = _gwas.load_gwas(gwas_path)
    print(f"      outcome: {sum(len(v) for v in store.values()):,} variants",
          flush=True)

    records: List[Dict[str, object]] = []
    for symbol, gene in genes.items():
        row = idx[idx["gene_id"] == gene]
        if row.empty:
            print(f"  {symbol:8s} not in the gene index, skipped", flush=True)
            continue
        r = row.iloc[0]

        q = _eqtl.qtl_region(str(r["chrom"]), int(r["lo"]), int(r["hi"]),
                             str(local), want)
        q = q[q["gene_id"] == gene].dropna(subset=["beta", "se", "pvalue"])
        q = q[(q["maf"] >= PARAMS.eqtl_maf_min)
              & (q["maf"] <= PARAMS.eqtl_maf_max) & (q["se"] > 0)]
        if len(q) < PARAMS.eqtl_min_snps:
            print(f"  {symbol:8s} too few eQTL variants ({len(q)})", flush=True)
            continue

        pairs = []
        for x in q.itertuples(index=False):
            hit = _gwas.lookup(store, str(r["chrom"]), int(x.position),
                               str(x.ref).upper(), str(x.alt).upper())
            if hit is None:
                continue
            beta_g, se_g, p_g = hit
            if np.isfinite(beta_g) and np.isfinite(se_g) and se_g > 0:
                pairs.append((x.beta, x.se, beta_g, se_g, p_g))
        if len(pairs) < PARAMS.eqtl_min_snps:
            print(f"  {symbol:8s} too few pairable variants ({len(pairs)})",
                  flush=True)
            continue

        d = pd.DataFrame(pairs, columns=["b_e", "se_e", "b_g", "se_g", "p_g"])
        res = _abf.coloc_abf(
            _abf.log_abf(d["b_e"].values, d["se_e"].values, 1.0),
            _abf.log_abf(d["b_g"].values, d["se_g"].values, 1.0),
            p1=PARAMS.p1, p2=PARAMS.p2, p12=PARAMS.p12,
        )
        records.append({
            "PP.H0": res["PP.H0"], "PP.H1": res["PP.H1"], "PP.H2": res["PP.H2"],
            "PP.H3": res["PP.H3"], "PP.H4": res["PP.H4"],
            "n_snps": len(d),
            "top_snp_share": res["top_snp_share"],
            "gene": symbol, "gene_id": gene, "tissue": tissue,
            "best_gwas_p": float(d["p_g"].min()),
        })
        print(f"  {symbol:8s} n={len(d):5d}  H4={res['PP.H4']:.3f}  "
              f"H3={res['PP.H3']:.3f}  best outcome p={d['p_g'].min():.1e}",
              flush=True)

    out = pd.DataFrame(records)
    if len(out):
        out = out.sort_values("PP.H4", ascending=False).reset_index(drop=True)
    return out


def load_power_region(tissue: str) -> pd.DataFrame:
    """Read the simulation window: every tested variant of one gene.

    Parameters
    ----------
    tissue : GTEx tissue label whose eQTL file supplies the variants.

    Returns
    -------
    DataFrame with pos, beta, se, pvalue, maf sorted by position.
    """
    qtd = TISSUES.get(tissue)
    if qtd is None:
        raise SystemExit(f"unknown tissue {tissue}; add it to config.TISSUES")
    local = local_eqtl(qtd)

    d = _eqtl.qtl_region(POWER_CHROM, POWER_LO, POWER_HI, str(local),
                         _eqtl.column_positions(local))
    d = d[d["gene_id"] == POWER_GENE].dropna(
        subset=["position", "beta", "se", "pvalue"])
    d = d[(d["maf"] >= PARAMS.eqtl_maf_min)
          & (d["maf"] <= PARAMS.eqtl_maf_max) & (d["se"] > 0)]
    return d.rename(columns={"position": "pos"}).sort_values("pos")


def one_trial(d: pd.DataFrame, z_eqtl: float, z_gwas: float, causal_idx: int,
              rng: np.random.Generator) -> Dict[str, float]:
    """One simulated pair of traits sharing a causal variant.

    Parameters
    ----------
    d : the real variant set; only ``se`` and the row count are used.
    z_eqtl, z_gwas : target z statistics at the causal variant.
    causal_idx : row of the shared causal variant.
    rng : generator, so that the caller controls reproducibility.

    Returns
    -------
    dict with PP.H4, PP.H3 and the number of variants.

    Notes
    -----
    Noise is drawn from the real per-variant standard errors, so the null
    distribution inherits the true allele-frequency and imputation profile of
    the window. Only the two effect sizes are simulated.
    """
    n = len(d)
    se_e = d["se"].to_numpy(dtype=float).copy()

    be = rng.normal(0.0, se_e)
    be[causal_idx] = z_eqtl * se_e[causal_idx]

    # The outcome's standard errors are not the expression's; a 0.7-1.3 multiple
    # spans the ratio observed between the real outcome and the eQTL file.
    se_g = se_e * rng.uniform(0.7, 1.3, n)
    bg = rng.normal(0.0, se_g)
    bg[causal_idx] = z_gwas * se_g[causal_idx]

    res = _abf.coloc_abf(_abf.log_abf(be, se_e, 1.0),
                         _abf.log_abf(bg, se_g, 1.0),
                         p1=PARAMS.p1, p2=PARAMS.p2, p12=PARAMS.p12)
    return {"PP.H4": res["PP.H4"], "PP.H3": res["PP.H3"], "n": n}


def power_curve(tissue: str, trials: int, seed: int) -> pd.DataFrame:
    """Build the eQTL-strength to PP.H4 curve.

    Parameters
    ----------
    tissue : GTEx tissue supplying the real variant set.
    trials : repetitions per cell.
    seed : RNG seed.

    Returns
    -------
    DataFrame over POWER_COLS.

    Notes
    -----
    The z grid is chosen to bracket the observed distribution of GTEx uterine
    eGene leads, whose median strongest eQTL z is about 5.3: the curve therefore
    answers "what would a signal of the strength we actually have achieve",
    which is the question a negative colocalisation result has to confront.
    """
    d = load_power_region(tissue)
    if len(d) < PARAMS.eqtl_min_snps:
        raise SystemExit(f"only {len(d)} variants in the power window; "
                         f"check {POWER_GENE} and the coordinates")

    # The shared causal variant is the true expression lead.
    causal_idx = int(np.argmin(d["pvalue"].to_numpy()))
    print(f"[power] window {POWER_GENE}: {len(d):,} variants, "
          f"shared causal variant at {int(d['pos'].iloc[causal_idx]):,}",
          flush=True)

    rng = np.random.default_rng(seed)
    rows: List[Dict[str, object]] = []
    for z_e in (4.0, 5.3, 6.5, 8.0, 10.0, 13.0):
        for z_g in (10.0, 15.0, 25.0):
            draws = [one_trial(d, z_e, z_g, causal_idx, rng)
                     for _ in range(trials)]
            h4 = np.array([x["PP.H4"] for x in draws])
            h3 = np.array([x["PP.H3"] for x in draws])
            rows.append({
                "eQTL_z": z_e,
                "GWAS_z": z_g,
                "H4_median": round(float(np.median(h4)), 3),
                "H4_p25": round(float(np.percentile(h4, 25)), 3),
                "H4_p75": round(float(np.percentile(h4, 75)), 3),
                "detect_rate_H4_gt_0.8": (
                    f"{100 * float(np.mean(h4 > PARAMS.pp_h4_high)):.0f}%"),
                "H3_median": round(float(np.median(h3)), 3),
            })
    return pd.DataFrame(rows)[POWER_COLS]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--section", default="both",
                    choices=["control", "power", "both"])
    ap.add_argument("--tissue", default="Pancreas",
                    help="GTEx tissue for the positive control")
    ap.add_argument("--gwas", default="finngen_R12_T2D.gz",
                    help="outcome file in data/ for the control")
    ap.add_argument("--set", default="canonical", choices=["canonical", "bmi"],
                    help="which gene panel the control runs")
    ap.add_argument("--index", default=None, help="pre-built gene window table")
    ap.add_argument("--trials", type=int, default=60,
                    help="Monte-Carlo repetitions per cell")
    ap.add_argument("--seed", type=int, default=SEED,
                    help="RNG seed for the simulation")
    args = ap.parse_args()

    PATHS.ensure()
    t0 = time.time()

    if args.section in ("control", "both"):
        print("[control] canonical gene panel", flush=True)
        gwas_path = resolve_input(args.gwas)
        if not gwas_path.exists():
            raise SystemExit(f"missing outcome file: {gwas_path}")
        genes = BMI_GENES if args.set == "bmi" else CANONICAL
        df = control_table(gwas_path, args.tissue, genes)
        if df.empty:
            raise SystemExit("no gene could be tested; the panel and the tissue "
                             "may not correspond to this eQTL dataset")
        dest = PATHS.tables / "positive_control_t2d.tsv"
        df.to_csv(dest, sep="\t", index=False)
        n_high = int((df["PP.H4"] > PARAMS.pp_h4_high).sum())
        print(f"      {len(df)} genes -> {dest.name}; "
              f"PP.H4>{PARAMS.pp_h4_high:g} in {n_high}", flush=True)
        if n_high == 0:
            print("      WARNING: the control recovered nothing. A negative "
                  "result elsewhere cannot be interpreted until this is fixed.",
                  flush=True)

    if args.section in ("power", "both"):
        print("[power] simulation", flush=True)
        out = power_curve(args.tissue, args.trials, args.seed)
        print(out.to_string(index=False), flush=True)
        dest = PATHS.tables / "coloc_power.tsv"
        out.to_csv(dest, sep="\t", index=False)
        print(f"      {len(out)} cells -> {dest.name}", flush=True)

    print(f"total {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
