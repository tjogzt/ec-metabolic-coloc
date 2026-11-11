"""
08_prior_sensitivity.py
=======================

Purpose
-------
Repeat every colocalisation test across a grid of the shared-causal-variant
prior, and report four diagnostics that the headline PP.H4 alone does not
expose:

    1. p12 sensitivity   PP.H4 recomputed over config.PARAMS.p12_sensitivity
                         with p1 and p2 held fixed
    2. conditional H4    PP.H4 / (PP.H3 + PP.H4): separates "one shared variant"
                         from "both traits have a signal but at distinct variants"
    3. credible set      variants carrying 95% of the shared-signal posterior,
                         with the size and the top posterior weight
    4. Bayesian FDR      expected number of false discoveries among the tests
                         declared positive, computed as sum(1 - PP.H4)

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
The locus-restricted eQTL subsets written by 07_run_colocalisation.py, cached
under PATHS.cache, plus the outcome GWAS.

Outputs
-------
results/tables/coloc_robustness.tsv, one row per test per p12 value, with
columns for all four diagnostics above.

Notes
-----
The credible set reported here follows from the single-causal-variant
assumption of the ABF model. It is therefore *not* a SuSiE credible set and must
not be described as one; the column names carry the distinction.

H3 is evaluated through the identity

    sum_{i!=j} exp(l1_i + l2_j) = (sum_i exp l1_i)(sum_j exp l2_j)
                                  - sum_i exp(l1_i + l2_i)

which avoids materialising an n by n matrix. With several hundred variants per
gene the dense form is both slower and needlessly memory hungry.
"""
from __future__ import annotations

import math
import pathlib
import sys
import time
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

from importlib import import_module

from config.config import PATHS, PARAMS, TISSUES, GWAS_PRIMARY_FILE

_abf = import_module("05_coloc_abf_engine")
_gwas = import_module("02_gwas_reader")
_loci = import_module("06_define_risk_loci")


def _logsum(x: np.ndarray) -> float:
    """Numerically stable log(sum(exp(x))); returns -inf for an empty input."""
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return -np.inf
    m = float(np.max(x))
    if not np.isfinite(m):
        return -np.inf
    return float(m + math.log(np.sum(np.exp(x - m))))


def posteriors(l1: np.ndarray, l2: np.ndarray, p12: float,
               p1: float = PARAMS.p1, p2: float = PARAMS.p2) -> Dict[str, float]:
    """Posterior probabilities with a freely varying p12.

    Parameters
    ----------
    l1, l2 : log ABF of the same variants in the two traits.
    p12 : prior probability of a shared causal variant; the quantity being
        varied across the sensitivity grid.
    p1, p2 : prior probabilities that each trait has a signal.

    Returns
    -------
    dict with keys H0 to H4.

    Notes
    -----
    Identical to 05_coloc_abf_engine.coloc_abf but with the H3 sum evaluated by
    the identity rather than by an outer product, and with the log(1 - p12)
    guarded so that an extreme prior cannot produce a log of a negative number.
    """
    lh0 = math.log(max(1e-300, 1.0 - p1 - p2 - p12))
    L1, L2, L4 = _logsum(l1), _logsum(l2), _logsum(l1 + l2)
    lh1 = L1 + math.log(p1) + math.log(1.0 - p2)
    lh2 = L2 + math.log(p2) + math.log(1.0 - p1)

    # delta = log( sum(e^{l1+l2}) / (sum e^{l1} * sum e^{l2}) ) <= 0.
    # log1p keeps precision when the shared support is a tiny fraction of the total.
    delta = min(0.0, L4 - L1 - L2)
    lh3 = (math.log(p1) + math.log(p2) + L1 + L2 + math.log1p(-math.exp(delta))
           if delta < 0 else -math.inf)
    lh4 = L4 + math.log(p12)

    v = np.array([lh0, lh1, lh2, lh3, lh4], dtype=float)
    v = np.exp(v - _logsum(v))
    return {"H0": float(v[0]), "H1": float(v[1]), "H2": float(v[2]),
            "H3": float(v[3]), "H4": float(v[4])}


def credible_set(l1: np.ndarray, l2: np.ndarray, indices: Sequence[int],
                 threshold: float = 0.95) -> Tuple[List[int], float, int]:
    """Credible set implied by the shared-signal posterior.

    Parameters
    ----------
    l1, l2 : log ABF of the two traits.
    indices : caller-supplied identifiers, returned for the set members.
    threshold : cumulative posterior mass required, conventionally 0.95.

    Returns
    -------
    (member_indices, top_posterior, top_index): the variants in the set, the
    posterior weight of the single best variant, and its index into ``indices``.

    Notes
    -----
    Variant posterior is proportional to exp(l1 + l2), so the set is built by
    ranking on the shared log ABF and accumulating until the threshold is met.
    A large set with a small top weight means the shared signal is diffuse,
    which is reported rather than hidden.
    """
    w = l1 + l2
    w = np.exp(w - np.max(w))
    w = w / w.sum()
    order = np.argsort(-w)
    cumulative, members = 0.0, []
    for j in order:
        members.append(int(indices[j]))
        cumulative += float(w[j])
        if cumulative >= threshold:
            break
    return members, float(w[order[0]]), int(order[0])


def main() -> None:
    PATHS.ensure()
    t0 = time.time()

    gwas_path = PATHS.data / GWAS_PRIMARY_FILE
    store = _gwas.load_gwas(gwas_path)
    print(f"[1/3] outcome loaded: {sum(len(v) for v in store.values()):,} variants", flush=True)

    loci = _gwas.peaks(store, PARAMS.gwas_p_thresh, PARAMS.locus_gap_bp)
    print(f"      {len(loci)} loci", flush=True)

    tar_path = PATHS.data / "GTEx_Analysis_v8_eQTL.tar"
    rows: List[Dict[str, object]] = []

    for tissue, qtl_id in TISSUES.items():
        subset = PATHS.cache / f"loci_{qtl_id}_{tissue}.tsv"
        if not subset.exists():
            print(f"[2/3] {tissue}: no cached subset ({subset.name}), skipped", flush=True)
            continue

        genes_by_locus = _loci.genes_near_loci(tar_path, tissue, loci,
                                               PARAMS.region_pad_bp)
        qtl = pd.read_csv(subset, sep="\t")
        qtl = qtl.dropna(subset=["beta", "se"])
        qtl = qtl[(qtl["maf"] >= PARAMS.eqtl_maf_min)
                  & (qtl["maf"] <= PARAMS.eqtl_maf_max) & (qtl["se"] > 0)]

        n_done = 0
        for i, row in loci.iterrows():
            genes = genes_by_locus.get(i, [])
            if not genes:
                continue
            lead = int(row.lead)
            near = qtl[(qtl["chromosome"].astype(str) == str(row.chrom))
                       & (qtl["position"].between(lead - PARAMS.region_pad_bp,
                                                  lead + PARAMS.region_pad_bp))]
            if near.empty:
                continue

            for gene in genes:
                sub = near[near["gene_id"] == gene]
                if len(sub) < PARAMS.eqtl_min_snps:
                    continue
                b_e, se_e, b_g, se_g, pos, rsid = [], [], [], [], [], []
                for x in sub.itertuples(index=False):
                    hit = _gwas.lookup(store, str(row.chrom), int(x.position),
                                       str(x.ref).upper(), str(x.alt).upper())
                    if hit is None:
                        continue
                    beta_g, se_gene, _p = hit
                    if np.isfinite(beta_g) and np.isfinite(se_gene) and se_gene > 0:
                        b_e.append(x.beta); se_e.append(x.se)
                        b_g.append(beta_g); se_g.append(se_gene)
                        pos.append(int(x.position)); rsid.append(str(getattr(x, "rsid", "")))
                if len(b_e) < PARAMS.eqtl_min_snps:
                    continue

                l1 = _abf.log_abf(np.array(b_e), np.array(se_e), 1.0)
                l2 = _abf.log_abf(np.array(b_g), np.array(se_g), 1.0)
                members, top_w, top_j = credible_set(l1, l2, range(len(b_e)))

                base = {
                    "tissue": tissue,
                    "locus": f"{row.chrom}:{lead}",
                    "gene_id": gene,
                    "n_snps": len(b_e),
                    "cred_set_size": len(members),
                    "cred_set_top_pos": pos[top_j],
                    "cred_set_top_rsid": rsid[top_j],
                    "cred_set_top_post": round(top_w, 4),
                    "cred_set_positions": ";".join(str(pos[j]) for j in members),
                }
                for p12 in PARAMS.p12_sensitivity:
                    po = posteriors(l1, l2, p12)
                    denom = po["H3"] + po["H4"]
                    rows.append(dict(base, p12=p12, PP_H3=po["H3"], PP_H4=po["H4"],
                                     cond_H4=(po["H4"] / denom) if denom > 0 else np.nan))
                n_done += 1
        print(f"[2/3] {tissue}: {n_done} tests ({time.time() - t0:.0f}s)", flush=True)

    if not rows:
        raise SystemExit("no tests produced; run 07_run_colocalisation.py first")

    out = pd.DataFrame(rows)
    dest = PATHS.tables / "coloc_robustness.tsv"
    out.to_csv(dest, sep="\t", index=False)

    # ---- summary -------------------------------------------------------
    print(f"[3/3] {len(out)} rows -> {dest.name}", flush=True)
    base = out[out["p12"] == PARAMS.p12]
    print(f"      baseline p12={PARAMS.p12:g}: {len(base)} tests, "
          f"PP.H4>{PARAMS.pp_h4_high}: {(base['PP_H4'] > PARAMS.pp_h4_high).sum()}, "
          f"PP.H3>{PARAMS.pp_h3_high}: {(base['PP_H3'] > PARAMS.pp_h3_high).sum()}",
          flush=True)

    print("      prior sensitivity:", flush=True)
    for p12 in PARAMS.p12_sensitivity:
        s = out[out["p12"] == p12]
        fdr = s.loc[s["PP_H4"] > PARAMS.pp_h4_high, "PP_H4"].rsub(1).sum()
        print(f"        p12={p12:<8g} H4>0.8: {(s['PP_H4'] > PARAMS.pp_h4_high).sum():3d}"
              f"  H4>0.5: {(s['PP_H4'] > 0.5).sum():3d}"
              f"  expected false discoveries: {fdr:.3f}", flush=True)


if __name__ == "__main__":
    main()
