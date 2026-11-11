"""
04_expression_matrix.py
=======================

Purpose
-------
Audit which GEO series can be analysed for expression at all, and inventory the
endometrial-cancer series that could enter the bulk layer.

Three questions are answered, one per step, because they cost very different
amounts of network time and are useful at different moments:

    matrix     does a series ship a pre-processed expression matrix (series
               matrix), or only CEL / RAW archives? Only the former can be used
               without a local array-processing stack, so this is checked before
               any cohort is chosen.
    inventory  search GDS for endometrial-cancer series, filter on
               "Expression profiling" and a minimum sample count, then spot-check
               the first few by reading their sample titles. The spot check is
               what separates a genuine tumour cohort from a study that merely
               mentions endometrial cancer in its abstract.
    census     recount samples per series by parsing the GEO text record
               directly. esummary is not used for this because it mixes nested
               Sample nodes into the count and inflates it.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
Network access to NCBI (GEO FTP listing, E-utilities, acc.cgi). No local files.

Outputs
-------
PATHS.tables/cohort_candidates.tsv, written by the census step, one row per
series that passed the filter:

    gse, n, pmid, title, types

The matrix and inventory steps print to stdout and write nothing.

Usage
-----
    python code/03_burden/04_expression_matrix.py --step matrix
    python code/03_burden/04_expression_matrix.py --step inventory
    python code/03_burden/04_expression_matrix.py --step census

Notes
-----
The default step is ``matrix`` because it is the cheapest and the most useful:
ten HEAD-level requests. ``inventory`` and ``census`` issue hundreds of requests
to NCBI and are therefore opt-in.

Requests are made through the ``curl`` binary rather than urllib so that a
stalled connection is bounded by an explicit per-request timeout instead of
hanging the whole audit; NCBI drops individual GEO connections routinely.
"""
from __future__ import annotations

import argparse
import pathlib
import re
import subprocess
import sys
import time
import urllib.parse
import xml.etree.ElementTree as ET
from typing import Dict, List

# Repository root, so that `config` resolves when the file is run as a script.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

from config.config import PATHS

# --------------------------------------------------------------------------
# Endpoints. These are service URLs, not filesystem paths: they are part of the
# protocol being spoken and live here rather than in config.
# --------------------------------------------------------------------------
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GEO_SERIES_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
GEO_ACC_CGI = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

# --------------------------------------------------------------------------
# Module constants
# --------------------------------------------------------------------------
# Series whose matrix availability was in question when the cohorts were being
# assembled: the four that later entered the analysis, plus the alternatives
# that were considered and rejected.
MATRIX_PROBE_SERIES = (
    "GSE120490", "GSE17025", "GSE106191", "GSE196033", "GSE56026",
    "GSE216872", "GSE14860", "GSE178671", "GSE32507", "GSE21882",
)

# Search strings. Both require the series to be human and to be a Series record;
# the census variant adds a fourth title synonym, which is why it finds a few
# more hits than the inventory search.
SEARCH_TERM_INVENTORY = (
    '(endometrial cancer[Title] OR endometrial carcinoma[Title] '
    'OR endometrioid[Title]) AND "Homo sapiens"[Organism] AND gse[Entry Type]'
)
SEARCH_TERM_CENSUS = (
    '(endometrial cancer[Title] OR endometrial carcinoma[Title] '
    'OR endometrioid[Title] OR endometrial neoplasms[Title]) '
    'AND "Homo sapiens"[Organism] AND gse[Entry Type]'
)

# A bulk cohort has to carry enough samples to support a two-group contrast;
# 40 is the floor below which the variance estimate stops being usable.
MIN_SAMPLES = 40

# Result caps. The inventory stops at 120 hits, the census at 200: beyond that
# the per-series verification dominates and the extra hits are off-topic.
RETMAX_INVENTORY = 120
RETMAX_CENSUS = 200

# Number of top hits whose sample titles are read in the inventory step, and how
# many titles are shown per series.
SPOT_CHECK_SERIES = 12
SPOT_CHECK_TITLES = 4

# Polite delay between requests, in seconds. NCBI throttles anonymous clients
# above three requests per second.
SLEEP_SERIES = 0.5
SLEEP_CENSUS = 0.3
SLEEP_SPOT = 0.6

# Per-request timeouts, in seconds.
TIMEOUT_LISTING = 60
TIMEOUT_EUTILS = 90
TIMEOUT_SPOT = 180

# Progress cadence of the census loop.
CENSUS_PROGRESS_EVERY = 40


# --------------------------------------------------------------------------
# Transport
# --------------------------------------------------------------------------
def curl(url: str, timeout: int = TIMEOUT_EUTILS) -> str:
    """Fetch a URL with ``curl`` and return the body as text.

    Parameters
    ----------
    url : absolute URL.
    timeout : per-request cap in seconds.

    Returns
    -------
    Response body, or an empty string if curl failed. Failures are not raised
    because the callers are audits that must survive one dead accession.
    """
    proc = subprocess.run(["curl", "-sS", "-m", str(timeout), url],
                          capture_output=True, text=True)
    return proc.stdout


# --------------------------------------------------------------------------
# Step 1: is a series matrix available?
# --------------------------------------------------------------------------
def matrix_availability(series=MATRIX_PROBE_SERIES) -> Dict[str, List[str]]:
    """List the series-matrix files published for each candidate series.

    Parameters
    ----------
    series : GEO accessions to probe.

    Returns
    -------
    dict mapping accession to a list of "file (size)" strings; an accession with
    no series matrix maps to an empty list.

    Notes
    -----
    GEO organises the directory tree by accession prefix, with the last three
    digits replaced by "nnn" (GSE17025 -> GSE17nnn). The listing is parsed with
    a regex rather than an HTML parser because the page is a plain Apache index.
    The date-and-time block before the size is optional in the pattern, so the
    captured token is the file size rather than the year of the listing.
    """
    found: Dict[str, List[str]] = {}
    pattern = re.compile(
        r'href="([^"]*series_matrix[^"]*)"[^<]*</a>'
        r'\s+(?:\d{4}-\d{2}-\d{2}\s+[\d:]{4,5}\s+)?'   # optional date and time
        r'([\d\.]+[KMG]?)')
    for acc in series:
        stub = acc[:-3] + "nnn"
        url = f"{GEO_SERIES_FTP}/{stub}/{acc}/matrix/"
        body = curl(url, TIMEOUT_LISTING)
        files = pattern.findall(body)
        found[acc] = [f"{name} ({size})" for name, size in files]
        # An empty body means the request failed, which is not the same finding
        # as an absent file: GEO drops individual connections routinely, and
        # reporting a failed request as "no matrix" would be a false negative.
        if not body.strip():
            report = "listing unavailable (request failed)"
        elif found[acc]:
            report = "; ".join(found[acc])
        else:
            report = "no series matrix"
        print(f"{acc:<12} {report}")
    return found


# --------------------------------------------------------------------------
# Step 2: inventory via E-utilities
# --------------------------------------------------------------------------
def esearch(term: str, retmax: int = RETMAX_INVENTORY) -> List[str]:
    """Return GDS UIDs matching a search term.

    Parameters
    ----------
    term : Entrez query, already in GDS syntax.
    retmax : maximum number of UIDs.

    Returns
    -------
    list of UID strings.
    """
    url = f"{EUTILS}/esearch.fcgi?db=gds&retmax={retmax}&term={urllib.parse.quote(term)}"
    return re.findall(r"<Id>(\d+)</Id>", curl(url))


def esummary(ids: List[str]) -> List[dict]:
    """Fetch DocSum records for a batch of UIDs.

    Parameters
    ----------
    ids : UIDs from :func:`esearch`.

    Returns
    -------
    list of dicts, one per record, keyed by the Item ``Name`` attribute.

    Notes
    -----
    Requests are batched at 100 ids, the E-utilities limit for esummary.
    """
    out: List[dict] = []
    for i in range(0, len(ids), 100):
        url = f"{EUTILS}/esummary.fcgi?db=gds&id={','.join(ids[i:i + 100])}"
        root = ET.fromstring(curl(url))
        for doc in root.iter("DocSum"):
            d = {}
            for item in doc.iter("Item"):
                d[item.get("Name")] = (item.text or "").strip()
            out.append(d)
    return out


def inventory() -> List[tuple]:
    """Inventory candidate series and spot-check the largest few.

    Returns
    -------
    list of (accession, n_samples, title, first type) tuples, sorted by sample
    count descending, restricted to expression-profiling series at or above
    ``MIN_SAMPLES``.

    Notes
    -----
    Two filters are applied, in this order: the series must be an expression
    profiling study (a methylation or ChIP series with the right title would
    otherwise slip through), and it must carry at least ``MIN_SAMPLES`` samples.
    Sample titles are then read for the largest few because the filters cannot
    tell an endometrial carcinoma cohort from, say, a decidua cohort with a
    similar title.
    """
    ids = esearch(SEARCH_TERM_INVENTORY, RETMAX_INVENTORY)
    rows = esummary(ids)
    time.sleep(SLEEP_SERIES)

    cand: List[tuple] = []
    for r in rows:
        acc = r.get("accession", "")
        try:
            n = int(r.get("n_samples", "0"))
        except ValueError:
            n = 0
        typ = r.get("gdstype", "")
        if "Expression profiling" in typ and n >= MIN_SAMPLES:
            cand.append((acc, n, r.get("title", "")[:88], typ.split(";")[0]))

    cand.sort(key=lambda x: -x[1])
    print(f"{len(rows)} records retrieved; {len(cand)} are expression profiling "
          f"with n>={MIN_SAMPLES}:\n")
    print(f"{'GSE':<14}{'n':>5}  title")
    for acc, n, title, _ in cand:
        print(f"{acc:<14}{n:>5}  {title}")

    # Sample-title spot check: four titles and the organism per series.
    print("\n" + "=" * 74)
    print(f"sample titles of the {SPOT_CHECK_SERIES} largest series "
          f"({SPOT_CHECK_TITLES} each)")
    print("=" * 74)
    for acc, _, _, _ in cand[:SPOT_CHECK_SERIES]:
        txt = curl(f"{GEO_ACC_CGI}?acc={acc}&targ=gsm&form=text&view=brief",
                   TIMEOUT_SPOT)
        titles = re.findall(r"!Sample_title = (.+)", txt)
        org = set(re.findall(r"!Sample_organism_ch1 = (.+)", txt))
        print(f"\n{acc} (n={len(titles)}, organism={','.join(sorted(org))[:24]})")
        print("   ", " | ".join(t.strip() for t in titles[:SPOT_CHECK_TITLES]))
        time.sleep(SLEEP_SPOT)
    return cand


# --------------------------------------------------------------------------
# Step 3: census by parsing the series record
# --------------------------------------------------------------------------
def series_brief(acc: str) -> dict:
    """Parse one series record in GEO text format.

    Parameters
    ----------
    acc : series accession.

    Returns
    -------
    dict with keys ``acc``, ``n``, ``type`` (list), ``title``, ``pubmed``.

    Notes
    -----
    ``n`` counts ``Series_sample_id`` lines, which is one line per sample in the
    self view of the record. This is the count esummary cannot give: its
    n_samples field is derived from a nested tree that includes child records.
    """
    url = f"{GEO_ACC_CGI}?acc={acc}&targ=self&form=text&view=brief"
    txt = curl(url)
    d = {"acc": acc, "n": 0, "type": [], "title": "", "pubmed": ""}
    for line in txt.splitlines():
        if not line.startswith("!"):
            continue
        key, _, value = line[1:].partition(" = ")
        key, value = key.strip(), value.strip()
        if key == "Series_sample_id":
            d["n"] += 1
        elif key == "Series_type":
            d["type"].append(value)
        elif key == "Series_title":
            d["title"] = value
        elif key == "Series_pubmed_id":
            d["pubmed"] = value
    return d


def census() -> "list[dict]":
    """Count samples for every hit of the search and write the candidate table.

    Returns
    -------
    list of candidate records, sorted by sample count descending.

    Notes
    -----
    The UID-to-accession conversion is arithmetic: a GDS UID is
    200000000 + series number, so GSE12345 has UID 200012345 and the accession
    is recovered by dropping the leading "200".
    """
    PATHS.ensure()
    ids = esearch(SEARCH_TERM_CENSUS, RETMAX_CENSUS)
    accs = [f"GSE{uid[3:]}" for uid in ids]
    print(f"{len(accs)} series matched; verifying each ...\n")

    rows: List[dict] = []
    for i, acc in enumerate(accs, 1):
        d = series_brief(acc)
        if any("Expression profiling" in t for t in d["type"]) and d["n"] >= MIN_SAMPLES:
            rows.append(d)
        if i % CENSUS_PROGRESS_EVERY == 0:
            print(f"  ... {i}/{len(accs)}, {len(rows)} kept")
    rows.sort(key=lambda x: -x["n"])

    print(f"\n{len(rows)} expression-profiling series with n>={MIN_SAMPLES}:\n")
    print(f"{'GSE':<13}{'n':>5}  {'PMID':<10}title")
    for d in rows:
        print(f"{d['acc']:<13}{d['n']:>5}  {d['pubmed']:<10}{d['title'][:80]}")

    dest = PATHS.tables / "cohort_candidates.tsv"
    with open(dest, "w") as fh:
        fh.write("gse\tn\tpmid\ttitle\ttypes\n")
        for d in rows:
            fh.write(f"{d['acc']}\t{d['n']}\t{d['pubmed']}\t{d['title']}\t"
                     f"{';'.join(d['type'])}\n")
    print(f"\nwritten {dest}")
    return rows


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------
def main() -> None:
    """Parse arguments and run the requested step."""
    ap = argparse.ArgumentParser(
        description="GEO availability audit for the bulk cohorts.")
    ap.add_argument("--step", choices=("matrix", "inventory", "census"),
                    default="matrix",
                    help="which audit to run (default: matrix; the other two "
                         "issue hundreds of NCBI requests)")
    args = ap.parse_args()

    if args.step == "matrix":
        matrix_availability()
    elif args.step == "inventory":
        inventory()
    else:
        census()


if __name__ == "__main__":
    main()
