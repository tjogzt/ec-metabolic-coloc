"""
09_curate_signature.py
======================

Purpose
-------
Curate the candidate GEO series by reading each series record and classifying
what the study actually is, so that the single-cell and bulk layers draw from
cohorts of the intended tissue and disease.

Classification is prioritised from the disease-specific labels to the
tissue-general ones:

    EC_tumor -> endometriosis -> organoid -> endometrium -> decidua -> off_topic

The order has two consequences that are easy to get wrong, and both are
deliberate:

    an endometrial-carcinoma organoid study is labelled ``EC_tumor``, not
    ``organoid``, because the tumour rule is matched first -- it is then treated
    as tumour tissue even though a model system has no microenvironment;

    a decidualisation study of endometrial stromal cells is labelled
    ``endometrium``, not ``decidua``, for the same reason.

The rules were fixed against the candidate list before any series was read, and
are reproduced here unchanged so that a re-run with the same list gives the same
labels.

Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)
Created  : 2026-11-11
Version  : 2.0.0

Inputs
------
PATHS.tables/<candidates>.tsv, default ``w3_candidates.tsv``, produced by the
GEO searches in stage 04 or assembled by hand. Required columns: ``gse`` and
``title``.

Network access to NCBI acc.cgi for each candidate.

Outputs
-------
PATHS.tables/w3_curated.tsv, one row per candidate, columns in this order:

    gse, n, cls, type, title

where ``n`` is the number of samples counted from the series record and ``cls``
is the classification label.

Usage
-----
    python code/03_burden/09_curate_signature.py
    python code/03_burden/09_curate_signature.py --candidates cohort_candidates.tsv

Notes
-----
The sample count is taken from the number of ``Series_sample_id`` values, not
from esummary's n_samples, which counts nested child records as samples and
inflates series that use sub-series.

A failed request returns an empty string and the series is classified as
off_topic. That is the conservative failure direction: a retrieval error cannot
promote a series into the tumour class.
"""
from __future__ import annotations

import argparse
import pathlib
import re
import sys
import time
import urllib.request
from typing import Dict, List, Tuple

import pandas as pd

# Repository root, so that `config` resolves when the file is run as a script.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

from config.config import PATHS

# --------------------------------------------------------------------------
# Module constants
# --------------------------------------------------------------------------
GEO_ACC_CGI = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

# Identifies the client to NCBI; not a secret and not a contact address.
USER_AGENT = {"User-Agent": "Mozilla/5.0 (research script)"}

# Retries per request and the per-attempt timeout in seconds. Backoff is
# linear, because the failures are rate limits rather than congestion.
REQUEST_TRIES = 3
REQUEST_TIMEOUT = 60
RETRY_BASE_SLEEP = 2

# Pause between series, to stay under NCBI's anonymous rate limit.
SLEEP_PER_SERIES = 0.3

# Classification rules, ordered. The first match wins.
RULES: List[Tuple[str, str]] = [
    ("EC_tumor",
     r"endometrial (cancer|carcinoma|tumou?r|neoplas)|endometrioid|"
     r"uterine (cancer|carcinoma|serous)|carcinosarcoma"),
    ("endometriosis",
     r"endometrio(sis|tic)|menstrual (blood|effluent|fluid)|retrograde menstruation"),
    ("organoid",
     r"organoid|assembloid|blastoid|cell line|cultured |E6/E7|TERT-1"),
    ("endometrium",
     r"endometri(um|al)|uterus|uterine|Asherman|recurrent implantation failure|"
     r"endometrial polyp|PCOS|polycystic"),
    ("decidua",
     r"decidual|trophoblast|placenta|implantation|blastocyst|embryo"),
]

# Studies that are clearly outside the topic; listed explicitly so that the
# off_topic label is a decision rather than a default.
OFF_TOPIC = (r"microglia|APOE4|hippocamp|pan-cancer|natural killer cells|"
             r"neutrophil|leukemia|triple-negative|colorectal|glioma")

# Fields of the series record that are read for classification.
FIELD_TITLE = "Series_title"
FIELD_SUMMARY = "Series_summary"
FIELD_DESIGN = "Series_overall_design"
FIELD_TYPE = "Series_type"
FIELD_SAMPLES = "Series_sample_id"

# Truncation lengths for console output and for the stored title.
TITLE_STORED = 95
TITLE_PRINTED = 58
MAX_COLWIDTH = 76

# Classes that are printed group by group at the end, in this order.
REPORT_CLASSES = ("EC_tumor", "endometrium", "endometriosis", "organoid")


# --------------------------------------------------------------------------
# Retrieval and parsing
# --------------------------------------------------------------------------
def get(url: str, tries: int = REQUEST_TRIES) -> str:
    """Fetch a URL with retries, returning an empty string on failure.

    Parameters
    ----------
    url : absolute URL.
    tries : number of attempts.

    Returns
    -------
    Response body, or "" when every attempt failed.
    """
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=USER_AGENT)
            with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:                       # noqa: BLE001
            time.sleep(RETRY_BASE_SLEEP + RETRY_BASE_SLEEP * i)
    return ""


def series_fields(gse: str) -> Dict[str, str]:
    """Parse a GEO series record into a field -> value mapping.

    Parameters
    ----------
    gse : series accession.

    Returns
    -------
    dict keyed by the field name without its leading ``!``.

    Notes
    -----
    GEO wraps long fields over several lines and repeats some keys. Continuation
    lines (no ``!`` prefix) are appended to the field being read, and a repeated
    key is joined with ";" rather than overwriting, because dropping the second
    occurrence is how a multi-line summary silently loses half its text.
    """
    txt = get(f"{GEO_ACC_CGI}?acc={gse}&targ=self&form=text&view=brief")
    out: Dict[str, str] = {}
    current = None
    for line in txt.splitlines():
        match = re.match(r"!(\w+) = (.*)", line)
        if match:
            current = match.group(1)
            value = match.group(2).strip()
            out[current] = (out[current] + ";" + value) if current in out else value
        elif line.startswith("!") or not line.strip():
            current = None
        elif current:
            out[current] += " " + line.strip()
    return out


def classify(title: str, summary: str, design: str) -> str:
    """Assign one series to a topic class.

    Parameters
    ----------
    title, summary, design : text fields of the series record.

    Returns
    -------
    Class label; ``off_topic`` when nothing matches.
    """
    blob = f"{title} {summary} {design}".lower()
    for name, pattern in RULES:
        if re.search(pattern, blob):
            return name
    if re.search(OFF_TOPIC, blob):
        return "off_topic"
    return "off_topic"


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------
def curate(candidates: str) -> pd.DataFrame:
    """Classify every candidate series and write the curated table.

    Parameters
    ----------
    candidates : file name inside ``PATHS.tables``.

    Returns
    -------
    The curated table as written.

    Raises
    ------
    SystemExit
        If the candidate file is missing or lacks a required column.
    """
    PATHS.ensure()
    src = PATHS.tables / candidates
    if not src.exists():
        raise SystemExit(f"missing candidate list: {src}\n"
                         f"  produce one with 04_expression_matrix.py --step census, "
                         f"or pass --candidates")
    cand = pd.read_csv(src, sep="\t")
    for column in ("gse", "title"):
        if column not in cand.columns:
            raise SystemExit(f"{src.name}: no '{column}' column; "
                             f"found {list(cand.columns)}")

    rows = []
    for _, r in cand.iterrows():
        f = series_fields(r["gse"])
        # Samples are counted from the id list, not from a numeric field.
        n = len([x for x in f.get(FIELD_SAMPLES, "").split(";") if x.strip()])
        cls = classify(r["title"], f.get(FIELD_SUMMARY, ""), f.get(FIELD_DESIGN, ""))
        rows.append({"gse": r["gse"], "n": n, "cls": cls,
                     "type": f.get(FIELD_TYPE, ""),
                     "title": r["title"][:TITLE_STORED]})
        print(f"{r['gse']}  n={n:>5}  {cls:<14} {r['title'][:TITLE_PRINTED]}")
        time.sleep(SLEEP_PER_SERIES)

    d = pd.DataFrame(rows).sort_values(["cls", "n"], ascending=[True, False])
    dest = PATHS.tables / "w3_curated.tsv"
    d.to_csv(dest, sep="\t", index=False)

    print("\n=== classes ===")
    print(d.cls.value_counts().to_string())
    for c in REPORT_CLASSES:
        sub = d[d.cls == c]
        print(f"\n--- {c} ({len(sub)} series, {int(sub.n.sum())} samples) ---")
        print(sub[["gse", "n", "title"]].to_string(index=False,
                                                   max_colwidth=MAX_COLWIDTH))
    print(f"\nwritten {dest}")
    return d


def main() -> None:
    """Parse arguments and curate the candidate list."""
    ap = argparse.ArgumentParser(
        description="Classify candidate GEO series by topic.")
    ap.add_argument("--candidates", default="w3_candidates.tsv",
                    help="candidate list inside the tables directory")
    args = ap.parse_args()
    curate(args.candidates)


if __name__ == "__main__":
    main()
