"""Path configuration for the replication package.

Every script resolves the repository root from this file. Override the root by
setting the environment variable ECIR_ROOT, otherwise the directory containing
this file is used.
"""
import os
import pathlib

ROOT = pathlib.Path(os.environ.get("ECIR_ROOT", pathlib.Path(__file__).resolve().parent))
DATA = ROOT / "data"
RESULTS = ROOT / "results"
TABLES = RESULTS / "tables"
FIGURES = RESULTS / "figures"
CODE = ROOT / "code"
