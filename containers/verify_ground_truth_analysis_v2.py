#!/usr/bin/env python3
"""Fail closed unless the analysis-v2 Python runtime is complete and pinned."""

import sys

import duckdb
import pandas


EXPECTED_PYTHON = (3, 11)
EXPECTED_DUCKDB = "1.0.0"

if sys.version_info[:2] != EXPECTED_PYTHON:
    raise SystemExit(
        f"Expected Python {EXPECTED_PYTHON[0]}.{EXPECTED_PYTHON[1]}; "
        f"found {sys.version_info.major}.{sys.version_info.minor}"
    )
if duckdb.__version__ != EXPECTED_DUCKDB:
    raise SystemExit(
        f"Expected DuckDB {EXPECTED_DUCKDB}; found {duckdb.__version__}"
    )

print("Python environment verification: PASSED")
print("Python version:", sys.version.split()[0])
print("DuckDB version:", duckdb.__version__)
print("pandas version:", pandas.__version__)
