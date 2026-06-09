#!/usr/bin/env python3
"""
paths.py
========
Single source of truth for filesystem locations used across the fetch/join
scripts.  Every path is derived from this file's own location, so nothing
depends on the process working directory — importing any module no longer needs
an os.chdir() side effect to make relative "data/..." paths resolve.
"""

from pathlib import Path

REPO_ROOT       = Path(__file__).resolve().parent.parent
DATA_DIR        = REPO_ROOT / "data"
RAW_DIR         = DATA_DIR / "raw"
PROC_DIR        = DATA_DIR / "processed"
INSTITUTION_MAP = DATA_DIR / "institution_map.json"
DOMAINS_JSON    = DATA_DIR / "domains.json"
