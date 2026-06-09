#!/usr/bin/env python3
"""
fetch_one_domain.py
===================
Fetch + filter a SINGLE EBI Search domain and write data/raw/<domain>/.

This is a thin CLI wrapper around fetch_ebi_data.fetch_and_save_domain().  It is
the per-domain entry point invoked once per Snakemake job (one matrix job per
domain).  All shared logic — domain config, HTTP, filtering, fetch/cache — lives
in fetch_ebi_data.py and the modules it imports, so there is a single source of
truth.

Usage
-----
  python scripts/fetch_one_domain.py <domain>

  where <domain> is one of the keys in DOMAINS, e.g.:
    python scripts/fetch_one_domain.py pride
    python scripts/fetch_one_domain.py sra-sample

Exit codes
----------
  0  – success (even if 0 Norwegian entries found)
  1  – unknown domain, or unrecoverable fetch error
"""

import logging
import sys
from pathlib import Path

# Ensure scripts/ is importable regardless of how this file is invoked.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from fetch_ebi_data import DOMAINS, fetch_and_save_domain, log

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)


def main() -> int:
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <domain>", file=sys.stderr)
        print(f"Known domains: {', '.join(DOMAINS)}", file=sys.stderr)
        return 1

    domain = sys.argv[1].strip()
    if domain not in DOMAINS:
        print(f"Unknown domain {domain!r}. Known: {', '.join(DOMAINS)}",
              file=sys.stderr)
        return 1

    try:
        fetch_and_save_domain(domain)
    except Exception as exc:
        log.error("fetch failed for %s: %s", domain, exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
