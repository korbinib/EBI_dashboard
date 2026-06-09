#!/usr/bin/env python3
"""
fetch_one_domain.py
===================
Fetches and filters EBI Search data for a SINGLE domain, then writes
the result to data/raw/<domain>/.

This script is designed to be called:
  - In PARALLEL by GitHub Actions (one matrix job per domain).
  - In SEQUENCE by fetch_ebi_data.py for local execution.

Usage
-----
  python scripts/fetch_one_domain.py <domain>

  where <domain> is one of the keys in DOMAINS, e.g.:
    python scripts/fetch_one_domain.py pride
    python scripts/fetch_one_domain.py sra-sample

All shared configuration (DOMAINS, filter helpers, HTTP helpers) is
imported directly from fetch_ebi_data.py so there is a single source
of truth for domain definitions and filtering logic.

Exit codes
----------
  0  – success (even if 0 Norwegian entries found)
  1  – unknown domain, or unrecoverable fetch error
"""

import json
import os
import sys
import logging
from pathlib import Path

# ── Ensure the repo root is on sys.path so we can import fetch_ebi_data ──────
_SCRIPTS_DIR = Path(__file__).resolve().parent
_REPO_ROOT   = _SCRIPTS_DIR.parent
sys.path.insert(0, str(_SCRIPTS_DIR))

# Change to repo root so relative paths (data/, etc.) resolve correctly
os.chdir(_REPO_ROOT)

from fetch_ebi_data import (
    DOMAINS,
    load_geo_tokens,
    load_institution_regexes,
    build_geo_regex,
    get_retrievable_fields,
    fetch_domain,
    save_domain,
    _partition_dir,
    log,
)


def load_domain_cfg(domain: str) -> dict:
    """
    Load the config for a single domain.

    The in-code DOMAINS dict (imported from fetch_ebi_data) is the single
    source of truth: it always reflects the running code version.  It is
    therefore preferred over data/domains.json, which is only a generated
    snapshot ("do not edit manually") that can drift — e.g. an older copy
    committed to git that lacks fields such as partition_date_field, which
    would silently disable partitioned fetching for large SRA domains.

    Preference order:
      1. DOMAINS dict       – authoritative; matches the running code.
      2. data/domains.json  – fallback only for a domain somehow absent
                              from DOMAINS.
    """
    if domain in DOMAINS:
        return DOMAINS[domain]

    domains_file = _REPO_ROOT / "data" / "domains.json"
    if domains_file.exists():
        try:
            data = json.loads(domains_file.read_text())
            cfg = data.get("domains", {}).get(domain)
            if cfg is not None:
                log.warning("Domain %r not in DOMAINS – falling back to data/domains.json",
                            domain)
                return cfg
        except Exception as exc:
            log.warning("Could not read data/domains.json (%s)", exc)

    return {}

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  [%(domain)s]  %(message)s"
    if False   # placeholder – real formatter set below
    else "%(asctime)s  %(levelname)-8s  %(message)s",
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

    cfg = load_domain_cfg(domain)

    log.info("=== fetch_one_domain: %s ===", domain)

    geo_tokens   = load_geo_tokens()
    inst_regexes = load_institution_regexes()
    geo_re       = build_geo_regex(geo_tokens)

    log.info("Filter ready: %d geo tokens, %d institution patterns",
             len(geo_tokens), len(inst_regexes))

    fields  = get_retrievable_fields(domain, cfg)
    entries = fetch_domain(domain, cfg, fields, geo_re, inst_regexes)
    save_domain(domain, entries, fields)

    # ── Partition checkpoint summary ─────────────────────────────────────────
    # If the domain used partitioned fetching, log a summary.
    # Checkpoint files are retained so that a re-run of this script for the
    # same domain resumes from where it left off rather than starting over.
    # To clear checkpoints after a successful run, uncomment the cleanup block.
    part_dir = _partition_dir(domain)
    if part_dir.exists():
        part_files = sorted(part_dir.glob("*.json"))
        log.info("  Partitions: %d checkpoint files retained in %s",
                 len(part_files), part_dir)
        # for pf in part_files:
        #     pf.unlink()
        # part_dir.rmdir()

    log.info("=== Done: %s (%d entries) ===", domain, len(entries))
    return 0


if __name__ == "__main__":
    sys.exit(main())
