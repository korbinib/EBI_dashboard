#!/usr/bin/env python3
"""
fetch_identifiers.py
====================
Cache the identifiers.org validation patterns for the *prefixes this project
actually uses* so the render step can turn raw accessions (PXD…, MTBLS…, EGAS…,
SAMEA…, …) into validated, resolvable links.

Rather than mirror the whole 860-namespace registry (most of which we never
touch), this fetches only the prefixes referenced by the domain definitions —
every `identifiers_prefix` in DOMAINS plus EXTRA_PREFIXES for domains handled
outside that dict (EGA studies).  The result is a ~1 KB file instead of ~108 KB.

What it writes
--------------
data/identifiers_namespaces.json — a compact map keyed by prefix:

    {
      "pride.project": {"prefix": "pride.project",
                        "pattern": "^P(X|R|A)D\\\\d{6}$",
                        "name": "PRIDE Project"},
      ...
    }

How it is used
--------------
R/plot_norwegian_data.R reads this file plus the per-domain `identifiers_prefix`
(from data/domains.json), checks each accession against the namespace `pattern`,
and builds `https://identifiers.org/<prefix>:<acc>` for the rows that match.

Source
------
identifiers.org Registry REST API — one lookup per prefix:
  GET https://registry.api.identifiers.org/restApi/namespaces/search/findByPrefix?prefix=<p>
See https://docs.identifiers.org/pages/api.html
"""

import json
import logging
import os
import sys
from pathlib import Path

# Ensure scripts/ is importable regardless of how this file is invoked.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from paths import IDENTIFIERS_JSON
from fetch_ebi_data import DOMAINS

try:
    import requests
    from tenacity import (retry, stop_after_attempt, wait_exponential,
                          retry_if_exception_type)
    _REQUESTS_AVAILABLE = True
except ImportError:
    _REQUESTS_AVAILABLE = False

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("fetch_identifiers")

FIND_BY_PREFIX = ("https://registry.api.identifiers.org/restApi/namespaces"
                  "/search/findByPrefix")

# Prefixes used by domains that aren't in the DOMAINS dict.  The EGA studies
# domain (handled by fetch_ega.py) resolves to ega.study; EGA samples (EGAN…)
# have no identifiers.org namespace, and ENA reuses sra-study's prefixes.
EXTRA_PREFIXES = ["ega.study"]


if _REQUESTS_AVAILABLE:
    @retry(
        stop=stop_after_attempt(5),
        wait=wait_exponential(multiplier=1, min=2, max=30),
        retry=retry_if_exception_type(requests.RequestException),
        reraise=True,
    )
    def _get(url: str, params: dict) -> dict:
        resp = requests.get(url, params=params, timeout=60,
                            headers={"Accept": "application/json"})
        resp.raise_for_status()
        return resp.json()
else:
    def _get(url: str, params: dict) -> dict:  # type: ignore[misc]
        raise RuntimeError("requests/tenacity not installed – cannot fetch registry")


def referenced_prefixes() -> list[str]:
    """Every identifiers.org prefix the project links against (DOMAINS + extras)."""
    prefixes: set[str] = set(EXTRA_PREFIXES)
    for cfg in DOMAINS.values():
        p = cfg.get("identifiers_prefix")
        if isinstance(p, str):
            prefixes.add(p)
        elif isinstance(p, (list, tuple)):
            prefixes.update(p)
    return sorted(prefixes)


def fetch_namespace(prefix: str) -> dict | None:
    """Look up one prefix; return {prefix, pattern, name} or None if absent."""
    try:
        data = _get(FIND_BY_PREFIX, {"prefix": prefix})
    except Exception as exc:
        log.warning("  findByPrefix(%s) failed: %s", prefix, exc)
        return None
    if not isinstance(data, dict) or not data.get("prefix"):
        log.warning("  prefix %r not found in registry", prefix)
        return None
    return {"prefix": data["prefix"],
            "pattern": data.get("pattern"),
            "name": data.get("name")}


def fetch_namespaces() -> dict[str, dict]:
    """Resolve every referenced prefix to its namespace record."""
    out: dict[str, dict] = {}
    for prefix in referenced_prefixes():
        ns = fetch_namespace(prefix)
        if ns:
            out[ns["prefix"]] = ns
    return out


def save(namespaces: dict[str, dict], path: Path = IDENTIFIERS_JSON) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"namespace_count": len(namespaces), "namespaces": namespaces}
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(payload, indent=2, ensure_ascii=False))
    os.replace(tmp, path)
    log.info("Wrote %s (%d namespaces)", path, len(namespaces))


def main() -> int:
    namespaces = fetch_namespaces()
    # If nothing came back (registry unreachable), keep any existing cache rather
    # than clobbering it — the link feature is optional, the pipeline goes on.
    if not namespaces and IDENTIFIERS_JSON.exists():
        log.warning("No namespaces fetched – keeping existing %s", IDENTIFIERS_JSON)
        return 0
    save(namespaces)
    return 0


if __name__ == "__main__":
    sys.exit(main())
