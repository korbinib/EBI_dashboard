#!/usr/bin/env python3
"""
ebi_api.py
==========
Low-level client for the EBI Search REST API.

This module knows only about HTTP: making requests, retrying transient
failures, discovering retrievable fields, and probing hit counts.  It holds no
domain configuration and no Norwegian-filtering logic, so it imports nothing
from the other scripts (keeping the import graph acyclic).
"""

import logging

try:
    import requests
    from tenacity import (retry, stop_after_attempt, wait_exponential,
                          retry_if_exception_type)
    _REQUESTS_AVAILABLE = True
except ImportError:
    _REQUESTS_AVAILABLE = False

log = logging.getLogger("ebi_api")

# ── HTTP constants ──────────────────────────────────────────────────────────────
BASE_URL        = "https://www.ebi.ac.uk/ebisearch/ws/rest"
PAGE_SIZE       = 500
RATE_SLEEP      = 0.4
CATCH_ALL_QUERY = "*:*"
FALLBACK_FIELDS = ["name", "title", "description"]


# ── Request helper (with retry) ─────────────────────────────────────────────────
if _REQUESTS_AVAILABLE:
    SESSION = requests.Session()
    SESSION.headers.update({"Accept": "application/json"})

    @retry(
        stop=stop_after_attempt(5),
        wait=wait_exponential(multiplier=1, min=2, max=30),
        retry=retry_if_exception_type(requests.RequestException),
        reraise=True,
    )
    def get_json(url: str, params: dict) -> dict:
        resp = SESSION.get(url, params=params, timeout=60)
        resp.raise_for_status()
        return resp.json()
else:
    def get_json(url: str, params: dict) -> dict:  # type: ignore[misc]
        raise RuntimeError(
            "requests/tenacity not installed – cannot make HTTP requests. "
            "Run: pip install requests tenacity"
        )


# ── Field discovery ─────────────────────────────────────────────────────────────
def get_retrievable_fields(domain: str, cfg: dict) -> list[str]:
    """Merge a domain's required fields with all retrievable fields the API reports."""
    required: list[str] = cfg.get("required_fields", FALLBACK_FIELDS)
    url = f"{BASE_URL}/{domain}"
    try:
        data = get_json(url, {"format": "json"})
        discovered = [
            f["id"]
            for f in data.get("fieldInfos", [])
            if f.get("retrievable", False)
        ]
        merged = list(dict.fromkeys(required + discovered))
        log.info("  %s: %d required + %d discovered = %d total fields",
                 domain, len(required), len(discovered), len(merged))
        return merged
    except Exception as exc:
        log.warning("  %s: metadata fetch failed (%s) – using required fields only",
                    domain, exc)
        return required


# ── Hit-count probe ─────────────────────────────────────────────────────────────
def get_hit_count(domain: str, fields: list[str], query: str) -> int:
    """Return the reported number of hits for a query (size=1 probe)."""
    url = f"{BASE_URL}/{domain}"
    try:
        data = get_json(url, {
            "query":  query,
            "fields": fields[0] if fields else "id",
            "format": "json",
            "size":   1,
            "start":  0,
        })
        return int(data.get("hitCount", 0))
    except Exception as exc:
        log.warning("  hit-count probe failed for %s (%s)", domain, exc)
        return 0
