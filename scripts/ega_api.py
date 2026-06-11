#!/usr/bin/env python3
"""
ega_api.py
==========
Low-level client for the EGA (European Genome-phenome Archive) Public Metadata
REST API — a DIFFERENT service from the EBI Search API used by every other
domain (see ebi_api.py).

Why a separate client
---------------------
The EBI Search index exposes EGA with only id/name/description populated — no
dates and no affiliation fields — so Norwegian detection is impossible there
(this is why "ega" is absent from the DOMAINS dict in fetch_ebi_data.py).  The
EGA Public Metadata API, by contrast, exposes the full object graph:

    DAC (EGAC, has contacts) ──< datasets (EGAD) ──< studies (EGAS, has dates)

This module knows only about HTTP: paginating list endpoints, fetching a single
object's sub-resources, and retrying transient failures.  It holds no
Norwegian-filtering logic and imports nothing from the other scripts.

Endpoints used by fetch_ega.py
------------------------------
  GET /dacs?limit=&offset=            → list DACs (accession_id, contacts, …)
  GET /dacs/{acc}/datasets            → datasets governed by a DAC
  GET /datasets/{acc}/studies         → studies a dataset belongs to (full objects)
  GET /datasets/{acc}/samples         → samples (EGAN) belonging to a dataset
"""

import logging

try:
    import requests
    from tenacity import (retry, stop_after_attempt, wait_exponential,
                          retry_if_exception_type)
    _REQUESTS_AVAILABLE = True
except ImportError:
    _REQUESTS_AVAILABLE = False

log = logging.getLogger("ega_api")

# ── HTTP constants ──────────────────────────────────────────────────────────────
BASE_URL   = "https://metadata.ega-archive.org"
PAGE_SIZE  = 500          # list endpoints accept limit/offset; 500 is comfortable
RATE_SLEEP = 0.2          # polite delay between paginated requests


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
    def get_json(url: str, params: dict | None = None) -> list | dict:
        resp = SESSION.get(url, params=params or {}, timeout=60)
        resp.raise_for_status()
        return resp.json()
else:
    def get_json(url: str, params: dict | None = None):  # type: ignore[misc]
        raise RuntimeError(
            "requests/tenacity not installed – cannot make HTTP requests. "
            "Run: pip install requests tenacity"
        )


# ── Pagination ──────────────────────────────────────────────────────────────────
def paginate(path: str, params: dict | None = None) -> list[dict]:
    """
    Page through a limit/offset list endpoint and return all rows.

    The EGA list endpoints return a bare JSON array; a short page (fewer than
    PAGE_SIZE rows) or an empty page signals the end.
    """
    import time

    url = f"{BASE_URL}/{path.lstrip('/')}"
    base_params = dict(params or {})
    rows: list[dict] = []
    offset = 0
    while True:
        page = get_json(url, {**base_params, "limit": PAGE_SIZE, "offset": offset})
        if not isinstance(page, list):
            log.warning("  %s: expected a list, got %s", url, type(page).__name__)
            break
        rows.extend(page)
        if len(page) < PAGE_SIZE:
            break
        offset += PAGE_SIZE
        time.sleep(RATE_SLEEP)
    return rows


def get_dacs() -> list[dict]:
    """Return every DAC the API exposes (each row already includes its contacts)."""
    return paginate("dacs")


def get_dac_datasets(dac_accession: str) -> list[dict]:
    """Return the datasets (EGAD) governed by a single DAC."""
    return paginate(f"dacs/{dac_accession}/datasets")


def get_dataset_studies(dataset_accession: str) -> list[dict]:
    """
    Return the studies (EGAS) a dataset belongs to, as full objects
    (accession_id, title, description, study_type, released_date, …).

    A dataset usually maps to one or a few studies, but we paginate anyway (the
    limit/offset params are harmless if the endpoint returns everything at once)
    so an unexpectedly large mapping is never silently truncated.
    """
    return paginate(f"datasets/{dataset_accession}/studies")


def get_dataset_samples(dataset_accession: str) -> list[dict]:
    """
    Return the samples (EGAN) belonging to a dataset.  Paginated, since a single
    dataset can hold thousands of samples (num_samples on the dataset object).
    Public sample fields are sparse (accession_id, biological_sex, phenotype, …)
    with no date or affiliation — those are inherited from the parent
    dataset/DAC by fetch_ega.py.
    """
    return paginate(f"datasets/{dataset_accession}/samples")
