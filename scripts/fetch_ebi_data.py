#!/usr/bin/env python3
"""
fetch_ebi_data.py
=================
Domain configuration and fetch/cache logic for the EBI Norway Dashboard.

Responsibilities
----------------
* DOMAINS         – the authoritative domain config (single source of truth;
                    imported by the Snakefile and fetch_one_domain.py).
* fetch_domain    – paginate a domain, routing large partitionable domains to
                    the incremental year/quarter/month fetch.
* incremental cache – per-window partition files + manifest with sha256 checks.
* save_domain     – write data/raw/<domain>/latest.json (+ dated copy).
* fetch_and_save_domain – per-domain orchestration used by both the CLI
                    (fetch_one_domain.py) and the local in-process run below.

The low-level HTTP client lives in ebi_api.py and the Norwegian filter in
norwegian_filter.py, so this module is free of those concerns.

Strategy
--------
1. DISCOVER  – ebi_api.get_retrievable_fields() lists retrievable field IDs.
2. FETCH     – *:* with all fields.  Domains with a partition_date_field that
               report ≥1M entries are split into year/quarter/month windows.
3. CACHE     – Only the last REFETCH_YEARS calendar years are re-fetched; older
               windows are served from sha256-verified partition files on disk.
4. FILTER    – sra-study is saved unfiltered (filter deferred to join_ena.py).
               sra-experiment and sra-sample are pre-filtered at page level
               (39–53 M rows unfiltered would exceed RAM).  Non-SRA domains are
               filtered at page level too.
5. SAVE      – data/raw/<domain>/latest.json (+ dated copy).

Output
------
  data/domains.json                       – domain config snapshot
  data/raw/<domain>/latest.json           – entries for this run
  data/raw/<domain>/manifest.json         – incremental-cache manifest
  data/raw/<domain>/partitions/<key>.json – per-window checkpoint files
"""

import hashlib
import json
import os
import sys
import logging
from datetime import date
from pathlib import Path

from paths import RAW_DIR, DOMAINS_JSON
from ebi_api import (
    BASE_URL, PAGE_SIZE, RATE_SLEEP, CATCH_ALL_QUERY,
    get_json, get_retrievable_fields, get_hit_count,
)
from norwegian_filter import (
    load_geo_tokens, load_institution_regexes, build_geo_regex,
    is_norwegian_entry,
)
import time
import re

# ──────────────────────────────────────────────────────────────────────────────
# Logging
# ──────────────────────────────────────────────────────────────────────────────

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("fetch_ebi")

# ──────────────────────────────────────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────────────────────────────────────

TODAY = date.today().isoformat()

# Hard pagination cap of the EBI Search API: it never returns more than this
# many entries for a single query, and it also CAPS the reported hitCount at
# this value.  A domain/window reporting exactly MAX_PAGEABLE may therefore hold
# far more — so any count that reaches the cap must be split into smaller
# (year → quarter → month) windows rather than paginated directly.
MAX_PAGEABLE         = 1_000_000
PARTITION_START_YEAR = 2007

# Number of trailing calendar years that are always re-fetched on every run.
# current_year and current_year-1 are refetched; everything older is immutable.
REFETCH_YEARS = 2

# sra-study: saved unfiltered so join_ena.py can find studies whose Norwegian
# signal is only visible after joining experiments or samples.  The domain
# has ~730 K entries total — well within the memory budget for a single pass.
#
# sra-experiment and sra-sample are NOT listed here: they have 39–53 M entries
# and are pre-filtered for Norwegian entries at page level inside _fetch_window.
# This trades a rare edge-case (Norwegian sample linked to a non-Norwegian
# experiment) for avoiding OOM on the fetch server.
#
# FILTER_VERSION must be bumped whenever the filtering strategy changes so that
# fetch_domain_partitioned() invalidates old unfiltered partition files.
SRA_DOMAINS    = frozenset({"sra-study"})
FILTER_VERSION = 2


# ──────────────────────────────────────────────────────────────────────────────
# Domain configuration
# ──────────────────────────────────────────────────────────────────────────────

DOMAINS: dict[str, dict] = {
    "bioimages": {
        "required_fields": [
            "acc", "attach_to", "author", "collection", "creation_date",
            "domain_source", "figure_sub", "figure_type", "id", "image_name",
            "journal_name", "legend", "license", "method", "modified_date",
            "name", "omics_type", "release_date", "repository", "species",
        ],
        # identifiers.org prefix(es) for this domain's accessions; consumed by
        # the R render step to build https://identifiers.org/<prefix>:<acc> links
        # (the accession is validated against the registry pattern first).
        "identifiers_prefix": "biostudies",
    },
    "biostudies-other": {
        "required_fields": [
            "abstract", "acc", "agency", "author", "collection",
            "creation_date", "data_source", "domain_source", "experiment_type",
            "grant_id", "id", "id_noversion", "journal", "method",
            "modified_date", "name", "omics_type", "organisation",
            "pagination", "pmcid", "project", "pub_date", "release_date",
            "repository", "species", "volume",
        ],
        "identifiers_prefix": "biostudies",
    },
    "metabolights": {
        "required_fields": [
            "description", "domain_source", "full_dataset_link",
            "ftp_download_link", "id", "instrument_platform", "name",
            "omics_type", "organism", "organism_group", "organism_part",
            "publication", "publication_date", "repository", "study",
            "study_design", "study_factor", "study_status", "submission_date",
            "submitter_affiliation", "submitter_email", "submitter_name",
            "technology_type",
        ],
        "identifiers_prefix": "metabolights",
    },
    "pride": {
        "required_fields": [
            "curator_keywords", "data_protocol", "description", "disease",
            "doi", "domain_source", "full_dataset_link", "id",
            "instrument_platform", "labhead", "labhead_affiliation",
            "labhead_mail", "modification", "name", "omics_type",
            "publication", "publication_date", "quantification_method",
            "repository", "sample_protocol", "software", "species",
            "submission_date", "submission_type", "submitter",
            "submitter_affiliation", "submitter_country", "submitter_keywords",
            "submitter_mail", "technology_type", "tissue",
        ],
        "identifiers_prefix": "pride.project",
    },
    "biomodels": {
        "required_fields": [
            "all_xrefs", "curationstatus", "description", "disease",
            "domain_source", "first_author", "full_dataset_link", "id",
            "isprivate", "last_modification_date", "levelversion", "modelflag",
            "modelformat", "modellingapproach", "name", "non_derived_xrefs",
            "omics_type", "publication", "publication_authors",
            "publication_date", "publication_doi", "publication_pubmed",
            "publication_title", "publication_url", "publication_year",
            "publicationid", "repository", "submission_date", "submissionid",
            "submitter", "submitter_affiliation", "submitter_keywords",
            "submitter_mail", "tokenised_name",
        ],
        "identifiers_prefix": "biomodels.db",
    },
    # EGA (European Genome-phenome Archive): NOT a DOMAINS entry because the EBI
    # Search index returns only id/description/name for it — no dates and no
    # affiliation fields — making Norwegian detection impossible here.  EGA is
    # instead fetched from the EGA Public Metadata API by scripts/fetch_ega.py
    # (DACs → datasets → studies) and written to data/raw/ega/latest.json.
    "sra-study": {
        "required_fields": [
            "abstract", "acc", "alias", "center_project_name", "description",
            "domain_source", "first_public_date", "id", "insdc-project",
            "study_keywords", "study_type", "tag",
        ],
        # first_public_date is not searchable in sra-study; total ≈730K (<1M)
        # so standard single-pass pagination handles it without partitioning.
        "join_key": "study_accession",
        # ENA/SRA study accessions are SRP/ERP/DRP (insdc.sra) or PRJ* (bioproject);
        # the render tries each prefix and links the one whose pattern matches.
        "identifiers_prefix": ["insdc.sra", "bioproject"],
    },
    "sra-sample": {
        "required_fields": [
            "acc", "alias", "broker_name", "center_name", "classification",
            "collection_date", "country", "description", "domain_source",
            "first_public_date", "host", "id", "isolate", "last_updated_date",
            "region", "sample_capture_status", "scientific_name", "strain",
            "submission_tool", "tag",
        ],
        "partition_date_field": "first_public_date",
        "join_key": "sample_accession",
        "identifiers_prefix": "biosample",
    },
    "sra-experiment": {
        "required_fields": [
            "abstract", "acc", "alias", "center_name", "classification",
            "collection_date", "country", "description", "domain_source",
            "first_public_date", "host", "id", "instrument_model",
            "instrument_platform", "last_updated_date", "library_layout",
            "library_name", "library_selection", "library_source",
            "library_strategy", "region", "scientific_name", "strain",
            "sub_species", "tag",
        ],
        "partition_date_field": "first_public_date",
        "join_key": "study_accession",
        "identifiers_prefix": ["insdc.sra", "bioproject"],
    },
}


# ──────────────────────────────────────────────────────────────────────────────
# Manifest (incremental cache bookkeeping)
# ──────────────────────────────────────────────────────────────────────────────

def _manifest_path(domain: str) -> Path:
    return RAW_DIR / domain / "manifest.json"


def _load_manifest(domain: str) -> dict:
    path = _manifest_path(domain)
    if path.exists():
        try:
            return json.loads(path.read_text())
        except Exception:
            pass
    return {"domain": domain, "partitions": {}}


def _save_manifest(domain: str, manifest: dict) -> None:
    path = _manifest_path(domain)
    path.parent.mkdir(parents=True, exist_ok=True)
    manifest["last_fetch_date"] = TODAY
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(manifest, indent=2, ensure_ascii=False))
    os.replace(tmp, path)


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


# ──────────────────────────────────────────────────────────────────────────────
# Partition helpers
# ──────────────────────────────────────────────────────────────────────────────

def _partition_dir(domain: str) -> Path:
    return RAW_DIR / domain / "partitions"


def _partition_path(domain: str, key: str) -> Path:
    return _partition_dir(domain) / f"{key}.json"


def _partition_ok(domain: str, key: str, manifest: dict) -> bool:
    """
    Return True if the partition file exists, is valid JSON, and its sha256
    matches the manifest record (if one exists).  On any mismatch returns False
    so the window is re-fetched.
    """
    path = _partition_path(domain, key)
    if not path.exists():
        return False
    try:
        json.loads(path.read_bytes())           # validity check
    except Exception:
        log.warning("  Corrupt partition %s – will re-fetch", path)
        return False
    recorded = manifest.get("partitions", {}).get(key, {}).get("sha256")
    if recorded:
        actual = _sha256_file(path)
        if actual != recorded:
            log.warning("  sha256 mismatch for partition %s – will re-fetch", key)
            return False
    return True


def _invalidate_partition(domain: str, key: str, manifest: dict) -> None:
    """Delete a partition file and remove it from the manifest."""
    path = _partition_path(domain, key)
    if path.exists():
        path.unlink()
    manifest.get("partitions", {}).pop(key, None)
    # Also delete sub-window files (Q and M keys) under the same year
    year = key.split("_")[0]
    part_dir = _partition_dir(domain)
    for child in list(part_dir.glob(f"{year}_*.json")):
        child.unlink()
        child_key = child.stem
        manifest.get("partitions", {}).pop(child_key, None)
    log.info("  Invalidated partition %s and sub-windows for year %s", key, year)


def _load_partition(domain: str, key: str) -> list[dict] | None:
    path = _partition_path(domain, key)
    if not path.exists():
        return None
    try:
        entries = json.loads(path.read_text())
        log.info("    ↩ loaded %s partition %s from disk (%d entries)",
                 domain, key, len(entries))
        return entries
    except Exception as exc:
        log.warning("    Corrupt partition %s – will re-fetch (%s)", path, exc)
        return None


def _save_partition(domain: str, key: str, entries: list[dict],
                    manifest: dict) -> None:
    """Write partition atomically and record its sha256 in the manifest."""
    path = _partition_path(domain, key)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = json.dumps(entries, ensure_ascii=False)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(raw)
    os.replace(tmp, path)
    manifest.setdefault("partitions", {})[key] = {
        "entries":    len(entries),
        "sha256":     _sha256_file(path),
        "fetch_date": TODAY,
    }


def _fetch_window(domain: str, fields: list[str], query: str,
                  geo_re: re.Pattern, inst_regexes: list[re.Pattern]) -> list[dict]:
    """Paginate through a single query window."""
    url    = f"{BASE_URL}/{domain}"
    params = {
        "query":  query,
        "fields": ",".join(fields),
        "format": "json",
        "size":   PAGE_SIZE,
        "start":  0,
    }
    entries:   list[dict] = []
    hit_count: int | None = None

    while True:
        try:
            data = get_json(url, params)
        except Exception as exc:
            log.error("    Window fetch failed at start=%d: %s", params["start"], exc)
            break

        if hit_count is None:
            hit_count = int(data.get("hitCount", 0))

        batch = data.get("entries", [])
        if domain in SRA_DOMAINS:
            entries.extend(batch)
        else:
            entries.extend(
                e for e in batch if is_norwegian_entry(e, geo_re, inst_regexes)
            )

        params["start"] += PAGE_SIZE
        if params["start"] >= (hit_count or 0) or not batch:
            break
        time.sleep(RATE_SLEEP)

    return entries


# ──────────────────────────────────────────────────────────────────────────────
# Incremental partitioned fetch
# ──────────────────────────────────────────────────────────────────────────────

def fetch_domain_partitioned(domain: str, cfg: dict, fields: list[str],
                              geo_re: re.Pattern,
                              inst_regexes: list[re.Pattern]) -> list[dict]:
    """
    Fetch a partitioned domain with incremental caching.

    Incremental rule
    ----------------
    immutable : year <= current_year - REFETCH_YEARS
        Loaded from disk if the partition file passes its sha256 check.
        Never re-fetched unless the file is absent or corrupt.

    refetch   : year >= current_year - REFETCH_YEARS + 1
        Deleted and re-fetched on every run (these windows grow as new
        records are added to EBI).

    Window splitting (recursive, year → quarter → month)
    --------------------------------------------------------
    If a window exceeds MAX_PAGEABLE it is split into sub-windows.
    Sub-windows follow the same immutable/refetch rule based on their year.
    """
    import calendar

    date_field   = cfg.get("partition_date_field", "first_public_date")
    current_year = date.today().year
    immutable_max = current_year - REFETCH_YEARS   # years ≤ this are immutable

    manifest = _load_manifest(domain)

    # If the filter strategy changed (e.g. domain moved from unfiltered to
    # pre-filtered), old partition files contain the wrong content.  Wipe them
    # so they are re-fetched with the current strategy.
    if manifest.get("filter_version", 1) != FILTER_VERSION:
        log.info("  %s: filter_version changed (%s→%d) – clearing cached partitions",
                 domain, manifest.get("filter_version", 1), FILTER_VERSION)
        part_dir = _partition_dir(domain)
        if part_dir.exists():
            for pf in sorted(part_dir.glob("*.json")):
                pf.unlink()
        manifest = {"domain": domain, "partitions": {}}
    manifest["filter_version"] = FILTER_VERSION

    all_entries: list[dict] = []
    seen_ids:    set[str]   = set()

    def _add(batch: list[dict]) -> None:
        for e in batch:
            # Entries with an id de-dup by id; id-less entries (rare) de-dup by
            # their content so the same record appearing in overlapping or
            # retried windows is not counted twice.
            eid = e.get("id", "")
            key = eid or json.dumps(e.get("fields", {}), sort_keys=True,
                                    ensure_ascii=False)
            if key in seen_ids:
                continue
            seen_ids.add(key)
            all_entries.append(e)

    def _date_range(year: int, month: int | None, quarter: int | None,
                    day: int | None = None) -> tuple[str, str]:
        if day is not None and month is not None:
            return (f"{year}-{month:02d}-{day:02d}",
                    f"{year}-{month:02d}-{day:02d}")
        if month is not None:
            last_day = calendar.monthrange(year, month)[1]
            return (f"{year}-{month:02d}-01", f"{year}-{month:02d}-{last_day:02d}")
        if quarter is not None:
            m_start = (quarter - 1) * 3 + 1
            m_end   = quarter * 3
            last_day = calendar.monthrange(year, m_end)[1]
            return (f"{year}-{m_start:02d}-01", f"{year}-{m_end:02d}-{last_day:02d}")
        return (f"{year}-01-01", f"{year}-12-31")

    def _fetch_or_load(key: str, year: int, quarter: int | None,
                       month: int | None, day: int | None = None,
                       indent: str = "  ") -> list[dict]:
        """
        Return entries for one date window, splitting recursively if needed
        (year → quarter → month → day).  Immutable windows are served from disk;
        refetch windows are re-downloaded.
        """
        # Immutable window: serve from disk if file is OK
        if year <= immutable_max:
            if _partition_ok(domain, key, manifest):
                return _load_partition(domain, key) or []
            # File missing or corrupt → fall through to re-fetch

        # Refetch window: always invalidate first
        else:
            _invalidate_partition(domain, key, manifest)

        d_start, d_end = _date_range(year, month, quarter, day)
        window_query   = f"{date_field}:[{d_start} TO {d_end}]"
        count = get_hit_count(domain, fields, window_query)
        log.info("%s%s  %s–%s  %d entries", indent, domain, d_start, d_end, count)

        if count == 0:
            _save_partition(domain, key, [], manifest)
            return []

        # < (not <=): a window reporting exactly the cap may hold more than it
        # reports, so only paginate directly when strictly below the cap.
        if count < MAX_PAGEABLE:
            entries = _fetch_window(domain, fields, window_query,
                                    geo_re, inst_regexes)
            _save_partition(domain, key, entries, manifest)
            log.info("%s→ fetched %d entries", indent + "  ", len(entries))
            time.sleep(RATE_SLEEP)
            return entries

        # Window too large: split to the next finer granularity.
        if day is not None:
            # A single day over the cap is unsplittable (the API has no sub-day
            # date field), so fetch up to the cap and warn.  Vanishingly rare.
            log.warning(
                "%s%s %s has %d entries > MAX_PAGEABLE=%d; fetching up to %d — "
                "some entries will be missed.",
                indent, domain, d_start, count, MAX_PAGEABLE, MAX_PAGEABLE,
            )
            entries = _fetch_window(domain, fields, window_query,
                                    geo_re, inst_regexes)
            _save_partition(domain, key, entries, manifest)
            time.sleep(RATE_SLEEP)
            return entries

        if month is not None:
            sub: list[dict] = []
            last_day = calendar.monthrange(year, month)[1]
            for d in range(1, last_day + 1):
                sub.extend(_fetch_or_load(
                    f"{year}_M{month:02d}_D{d:02d}", year, None, month, d,
                    indent + "  "))
            _save_partition(domain, key, sub, manifest)
            return sub

        if quarter is not None:
            sub = []
            m0 = (quarter - 1) * 3 + 1
            for m in range(m0, m0 + 3):
                sub.extend(_fetch_or_load(
                    f"{year}_M{m:02d}", year, None, m, None, indent + "  "))
            _save_partition(domain, key, sub, manifest)
            return sub

        sub = []
        for q in range(1, 5):
            sub.extend(_fetch_or_load(
                f"{year}_Q{q}", year, q, None, None, indent + "  "))
        _save_partition(domain, key, sub, manifest)
        return sub

    log.info("  %s: incremental partitioned fetch  date_field=%s"
             "  years=%d–%d  immutable_up_to=%d",
             domain, date_field, PARTITION_START_YEAR,
             current_year, immutable_max)

    for year in range(PARTITION_START_YEAR, current_year + 1):
        year_entries = _fetch_or_load(str(year), year, None, None)
        _add(year_entries)

    _save_manifest(domain, manifest)
    log.info("  %s: done – %d unique entries (%d immutable years cached)",
             domain, len(all_entries), max(0, immutable_max - PARTITION_START_YEAR + 1))
    return all_entries


# ──────────────────────────────────────────────────────────────────────────────
# Core fetch dispatcher
# ──────────────────────────────────────────────────────────────────────────────

def fetch_domain(domain: str, cfg: dict, fields: list[str],
                 geo_re: re.Pattern, inst_regexes: list[re.Pattern]) -> list[dict]:
    """
    Fetch all entries for a domain.

    Domains with partition_date_field are routed to fetch_domain_partitioned()
    which implements both the year/quarter/month window splitting and the
    incremental cache.  Small domains use standard single-pass pagination.
    """
    if cfg.get("partition_date_field"):
        hit_count = get_hit_count(domain, fields, CATCH_ALL_QUERY)
        log.info("  %s → %d total entries in domain", domain, hit_count)
        # >= (not >): the API caps hitCount at MAX_PAGEABLE, so a domain reporting
        # exactly the cap may actually hold many more and must be partitioned.
        if hit_count >= MAX_PAGEABLE:
            log.info("  %s at/above MAX_PAGEABLE – using incremental partitioned fetch",
                     domain)
            return fetch_domain_partitioned(domain, cfg, fields, geo_re, inst_regexes)

    # Standard single-pass pagination for small domains
    url    = f"{BASE_URL}/{domain}"
    params = {
        "query":  CATCH_ALL_QUERY,
        "fields": ",".join(fields),
        "format": "json",
        "size":   PAGE_SIZE,
        "start":  0,
    }
    entries:    list[dict] = []
    hit_count:  int | None = None
    total_seen: int        = 0

    while True:
        log.info("  GET %s  start=%d", domain, params["start"])
        try:
            data = get_json(url, params)
        except Exception as exc:
            log.error("  Failed fetching %s at start=%d: %s",
                      domain, params["start"], exc)
            break

        if hit_count is None:
            hit_count = int(data.get("hitCount", 0))
            log.info("  %s → %d total entries in domain", domain, hit_count)

        batch = data.get("entries", [])
        total_seen += len(batch)

        if domain in SRA_DOMAINS:
            entries.extend(batch)
        else:
            entries.extend(
                e for e in batch if is_norwegian_entry(e, geo_re, inst_regexes)
            )

        params["start"] += PAGE_SIZE
        if params["start"] >= (hit_count or 0) or not batch:
            break
        time.sleep(RATE_SLEEP)

    if domain in SRA_DOMAINS:
        log.info("  %s → saved all %d entries unfiltered (filter deferred to join_ena.py)",
                 domain, len(entries))
    else:
        log.info("  %s → kept %d / %d after Norwegian filter",
                 domain, len(entries), total_seen)
    return entries


# ──────────────────────────────────────────────────────────────────────────────
# Save
# ──────────────────────────────────────────────────────────────────────────────

def save_domain(domain: str, entries: list[dict], fields: list[str]) -> Path:
    """Write entries to data/raw/<domain>/latest.json (and a dated copy)."""
    raw_dir = RAW_DIR / domain
    raw_dir.mkdir(parents=True, exist_ok=True)

    payload = {
        "domain":      domain,
        "fetch_date":  TODAY,
        "query":       CATCH_ALL_QUERY,
        "fields_used": fields,
        "entry_count": len(entries),
        "entries":     entries,
    }
    body = json.dumps(payload, indent=2, ensure_ascii=False)

    dated_path = raw_dir / f"{TODAY}.json"
    dated_path.write_text(body)
    log.info("  Saved %d entries → %s", len(entries), dated_path)

    latest_path = raw_dir / "latest.json"
    latest_path.write_text(body)
    return dated_path


def save_domains_json(path: Path = DOMAINS_JSON) -> None:
    """Write DOMAINS config to data/domains.json for external consumers."""
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "_comment": (
            "Auto-generated by fetch_ebi_data.py – do not edit manually. "
            "Edit DOMAINS in scripts/fetch_ebi_data.py instead."
        ),
        "generated":    TODAY,
        "domain_count": len(DOMAINS),
        "domains":      DOMAINS,
    }
    tmp = out_path.with_suffix(".tmp")
    tmp.write_text(json.dumps(payload, indent=2, ensure_ascii=False))
    os.replace(tmp, out_path)
    log.info("Wrote %s (%d domains)", path, len(DOMAINS))


# ──────────────────────────────────────────────────────────────────────────────
# Per-domain orchestration  (single source of truth for one domain's fetch)
# ──────────────────────────────────────────────────────────────────────────────

def fetch_and_save_domain(domain: str,
                          geo_re: re.Pattern | None = None,
                          inst_regexes: list[re.Pattern] | None = None) -> int:
    """
    Fetch + filter + save a single domain.  Used by fetch_one_domain.py (one
    domain per Snakemake job) and by main()'s local in-process loop.

    The Norwegian filter is built on demand; callers running many domains can
    pass a prebuilt geo_re / inst_regexes to avoid recompiling per domain.
    Returns the number of entries saved.
    """
    cfg = DOMAINS[domain]

    if geo_re is None or inst_regexes is None:
        geo_tokens   = load_geo_tokens()
        inst_regexes = load_institution_regexes()
        geo_re       = build_geo_regex(geo_tokens)
        log.info("Filter ready: %d geo tokens, %d institution patterns",
                 len(geo_tokens), len(inst_regexes))

    log.info("=== fetch_and_save_domain: %s ===", domain)
    fields  = get_retrievable_fields(domain, cfg)
    entries = fetch_domain(domain, cfg, fields, geo_re, inst_regexes)
    save_domain(domain, entries, fields)

    # Partition checkpoint summary (files are retained so a re-run resumes)
    part_dir = _partition_dir(domain)
    if part_dir.exists():
        part_files = sorted(part_dir.glob("*.json"))
        log.info("  Partitions: %d checkpoint files retained in %s",
                 len(part_files), part_dir)

    log.info("=== Done: %s (%d entries) ===", domain, len(entries))
    return len(entries)


# ──────────────────────────────────────────────────────────────────────────────
# Main  –  sequential local orchestrator (Snakemake is the production driver)
# ──────────────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) == 2 and sys.argv[1] == "--list-domains":
        save_domains_json()
        print(json.dumps(list(DOMAINS.keys())))
        return

    log.info("Starting EBI fetch (sequential / local)  date=%s", TODAY)
    save_domains_json()

    # Cache the identifiers.org namespace registry for the link feature.  Failure
    # is non-fatal (the render step simply skips links it can't build).
    log.info("─── Caching identifiers.org namespaces ───")
    try:
        import fetch_identifiers
        fetch_identifiers.main()
    except Exception as exc:
        log.error("fetch_identifiers failed: %s", exc)

    # Build the Norwegian filter once and reuse it across all domains.
    geo_tokens   = load_geo_tokens()
    inst_regexes = load_institution_regexes()
    geo_re       = build_geo_regex(geo_tokens)
    log.info("Filter ready: %d geo tokens, %d institution patterns",
             len(geo_tokens), len(inst_regexes))

    failed: list[str] = []
    for domain in DOMAINS:
        log.info("─── Domain: %s ───", domain)
        try:
            fetch_and_save_domain(domain, geo_re, inst_regexes)
        except Exception as exc:
            log.error("fetch_and_save_domain failed for %s: %s", domain, exc)
            failed.append(domain)

    if failed:
        log.error("Failed domains: %s", ", ".join(failed))
    else:
        log.info("All domains fetched ✓")

    # Run the SRA join in-process (no subprocess indirection).
    log.info("─── Running join_ena ───")
    try:
        import join_ena
        join_ena.main()
        log.info("join_ena ✓")
    except Exception as exc:
        log.error("join_ena failed: %s", exc)

    # Fetch EGA studies via the EGA Public Metadata API (separate service; not a
    # DOMAINS entry — see scripts/fetch_ega.py for why).
    log.info("─── Running fetch_ega ───")
    try:
        import fetch_ega
        fetch_ega.main()
        log.info("fetch_ega ✓")
    except Exception as exc:
        log.error("fetch_ega failed: %s", exc)


if __name__ == "__main__":
    main()
