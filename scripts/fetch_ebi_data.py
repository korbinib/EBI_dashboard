#!/usr/bin/env python3
"""
fetch_ebi_data.py
=================
Queries the EBI Search REST API for all configured domains and saves
confirmed Norwegian entries to data/raw/<domain>/.

Strategy
--------
1. DISCOVER  – For each domain call GET /ws/rest/{domain} to retrieve the
               complete list of retrievable field IDs.

2. FETCH     – Send *:* (catch-all) with all retrievable fields requested.
               Large SRA domains (>1M entries) are split into year/quarter/month
               date-range windows.

3. INCREMENTAL CACHE – Only the last 2 calendar years are re-fetched on each
               run.  Windows for year ≤ (current_year - 2) are immutable: once
               a partition file exists on disk it is loaded from disk and never
               re-downloaded.  This applies to all domains that have a
               partition_date_field.

               Manifest files (data/raw/<domain>/manifest.json) record the
               sha256 of each immutable partition so corruption is detected
               and the affected window is re-fetched.

4. FILTER    – Concatenate every field value into a text blob and test against:
                 a) Geographic indicators  (Norway, Oslo, Bergen, Tromsø, \\bNO\\b …)
                 b) Institution regexes    (NTNU, Folkehelseinstituttet, …)
                 c) Norwegian TLD email    (@*.no)
               sra-study is saved unfiltered (730 K entries, filter deferred to
               join_ena.py so studies detectable only via joined experiment/sample
               signals are not lost).  sra-experiment and sra-sample are pre-
               filtered at page level — fetching 39–53 M rows unfiltered exceeds
               available RAM.

5. SAVE      – Write entries to data/raw/<domain>/latest.json (and a dated copy).

Output
------
  data/domains.json                      – domain config snapshot
  data/raw/<domain>/latest.json          – filtered entries for this run
  data/raw/<domain>/manifest.json        – incremental-cache manifest
  data/raw/<domain>/partitions/<key>.json – per-window checkpoint files
"""

import hashlib
import json
import os
import re
import time
import logging
from datetime import date
from pathlib import Path

try:
    import requests
    from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type
    _REQUESTS_AVAILABLE = True
except ImportError:
    _REQUESTS_AVAILABLE = False

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

BASE_URL   = "https://www.ebi.ac.uk/ebisearch/ws/rest"
PAGE_SIZE  = 500
RATE_SLEEP = 0.4
TODAY      = date.today().isoformat()

CATCH_ALL_QUERY     = "*:*"
# Hard pagination cap of the EBI Search API: it never returns more than this
# many entries for a single query, and it also CAPS the reported hitCount at
# this value.  A domain/window reporting exactly MAX_PAGEABLE may therefore hold
# far more — so any count that reaches the cap must be split into smaller
# (year → quarter → month) windows rather than paginated directly.
MAX_PAGEABLE        = 1_000_000
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
SRA_DOMAINS     = frozenset({"sra-study"})
FILTER_VERSION  = 2

FALLBACK_FIELDS = ["name", "title", "description"]


# ──────────────────────────────────────────────────────────────────────────────
# Norwegian filter helpers
# ──────────────────────────────────────────────────────────────────────────────

def load_geo_tokens(path: str = "data/institution_map.json") -> list[str]:
    try:
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
    except FileNotFoundError:
        log.warning("Institution map not found: %s – using fallback geo tokens", path)
        return ["Norway", "Norge", "Norwegian", "Norsk", "Oslo", "Bergen",
                "Trondheim", "Tromsø", "Stavanger"]
    result: list[str] = []
    for t in data.get("norway_indicators", []):
        t = t.strip()
        if not t or any(c in t for c in r"\()[]{}?*+^$|") or len(t) <= 2:
            continue
        result.append(t)
    return sorted(set(result))


def load_institution_regexes(path: str = "data/institution_map.json") -> list[re.Pattern]:
    try:
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
    except FileNotFoundError:
        log.warning("Institution map not found: %s – institution filter disabled", path)
        return []
    compiled: list[re.Pattern] = []
    for inst in data.get("institutions", []):
        for p in inst.get("patterns", []):
            try:
                compiled.append(re.compile(p, re.IGNORECASE))
            except re.error as exc:
                log.debug("Skipping invalid pattern %r: %s", p, exc)
    log.info("Loaded %d institution regex patterns", len(compiled))
    return compiled


def build_geo_regex(geo_tokens: list[str]) -> re.Pattern:
    parts = [re.escape(t) for t in geo_tokens] + [r"\bNO\b"]
    return re.compile("|".join(parts), re.IGNORECASE)


_EMAIL_NO_RE = re.compile(r"@[\w.\-]+\.no\b", re.IGNORECASE)


def is_norwegian_entry(entry: dict, geo_re: re.Pattern,
                       inst_regexes: list[re.Pattern]) -> bool:
    fields = entry.get("fields", {})
    parts: list[str] = []
    for vals in fields.values():
        if isinstance(vals, list):
            parts.extend(str(v) for v in vals if v is not None and str(v).strip())
        elif vals is not None and str(vals).strip():
            parts.append(str(vals))
    combined = " ".join(parts)
    if not combined.strip():
        return False
    if geo_re.search(combined):
        return True
    for pat in inst_regexes:
        if pat.search(combined):
            return True
    return bool(_EMAIL_NO_RE.search(combined))


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
    },
    # EGA (European Genome-phenome Archive): omitted because the EBI Search
    # index returns only id, description, and name — no date fields and no
    # affiliation fields are populated, making Norwegian detection impossible.
    "sra-study": {
        "required_fields": [
            "abstract", "acc", "alias", "center_project_name", "description",
            "domain_source", "first_public_date", "id", "insdc-project",
            "study_keywords", "study_type", "tag",
        ],
        # first_public_date is not searchable in sra-study; total ≈730K (<1M)
        # so standard single-pass pagination handles it without partitioning.
        "join_key": "study_accession",
    },
    "sra-sample": {
        "required_fields": [
            "acc", "alias", "broker_name", "center_name", "classification",
            "collection_date", "country", "description", "domain_source",
            "first_public_date", "host", "id", "isolate", "last_updated_date",
            "region", "sample_capture_status", "scientific_name", "strain",
            "submission_tool", "tag",
        ],
        # country:Norway filters at query time: ~150K entries vs 53M total.
        # This keeps the result well below MAX_PAGEABLE so no year-partitioning
        # is needed and no page-level regex filter is required.
        "query": "country:Norway",
        "join_key": "sample_accession",
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
            # XREF fields: these are the actual join keys in the EBI Search API.
            # study_accession and sample_accession are never populated; the
            # cross-domain links live here instead.
            "SRA-STUDY",   # ERP/SRP study accession → sra-study.id
            "SAMPLE",      # SAMEA/SAMN BioSample accession → sra-sample.id/acc
            "SRA-SAMPLE",  # ERS accession → sra-sample.id (ENA-native fallback)
        ],
        "partition_date_field": "first_public_date",
        "join_key": "SRA-STUDY",
    },
}


# ──────────────────────────────────────────────────────────────────────────────
# HTTP helpers
# ──────────────────────────────────────────────────────────────────────────────

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


# ──────────────────────────────────────────────────────────────────────────────
# Field discovery
# ──────────────────────────────────────────────────────────────────────────────

def get_retrievable_fields(domain: str, cfg: dict) -> list[str]:
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


# ──────────────────────────────────────────────────────────────────────────────
# Manifest (incremental cache bookkeeping)
# ──────────────────────────────────────────────────────────────────────────────

def _manifest_path(domain: str) -> Path:
    return Path("data/raw") / domain / "manifest.json"


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
    return Path("data/raw") / domain / "partitions"


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


def _get_hit_count(domain: str, fields: list[str], query: str) -> int:
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
            eid = e.get("id", "")
            if eid and eid in seen_ids:
                continue
            if eid:
                seen_ids.add(eid)
            all_entries.append(e)

    def _date_range(year: int, month: int | None,
                    quarter: int | None) -> tuple[str, str]:
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
                       month: int | None, indent: str = "  ") -> list[dict]:
        """
        Return entries for one date window, splitting recursively if needed.
        Immutable windows are served from disk; refetch windows are re-downloaded.
        """
        # Immutable window: serve from disk if file is OK
        if year <= immutable_max:
            if _partition_ok(domain, key, manifest):
                return _load_partition(domain, key) or []
            # File missing or corrupt → fall through to re-fetch

        # Refetch window: always invalidate first
        else:
            _invalidate_partition(domain, key, manifest)

        d_start, d_end = _date_range(year, month, quarter)
        window_query   = f"{date_field}:[{d_start} TO {d_end}]"
        count = _get_hit_count(domain, fields, window_query)
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

        # Window too large: split recursively
        if month is not None:
            log.warning(
                "%s%s %s has %d entries > MAX_PAGEABLE=%d; "
                "fetching up to %d — some entries will be missed.",
                indent, domain, d_start[:7], count, MAX_PAGEABLE, MAX_PAGEABLE,
            )
            entries = _fetch_window(domain, fields, window_query,
                                    geo_re, inst_regexes)
            _save_partition(domain, key, entries, manifest)
            time.sleep(RATE_SLEEP)
            return entries

        if quarter is not None:
            sub: list[dict] = []
            m0 = (quarter - 1) * 3 + 1
            for m in range(m0, m0 + 3):
                sub.extend(_fetch_or_load(
                    f"{year}_M{m:02d}", year, None, m, indent + "  "))
            _save_partition(domain, key, sub, manifest)
            return sub

        sub = []
        for q in range(1, 5):
            sub.extend(_fetch_or_load(
                f"{year}_Q{q}", year, q, None, indent + "  "))
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
    incremental cache.  Small domains (or domains with a narrow query) use
    standard single-pass pagination.

    If cfg["query"] is set it replaces the default catch-all (*:*) at the API
    level, pre-filtering results before any page-level Norwegian regex runs.
    """
    domain_query = cfg.get("query", CATCH_ALL_QUERY)

    if cfg.get("partition_date_field"):
        hit_count = _get_hit_count(domain, fields, domain_query)
        log.info("  %s → %d entries for query %r", domain, hit_count,
                 domain_query if domain_query != CATCH_ALL_QUERY else "*:*")
        # >= (not >): the API caps hitCount at MAX_PAGEABLE, so a domain reporting
        # exactly the cap may actually hold many more and must be partitioned.
        if hit_count >= MAX_PAGEABLE:
            log.info("  %s at/above MAX_PAGEABLE – using incremental partitioned fetch",
                     domain)
            return fetch_domain_partitioned(domain, cfg, fields, geo_re, inst_regexes)

    # Standard single-pass pagination
    url    = f"{BASE_URL}/{domain}"
    params = {
        "query":  domain_query,
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
    raw_dir = Path("data/raw") / domain
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


def save_domains_json(path: str = "data/domains.json") -> None:
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
# Main  –  sequential local orchestrator
# ──────────────────────────────────────────────────────────────────────────────

import subprocess
import sys


def main():
    if len(sys.argv) == 2 and sys.argv[1] == "--list-domains":
        save_domains_json()
        data = json.loads(Path("data/domains.json").read_text())
        print(json.dumps(list(data["domains"].keys())))
        return

    log.info("Starting EBI fetch (sequential / local)  date=%s", TODAY)
    save_domains_json()

    scripts_dir = Path(__file__).resolve().parent
    fetch_one   = scripts_dir / "fetch_one_domain.py"

    failed: list[str] = []
    for domain in DOMAINS:
        log.info("─── Domain: %s ───", domain)
        result = subprocess.run([sys.executable, str(fetch_one), domain], check=False)
        if result.returncode != 0:
            log.error("fetch_one_domain.py failed for %s (exit %d)",
                      domain, result.returncode)
            failed.append(domain)

    if failed:
        log.error("Failed domains: %s", ", ".join(failed))
    else:
        log.info("All domains fetched ✓")

    join_ena = scripts_dir / "join_ena.py"
    if join_ena.exists():
        log.info("─── Running join_ena.py ───")
        result = subprocess.run([sys.executable, str(join_ena)], check=False)
        if result.returncode != 0:
            log.error("join_ena.py failed (exit %d)", result.returncode)
        else:
            log.info("join_ena.py ✓")


if __name__ == "__main__":
    os.chdir(Path(__file__).parent.parent)
    main()
