#!/usr/bin/env python3
"""
join_ena.py
===========
Joins the three SRA/ENA sub-tables fetched by fetch_ebi_data.py into a single
flat table `data/processed/ena_joined.json`, filtered to Norwegian entries.

Sub-tables used
---------------
sra-study      (primary key: study_accession)
sra-experiment (study_accession → study)
sra-sample     (sample_accession, via experiment)

Dropped vs. earlier version
----------------------------
sra-submission  No searchable date field; 31M entries not partitionable.
sra-analysis    Date field mismatch (last_updated_date ≠ first_public_date).
sra-run         No searchable date field; 42M entries; only used for run counts.

Filtering strategy
------------------
sra-study      Saved unfiltered by the fetch step (~730 K entries); the
               Norwegian filter is applied here post-join so that studies
               detectable only via joined experiment/sample signals are kept.

sra-experiment Pre-filtered for Norwegian entries at page level during the
sra-sample     fetch step (memory constraint: 39–53 M rows each).  The
               partition cache is versioned (FILTER_VERSION=2) so stale
               unfiltered partition files from older runs are discarded.

Recovery for uncovered samples
------------------------------
After the initial join, any Norwegian sample whose linked experiment was
filtered out (non-Norwegian center_name/country) is detected and a targeted
API query fetches the missing experiment→study links from sra-experiment.
This uses the SAMPLE XREF field (SAMEA/SAMN → ERS fallback) in batches of
50.  Network access is required for this step; failures are logged and the
join continues with the links already available.

EBI API join-key note: study_accession and sample_accession fields in
sra-experiment are always empty.  The actual links are in the XREF fields
SRA-STUDY (ERP/SRP accession) and SAMPLE (SAMEA/SAMN BioSample accession).

Output
------
  data/processed/ena_joined.json
  data/processed/ena_joined_<date>.json
"""

import json
import logging
import re
import time
import sys
from pathlib import Path
from datetime import date
from typing import Iterator

import pandas as pd

# Ensure scripts/ is importable regardless of how this file is invoked.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from norwegian_filter import (
    load_geo_tokens, load_institution_regexes, build_geo_regex,
    make_combined_filter,
)
from paths import RAW_DIR, PROC_DIR

try:
    import requests as _requests
    _REQUESTS_AVAILABLE = True
except ImportError:
    _REQUESTS_AVAILABLE = False

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("join_ena")

TODAY = date.today().isoformat()
PROC_DIR.mkdir(parents=True, exist_ok=True)


# ──────────────────────────────────────────────────────────────────────────────
# Projected loading helpers
# ──────────────────────────────────────────────────────────────────────────────

def _fv(fields: dict, key: str) -> str:
    v = fields.get(key, [])
    if isinstance(v, list):
        return v[0] if v else ""
    return str(v) if v else ""


def _fvlist(fields: dict, key: str) -> list[str]:
    v = fields.get(key, [])
    items = v if isinstance(v, list) else [v]
    return [str(x) for x in items if x]


def iter_entries(domain: str) -> Iterator[dict]:
    path = RAW_DIR / domain / "latest.json"
    if not path.exists():
        log.warning("No latest.json for %s – skipping", domain)
        return
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    for entry in data.get("entries", []):
        yield entry


def load_studies() -> pd.DataFrame:
    """
    Load all sra-study entries (unfiltered).  Norwegian filter is applied post-join
    so that studies detectable only via their joined experiment signals are kept.
    """
    rows = []
    for e in iter_entries("sra-study"):
        f = e.get("fields", {})
        rows.append({
            "study_acc":   e.get("id", "") or _fv(f, "acc"),
            "title":       _fv(f, "abstract") or _fv(f, "description"),
            "center_name": _fv(f, "center_project_name"),
            "description": _fv(f, "description"),
            "study_text":  " ".join(filter(None, [
                _fv(f, "abstract"), _fv(f, "description"),
                _fv(f, "center_project_name"), _fv(f, "alias"),
                _fv(f, "study_keywords"), _fv(f, "study_type"),
            ])),
        })
    df = pd.DataFrame(rows, dtype="string")
    log.info("  sra-study:      %d rows loaded", len(df))
    return df


def load_experiments() -> pd.DataFrame:
    """Load pre-filtered Norwegian sra-experiment entries; project to join keys + signal columns."""
    rows = []
    for e in iter_entries("sra-experiment"):
        f = e.get("fields", {})
        # SRA-STUDY and SAMPLE are the actual join-key XREF fields in the EBI Search API.
        # study_accession / sample_accession are never populated.
        rows.append({
            "exp_acc":           e.get("id", "") or _fv(f, "acc"),
            "study_acc":         _fv(f, "SRA-STUDY"),
            "sample_acc":        _fv(f, "SAMPLE") or _fv(f, "SRA-SAMPLE"),
            "first_public_date": _fv(f, "first_public_date"),
            "exp_country":       _fv(f, "country"),
            "exp_center":        _fv(f, "center_name"),
            "exp_text": " ".join(filter(None, [
                _fv(f, "abstract"),
                _fv(f, "alias"),
                _fv(f, "country"),
                _fv(f, "center_name"),
                _fv(f, "description"),
                _fv(f, "region"),
            ])),
        })
    df = pd.DataFrame(rows, dtype="string")
    log.info("  sra-experiment: %d rows loaded", len(df))
    return df


def load_samples() -> pd.DataFrame:
    """Load pre-filtered Norwegian sra-sample entries; project to acc + signal columns."""
    rows = []
    for e in iter_entries("sra-sample"):
        f = e.get("fields", {})
        rows.append({
            "sample_acc":     e.get("id", "") or _fv(f, "acc"),
            "sample_country": _fv(f, "country"),
            "sample_center":  _fv(f, "center_name"),
            "sample_region":  _fv(f, "region"),
            # broker_name, alias, description are fetched but were previously
            # dropped before the Norwegian filter; include them now.
            "sample_text":    " ".join(filter(None, [
                _fv(f, "broker_name"),
                _fv(f, "alias"),
                _fv(f, "description"),
            ])),
        })
    df = pd.DataFrame(rows, dtype="string")
    log.info("  sra-sample:     %d rows loaded", len(df))
    return df


# ──────────────────────────────────────────────────────────────────────────────
# Targeted experiment-link recovery for uncovered Norwegian samples
# ──────────────────────────────────────────────────────────────────────────────

_EBI_URL = "https://www.ebi.ac.uk/ebisearch/ws/rest"
_LINK_BATCH = 50   # Lucene OR-clause limit for XREF queries


def _fetch_experiment_links(sample_accs: list[str]) -> list[dict]:
    """
    For Norwegian samples not covered by any experiment in df_exps (i.e. their
    experiment was filtered out because it has no Norwegian center/country),
    fetch the minimal (exp_acc, study_acc, sample_acc) link rows directly from
    the sra-experiment API using SAMPLE XREF queries.

    Returns a list of bare link rows (no Norwegian-signal fields populated)
    suitable for appending to df_exps before the join aggregation.
    """
    if not _REQUESTS_AVAILABLE:
        log.warning("requests not installed – cannot recover uncovered sample links")
        return []

    rows: list[dict] = []
    for i in range(0, len(sample_accs), _LINK_BATCH):
        batch = sample_accs[i : i + _LINK_BATCH]
        query = " OR ".join(f"SAMPLE:{acc}" for acc in batch)
        start = 0
        while True:
            try:
                resp = _requests.get(
                    f"{_EBI_URL}/sra-experiment",
                    params={
                        "query":  query,
                        "fields": "acc,SRA-STUDY,SAMPLE,SRA-SAMPLE",
                        "format": "json",
                        "size":   500,
                        "start":  start,
                    },
                    timeout=60,
                )
                resp.raise_for_status()
                data = resp.json()
            except Exception as exc:
                log.warning("    Link-fetch batch %d failed: %s", i // _LINK_BATCH, exc)
                break

            for entry in data.get("entries", []):
                f = entry.get("fields", {})
                rows.append({
                    "exp_acc":           entry.get("id", ""),
                    "study_acc":         _fv(f, "SRA-STUDY"),
                    "sample_acc":        _fv(f, "SAMPLE") or _fv(f, "SRA-SAMPLE"),
                    "first_public_date": "",
                    "exp_country":       "",
                    "exp_center":        "",
                    "exp_text":          "",
                })

            hit_count = data.get("hitCount", 0)
            start += 500
            if start >= hit_count or not data.get("entries"):
                break
            time.sleep(0.4)

    return rows


# ──────────────────────────────────────────────────────────────────────────────
# Main join
# ──────────────────────────────────────────────────────────────────────────────

def main():
    geo_tokens   = load_geo_tokens()
    inst_regexes = load_institution_regexes()
    geo_re       = build_geo_regex(geo_tokens)
    combined_re  = make_combined_filter(geo_re, inst_regexes)
    log.info("Filter: %d geo tokens, %d institution patterns",
             len(geo_tokens), len(inst_regexes))

    log.info("Loading SRA sub-tables …")
    df_studies  = load_studies()
    df_exps     = load_experiments()
    df_samples  = load_samples()

    # ── Recover Norwegian samples whose experiment was filtered out ────────────
    # sra-experiment is pre-filtered for Norwegian entries, so samples linked
    # via non-Norwegian experiments have no path to their study.  Find those
    # samples and fetch the missing experiment→study link rows from the API.
    covered_sample_accs = set(df_exps["sample_acc"].dropna())
    covered_sample_accs.discard("")
    uncovered = [
        acc for acc in df_samples["sample_acc"].dropna()
        if acc and acc not in covered_sample_accs
    ]
    if uncovered:
        log.info("  %d Norwegian samples lack a Norwegian experiment – fetching links …",
                 len(uncovered))
        link_rows = _fetch_experiment_links(uncovered)
        if link_rows:
            df_links = pd.DataFrame(link_rows, dtype="string")
            df_exps = pd.concat([df_exps, df_links], ignore_index=True).drop_duplicates(
                subset=["exp_acc"]
            )
            log.info("  df_exps after link recovery: %d rows", len(df_exps))
    else:
        log.info("  All Norwegian samples covered by a Norwegian experiment ✓")

    # ── Join experiments → samples ────────────────────────────────────────────
    df_exp_sample = df_exps.merge(df_samples, on="sample_acc", how="left")

    # ── Aggregate experiment+sample signals per study ─────────────────────────
    def join_unique(series: pd.Series) -> str:
        vals = series.dropna()
        vals = vals[vals != ""]
        return " | ".join(sorted(set(vals)))

    exp_agg = df_exp_sample.groupby("study_acc", as_index=False).agg(
        n_experiments     = ("exp_acc",             "nunique"),
        first_public_date = ("first_public_date",   lambda s: min((v for v in s if v), default="")),
        exp_countries     = ("exp_country",         join_unique),
        exp_centers       = ("exp_center",          join_unique),
        sample_countries  = ("sample_country",      join_unique),
        sample_centers    = ("sample_center",       join_unique),
        sample_regions    = ("sample_region",       join_unique),
        exp_text_blob     = ("exp_text",            lambda s: " ".join(s.dropna())),
        sample_text_blob  = ("sample_text",         lambda s: " ".join(s.dropna())),
    )

    # ── Assemble master join ──────────────────────────────────────────────────
    df = df_studies.merge(exp_agg, on="study_acc", how="left")

    for col in ("n_experiments",):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)

    log.info("Master join: %d studies × %d columns", len(df), len(df.columns))

    # ── Norwegian filter ──────────────────────────────────────────────────────
    signal_cols = [
        "study_text",
        "exp_text_blob",
        "exp_countries", "exp_centers",
        "sample_countries", "sample_centers", "sample_regions",
        "sample_text_blob",
    ]
    signal_cols = [c for c in signal_cols if c in df.columns]

    df["_text_blob"] = (
        df[signal_cols]
        .fillna("")
        .astype(str)
        .agg(" ".join, axis=1)
    )
    mask = df["_text_blob"].str.contains(
        combined_re.pattern, regex=True, flags=re.IGNORECASE, na=False,
    )
    df_nor = df[mask].drop(columns=["_text_blob"])

    log.info("Norwegian filter: kept %d / %d studies", len(df_nor), len(df))

    # ── Serialise ─────────────────────────────────────────────────────────────
    output_cols = {
        "study_acc":          "accession",
        "title":              "title",
        "description":        "description",
        "center_name":        "center_name",
        "first_public_date":  "first_public_date",
        "sample_countries":   "sample_countries_str",
        "sample_centers":     "sample_centers_str",
        "n_experiments":      "n_experiments",
    }
    present = {k: v for k, v in output_cols.items() if k in df_nor.columns}
    df_out = df_nor[list(present.keys())].rename(columns=present)

    def pipe_to_list(s) -> list:
        if pd.isna(s) or s == "":
            return []
        return [x.strip() for x in str(s).split("|") if x.strip()]

    entries: list[dict] = []
    for row in df_out.to_dict(orient="records"):
        row["source"] = "ENA"
        row["domain"] = "sra-study"
        row["sample_countries"] = pipe_to_list(row.pop("sample_countries_str", ""))
        row["sample_centers"]   = pipe_to_list(row.pop("sample_centers_str",   ""))
        if "n_experiments" in row:
            row["n_experiments"] = int(row["n_experiments"]) if pd.notna(row["n_experiments"]) else 0
        entries.append(row)

    out = {
        "join_date":   TODAY,
        "study_count": len(entries),
        "entries":     entries,
    }

    payload = json.dumps(out, indent=2, ensure_ascii=False, default=str)
    (PROC_DIR / f"ena_joined_{TODAY}.json").write_text(payload)
    latest = PROC_DIR / "ena_joined.json"
    latest.write_text(payload)
    log.info("Wrote %s (%d Norwegian studies) ✓", latest, len(entries))


if __name__ == "__main__":
    main()
