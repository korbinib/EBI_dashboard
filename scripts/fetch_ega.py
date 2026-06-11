#!/usr/bin/env python3
"""
fetch_ega.py
============
Fetch Norwegian EGA (European Genome-phenome Archive) studies and write them to
data/raw/ega/latest.json in the same {id, fields} shape every other domain uses,
so R/plot_norwegian_data.R picks them up through its generic parse_entry() path
(only a DOMAIN_LABELS entry is needed on the R side).

Why EGA is a separate fetch (not a DOMAINS entry)
-------------------------------------------------
EGA is invisible to the EBI Search API used by every other domain: that index
returns only id/name/description for EGA, with no dates and no affiliations, so
Norwegian detection is impossible there.  Instead this script talks to the EGA
Public Metadata API (ega_api.py), which exposes the full object graph.

Strategy (as requested: DACs → filter → matching studies)
---------------------------------------------------------
1. DACS     – fetch every Data Access Committee (EGAC).  Each carries a
              `contacts` array with institution_name + email.
2. FILTER   – keep a DAC when a contact looks Norwegian, judged precisely:
                • institution_name matches a geo token or institution regex, OR
                • the email DOMAIN ends in .no or a known institution web_domain
                  (institution_map.json).  Matching the email *domain* — not the
                  whole address — avoids false hits from name-like local parts
                  (e.g. "nina.gasparoni@uni-saarland.de" must not match NINA).
3. RECORDS  – for each Norwegian DAC walk DAC → datasets (EGAD), and from each
              dataset collect both its studies (EGAS) and its samples (EGAN),
              de-duplicating records that several datasets/DACs share.
4. SHAPE    – emit each study/sample as an EBI-Search-style entry whose fields
              carry the matched Norwegian institution_name(s) and email(s), so
              the R institution normaliser (regex + email web_domain lookup)
              resolves the right institution.  Samples have no public date or
              affiliation of their own, so they inherit the parent dataset's
              released_date and the parent DAC's institution/email — mirroring
              how the "ENA Samples" (sra-sample) domain is plotted.
5. SAVE     – data/raw/ega/latest.json and data/raw/ega-sample/latest.json
              (+ dated copies) via save_domain().

Output
------
  data/raw/ega/latest.json          – Norwegian EGA studies  ({id, fields})
  data/raw/ega-sample/latest.json   – Norwegian EGA samples  ({id, fields})
  data/raw/<domain>/<date>.json     – dated snapshots
"""

import logging
import sys
import time
from pathlib import Path

# Ensure scripts/ is importable regardless of how this file is invoked.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from ega_api import (
    RATE_SLEEP, get_dacs, get_dac_datasets, get_dataset_studies,
    get_dataset_samples,
)
from norwegian_filter import (
    load_geo_tokens, load_institution_regexes, build_geo_regex,
    load_web_domains, email_domain_is_norwegian,
)
from fetch_ebi_data import save_domain

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("fetch_ega")

DOMAIN        = "ega"          # EGA studies (EGAS)
DOMAIN_SAMPLE = "ega-sample"  # EGA samples (EGAN), mirrors the sra-sample domain

# Fields populated on each emitted entry (recorded in latest.json metadata).
FIELDS_USED = [
    "name", "title", "description", "study_type",
    "submitter_affiliation", "submitter_mail", "email",
    "release_date", "dac_accession",
]


# ──────────────────────────────────────────────────────────────────────────────
# DAC Norwegian-signal extraction
# ──────────────────────────────────────────────────────────────────────────────

def _institution_is_norwegian(name: str, geo_re, inst_regexes) -> bool:
    """True if a free-text institution name carries a Norwegian signal."""
    if not name or not name.strip():
        return False
    if geo_re.search(name):
        return True
    return any(p.search(name) for p in inst_regexes)


def dac_norwegian_signal(dac: dict, geo_re, inst_regexes,
                         web_domains: set[str]) -> tuple[list[str], list[str]]:
    """
    Inspect a DAC's contacts and return (institutions, emails) carrying a
    Norwegian signal.  An empty-and-empty return means the DAC is not Norwegian.

    institutions – institution_name values matched by geo/institution regex
    emails        – contact emails whose DOMAIN is Norwegian (.no or web_domain)
    """
    institutions: list[str] = []
    emails:       list[str] = []
    for c in (dac.get("contacts") or []):
        inst  = (c.get("institution_name") or "").strip()
        email = (c.get("email") or "").strip()
        if inst and _institution_is_norwegian(inst, geo_re, inst_regexes):
            institutions.append(inst)
        if email and email_domain_is_norwegian(email, web_domains):
            emails.append(email)
    # De-duplicate while preserving order.
    return list(dict.fromkeys(institutions)), list(dict.fromkeys(emails))


# ──────────────────────────────────────────────────────────────────────────────
# Study collection
# ──────────────────────────────────────────────────────────────────────────────

def _iso_to_date(value) -> str:
    """'2012-07-18T17:56:23+02:00' → '2012-07-18'; None/'' → ''."""
    if not value:
        return ""
    return str(value).split("T", 1)[0]


def _accumulate(store: dict[str, dict], acc: str, base: dict,
                insts: list[str], emails: list[str], dac_acc: str) -> None:
    """Insert/merge one record into `store`, OR-ing in this path's signals.

    A study or sample reached via several datasets/DACs keeps every Norwegian
    institution_name, email, and source DAC seen along the way.
    """
    rec = store.get(acc)
    if rec is None:
        rec = {**base, "institutions": [], "emails": [], "dacs": []}
        store[acc] = rec
    for x in insts:
        if x not in rec["institutions"]:
            rec["institutions"].append(x)
    for x in emails:
        if x not in rec["emails"]:
            rec["emails"].append(x)
    if dac_acc not in rec["dacs"]:
        rec["dacs"].append(dac_acc)


def collect_norwegian_records(geo_re, inst_regexes,
                              web_domains: set[str]) -> tuple[dict, dict]:
    """
    Walk Norwegian DACs → datasets and, from each dataset, collect both its
    studies (EGAS) and its samples (EGAN).  Returns (studies, samples), each a
    dict keyed by accession.

    Samples carry no public date or affiliation, so they inherit the parent
    dataset's released_date and the parent DAC's Norwegian institution/email —
    the same shape the studies use, so both flow through one R parser.
    """
    log.info("Fetching all DACs …")
    dacs = get_dacs()
    log.info("  %d DACs returned by the API", len(dacs))

    nor_dacs = []
    for dac in dacs:
        acc = str(dac.get("accession_id") or "")
        if not acc.startswith("EGAC"):
            continue
        insts, emails = dac_norwegian_signal(dac, geo_re, inst_regexes, web_domains)
        if insts or emails:
            nor_dacs.append((acc, insts, emails))
    log.info("  %d Norwegian DACs after filtering", len(nor_dacs))

    studies: dict[str, dict] = {}
    samples: dict[str, dict] = {}
    for i, (dac_acc, insts, emails) in enumerate(nor_dacs, 1):
        try:
            datasets = get_dac_datasets(dac_acc)
        except Exception as exc:
            log.warning("  [%d/%d] %s: dataset fetch failed: %s",
                        i, len(nor_dacs), dac_acc, exc)
            continue

        seen_studies: set[str] = set()
        n_samples = 0
        for ds in datasets:
            if ds.get("is_deprecated"):
                continue
            ds_acc = str(ds.get("accession_id") or "")
            if not ds_acc:
                continue
            ds_title = ds.get("title") or ""
            ds_date  = _iso_to_date(ds.get("released_date"))

            # ── studies (EGAS) ────────────────────────────────────────────────
            try:
                ds_studies = get_dataset_studies(ds_acc)
            except Exception as exc:
                log.warning("    %s/%s: study fetch failed: %s", dac_acc, ds_acc, exc)
                ds_studies = []
            time.sleep(RATE_SLEEP)
            for st in ds_studies:
                if st.get("is_deprecated"):
                    continue
                s_acc = str(st.get("accession_id") or "")
                if not s_acc:
                    continue
                seen_studies.add(s_acc)
                _accumulate(studies, s_acc, {
                    "accession":    s_acc,
                    "title":        st.get("title") or "",
                    "description":  st.get("description") or "",
                    "study_type":   st.get("study_type") or "",
                    "release_date": _iso_to_date(st.get("released_date")),
                }, insts, emails, dac_acc)

            # ── samples (EGAN) ────────────────────────────────────────────────
            try:
                ds_samples = get_dataset_samples(ds_acc)
            except Exception as exc:
                log.warning("    %s/%s: sample fetch failed: %s", dac_acc, ds_acc, exc)
                ds_samples = []
            time.sleep(RATE_SLEEP)
            for sm in ds_samples:
                if sm.get("is_deprecated"):
                    continue
                m_acc = str(sm.get("accession_id") or "")
                if not m_acc:
                    continue
                n_samples += 1
                _accumulate(samples, m_acc, {
                    # Sample title/description are usually null; fall back to the
                    # dataset title so the entry has a meaningful display label.
                    "accession":    m_acc,
                    "title":        sm.get("title") or ds_title,
                    "description":  sm.get("description") or "",
                    "study_type":   "",
                    "release_date": ds_date,
                }, insts, emails, dac_acc)

        log.info("  [%d/%d] %s → %d datasets, %d studies, %d samples",
                 i, len(nor_dacs), dac_acc, len(datasets),
                 len(seen_studies), n_samples)

    log.info("Collected %d unique Norwegian EGA studies, %d samples",
             len(studies), len(samples))
    return studies, samples


# ──────────────────────────────────────────────────────────────────────────────
# Entry shaping (EBI-Search {id, fields} format)
# ──────────────────────────────────────────────────────────────────────────────

def record_to_entry(rec: dict) -> dict:
    """Shape one accumulated study/sample record as an EBI-Search-style entry.

    Field values are lists to match the EBI Search format that R's parse_entry()
    expects.  institution_name(s) go in submitter_affiliation and email(s) in
    submitter_mail/email so the R normaliser resolves the institution via both
    its regex patterns and its email web_domain lookup.
    """
    fields: dict[str, list[str]] = {
        "name":        [rec["title"]] if rec["title"] else [],
        "title":       [rec["title"]] if rec["title"] else [],
        "description": [rec["description"]] if rec["description"] else [],
        "release_date": [rec["release_date"]] if rec["release_date"] else [],
        "dac_accession": list(rec["dacs"]),
    }
    if rec["study_type"]:
        fields["study_type"] = [rec["study_type"]]
    if rec["institutions"]:
        fields["submitter_affiliation"] = list(rec["institutions"])
    if rec["emails"]:
        fields["submitter_mail"] = list(rec["emails"])
        fields["email"] = list(rec["emails"])
    # Drop empty lists so the entry mirrors a real (sparsely populated) EBI entry.
    fields = {k: v for k, v in fields.items() if v}
    return {"id": rec["accession"], "fields": fields}


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def main() -> int:
    geo_re       = build_geo_regex(load_geo_tokens())
    inst_regexes = load_institution_regexes()
    web_domains  = load_web_domains()
    log.info("Filter ready: %d institution patterns, %d web domains",
             len(inst_regexes), len(web_domains))

    studies, samples = collect_norwegian_records(geo_re, inst_regexes, web_domains)

    study_entries  = [record_to_entry(rec) for rec in studies.values()]
    sample_entries = [record_to_entry(rec) for rec in samples.values()]

    save_domain(DOMAIN,        study_entries,  FIELDS_USED)
    save_domain(DOMAIN_SAMPLE, sample_entries, FIELDS_USED)
    log.info("=== Done: ega (%d studies), ega-sample (%d samples) ===",
             len(study_entries), len(sample_entries))
    return 0


if __name__ == "__main__":
    sys.exit(main())
