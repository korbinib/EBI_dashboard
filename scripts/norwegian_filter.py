#!/usr/bin/env python3
"""
norwegian_filter.py
===================
Norwegian-entry detection shared by the fetch and join stages.

A single copy of the geo/institution/email filter lives here so the two
consumers (fetch_ebi_data.py and join_ena.py) cannot drift apart.

Detection signals
-----------------
  a) Geographic indicators  (Norway, Oslo, Bergen, Tromsø, \\bNO\\b …)
  b) Institution regexes    (NTNU, Folkehelseinstituttet, …)
  c) Norwegian TLD email    (@*.no)

All three are sourced from data/institution_map.json.
"""

import json
import logging
import re

from paths import INSTITUTION_MAP

log = logging.getLogger("norwegian_filter")

_EMAIL_NO_RE = re.compile(r"@[\w.\-]+\.no\b", re.IGNORECASE)


def load_geo_tokens(path=INSTITUTION_MAP) -> list[str]:
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


def load_institution_regexes(path=INSTITUTION_MAP) -> list[re.Pattern]:
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


def make_combined_filter(geo_re: re.Pattern,
                         inst_regexes: list[re.Pattern]) -> re.Pattern:
    """A single regex OR-ing geo, institution, and Norwegian-email patterns."""
    all_patterns = [geo_re.pattern] + [p.pattern for p in inst_regexes] + \
                   [_EMAIL_NO_RE.pattern]
    return re.compile("|".join(f"(?:{p})" for p in all_patterns), re.IGNORECASE)


def is_norwegian_entry(entry: dict, geo_re: re.Pattern,
                       inst_regexes: list[re.Pattern]) -> bool:
    """True if any field value in an EBI Search entry carries a Norwegian signal."""
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
