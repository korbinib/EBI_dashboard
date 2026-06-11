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


def load_web_domains(path=INSTITUTION_MAP) -> set[str]:
    """
    Return the set of lowercase Norwegian institution web-domains declared in the
    institution map (e.g. {"uib.no", "uio.no", …}).  Used by
    email_domain_is_norwegian() to recognise affiliations from a contact's email
    address even when no free-text affiliation string is present.
    """
    try:
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
    except FileNotFoundError:
        log.warning("Institution map not found: %s – no web domains loaded", path)
        return set()
    result: set[str] = set()
    for inst in data.get("institutions", []):
        d = inst.get("web_domain")
        if d and isinstance(d, str) and d.strip():
            result.add(d.strip().lower())
    return result


def email_domain_is_norwegian(email: str, web_domains: set[str]) -> bool:
    """
    True if an email address carries a Norwegian signal, judged **only on its
    domain part** (everything after the last '@').

    Two ways to qualify:
      a) the domain ends in the Norwegian ccTLD '.no'         (e.g. *@uib.no)
      b) the domain equals, or is a sub-domain of, a known institution
         web_domain from the institution map                 (e.g. *@ous-research.no)

    Matching the domain rather than the whole address avoids false positives
    from name-like local parts (e.g. "nina.gasparoni@uni-saarland.de" must NOT
    match the NINA institution pattern).
    """
    if not email or "@" not in email:
        return False
    domain = email.rsplit("@", 1)[1].strip().lower().rstrip(".")
    if not domain:
        return False
    if domain == "no" or domain.endswith(".no"):
        return True
    return any(domain == wd or domain.endswith("." + wd) for wd in web_domains)


def load_geo_tokens(path=INSTITUTION_MAP) -> list[str]:
    try:
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
    except FileNotFoundError:
        log.warning("Institution map not found: %s – using fallback geo tokens", path)
        return ["Norway", "Norge", "Norwegian", "Norsk", "Oslo", "Bergen",
                "Trondheim", "Tromsø", "Stavanger"]
    # Indicators are treated as regex patterns (matching the R side's NORWAY_RE),
    # so entries like "Troms[øo]" are kept rather than dropped.  We only skip
    # blanks and patterns that fail to compile, keeping Python and R detection in
    # sync instead of silently weaker on the fetch side.
    result: list[str] = []
    for t in data.get("norway_indicators", []):
        t = t.strip()
        if not t:
            continue
        try:
            re.compile(t)
        except re.error as exc:
            log.debug("Skipping invalid norway_indicator %r: %s", t, exc)
            continue
        result.append(t)
    return sorted(set(result))


def _institution_name_patterns(inst: dict) -> list[str]:
    """
    Escaped, word-bounded regexes for an institution's identifying names so an
    entry mentioning any of them counts as Norwegian, even when the curated
    `patterns` list doesn't spell that variant out:

      canonical      English name      "University of Bergen"
      canonical_no   Norwegian name    "Universitetet i Bergen"
      abbrev         abbreviation      "UiB"
      ror            ROR id            "03zga2b32" (matched, not the full URL)

    Escaping + \\b boundaries keep these literal (no accidental regex meaning,
    and "UiB" won't match inside another word).
    """
    out: list[str] = []
    for key in ("canonical", "canonical_no", "abbrev"):
        v = inst.get(key)
        if isinstance(v, str) and v.strip():
            out.append(r"\b" + re.escape(v.strip()) + r"\b")
    ror = inst.get("ror")
    if isinstance(ror, str) and ror.strip():
        ror_id = ror.strip().rstrip("/").rsplit("/", 1)[-1]   # id after last '/'
        if ror_id:
            out.append(r"\b" + re.escape(ror_id) + r"\b")
    return out


def load_institution_regexes(path=INSTITUTION_MAP) -> list[re.Pattern]:
    try:
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
    except FileNotFoundError:
        log.warning("Institution map not found: %s – institution filter disabled", path)
        return []
    compiled: list[re.Pattern] = []
    for inst in data.get("institutions", []):
        # Curated regex patterns first, then literal name/abbrev/ROR fallbacks.
        for p in list(inst.get("patterns", [])) + _institution_name_patterns(inst):
            try:
                compiled.append(re.compile(p, re.IGNORECASE))
            except re.error as exc:
                log.debug("Skipping invalid pattern %r: %s", p, exc)
    log.info("Loaded %d institution regex patterns", len(compiled))
    return compiled


def build_geo_regex(geo_tokens: list[str]) -> re.Pattern:
    # geo_tokens are regex patterns (see load_geo_tokens), joined as a case-
    # insensitive alternation — same as R's NORWAY_RE.  The bare country code
    # "NO" is matched *case-sensitively* via a scoped (?-i:) group so the English
    # word "no"/"No" in free text does not produce a false Norwegian hit.
    parts = list(geo_tokens) + [r"(?-i:\bNO\b)"]
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
