# Norwegian EBI Submissions Dashboard

> **Weekly-updated tracker of Norwegian research data deposited in EBI repositories.**  
> Data is fetched automatically via GitHub Actions, institution names are normalised
> against a curated list, and bar plots are rendered with R/ggplot2 (Shiny-ready).

---

## Repositories tracked

| Repository | EBI domain | Primary affiliation fields |
|---|---|---|
| BioImages | `bioimages` | `author`
| BioStudies | `biostudies-other` | `organisation` |
| MetaboLights | `metabolights` | `submitter_affiliation` |
| PRIDE | `pride` | `submitter_affiliation`, `submitter_country` |
| BioModels | `biomodels` | `submitter_affiliation` |
| ENA/SRA | `sra-study`, `sra-experiment`, `sra-sample` | `center_name`, `country` (sample) |
| EGA | `ega`, `ega-sample` (EGA Metadata API, not EBI Search) | DAC contact `institution_name`, `email` |

> **EGA** is fetched from the separate [EGA Public Metadata API](https://metadata.ega-archive.org)
> rather than EBI Search (which exposes no dates/affiliations for it).  Data Access
> Committees (DACs) are filtered to Norwegian ones — by contact `institution_name`
> or by email **domain** (`.no` or a known institution `web_domain`) — and the
> records below the matching DACs are plotted:
> - `ega` — **studies** (DAC → datasets → studies).
> - `ega-sample` — **samples** (DAC → datasets → samples), mirroring `ENA Samples`.
>   Samples have no public date/affiliation, so each inherits the parent dataset's
>   release date and the parent DAC's institution/email.

---

## Architecture

```
GitHub Actions (cron 02:30 UTC)
│
├─ scripts/fetch_ebi_data.py
│    ├─ Queries EBI Search REST API  (GET /ws/rest/{domain}?query=Norway…)
│    ├─ Paginates until all hits retrieved (PAGE_SIZE = 500)
│    └─ Writes data/raw/{domain}/latest.json  (+ dated snapshot)
│
├─ scripts/join_ena.py
│    ├─ Loads sra-study, sra-experiment, sra-sample
│    ├─ Joins studies → experiments → samples via shared accession keys
│    └─ Writes data/processed/ena_joined.json
│
├─ scripts/fetch_ega.py
│    ├─ Queries EGA Public Metadata API  (GET metadata.ega-archive.org/dacs…)
│    ├─ Keeps Norwegian DACs (institution_name regex / email .no|web_domain)
│    ├─ Walks DAC → datasets → studies + samples
│    └─ Writes data/raw/ega/latest.json + data/raw/ega-sample/latest.json
│       (same {id, fields} shape as the EBI Search domains)
│
├─ scripts/fetch_identifiers.py
│    ├─ Caches the identifiers.org namespace registry (prefix → pattern)
│    └─ Writes data/identifiers_namespaces.json
│
└─ R/plot_norwegian_data.R
     ├─ Loads all latest.json + ena_joined.json
     ├─ Filters rows mentioning Norway/Norge/Norwegian in affil/country fields
     ├─ Normalises institution names (regex → fuzzy fallback → "Other Norway")
     ├─ Builds identifiers.org links (validated against the registry pattern)
     └─ Saves PNG plots + norwegian_entries.csv to output/
```

### Accession links (identifiers.org)

Each domain declares an `identifiers_prefix` in its definition (e.g. `pride` →
`pride.project`, `sra-study` → `["insdc.sra", "bioproject"]`).  The render step
checks every accession against that prefix's [identifiers.org](https://identifiers.org)
registry pattern and, on a match, builds a `https://identifiers.org/<prefix>:<acc>`
link (an unmatched or mis-mapped accession simply gets no link).  These render as
clickable accessions in the Shiny dashboard's table.  EGA samples (`EGAN…`) have
no identifiers.org namespace, so they are left unlinked.

---

## Output plots

| File | Description |
|---|---|
| `output/01_overview_by_year.png` | Stacked bar: entries per year, coloured by repository |
| `output/02_by_institution_all.png` | Faceted bar: entries per institution, one panel per repository |
| `output/03_time_year_by_domain.png` | Grouped bars: year × domain |
| `output/04_time_quarter_by_inst.png` | Quarterly trend, coloured by institution |
| `output/norwegian_entries.csv` | Flat export of all matched entries |

---

## Institution normalisation

Institution names in EBI metadata are free text – a single lab may appear as
`"UiO"`, `"University of Oslo"`, `"Universitetet i Oslo"`, `"Univ. of Oslo"`, etc.

Normalisation happens in three stages (see `data/institution_map.json`):

1. **Regex matching** – each institution has a list of case-insensitive Perl
   regex patterns covering abbreviations, Norwegian/English spellings, and
   common misspellings.
2. **Fuzzy fallback** – if no regex matches, the affiliation string is compared
   against all canonical names and abbreviations using Jaro-Winkler distance
   (`stringdist`).
3. **"Other Norway"** – entries where Norway can be detected (city names,
   country field = NO, etc.) but no institution is matched.

### Adding an institution

Edit `data/institution_map.json` and add an entry:

```json
{
  "canonical": "My Institute",
  "canonical_no": "Mitt Institutt",
  "abbrev": "MI",
  "ror": "https://ror.org/...",
  "patterns": [
    "My Institute",
    "Mitt Institutt",
    "\\bMI\\b"
  ]
}
```

---
## OS Dependencies

On Ubuntu 26.04 LTS

```bash

sudo apt-get install r-base build-essential pkg-config git make \
libcurl4-openssl-dev libssl-dev libxml2-dev zlib1g-dev libicu-dev \
libfontconfig1-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
libharfbuzz-dev libfribidi-dev libcairo2-dev libpango1.0-dev libx11-dev libxt-dev \
libblas-dev liblapack-dev libopenblas-dev gfortran libgfortran5 \
liblzma-dev libbz2-dev libreadline-dev libsqlite3-dev libpq-dev \
libgit2-dev libgmp-dev libglpk-dev imagemagick

```



## Running locally

### 1. Fetch data

```bash
pip install requests tenacity
python scripts/fetch_ebi_data.py
python scripts/join_ena.py
```

### 2. Render plots (static)

```r
Rscript R/plot_norwegian_data.R
# → output/*.png  +  output/norwegian_entries.csv
```

### 3. Launch Shiny app

Two ways to run the interactive dashboard:

**A. Quick / debugging** — loads directly from the raw JSON files, so changes to
parsing or institution normalisation are reflected immediately on reload:

```bash
SHINY=1 Rscript R/plot_norwegian_data.R
```

**B. Production / shinylive** — reads the pre-built CSV, identical to what is
deployed via shinylive. Run the render step first:

```bash
cp output/norwegian_entries.csv shiny/data/norwegian_entries.csv
Rscript -e 'shiny::runApp("shiny")'
```

Both apps expose the same controls:
- Time granularity (year / quarter / month)
- Year range (default: last 10 years)
- Top N institutions (default: 8), with per-institution checkboxes
- Repository filter (checkboxes) and minimum entries per repository

---

## GitHub Actions workflows

| Workflow | Trigger | Description |
|---|---|---|
| `01_fetch_data.yml` | cron 02:30 UTC + manual | Runs fetch + join scripts, commits data, triggers plot workflow |
| `02_render_plots.yml` | called by #01 + manual | Installs R deps, renders plots, commits output/ |

Both workflows require the default `GITHUB_TOKEN` with **write** access to
`contents` (enabled automatically in public repos; check repo Settings →
Actions → General if you have issues).

---

## ENA join logic

ENA data is fetched from three EBI Search sub-tables and joined by
`join_ena.py`:

```
sra-study  ←─ study_accession ─→  sra-experiment
sra-experiment ←─ sample_accession ─→  sra-sample
```

The joined row uses the **study** as its unit; sample-level country
information (best source for Norway detection) is propagated upward.

Three sub-tables are intentionally omitted:

- **sra-submission** (31M entries, no searchable date field — cannot be year-partitioned)
- **sra-analysis** (24M entries, date field mismatch)
- **sra-run** (42M entries, no searchable date field — was only used for run counts)

---

## Data refresh cadence

By default the workflow runs once per day at 02:30 UTC. To change, edit the
`cron` line in `.github/workflows/01_fetch_data.yml`.

The dated raw snapshots (e.g. `data/raw/pride/2026-05-18.json`) are excluded
from git (see `.gitignore`) to keep the repository lean; only `latest.json`
per domain is committed.
