#!/usr/bin/env Rscript
# =============================================================================
# plot_norwegian_data.R
# =============================================================================
# Reads all EBI Search raw JSON files (and the joined ENA file), filters for
# Norwegian entries, normalises institution names, and produces static ggplot2
# bar charts + norwegian_entries.csv under output/.
#
# The interactive dashboard is a SEPARATE app: shiny/app.R (run with
# `Rscript -e 'shiny::runApp("shiny")'`).  It reads the CSV this script writes.
#
# NOTE: the plot style (make_inst_palette, theme_nor, the grouped/dodged geom)
# is intentionally duplicated in shiny/app.R — that file must stay self-contained
# for the WebR/shinylive export.  Keep the two copies in sync.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(tibble)
  library(stringr)
  library(forcats)
  library(lubridate)
  library(jsonlite)
  library(stringdist)
  library(scales)
  library(ggtext)
  library(patchwork)
  library(glue)
})

# Null-coalescing operator: fall back to `b` only when `a` is NULL or empty.
# (Pure container-safe coalesce — does NOT inspect a[[1]], so passing a list of
# fields whose first element happens to be empty returns the list unchanged.)
`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0) return(b)
  a
}

#' Return the first present, non-empty, non-NA field value from `fields`,
#' trying `keys` in priority order; `default` if none match.
#' Used where several field names may carry the same information (e.g. a title
#' that may live in name / title / abstract depending on the domain).
pick_field <- function(fields, keys, default = NA_character_) {
  for (k in keys) {
    v <- fields[[k]]
    if (!is.null(v) && length(v) > 0) {
      v1 <- as.character(v[[1]])
      if (!is.na(v1) && nzchar(v1)) return(v1)
    }
  }
  default
}

# ── Paths ─────────────────────────────────────────────────────────────────────
if (requireNamespace("here", quietly = TRUE)) {
  ROOT <- here::here()
} else {
  ROOT <- "."
}

RAW_DIR   <- file.path(ROOT, "data", "raw")
PROC_DIR  <- file.path(ROOT, "data", "processed")
INST_MAP  <- file.path(ROOT, "data", "institution_map.json")
OUT_DIR   <- file.path(ROOT, "output")

# create OUT_DIR if it doesn't exist (recursively), avoid warning if it already exists
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
}

# ── Domain labels (pretty names for plots) ────────────────────────────────────
DOMAIN_LABELS <- c(
  "bioimages"        = "BioImages",
  "biostudies-other" = "BioStudies",
  "metabolights"     = "MetaboLights",
  "pride"            = "PRIDE",
  "biomodels"        = "BioModels",
  "ENA"              = "ENA Studies",
  "sra-sample"       = "ENA Samples"
)

# Non-SRA domains to read from raw/
STANDARD_DOMAINS <- names(DOMAIN_LABELS)[!names(DOMAIN_LABELS) %in% c("ENA", "sra-sample")]

# Domains with fewer total entries than this threshold are excluded from plots.
# Applied before faceting so small domains don't produce near-empty panels.
MIN_DOMAIN_ENTRIES <- 10L

# =============================================================================
# 1.  Institution normalisation
# =============================================================================

inst_map <- fromJSON(INST_MAP, simplifyDataFrame = FALSE)

# Build a lookup: pattern -> canonical name
build_pattern_df <- function(inst_list) {
  rows <- lapply(inst_list, function(i) {
    tibble(
      canonical = i$canonical,
      pattern   = i$patterns
    )
  })
  bind_rows(rows)
}

PATTERN_DF <- build_pattern_df(inst_map$institutions)
NORWAY_RE  <- paste(inst_map$norway_indicators, collapse = "|")

# Email-domain → canonical lookup built once from web_domain fields.
# Keys are base domains (e.g. "uio.no"); matching also handles sub-domains.
DOMAIN_LU <- local({
  cans <- sapply(inst_map$institutions, `[[`, "canonical")
  doms <- sapply(inst_map$institutions, function(i) {
    d <- i[["web_domain"]]
    if (is.null(d) || is.na(d) || !nzchar(d)) NA_character_ else tolower(trimws(d))
  })
  mask <- !is.na(doms)
  setNames(cans[mask], doms[mask])
})

# Canonical name → abbreviation lookup (e.g. "University of Oslo" → "UiO").
ABBREV_LU <- local({
  cans  <- sapply(inst_map$institutions, `[[`, "canonical")
  abbrs <- sapply(inst_map$institutions, `[[`, "abbrev")
  setNames(abbrs, cans)
})

# Map a canonical institution name to its abbreviation; leave unrecognised
# values (including "Other Norway") unchanged.
to_abbrev <- function(canonical) {
  ab <- ABBREV_LU[[canonical]]
  if (!is.null(ab) && !is.na(ab) && nzchar(ab)) ab else canonical
}

#' Robust EBI date parser.
#' Handles: "20230115" (YYYYMMDD compact, most EBI fields),
#'          "2023-01-15", "2023-01-15T00:00:00Z", "2023-01",
#'          "2019 Jan" (biostudies pub_date), "2023", Unix ms integers, NA/empty.
parse_ebi_date <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) return(NA_Date_)
  x <- trimws(as.character(x))
  # Unix milliseconds (13-digit number)
  if (grepl("^\\d{13}$", x)) return(as.Date(as.POSIXct(as.numeric(x) / 1000,
                                                       origin = "1970-01-01")))
  # Try lubridate with progressively looser formats.
  # "Ymd" covers both compact 20230115 and dash-separated 2023-01-15.
  # "Y b" covers "2019 Jan" style returned by biostudies pub_date.
  fmts <- c("Ymd HMS", "Ymd HM", "Ymd", "Y b", "Y-b", "Y-m", "Y/m/d",
            "d/m/Y", "d-m-Y", "d b Y", "b d Y", "Y")
  d <- lubridate::parse_date_time(x, orders = fmts, quiet = TRUE)
  if (is.na(d)) return(NA_Date_)
  d <- as.Date(d)
  if (d > Sys.Date() + lubridate::years(2)) return(NA_Date_)
  d
}

#' Normalise affiliation + email signals to a canonical institution name.
#'
#' Matching priority:
#'   1. Email domain lookup  (@uio.no → "University of Oslo", decisive signal)
#'   2. Regex pattern matching on combined affiliation text
#'   3. Per-token Jaro-Winkler against abbreviations (catches "UiO"/"NTNU" in noisy strings)
#'   4. Full-string Jaro-Winkler against canonical names (informal long-form names)
#'
#' Returns the canonical institution name, or "Other Norway" if nothing matches.
normalise_institution <- function(affil_vec, email_vec = character(0)) {

  # 1. Email domain lookup — highest confidence, very few false positives.
  valid_emails <- email_vec[!is.na(email_vec) & nzchar(email_vec)]
  for (em in valid_emails) {
    dom <- tolower(sub(".*@", "", trimws(em)))
    for (d in names(DOMAIN_LU)) {
      if (dom == d || endsWith(dom, paste0(".", d))) return(DOMAIN_LU[[d]])
    }
  }

  # 2. Regex pattern matching on affiliation text.
  affil <- paste(affil_vec[!is.na(affil_vec)], collapse = " ")
  if (!nzchar(trimws(affil))) return("Other Norway")
  for (i in seq_len(nrow(PATTERN_DF))) {
    if (grepl(PATTERN_DF$pattern[i], affil, ignore.case = TRUE, perl = TRUE)) {
      return(PATTERN_DF$canonical[i])
    }
  }

  # 3 & 4. Fuzzy fallback.
  all_canonical <- sapply(inst_map$institutions, `[[`, "canonical")
  all_abbrevs   <- sapply(inst_map$institutions, `[[`, "abbrev")

  # 3. Per-token JW against abbreviations (abbreviation-length tokens only).
  tokens <- unlist(strsplit(affil, "[,;/()\\s]+"))
  tokens <- tokens[nchar(tokens) >= 2L & nchar(tokens) <= 8L]
  for (tok in tokens) {
    dist_abbr <- stringdist(tolower(tok), tolower(all_abbrevs), method = "jw")
    best <- which.min(dist_abbr)
    if (dist_abbr[best] < 0.15) return(all_canonical[best])
  }

  # 4. Full-string JW against canonical names (run once, outside the token loop).
  dist_full <- stringdist(tolower(affil), tolower(all_canonical), method = "jw")
  best_full <- which.min(dist_full)
  if (dist_full[best_full] < 0.22) return(all_canonical[best_full])

  "Other Norway"
}

#' Return the single affiliation string most likely to have triggered the
#' institution match, for display in the affiliation column.
#'
#' Priority:
#'   1. Email address whose domain matched DOMAIN_LU
#'   2. First affil string that normalise_institution() resolves to a known institution
#'   3. First non-empty affil string (when only the fuzzy fallback fires or nothing matches)
#'
#' This avoids the old pipe-joined blob and lets the user see which piece of
#' metadata actually drove the institution assignment.
pick_affiliation <- function(affil_vec, email_vec = character(0)) {
  # 1. Email domain lookup — return the canonical institution name if it matched.
  for (em in email_vec[!is.na(email_vec) & nzchar(email_vec)]) {
    dom <- tolower(sub(".*@", "", trimws(em)))
    for (d in names(DOMAIN_LU)) {
      if (dom == d || endsWith(dom, paste0(".", d))) return(DOMAIN_LU[[d]])
    }
  }

  # 2. Try each affil string individually; return the first that matches.
  valid <- affil_vec[!is.na(affil_vec) & nzchar(affil_vec)]
  for (s in valid) {
    if (normalise_institution(s) != "Other Norway") return(s)
  }

  # 3. Nothing matched — return the first available string as raw fallback.
  if (length(valid) > 0) valid[[1L]] else NA_character_
}

# =============================================================================
# 2.  Load and flatten all domains
# =============================================================================
 
is_norwegian <- function(values) {
  any(grepl(NORWAY_RE, unlist(values), ignore.case = TRUE))
}
 
#' Parse a single entry from a standard (non-ENA) domain.
#'
#' Every entry reaching this function has already been confirmed Norwegian
#' by fetch_ebi_data.py, so no is_norwegian() guard is needed here.
#' Removing it prevents silent drops when the only Norwegian signal is in
#' a field like labhead_affiliation or submitter_email that the old check
#' did not include in all_text.
parse_entry <- function(entry, domain) {
  fields <- entry$fields %||% list()

  affil_field_names <- c(
    "affiliation", "submitter_affiliation",
    "labhead_affiliation",       # PRIDE lab head
    "labhead",                   # PRIDE lab head name (may carry affil)
    "organisation",              # BioStudies
    "center_name",
    "author",                    # BioImages, BioStudies
    "submitter",
    "first_author",              # BioModels, EGA
    "publication_authors"        # BioModels, EGA
  )

  country_field_names <- c("country", "submitter_country")

  email_field_names <- c(
    "submitter_mail",            # PRIDE, BioModels, EGA
    "submitter_email",           # MetaboLights
    "labhead_mail",              # PRIDE lab head
    "email"                      # EGA
  )

  affil_vals   <- unlist(fields[names(fields) %in% affil_field_names])
  country_vals <- unlist(fields[names(fields) %in% country_field_names])
  email_vals   <- unlist(fields[names(fields) %in% email_field_names])
  email_vals   <- email_vals[grepl("@", email_vals, fixed = TRUE)]

  # ── Date: walk known EBI date field names in priority order ──────────────
  # EBI returns dates as YYYYMMDD (compact) in most fields; parse_ebi_date()
  # handles that format along with ISO and other variants.
  date_raw <- NA_character_
  for (.df in c("submission_date", "creation_date", "pub_date",
                "publication_date", "last_modification_date",
                "collection_date", "release_date", "modified_date",
                "updated_date", "first_public_date")) {
    .v <- fields[[.df]]
    if (!is.null(.v) && length(.v) > 0 && nzchar(as.character(.v[[1]]))) {
      date_raw <- as.character(.v[[1]]); break
    }
  }
  parsed_date <- parse_ebi_date(date_raw)

  tibble(
    domain      = domain,
    accession   = entry$id %||% NA_character_,
    title       = pick_field(fields, c("name", "title", "abstract")),
    affiliation = pick_affiliation(affil_vals, email_vec = email_vals),
    country     = paste(country_vals, collapse = " | "),
    email       = if (length(email_vals) > 0) email_vals[[1]] else NA_character_,
    date        = parsed_date,
    year        = year(parsed_date),
    quarter     = quarter(parsed_date),
    month       = month(parsed_date),
    institution = to_abbrev(as.character(normalise_institution(affil_vals, email_vec = email_vals))[1L]),
    broker      = NA_character_
  )
}

load_standard_domain <- function(domain) {
  path <- file.path(RAW_DIR, domain, "latest.json")
  if (!file.exists(path)) {
    message("  Skipping ", domain, " (no latest.json)")
    return(NULL)
  }
  message("Loading ", domain, " …")
  raw  <- fromJSON(path, simplifyDataFrame = FALSE)
  entries <- raw$entries %||% list()
  rows <- lapply(entries, parse_entry, domain = domain)
  bind_rows(Filter(Negate(is.null), rows))
}
 
#' Parse a joined ENA row from ena_joined.json.
#'
#' join_ena.py now filters for Norwegian entries post-join, so every row
#' here is already confirmed Norwegian.  The is_norwegian() guard is removed
#' to avoid silent drops.  email is NA for ENA rows (not available post-join).
parse_ena_row <- function(row) {
  # join_ena.py outputs: accession, title, center_name, first_public_date,
  # sample_countries (list), sample_centers (list), n_experiments.
  affil_vals   <- c(row$center_name, unlist(row$sample_centers))
  country_vals <- unlist(row$sample_countries)

  # first_public_date from sra-experiment (earliest across experiments for the study).
  # Format is YYYYMMDD compact — handled by parse_ebi_date().
  parsed_date <- parse_ebi_date(as.character(row$first_public_date %||% NA_character_))

  tibble(
    domain      = "ENA",
    accession   = row$accession    %||% NA_character_,
    title       = row$title        %||% NA_character_,
    affiliation = pick_affiliation(affil_vals),
    country     = paste(country_vals, collapse = " | "),
    email       = NA_character_,
    date        = parsed_date,
    year        = year(parsed_date),
    quarter     = quarter(parsed_date),
    month       = month(parsed_date),
    institution = to_abbrev(as.character(normalise_institution(affil_vals))[1L]),
    broker      = row$center_name %||% NA_character_
  )
}

load_ena <- function() {
  path <- file.path(PROC_DIR, "ena_joined.json")
  if (!file.exists(path)) {
    message("  Skipping ENA (no ena_joined.json – has join_ena.py been run?)")
    return(NULL)
  }
  message("Loading ENA joined …")
  raw  <- fromJSON(path, simplifyDataFrame = FALSE)
  rows <- lapply(raw$entries %||% list(), parse_ena_row)
  df   <- bind_rows(Filter(Negate(is.null), rows))

  if (nrow(df) == 0L) return(df)

  # ── Collapse to one row per study accession ──────────────────────────────
  # join_ena.py already produces one row per study, but partitioned fetching
  # can introduce duplicates (the same study appearing in multiple year
  # windows).  This step guarantees one row per accession regardless.
  #
  # Per-column strategy:
  #   accession   – identity (the grouping key)
  #   date        – earliest non-NA date (first public)
  #   title       – first non-NA, non-empty value
  #   country     – union of all pipe-separated values across rows
  #   affiliation – union of all pipe-separated values across rows
  #   email       – first non-NA (ENA rows are always NA here)
  #   institution – majority vote across rows; ties broken by first occurrence
  #   domain      – always "ENA", kept as-is
  #   year/quarter/month – recomputed from the kept date

  collapse_pipe <- function(x) {
    vals <- unlist(strsplit(x[!is.na(x) & nzchar(x)], " | ", fixed = TRUE))
    paste(unique(trimws(vals[nzchar(trimws(vals))])), collapse = " | ")
  }

  majority <- function(x) {
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0L) return(NA_character_)
    names(sort(table(x), decreasing = TRUE))[1L]
  }

  df <- df %>%
    group_by(accession) %>%
    summarise(
      domain      = first(domain),
      title       = first(title[!is.na(title) & nzchar(title)]) %||% NA_character_,
      affiliation = collapse_pipe(affiliation),
      country     = collapse_pipe(country),
      email       = first(email[!is.na(email)]) %||% NA_character_,
      date        = suppressWarnings(min(date, na.rm = TRUE)),
      institution = majority(institution),
      broker      = first(broker[!is.na(broker) & nzchar(broker)]) %||% NA_character_,
      .groups     = "drop"
    ) %>%
    mutate(
      date        = if_else(is.infinite(date), as.Date(NA), date),
      institution = if_else(is.na(institution), "Other Norway", institution),
      year        = year(date),
      quarter     = quarter(date),
      month       = month(date),
    )

  message("  ENA: ", nrow(df), " unique studies after deduplication")
  df
}

#' Parse a single pre-filtered Norwegian sra-sample entry.
#' Produces a row comparable to parse_ena_row but at sample granularity.
parse_sra_sample <- function(entry) {
  f <- entry$fields %||% list()

  # broker_name: the ENA submitter / data broker (fetched field in sra-sample config).
  # center_name: the submitting center / research institution.
  # Both are tried for institution guessing; broker_name also drives the broker column.
  broker_name <- pick_field(f, c("broker_name"))
  center_name <- pick_field(f, c("center_name"))
  country_val <- pick_field(f, c("country"))
  affil_vals  <- Filter(function(x) !is.na(x) && nzchar(x),
                        c(center_name, broker_name))

  date_raw <- NA_character_
  for (.df in c("first_public_date", "collection_date", "last_updated_date")) {
    .v <- f[[.df]]
    if (!is.null(.v) && length(.v) > 0 && nzchar(as.character(.v[[1]]))) {
      date_raw <- as.character(.v[[1]]); break
    }
  }
  parsed_date <- parse_ebi_date(date_raw)

  tibble(
    domain      = "sra-sample",
    accession   = entry$id %||% NA_character_,
    title       = pick_field(f, c("alias", "description")),
    affiliation = pick_affiliation(affil_vals),
    country     = country_val %||% NA_character_,
    email       = NA_character_,
    date        = parsed_date,
    year        = year(parsed_date),
    quarter     = quarter(parsed_date),
    month       = month(parsed_date),
    institution = to_abbrev(as.character(normalise_institution(affil_vals))[1L]),
    broker      = broker_name %||% center_name %||% NA_character_
  )
}

load_sra_samples <- function() {
  path <- file.path(RAW_DIR, "sra-sample", "latest.json")
  if (!file.exists(path)) {
    message("  Skipping ENA Samples (no sra-sample/latest.json – run fetch_ebi_data.py first)")
    return(NULL)
  }
  message("Loading ENA Samples (sra-sample) …")
  raw     <- fromJSON(path, simplifyDataFrame = FALSE)
  entries <- raw$entries %||% list()
  rows    <- lapply(entries, parse_sra_sample)
  df      <- bind_rows(Filter(Negate(is.null), rows))
  message("  ENA Samples: ", nrow(df), " entries")
  df
}

# ── Combine ───────────────────────────────────────────────────────────────────

load_all_data <- function() {
  standard_rows <- lapply(STANDARD_DOMAINS, load_standard_domain)
  ena_rows      <- load_ena()
  sample_rows   <- load_sra_samples()

  df <- bind_rows(c(standard_rows, list(ena_rows), list(sample_rows)))
  df <- df %>%
    mutate(
      domain_label = DOMAIN_LABELS[domain],
      domain_label = if_else(is.na(domain_label), domain, domain_label),
      # Entries without a broker (non-ENA domains) get a label so the
      # broker colour mode can include them with a neutral category.
      broker       = if_else(is.na(broker) | !nzchar(broker), "Non-ENA", broker)
    ) %>%
    filter(!is.na(date))
  df
}
 
 
# =============================================================================
# 3.  Plotting functions
# =============================================================================
 
# Distinct ColorBrewer palette, extended via interpolation when n > 12.
# "Other Norway", "Other", and "Non-ENA" are pinned to grey so they recede.
make_inst_palette <- function(values) {
  values     <- as.character(values)
  grey_keys  <- c("Other Norway", "Other", "Non-ENA")
  non_grey   <- sort(setdiff(values, grey_keys))
  n          <- length(non_grey)

  base_pal <- if (n <= 8) {
    RColorBrewer::brewer.pal(max(3L, n), "Set2")
  } else if (n <= 12) {
    RColorBrewer::brewer.pal(12L, "Set3")
  } else {
    colorRampPalette(
      c(RColorBrewer::brewer.pal(9,  "Set1"),
        RColorBrewer::brewer.pal(8,  "Set2"),
        RColorBrewer::brewer.pal(8,  "Dark2"))
    )(n)
  }

  pal <- setNames(base_pal[seq_along(non_grey)], non_grey)
  for (g in intersect(grey_keys, values)) pal[g] <- "#AAAAAA"
  pal
}
 
theme_nor <- function() {
  theme_classic(base_size = 13) +
    theme(
      plot.title      = element_markdown(face = "bold", size = 15),
      plot.subtitle   = element_markdown(colour = "grey40"),
      axis.text.x     = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title    = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      strip.text      = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey92", colour = NA)
    )
}
 
#' Grouped bar chart, faceted by domain, X axis = time at chosen granularity.
#'
#' @param df               data frame from load_all_data()
#' @param granularity      "year" | "quarter" | "month"
#' @param top_n_inst       keep this many fill values individually; rest → residual category
#' @param domains          character vector of domain_label values to include (NULL = all)
#' @param min_year         drop entries before this year
#' @param min_domain_entries exclude domains with fewer total entries than this threshold
#' @param color_by         "institution" (default) or "broker" (ENA center_name)
plot_time_by_domain <- function(df,
                                granularity        = c("year", "quarter", "month"),
                                top_n_inst         = 12L,
                                domains            = NULL,
                                min_year           = 2000L,
                                min_domain_entries = MIN_DOMAIN_ENTRIES,
                                color_by           = c("institution", "broker")) {
  granularity <- match.arg(granularity)
  color_by    <- match.arg(color_by)

  fill_col    <- color_by                          # column name in df
  other_label <- if (color_by == "institution") "Other Norway" else "Other"
  legend_name <- if (color_by == "institution") "Institution" else "ENA Broker / Center"

  d <- df
  if (!is.null(domains)) d <- d %>% filter(domain_label %in% domains)
  d <- d %>% filter(!is.na(date), year >= min_year, year <= year(Sys.Date()))

  # ── Drop domains below the entry threshold ────────────────────────────────
  domain_counts <- d %>% count(domain_label, name = "total")
  keep_domains  <- domain_counts %>%
    filter(total >= min_domain_entries) %>%
    pull(domain_label)
  dropped <- setdiff(domain_counts$domain_label, keep_domains)
  if (length(dropped) > 0)
    message("  Excluded (< ", min_domain_entries, " entries): ",
            paste(dropped, collapse = ", "))
  d <- d %>% filter(domain_label %in% keep_domains)

  if (nrow(d) == 0) {
    return(ggplot() +
             labs(title = "No domains meet the minimum entry threshold") +
             theme_void())
  }

  # ── Time axis ────────────────────────────────────────────────────────────
  d <- switch(granularity,
    year = d %>% mutate(
      time_val   = as.integer(year),
      time_label = as.character(year)
    ),
    quarter = d %>% mutate(
      time_val   = year + (quarter - 1) / 4,
      time_label = paste0(year, "\nQ", quarter)
    ),
    month = d %>% mutate(
      time_val   = year + (month - 1) / 12,
      time_label = format(date, "%Y\n%b")
    )
  )

  # ── Top-N lumping on the chosen fill column ──────────────────────────────
  top_vals <- d %>%
    count(.data[[fill_col]], sort = TRUE) %>%
    slice_head(n = top_n_inst) %>%
    pull(.data[[fill_col]])

  d <- d %>%
    mutate(fill_val = if_else(.data[[fill_col]] %in% top_vals,
                              .data[[fill_col]], other_label))

  # ── Aggregate ────────────────────────────────────────────────────────────
  counts <- d %>%
    count(domain_label, time_val, time_label, fill_val, name = "n") %>%
    mutate(fill_val = fct_reorder(fill_val, n, .fun = sum) %>%
             fct_relevel(other_label, after = 0L))

  pal      <- make_inst_palette(levels(counts$fill_val))
  x_labels <- counts %>% distinct(time_val, time_label) %>% arrange(time_val)

  # Grouped (dodged) bars — kept in sync with shiny/app.R (see header note).
  ggplot(counts, aes(x = time_val, y = n, fill = fill_val)) +
    geom_col(
      position = position_dodge2(padding = 0.1, preserve = "single"),
      width    = if (granularity == "year") 0.8 else
                 if (granularity == "quarter") 0.22 else 0.07
    ) +
    facet_wrap(~domain_label, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = pal, name = legend_name) +
    scale_x_continuous(
      breaks = x_labels$time_val,
      labels = x_labels$time_label
    ) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .05))) +
    guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
    labs(
      title    = "**Norwegian submissions to EBI repositories**",
      subtitle = glue(
        "Coloured by {color_by} · faceted by repository · granularity: {granularity}"
      ),
      x = NULL, y = "Number of entries"
    ) +
    theme_nor()
}
 
 
# =============================================================================
# 4.  Static output (used by GitHub Actions)
# =============================================================================
 
save_plots <- function(df) {
  message("Saving plots → ", OUT_DIR)
 
  for (gran in c("year", "quarter", "month")) {
    fname <- glue("norwegian_ebi_{gran}.png")
    ggsave(
      file.path(OUT_DIR, fname),
      plot_time_by_domain(df, granularity = gran),
      width = 18, height = 14, dpi = 150
    )
    message("  Saved ", fname)
  }
 
  df %>%
    select(domain, domain_label, accession, title, date, year,
           quarter, month, institution, broker, affiliation, country, email) %>%
    readr::write_csv(file.path(OUT_DIR, "norwegian_entries.csv"))
 
  message("Done ✓  Files in ", OUT_DIR)
}
 
 
# =============================================================================
# 5.  Interactive Shiny app (local debugging)
# =============================================================================
# Loads data live from JSON files via load_all_data() — exercises the full
# parse/normalise pipeline, useful for debugging date parsing, institution
# matching, filter gaps, etc.  Unlike shiny/app.R (which reads the pre-built
# CSV), changes to parse_entry() / normalise_institution() are reflected
# immediately on reload.
#
# Launch:
#   SHINY=1 Rscript R/plot_norwegian_data.R
#   Rscript -e 'source("R/plot_norwegian_data.R"); shiny_app(load_all_data())'

shiny_app <- function(df) {
  for (pkg in c("shiny", "shinythemes", "DT")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required to run the interactive app: ",
           "install.packages('", pkg, "')")
  }
  library(shiny)
  library(shinythemes)
  library(DT)

  domain_choices  <- sort(unique(df$domain_label))
  year_range_data <- range(df$year, na.rm = TRUE)
  latest_date     <- max(df$date, na.rm = TRUE)

  ui <- fluidPage(
    theme = shinytheme("flatly"),
    titlePanel("Norwegian EBI Submissions (live data)"),

    sidebarLayout(
      sidebarPanel(
        width = 3,

        selectInput(
          "granularity", "Time granularity",
          choices  = c("Year" = "year", "Quarter" = "quarter", "Month" = "month"),
          selected = "year"
        ),

        sliderInput(
          "year_range", "Year range",
          min   = year_range_data[1],
          max   = year_range_data[2],
          value = c(max(year_range_data[1], year_range_data[2] - 10L),
                    year_range_data[2]),
          step  = 1L, sep = ""
        ),

        radioButtons(
          "color_by", "Colour bars by",
          choices  = c("Institution" = "institution", "ENA Broker / Center" = "broker"),
          selected = "institution", inline = TRUE
        ),

        sliderInput(
          "top_n_inst", "Top N values to show",
          min = 3, max = 30, value = 8L, step = 1
        ),

        sliderInput(
          "min_domain_entries", "Min entries per repository",
          min = 1, max = 100, value = MIN_DOMAIN_ENTRIES, step = 1
        ),

        checkboxGroupInput(
          "domains", "Repositories",
          choices  = domain_choices,
          selected = domain_choices
        ),

        hr(),
        # Dynamic fill-value checkboxes: rebuilt when color_by / top_n / year /
        # domain filters change.  All top-N values are pre-ticked by default.
        uiOutput("fill_checkbox"),

        hr(),
        p(em(paste("Latest entry:", format(latest_date, "%Y-%m-%d")))),
        p(em(paste(nrow(df), "Norwegian entries loaded from JSON")))
      ),

      mainPanel(
        width = 9,
        plotOutput("main_plot", height = "700px"),
        hr(),
        DT::dataTableOutput("entry_table")
      )
    )
  )

  server <- function(input, output, session) {

    # Filtered base (year + domain); institution lumping applied on top.
    base_df <- reactive({
      req(input$year_range, input$domains)
      df |>
        filter(domain_label %in% input$domains,
               !is.na(year),
               year >= input$year_range[1],
               year <= input$year_range[2])
    })

    # Top-N fill values for the current color_by mode, year, and domain window.
    # The residual sentinel ("Other Norway" / "Other") is always appended.
    top_fill_vals <- reactive({
      req(input$top_n_inst, input$color_by)
      col        <- input$color_by
      other_lbl  <- if (col == "institution") "Other Norway" else "Other"
      top_vals   <- base_df() |>
        count(.data[[col]], sort = TRUE) |>
        slice_head(n = input$top_n_inst) |>
        pull(.data[[col]])
      unique(c(top_vals, other_lbl))
    })

    output$fill_checkbox <- renderUI({
      vals  <- top_fill_vals()
      label <- if (input$color_by == "institution") "Show institutions"
               else "Show brokers / centers"
      checkboxGroupInput("selected_fill", label,
                         choices = vals, selected = vals)
    })

    output$main_plot <- renderPlot({
      req(input$year_range, input$top_n_inst, input$color_by)

      col       <- input$color_by
      other_lbl <- if (col == "institution") "Other Norway" else "Other"

      # Lump the fill column to top-N, then apply checkbox filter.
      # plot_time_by_domain is called with top_n_inst = 999 to skip re-lumping.
      top_vals <- setdiff(top_fill_vals(), other_lbl)
      d <- base_df() |>
        mutate(across(all_of(col),
                      ~ if_else(.x %in% top_vals, .x, other_lbl)))

      if (!is.null(input$selected_fill) && length(input$selected_fill) > 0)
        d <- d |> filter(.data[[col]] %in% input$selected_fill)

      plot_time_by_domain(
        d,
        granularity        = input$granularity,
        top_n_inst         = 999L,
        color_by           = col,
        min_domain_entries = input$min_domain_entries
      )
    }, res = 120)

    output$entry_table <- DT::renderDataTable({
      req(input$year_range, input$domains)
      keep <- base_df() |>
        count(domain_label, name = "total") |>
        filter(total >= input$min_domain_entries) |>
        pull(domain_label)
      base_df() |>
        filter(domain_label %in% keep) |>
        select(
          Repository  = domain_label,
          Accession   = accession,
          Title       = title,
          Date        = date,
          Institution = institution,
          Broker      = broker,
          Email       = email
        ) |>
        arrange(desc(Date))
    }, options = list(pageLength = 10, scrollX = TRUE), filter = "top")
  }

  shinyApp(ui, server)
}


# =============================================================================
# 6.  Entry point
# =============================================================================
# Default:   Rscript R/plot_norwegian_data.R          → saves static PNGs + CSV
# Shiny:     SHINY=1 Rscript R/plot_norwegian_data.R  → launches interactive app

main <- function() {
  df <- load_all_data()
  message(glue("Loaded {nrow(df)} Norwegian entries across {n_distinct(df$domain)} domains"))

  if (nrow(df) == 0) {
    message("No data found – have you run fetch_ebi_data.py yet?")
    return(invisible(NULL))
  }

  if (identical(Sys.getenv("SHINY"), "1")) {
    shiny_app(df)
  } else {
    save_plots(df)
  }
}

main()
