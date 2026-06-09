#!/usr/bin/env Rscript
# =============================================================================
# plot_norwegian_data.R
# =============================================================================
# Reads all EBI Search raw JSON files (and the joined ENA file), filters for
# Norwegian entries, normalises institution names, and produces ggplot2 bar
# charts saved under output/.
#
# The script is also structured to run as a Shiny app – set SHINY=TRUE via
# environment variable or call it directly with `shiny::runApp("R/")`.
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

RUN_AS_SHINY <- nzchar(Sys.getenv("SHINY"))  # set env var SHINY=1 to run app

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
  "ENA"              = "ENA/SRA"
)

# Non-SRA domains to read from raw/
STANDARD_DOMAINS <- names(DOMAIN_LABELS)[names(DOMAIN_LABELS) != "ENA"]

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
  if (!is.na(d)) return(as.Date(d))
  NA_Date_
}

#' Normalise a free-text affiliation string to a canonical institution name.
#' Returns the canonical name, or "Other Norway" if no match.
normalise_institution <- function(affil_vec) {
  if (is.null(affil_vec) || length(affil_vec) == 0) return("Other Norway")
  affil <- paste(affil_vec, collapse = " ")
  if (is.na(affil) || affil == "") return("Other Norway")

  for (i in seq_len(nrow(PATTERN_DF))) {
    if (grepl(PATTERN_DF$pattern[i], affil, ignore.case = TRUE, perl = TRUE)) {
      return(PATTERN_DF$canonical[i])
    }
  }

  # Fuzzy fallback: compare against all canonical names (+ Norwegian names)
  all_canonical  <- sapply(inst_map$institutions, `[[`, "canonical")
  all_canonical_no <- sapply(inst_map$institutions, `[[`, "canonical_no")
  all_abbrevs    <- sapply(inst_map$institutions, `[[`, "abbrev")

  # Tokenise affil into words and try each against canonicals
  tokens <- unlist(strsplit(affil, "[,;/()\\s]+"))
  tokens <- tokens[nchar(tokens) > 3]

  for (tok in tokens) {
    distances <- stringdist(tolower(tok), tolower(all_abbrevs), method = "jw")
    best_abbrev <- which.min(distances)
    if (distances[best_abbrev] < 0.12) {
      return(all_canonical[best_abbrev])
    }
    distances2 <- stringdist(tolower(affil), tolower(all_canonical), method = "jw")
    best2 <- which.min(distances2)
    if (distances2[best2] < 0.18) {
      return(all_canonical[best2])
    }
  }

  "Other Norway"
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
    title       = (fields[["name"]] %||% fields[["title"]] %||%
                   fields[["abstract"]] %||% list(NA_character_))[[1]],
    affiliation = paste(affil_vals, collapse = " | "),
    country     = paste(country_vals, collapse = " | "),
    email       = if (length(email_vals) > 0) email_vals[[1]] else NA_character_,
    date        = parsed_date,
    year        = year(parsed_date),
    quarter     = quarter(parsed_date),
    month       = month(parsed_date),
    institution = as.character(normalise_institution(affil_vals))[1L]
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
    affiliation = paste(affil_vals, collapse = " | "),
    country     = paste(country_vals, collapse = " | "),
    email       = NA_character_,
    date        = parsed_date,
    year        = year(parsed_date),
    quarter     = quarter(parsed_date),
    month       = month(parsed_date),
    institution = as.character(normalise_institution(affil_vals))[1L]
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

# ── Combine ───────────────────────────────────────────────────────────────────

load_all_data <- function() {
  standard_rows <- lapply(STANDARD_DOMAINS, load_standard_domain)
  ena_rows      <- load_ena()

  df <- bind_rows(c(standard_rows, list(ena_rows)))
  df <- df %>%
    mutate(
      domain_label = DOMAIN_LABELS[domain],
      domain_label = if_else(is.na(domain_label), domain, domain_label)
    )
  df
}
 
 
# =============================================================================
# 3.  Plotting functions
# =============================================================================
 
# Distinct ColorBrewer palette, extended via interpolation when n > 12.
# "Other Norway" is always pinned to mid-grey so it recedes visually.
make_inst_palette <- function(institutions) {
  institutions <- as.character(institutions)
  non_other    <- sort(setdiff(institutions, "Other Norway"))
  n            <- length(non_other)
 
  base_pal <- if (n <= 8) {
    RColorBrewer::brewer.pal(max(3L, n), "Set2")
  } else if (n <= 12) {
    RColorBrewer::brewer.pal(12L, "Set3")
  } else {
    # Blend Set1 + Set2 + Dark2 for maximum distinctiveness
    colorRampPalette(
      c(RColorBrewer::brewer.pal(9,  "Set1"),
        RColorBrewer::brewer.pal(8,  "Set2"),
        RColorBrewer::brewer.pal(8,  "Dark2"))
    )(n)
  }
 
  pal <- setNames(base_pal[seq_along(non_other)], non_other)
  if ("Other Norway" %in% institutions) pal["Other Norway"] <- "#AAAAAA"
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
 
#' Primary plot: stacked bars coloured by institution, faceted by domain,
#' X axis = time at chosen granularity.
#'
#' @param df               data frame from load_all_data()
#' @param granularity      "year" | "quarter" | "month"
#' @param top_n_inst       keep this many institutions individually; rest → "Other Norway"
#' @param domains          character vector of domain_label values to include (NULL = all)
#' @param min_year         drop entries before this year
#' @param min_domain_entries exclude domains with fewer total entries than this threshold
plot_time_by_domain <- function(df,
                                granularity        = c("year", "quarter", "month"),
                                top_n_inst         = 12L,
                                domains            = NULL,
                                min_year           = 2000L,
                                min_domain_entries = MIN_DOMAIN_ENTRIES) {
  granularity <- match.arg(granularity)
 
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
 
  # ── Institution lumping ──────────────────────────────────────────────────
  top_inst <- d %>%
    count(institution, sort = TRUE) %>%
    slice_head(n = top_n_inst) %>%
    pull(institution)
 
  d <- d %>%
    mutate(institution = if_else(institution %in% top_inst,
                                 institution, "Other Norway"))
 
  # ── Aggregate ────────────────────────────────────────────────────────────
  counts <- d %>%
    count(domain_label, time_val, time_label, institution, name = "n") %>%
    # Order institutions: named first (by total), "Other Norway" last
    mutate(institution = fct_reorder(institution, n, .fun = sum) %>%
             fct_relevel("Other Norway", after = 0L))
 
  # Build palette keyed on the institutions present in this subset
  pal <- make_inst_palette(levels(counts$institution))
 
  # ── One representative x-axis label per major tick ───────────────────────
  x_labels <- counts %>%
    distinct(time_val, time_label) %>%
    arrange(time_val)
 
  ggplot(counts, aes(x = time_val, y = n, fill = institution)) +
    geom_col(width = if (granularity == "year") 0.7 else
                     if (granularity == "quarter") 0.22 else 0.07) +
    facet_wrap(~domain_label, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = pal, name = "Institution") +
    scale_x_continuous(
      breaks = x_labels$time_val,
      labels = x_labels$time_label
    ) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .05))) +
    guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
    labs(
      title    = "**Norwegian submissions to EBI repositories**",
      subtitle = glue(
        "Coloured by institution · faceted by repository · granularity: {granularity}"
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
           quarter, month, institution, affiliation, country, email) %>%
    readr::write_csv(file.path(OUT_DIR, "norwegian_entries.csv"))
 
  message("Done ✓  Files in ", OUT_DIR)
}
 
 
# =============================================================================
# 5.  Shiny app wrapper
# =============================================================================
 
shiny_app <- function(df) {
  require(shiny)
  require(shinythemes)
  require(DT)
 
  domain_choices <- sort(unique(df$domain_label))
  year_range     <- range(df$year, na.rm = TRUE)
 
  ui <- fluidPage(
    theme = shinytheme("flatly"),
    titlePanel("Norwegian EBI Submissions Dashboard"),
    sidebarLayout(
      sidebarPanel(
        width = 3,
 
        selectInput("granularity", "Time granularity",
                    choices  = c("Year" = "year", "Quarter" = "quarter", "Month" = "month"),
                    selected = "year"),
 
        sliderInput("top_n_inst", "Top N institutions (rest → Other Norway)",
                    min = 3, max = 30, value = 12, step = 1),
 
        sliderInput("min_year", "From year",
                    min   = year_range[1],
                    max   = year_range[2],
                    value = max(year_range[1], year_range[2] - 10L),
                    step  = 1L,
                    sep   = ""),

        sliderInput("min_domain_entries",
                    "Min entries per repository (hide smaller ones)",
                    min = 1, max = 100, value = MIN_DOMAIN_ENTRIES, step = 1),
 
        checkboxGroupInput("domains", "Repositories",
                           choices  = domain_choices,
                           selected = domain_choices),
 
        hr(),
        p(em(paste("Data fetched:", Sys.Date()))),
        p(em(paste(nrow(df), "Norwegian entries total")))
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
    output$main_plot <- renderPlot({
      plot_time_by_domain(
        df,
        granularity        = input$granularity,
        top_n_inst         = input$top_n_inst,
        domains            = input$domains,
        min_year           = input$min_year,
        min_domain_entries = input$min_domain_entries
      )
    }, res = 120)
 
    output$entry_table <- DT::renderDataTable({
      # Respect the same domain threshold in the table
      keep <- df %>%
        filter(domain_label %in% input$domains,
               !is.na(year), year >= input$min_year) %>%
        count(domain_label, name = "total") %>%
        filter(total >= input$min_domain_entries) %>%
        pull(domain_label)
      df %>%
        filter(domain_label %in% keep,
               !is.na(year), year >= input$min_year) %>%
        select(Repository = domain_label, Accession = accession,
               Title = title, Date = date, Institution = institution) %>%
        arrange(desc(Date))
    }, options = list(pageLength = 10, scrollX = TRUE), filter = "top")
  }
 
  shinyApp(ui, server)
}
 
 
# =============================================================================
# 6.  Entry point
# =============================================================================
 
main <- function() {
  df <- load_all_data()
  message(glue("Loaded {nrow(df)} Norwegian entries across {n_distinct(df$domain)} domains"))
 
  if (nrow(df) == 0) {
    message("No data found – have you run fetch_ebi_data.py yet?")
    return(invisible(NULL))
  }
 
  if (RUN_AS_SHINY) {
    require(DT)
    shiny_app(df)
  } else {
    save_plots(df)
  }
}
 
# Null-coalescing operator.
# Safe when a[[1]] is a vector of length > 1 (multi-value API fields).
`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0) return(b)
  first <- a[[1]]
  if (length(first) == 0 || all(is.na(first))) return(b)
  a
}
 
main()
