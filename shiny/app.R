#!/usr/bin/env Rscript
# =============================================================================
# shiny/app.R  –  Shinylive-compatible standalone dashboard
# =============================================================================
# This file is the entry point for the static Shinylive / GitHub Pages build.
# It is deliberately self-contained and does NOT share code with
# R/plot_norwegian_data.R.
#
# Data source
# -----------
# At CI build time (before shinylive::export() is called) the GitHub Actions
# workflow copies output/norwegian_entries.csv → shiny/data/norwegian_entries.csv
# so the file is bundled inside the exported static site.
# The CSV is produced by R/plot_norwegian_data.R and already contains a fully
# resolved `institution` column – no JSON parsing or fuzzy matching needed here.
#
# WebR / Shinylive package notes
# ------------------------------
# - ggtext is intentionally excluded: element_markdown() is replaced with
#   element_text() because ggtext's C dependencies (systemfonts / textshaping)
#   are not reliably available in WebR.
# - stringdist and jsonlite are not needed (no data processing at runtime).
# - All remaining packages are available as pre-compiled WebR binaries.
# =============================================================================

library(shiny)
library(shinythemes)
library(ggplot2)
library(dplyr)
library(forcats)
library(scales)
library(RColorBrewer)
library(glue)
library(DT)
library(readr)
library(lubridate)

# ── Load data (bundled at export time) ────────────────────────────────────────
df <- readr::read_csv(
  "data/norwegian_entries.csv",
  col_types = cols(
    domain        = col_character(),
    domain_label  = col_character(),
    accession     = col_character(),
    title         = col_character(),
    date          = col_date(),
    year          = col_integer(),
    quarter       = col_integer(),
    month         = col_integer(),
    institution   = col_character(),
    affiliation   = col_character(),
    country       = col_character(),
    email         = col_character()
  ),
  show_col_types = FALSE
)

# ── Helpers ───────────────────────────────────────────────────────────────────

#' Distinct, perceptually-separated palette; "Other Norway" always grey.
make_inst_palette <- function(institutions) {
  institutions <- as.character(institutions)
  non_other    <- sort(setdiff(institutions, "Other Norway"))
  n            <- length(non_other)

  base_pal <- if (n <= 8) {
    RColorBrewer::brewer.pal(max(3L, n), "Set2")
  } else if (n <= 12) {
    RColorBrewer::brewer.pal(12L, "Set3")
  } else {
    colorRampPalette(
      c(RColorBrewer::brewer.pal(9L, "Set1"),
        RColorBrewer::brewer.pal(8L, "Set2"),
        RColorBrewer::brewer.pal(8L, "Dark2"))
    )(n)
  }

  pal <- setNames(base_pal[seq_along(non_other)], non_other)
  if ("Other Norway" %in% institutions) pal["Other Norway"] <- "#AAAAAA"
  pal
}

#' Shared ggplot2 theme.
#' Uses element_text (not element_markdown) for WebR compatibility.
theme_nor <- function() {
  theme_classic(base_size = 13) +
    theme(
      plot.title       = element_text(face = "bold", size = 15),
      plot.subtitle    = element_text(colour = "grey40"),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey92", colour = NA)
    )
}

MIN_DOMAIN_ENTRIES <- 10L   # configurable default; exposed as a Shiny slider

#' Stacked-bar chart coloured by institution, faceted by repository.
plot_time_by_domain <- function(df,
                                granularity        = c("year", "quarter", "month"),
                                top_n_inst         = 12L,
                                domains            = NULL,
                                min_year           = 2000L,
                                min_domain_entries = MIN_DOMAIN_ENTRIES) {
  granularity <- match.arg(granularity)

  d <- df
  if (!is.null(domains)) d <- d |> filter(domain_label %in% domains)
  d <- d |> filter(!is.na(date), year >= min_year, year <= year(Sys.Date()))

  # ── Drop domains below the entry threshold ────────────────────────────────
  keep_domains <- d |>
    count(domain_label, name = "total") |>
    filter(total >= min_domain_entries) |>
    pull(domain_label)
  d <- d |> filter(domain_label %in% keep_domains)

  if (nrow(d) == 0) {
    return(
      ggplot() +
        labs(title = "No domains meet the minimum entry threshold") +
        theme_void()
    )
  }

  # ── Time axis ────────────────────────────────────────────────────────────
  d <- switch(granularity,
    year    = d |> mutate(time_val   = as.integer(year),
                          time_label = as.character(year)),
    quarter = d |> mutate(time_val   = year + (quarter - 1) / 4,
                          time_label = paste0(year, "\nQ", quarter)),
    month   = d |> mutate(time_val   = year + (month - 1) / 12,
                          time_label = format(date, "%Y\n%b"))
  )

  # ── Institution lumping ──────────────────────────────────────────────────
  top_inst <- d |> count(institution, sort = TRUE) |>
    slice_head(n = top_n_inst) |> pull(institution)

  d <- d |>
    mutate(institution = if_else(institution %in% top_inst,
                                 institution, "Other Norway"))

  # ── Aggregate ────────────────────────────────────────────────────────────
  counts <- d |>
    count(domain_label, time_val, time_label, institution, name = "n") |>
    mutate(institution = fct_reorder(institution, n, .fun = sum) |>
             fct_relevel("Other Norway", after = 0L))

  pal      <- make_inst_palette(levels(counts$institution))
  x_labels <- counts |> distinct(time_val, time_label) |> arrange(time_val)

  ggplot(counts, aes(x = time_val, y = n, fill = institution)) +
    geom_col(width = if (granularity == "year")    0.70 else
                     if (granularity == "quarter") 0.22 else 0.07) +
    facet_wrap(~domain_label, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = pal, name = "Institution") +
    scale_x_continuous(
      breaks = x_labels$time_val,
      labels = x_labels$time_label
    ) +
    scale_y_continuous(
      labels = scales::comma,
      expand = expansion(mult = c(0, .05))
    ) +
    guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
    labs(
      title    = "Norwegian submissions to EBI repositories",
      subtitle = glue(
        "Coloured by institution \u00b7 faceted by repository \u00b7 granularity: {granularity}"
      ),
      x = NULL,
      y = "Number of entries"
    ) +
    theme_nor()
}

# ── Derived constants ─────────────────────────────────────────────────────────
domain_choices <- sort(unique(df$domain_label))
year_range     <- range(df$year, na.rm = TRUE)
latest_date    <- max(df$date,   na.rm = TRUE)

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel("Norwegian EBI Submissions Dashboard"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      selectInput(
        "granularity", "Time granularity",
        choices  = c("Year" = "year", "Quarter" = "quarter", "Month" = "month"),
        selected = "year"
      ),

      sliderInput(
        "top_n_inst", "Top N institutions (rest \u2192 Other Norway)",
        min = 3, max = 30, value = 12, step = 1
      ),

      sliderInput(
        "min_year", "From year",
        min   = year_range[1],
        max   = year_range[2],
        value = max(year_range[1], year_range[2] - 10L),
        step  = 1L,
        sep   = ""
      ),

      sliderInput(
        "min_domain_entries",
        "Min entries per repository (hide smaller ones)",
        min = 1, max = 100, value = MIN_DOMAIN_ENTRIES, step = 1
      ),

      checkboxGroupInput(
        "domains", "Repositories",
        choices  = domain_choices,
        selected = domain_choices
      ),

      hr(),
      p(em(paste("Latest entry:", format(latest_date, "%Y-%m-%d")))),
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

# ── Server ────────────────────────────────────────────────────────────────────
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
    keep <- df |>
      filter(domain_label %in% input$domains,
             !is.na(year), year >= input$min_year) |>
      count(domain_label, name = "total") |>
      filter(total >= input$min_domain_entries) |>
      pull(domain_label)
    df |>
      filter(domain_label %in% keep,
             !is.na(year), year >= input$min_year) |>
      select(
        Repository   = domain_label,
        Accession    = accession,
        Title        = title,
        Date         = date,
        Institution  = institution,
        Email        = email
      ) |>
      arrange(desc(Date))
  }, options = list(pageLength = 10, scrollX = TRUE), filter = "top")

}

# ── Launch ────────────────────────────────────────────────────────────────────
shinyApp(ui, server)
