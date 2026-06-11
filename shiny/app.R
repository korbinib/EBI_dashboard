#!/usr/bin/env Rscript
# =============================================================================
# shiny/app.R  –  Shinylive-compatible standalone dashboard
# =============================================================================
# Data source: shiny/data/norwegian_entries.csv bundled at CI export time.
# The CSV is produced by R/plot_norwegian_data.R (institution already resolved).
#
# WebR / Shinylive package notes
# - ggtext excluded (C deps not reliably in WebR); element_text() used instead.
# - stringdist and jsonlite not needed here.
#
# NOTE: the plot style (make_inst_palette, theme_nor, the grouped/dodged geom)
# is intentionally duplicated from R/plot_norwegian_data.R so this app stays
# self-contained for the shinylive export.  Keep the two copies in sync.
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
    broker        = col_character(),
    affiliation   = col_character(),
    country       = col_character(),
    email         = col_character()
  ),
  show_col_types = FALSE
)

# ── Helpers ───────────────────────────────────────────────────────────────────

# "Other Norway", "Other", and "Non-ENA" are pinned to grey so they recede.
make_inst_palette <- function(values) {
  values    <- as.character(values)
  grey_keys <- c("Other Norway", "Other", "Non-ENA")
  non_grey  <- sort(setdiff(values, grey_keys))
  n         <- length(non_grey)

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

  pal <- setNames(base_pal[seq_along(non_grey)], non_grey)
  for (g in intersect(grey_keys, values)) pal[g] <- "#AAAAAA"
  pal
}

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

MIN_DOMAIN_ENTRIES <- 10L

#' Grouped bar chart, faceted by domain.  Colour dimension is controlled by
#' `color_by` ("institution" or "broker").  Pre-lumped data is accepted when
#' top_n_inst = 999 (lumping already done in the server reactive).
plot_time_by_domain <- function(df,
                                granularity        = c("year", "quarter", "month"),
                                top_n_inst         = 8L,
                                domains            = NULL,
                                year_range         = NULL,
                                selected_fill      = NULL,
                                min_domain_entries = MIN_DOMAIN_ENTRIES,
                                color_by           = c("institution", "broker")) {
  granularity <- match.arg(granularity)
  color_by    <- match.arg(color_by)

  fill_col    <- color_by
  other_label <- if (color_by == "institution") "Other Norway" else "Other"
  legend_name <- if (color_by == "institution") "Institution" else "ENA Broker / Center"

  d <- df
  if (!is.null(domains))
    d <- d |> filter(domain_label %in% domains)
  if (!is.null(year_range) && length(year_range) == 2)
    d <- d |> filter(year >= year_range[1], year <= year_range[2])
  d <- d |> filter(!is.na(date), year <= year(Sys.Date()))

  keep_domains <- d |>
    count(domain_label, name = "total") |>
    filter(total >= min_domain_entries) |>
    pull(domain_label)
  d <- d |> filter(domain_label %in% keep_domains)

  if (nrow(d) == 0)
    return(ggplot() + labs(title = "No data for current filters") + theme_void())

  d <- switch(granularity,
    year    = d |> mutate(time_val   = as.integer(year),
                          time_label = as.character(year)),
    quarter = d |> mutate(time_val   = year + (quarter - 1) / 4,
                          time_label = paste0(year, "\nQ", quarter)),
    month   = d |> mutate(time_val   = year + (month - 1) / 12,
                          time_label = format(date, "%Y\n%b"))
  )

  # Lump fill column to top-N (skipped when top_n_inst = 999 — pre-lumped).
  if (top_n_inst < 999L) {
    top_vals <- d |> count(.data[[fill_col]], sort = TRUE) |>
      slice_head(n = top_n_inst) |> pull(.data[[fill_col]])
    d <- d |> mutate(across(all_of(fill_col),
                             ~ if_else(.x %in% top_vals, .x, other_label)))
  }

  if (!is.null(selected_fill) && length(selected_fill) > 0)
    d <- d |> filter(.data[[fill_col]] %in% selected_fill)

  if (nrow(d) == 0)
    return(ggplot() + labs(title = "No values selected") + theme_void())

  counts <- d |>
    count(domain_label, time_val, time_label, .data[[fill_col]], name = "n") |>
    rename(fill_val = all_of(fill_col)) |>
    mutate(fill_val = fct_reorder(fill_val, n, .fun = sum) |>
             fct_relevel(other_label, after = 0L))

  pal      <- make_inst_palette(levels(counts$fill_val))
  x_labels <- counts |> distinct(time_val, time_label) |> arrange(time_val)
  bar_width <- switch(granularity, year = 0.8, quarter = 0.22, month = 0.07)

  ggplot(counts, aes(x = time_val, y = n, fill = fill_val)) +
    geom_col(
      position = position_dodge2(padding = 0.1, preserve = "single"),
      width    = bar_width
    ) +
    facet_wrap(~domain_label, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = pal, name = legend_name) +
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
        "Coloured by {color_by} · faceted by repository · {granularity}"
      ),
      x = NULL, y = "Number of entries"
    ) +
    theme_nor()
}

# ── Derived constants ─────────────────────────────────────────────────────────
# Guard against an empty / all-NA CSV: range()/max() on no values return
# Inf/-Inf, which would feed Inf into the year sliderInput and the "Latest entry"
# label.  Fall back to a sane window when there is no dated data yet.
domain_choices  <- sort(unique(df$domain_label))
valid_years     <- df$year[!is.na(df$year)]
this_year       <- as.integer(format(Sys.Date(), "%Y"))
year_range_data <- if (length(valid_years)) range(valid_years) else c(2000L, this_year)
valid_dates     <- df$date[!is.na(df$date)]
latest_date     <- if (length(valid_dates)) max(valid_dates) else Sys.Date()

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

      # Two-handle year range; default: last 10 years to current year
      sliderInput(
        "year_range", "Year range",
        min   = year_range_data[1],
        max   = year_range_data[2],
        value = c(max(year_range_data[1], year_range_data[2] - 10L),
                  year_range_data[2]),
        step  = 1L,
        sep   = ""
      ),

      radioButtons(
        "color_by", "Colour bars by",
        choices  = c("Institution" = "institution",
                     "ENA Broker / Center" = "broker"),
        selected = "institution", inline = TRUE
      ),

      sliderInput(
        "top_n_inst", "Top N values to show",
        min = 3, max = 30, value = 8L, step = 1
      ),

      sliderInput(
        "min_domain_entries",
        "Min entries per repository",
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

  base_df <- reactive({
    req(input$year_range, input$domains)
    df |>
      filter(domain_label %in% input$domains,
             !is.na(year),
             year >= input$year_range[1],
             year <= input$year_range[2])
  })

  # Top-N fill values for the current color_by mode, domain, and year window.
  top_fill_vals <- reactive({
    req(input$top_n_inst, input$color_by)
    col       <- input$color_by
    other_lbl <- if (col == "institution") "Other Norway" else "Other"
    top_vals  <- base_df() |>
      count(.data[[col]], sort = TRUE) |>
      slice_head(n = input$top_n_inst) |>
      pull(.data[[col]])
    unique(c(top_vals, other_lbl))
  })

  # Rebuild checkboxes whenever the fill list changes; all pre-ticked.
  output$fill_checkbox <- renderUI({
    vals  <- top_fill_vals()
    label <- if (input$color_by == "institution") "Show institutions"
             else "Show brokers / centers"
    checkboxGroupInput("selected_fill", label,
                       choices = vals, selected = vals)
  })

  output$main_plot <- renderPlot({
    req(input$year_range, input$color_by)

    col       <- input$color_by
    other_lbl <- if (col == "institution") "Other Norway" else "Other"

    # Pre-lump then apply checkbox filter; pass top_n_inst=999 to skip re-lump.
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
      selected_fill      = NULL,   # already filtered above
      min_domain_entries = input$min_domain_entries
    )
  }, res = 120)

  output$entry_table <- DT::renderDataTable({
    req(input$year_range)
    keep <- df |>
      filter(domain_label %in% input$domains,
             !is.na(year),
             year >= input$year_range[1],
             year <= input$year_range[2]) |>
      count(domain_label, name = "total") |>
      filter(total >= input$min_domain_entries) |>
      pull(domain_label)
    df |>
      filter(domain_label %in% keep,
             !is.na(year),
             year >= input$year_range[1],
             year <= input$year_range[2]) |>
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
