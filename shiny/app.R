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

#' Grouped bar chart, one bar per institution per time period, faceted by domain.
#'
#' @param selected_inst  Character vector of institution names to show after
#'   lumping.  NULL means show all (used while the checkbox UI is re-rendering).
plot_time_by_domain <- function(df,
                                granularity        = c("year", "quarter", "month"),
                                top_n_inst         = 8L,
                                domains            = NULL,
                                year_range         = NULL,
                                selected_inst      = NULL,
                                min_domain_entries = MIN_DOMAIN_ENTRIES) {
  granularity <- match.arg(granularity)

  d <- df
  if (!is.null(domains))                d <- d |> filter(domain_label %in% domains)
  if (!is.null(year_range) && length(year_range) == 2)
    d <- d |> filter(year >= year_range[1], year <= year_range[2])
  d <- d |> filter(!is.na(date), year <= year(Sys.Date()))

  # Drop repositories below entry threshold
  keep_domains <- d |>
    count(domain_label, name = "total") |>
    filter(total >= min_domain_entries) |>
    pull(domain_label)
  d <- d |> filter(domain_label %in% keep_domains)

  if (nrow(d) == 0)
    return(ggplot() + labs(title = "No data for current filters") + theme_void())

  # Time axis
  d <- switch(granularity,
    year    = d |> mutate(time_val   = as.integer(year),
                          time_label = as.character(year)),
    quarter = d |> mutate(time_val   = year + (quarter - 1) / 4,
                          time_label = paste0(year, "\nQ", quarter)),
    month   = d |> mutate(time_val   = year + (month - 1) / 12,
                          time_label = format(date, "%Y\n%b"))
  )

  # Lump institutions outside top N into "Other Norway"
  top_inst <- d |> count(institution, sort = TRUE) |>
    slice_head(n = top_n_inst) |> pull(institution)
  d <- d |> mutate(institution = if_else(institution %in% top_inst,
                                         institution, "Other Norway"))

  # Apply institution checkbox filter (NULL = show all, used during re-render)
  if (!is.null(selected_inst) && length(selected_inst) > 0)
    d <- d |> filter(institution %in% selected_inst)

  if (nrow(d) == 0)
    return(ggplot() + labs(title = "No institutions selected") + theme_void())

  counts <- d |>
    count(domain_label, time_val, time_label, institution, name = "n") |>
    mutate(institution = fct_reorder(institution, n, .fun = sum) |>
             fct_relevel("Other Norway", after = 0L))

  pal      <- make_inst_palette(levels(counts$institution))
  x_labels <- counts |> distinct(time_val, time_label) |> arrange(time_val)

  # Per-granularity group width: bars within a time unit span this fraction
  bar_width <- switch(granularity, year = 0.8, quarter = 0.22, month = 0.07)

  ggplot(counts, aes(x = time_val, y = n, fill = institution)) +
    geom_col(
      position = position_dodge2(padding = 0.1, preserve = "single"),
      width    = bar_width
    ) +
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
        "Grouped by institution · faceted by repository · {granularity}"
      ),
      x = NULL, y = "Number of entries"
    ) +
    theme_nor()
}

# ── Derived constants ─────────────────────────────────────────────────────────
domain_choices <- sort(unique(df$domain_label))
year_range_data <- range(df$year, na.rm = TRUE)
latest_date     <- max(df$date, na.rm = TRUE)

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

      sliderInput(
        "top_n_inst", "Top N institutions",
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
      # Dynamic institution list: rebuilt when top_n or domain/year filters change.
      # All top-N institutions are ticked by default; user may untick to hide.
      uiOutput("inst_checkbox"),

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

  # Top N institutions for the current domain + year filter.
  # Used both to populate the checkbox and to drive the lump in the plot.
  top_institutions <- reactive({
    req(input$top_n_inst, input$year_range)
    d <- df
    if (!is.null(input$domains))
      d <- d |> filter(domain_label %in% input$domains)
    d <- d |> filter(year >= input$year_range[1], year <= input$year_range[2])

    top_inst <- d |>
      count(institution, sort = TRUE) |>
      slice_head(n = input$top_n_inst) |>
      pull(institution)

    # Always include "Other Norway" so the user can toggle it
    unique(c(top_inst, "Other Norway"))
  })

  # Rebuild the checkbox whenever the institution list changes.
  # All entries ticked by default; selections reset when the list changes.
  output$inst_checkbox <- renderUI({
    inst <- top_institutions()
    checkboxGroupInput(
      "selected_inst", "Show institutions",
      choices  = inst,
      selected = inst
    )
  })

  output$main_plot <- renderPlot({
    req(input$year_range)
    plot_time_by_domain(
      df,
      granularity        = input$granularity,
      top_n_inst         = input$top_n_inst,
      domains            = input$domains,
      year_range         = input$year_range,
      selected_inst      = input$selected_inst,
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
        Email       = email
      ) |>
      arrange(desc(Date))
  }, options = list(pageLength = 10, scrollX = TRUE), filter = "top")
}

shinyApp(ui, server)
