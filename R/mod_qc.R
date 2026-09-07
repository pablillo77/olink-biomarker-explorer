# R/mod_qc.R
#
# SHINY MODULES, explained: a module is just a pair of functions —
# `mod_qc_ui(id)` and `mod_qc_server(id, npx)` — that share a private
# namespace created by `NS(id)`. Every input/output id gets wrapped in
# ns(...), so `ns("threshold")` might resolve to "qc-threshold" internally.
# That means you could mount two QC modules side by side (e.g. comparing
# two datasets) and their sliders/tables would never collide.
#
# The server function *returns* a reactive (here, the vector of excluded
# SampleIDs) so app.R — or a future PCA/differential-abundance module —
# can consume "samples the user excluded" without knowing anything about
# how the QC tab is built. That's the module boundary: inputs in via
# arguments, outputs out via a returned reactive.

mod_qc_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "QC settings",
      sliderInput(
        ns("threshold"), "Flag samples with QC-warning rate above",
        min = 0, max = 1, value = 0.10, step = 0.01
      ),
      checkboxInput(ns("hide_controls"), "Hide Olink control samples", TRUE),
      hr(),
      uiOutput(ns("exclude_ui")),
      hr(),
      helpText("Excluding samples here removes them from every other tab.")
    ),
    bslib::layout_columns(
      col_widths = c(4, 4, 4),
      bslib::value_box(
        title = "Samples flagged",
        value = textOutput(ns("n_flagged")),
        showcase = bsicons::bs_icon("exclamation-triangle"),
        theme = "warning"
      ),
      bslib::value_box(
        title = "NPX values below LOD",
        value = textOutput(ns("pct_below_lod")),
        showcase = bsicons::bs_icon("graph-down"),
        theme = "secondary"
      ),
      bslib::value_box(
        title = "Samples excluded",
        value = textOutput(ns("n_excluded")),
        showcase = bsicons::bs_icon("x-circle"),
        theme = "danger"
      )
    ),
    bslib::card(
      bslib::card_header("NPX distribution per sample (first 40 shown)"),
      plotly::plotlyOutput(ns("dist_plot"), height = "380px")
    ),
    bslib::card(
      bslib::card_header("Per-sample QC summary"),
      DT::DTOutput(ns("qc_table"))
    )
  )
}

mod_qc_server <- function(id, npx) {
  # `npx` is a *reactive* passed in from app.R (e.g. the dataset the user
  # picked). Modules receive reactives as plain arguments and call them,
  # npx(), wherever the current value is needed — same rule as anywhere else.
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    qc_summary <- reactive({
      s <- sample_qc_summary(npx())
      if (input$hide_controls) s <- dplyr::filter(s, !.data$is_control)
      s
    })

    flagged_ids <- reactive({
      flag_failing_samples(npx(), threshold = input$threshold)
    })

    # reactiveVal holds state that persists across reactive ticks but isn't
    # derived from other reactives — here, the user's manual exclusion list.
    excluded_ids <- reactiveVal(character(0))

    output$exclude_ui <- renderUI({
      choices <- qc_summary()$SampleID
      selectizeInput(
        ns("manual_exclude"), "Manually exclude samples",
        choices = choices, multiple = TRUE,
        selected = intersect(isolate(excluded_ids()), choices)
      )
    })

    observeEvent(input$manual_exclude, {
      excluded_ids(input$manual_exclude)
    }, ignoreNULL = FALSE)

    output$n_flagged     <- renderText(length(flagged_ids()))
    output$n_excluded    <- renderText(length(excluded_ids()))
    output$pct_below_lod <- renderText({
      scales::percent(missingness_overview(npx())$pct_below_lod, accuracy = 0.1)
    })

    output$dist_plot <- plotly::renderPlotly({
      d <- npx()
      keep <- unique(d$SampleID)[1:min(40, dplyr::n_distinct(d$SampleID))]
      p <- d |>
        dplyr::filter(.data$SampleID %in% keep) |>
        ggplot2::ggplot(ggplot2::aes(.data$SampleID, .data$NPX,
                                      fill = .data$QC_Warning)) +
        ggplot2::geom_boxplot(outlier.size = 0.5) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, size = 6)) +
        ggplot2::labs(x = NULL, y = "NPX", fill = "QC")
      plotly::ggplotly(p) |> plotly::layout(legend = list(orientation = "h"))
    })

    output$qc_table <- DT::renderDT({
      qc_summary() |>
        dplyr::mutate(
          pct_warning   = scales::percent(.data$pct_warning, accuracy = 0.1),
          pct_below_lod = scales::percent(.data$pct_below_lod, accuracy = 0.1),
          flagged       = .data$SampleID %in% flagged_ids()
        ) |>
        DT::datatable(
          options = list(pageLength = 10, order = list(list(2, "desc"))),
          rownames = FALSE
        ) |>
        DT::formatStyle("flagged", target = "row",
                         backgroundColor = DT::styleEqual(TRUE, "#fff3cd"))
    })

    # The module's public contract: whoever calls mod_qc_server() gets this
    # reactive back and can plug it straight into filtering for other tabs.
    reactive(excluded_ids())
  })
}
