# R/mod_diffab.R
# Differential abundance tab: mixed-model test per protein, BH-corrected,
# shown as a volcano plot, with a downloadable results table.

mod_diffab_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "Differential abundance",
      sliderInput(ns("fdr_threshold"), "FDR (adjusted p-value) cutoff",
                  min = 0.01, max = 0.25, value = 0.05, step = 0.01),
      helpText("Model per protein: NPX ~ Treatment + (1 | Subject) — a ",
               "random intercept per Subject accounts for repeated ",
               "measures (each Subject contributes samples at multiple ",
               "timepoints), avoiding pseudoreplication."),
      hr(),
      downloadButton(ns("download_results"), "Download results (CSV)")
    ),
    bslib::value_box(
      title = "Significant proteins (FDR-corrected)",
      value = textOutput(ns("n_significant")),
      theme = "success"
    ),
    bslib::card(
      bslib::card_header("Volcano plot"),
      plotly::plotlyOutput(ns("volcano"), height = "460px")
    ),
    bslib::card(
      bslib::card_header("Results table"),
      DT::DTOutput(ns("results_table"))
    )
  )
}

mod_diffab_server <- function(id, npx, excluded) {
  moduleServer(id, function(input, output, session) {

    diffab_result <- reactive({
      run_differential_abundance(npx(), excluded = excluded())
    })

    flagged <- reactive({
      diffab_result() |>
        dplyr::mutate(significant = .data$Adjusted_pval < input$fdr_threshold)
    })

    output$n_significant <- renderText({
      sum(flagged()$significant, na.rm = TRUE)
    })

    output$volcano <- plotly::renderPlotly({
      d <- flagged() |>
        dplyr::mutate(neg_log10_p = -log10(.data$p.value))
      p <- ggplot2::ggplot(d, ggplot2::aes(.data$estimate, .data$neg_log10_p,
                                            color = .data$significant,
                                            text = .data$Assay)) +
        ggplot2::geom_point(alpha = 0.8) +
        ggplot2::scale_color_manual(values = c(`TRUE` = "#d9534f", `FALSE` = "grey70")) +
        ggplot2::theme_minimal() +
        ggplot2::labs(x = "Effect size (Treatment)", y = "-log10(p-value)",
                       color = "FDR-significant")
      plotly::ggplotly(p, tooltip = c("text", "x", "y"))
    })

    output$results_table <- DT::renderDT({
      flagged() |>
        dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) signif(x, 4))) |>
        DT::datatable(options = list(pageLength = 10, order = list(list(3, "asc"))),
                       rownames = FALSE) |>
        DT::formatStyle("significant", target = "row",
                         backgroundColor = DT::styleEqual(TRUE, "#e6f4ea"))
    })

    output$download_results <- downloadHandler(
      filename = function() paste0("differential_abundance_", Sys.Date(), ".csv"),
      content = function(file) {
        utils::write.csv(flagged(), file, row.names = FALSE)
      }
    )
  })
}
