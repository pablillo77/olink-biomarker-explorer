# R/mod_pca.R
# Same module pattern as mod_qc.R: ui + server pair sharing a namespace.
# This module takes TWO reactives in: the dataset and the excluded-sample
# list returned by the QC module — showing how modules compose.

mod_pca_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "PCA settings",
      selectInput(ns("color_by"), "Color by",
                  choices = c("Treatment", "Site", "Time"), selected = "Treatment"),
      helpText("Each protein is scaled (mean 0, sd 1) before PCA so that ",
               "no single high-variance protein dominates the components.")
    ),
    bslib::layout_columns(
      col_widths = c(6, 6),
      bslib::value_box(title = "Variance explained (PC1)",
                        value = textOutput(ns("var_pc1")), theme = "primary"),
      bslib::value_box(title = "Variance explained (PC2)",
                        value = textOutput(ns("var_pc2")), theme = "primary")
    ),
    bslib::card(
      bslib::card_header("PC1 vs PC2"),
      plotly::plotlyOutput(ns("pca_plot"), height = "460px")
    ),
    bslib::card(
      bslib::card_header("Top proteins loading on PC1 / PC2"),
      DT::DTOutput(ns("loadings_table"))
    )
  )
}

mod_pca_server <- function(id, npx, excluded) {
  moduleServer(id, function(input, output, session) {

    # Both npx() and excluded() are reactives from other parts of the app;
    # this reactive recomputes whenever EITHER changes.
    pca_result <- reactive({
      run_pca(npx(), excluded = excluded())
    })

    output$var_pc1 <- renderText({
      scales::percent(pca_result()$variance_explained[1], accuracy = 0.1)
    })
    output$var_pc2 <- renderText({
      scales::percent(pca_result()$variance_explained[2], accuracy = 0.1)
    })

    output$pca_plot <- plotly::renderPlotly({
      s <- pca_result()$scores
      p <- ggplot2::ggplot(s, ggplot2::aes(.data$PC1, .data$PC2,
                                            color = .data[[input$color_by]],
                                            text = .data$SampleID)) +
        ggplot2::geom_point(size = 2, alpha = 0.8) +
        ggplot2::theme_minimal() +
        ggplot2::labs(x = "PC1", y = "PC2", color = input$color_by)
      plotly::ggplotly(p, tooltip = c("text", "colour"))
    })

    output$loadings_table <- DT::renderDT({
      pca_result()$loadings |>
        dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) round(x, 3))) |>
        DT::datatable(options = list(pageLength = 8), rownames = FALSE)
    })
  })
}
