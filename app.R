# app.R — Interactive Proteomic Biomarker Explorer
#
# Shiny in one sentence: `ui` describes the page; `server` describes how
# outputs are computed from inputs; Shiny re-runs only the pieces whose
# inputs changed (reactivity). Compare to a Python dashboard where you'd
# manually attach callbacks — here the dependency graph is inferred.
#
# app.R stays a thin shell: it loads data and mounts modules. All the real
# logic lives in R/ so it's independently testable (see tests/testthat/).

library(shiny)
library(bslib)
library(dplyr)

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

ui <- page_navbar(
  title = "Olink NPX Biomarker Explorer",
  theme = bs_theme(version = 5, bootswatch = "flatly", base_font = font_google("Inter")),
  sidebar = sidebar(
    title = "Data",
    selectInput(
      "dataset", "Dataset",
      choices = c("Olink demo 1 (npx_data1)" = "npx_data1",
                  "Olink demo 2 (npx_data2)" = "npx_data2")
    ),
    helpText("Synthetic example data shipped with OlinkAnalyze ",
             "(randomly generated, not patient data).")
  ),
  nav_panel("Overview", verbatimTextOutput("shape")),
  nav_panel("Quality Control", mod_qc_ui("qc"))
)

server <- function(input, output, session) {

  # reactive(): a cached expression that re-evaluates only when input$dataset
  # changes. Call it like a function — npx() — wherever you need the data.
  npx <- reactive({
    load_npx_demo(input$dataset)
  })

  output$shape <- renderPrint({
    str(npx_summary(npx()))
  })

  # Mount the QC module under namespace "qc". It gets the raw dataset in,
  # and hands back a reactive of excluded SampleIDs — wired up to nothing
  # yet, but this is exactly what the PCA/differential-abundance tabs will
  # consume once they exist.
  excluded_samples <- mod_qc_server("qc", npx = npx)
}

shinyApp(ui, server)
