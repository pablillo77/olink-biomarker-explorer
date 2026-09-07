# app.R — Interactive Proteomic Biomarker Explorer (v0: load data, show shape)
#
# Shiny in one sentence: `ui` describes the page; `server` describes how
# outputs are computed from inputs; Shiny re-runs only the pieces whose
# inputs changed (reactivity). Compare to a Python dashboard where you'd
# manually attach callbacks — here the dependency graph is inferred.

library(shiny)
library(dplyr)

# Source every helper in /R so the app and the notebook share one code base.
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

ui <- fluidPage(
  titlePanel("Olink NPX Biomarker Explorer"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "dataset", "Dataset",
        choices = c("Olink demo 1 (npx_data1)" = "npx_data1",
                    "Olink demo 2 (npx_data2)" = "npx_data2")
      ),
      helpText("Synthetic example data shipped with OlinkAnalyze ",
               "(randomly generated, not patient data).")
    ),
    mainPanel(
      h4("Dataset shape"),
      verbatimTextOutput("shape"),
      h4("First rows"),
      tableOutput("head")
    )
  )
)

server <- function(input, output, session) {

  # reactive(): a cached expression that re-evaluates only when input$dataset
  # changes. Call it like a function — npx() — wherever you need the data.
  npx <- reactive({
    load_npx_demo(input$dataset)
  })

  # render*(): each output declares what it depends on simply by using it.
  output$shape <- renderPrint({
    str(npx_summary(npx()))
  })

  output$head <- renderTable({
    npx() |>
      select(SampleID, OlinkID, Assay, Panel, NPX, LOD, QC_Warning,
             Treatment, Site, Time) |>
      head(10)
  })
}

shinyApp(ui, server)
