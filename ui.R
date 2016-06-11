library(shiny)

shinyUI(fluidPage(

  titlePanel("Mouse brain cell type data"),

  sidebarLayout(
    sidebarPanel(
      textInput(inputId = "gene", label = h3("Enter mouse gene symbol:"), value = "Plp1"),
      selectInput("data_type", h3("Data Source:"),
            c("Sharma et al. RNA (RPKM)" = "Sharma_RNA",
              "Sharma et al. Proteomics (LFQ)" = "Sharma_Proteomics")),
      checkboxInput(inputId = "log", label = "Log2 transform?", value = FALSE)
    ),

    mainPanel(
      plotOutput("cellPlot")
    )
  ),

  h5("Reference:"),
  h5("Sharma K, Schmitt S, Bergner CG, et al."),
  h5("Cell type- and brain region-resolved mouse brain proteome. "),
  h5("Nat Neurosci. 2015;18(12):1819-31."),
  a("doi:10.1038/nn.4160", href = "http://www.nature.com/neuro/journal/v18/n12/abs/nn.4160.html")
  a("Source Code", href = "http://www.nature.com/neuro/journal/v18/n12/abs/nn.4160.html")

))
