library(shiny)

shinyUI(fluidPage(

  titlePanel("Brain cell gene expression data"),

  sidebarLayout(
    sidebarPanel(
      textInput(inputId = "gene", label = h3("Enter gene symbol:"), value = "Plp1"),
      selectInput("data_type", h3("Data Source:"),
            c("Sharma et al. RNA Mouse (RPKM)" = "Sharma_RNA",
              "Sharma et al. Proteomics Mouse (LFQ)" = "Sharma_Proteomics",
              "Darmanis et al. Single Cell RNA Human Cortex (Norm. CPM)" = "Darmanis_RNA",
              "Tasic et al. Single Cell RNA V1 Mouse (Norm. CPM)" = "Tasic_RNA",
              "Zeisel et al. Single Cell RNA Mouse (Norm. CPM)" = "Zeisel_RNA",
              "Marques et al. Single Cell Oligodendrocyte RNA Mouse (Norm. CPM)" = "Marques_RNA")),
      checkboxInput(inputId = "log", label = "Log2 transform?", value = FALSE)
    ),

    mainPanel(
      plotOutput("cellPlot")
    )
  ),

  h5("References:"),
  a("Sharma K, et al. Nat Neurosci 2015; doi:10.1038/nn.4160",
    href = "http://www.nature.com/neuro/journal/v18/n12/abs/nn.4160.html"),
  br(),
  a("Darmanis S, et al. PNAS 2015; doi:10.1073/pnas.1507125112",
    href = "http://www.pnas.org/content/112/23/7285"),
  br(),
  a("Marques S, et al. Science 2016; doi:10.1126/science.aaf6463",
    href = "http://science.sciencemag.org/content/352/6291/1326.full"),
  br(),
  a("Tasic B, et al. Nat Neurosci 2016; doi:10.1038/nn.4216",
    href = "http://www.nature.com/neuro/journal/v19/n2/full/nn.4216.html"),
  br(),
  a("Zeisel A, et al. Science 2015; doi:10.1126/science.aaa1934",
    href = "http://science.sciencemag.org/content/347/6226/1138.full"),
  h5("Notes:"),
  div("- Error bars represent the standard error of the mean."),
  div("- Certain sub-cell types are averaged for visualization in Zeisel et al."),
  div("- At least 15 samples per cell type were required for the single cell RNA-seq data."),
  div("- Data is averaged across regions for some cell types in Marques et al. and Zeisel et al."),
  div("- We attempt to convert human to mouse gene symbols (and vv.) through capitilization."),
  h5("All source code is available:"),
  a("GitHub", href = "https://github.com/andymckenzie/brain_gene_expression")

))
