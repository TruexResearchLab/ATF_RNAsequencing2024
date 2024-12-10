#The purpose of this RStudio ShinyApp is to take two input files, a count sheet and a metadata sheet, and generate a PCA plot, as well as offer to download the DESeq2 results from the ATF vs Shock

# Load Libraries
library(shiny)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(plotly)
library(matrixStats)

# Set the maximum upload file size to 100 MB
options(shiny.maxRequestSize = 100 * 1024^2)

# UI
ui <- fluidPage(
  titlePanel("DESeq2 Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("count_file", "Upload Count Data (CSV)", accept = c(".csv")),
      fileInput("metadata_file", "Upload Metadata (CSV)", accept = c(".csv")),
      actionButton("run_analysis", "Run Analysis"),
      downloadButton("download_results", "Download DESeq2 Results as CSV")
    ),
    mainPanel(
      plotlyOutput("pca_plot")
    )
  )
)

server <- function(input, output, session) {
  
  # Reactive values to store data 
  data <- reactiveValues(
    raw_counts = NULL, 
    metadata = NULL, 
    dds = NULL,
    deseq_results = NULL
  )
  
  observeEvent(input$run_analysis, {
    req(input$count_file, input$metadata_file)
    withProgress(message = "Running DESeq2 analysis...", value = 0, {
      incProgress(0.2, detail = "Loading data...")
      
      # Load metadata and count data
      metadata <- read.csv(input$metadata_file$datapath, header = TRUE)
      counts <- read.csv(input$count_file$datapath, header = TRUE)
      
      # Validate that the columns match between count data and metadata
      if (!all(metadata$name %in% colnames(counts))) {
        showModal(modalDialog(
          title = "Error",
          "Metadata samples do not match the count data columns.",
          easyClose = TRUE,
          footer = NULL
        ))
        return()
      }
      
      # Match count data columns with metadata 'name' field
      matched_counts <- counts[, intersect(colnames(counts), metadata$name)]
      rownames(matched_counts) <- counts$Geneid
      
      data$metadata <- metadata
      data$raw_counts <- matched_counts
      
      incProgress(0.3, detail = "Data loaded. Running DESeq2...")
      
      # Run DESeq2
      dds <- DESeqDataSetFromMatrix(countData = data$raw_counts, colData = metadata, design = ~ treatment)
      dds <- DESeq(dds)
      data$dds <- dds
      
      incProgress(0.5, detail = "DESeq2 analysis completed.")
      
      # Get DE results
      results_atf3_vs_shock <- results(dds, contrast = c("treatment", "ATF3", "Shock"))
      data$deseq_results <- results_atf3_vs_shock
      
      incProgress(1, detail = "Analysis done.")
    })
  })
  
  # PCA Plot
  output$pca_plot <- renderPlotly({
    req(data$deseq_results, data$raw_counts, data$metadata)
    
    counts_matrix <- as.matrix(data$raw_counts)
    
    # Filter out genes with zero or constant variance across samples
    gene_variances <- rowVars(counts_matrix)
    counts_nonzero_var <- counts_matrix[gene_variances > 0, , drop = FALSE]
    
    # Perform PCA on the transposed data
    pca_data <- prcomp(t(counts_nonzero_var), scale. = TRUE)
    
    # Create a data frame with the first three PCs and sample info
    pca_df <- data.frame(
      PC1 = pca_data$x[,1],
      PC2 = pca_data$x[,2],
      PC3 = pca_data$x[,3],
      sample = colnames(counts_matrix),  # Each column is a sample/replicate
      treatment = data$metadata$treatment[match(colnames(counts_matrix), data$metadata$name)],
      stringsAsFactors = FALSE
    )
    
    # Use the sample name as the label for each replicate
    plot_ly(pca_df, x = ~PC1, y = ~PC2, z = ~PC3,
            text = ~sample, hoverinfo = 'text',  # Hover text shows replicate name
            color = ~treatment, type = 'scatter3d', mode = 'markers') %>%
      layout(title = '3D PCA Plot',
             scene = list(
               xaxis = list(title = 'PC1'),
               yaxis = list(title = 'PC2'),
               zaxis = list(title = 'PC3')
             ))
  })
  
  # Download DESeq2 Results
  output$download_results <- downloadHandler(
    filename = function() {
      paste("DESeq2_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(data$deseq_results)
      write.csv(as.data.frame(data$deseq_results), file, row.names = TRUE)
    }
  )
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
