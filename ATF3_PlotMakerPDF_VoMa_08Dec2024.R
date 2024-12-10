# Load Libraries
library(shiny)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Set the maximum upload file size to 100 MB
options(shiny.maxRequestSize = 100 * 1024^2)

# UI
ui <- fluidPage(
  titlePanel("DESeq2 Plot Generator"),
  sidebarLayout(
    sidebarPanel(
      fileInput("deseq_results", "Upload DESeq2 Results (CSV)", accept = c(".csv"))
    ),
    mainPanel(
      plotOutput("volcano_plot"),
      downloadButton("download_volcano", "Download Volcano Plot as PDF"),
      br(), br(),
      plotOutput("ma_plot"),
      downloadButton("download_ma", "Download MA Plot as PDF")
    )
  )
)

server <- function(input, output, session) {
  # Reactive values to store uploaded data
  data <- reactiveValues(
    deseq_results = NULL
  )
  
  # Load and process DESeq2 results
  observeEvent(input$deseq_results, {
    req(input$deseq_results)  # Ensure a file is uploaded
    tryCatch({
      # Read the uploaded CSV file
      deseq_results <- read.csv(input$deseq_results$datapath, row.names = 1)
      
      # Validate necessary columns in the input
      required_cols <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
      if (!all(required_cols %in% colnames(deseq_results))) {
        showModal(modalDialog(
          title = "Error",
          "The uploaded file does not contain the required columns for DESeq2 results.",
          easyClose = TRUE,
          footer = NULL
        ))
        return()
      }
      
      # Store the DESeq2 results in reactive values
      data$deseq_results <- deseq_results
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("An error occurred while reading the file:", e$message),
        easyClose = TRUE,
        footer = NULL
      ))
    })
  })
  
  volcanoPlot <- reactive({
    req(data$deseq_results)
    
    # Extract DESeq2 results
    volcano_data <- as.data.frame(data$deseq_results)
    volcano_data$Gene <- rownames(volcano_data)  # Add Gene names as a column
    
    # Identify rows where padj =0 or NaN
    zero_nan_rows <- which(is.na(volcano_data$padj) | volcano_data$padj == 0)
    
    # For non-zero and non-NaN padj, calculate y-axis as -log10(padj)
    volcano_data$y_val <- -log10(volcano_data$padj)
    
    # For zero/NaN padj values, set y-values around 350 with a sd of 5
    if (length(zero_nan_rows) > 0) {
      # Jitter values
      jittered_values <- round(rnorm(length(zero_nan_rows), mean = 350, sd = 5), 0)
      volcano_data$y_val[zero_nan_rows] <- jittered_values
    }
    
    # Define the genes to label
    label_genes <- c("IFNG", "HLA-A", "HLA-B", "HLA-C", "HLA-F",
                     "IL11", "HBA1", "AIF1L", "SLCO3A1", "MAF", "CRABP2")
    
    # Assign colors to points
    volcano_data$Color <- ifelse(volcano_data$Gene %in% label_genes, "black",
                                 ifelse(volcano_data$log2FoldChange < -1, "blue",
                                        ifelse(volcano_data$log2FoldChange > 1, "red", "grey")))
    
    # Separate data for labeled and non-labeled genes
    labeled_data <- volcano_data[volcano_data$Gene %in% label_genes, ]
    other_data <- volcano_data[!volcano_data$Gene %in% label_genes, ]
    
    # Create the plot
    p <- ggplot() +
      # Plot non-labeled, non-red/blue points first
      geom_point(data = other_data %>% filter(!Color %in% c("red", "blue")), 
                 aes(x = log2FoldChange, y = y_val, color = Color), 
                 alpha = 0.7, size = 0.25) +
      # Plot red and blue points next
      geom_point(data = other_data %>% filter(Color %in% c("red", "blue")), 
                 aes(x = log2FoldChange, y = y_val, color = Color), 
                 alpha = 0.7, size = 2) +
      # Plot labeled points last
      geom_point(data = labeled_data, aes(x = log2FoldChange, y = y_val), 
                 color = "black", size = 2) +
      # Add labels for specific genes
      geom_text_repel(data = labeled_data, aes(x = log2FoldChange, y = y_val, label = Gene),
                      size = 6, fontface = "italic", color = "black",
                      max.overlaps = 30, force = 2, nudge_x = 0.1, nudge_y = 0.1) +
      scale_color_manual(values = c("red" = "red", "blue" = "blue",
                                    "grey" = "grey", "black" = "black")) +
      labs(x = "Log2 Fold Change", y = "-Log10 Adjusted p-value", title = "Volcano Plot") +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17)
      ) +
      # Add vertical lines at log2FoldChange thresholds
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
      coord_cartesian(xlim = c(-10, 10), ylim = c(0, 400))
    
    return(p)
  })
  
  
  
  
  # Render Volcano Plot in UI
  output$volcano_plot <- renderPlot({
    volcanoPlot()
  })
  
  # Download handler for Volcano Plot
  output$download_volcano <- downloadHandler(
    filename = function() {
      paste("Volcano_Plot_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = volcanoPlot(), device = "pdf", width = 10, height = 8)
    }
  )
  
  # Reactive expression for MA Plot
  maPlot <- reactive({
    req(data$deseq_results)  # Ensure DESeq2 results are available
    
    # Extract results and add baseMean for x-axis
    results_df <- as.data.frame(data$deseq_results)
    results_df$baseMean_log <- log10(results_df$baseMean + 1)
    results_df$Gene <- rownames(results_df)  # Add Gene names as a column
    
    # Define the genes to label
    label_genes <- c("IFNG", "HLA-A", "HLA-B", "HLA-C", "HLA-F",
                     "IL11", "HBA1", "AIF1L", "SLCO3A1", "MAF", "CRABP2")
    
    # Assign base colors to points
    results_df$Color <- ifelse(results_df$Gene %in% label_genes, "black",
                               ifelse(results_df$log2FoldChange < -1, "blue",
                                      ifelse(results_df$log2FoldChange > 1, "red", "grey")))
    
    # Separate data for labeled and non-labeled genes
    labeled_data <- results_df[results_df$Gene %in% label_genes, ]
    other_data <- results_df[!results_df$Gene %in% label_genes, ]
    
    # Create the plot
    p <- ggplot() +
      # Plot non-labeled points first
      geom_point(data = other_data %>% filter(!is.na(log2FoldChange)), 
                 aes(x = baseMean_log, y = log2FoldChange, color = Color), 
                 alpha = 0.7, size = 0.25) +
      # Plot red and blue points with size 1
      geom_point(data = other_data %>% filter(Color %in% c("red", "blue")), 
                 aes(x = baseMean_log, y = log2FoldChange, color = Color), 
                 alpha = 0.7, size = 2) +
      # Plot labeled points with larger size
      geom_point(data = labeled_data, aes(x = baseMean_log, y = log2FoldChange), 
                 color = "black", size = 2) +
      # Add labels for specific genes
      geom_text_repel(data = labeled_data, aes(x = baseMean_log, y = log2FoldChange, label = Gene),
                      size = 6, fontface = "italic", color = "black",
                      max.overlaps = 30, force = 2, nudge_x = 0.1, nudge_y = 0.1) +
      # Define color scale for non-labeled points
      scale_color_manual(values = c("red" = "red", "blue" = "blue", "grey" = "grey", "black" = "black")) +
      labs(x = "Log10 Base Mean", y = "Log2 Fold Change", title = "MA Plot") +
      theme_minimal() +
      theme(legend.position = "none",
            axis.title = element_text(size = 17),
            axis.text = element_text(size = 17)
      ) +
      # Add horizontal lines at thresholds
      geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "black") +
      # Instead, use coord_cartesian for visible range:
      coord_cartesian(ylim = c(-10, 25))
    
    return(p)
  })
  
  # Render MA Plot in UI
  output$ma_plot <- renderPlot({
    maPlot()
  })
  
  # Download handler for MA Plot
  output$download_ma <- downloadHandler(
    filename = function() {
      paste("MA_Plot_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = maPlot(), device = "pdf", width = 10, height = 8)
    }
  )
}

# Run the Shiny app
shinyApp(ui = ui, server = server)