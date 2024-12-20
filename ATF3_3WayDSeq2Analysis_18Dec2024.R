# Load Libraries
library(shiny)
library(DESeq2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

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
      downloadButton("download_atf0_shock_results", "Download ATF0vsShock DESeq2 Results as CSV"),
      downloadButton("download_atf0_aart6_results", "Download ATF0vsAart6 DESeq2 Results as CSV"),
      downloadButton("download_aart6_shock_results", "Download Aart6vsShock DESeq2 Results as CSV"),
      downloadButton("download_commonDEG_results", "Download Common Gene DESeq2 Results as CSV"),
      downloadButton("download_heatmap", "Download Heatmap as PDF"),
      
      # Text input for additional genes for the heatmap
      textInput("heatmap_genes", "Additional Genes for Heatmap (space-separated)", 
                placeholder = "e.g. GENE_X GENE_Y GENE_Z"),
      
      uiOutput("boxplot_gene_selector"),
      downloadButton("download_boxplots", "Download Boxplots as PDF")
    ),
    mainPanel(
      h4("Top 10 Rows of DESeq2 Results:"),
      tableOutput("deseq_table"),
      plotOutput("heatmap_plot", height = "1200px"),
      plotOutput("boxplot_plot")
    )
  )
)

server <- function(input, output, session) {

  #Set Data as Reactive
  data <- reactiveValues(
    raw_counts = NULL, 
    metadata = NULL, 
    deseq_results_atf0_vs_shock = NULL, 
    deseq_results_atf0_vs_aart6 = NULL, 
    deseq_results_aart6_vs_shock = NULL,
    common_genes_results = NULL, 
    dds = NULL
  )
  
  # Run DESeq2 analysis
  observeEvent(input$run_analysis, {
    req(input$count_file, input$metadata_file)  # Ensure both files are uploaded
    withProgress(message = "Running DESeq2 analysis...", value = 0, {
      incProgress(0.2, detail = "Loading data...")
      
      # Load metadata and count data
      metadata <- read.csv(input$metadata_file$datapath, header = TRUE)
      counts <- read.csv(input$count_file$datapath, header = TRUE)
      
      # Ensure that the columns match between count data and metadata
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
      
      # DESeq2 Analysis
      dds <- DESeqDataSetFromMatrix(countData = data$raw_counts, colData = metadata, design = ~ treatment)
      dds <- DESeq(dds)
      data$dds <- dds
      
      incProgress(0.5, detail = "DESeq2 analysis completed.")
      
      # Comparisons
      results_atf0_vs_shock <- results(dds, contrast = c("treatment", "ATF0", "Shock"))
      results_atf0_vs_aart6 <- results(dds, contrast = c("treatment", "ATF0", "Aart6"))
      results_aart6_vs_shock <- results(dds, contrast = c("treatment", "Aart6", "Shock"))
      
      data$deseq_results_atf0_vs_shock <- results_atf0_vs_shock
      data$deseq_results_atf0_vs_aart6 <- results_atf0_vs_aart6
      data$deseq_results_aart6_vs_shock <- results_aart6_vs_shock
      
      # Find common genes
      significant_atf0_vs_shock <- as.data.frame(na.omit(results_atf0_vs_shock))
      significant_atf0_vs_aart6 <- as.data.frame(na.omit(results_atf0_vs_aart6))
      
      significant_atf0_vs_shock <- significant_atf0_vs_shock[significant_atf0_vs_shock$padj < 0.05, ]
      significant_atf0_vs_aart6 <- significant_atf0_vs_aart6[significant_atf0_vs_aart6$padj < 0.05, ]
      
      common_genes <- intersect(rownames(significant_atf0_vs_shock), rownames(significant_atf0_vs_aart6))
      
      common_genes_results <- list(
        ATF0_vs_Shock = significant_atf0_vs_shock[common_genes, ],
        ATF0_vs_Aart6 = significant_atf0_vs_aart6[common_genes, ]
      )
      data$common_genes_results <- common_genes_results
      
      # Create rlog transformed data
      data$rld <- rlog(data$dds, blind = FALSE)
    })
  })
  
  # Reactive expression to parse UI input for additional heatmap genes
  user_heatmap_genes <- reactive({
    if (is.null(input$heatmap_genes) || input$heatmap_genes == "") {
      return(character(0))
    }
    genes <- unlist(strsplit(input$heatmap_genes, " "))
    genes <- genes[genes != ""]
    return(genes)
  })
  
  # Reactive expression for the heatmap
  heatmap_matrix_reactive <- reactive({
    req(data$raw_counts, data$metadata)
    
    # Log-transform the raw counts
    heatmap_data_log <- log2(data$raw_counts + 1)
    
    # Define selected genes (always present)
    selected_genes <- c("IFNG", "IFNGR1", "IFNGR2", "TNF", "HLA-A", "HLA-B", "HLA-C", "HLA-F", "HLA-J", 
                         "HLA-DPB2", "HLA-DQA2", "HLA-DRA", "IL4", "IL11", "SLCO3A1", "MAF", "AIF1L", "CRABP2", 
                         "IRF8", "IRF9", "TAP1", "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6", "REL",
                         "RELA", "RELB", "NFKB1", "NFKB2")
    
    extra_genes <- user_heatmap_genes()
    all_genes <- c(selected_genes, extra_genes)
    
    # Subset the data
    heatmap_data <- heatmap_data_log[rownames(heatmap_data_log) %in% all_genes, , drop = FALSE]
    # Reorder rows to match the order in all_genes
    heatmap_data <- heatmap_data[match(all_genes, rownames(heatmap_data)), , drop = FALSE]
    
    # Remove rows that aren't found (NAs)
    heatmap_data <- heatmap_data[!is.na(rownames(heatmap_data)), , drop = FALSE]
    
    if (nrow(heatmap_data) == 0) {
      # Return NULL if no genes found
      return(NULL)
    }
    
    # Reorder columns based on metadata
    metadata <- data$metadata
    heatmap_data <- heatmap_data[, metadata$name, drop = FALSE]
    
    # Average the data by treatment group
    grouped_data <- t(heatmap_data) %>%
      as.data.frame() %>%
      mutate(SampleGroup = metadata$treatment) %>%
      group_by(SampleGroup) %>%
      summarise(across(everything(), mean)) %>%
      column_to_rownames("SampleGroup")
    
    # Convert back to matrix
    heatmap_matrix <- t(as.matrix(grouped_data))
    return(heatmap_matrix)
  })
  
  # Heatmap Function
  draw_heatmap <- function(heatmap_matrix) {
    pheatmap(
      heatmap_matrix,
      scale = "row",
      show_rownames = TRUE,
      show_colnames = TRUE,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      cellwidth = 40,
      cellheight = 15,
      fontsize_row = 12,
      fontsize_col = 12,
      angle_col = 0,
      color = colorRampPalette(c("blue", "white", "red"))(100),
      labels_row = parse(text = paste0("italic('", rownames(heatmap_matrix), "')"))
    )
  }
  
  # Plot Heatmap
  output$heatmap_plot <- renderPlot({
    req(data$dds) # Ensure DESeq2 is done
    heatmap_matrix <- heatmap_matrix_reactive()
    req(!is.null(heatmap_matrix))  # Ensure the matrix is valid
    draw_heatmap(heatmap_matrix)
  })
  
  # Download heatmap as PDF
  output$download_heatmap <- downloadHandler(
    filename = function() {
      paste("Heatmap_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width = 10, height = 8)
      heatmap_matrix <- heatmap_matrix_reactive()
      req(!is.null(heatmap_matrix))
      draw_heatmap(heatmap_matrix)
      dev.off()
    }
  )
  
  
  # UI output for gene selection
  output$boxplot_gene_selector <- renderUI({
    req(data$dds)  # Ensure DESeq2 analysis is done
    
    selected_genes <- c("IFNG", "IFNGR1", "IFNGR2", "TNF", "HLA-A", "HLA-B", "HLA-C", "HLA-F", "HLA-J", 
                        "HLA-DPB2", "HLA-DQA2", "HLA-DRA", "IL4", "IL11", "SLCO3A1", "MAF", "AIF1L", "CRABP2", 
                        "IRF8", "IRF9", "TAP1", "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6", "REL",
                        "RELA", "RELB", "NFKB1", "NFKB2")
    selectInput("selected_genes_boxplot", "Select Genes for Boxplot:", 
                choices = selected_genes, multiple = TRUE, selected = selected_genes[1:25])
  })
  
  # Reactive Boxplot Data
  boxplot_data_reactive <- reactive({
    req(data$rld, input$selected_genes_boxplot)
    rld_mat <- assay(data$rld)
    metadata <- data$metadata
    
    sel_genes <- input$selected_genes_boxplot
    rld_sub <- rld_mat[rownames(rld_mat) %in% sel_genes, , drop = FALSE]
    
    if (nrow(rld_sub) == 0) {
      return(NULL)
    }
    
    df <- as.data.frame(t(rld_sub))
    df$Sample <- rownames(df)
    df <- df %>%
      tidyr::gather(key = "Gene", value = "Expression", -Sample) %>%
      left_join(metadata, by = c("Sample" = "name"))
    return(df)
  })
  
  # Plot Boxplot
  boxplot_ggplot_reactive <- reactive({
    df <- boxplot_data_reactive()
    req(df)
    
    p <- ggplot(df, aes(x = treatment, y = Expression)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      facet_wrap(~ Gene, scales = "free_y") +
      labs(title = "Gene Expression by Treatment", x = "Treatment", y = "rlog Expression") +
      theme_minimal() +
      theme(
        strip.text = element_text(size = 14),  # Gene labels (facet titles)
        axis.title = element_text(size = 14),  # Axis titles
        axis.text = element_text(size = 14),   # Axis text
        plot.title = element_text(size = 14)   # Plot title
      )
    return(p)
  })
  
  output$boxplot_plot <- renderPlot({
    req(data$dds) # Ensure analysis is done
    boxplot_ggplot_reactive()
  })
  
  output$download_boxplots <- downloadHandler(
    filename = function() {
      paste("Boxplots_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      p <- boxplot_ggplot_reactive()
      ggsave(file, plot = p, device = "pdf", width = 12, height = 8)
    }
  )
  
  
  # Display the top 10 rows of DESeq2 results
  output$deseq_table <- renderTable({
    req(data$deseq_results_atf0_vs_shock)
    head(data$deseq_results_atf0_vs_shock, 10)
  }, rownames = TRUE)
  
  # Download DESeq2 Results
  output$download_atf0_shock_results <- downloadHandler(
    filename = function() {
      paste("DESeq2results_ATF0vsShock", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(data$deseq_results_atf0_vs_shock)
      write.csv(as.data.frame(data$deseq_results_atf0_vs_shock), file, row.names = TRUE)
    }
  )
  
  output$download_atf0_aart6_results <- downloadHandler(
    filename = function() {
      paste("DESeq2results_ATF0vsAart6_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(data$deseq_results_atf0_vs_aart6)
      write.csv(as.data.frame(data$deseq_results_atf0_vs_aart6), file, row.names = TRUE)
    }
  )
  
  output$download_aart6_shock_results <- downloadHandler(
    filename = function() {
      paste("DESeq2results_Aart6vsShock_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(data$deseq_results_aart6_vs_shock)
      write.csv(as.data.frame(data$deseq_results_aart6_vs_shock), file, row.names = TRUE)
    }
  )
  
  output$download_commonDEG_results <- downloadHandler(
    filename = function() {
      paste("DESeq2results_ATF0vsShockandAart6_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(data$common_genes_results)
      # Convert the list to a data frame
      common_df <- do.call(rbind, lapply(data$common_genes_results, as.data.frame))
      write.csv(common_df, file, row.names = TRUE)
    }
  )
}

# Run the Shiny app
shinyApp(ui = ui, server = server)


