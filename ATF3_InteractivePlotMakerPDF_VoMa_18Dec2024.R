# Load Libraries
library(shiny)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(plotly)

# Set the maximum upload file size to 100 MB
options(shiny.maxRequestSize = 100 * 1024^2)

# Predefined label genes
selected_gene <- c("IFNG", "IFNGR1", "IFNGR2", "TNF", "HLA-A", "HLA-B", "HLA-C", "HLA-F", "HLA-J", 
                   "HLA-DPB2", "HLA-DQA2", "HLA-DRA", "IL4", "IL11", "SLCO3A1", "MAF", "AIF1L", "CRABP2", 
                   "IRF8", "IRF9", "TAP1", "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6", "REL",
                   "RELA", "RELB", "NFKB1", "NFKB2")

ui <- fluidPage(
  titlePanel("DESeq2 Plot Generator"),
  sidebarLayout(
    sidebarPanel(
      fileInput("deseq_results", "Upload DESeq2 Results (CSV)", accept = c(".csv")),
      textAreaInput("extra_genes", "Additional Genes to Label (separated by spaces):",
                    value = paste(selected_gene, collapse = " "), rows = 3, width = "100%"),
      checkboxInput("show_gridlines", "Show Gridlines", value = FALSE),
      checkboxInput("interactive_plots", "Interactive Plots", value = FALSE),
      h4("Volcano Plot Axis Limits"),
      numericInput("volcano_xmin", "X min:", -10),
      numericInput("volcano_xmax", "X max:", 10),
      numericInput("volcano_ymin", "Y min:", 0),
      numericInput("volcano_ymax", "Y max:", 400),
      h4("MA Plot Axis Limits"),
      numericInput("ma_xmin", "X min:", 0),
      numericInput("ma_xmax", "X max:", 5),
      numericInput("ma_ymin", "Y min:", -10),
      numericInput("ma_ymax", "Y max:", 25)
    ),
    mainPanel(
      conditionalPanel(
        condition = "!input.interactive_plots",
        plotOutput("volcano_plot"),
        downloadButton("download_volcano", "Download Volcano Plot as PDF"),
        br(), br(),
        plotOutput("ma_plot"),
        downloadButton("download_ma", "Download MA Plot as PDF")
      ),
      conditionalPanel(
        condition = "input.interactive_plots",
        plotlyOutput("volcano_plotly"),
        downloadButton("download_volcano", "Download Volcano Plot as PDF"),
        br(), br(),
        plotlyOutput("ma_plotly"),
        downloadButton("download_ma", "Download MA Plot as PDF")
      )
    )
  )
)

server <- function(input, output, session) {
  data <- reactiveValues(deseq_results = NULL)
  
  # Load and process DESeq2 results
  observeEvent(input$deseq_results, {
    req(input$deseq_results)  
    tryCatch({
      deseq_results <- read.csv(input$deseq_results$datapath, row.names = 1)
      
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
      
      data$deseq_results <- deseq_results

      # Volcano defaults
      volcano_data <- as.data.frame(deseq_results)
      volcano_data$Gene <- rownames(volcano_data)
      volcano_data$y_val <- -log10(volcano_data$padj)
      # Handle padj =0/NA
      zero_nan_rows <- which(is.na(volcano_data$padj) | volcano_data$padj == 0)
      if (length(zero_nan_rows) > 0) {
        jittered_values <- round(rnorm(length(zero_nan_rows), mean = 350, sd = 5), 0)
        volcano_data$y_val[zero_nan_rows] <- jittered_values
      }
      v_xmin <- floor(min(volcano_data$log2FoldChange, na.rm = TRUE))
      v_xmax <- ceiling(max(volcano_data$log2FoldChange, na.rm = TRUE))
      v_ymin <- floor(min(volcano_data$y_val, na.rm = TRUE))
      v_ymax <- ceiling(max(volcano_data$y_val, na.rm = TRUE))
      
      # Update Volcano axis numericInputs
      updateNumericInput(session, "volcano_xmin", value = v_xmin)
      updateNumericInput(session, "volcano_xmax", value = v_xmax)
      updateNumericInput(session, "volcano_ymin", value = ifelse(v_ymin < 0, 0, v_ymin))
      updateNumericInput(session, "volcano_ymax", value = v_ymax)
      
      # MA defaults
      results_df <- as.data.frame(deseq_results)
      results_df$baseMean_log <- log10(results_df$baseMean + 1)
      
      m_xmin <- floor(min(results_df$baseMean_log, na.rm = TRUE))
      m_xmax <- ceiling(max(results_df$baseMean_log, na.rm = TRUE))
      m_ymin <- floor(min(results_df$log2FoldChange, na.rm = TRUE))
      m_ymax <- ceiling(max(results_df$log2FoldChange, na.rm = TRUE))
      
      # Update MA axis numericInputs
      updateNumericInput(session, "ma_xmin", value = m_xmin)
      updateNumericInput(session, "ma_xmax", value = m_xmax)
      updateNumericInput(session, "ma_ymin", value = m_ymin)
      updateNumericInput(session, "ma_ymax", value = m_ymax)
      
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("An error occurred while reading the file:", e$message),
        easyClose = TRUE,
        footer = NULL
      ))
    })
  })
  
  volcanoPlotData <- reactive({
    req(data$deseq_results)
    volcano_data <- as.data.frame(data$deseq_results)
    volcano_data$Gene <- rownames(volcano_data)
    zero_nan_rows <- which(is.na(volcano_data$padj) | volcano_data$padj == 0)
    volcano_data$y_val <- -log10(volcano_data$padj)
    if (length(zero_nan_rows) > 0) {
      jittered_values <- round(rnorm(length(zero_nan_rows), mean = 350, sd = 5), 0)
      volcano_data$y_val[zero_nan_rows] <- jittered_values
    }
    volcano_data
  })
  
  volcanoPlot <- reactive({
    volcano_data <- volcanoPlotData()
    
    # Combine UI-defined additional genes with the pre-defined selected_gene group
    extra_genes_vec <- unlist(strsplit(input$extra_genes, "\\s+"))
    extra_genes_vec <- extra_genes_vec[extra_genes_vec != ""]
    label_genes <- unique(c(selected_gene, extra_genes_vec))
    
    volcano_data$Color <- ifelse(volcano_data$Gene %in% label_genes, "black",
                                 ifelse(volcano_data$log2FoldChange < -1, "blue",
                                        ifelse(volcano_data$log2FoldChange > 1, "red", "grey")))
    
    x_min <- input$volcano_xmin
    x_max <- input$volcano_xmax
    y_min <- input$volcano_ymin
    y_max <- input$volcano_ymax
    
    grid_theme <- if (input$show_gridlines) {
      theme(panel.grid.major = element_line(color = "grey80"),
            panel.grid.minor = element_line(color = "grey90"))
    } else {
      theme(panel.grid = element_blank())
    }
    
    x_breaks <- pretty(c(x_min, x_max), n = 5)
    y_breaks <- pretty(c(y_min, y_max), n = 5)
    
    #If interactive...
    volcano_data$hover_text <- paste("Gene:", volcano_data$Gene)
    
    labeled_data <- volcano_data[volcano_data$Gene %in% label_genes, ]
    other_data <- volcano_data[!volcano_data$Gene %in% label_genes, ]
    
    #Plot Volcano
    p <- ggplot() +
      geom_point(data = other_data %>% filter(!Color %in% c("red", "blue")), 
                 aes(x = log2FoldChange, y = y_val, color = Color, text=hover_text), 
                 alpha = 0.7, size = 0.25, shape = 16) +
      geom_point(data = other_data %>% filter(Color %in% c("red", "blue")), 
                 aes(x = log2FoldChange, y = y_val, color = Color, text=hover_text), 
                 alpha = 0.7, size = 2, shape = 16) +
      geom_point(data = labeled_data, aes(x = log2FoldChange, y = y_val, text=hover_text),
                 color = "black", size = 2, shape = 16) +
      geom_text_repel(data = labeled_data, aes(x = log2FoldChange, y = y_val, label = Gene),
                      size = 6, fontface = "italic", color = "black",
                      max.overlaps = 30, force = 2, nudge_x = 0.1, nudge_y = 0.1) +
      scale_color_manual(values = c("red" = "red", "blue" = "blue",
                                    "grey" = "grey", "black" = "black")) +
      labs(x = "Log2 Fold Change", y = "-Log10 Adjusted p-value", title = "Volcano Plot") +
      theme_bw() +
      theme(
        legend.position = "none",
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(0.25,"cm")
      ) +
      grid_theme +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
      scale_x_continuous(breaks = x_breaks, expand = c(0,0)) +
      scale_y_continuous(breaks = y_breaks, expand = c(0,0)) +
      coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max))
    
    p
  })
  
  output$volcano_plot <- renderPlot({
    volcanoPlot()
  })
  
  output$volcano_plotly <- renderPlotly({
    ggplotly(volcanoPlot(), tooltip = "text")
  })
  
  output$download_volcano <- downloadHandler(
    filename = function() {
      paste("Volcano_Plot_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = volcanoPlot(), device = "pdf", width = 10, height = 8)
    }
  )
  
  maPlotData <- reactive({
    req(data$deseq_results)
    results_df <- as.data.frame(data$deseq_results)
    results_df$Gene <- rownames(results_df)
    results_df$baseMean_log <- log10(results_df$baseMean + 1)
    results_df
  })
  
  maPlot <- reactive({
    results_df <- maPlotData()
    
    # Combine UI-defined additional genes with the pre-defined selected_gene group
    extra_genes_vec <- unlist(strsplit(input$extra_genes, "\\s+"))
    extra_genes_vec <- extra_genes_vec[extra_genes_vec != ""]
    label_genes <- unique(c(selected_gene, extra_genes_vec))
    
    results_df$Color <- ifelse(results_df$Gene %in% label_genes, "black",
                               ifelse(results_df$log2FoldChange < -1, "blue",
                                      ifelse(results_df$log2FoldChange > 1, "red", "grey")))
    
    x_min <- input$ma_xmin
    x_max <- input$ma_xmax
    y_min <- input$ma_ymin
    y_max <- input$ma_ymax
    
    grid_theme <- if (input$show_gridlines) {
      theme(panel.grid.major = element_line(color = "grey80"),
            panel.grid.minor = element_line(color = "grey90"))
    } else {
      theme(panel.grid = element_blank())
    }
    
    x_breaks <- pretty(c(x_min, x_max), n = 5)
    y_breaks <- pretty(c(y_min, y_max), n = 5)
    
    results_df$hover_text <- paste("Gene:", results_df$Gene)
    
    labeled_data <- results_df[results_df$Gene %in% label_genes, ]
    other_data <- results_df[!results_df$Gene %in% label_genes, ]
    
    #Plot MA
    p <- ggplot() +
      geom_point(data = other_data %>% filter(!is.na(log2FoldChange)), 
                 aes(x = baseMean_log, y = log2FoldChange, color = Color, text=hover_text), 
                 alpha = 0.7, size = 0.25, shape = 16) +
      geom_point(data = other_data %>% filter(Color %in% c("red", "blue")), 
                 aes(x = baseMean_log, y = log2FoldChange, color = Color, text=hover_text), 
                 alpha = 0.7, size = 2, shape = 16) +
      geom_point(data = labeled_data, aes(x = baseMean_log, y = log2FoldChange, text=hover_text), 
                 color = "black", size = 2, shape = 16) +
      geom_text_repel(data = labeled_data, aes(x = baseMean_log, y = log2FoldChange, label = Gene),
                      size = 6, fontface = "italic", color = "black",
                      max.overlaps = 30, force = 2, nudge_x = 0.1, nudge_y = 0.1) +
      scale_color_manual(values = c("red" = "red", "blue" = "blue", "grey" = "grey", "black" = "black")) +
      labs(x = "Log10 Base Mean", y = "Log2 Fold Change", title = "MA Plot") +
      theme_bw() +
      theme(
        legend.position = "none",
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(0.25,"cm")
      ) +
      grid_theme +
      geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "black") +
      scale_x_continuous(breaks = x_breaks, expand = c(0,0)) +
      scale_y_continuous(breaks = y_breaks, expand = c(0,0)) +
      coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max))
    
    p
  })
  
  output$ma_plot <- renderPlot({
    maPlot()
  })
  
  output$ma_plotly <- renderPlotly({
    ggplotly(maPlot(), tooltip = "text")
  })
  
  output$download_ma <- downloadHandler(
    filename = function() {
      paste("MA_Plot_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = maPlot(), device = "pdf", width = 10, height = 8)
    }
  )
}

shinyApp(ui = ui, server = server)
