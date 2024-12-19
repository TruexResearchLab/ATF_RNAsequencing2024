# Load Libraries
library(shiny)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(plotly)

# Set the maximum upload file size to 100 MB
options(shiny.maxRequestSize = 100 * 1024^2)

# UI
ui <- fluidPage(
  titlePanel("DESeq2 Plot Generator"),
  sidebarLayout(
    sidebarPanel(
      fileInput("deseq_results_atf_shock", "Upload DESeq2 Results (ATF vs Shock) CSV", accept = c(".csv")),
      fileInput("deseq_results_atf_aart6", "Upload DESeq2 Results (ATF vs Aart6) CSV", accept = c(".csv")),
      
      textInput("user_genes", "Additional Genes to Label (space-separated)",
                placeholder = "e.g. GENE1 GENE2 GENE3"),
      
      checkboxInput("show_gridlines", "Show Gridlines", value = FALSE),
      
      # Checkbox for making Scatter1 interactive
      checkboxInput("interactive_scatter1", "Interactive Scatter Plot 1", value = FALSE),
      
      h4("Scatter Plot 1 Axis Limits"),
      numericInput("scatter1_xmin", "X min:", -50),
      numericInput("scatter1_xmax", "X max:", 90),
      numericInput("scatter1_ymin", "Y min:", -100),
      numericInput("scatter1_ymax", "Y max:", 80),
      
      h4("Scatter Plot 2 Axis Limits"),
      numericInput("scatter2_xmin", "X min:", -50),
      numericInput("scatter2_xmax", "X max:", 90),
      numericInput("scatter2_ymin", "Y min:", -100),
      numericInput("scatter2_ymax", "Y max:", 80),
      
      checkboxInput("auto_sig_label", "Automatic Significance Labeling", value = FALSE)
    ),
    mainPanel(
      # Conditional panels for Scatter1 (interactive vs static)
      conditionalPanel(
        condition = "!input.interactive_scatter1",
        plotOutput("scatter_plot_static"),
        downloadButton("download_scatter", "Download Scatter Plot as PDF")
      ),
      conditionalPanel(
        condition = "input.interactive_scatter1",
        plotlyOutput("scatter_plot_interactive"),
        downloadButton("download_scatter", "Download Scatter Plot as PDF")
      ),
      
      br(), br(),
      
      plotOutput("scatter_plot2"),
      downloadButton("download_scatter2", "Download Scatter2 Plot as PDF")
    )
  )
)

server <- function(input, output, session) {
  # Reactive values to store uploaded data
  data <- reactiveValues(
    deseq_results_atf_shock = NULL,
    deseq_results_atf_aart6 = NULL
  )
  
  # Load and process ATF vs Shock DESeq2 results
  observeEvent(input$deseq_results_atf_shock, {
    req(input$deseq_results_atf_shock)
    tryCatch({
      # Read the CSV file
      res_atf_shock <- read.csv(input$deseq_results_atf_shock$datapath, row.names = 1)
      
      # Validate columns
      required_cols <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
      if (!all(required_cols %in% colnames(res_atf_shock))) {
        showModal(modalDialog(
          title = "Error",
          "The uploaded file (ATF vs Shock) does not contain the required columns.",
          easyClose = TRUE,
          footer = NULL
        ))
        return()
      }
      
      data$deseq_results_atf_shock <- res_atf_shock
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("Error reading ATF vs Shock file:", e$message),
        easyClose = TRUE,
        footer = NULL
      ))
    })
  })
  
  # Load/process ATF vs Aart6 DESeq2 results
  observeEvent(input$deseq_results_atf_aart6, {
    req(input$deseq_results_atf_aart6)
    tryCatch({
      # Read the CSV file
      res_atf_aart6 <- read.csv(input$deseq_results_atf_aart6$datapath, row.names = 1)
      
      # Validate columns
      required_cols <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
      if (!all(required_cols %in% colnames(res_atf_aart6))) {
        showModal(modalDialog(
          title = "Error",
          "The uploaded file (ATF vs Aart6) does not contain the required columns.",
          easyClose = TRUE,
          footer = NULL
        ))
        return()
      }
      
      data$deseq_results_atf_aart6 <- res_atf_aart6
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("Error reading ATF vs Aart6 file:", e$message),
        easyClose = TRUE,
        footer = NULL
      ))
    })
  })
  
  # Reactive expression for defined genes
  user_genes <- reactive({
    genes <- unlist(strsplit(input$user_genes, " "))
    genes <- genes[genes != ""]
    return(genes)
  })
  
  # Reactive expression to create data frame for the scatter plot
  scatterData <- reactive({
    req(data$deseq_results_atf_shock, data$deseq_results_atf_aart6)
    
    df_shock <- as.data.frame(data$deseq_results_atf_shock)
    df_aart6 <- as.data.frame(data$deseq_results_atf_aart6)
    
    df_shock$Gene <- rownames(df_shock)
    df_aart6$Gene <- rownames(df_aart6)
    
    merged_df <- inner_join(df_shock %>% select(Gene, stat),
                            df_aart6 %>% select(Gene, stat),
                            by="Gene",
                            suffix=c("_ATF_Shock","_ATF_Aart6"))
    return(merged_df)
  })
  
  # Once both datasets are available, update axis limits
  observeEvent(scatterData(), {
    req(scatterData())
    df <- scatterData()
    # Compute min and max
    x_min <- floor(min(df$stat_ATF_Shock, na.rm=TRUE))
    x_max <- ceiling(max(df$stat_ATF_Shock, na.rm=TRUE))
    y_min <- floor(min(df$stat_ATF_Aart6, na.rm=TRUE))
    y_max <- ceiling(max(df$stat_ATF_Aart6, na.rm=TRUE))
    
    # Update axis inputs
    updateNumericInput(session, "scatter1_xmin", value = x_min)
    updateNumericInput(session, "scatter1_xmax", value = x_max)
    updateNumericInput(session, "scatter1_ymin", value = y_min)
    updateNumericInput(session, "scatter1_ymax", value = y_max)
    
    updateNumericInput(session, "scatter2_xmin", value = x_min)
    updateNumericInput(session, "scatter2_xmax", value = x_max)
    updateNumericInput(session, "scatter2_ymin", value = y_min)
    updateNumericInput(session, "scatter2_ymax", value = y_max)
  })
  
  #Make Gridlines
  grid_theme <- reactive({
    base_theme <- theme_bw() +
      theme(
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(0.25,"cm")
      )
    if (input$show_gridlines) {
      base_theme
    } else {
      base_theme + theme(panel.grid = element_blank())
    }
  })
  
  # First scatter plot 
  scatterPlot <- reactive({
    req(scatterData())
    df <- scatterData()
    
    # Determine colors
    df$color <- "grey"
    df$color[df$stat_ATF_Shock <= -10 & df$stat_ATF_Aart6 <= -10] <- "blue"
    df$color[df$stat_ATF_Shock >= 10 & df$stat_ATF_Aart6 >= 10] <- "red"
    
    p <- ggplot(df, aes(x = stat_ATF_Shock, y = stat_ATF_Aart6,
                        text = paste("Gene:", Gene))) +
      geom_point(aes(color = color), size = 1.5, alpha = 0.7) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      labs(title = "Comparison of T-values",
           x = "T (ATF vs Shock)",
           y = "T (ATF vs Aart6)") +
      grid_theme() +
      scale_color_identity() +
      coord_cartesian(xlim = c(input$scatter1_xmin, input$scatter1_xmax),
                      ylim = c(input$scatter1_ymin, input$scatter1_ymax)) +
      theme(legend.position = "none")
    
    p
  })
  
  # Render static plot if not interactive
  output$scatter_plot_static <- renderPlot({
    scatterPlot()
  })
  
  # Render interactive plot if interactive
  output$scatter_plot_interactive <- renderPlotly({
    ggplotly(scatterPlot(), tooltip = "text")
  })
  
  output$download_scatter <- downloadHandler(
    filename = function(){
      paste("Scatter_Plot_Tvalues_", Sys.Date(), ".pdf", sep="")
    },
    content = function(file){
      df <- scatterData()
      df$color <- "grey"
      df$color[df$stat_ATF_Shock <= -10 & df$stat_ATF_Aart6 <= -10] <- "blue"
      df$color[df$stat_ATF_Shock >= 10 & df$stat_ATF_Aart6 >= 10] <- "red"
      
      p <- ggplot(df, aes(x = stat_ATF_Shock, y = stat_ATF_Aart6)) +
        geom_point(aes(color = color), size = 1.5, alpha = 0.7) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        labs(title = "Comparison of T-values",
             x = "T (ATF vs Shock)",
             y = "T (ATF vs Aart6)") +
        grid_theme() +
        scale_color_identity() +
        coord_cartesian(xlim = c(input$scatter1_xmin, input$scatter1_xmax),
                        ylim = c(input$scatter1_ymin, input$scatter1_ymax))
      
      ggsave(file, plot=p, device="pdf", width=8, height=6)
    }
  )
  
  # Second scatter plot
  scatterPlot2 <- reactive({
    req(scatterData())
    df <- scatterData()
    
    df$color <- "grey"
    df$color[df$stat_ATF_Shock <= -10 & df$stat_ATF_Aart6 <= -10] <- "blue"
    df$color[df$stat_ATF_Shock >= 10 & df$stat_ATF_Aart6 >= 10] <- "red"
    
    user_genes_vec <- user_genes()
    df$UserLabel <- df$Gene %in% user_genes_vec
    
    # Separate data by color groups
    rb_df <- df[df$color %in% c("red","blue"), ]
    grey_df <- df[df$color == "grey", ]
    
    # Condition for automatic labeling
    auto_label <- input$auto_sig_label
    
    # Label sets:
    rb_to_label <- if (auto_label) {
      rb_df  # all red/blue
    } else {
      rb_df[rb_df$UserLabel, ] # only user-labeled red/blue
    }
    
    # Grey to label: only user-labeled
    grey_to_label <- grey_df[grey_df$UserLabel, ]
    
    # Points not labeled
    rb_non_label <- rb_df[!rb_df$Gene %in% rb_to_label$Gene, ]
    grey_non_label <- grey_df[!grey_df$Gene %in% grey_to_label$Gene, ]
    
    p <- ggplot() +
      # Grey non-labeled
      geom_point(data = grey_non_label,
                 aes(x = stat_ATF_Shock, y = stat_ATF_Aart6),
                 color = "grey", size = 0.25, alpha = 0.7) +
      
      # Red/Blue non-labeled
      geom_point(data = rb_non_label,
                 aes(x = stat_ATF_Shock, y = stat_ATF_Aart6, color = color),
                 size = 2, alpha = 0.7) +
      
      # Red/Blue labeled
      geom_point(data = rb_to_label,
                 aes(x = stat_ATF_Shock, y = stat_ATF_Aart6, color = color),
                 size = 2, alpha = 0.7) +
      geom_text_repel(data = rb_to_label,
                      aes(x = stat_ATF_Shock, y = stat_ATF_Aart6, label = Gene, color = color),
                      size = 6, fontface = "italic", max.overlaps = Inf) +
      
      # Grey labeled (black)
      geom_point(data = grey_to_label,
                 aes(x = stat_ATF_Shock, y = stat_ATF_Aart6),
                 color = "black", size = 2, alpha = 0.7) +
      geom_text_repel(data = grey_to_label,
                      aes(x = stat_ATF_Shock, y = stat_ATF_Aart6, label = Gene),
                      size = 6, fontface = "italic", color = "black", max.overlaps = Inf) +
      
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      labs(title = "Comparison of T-values",
           x = "T (ATF vs Shock)",
           y = "T (ATF vs Aart6)") +
      grid_theme() +
      scale_color_identity() +
      coord_cartesian(xlim = c(input$scatter2_xmin, input$scatter2_xmax),
                      ylim = c(input$scatter2_ymin, input$scatter2_ymax))
    
    return(p)
  })
  
  output$scatter_plot2 <- renderPlot({
    scatterPlot2()
  })
  
  output$download_scatter2 <- downloadHandler(
    filename = function(){
      paste("Scatter_Plot2_Tvalues_", Sys.Date(), ".pdf", sep="")
    },
    content = function(file){
      ggsave(file, plot=scatterPlot2(), device="pdf", width=8, height=6)
    }
  )
}

# Run the Shiny app
shinyApp(ui=ui, server=server)
