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
      fileInput("deseq_results_atf_shock", "Upload DESeq2 Results (ATF vs Shock) CSV", accept = c(".csv")),
      fileInput("deseq_results_atf_aart6", "Upload DESeq2 Results (ATF vs Aart6) CSV", accept = c(".csv"))
    ),
    mainPanel(
      plotOutput("scatter_plot"),
      downloadButton("download_scatter", "Download Scatter Plot as PDF"),
      plotOutput("scatter_plot2"),
      downloadButton("download_scatter2", "Download Scatter2 Plot as PDF"),
      plotOutput("scatter_plot3"),
      downloadButton("download_scatter3", "Download Scatter3 Plot as PDF")
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
      # Read the uploaded CSV file
      res_atf_shock <- read.csv(input$deseq_results_atf_shock$datapath, row.names = 1)
      
      # Validate necessary columns in the input
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
  
  # load and process ATF vs Aart6 DESeq2 results
  observeEvent(input$deseq_results_atf_aart6, {
    req(input$deseq_results_atf_aart6)
    tryCatch({
      # Read the uploaded CSV file
      res_atf_aart6 <- read.csv(input$deseq_results_atf_aart6$datapath, row.names = 1)
      
      # Validate necessary columns
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
  
  
  
  # Reactive expression to create data frame for the scatter plot, requires two sets of results to be uploaded
  scatterData <- reactive({
    req(data$deseq_results_atf_shock, data$deseq_results_atf_aart6)
    
    df_shock <- as.data.frame(data$deseq_results_atf_shock)
    df_aart6 <- as.data.frame(data$deseq_results_atf_aart6)
    
    # Extract stat columns
    df_shock$Gene <- rownames(df_shock)
    df_aart6$Gene <- rownames(df_aart6)
    
    merged_df <- inner_join(df_shock %>% select(Gene, stat),
                            df_aart6 %>% select(Gene, stat),
                            by="Gene",
                            suffix=c("_ATF_Shock","_ATF_Aart6"))
    return(merged_df)
  })
  
  scatterPlot <- reactive({
    req(scatterData())
    df <- scatterData()
    
    p <- ggplot(df, aes(x = stat_ATF_Shock, y = stat_ATF_Aart6)) +
      geom_point(alpha = 0.5) +
      # Add axis lines at x=0 and y=0
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      
      labs(title = "Comparison of T-values",
           x = "T (ATF vs Shock)",
           y = "T (ATF vs Aart6)") +
      theme_minimal() +
      theme(
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14)
      ) +
      coord_cartesian(xlim = c(-50,90), ylim = c(-100,80))
    
    return(p)
  })
  
  # Render scatter plot
  output$scatter_plot <- renderPlot({
    scatterPlot()
  })
  
  # Download scatter plot as PDF
  output$download_scatter <- downloadHandler(
    filename = function(){
      paste("Scatter_Plot_Tvalues_",Sys.Date(),".pdf",sep="")
    },
    content = function(file){
      ggsave(file, plot=scatterPlot(), device="pdf", width=8, height=6)
    }
  )
  
  scatterPlot2 <- reactive({
    req(scatterData())
    df <- scatterData()
    
    # Assign colors based on conditions
    df$color <- "black"
    df$color[df$stat_ATF_Shock <= -10 & df$stat_ATF_Aart6 <= -10] <- "blue"
    df$color[df$stat_ATF_Shock >= 10 & df$stat_ATF_Aart6 >= 10] <- "red"
    
    # Separate subsets for labeling
    highlight_blue <- df[df$color == "blue", ]
    highlight_red <- df[df$color == "red", ]
    
    p <- ggplot(df, aes(x = stat_ATF_Shock, y = stat_ATF_Aart6)) +
      geom_point(aes(color = color), alpha = 0.5) +
      # Add axis lines at x=0 and y=0
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      # Label the highlighted genes
      geom_text_repel(data = highlight_blue, aes(label = Gene),
                      color = "blue", size = 6, fontface = "italic",
                      max.overlaps = Inf) +
      geom_text_repel(data = highlight_red, aes(label = Gene),
                      color = "red", size = 6, fontface = "italic",
                      max.overlaps = Inf) +
      labs(title = "Comparison of T-values",
           x = "T (ATF vs Shock)",
           y = "T (ATF vs Aart6)") +
      theme_minimal() +
      theme(
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14)
      ) +
      # Add extra space around the data on both axes
      coord_cartesian(xlim = c(-50,90), ylim = c(-100,80)) + 
      # Interpret color values directly
      scale_color_identity()
    
    return(p)
  })
  
  
  # Render scatter plot
  output$scatter_plot2 <- renderPlot({
    scatterPlot2()
  })
  
  # Download scatter plot as PDF
  output$download_scatter2 <- downloadHandler(
    filename = function(){
      paste("Scatter_Plot2_Tvalues_",Sys.Date(),".pdf",sep="")
    },
    content = function(file){
      ggsave(file, plot=scatterPlot2(), device="pdf", width=8, height=6)
    }
  )
  
  scatterPlot3 <- reactive({
    req(scatterData())
    df <- scatterData()
    
    # Assign colors based on conditions
    df$color <- "black"
    df$color[df$stat_ATF_Shock <= -10 & df$stat_ATF_Aart6 <= -10] <- "blue"
    df$color[df$stat_ATF_Shock >= 10 & df$stat_ATF_Aart6 >= 10] <- "red"
    
    # Separate subsets for labeling
    highlight_blue <- df[df$color == "blue", ]
    highlight_red <- df[df$color == "red", ]
    
    p <- ggplot(df, aes(x = stat_ATF_Shock, y = stat_ATF_Aart6)) +
      geom_point(aes(color = color), alpha = 0.5) +
      # Add axis lines at x=0 and y=0
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      # Label the highlighted genes
      geom_text_repel(data = highlight_blue, aes(label = Gene),
                      color = "blue", size = 6, fontface = "italic",
                      max.overlaps = Inf) +
      geom_text_repel(data = highlight_red, aes(label = Gene),
                      color = "red", size = 6, fontface = "italic",
                      max.overlaps = Inf) +
      labs(title = "Comparison of T-values",
           x = "T (ATF vs Shock)",
           y = "T (ATF vs Aart6)") +
      theme_minimal() +
      theme(
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14)
      ) +
      # Add extra space around the data on both axes
      scale_x_continuous(expand = expansion(mult = 0.1)) +
      scale_y_continuous(expand = expansion(mult = 0.1)) +
      # Interpret color values directly
      scale_color_identity()
    
    return(p)
  })
  
  
  # Render scatter plot
  output$scatter_plot3 <- renderPlot({
    scatterPlot3()
  })
  
  # Download scatter plot as PDF
  output$download_scatter3 <- downloadHandler(
    filename = function(){
      paste("Scatter_Plot3_Tvalues_",Sys.Date(),".pdf",sep="")
    },
    content = function(file){
      ggsave(file, plot=scatterPlot3(), device="pdf", width=8, height=6)
    }
  )
  
}

# Run the Shiny app
shinyApp(ui=ui, server=server)