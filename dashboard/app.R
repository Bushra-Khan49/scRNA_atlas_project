library(shiny)
library(bslib)
library(Seurat)
library(plotly)
library(DT)

# Load data
pbmc <- readRDS("data/03_pbmc_slim.rds")

# Curated, vibrant color palette ensuring no dull greys
custom_colors <- c(
  "T_cells" = "#E63946",
  "B_cell" = "#457B9D",
  "Monocyte" = "#2A9D8F",
  "NK_cell" = "#F4A261",
  "Platelets" = "#E9C46A",
  "CMP" = "#8338EC",
  "Pre-B_cell_CD34-" = "#3A86FF",
  "Pro-B_cell_CD34+" = "#FB5607"
)

premium_theme <- bs_theme(
  version = 5,
  bootswatch = "darkly", 
  primary = "#3A86FF",
  secondary = "#8338EC",
  success = "#2A9D8F",
  base_font = font_google("Outfit"),
  heading_font = font_google("Outfit")
)

# Custom Mac Toggle UI Component
mac_toggle_ui <- div(
  class = "mac-switch-container",
  span("Light", class = "mac-switch-label"),
  tags$label(class = "mac-switch",
        tags$input(type = "checkbox", id = "mac_toggle", checked = "checked"),
        span(class = "mac-slider")
  ),
  span("Dark", class = "mac-switch-label")
)

# Custom Metric Card HTML Generator
custom_card <- function(title, value, info, icon_class, bg_class) {
  div(class = paste("custom-metric-card", bg_class),
      div(class = "metric-left",
          div(class = "metric-icon", tags$i(class = icon_class)),
          div(class = "metric-value", value)
      ),
      div(class = "metric-right",
          div(class = "metric-title", title),
          div(class = "metric-info", info)
      )
  )
}

ui <- page_navbar(
  theme = premium_theme,
  title = "Single-Cell Transcriptomics Explorer",
  fillable = FALSE, 
  
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
    tags$link(rel = "stylesheet", href = "https://cdn.jsdelivr.net/npm/bootstrap-icons@1.10.5/font/bootstrap-icons.css"),
    tags$script(HTML("
      $(document).on('shiny:connected', function() {
        // Broadcast the initial dark mode state to the server on load
        Shiny.setInputValue('dark_mode', $('#mac_toggle').is(':checked') ? 'dark' : 'light');
      });
      $(document).ready(function() {
        $(document).on('change', '#mac_toggle', function() {
           var theme = $(this).is(':checked') ? 'dark' : 'light';
           document.documentElement.setAttribute('data-bs-theme', theme);
           Shiny.setInputValue('dark_mode', theme);
        });
      });
    "))
  ),
  
  nav_spacer(),
  nav_item(mac_toggle_ui),
  
  nav_panel("Dashboard", 
    layout_sidebar(
      sidebar = sidebar(
        width = 380,
        h4("Controls"),
        
        selectInput("group_var", "Color Map By:", 
                    choices = c("Biological Cell Type" = "cell_type", 
                                "Mathematical Clusters" = "seurat_clusters")),
        
        selectizeInput("gene", "Examine Specific Gene:", 
                       choices = NULL, 
                       selected = "CD3D"),
                       
        selectInput("chart_type", "Expression Plot Type:",
                    choices = c("Violin Plot", "Box Plot", "Strip Plot"),
                    selected = "Violin Plot"),
        
        hr(),
        
        div(
          class = "sidebar-explanation",
          h4("Understanding the Data"),
          
          h5("UMAP (Dimensionality Reduction)"),
          p("Blood contains many different cell types. The UMAP algorithm takes 2,000+ genes and compresses them into 2D space. Cells that are close together are biologically similar. This allows us to instantly identify clusters of distinct immune cells."),
          
          h5("Expression Plot"),
          p("Shows the statistical distribution of the specific gene you selected. For example, if you select 'CD3D', you will see it is highly expressed almost exclusively in T-cells, proving it is a T-cell biomarker."),
          
          h5("Biomarker Data"),
          p("The table at the bottom shows the raw differential expression statistics (p-values and log fold changes). It mathematically proves which genes define each cluster.")
        )
      ),
      
      layout_columns(
        fill = FALSE,
        custom_card(
          title = "Total Cells Analyzed",
          value = ncol(pbmc),
          info = "Quality-controlled single cells passing mitochondrial thresholds.",
          icon_class = "bi-diagram-3",
          bg_class = "card-primary"
        ),
        custom_card(
          title = "Identified Cell Types",
          value = length(unique(pbmc$cell_type)),
          info = "Distinct biological identities discovered via reference annotation.",
          icon_class = "bi-tags",
          bg_class = "card-secondary"
        ),
        custom_card(
          title = "Total Genes Tracked",
          value = nrow(pbmc),
          info = "Expressed transcripts measured across all cells.",
          icon_class = "bi-activity",
          bg_class = "card-success"
        )
      ),
      
      layout_columns(
        col_widths = c(6, 6),
        
        card(
          full_screen = TRUE,
          card_header("Cellular Map (UMAP Projection)"),
          card_body(class="plotly-container", plotlyOutput("umap_plot", height = "550px"))
        ),
        
        card(
          full_screen = TRUE,
          card_header("Gene Expression Distribution"),
          card_body(class="plotly-container", plotlyOutput("expr_plot", height = "550px"))
        )
      ),
      
      card(
        full_screen = TRUE,
        card_header("Top Biomarker Genes (Differential Expression Statistics)"),
        card_body(DTOutput("marker_table"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  updateSelectizeInput(session, "gene", choices = sort(rownames(pbmc)), server = TRUE, selected = "CD3D")
  
  output$umap_plot <- renderPlotly({
    # CRITICAL: DO NOT use req(input$dark_mode) here, or the plot will be blank on initial load!
    
    df <- FetchData(pbmc, vars = c("umap_1", "umap_2", input$group_var))
    color_mapping <- if(input$group_var == "cell_type") custom_colors else "Paired"
    
    is_dark <- isTRUE(input$dark_mode == "dark") || is.null(input$dark_mode)
    text_color <- if(is_dark) "#ffffff" else "#222222"
    grid_color <- if(is_dark) "#444444" else "#dddddd"
    
    p <- plot_ly(df, x = ~umap_1, y = ~umap_2, 
                 color = as.formula(paste0("~", input$group_var)),
                 colors = color_mapping,
                 type = "scatter", mode = "markers", 
                 marker = list(size = 4, opacity = 0.9, line = list(width = 0)),
                 hoverinfo = "text",
                 text = ~paste("Identity:", get(input$group_var))) %>%
      layout(
        plot_bgcolor = "rgba(0,0,0,0)",
        paper_bgcolor = "rgba(0,0,0,0)",
        font = list(color = text_color),
        xaxis = list(
            title = "UMAP_1 (Component 1)", 
            zeroline = FALSE, showgrid = FALSE, 
            color = text_color, 
            tickfont = list(color = text_color, size = 12),
            showline = TRUE, linecolor = text_color, linewidth = 2,
            ticks = "outside", tickcolor = text_color, tickwidth = 2
        ),
        yaxis = list(
            title = "UMAP_2 (Component 2)", 
            zeroline = FALSE, showgrid = FALSE, 
            color = text_color, 
            tickfont = list(color = text_color, size = 12),
            showline = TRUE, linecolor = text_color, linewidth = 2,
            ticks = "outside", tickcolor = text_color, tickwidth = 2
        ),
        legend = list(font = list(color = text_color)),
        annotations = list(
            list(
              x = 1, y = 0, xref = "paper", yref = "paper",
              ax = -20, ay = 0, showarrow = TRUE, arrowhead = 2, arrowsize = 1.2, arrowwidth=2,
              arrowcolor = text_color, text = ""
            ),
            list(
              x = 0, y = 1, xref = "paper", yref = "paper",
              ax = 0, ay = 20, showarrow = TRUE, arrowhead = 2, arrowsize = 1.2, arrowwidth=2,
              arrowcolor = text_color, text = ""
            )
        )
      )
    p
  })

  output$expr_plot <- renderPlotly({
    req(input$gene) 
    req(input$gene %in% rownames(pbmc))
    
    df <- FetchData(pbmc, vars = c(input$gene, "cell_type"))
    colnames(df) <- c("expression", "cell_type")
    
    is_dark <- isTRUE(input$dark_mode == "dark") || is.null(input$dark_mode)
    text_color <- if(is_dark) "#ffffff" else "#222222"
    grid_color <- if(is_dark) "#444444" else "#cccccc"
    
    if (input$chart_type == "Violin Plot") {
      p <- plot_ly(df, x = ~cell_type, y = ~expression, 
                   color = ~cell_type, colors = custom_colors,
                   type = "violin", box = list(visible = TRUE),
                   meanline = list(visible = TRUE),
                   points = FALSE)
    } else if (input$chart_type == "Box Plot") {
      p <- plot_ly(df, x = ~cell_type, y = ~expression, 
                   color = ~cell_type, colors = custom_colors,
                   type = "box", boxpoints = FALSE)
    } else {
      p <- plot_ly(df, x = ~cell_type, y = ~expression, 
                   color = ~cell_type, colors = custom_colors,
                   type = "scatter", mode = "markers",
                   marker = list(size = 3, opacity = 0.6, line = list(width = 0)),
                   x0 = 0)
    }
    
    p %>% layout(
        plot_bgcolor = "rgba(0,0,0,0)",
        paper_bgcolor = "rgba(0,0,0,0)",
        font = list(color = text_color),
        xaxis = list(
            title = "Biological Cell Type", tickangle = 45, 
            color = text_color, tickfont = list(color = text_color, size = 12),
            showline = TRUE, linecolor = text_color, linewidth = 2,
            ticks = "outside", tickcolor = text_color, tickwidth = 2
        ),
        yaxis = list(
            title = paste(input$gene, "Log-Normalized Expression"), 
            color = text_color, tickfont = list(color = text_color, size = 12), 
            gridcolor = grid_color, gridwidth = 1,
            showline = TRUE, linecolor = text_color, linewidth = 2,
            ticks = "outside", tickcolor = text_color, tickwidth = 2
        ),
        showlegend = FALSE,
        margin = list(b = 100),
        annotations = list(
            list(
              x = 1, y = 0, xref = "paper", yref = "paper",
              ax = -20, ay = 0, showarrow = TRUE, arrowhead = 2, arrowsize = 1.2, arrowwidth=2,
              arrowcolor = text_color, text = ""
            ),
            list(
              x = 0, y = 1, xref = "paper", yref = "paper",
              ax = 0, ay = 20, showarrow = TRUE, arrowhead = 2, arrowsize = 1.2, arrowwidth=2,
              arrowcolor = text_color, text = ""
            )
        )
      )
  })

  output$marker_table <- renderDT({
    table_class <- 'cell-border stripe hover dark-table'
    
    if(file.exists("data/03_marker_genes.csv")) {
      data <- read.csv("data/03_marker_genes.csv")
      data$p_val <- signif(data$p_val, 3)
      data$p_val_adj <- signif(data$p_val_adj, 3)
      data$avg_log2FC <- round(data$avg_log2FC, 2)
      
      datatable(data, 
                options = list(pageLength = 5, scrollX = TRUE, dom = 'Bfrtip'),
                rownames = FALSE,
                class = table_class)
    }
  })
}

shinyApp(ui, server)
