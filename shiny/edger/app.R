### EDGER EXPLORER

libDir <- "/cluster/gjb_lab/mgierlinski/R_shiny/library/4.1"
if (dir.exists(libDir)) .libPaths(libDir)

library(scales)
library(shiny)
library(tidyverse)
library(DT)
library(fenr)
source("../shiny_func.R")

css <- "table{font-size: 11px; background-color: #EAF5FF}"

### Read data ###

data <- sh_read_edger_data("../data")
contrasts <- unique(data$edger$contrast) |> as.character()
ontologies <- names(data$fterms)

max_points <- 1000


#######################################################################

ui <- shinyUI(fluidPage(

  tags$style(css),

  titlePanel("Analysis of Ssp4 and Ssp6 effects on Pseudomonas fluorescens using Tn-seq"),

  fluidRow(
    column(12,
      fluidRow(
        column(4,
          selectInput("contrast", "Contrast:", choices = contrasts),
          radioButtons("plot_type", "Plot type:", choices = c("Volcano" = "vol", "MA" = "ma"), inline = TRUE),
          plotOutput("main_plot", height = "480px", width = "100%", brush = "plot_brush", hover = "plot_hover")
        ),
        column(3,
          plotOutput("gene_plot", height = "400px",width = "100%")
        ),
        column(5,
          p("Gene list"),
          div(style = 'height: 200px; overflow-y: scroll', tableOutput("gene_info")),
          br(),
          div(style = 'height: 400px; overflow-y: scroll', tableOutput("enrichment")),
        )
      ),
      fluidRow(
        DT::dataTableOutput("all_gene_table")
      )
    )
  )
)
)


########################################################################################


# Define server logic required to draw a histogram
server <- function(input, output, session) {

  # Prevents RStudio from crashing when Shiny window closed manually
  session$onSessionEnded(function() {
    stopApp()
  })

  get_de <- function() {
    ctr <- input$contrast
    data$edger |>
      filter(contrast == ctr)
  }

  get_xy_data <- function() {
    de <- get_de()
    if (input$plot_type == "vol") {
      xy_data <- de |>
        mutate(x = logFC, y = -log10(PValue))
    } else if (input$plot_type == "ma") {
      xy_data <- de  |>
        mutate(x = logCPM, y = logFC)
    }
    xy_data
  }

  select_gene <- function(max_hover = 1) {
    xy_data <- get_xy_data()
    sel <- NULL
    tab_idx <- as.numeric(input$all_gene_table_rows_selected)
    if (!is.null(input$plot_brush)) {
      brushed <- brushedPoints(xy_data, input$plot_brush)
      sel <- brushed$feature_id
    } else if (!is.null(input$plot_hover)) {
      near <- nearPoints(xy_data, input$plot_hover, threshold = 20, maxpoints = max_hover)
      sel <- near$feature_id
    } else if (length(tab_idx) > 0) {
      sel <- xy_data[tab_idx, ] |> pull(feature_id)
    }
    return(sel)
  }

  output$gene_info <- renderTable({
    xy_data <- get_xy_data()
    sel <- select_gene()
    df <- NULL
    if (!is.null(sel) && length(sel) >= 1 && length(sel) <= max_points) {
      df <- xy_data |>
        filter(feature_id %in% sel) |>
        arrange(feature_id) |>
        select(feature_id, description, FDR)
    } else if (length(sel) > max_points) {
      df <- data.frame(Error = paste0('only ',max_points,' points can be selected.'))
    }
    df
  })

  enrichment_table <- function(terms) {
    xy_data <- get_xy_data()
    sel <- NULL
    fe <- NULL
    if (!is.null(input$plot_brush)) {
      brushed <- brushedPoints(xy_data, input$plot_brush)
      sel <- brushed$feature_id
      n <- length(sel)
      if (n > 0 && n <= max_points) {
        fe <- fenr::functional_enrichment(data$all_genes, sel, terms)
        if(!is.null(fe))
          fe <- filter(fe, p_adjust < 0.05)
      } else if (n > 0) {
        fe <- data.frame(Error = paste0('only ',max_points,' points can be selected.'))
      }
    }
    fe
  }

  output$enrichment <- renderTable({
    map_dfr(ontologies, function(ont) {
      enrichment_table(data$fterms[[ont]])
    })
  })


  output$gene_plot <- renderPlot({
    sel <- select_gene()
    if (!is.null(sel) && length(sel) > 0 && length(sel) <= max_points) {
      sh_plot_genes(data$cnts, sel)
    }
  })

  output$main_plot <- renderPlot({
    xy_data <- get_xy_data()
    tab_idx <- as.numeric(input$all_gene_table_rows_selected)

    if (input$plot_type == "vol") {
      g <- sh_plot_volcano(xy_data)
    } else if (input$plot_type == "ma") {
      g <- sh_plot_ma(xy_data)
    }
    if (length(tab_idx) >= 1) {
      g <- g + geom_point(data = xy_data[tab_idx, ], colour = "red", size = 2)
    }
    g
  })

  output$all_gene_table <- DT::renderDataTable({
    d <- get_xy_data() |>
      select(feature_id, logFC, FDR, description) |>
      mutate(across(where(is.numeric), ~signif(.x, 3)))
    DT::datatable(d, class = 'cell-border strip hover', selection = "single", rownames = FALSE)
  })
}

# Run the application
shinyApp(ui = ui, server = server)


# For testing
# input <- list(method = "fi", plot_type = "vol", contrast = "strainTfe2", enrichment = "go", intensity_scale = "lin")
