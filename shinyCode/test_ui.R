# Load packages 
library(shiny)
library(ggplot2)
library(TSstudio)
library(readr)
library(plotly)
library(dplyr)
library(tidyr)

# import modules 
source("sim_modules.R")

# Define UI 
ui <- fluidPage(title="Demo", 
                sidebarLayout(
                  sidebarPanel(
                    # user uploads csv file 
                    fileInput("file", "Choose CSV File",
                              accept = c("text/csv",
                                         "text/comma-separated-values,
                              .csv")),
                    
                    uiOutput('action'), 
                    uiOutput('size'), 
                    uiOutput('steps'),
                    uiOutput('ampContemp'),
                    uiOutput('ampLagged'),
                    uiOutput('ampMeasure'),
                    uiOutput('covContemptName'),
                    uiOutput('cogLaggedName'),
                    uiOutput('mask'),
                    uiOutput('reset'), 
                    uiOutput('variables')
                    
                  ),
                  
                  mainPanel(
                    
                    h1("Plotted Data"),
                    
                    plotOutput(outputId = "plot1",
                               dblclick = "plot1_dblclick",
                               brush = brushOpts(
                                 id = "plot1_brush",
                                 resetOnNew = TRUE),
                    ),
                    
                    plotOutput(outputId = "plot2"),
                             
                    downloadButton(
                      outputId = "download_plot",
                      label = "Download Plot Image"
                    ),
                    
                    downloadButton(
                      outputId = "unclipped",
                      label = "Download Unclipped Matrices"
                    ),
                    
                    downloadButton(
                      outputId = "clipped",
                      label = "Download Clipped Matrices"
                    )
                  )
                )
)

# Define server ----------------------------------------------------------------

server <- function(input, output, session) {
  
  filedata <- reactive({
    infile <- input$file
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath, header = TRUE)
  })
  
  output$action <- renderUI({
    if (is.null(filedata())) return()
    actionButton("analyze", "Analyze CSV file")
  })
  
  output$size <- renderUI({
    df <- filedata()
    if (is.null(df)) return(NULL)
    numericInput(inputId = 'size',
                 label = 'Number of variables',
                 value = 6)  # Default value
  })
  
  output$steps <- renderUI({
    df <- filedata()
    if (is.null(df)) return(NULL)
    numericInput(inputId = 'steps',
                 label = 'Number of time steps',
                 value = 100) 
  })
  
  output$ampContemp <- renderUI({
    df <- filedata()
    if (is.null(df)) return(NULL)
    numericInput(inputId = 'ampContemp', 
                 label = 'Amplitude of contemp matrix',
                 value = 0.1)  
  })
  
  output$ampLagged <- renderUI({
    df <- filedata()
    if (is.null(df)) return(NULL)
    numericInput(inputId = 'ampLagged', 
                 label = 'Amplitude of lagged matrix',
                 value = 0.1)
  })
  
  output$ampMeasure <- renderUI({
    df <- filedata()
    if (is.null(df)) return(NULL)
    numericInput(inputId = 'ampMeasure', 
                 label = 'Measure amplitude',
                 value = 1.0)
  })
  
  output$covContempName <- renderUI({
    df <- filedata()
    if (is.null(df)) return(NULL)
    selectInput(inputId = 'covContempName',
                label = 'Select contemp covariance type',
                choices = list('random' = 1, 'heterogeneous' = 2, 'compound' = 3), selected = 1)
  })
  
  output$covLaggedName <- renderUI({
    df <- filedata()
    if (is.null(df)) return(NULL)
    selectInput(inputId = 'covLaggedName',
                label = 'Select lagged covariance type',
                choices = list('random' = 1, 'heterogeneous' = 2, 'compound' = 3), selected = 1)
  })
  
  output$mask <- renderUI({
    df <- filedata()
    if (is.null(df)) return(NULL)
    selectInput(inputId = 'mask',
                label = 'Apply mask',
                choices = list('Yes' = 1, 'No' = 2), selected = 1)
  })
  
  output$reset <- renderUI({
    df <- filedata()
    if (is.null(df)) return(NULL)
    actionButton("reset", "Reset form")
    
  })
  
  observeEvent(input$reset, {
    session$reload()
  })
  
  # load in contemp and lagged matrices
  observeEvent(input$analyze, {
    data <- filedata()
    size <- input$size
    start <- rep(0, size)
    steps <- 100
    start_index <- size+1
    stop <- length(data)
    matContemp <- data[, start_index:stop]
    matLagged <- data[, 1:size]
    
    # covariances
    covContemp <- matrix(rnorm(size**2), nrow=size, ncol=size)
    covLagged <- matrix(rnorm(size**2), nrow=size, ncol=size)
    ampContemp = 0.1 # input$ampContemp
    ampLagged = 0.1 # input$ampLagged
    ampMeasure = 1.0 #input$ampMeasure
    measureCov <- ampMeasure * diag(size)
    covContempName <- input$covContempName
    covLaggedName <- input$covLaggedName
    
    # masks
    mask <- input$mask
    maskContemp =  make_mask(matContemp, contemp=T)
    maskLagged = make_mask(matLagged, contemp=F)
    
    # clipping
    clip_indices <- c(0,1)
    clip_mins <- c(0.5, 0.7)
    clip_maxs <- c(0.8,1.3)
    
    # clipping outliers
    clip_sigma = 2
    
    ts <- generate_timeseries(start, steps, ampContemp, matContemp, covContemp, ampLagged, matLagged, covLagged, measureCov, debug=F)
    
    ts_clipped <- clip_outliers(ts, clip_sigma, ampMeasure, debug=T, start, ampContemp, matContemp, covContemp, ampLagged, matLagged, covLagged, measureCov)
    
    output$variables <- renderUI({
      items = names(as.data.frame(ts))
      names(items) = items
      checkboxGroupInput(inputId = "variables",
                         label = "Choose variables you want shown on plot", 
                         choices = items, 
                         selected = items)
    })
    
    data <- data.frame(days=1:(steps-1),
                       gather(select(as.data.frame(ts), input$variables), variable_o, value_o),
                       gather(select(as.data.frame(ts_clipped), input$variables), variable_c, value_c)
    )
    
    plotInput1 <- function() {
      data <- data.frame(days=1:(steps-1),
                         gather(select(as.data.frame(ts), input$variables), variable_o, value_o),
                         gather(select(as.data.frame(ts_clipped), input$variables), variable_c, value_c)
      )
      p <- ggplot(data, aes(x=days)) +
        geom_line(aes(y = value_o, color = variable_o, linetype = "Original")) +
        geom_line(aes(y = value_c, color = variable_c, linetype = "Clipped")) +
        xlab("time [days]") +
        ylab("value") +
        scale_linetype_manual(values = c(Original = "solid", Clipped = "dashed")) +
        labs(color = "Timeseries")
    }
    
    plotInput2 <- function() {
      data <- data.frame(days=1:(steps-1),
                         gather(select(as.data.frame(ts), input$variables), variable_o, value_o),
                         gather(select(as.data.frame(ts_clipped), input$variables), variable_c, value_c)
      )
    p <- ggplot(data, aes(x=days)) +
      geom_line(aes(y = value_o, color = variable_o, linetype = "Original")) +
      geom_line(aes(y = value_c, color = variable_c, linetype = "Clipped")) +
      xlab("time [days]") +
      ylab("value") +
      scale_linetype_manual(values = c(Original = "solid", Clipped = "dashed")) +
      labs(color = "Timeseries") +
      coord_cartesian(xlim = ranges2$x, ylim = ranges2$y, expand = FALSE)
    }
    
    output$plot1 <- renderPlot({
      validate(need(!is.null(input$variables), 'Please tick a box to show a plot.'))
      print(plotInput1())
    })
    
    output$plot2 <- renderPlot({
      print(plotInput2())
    })
    
    # https://stackoverflow.com/questions/70391083/shiny-double-click-zoom-plot-interaction-not-working-properly
    
    ranges2 <- reactiveValues(x = NULL, y = NULL)

    observe({
      brush <- input$plot2_brush
      if (!is.null(brush)) {
        ranges2$x <- c(brush$xmin, brush$xmax)
        ranges2$y <- c(brush$ymin, brush$ymax)
        
      } else {
        ranges2$x <- NULL
        ranges2$y <- NULL
      }
    })
    
    
    output$download_plot <- downloadHandler(
      filename = "plot.png",
      content = function(file) {
        ggsave(file, plotInput1())
      }
    )
    
    output$unclipped <- downloadHandler(
      filename = "unclipped_matrices.csv",
      content = function(file) {
        write_csv(as.data.frame(ts), path = file)
      }
    )
    
    output$clipped <- downloadHandler(
      filename = "clipped_matrices.csv",
      content = function(file) {
        write_csv(as.data.frame(ts_clipped), path = file)
      }
    )
  })
  
}

# Create the Shiny app object --------------------------------------------------

shinyApp(ui = ui, server = server)
