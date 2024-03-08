# Load packages 
library(shiny)
library(ggplot2)
library(TSstudio)
library(readr)
library(plotly)
library(dplyr)

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
                    
                    plotOutput(outputId = "plot"),
                    
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
                 value = 1)  
  })
  
  output$ampLagged <- renderUI({
    df <- filedata()
    if (is.null(df)) return(NULL)
    numericInput(inputId = 'ampLagged', 
                 label = 'Amplitude of lagged matrix',
                 value = 1)
  })
  
  output$ampMeasure <- renderUI({
    df <- filedata()
    if (is.null(df)) return(NULL)
    numericInput(inputId = 'ampMeasure', 
                 label = 'Measure amplitude',
                 value = 1)
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
    steps <- input$steps
    start_index <- size+1
    stop <- ncol(data)
    matContemp <- data[, start_index:stop]
    matLagged <- data[, 1:size]
    
    # covariances
    covContemp <- matrix(rnorm(size**2), nrow=size, ncol=size)
    covLagged <- matrix(rnorm(size**2), nrow=size, ncol=size)
    ampContemp = input$ampContemp
    ampLagged = input$ampLagged
    ampMeasure = input$ampMeasure
    covMeasure <- ampMeasure * diag(size)
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
    
    ts = generate_timeseries(start, steps, ampContemp, matContemp, covContemp, ampLagged, matLagged, covLagged, covMeasure, debug=F)
    
    ts_clipped = clip_outliers(ts, clip_sigma, ampMeasure, debug=T, start, ampContemp, matContemp, covContemp, ampLagged, matLagged, covLagged, covMeasure)
    
    # current thought process: we have separated the columns of ts and clipped into parameters 
    # now need to work on making the check box appear AFTER file is uploaded and we know how many variables 
    # then allow checking boxing specific parameters and them showing on the plot
    # https://stackoverflow.com/questions/67595612/shiny-app-which-uploads-a-datafile-and-create-a-checkbox-and-textinput-and-dropd
    
    output$variables <- renderUI({
      df <- filedata()
      if (is.null(df)) return(NULL)
      items = names(df)
      names(items) = items
      checkboxGroupInput(inputId = "variables",
                         label = "Choose variables you want shown on plot", 
                         choices = items, 
                         selected = items)
    })
    
    plotInputs = numeric(length(input$variables))
    for (i in 1:length(input$variables)) {
      plotInputs[i] <- function() {
        df <- data.frame(
          days = 1:(steps-1),
          ts[, i],
          ts_clipped[, i]
        )
        p <- ggplot(df, aes(x = days)) +
          geom_line(aes(y = ts), linetype = "solid") +
          geom_line(aes(y = ts_clipped), linetype = "dashed") +
          xlab("time [days]") +
          ylab("value") +
          ggtitle("original vs. clipped timeseries") +
          labs(color = "Timeseries")
      }
      i = i + 1
    }

    
    output$plot <- renderPlot({
      for (input in plotInputs) {
        print(input())
      }
    })
    
    output$download_plot <- downloadHandler(
      filename = "plot.png",
      content = function(file) {
        ggsave(file, plotInput())
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
