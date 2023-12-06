# Load packages 
library(shiny)
library(ggplot2)
library(TSstudio)
library(readr)

# import modules 
source("~/ShinyApp/CIFAR/sim_modules.R")

# Define UI 
ui <- fluidPage(
  titlePanel("Demo"),
  
  sidebarLayout(
    sidebarPanel(
      # user uploads csv file 
      fileInput("file1", "Choose CSV File",
                accept = c("text/csv",
                           "text/comma-separated-values,
                       .csv")),
      
      numericInput(inputId = 'size',
                   label = 'Number of variables',
                   value = 6),  # Default value
      
      numericInput(inputId = 'steps',
                   label = 'Number of time steps',
                   value = 100),  
      
      numericInput(inputId = 'ampContemp', 
                   label = 'Amplitude of contemp matrix',
                   value = 1),  
      
      numericInput(inputId = 'ampLagged', 
                   label = 'Amplitude of lagged matrix',
                   value = 1),  
      
      numericInput(inputId = 'ampMeasure', 
                   label = 'Measure amplitude',
                   value = 1),
      
      selectInput(inputId = 'covContempName',
                  label = 'Select contemp covariance type',
                  choices = list('random' = 1, 'heterogeneous' = 2, 'compound' = 3), selected = 1),
      
      selectInput(inputId = 'covLaggedName',
                  label = 'Select lagged covariance type',
                  choices = list('random' = 1, 'heterogeneous' = 2, 'compound' = 3), selected = 1),
      
      selectInput(inputId = 'mask',
                  label = 'Apply mask',
                  choices = list('Yes' = 1, 'No' = 2), selected = 1),
      
      actionButton(inputId = "analyze", label = "Analyze CSV file")
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
  

  data_input <- reactive({
    infile <- input$file1
    if (is.null(infile)) {
      return(NULL)
    }
    read.csv(infile$datapath, header = TRUE)
  })
  
  # load in contemp and lagged matrices
  observeEvent(input$analyze, {
    data <- data_input()
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


  plotInput <- function() {
    df <- data.frame(
      days = 1:(steps-1),
      original = ts[, 1],  # only param 1
      clipped = ts_clipped[, 1]
    )
    p <- ggplot(df, aes(x = days)) +
      geom_line(aes(y = original, color = "Original"), linetype = "solid") +
      geom_line(aes(y = clipped, color = "Clipped"), linetype = "dashed") +
      xlab("time [days]") +
      ylab("value") +
      ggtitle("original vs. clipped timeseries") +
      scale_color_manual(values = c(Original = "blue", Clipped = "red")) +
      scale_linetype_manual(values = c(Original = "solid", Clipped = "dashed")) +
      labs(color = "Timeseries")
  }

  output$plot <- renderPlot({
    print(plotInput())
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
