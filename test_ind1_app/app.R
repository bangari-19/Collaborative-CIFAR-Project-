library(shiny)
library(ggplot2)

# Load dataset
data <- read.csv("~/Desktop/CIFAR/220429.FinalMatrices/1.csv") 

# Define UI ----
ui <- fluidPage(
  titlePanel("Timeseries for an individual"),
  
  sidebarLayout(
    sidebarPanel(
      # INPUTS
      numericInput(inputId = 'Size',
                   label = 'Number of variables',
                   value = 6),  # Default value
      
      numericInput(inputId = 'Steps',
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
      
      selectInput(inputId = 'Mask',
                  label = 'Apply mask',
                  choices = list('Yes' = 1, 'No' = 2), selected = 1)
    ),
    
    mainPanel(
      h1("Timeseries data"),
      plotOutput("timeseries_plot"),  # Output plot
      actionButton("generate", "Generate Data")  
    )
  )
)

# Define server logic --- 
server <- function(input, output) {
  data <- read.csv("~/Desktop/CIFAR/220429.FinalMatrices/1.csv")
  
  # Access inputs
  observeEvent(input$generate, {
    size <- input$Size
    steps <- input$Steps
    ampContemp <- input$ampContemp
    ampLagged <- input$ampLagged
    ampMeasure <- input$ampMeasure
    covContempName <- input$covContempName
    covLaggedName <- input$covLaggedName
    Mask <- input$Mask
    
    start_index <- size + 1
    stop <- ncol(data)
    
    # Generate time series data
    simulated_data <- generate_timeseries(
      start = rep(0, size),  
      len = steps,
      contempamp = ampContemp,
      contempmat = data[, start_index:stop],  
      contempcov = matrix(rnorm(size**2), nrow=size, ncol=size),
      lagamp = ampLagged,
      lagmat = data[, 1:size],  
      lagcov = matrix(rnorm(size**2), nrow=size, ncol=size),
      measurecov = ampMeasure * diag(size),  
      save = FALSE,
      debug = FALSE
    )
    
    # Create a data frame for ggplot
    df <- data.frame(
      days = 1:(steps-1),
      original = simulated_data[, 1]
    )
    
    # Plot the generated time series
    output$timeseries_plot <- renderPlot({
      ggplot(df, aes(x = days, y = original)) +
        geom_line(color = "blue", linetype = "solid") +
        xlab("time [days]") +
        ylab("value") +
        ggtitle("Original Timeseries") +
        theme_minimal()
    })
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)