library(shiny)
library(DT)
library(ggplot2)
library(reshape2)

ui <- fluidPage(
  titlePanel("Dynamic Connection Matrix Generator"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("var_name", "Enter Variable Name:", ""),
      actionButton("add_var", "Add Variable"),
      hr(),
      uiOutput("matrix_inputs")
    ),
    
    mainPanel(
      h3("Generated Matrix"),
      DTOutput("matrix_table"),
      h3("Heatmap"),
      plotOutput("heatmap_plot")
    )
  )
)

server <- function(input, output, session) {
  variables <- reactiveVal(character())  # Stores variable names
  
  observeEvent(input$add_var, {
    if (input$var_name != "" && !(input$var_name %in% variables())) {
      variables(c(variables(), input$var_name))
    }
    updateTextInput(session, "var_name", value = "")
  })
  
  output$matrix_inputs <- renderUI({
    var_list <- variables()
    if (length(var_list) < 2) return(NULL)
    
    matrix_inputs <- lapply(seq_along(var_list), function(i) {
      lapply(seq_along(var_list), function(j) {
        if (i != j) {
          numericInput(
            inputId = paste0("conn_", i, "_", j),
            label = paste(var_list[i], "<->", var_list[j]),
            value = 0, min = -1, max = 1, step = 0.1
          )
        }
      })
    })
    do.call(tagList, unlist(matrix_inputs, recursive = FALSE))
  })
  
  output$matrix_table <- renderDT({
    var_list <- variables()
    if (length(var_list) < 2) return(NULL)
    
    mat <- matrix(0, nrow = length(var_list), ncol = length(var_list), 
                  dimnames = list(var_list, var_list))
    
    for (i in seq_along(var_list)) {
      for (j in seq_along(var_list)) {
        if (i != j) {
          input_id <- paste0("conn_", i, "_", j)
          if (!is.null(input[[input_id]])) {
            mat[i, j] <- input[[input_id]]
          }
        }
      }
    }
    datatable(mat, options = list(pageLength = 10))
  })
  
  output$heatmap_plot <- renderPlot({
    var_list <- variables()
    if (length(var_list) < 2) return(NULL)
    
    mat <- matrix(0, nrow = length(var_list), ncol = length(var_list), 
                  dimnames = list(var_list, var_list))
    
    for (i in seq_along(var_list)) {
      for (j in seq_along(var_list)) {
        if (i != j) {
          input_id <- paste0("conn_", i, "_", j)
          if (!is.null(input[[input_id]])) {
            mat[i, j] <- input[[input_id]]
          }
        }
      }
    }
    
    melted_mat <- melt(mat)
    ggplot(melted_mat, aes(Var1, Var2, fill = value)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                           midpoint = 0, limit = c(-1, 1), space = "Lab") +
      theme_minimal() +
      labs(title = "Connection Strength Heatmap", fill = "Strength") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
}

shinyApp(ui, server)
