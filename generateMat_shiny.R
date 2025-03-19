library(shiny)
library(DT)
library(ggplot2)
library(reshape2)

ui <- fluidPage(
  titlePanel("Randomized Connection Matrix Generator"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("var_name", "Enter Variable Name:", ""),
      actionButton("add_var", "Add Variable"),
      hr(),
      textInput("vector_positions", "Enter Fixed Connections (e.g., A-B, C-D):", "A-B, C-D"),
      numericInput("fixed_value", "Set Fixed Value for Connections:", 0.8, step = 0.1),
      actionButton("apply_fixed", "Apply Fixed Values"),
      hr(),
      numericInput("rand_mean", "Mean for Random Values:", 0.2, step = 0.1),
      numericInput("rand_sd", "Standard Deviation for Random Values:", 0.1, step = 0.1),
      hr(),
      actionButton("generate_matrix", "Generate Matrix & Heatmap", class = "btn-primary")
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
  # Default variables
  variables <- reactiveVal(c("A", "B", "C", "D"))
  connections <- reactiveVal(data.frame(from = c("A", "C"), to = c("B", "D"), value = c(0.8, 0.8), stringsAsFactors = FALSE))
  
  # Matrix generation logic (asymmetric)
  matrix_result <- eventReactive(input$generate_matrix, {
    var_list <- variables()
    if (length(var_list) < 2) return(NULL)
    
    n <- length(var_list)
    mat <- matrix(rnorm(n^2, mean = input$rand_mean, sd = input$rand_sd), nrow = n, ncol = n)
    diag(mat) <- 0  # Set diagonal to zero (no self-connections)
    
    # Apply fixed connections
    conn_df <- connections()
    if (nrow(conn_df) > 0) {
      for (i in 1:nrow(conn_df)) {
        from <- conn_df$from[i]
        to <- conn_df$to[i]
        value <- conn_df$value[i]
        row_idx <- which(var_list == from)
        col_idx <- which(var_list == to)
        if (length(row_idx) > 0 && length(col_idx) > 0) {
          mat[row_idx, col_idx] <- value  # Asymmetric: only (i, j), not (j, i)
        }
      }
    }
    
    colnames(mat) <- rownames(mat) <- var_list
    return(mat)
  })
  
  observeEvent(input$add_var, {
    if (input$var_name != "" && !(input$var_name %in% variables())) {
      variables(c(variables(), input$var_name))
    }
    updateTextInput(session, "var_name", value = "")
  })
  
  observeEvent(input$apply_fixed, {
    var_list <- variables()
    if (length(var_list) < 2) return(NULL)
    
    conn_text <- gsub("\\s+", "", input$vector_positions)  # Remove spaces
    conn_pairs <- unlist(strsplit(conn_text, ","))  # Split into pairs
    
    new_connections <- do.call(rbind, lapply(conn_pairs, function(pair) {
      nodes <- unlist(strsplit(pair, "-"))
      if (length(nodes) == 2 && all(nodes %in% var_list)) {
        data.frame(from = nodes[1], to = nodes[2], value = input$fixed_value, stringsAsFactors = FALSE)
      } else {
        NULL
      }
    }))
    
    if (!is.null(new_connections) && nrow(new_connections) > 0) {
      connections(rbind(connections(), new_connections))
    }
  })
  
  output$matrix_table <- renderDT({
    mat <- matrix_result()
    if (is.null(mat)) return(NULL)
    datatable(round(mat, 2), options = list(pageLength = 10))
  })
  
  output$heatmap_plot <- renderPlot({
    mat <- matrix_result()
    if (is.null(mat)) return(NULL)
    
    melted_mat <- melt(mat)
    colnames(melted_mat) <- c("Var1", "Var2", "value")
    
    # Get min, max, and midpoint dynamically
    min_val <- min(mat)
    max_val <- max(mat)
    mid_val <- mean(mat)
    
    ggplot(melted_mat, aes(Var1, Var2, fill = value)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                           midpoint = mid_val, limit = c(min_val, max_val), space = "Lab") +
      theme_minimal() +
      labs(title = "Asymmetric Connection Strength Heatmap", fill = "Strength") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Auto-generate matrix on startup
  observe({
    isolate({
      input$generate_matrix
      matrix_result()
    })
  })
}

shinyApp(ui, server)
