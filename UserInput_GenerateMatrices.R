library(shiny)
library(ggplot2)
library(reshape2)

ui <- fluidPage(
  titlePanel("Group & AR Path Matrix Generator"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("nvar", "Number of Variables (nvar):", value = 6, min = 1, step = 1),
      numericInput("pcon", "Proportion Contemporaneous (0-1):", value = 0.5, min = 0, max = 1, step = 0.1),
      numericInput("AR", "Auto-Regressive Coefficient:", value = 0.5),
      numericInput("dens", "Density (0-1):", value = 0.25, min = 0, max = 1, step = 0.05),
      numericInput("pgroup", "Proportion Group-Level (0-1):", value = 1, min = 0, max = 1, step = 0.1),
      numericInput("conb", "Contemporaneous Beta:", value = 0.3),
      numericInput("lagb", "Lagged Beta:", value = 0.3),
      textAreaInput("group_paths", "Group-Level Paths (one per line, format: row,col):", 
                    value = "3,5\n4,6\n6,1\n5,4", rows = 5),
      textAreaInput("ar_paths", "AR Paths (one per line, format: row,col):", 
                    value = "2,2\n4,4\n6,6", rows = 3),
      actionButton("generate", "Generate Matrices")
    ),
    
    mainPanel(
      h4("Lagged Matrix (Phi)"),
      plotOutput("phi_plot"),
      h4("Contemporaneous Matrix (A)"),
      plotOutput("a_plot")
    )
  )
)

# Define the Server
server <- function(input, output) {
  generate_matrices <- eventReactive(input$generate, {
    nvar <- input$nvar
    
    # Parse group-level paths from multiline input
    group_list <- unlist(strsplit(input$group_paths, "\n"))
    group_coords <- do.call(rbind, lapply(group_list, function(x) as.integer(unlist(strsplit(x, ",")))))
    
    # Parse AR paths from multiline input
    ar_list <- unlist(strsplit(input$ar_paths, "\n"))
    ar_coords <- do.call(rbind, lapply(ar_list, function(x) as.integer(unlist(strsplit(x, ",")))))
    
    # Create empty matrices
    A <- matrix(0, nvar, nvar)
    Phi <- matrix(0, nvar, nvar)
    
    # Split group paths into contemporaneous and lagged
    if (nrow(group_coords) > 0) {
      split_point <- round(input$pcon * nrow(group_coords))
      grp_con <- if (split_point > 0) group_coords[1:split_point, , drop = FALSE] else matrix(nrow = 0, ncol = 2)
      grp_lag <- if (split_point < nrow(group_coords)) group_coords[(split_point + 1):nrow(group_coords), , drop = FALSE] else matrix(nrow = 0, ncol = 2)
      
      # Fill matrices
      if (nrow(grp_con) > 0) {
        for (i in 1:nrow(grp_con)) A[grp_con[i,1], grp_con[i,2]] <- input$conb
      }
      if (nrow(grp_lag) > 0) {
        for (i in 1:nrow(grp_lag)) Phi[grp_lag[i,1], grp_lag[i,2]] <- input$lagb
      }
    }
    
    if (nrow(ar_coords) > 0) {
      for (i in 1:nrow(ar_coords)) Phi[ar_coords[i,1], ar_coords[i,2]] <- input$AR
    }
    
    list(A = A, Phi = Phi)
  })
  
  draw_matrix_plot <- function(mat, title, lag_labels = FALSE) {
    n <- nrow(mat)
    rownames(mat) <- paste0("V", 1:n)
    colnames(mat) <- if (lag_labels) paste0("V", 1:n, "_lag") else paste0("V", 1:n)
    melted <- melt(mat)
    melted$label <- ifelse(melted$value != 0, round(melted$value, 2), "")
    
    ggplot(melted, aes(Var2, Var1, fill = value)) +
      geom_tile(color = "black") +
      geom_text(aes(label = label), size = 4) +
      scale_fill_gradient(low = "white", high = if (lag_labels) "darkred" else "steelblue") +
      labs(title = title) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(hjust = 0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.text.x.top = element_text(),
        axis.ticks.length = unit(0, "pt")
      ) +
      scale_x_discrete(position = "top")
  }
  
  output$a_plot <- renderPlot({
    req(generate_matrices())
    draw_matrix_plot(generate_matrices()$A, "Contemporaneous Matrix (A)", lag_labels = FALSE)
  })
  
  output$phi_plot <- renderPlot({
    req(generate_matrices())
    draw_matrix_plot(generate_matrices()$Phi, "Lagged Matrix (Phi)", lag_labels = TRUE)
  })
}

shinyApp(ui = ui, server = server)
