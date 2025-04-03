library(shiny)
library(ggplot2)
library(reshape2)

ui <- fluidPage(
  titlePanel("Group & AR Path Matrix Generator with Individual Paths and Noise"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("nvar", "Number of Variables (nvar):", value = 6, min = 1, step = 1),
      numericInput("AR", "Auto-Regressive Coefficient:", value = 0.5),
      numericInput("conb", "Contemporaneous Beta:", value = 0.3),
      numericInput("lagb", "Lagged Beta:", value = 0.3),
      textAreaInput("ar_paths", "AR Paths (row,col):", 
                    value = "2,2\n4,4\n6,6", rows = 3),
      
      fluidRow(
        column(6, textAreaInput("group_con_paths", "Contemporaneous Group-Level Paths (row,col):", 
                                value = "3,5\n4,6", rows = 4)),
        column(6, textAreaInput("group_lag_paths", "Lagged Group-Level Paths (row,col):", 
                                value = "6,1\n5,4", rows = 4))
      ),
      
      fluidRow(
        column(6, textAreaInput("indiv_con_paths", "Contemporaneous Individual-Level Paths (row,col):", 
                                value = "1,2", rows = 3)),
        column(6, textAreaInput("indiv_lag_paths", "Lagged Individual-Level Paths (row,col):", 
                                value = "2,4", rows = 3))
      ),
      
      numericInput("noise_sd", "Standard Deviation for Gaussian Noise:", value = 0.01, min = 0, step = 0.005),
      numericInput("ts_length", "Length of Time Series (t):", value = 100, min = 1, step = 1),
      actionButton("generate", "Generate Matrices and Simulate Time Series")
    ),
    
    mainPanel(
      h4("Lagged Matrix (Phi)"),
      plotOutput("phi_plot"),
      h4("Contemporaneous Matrix (A)"),
      plotOutput("a_plot"),
      h4("Simulated Time Series Plot"),
      plotOutput("ts_plot")
    )
  )
)

server <- function(input, output) {
  generate_matrices <- eventReactive(input$generate, {
    nvar <- input$nvar
    
    parse_coords <- function(txt) {
      if (nzchar(txt)) {
        rows <- unlist(strsplit(txt, "\n"))
        do.call(rbind, lapply(rows, function(x) as.integer(unlist(strsplit(x, ",")))))
      } else {
        matrix(nrow = 0, ncol = 2)
      }
    }
    
    group_con_coords <- parse_coords(input$group_con_paths)
    group_lag_coords <- parse_coords(input$group_lag_paths)
    ar_coords <- parse_coords(input$ar_paths)
    indiv_con_coords <- parse_coords(input$indiv_con_paths)
    indiv_lag_coords <- parse_coords(input$indiv_lag_paths)
    
    A <- matrix(0, nvar, nvar)
    Phi <- matrix(0, nvar, nvar)
    
    if (nrow(group_con_coords) > 0) {
      for (i in 1:nrow(group_con_coords)) A[group_con_coords[i,1], group_con_coords[i,2]] <- input$conb
    }
    if (nrow(group_lag_coords) > 0) {
      for (i in 1:nrow(group_lag_coords)) Phi[group_lag_coords[i,1], group_lag_coords[i,2]] <- input$lagb
    }
    if (nrow(ar_coords) > 0) {
      for (i in 1:nrow(ar_coords)) Phi[ar_coords[i,1], ar_coords[i,2]] <- input$AR
    }
    if (nrow(indiv_con_coords) > 0) {
      for (i in 1:nrow(indiv_con_coords)) A[indiv_con_coords[i,1], indiv_con_coords[i,2]] <- input$conb
    }
    if (nrow(indiv_lag_coords) > 0) {
      for (i in 1:nrow(indiv_lag_coords)) Phi[indiv_lag_coords[i,1], indiv_lag_coords[i,2]] <- input$lagb
    }
    
    A[which(A != 0, arr.ind = TRUE)] <- A[which(A != 0, arr.ind = TRUE)] + 
      rnorm(n = length(which(A != 0)), mean = 0, sd = input$noise_sd)
    
    Phi[which(Phi != 0, arr.ind = TRUE)] <- Phi[which(Phi != 0, arr.ind = TRUE)] + 
      rnorm(n = length(which(Phi != 0)), mean = 0, sd = input$noise_sd)
    
    list(A = A, Phi = Phi)
  })
  
  simulate_timeseries <- eventReactive(input$generate, {
    mat <- generate_matrices()
    A <- mat$A
    Phi <- mat$Phi
    v <- nrow(A)
    t <- input$ts_length
    st <- t + 50
    
    noise <- matrix(rnorm(v * st, 0, 1), v)
    I <- diag(v)
    time <- matrix(0, nrow = v, ncol = (st + 1))
    time1 <- matrix(0, nrow = v, ncol = st)
    
    for (i in 1:st) {
      time1[, i] <- solve(I - A) %*% (Phi %*% time[, i] + noise[, i])
      time[, i + 1] <- time1[, i]
    }
    
    time1 <- time1[, 51:(50 + t)]
    series <- t(time1)
    colnames(series) <- paste0("V", 1:v)
    series_df <- as.data.frame(series)
    series_df$Time <- 1:t
    series_df
  })
  
  draw_matrix_plot <- function(mat, title, lag_labels = FALSE) {
    n <- nrow(mat)
    rownames(mat) <- paste0("V", 1:n)
    colnames(mat) <- if (lag_labels) paste0("V", 1:n, "_lag") else paste0("V", 1:n)
    melted <- melt(mat)
    melted$label <- ifelse(melted$value != 0, round(melted$value, 2), "")
    melted$Var1 <- factor(melted$Var1, levels = rev(levels(melted$Var1)))
    
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
  
  output$ts_plot <- renderPlot({
    req(simulate_timeseries())
    df <- simulate_timeseries()
    df_melt <- melt(df, id.vars = "Time")
    ggplot(df_melt, aes(x = Time, y = value, color = variable)) +
      geom_line() +
      labs(title = "Simulated Time Series", y = "Value", x = "Time") +
      theme_minimal()
  })
}

shinyApp(ui = ui, server = server)
