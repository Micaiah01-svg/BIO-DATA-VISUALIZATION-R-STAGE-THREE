# ======================================================
#  Cell Type Snapshot Explorer (Shiny App)
#  Author: Micaiah Adeoluwa Adedeji
#  Description: Minimal reproducible scRNA-seq explorer
# ======================================================

library(shiny)
library(ggplot2)
library(dplyr)
library(viridis)


# ----------------------------
# 1. Helper Functions
# ----------------------------

# Transparent color helper
t_col <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(
    rgb.val[1], rgb.val[2], rgb.val[3],
    max = 255, alpha = (100 - percent) * 255 / 100,
    names = name
  )
  invisible(t.col)
}

# Scale values to 0–100
scale_0_100 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(50, length(x)))
  (x - rng[1]) / (rng[2] - rng[1]) * 100
}

# Compute per-gene statistics
compute_gene_stats <- function(expr_mat, meta_df, target_ct) {
  in_cells <- meta_df$cell_id[meta_df$cell_type == target_ct]
  out_cells <- meta_df$cell_id[meta_df$cell_type != target_ct]
  xin <- expr_mat[in_cells, , drop = FALSE]
  xout <- expr_mat[out_cells, , drop = FALSE]
  det_in <- colMeans(xin > 0)
  det_out <- colMeans(xout > 0)
  mean_in <- colMeans(xin)
  mean_out <- colMeans(xout)
  diff <- mean_in - mean_out
  
  data.frame(
    gene = colnames(expr_mat),
    det_in = det_in,
    det_out = det_out,
    mean_in = mean_in,
    mean_out = mean_out,
    diff = diff,
    stringsAsFactors = FALSE
  )
}

# Select marker gene: max diff then max det_in
pick_marker_gene <- function(gene_stats_df) {
  gene_stats_df <- gene_stats_df %>%
    arrange(desc(diff), desc(det_in))
  gene_stats_df$gene[1]
}

# ----------------------------
# 2. Load Data (Local or from URL)
# ----------------------------

# Define local paths and GitHub URLs
expr_path <- "expression_matrix.csv"
meta_path <- "cell_metadata.csv"
umap_path <- "umap_coordinates.csv"

expr_url <- "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/sc_synthetic/expression_matrix.csv"
meta_url <- "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/sc_synthetic/cell_metadata.csv"
umap_url <- "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/sc_synthetic/umap_coordinates.csv"

# Smart loader: reads local if exists, otherwise from GitHub
safe_read <- function(local_path, url, ...) {
  if (file.exists(local_path)) {
    message("✅ Reading local file: ", local_path)
    read.csv(local_path, ...)
  } else {
    message("🌐 Local file not found — downloading from GitHub: ", basename(url))
    read.csv(url, ...)
  }
}

# Load files
expr_mat <- safe_read(expr_path, expr_url, row.names = 1, check.names = FALSE)
meta_df <- safe_read(meta_path, meta_url, stringsAsFactors = FALSE)
umap.coord <- safe_read(umap_path, umap_url, stringsAsFactors = FALSE)

# Ensure all cell IDs align
common_ids <- Reduce(intersect, list(rownames(expr_mat), meta_df$cell_id, umap.coord$cell_id))
expr_mat <- expr_mat[common_ids, , drop = FALSE]
meta_df <- meta_df %>% filter(cell_id %in% common_ids)
umap.coord <- umap.coord %>% filter(cell_id %in% common_ids)

# Reorder to align cells
expr_mat <- expr_mat[meta_df$cell_id, , drop = FALSE]

# ----------------------------
# 3. Shiny UI
# ----------------------------

ui <- fluidPage(
  titlePanel("🧬 Cell Type Snapshot Explorer (scRNA-seq)"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("cell_type", "Select cell type:",
                  choices = sort(unique(meta_df$cell_type))),
      br(),
      h5("Outputs: UMAP visualization, marker gene, and gene statistics."),
      width = 3
    ),
    mainPanel(
      plotOutput("umapPlot", height = "500px"),
      br(),
      textOutput("markerInfo"),
      br(),
      tableOutput("geneStats")
    )
  )
)

# ----------------------------
# 4. Shiny Server
# ----------------------------

server <- function(input, output, session) {
  
  # Reactive: compute gene stats for selected cell type
  gene_stats <- reactive({
    req(input$cell_type)
    compute_gene_stats(expr_mat, meta_df, input$cell_type)
  })
  
  # Reactive: identify marker gene
  marker_gene <- reactive({
    req(gene_stats())
    pick_marker_gene(gene_stats())
  })
  
  # Reactive: merge UMAP + metadata + expression
  umap_data <- reactive({
    req(marker_gene())
    df <- meta_df %>%
      left_join(umap.coord, by = "cell_id")
    df$expr <- expr_mat[df$cell_id, marker_gene()]
    df$expr_scaled <- scale_0_100(df$expr)
    df
  })
  
  # Marker info text
  output$markerInfo <- renderText({
    req(marker_gene(), gene_stats())
    gene <- marker_gene()
    diff_val <- gene_stats() %>% filter(gene == !!gene) %>% pull(diff)
    paste0("📍 Selected Cell Type: ", input$cell_type,
           " | Marker Gene: ", gene,
           " | Diff Score: ", round(diff_val, 3))
  })
  
  # UMAP plot: cells colored by marker expression
  output$umapPlot <- renderPlot({
    df <- umap_data()
    ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(aes(color = expr_scaled), size = 1.8) +
      scale_color_viridis_c(option = "plasma", name = paste0("Expr(", marker_gene(), ")")) +
      theme_minimal() +
      labs(title = paste("UMAP —", input$cell_type),
           subtitle = paste("Marker gene:", marker_gene())) +
      theme(
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10)
      )
  })
  
  # Gene statistics table (top 15 by diff)
  output$geneStats <- renderTable({
    head(gene_stats() %>% arrange(desc(diff)), 15)
  })
}

# ----------------------------
# 5. Run the App
# ----------------------------
shinyApp(ui = ui, server = server)




