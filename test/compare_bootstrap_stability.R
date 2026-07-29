# test/compare_bootstrap_stability.R
# ==============================================================================
# !! STALE — DOES NOT RUN (verified 2026-07-29) !!
# This script reads `Multi_Scenario_Analysis_Report.xlsx`, which the pipeline no longer
# produces; src/05_reporting.R now writes `Master_Analytical_Report.xlsx`. It also points
# at `results_longsurvivors_onecore/`, which does not exist. Kept only as a sketch of the
# stability-comparison approach — rewrite the paths before using it.
# ==============================================================================
# BOOTSTRAP STABILITY DIAGNOSTICS
# Description: Compares two Multi_Scenario_Analysis_Report.xlsx outputs 
#              to evaluate RNG/Parallelization variance on edge stability.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(ggplot2)
  library(gridExtra)
})

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================
# Define the paths to your two execution reports
# (Adjust these paths according to your local directory structure)
file_run_A <- "results_longsurvivors/results_analysis/Multi_Scenario_Analysis_Report.xlsx"
file_run_B <- "results_longsurvivors_onecore/results_analysis/Multi_Scenario_Analysis_Report.xlsx"

label_A <- "Run_1_Core"
label_B <- "Run_15_Cores"

out_pdf <- "Bootstrap_Stability_Diagnostics.pdf"
out_csv <- "Bootstrap_Metrics_Summary.csv"

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

#' Create undirected, alphabetically sorted Edge IDs
create_edge_id <- function(df) {
  # Handle both 'Node1/Node2' and 'Source/Target' naming conventions
  col_n1 <- if("Node1" %in% names(df)) "Node1" else "Source"
  col_n2 <- if("Node2" %in% names(df)) "Node2" else "Target"
  
  df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Edge_ID = paste(sort(c(!!sym(col_n1), !!sym(col_n2))), collapse = "~")) %>%
    dplyr::ungroup()
}

# ==============================================================================
# 3. EXECUTION
# ==============================================================================

if (!file.exists(file_run_A)) stop(sprintf("File A not found: %s", file_run_A))
if (!file.exists(file_run_B)) stop(sprintf("File B not found: %s", file_run_B))

sheets_A <- getSheetNames(file_run_A)
sheets_B <- getSheetNames(file_run_B)

# Isolate Network sheets (they contain the full edge tables and StabFreq)
net_sheets <- intersect(
  sheets_A[grepl("_Net$", sheets_A)],
  sheets_B[grepl("_Net$", sheets_B)]
)

if (length(net_sheets) == 0) {
  stop("No common '_Net' sheets found between the two reports.")
}

message(sprintf("[Diagnostics] Found %d matching network scenarios.", length(net_sheets)))

pdf(out_pdf, width = 12, height = 8)
summary_metrics <- list()

for (sheet in net_sheets) {
  message(sprintf("   -> Analyzing Scenario: %s", sheet))
  
  # Read Data
  df_A <- read.xlsx(file_run_A, sheet = sheet)
  df_B <- read.xlsx(file_run_B, sheet = sheet)
  
  # Standardize IDs
  df_A <- create_edge_id(df_A)
  df_B <- create_edge_id(df_B)
  
  # Identify Continuous Stability Columns dynamically
  stab_cols <- grep("^StabFreq_", names(df_A), value = TRUE)
  
  if (length(stab_cols) == 0) {
    message("      [Skip] No continuous stability columns (StabFreq_) found.")
    next
  }
  
  # Select necessary columns and prefix them
  cols_to_keep <- c("Edge_ID", "Significant", "Edge_Category", stab_cols)
  
  sub_A <- df_A[, cols_to_keep] %>% dplyr::rename_with(~paste0(., "_A"), -Edge_ID)
  sub_B <- df_B[, cols_to_keep] %>% dplyr::rename_with(~paste0(., "_B"), -Edge_ID)
  
  # Outer Join to compare identical edges
  merged_df <- dplyr::inner_join(sub_A, sub_B, by = "Edge_ID")
  
  if (nrow(merged_df) == 0) {
    message("      [Skip] No matching edges found after join.")
    next
  }
  
  # Calculate Jaccard for Significant Edges
  sig_A <- merged_df$Edge_ID[merged_df$Significant_A == TRUE]
  sig_B <- merged_df$Edge_ID[merged_df$Significant_B == TRUE]
  intersection_n <- length(intersect(sig_A, sig_B))
  union_n <- length(union(sig_A, sig_B))
  jaccard_sig <- if(union_n > 0) intersection_n / union_n else NA
  
  plots <- list()
  
  for (scol in stab_cols) {
    col_A <- paste0(scol, "_A")
    col_B <- paste0(scol, "_B")
    
    # Calculate Continuous Error Metrics
    mae <- mean(abs(merged_df[[col_A]] - merged_df[[col_B]]), na.rm = TRUE)
    cor_val <- suppressWarnings(cor(merged_df[[col_A]], merged_df[[col_B]], method = "pearson", use = "complete.obs"))
    
    # Save Metrics
    summary_metrics[[paste(sheet, scol, sep="_")]] <- data.frame(
      Scenario = sheet,
      Metric_Column = scol,
      Total_Edges = nrow(merged_df),
      MAE = mae,
      Pearson_R = cor_val,
      Jaccard_Significant_Edges = jaccard_sig,
      stringsAsFactors = FALSE
    )
    
    # Plot 1: Scatterplot
    p_scatter <- ggplot(merged_df, aes_string(x = col_A, y = col_B)) +
      geom_point(alpha = 0.3, color = "steelblue") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      labs(
        title = paste("Scatter:", scol),
        subtitle = sprintf("MAE: %.4f | R: %.4f", mae, cor_val),
        x = paste(scol, label_A),
        y = paste(scol, label_B)
      ) +
      theme_bw()
    
    # Plot 2: Density of Differences
    merged_df$Delta <- merged_df[[col_A]] - merged_df[[col_B]]
    p_density <- ggplot(merged_df, aes(x = Delta)) +
      geom_density(fill = "darkorange", alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      labs(
        title = expression(Delta ~ "Distribution (A - B)"),
        subtitle = "Zero-centered implies unbiased algorithmic noise",
        x = paste("Difference in", scol),
        y = "Density"
      ) +
      theme_bw()
    
    plots[[paste0(scol, "_scatter")]] <- p_scatter
    plots[[paste0(scol, "_density")]] <- p_density
  }
  
  # Plot 3: Categorical Shift (Confusion Matrix / Barplot)
  cat_A <- merged_df %>% count(Edge_Category_A) %>% rename(Category = Edge_Category_A, N_A = n)
  cat_B <- merged_df %>% count(Edge_Category_B) %>% rename(Category = Edge_Category_B, N_B = n)
  cat_merged <- full_join(cat_A, cat_B, by = "Category") %>% replace(is.na(.), 0)
  cat_long <- cat_merged %>% pivot_longer(cols = c("N_A", "N_B"), names_to = "Run", values_to = "Count")
  
  p_bars <- ggplot(cat_long, aes(x = Category, y = Count, fill = Run)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("N_A" = "#4682B4", "N_B" = "#CD5C5C"), labels = c(label_A, label_B)) +
    labs(title = "Edge Classification Stability", x = "Category", y = "Count") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  plots[["categorical_bars"]] <- p_bars
  
  # Arrange and Print Page
  title_grob <- grid::textGrob(sprintf("Diagnostics for Scenario: %s", sheet), 
                               gp = grid::gpar(fontsize = 16, fontface = "bold"))
  grid.arrange(grobs = plots, ncol = 2, top = title_grob)
}

dev.off()

# Save Summary Table
if (length(summary_metrics) > 0) {
  df_summary <- do.call(rbind, summary_metrics)
  write.csv(df_summary, out_csv, row.names = FALSE)
  message(sprintf("\n[Success] Diagnostics completed. PDF saved to %s", out_pdf))
  message(sprintf("          Metrics summary saved to %s", out_csv))
} else {
  message("\n[Warn] No valid metrics extracted.")
}