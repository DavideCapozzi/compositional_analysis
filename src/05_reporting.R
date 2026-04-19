# src/05_reporting.R
# ==============================================================================
# STEP 05: REPORTING
# Description: Compiles all analytical outputs into a single Master Excel 
#              report and dedicated Meta-Analysis summaries.
# ==============================================================================

source("R/utils_io.R")
source("R/modules_qc.R")
source("R/modules_viz.R")

message("\n=== PIPELINE STEP 5: REPORTING ===")

# 1. Initialization
# ------------------------------------------------------------------------------
config <- load_config("config/global_params.yml")
results_dir <- file.path(config$output_root, "results_analysis")
meta_dir <- file.path(results_dir, "Meta_Analysis")

if (!dir.exists(meta_dir)) dir.create(meta_dir, recursive = TRUE)

# Loading Payloads and Data
data_file <- file.path(config$output_root, "01_data_processing", "data_processed.rds")
if (!file.exists(data_file)) stop("[Fatal] Step 01 output not found.")
DATA <- readRDS(data_file)

payload03 <- readRDS(file.path(results_dir, "step03_payload.rds"))
payload04 <- readRDS(file.path(results_dir, "step04_payload.rds"))

# Initialize Master Workbook
wb_master <- createWorkbook()

# 2. MASTER ANALYTICAL REPORT (Global Level)
# ------------------------------------------------------------------------------
addWorksheet(wb_master, "QC_Overview")
writeData(wb_master, "QC_Overview", "Quality Control Summary", startRow = 1)

# Global PERMANOVA results
if (!is.null(payload03$global$permanova)) {
  addWorksheet(wb_master, "Global_PERMANOVA")
  writeData(wb_master, "Global_PERMANOVA", payload03$global$permanova)
}

# sPLS-DA Model Performance (BER/Validation)
if (!is.null(payload03$global$splsda_perf)) {
  addWorksheet(wb_master, "Global_sPLSDA_Performance")
  writeData(wb_master, "Global_sPLSDA_Performance", payload03$global$splsda_perf)
}

# sPLS-DA Discriminant Markers
if (!is.null(payload03$global$splsda_drivers)) {
  addWorksheet(wb_master, "Global_sPLSDA_Drivers")
  writeData(wb_master, "Global_sPLSDA_Drivers", payload03$global$splsda_drivers)
}

# 3. META-ANALYSIS: EDGE INTERSECTIONS
# ------------------------------------------------------------------------------
if (length(payload04$meta_analysis$diff_edges_list) >= 2) {
  message("   [Report] Generating Meta-Analysis Edge Intersections...")
  wb_meta <- createWorkbook()
  
  diff_list <- payload04$meta_analysis$diff_edges_list
  all_edges <- unique(unlist(diff_list))
  
  intersection_mat <- data.frame(Edge_ID = all_edges, stringsAsFactors = FALSE)
  for (scen_id in names(diff_list)) {
    intersection_mat[[scen_id]] <- intersection_mat$Edge_ID %in% diff_list[[scen_id]]
  }
  
  intersection_mat$Sum_Hits <- rowSums(intersection_mat[, -1])
  intersection_mat <- intersection_mat[order(intersection_mat$Sum_Hits, decreasing = TRUE), ]
  
  addWorksheet(wb_meta, "Edge_Intersections")
  writeData(wb_meta, "Edge_Intersections", intersection_mat)
  
  saveWorkbook(wb_meta, file.path(meta_dir, "Differential_Intersections_List.xlsx"), overwrite = TRUE)
}

# 4. META-ANALYSIS: INVERTED EDGES (Phenotypic Rewiring)
# ------------------------------------------------------------------------------
message("   [Report] Checking for Inverted Edges across all scenarios...")
wb_inverted <- createWorkbook()
inverted_master_list <- list()

for (scen_id in names(payload04$scenarios)) {
  net_data <- payload04$scenarios[[scen_id]]
  
  scen_conf <- purrr::keep(config$analysis_scenarios, ~ .x$id == scen_id)
  if (length(scen_conf) == 0) next
  scen_conf <- scen_conf[[1]]
  
  if (!is.null(net_data$edges_table)) {
    # Filter edges with sign inversion
    inv_edges <- net_data$edges_table %>% 
      dplyr::filter(Edge_Category == "Inverted")
    
    if (nrow(inv_edges) > 0) {
      message(sprintf("      -> Found %d inverted edges in scenario: %s", nrow(inv_edges), scen_id))
      
      # Scenario-specific formatting (Keeping Node1 and Node2 here for detail)
      inv_formatted <- inv_edges %>%
        dplyr::rename(
          Cor_Ctrl = dplyr::matches(paste0("^Spearman_", scen_conf$control_label, "$")),
          Weight_Ctrl = dplyr::matches(paste0("^Pcor_", scen_conf$control_label, "$")),
          Cor_Case = dplyr::matches(paste0("^Spearman_", scen_conf$case_label, "$")),
          Weight_Case = dplyr::matches(paste0("^Pcor_", scen_conf$case_label, "$"))
        ) %>%
        dplyr::select(Node1, Node2, Cor_Ctrl, Weight_Ctrl, Cor_Case, Weight_Case, 
                      Diff_Score, P_Value, FDR, 
                      dplyr::matches("^Mech_"))
      
      # Add detailed sheet for this scenario
      sh_inv_name <- safe_sheet_name(scen_id)
      addWorksheet(wb_inverted, sh_inv_name)
      writeData(wb_inverted, sh_inv_name, inv_formatted)
      
      # Prepare data for Master Summary collection
      inv_summary_entry <- inv_formatted %>%
        dplyr::mutate(
          Scenario = scen_id,
          Edge_ID = paste(pmin(Node1, Node2), pmax(Node1, Node2), sep = "~")
        ) %>%
        dplyr::select(Edge_ID, Node1, Node2, Scenario)
      
      inverted_master_list[[scen_id]] <- inv_summary_entry
    }
  }
}

# Finalize Inverted Report
if (length(inverted_master_list) > 0) {
  # Master Summary: Removing Node1 and Node2 as they are redundant with Edge_ID
  master_inv_df <- dplyr::bind_rows(inverted_master_list) %>%
    dplyr::group_by(Edge_ID) %>%
    dplyr::summarise(
      Scenarios_Inverted = paste(Scenario, collapse = ", "),
      Count_Inverted = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(Count_Inverted), Edge_ID)
  
  # Create the Master Summary sheet
  addWorksheet(wb_inverted, "Master_Inverted_Summary")
  writeData(wb_inverted, "Master_Inverted_Summary", master_inv_df)
  
  # Reorder to place the Master summary as the first sheet
  n_sheets <- length(names(wb_inverted))
  if (n_sheets > 1) {
    worksheetOrder(wb_inverted) <- c(n_sheets, 1:(n_sheets - 1))
  }
  
  inv_report_path <- file.path(meta_dir, "Inverted_Edges_MetaAnalysis_Report.xlsx")
  saveWorkbook(wb_inverted, inv_report_path, overwrite = TRUE)
  message(sprintf("   [Output] Inverted Edges Report saved: %s", inv_report_path))
} else {
  message("   [Info] No inverted edges found. Skipping specific report.")
}

# 5. Save Final Master Analytical Report
# ------------------------------------------------------------------------------
master_report_path <- file.path(results_dir, "Master_Analytical_Report.xlsx")
saveWorkbook(wb_master, master_report_path, overwrite = TRUE)

message(sprintf("\n[Final] Reporting complete. Master file: %s", master_report_path))
message("=== STEP 5 COMPLETE ===\n")