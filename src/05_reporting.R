# src/05_reporting.R
# ==============================================================================
# STEP 05: MASTER REPORTING
# Description: Consolidates outputs from Step 03 and Step 04 payloads into
#              a unified Excel Workbook, isolating I/O formatting.
# ==============================================================================

source("R/utils_io.R")

message("\n=== PIPELINE STEP 5: MASTER REPORT GENERATION ===")

config <- load_config("config/global_params.yml")
results_dir <- file.path(config$output_root, "results_analysis")

payload03_path <- file.path(results_dir, "step03_payload.rds")
payload04_path <- file.path(results_dir, "step04_payload.rds")

if (!file.exists(payload03_path) || !file.exists(payload04_path)) {
  stop("[Fatal] Required payloads missing. Ensure Steps 03 and 04 have completed.")
}

payload03 <- readRDS(payload03_path)
payload04 <- readRDS(payload04_path)

wb_master <- createWorkbook()

# 1. WRITE GLOBAL STATISTICS
message("   [Report] Compiling Global Statistics...")
if (!is.null(payload03$global$permanova)) {
  addWorksheet(wb_master, "Global_PERMANOVA")
  writeData(wb_master, "Global_PERMANOVA", payload03$global$permanova, rowNames = TRUE)
}
if (!is.null(payload03$global$dispersion)) {
  addWorksheet(wb_master, "Global_Dispersion")
  writeData(wb_master, "Global_Dispersion", payload03$global$dispersion, rowNames = TRUE)
}
if (!is.null(payload03$global$pairwise_dispersion)) {
  addWorksheet(wb_master, "All_Pairs_Dispersion")
  writeData(wb_master, "All_Pairs_Dispersion", payload03$global$pairwise_dispersion)
}
if (!is.null(payload03$global$splsda_drivers)) {
  addWorksheet(wb_master, "Global_sPLSDA_Drivers")
  writeData(wb_master, "Global_sPLSDA_Drivers", payload03$global$splsda_drivers)
}

# 2. WRITE SCENARIO METRICS
message("   [Report] Compiling Scenario Metrics...")
combined_drivers_list <- list()

for (scen_id in names(payload03$scenarios)) {
  scen_data <- payload03$scenarios[[scen_id]]
  net_data  <- payload04$scenarios[[scen_id]]
  
  if (!is.null(scen_data$dispersion)) {
    disp_sheet <- safe_sheet_name(paste0(scen_id, "_Disp"))
    addWorksheet(wb_master, disp_sheet); writeData(wb_master, disp_sheet, scen_data$dispersion, rowNames = TRUE)
  }
  if (!is.null(scen_data$permanova)) {
    perm_sheet <- safe_sheet_name(paste0(scen_id, "_Perm"))
    addWorksheet(wb_master, perm_sheet); writeData(wb_master, perm_sheet, scen_data$permanova, rowNames = TRUE)
  }
  if (!is.null(scen_data$drivers)) {
    drv_sheet <- safe_sheet_name(paste0(scen_id, "_Drv"))
    addWorksheet(wb_master, drv_sheet); writeData(wb_master, drv_sheet, scen_data$drivers)
    
    df <- scen_data$drivers
    df$Scenario <- scen_id 
    df_clean <- df %>%
      dplyr::rename_with(.fn = ~ "Weight_Case_VS_Control_PC1", .cols = dplyr::matches("^Weight_.*_PC1$")) %>%
      dplyr::rename_with(.fn = ~ "Weight_Case_VS_Control_PC2", .cols = dplyr::matches("^Weight_.*_PC2$")) %>%
      dplyr::select(Scenario, Marker, dplyr::matches("^Weight_Case_VS_Control_PC[12]$"))
    combined_drivers_list[[scen_id]] <- df_clean
  }
  
  if (!is.null(net_data$edges_table)) {
    net_sheet <- safe_sheet_name(paste0(scen_id, "_Net"))
    addWorksheet(wb_master, net_sheet); writeData(wb_master, net_sheet, net_data$edges_table)
  }
  if (!is.null(net_data$sig_edges) && nrow(net_data$sig_edges) > 0) {
    sig_sheet <- safe_sheet_name(paste0(scen_id, "_SigNet"))
    addWorksheet(wb_master, sig_sheet); writeData(wb_master, sig_sheet, net_data$sig_edges)
  }
  
  # --- Generate Local Topology Metrics Report ---
  if (!is.null(net_data$topo_ctrl) || !is.null(net_data$topo_case)) {
    # Extract labels safely from config
    scen_conf <- purrr::keep(config$analysis_scenarios, ~ .x$id == scen_id)[[1]]
    
    topo_xlsx_path <- file.path(results_dir, scen_id, paste0(scen_id, "_Topology_Metrics.xlsx"))
    wb_topo <- createWorkbook()
    
    if (!is.null(net_data$topo_ctrl)) {
      sh_name <- safe_sheet_name(paste0("Topology_", scen_conf$control_label))
      addWorksheet(wb_topo, sh_name); writeData(wb_topo, sh_name, net_data$topo_ctrl)
    }
    if (!is.null(net_data$topo_case)) {
      sh_name <- safe_sheet_name(paste0("Topology_", scen_conf$case_label))
      addWorksheet(wb_topo, sh_name); writeData(wb_topo, sh_name, net_data$topo_case)
    }
    if (!is.null(net_data$rewiring)) {
      addWorksheet(wb_topo, "Rewiring_Analysis"); writeData(wb_topo, "Rewiring_Analysis", net_data$rewiring)
    }
    if (!is.null(net_data$edges_table)) {
      addWorksheet(wb_topo, "Differential_Edges"); writeData(wb_topo, "Differential_Edges", net_data$edges_table)
    }
    
    saveWorkbook(wb_topo, topo_xlsx_path, overwrite = TRUE)
  }
}

if (length(combined_drivers_list) > 0) {
  combined_drivers <- dplyr::bind_rows(combined_drivers_list)
  if (nrow(combined_drivers) > 0) {
    addWorksheet(wb_master, "All_Drivers_Summary")
    writeData(wb_master, "All_Drivers_Summary", combined_drivers)
  }
}

# 3. META-ANALYSIS SUMMARY
message("   [Report] Compiling Meta-Analysis Summary...")
diff_edges_list <- payload04$meta_analysis$diff_edges_list

if (!is.null(diff_edges_list) && length(diff_edges_list) >= 2) {
  wb_meta <- createWorkbook()
  addWorksheet(wb_meta, "Summary")
  m_comb <- ComplexHeatmap::make_comb_mat(diff_edges_list)
  comb_names <- ComplexHeatmap::comb_name(m_comb)
  summary_df <- data.frame()
  
  for (nm in comb_names) {
    edges_in_comb <- ComplexHeatmap::extract_comb(m_comb, nm)
    if (length(edges_in_comb) > 0) {
      set_indices <- as.numeric(strsplit(nm, "")[[1]])
      involved_sets <- names(diff_edges_list)[which(set_indices == 1)]
      set_label <- paste(involved_sets, collapse = " & ")
      sheet_id <- paste0("Int_", nrow(summary_df) + 1)
      
      summary_df <- rbind(summary_df, data.frame(Intersection = set_label, Count = length(edges_in_comb), Sheet_Link = sheet_id))
      addWorksheet(wb_meta, sheet_id)
      details <- data.frame(Edge_ID = edges_in_comb) %>% separate(Edge_ID, into = c("Node_A", "Node_B"), sep = "~")
      writeData(wb_meta, sheet_id, details)
    }
  }
  writeData(wb_meta, "Summary", summary_df)
  meta_dir <- file.path(results_dir, "Meta_Analysis")
  saveWorkbook(wb_meta, file.path(meta_dir, "Differential_Intersections_List.xlsx"), overwrite = TRUE)
}

# 4. FINALIZE MASTER WORKBOOK
sheet_names <- names(wb_master)
if ("All_Drivers_Summary" %in% sheet_names) {
  target_idx <- which(sheet_names == "All_Drivers_Summary")
  other_idx <- setdiff(seq_along(sheet_names), target_idx)
  openxlsx::worksheetOrder(wb_master) <- c(target_idx, other_idx)
}

master_report_path <- file.path(results_dir, "Multi_Scenario_Analysis_Report.xlsx")
saveWorkbook(wb_master, master_report_path, overwrite = TRUE)
message(sprintf("   [Output] Master Report saved to: %s", master_report_path))
message("=== STEP 5 COMPLETE ===\n")