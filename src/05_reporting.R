# src/05_reporting.R
# ==============================================================================
# STEP 05: FINAL REPORT GENERATION
# Description: Consolidates results from all steps into Excel reports.
# ==============================================================================

source("R/utils_io.R")
source("R/workflows.R")

message("\n=== PIPELINE STEP 5: REPORTING ===")

config <- load_config("config/global_params.yml")
results_dir <- file.path(config$output_root, "results_analysis")

# 1. LOAD PAYLOADS
# ------------------------------------------------------------------------------
payload03_path <- file.path(results_dir, "step03_payload.rds")
payload04_path <- file.path(results_dir, "step04_payload.rds")

if (!file.exists(payload03_path)) stop("Step 03 payload not found.")
payload03 <- readRDS(payload03_path)

payload04 <- if (file.exists(payload04_path)) readRDS(payload04_path) else NULL

# 2. GENERATE SCENARIO REPORTS
# ------------------------------------------------------------------------------
message("   [Report] Generating Multi-Scenario Analysis Report...")

wb_scenarios <- createWorkbook()

# Sheet: Global Stats
if (!is.null(payload03$global)) {
  addWorksheet(wb_scenarios, "Global_Stats")
  curr_row <- 1
  
  if (!is.null(payload03$global$permanova)) {
    writeData(wb_scenarios, "Global_Stats", "GLOBAL PERMANOVA:", startRow = curr_row)
    writeData(wb_scenarios, "Global_Stats", payload03$global$permanova, startRow = curr_row + 1)
    curr_row <- curr_row + nrow(payload03$global$permanova) + 3
  }
  
  if (!is.null(payload03$global$splsda_drivers)) {
    writeData(wb_scenarios, "Global_Stats", "GLOBAL sPLS-DA TOP DRIVERS:", startRow = curr_row)
    writeData(wb_scenarios, "Global_Stats", payload03$global$splsda_drivers, startRow = curr_row + 1)
  }
}

# Sheets: Scenario Details (Hypothesis & Drivers)
if (!is.null(payload03$scenarios)) {
  for (scen_id in names(payload03$scenarios)) {
    scen_dat <- payload03$scenarios[[scen_id]]
    sh_name <- substr(scen_id, 1, 31)
    addWorksheet(wb_scenarios, sh_name)
    
    curr_row <- 1
    if (!is.null(scen_dat$permanova)) {
      writeData(wb_scenarios, sh_name, "PERMANOVA:", startRow = curr_row)
      writeData(wb_scenarios, sh_name, scen_dat$permanova, startRow = curr_row + 1)
      curr_row <- curr_row + nrow(scen_dat$permanova) + 3
    }
    
    if (!is.null(scen_dat$drivers)) {
      writeData(wb_scenarios, sh_name, "TOP DISCRIMINANT MARKERS (sPLS-DA):", startRow = curr_row)
      writeData(wb_scenarios, sh_name, scen_dat$drivers, startRow = curr_row + 1)
    }
  }
}

saveWorkbook(wb_scenarios, file.path(results_dir, "Multi_Scenario_Analysis_Report.xlsx"), overwrite = TRUE)

# 3. META-ANALYSIS: FOCUSED INVERTED EDGES
# ------------------------------------------------------------------------------
# We replace the old intersection list with a descriptive Inverted Switches report
if (!is.null(payload04$scenarios)) {
  message("   [Report] Generating Focused Inverted Edges Summary...")
  
  inverted_report_data <- create_inverted_edges_report(payload04$scenarios, config)
  
  if (!is.null(inverted_report_data)) {
    wb_inverted <- createWorkbook()
    
    for (sheet_name in names(inverted_report_data)) {
      addWorksheet(wb_inverted, sheet_name)
      writeData(wb_inverted, sheet_name, inverted_report_data[[sheet_name]])
      
      # Formatting for readability
      addStyle(wb_inverted, sheet_name, createStyle(textDecoration = "bold"), rows = 1, cols = 1:ncol(inverted_report_data[[sheet_name]]))
      setColWidths(wb_inverted, sheet_name, cols = 1:ncol(inverted_report_data[[sheet_name]]), widths = "auto")
    }
    
    meta_dir <- file.path(results_dir, "Meta_Analysis")
    if (!dir.exists(meta_dir)) dir.create(meta_dir, recursive = TRUE)
    
    inv_file <- file.path(meta_dir, "Inverted_Edges_MetaAnalysis_Report.xlsx")
    saveWorkbook(wb_inverted, inv_file, overwrite = TRUE)
    message(sprintf("   [Output] Focused Inverted Analysis saved: %s", basename(inv_file)))
  }
}

message("=== STEP 5 COMPLETE ===\n")