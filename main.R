# main.R
# ==============================================================================
# MAIN PIPELINE ORCHESTRATOR
# Description: End-to-end execution of the Immunological Analysis Pipeline.
#              Orchestrates ETL, Visualization, and Multi-Scenario Analysis.
# ==============================================================================

# 1. Environment Setup
rm(list = ls())
graphics.off()

suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(here)
  library(openxlsx)
})

# Disable ANSI colors for clean log files
options(crayon.enabled = FALSE)

# Load Configuration
config_path <- here("config/global_params.yml")
if (!file.exists(config_path)) stop("Config file not found!")
config <- yaml::read_yaml(config_path)

# 2. Logging Setup
# ------------------------------------------------------------------------------
out_root <- here(config$output_root)
if (!is.null(config$project_name) && config$project_name != "") {
  out_root <- paste0(out_root, "_", config$project_name)
}
if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(out_root, paste0("pipeline_log_", timestamp, ".txt"))

# Initialize Log
cat(sprintf("=== PIPELINE STARTED: %s ===\n", Sys.time()), file = log_file)
message(sprintf("[System] Saving full log to: %s", log_file))

# Helper to log messages to file AND console
log_handler <- function(m) {
  cat(conditionMessage(m), file = log_file, append = TRUE, sep = "\n")
  # Message bubbles up to console automatically
}

warn_handler <- function(w) {
  cat(paste0("WARNING: ", conditionMessage(w)), file = log_file, append = TRUE, sep = "\n")
  # Warning bubbles up to console automatically
}

# 3. Pipeline Execution
# ------------------------------------------------------------------------------
# We wrap the entire execution in withCallingHandlers to capture output reliably
tryCatch({
  withCallingHandlers({
    
    # --- Module Loading ---
    message("\n>>> LOADING MODULES <<<")
    list.files(here("R"), pattern = "\\.R$", full.names = TRUE) %>% walk(source)
    message("[System] Modules loaded successfully.")
    
    # --- PHASE 1: DATA PROCESSING & VIZ ---
    message("\n>>> RUNNING PHASE 1: DATA PROCESSING & VIZ <<<")
    
    # Sourcing these files executes the steps immediately
    source(here("src/01_data_processing.R"), echo = FALSE)
    source(here("src/02_visualization.R"), echo = FALSE)
    
    # --- PHASE 2: STATISTICAL ANALYSIS & SCENARIOS ---
    message("\n>>> RUNNING PHASE 2: STATISTICAL ANALYSIS & SCENARIOS <<<")
    
    source(here("src/03_statistical_analysis.R"), echo = FALSE)
    
    # --- PHASE 3: NETWORK ANALYSIS & META-ANALYSIS ---
    message("\n>>> RUNNING PHASE 3: NETWORK ANALYSIS & META-ANALYSIS <<<")
    
    source(here("src/04_network_analysis.R"), echo = FALSE)
    
    # --- PHASE 4: REPORT GENERATION ---
    message("\n>>> RUNNING PHASE 4: REPORT GENERATION <<<")
    source(here("src/05_reporting.R"), echo = FALSE)

    # --- PHASE 5: TRUSTWORTHY / RESPONSIBLE-AI MODULE ---
    message("\n>>> RUNNING PHASE 5: TRUSTWORTHY AI <<<")
    source(here("src/06_trustworthy_ai.R"), echo = FALSE)

    final_msg <- sprintf("\n=== PIPELINE FINISHED SUCCESSFULLY: %s ===", Sys.time())
    message(final_msg)
    
  }, message = log_handler, warning = warn_handler)
  
}, error = function(e) {
  # Fatal Error Handler
  err_msg <- paste0("\n[FATAL ERROR] Pipeline stopped: ", e$message)
  cat(err_msg, file = log_file, append = TRUE, sep = "\n")
  stop(e$message)
  
}, finally = {
  options(crayon.enabled = TRUE)
})