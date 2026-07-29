# src/03_statistical_analysis.R
# ==============================================================================
# STEP 03: STATISTICAL ANALYSIS & SCENARIO ORCHESTRATION
# Description: Global Hypothesis testing (PERMANOVA, sPLS-DA) and 
#              Iterative Multi-Scenario execution. Outputs a binary payload.
# ==============================================================================

source("R/utils_io.R")          
source("R/modules_hypothesis.R") 
source("R/modules_multivariate.R")
source("R/modules_viz.R") 
source("R/workflows.R")

message("\n=== PIPELINE STEP 3: STATISTICAL ANALYSIS & SCENARIOS ===")

# 1. INITIALIZATION & DATA SETUP
config <- load_config("config/global_params.yml")
input_file <- file.path(config$output_root, "01_data_processing", "data_processed.rds")

if (!file.exists(input_file)) stop("[Fatal] Step 01 output not found. Run src/01_data_processing.R first.")
DATA <- readRDS(input_file)

results_dir <- file.path(config$output_root, "results_analysis")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

global_dir <- file.path(results_dir, "Global_Stats")
if (!dir.exists(global_dir)) dir.create(global_dir, recursive = TRUE)

df_global <- DATA$hybrid_data_z
safe_markers <- DATA$hybrid_markers
meta_viz <- DATA$metadata 

strat_col <- config$stratification$column
if (!strat_col %in% names(df_global)) {
  stop(sprintf("[Fatal] Column '%s' not found in processed data.", strat_col))
}

df_global$Group <- df_global[[strat_col]]
meta_viz$Group  <- meta_viz[[strat_col]]

# Initialize Payload
payload03 <- list(global = list(), scenarios = list())

# 2. GLOBAL STATISTICAL ANALYSIS
message("\n--- RUNNING GLOBAL STATISTICS ---")

target_control <- as.vector(config$control_group)
target_cases <- as.vector(config$case_groups)
target_all <- unique(c(target_control, target_cases))

df_stats_global <- df_global %>% filter(Group %in% target_all)
meta_stats_global <- meta_viz %>% filter(Patient_ID %in% df_stats_global$Patient_ID)
df_stats_global$Group <- factor(df_stats_global$Group, levels = target_all)
meta_stats_global$Group <- factor(meta_stats_global$Group, levels = target_all)

clean_targets <- unique(as.character(na.omit(target_all)))
current_palette <- get_palette(config, match_groups = clean_targets)[clean_targets]
if (any(is.na(current_palette))) current_palette[is.na(current_palette)] <- "grey50"

# 2.1 Global PERMANOVA
# Guarded like every other analysis in this step: test_coda_permanova() hard-stops on NA
# input and adonis2 errors on a single-level factor, either of which would otherwise abort
# steps 04/05/06 via the re-stop() in main.R.
message("   [Stats] Running Global PERMANOVA (Multiclass)...")
tryCatch({
  perm_global <- test_coda_permanova(
    data_input = df_stats_global[, c("Group", safe_markers)],
    group_col = "Group", n_perm = config$stats$n_perm
  )
  payload03$global$permanova <- as.data.frame(perm_global)
}, error = function(e) {
  message(sprintf("      [WARN] Global PERMANOVA failed: %s", conditionMessage(e)))
})

# 2.2 Global Dispersion
message("   [Stats] Running Global Dispersion Check...")
tryCatch({
  disp_global <- test_coda_dispersion(
    data_input = df_stats_global[, c("Group", safe_markers)], 
    group_col = "Group", n_perm = config$stats$n_perm
  )
  payload03$global$dispersion <- disp_global$anova_table
  
  disp_pairwise <- run_pairwise_betadisper(
    data_input = df_stats_global[, c("Group", safe_markers)], 
    group_col = "Group", metadata_cols = c("Patient_ID"),
    n_perm = config$stats$n_perm, min_n = if(!is.null(config$stats$min_sample_size)) config$stats$min_sample_size else 4
  )
  if (!is.null(disp_pairwise) && nrow(disp_pairwise) > 0) {
    payload03$global$pairwise_dispersion <- disp_pairwise
  }
}, error = function(e) message(paste("   [ERROR] Beta-Dispersion failed:", e$message)))

# 2.3 Global sPLS-DA
if (config$multivariate$run_plsda) {
  message("   [sPLS-DA] Fitting STRATIFIED model (Multiclass)...")
  set.seed(config$stats$seed) 
  tryCatch({
    pls_res <- run_splsda_model(
      data_z = df_stats_global[, safe_markers], metadata = meta_stats_global, 
      group_col = "Group", n_comp = config$multivariate$n_comp,
      folds = config$multivariate$validation_folds, n_repeat = if(!is.null(config$multivariate$n_repeat_cv)) config$multivariate$n_repeat_cv else 50
    ) 
    top_drivers <- extract_plsda_loadings(pls_res)
    payload03$global$splsda_drivers <- top_drivers
    payload03$global$splsda_perf <- extract_plsda_performance(pls_res)
    
    viz_report_plsda(
      pls_res = pls_res, drivers_df = top_drivers, metadata_viz = meta_stats_global, 
      colors_viz = current_palette, out_path = file.path(global_dir, "Global_sPLSDA_Results.pdf"),
      group_col = "Group", n_top_boxplots = if(!is.null(config$multivariate$top_n_loadings)) config$multivariate$top_n_loadings else 9
    )
  }, error = function(e) message(paste("   [ERROR] Stratified sPLS-DA Failed:", e$message)))
}

# 3. SCENARIO ORCHESTRATION
message("\n--- RUNNING SCENARIO-SPECIFIC STATISTICS ---")

if (!is.null(config$analysis_scenarios) && length(config$analysis_scenarios) > 0) {
  full_meta <- DATA$metadata
  if (strat_col %in% colnames(full_meta)) full_meta$Group <- full_meta[[strat_col]]
  
  full_mat <- DATA$hybrid_data_z %>% dplyr::select(dplyr::all_of(safe_markers)) %>% as.matrix()
  rownames(full_mat) <- DATA$metadata$Patient_ID
  
  for (scenario in config$analysis_scenarios) {
    message(sprintf("\n   [Scenario] Executing: %s", scenario$id))
    scen_dir <- file.path(results_dir, scenario$id)
    if (!dir.exists(scen_dir)) dir.create(scen_dir, recursive = TRUE)
    
    tryCatch({
      scen_res <- run_scenario_pipeline(
        scenario = scenario, full_mat = full_mat, full_meta = full_meta, config = config, out_dir = scen_dir
      )
      if (!is.null(scen_res)) {
        payload03$scenarios[[scenario$id]] <- list(
          dispersion = scen_res$dispersion, permanova = scen_res$permanova, drivers = scen_res$drivers
        )
      }
    }, error = function(e) {
      message(sprintf("!!! [Error] Scenario '%s' failed: %s", scenario$id, e$message))
    })
  }
}

payload_path <- file.path(results_dir, "step03_payload.rds")
saveRDS(payload03, payload_path)
message(sprintf("\n   [Output] Step 03 Payload saved: %s", payload_path))

message("=== STEP 3 COMPLETE ===\n")