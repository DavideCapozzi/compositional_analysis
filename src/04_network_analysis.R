# src/04_network_analysis.R
# ==============================================================================
# STEP 04: NETWORK ANALYSIS & META-ANALYSIS
# Description: Orchestrates Universal Baselines and iterates scenarios 
#              via workflows. Outputs binary payload.
# ==============================================================================

source("R/utils_io.R")
source("R/modules_viz.R")
source("R/modules_network.R")
source("R/workflows.R") 

message("\n=== PIPELINE STEP 4: NETWORK ANALYSIS & META-ANALYSIS ===")

config <- load_config("config/global_params.yml")
results_dir <- file.path(config$output_root, "results_analysis")

data_file <- file.path(config$output_root, "01_data_processing", "data_processed.rds")
if (!file.exists(data_file)) stop("[Fatal] Step 01 output not found.")
DATA <- readRDS(data_file)

payload03_path <- file.path(results_dir, "step03_payload.rds")
if (!file.exists(payload03_path)) stop("[Fatal] Step 03 payload not found. Run Step 03 first.")
payload03 <- readRDS(payload03_path)

payload04 <- list(scenarios = list(), meta_analysis = list())

pad_matrix <- function(small_mat, full_names, default_val = 0) {
  p <- length(full_names)
  big_mat <- matrix(default_val, nrow = p, ncol = p, dimnames = list(full_names, full_names))
  if (!is.null(small_mat) && nrow(small_mat) > 0) {
    valid_n <- rownames(small_mat)
    big_mat[valid_n, valid_n] <- small_mat
  }
  return(big_mat)
}

# 1. COMPUTE UNIVERSAL BASELINES
message("\n--- PHASE 1: GENERATING UNIVERSAL BASELINES ---")
cyto_base_dir <- file.path(results_dir, "cytoscape_baselines")
if (!dir.exists(cyto_base_dir)) dir.create(cyto_base_dir, recursive = TRUE)

macro_groups <- list()
for (scen in config$analysis_scenarios) {
  macro_groups[[scen$control_label]] <- unique(c(macro_groups[[scen$control_label]], scen$control_groups))
  macro_groups[[scen$case_label]]    <- unique(c(macro_groups[[scen$case_label]], scen$case_groups))
}

baselines <- list()
strat_col <- config$stratification$column

for (label in names(macro_groups)) {
  sub_meta <- DATA$metadata %>% dplyr::filter(.data[[strat_col]] %in% macro_groups[[label]])
  sub_mat <- DATA$hybrid_data_z[sub_meta$Patient_ID, DATA$hybrid_markers, drop=FALSE]
  
  if (nrow(sub_mat) < if(!is.null(config$stats$min_network_n)) config$stats$min_network_n else 5) next
  
  var_thresh <- if(!is.null(config$stats$min_variance)) config$stats$min_variance else 1e-6
  vars <- apply(sub_mat, 2, var, na.rm = TRUE)
  valid_cols <- names(vars)[!(vars < var_thresh | is.na(vars))]
  if (length(valid_cols) < 2) next
  
  auto_cores <- max(1, parallel::detectCores(logical = TRUE) - 1, na.rm = TRUE)
  target_cores <- if(!is.null(config$stats$n_cores) && config$stats$n_cores != "auto") as.numeric(config$stats$n_cores) else auto_cores
  
  base_obj_raw <- compute_universal_baseline(
    mat = sub_mat[, valid_cols, drop = FALSE], label = label,
    n_boot = if(!is.null(config$stats$n_boot)) config$stats$n_boot else 100, seed = if(!is.null(config$stats$seed)) config$stats$seed else 123,
    n_cores = target_cores,
    threshold_type = if(!is.null(config$stats$network_threshold_type)) config$stats$network_threshold_type else "percentile",
    threshold_value = if(!is.null(config$stats$network_edge_threshold)) config$stats$network_edge_threshold else 0.15
  )
  
  all_features <- colnames(sub_mat)
  baselines[[label]] <- list(
    label = base_obj_raw$label, mat = sub_mat, 
    pcor = pad_matrix(base_obj_raw$pcor, all_features, 0), raw_cor = pad_matrix(base_obj_raw$raw_cor, all_features, 0),
    stability = pad_matrix(base_obj_raw$stability, all_features, 0), 
    stability_freq = pad_matrix(base_obj_raw$stability_freq, all_features, 0),
    adj_final = pad_matrix(base_obj_raw$adj_final, all_features, FALSE),
    applied_threshold = base_obj_raw$applied_threshold
  )
  
  # --- CYTOSCAPE BASELINE EXPORT (Dynamic Filtering) ---
  pcor_mat <- base_obj_raw$pcor
  stab_mat <- base_obj_raw$stability_freq
  edges_idx <- which(upper.tri(pcor_mat), arr.ind = TRUE)
  
  if (nrow(edges_idx) > 0) {
    nodes <- colnames(pcor_mat)
    
    # Build complete raw edge list
    base_edges_all <- data.frame(
      Source = nodes[edges_idx[,1]],
      Target = nodes[edges_idx[,2]],
      Pcor = pcor_mat[edges_idx],
      StabFreq = stab_mat[edges_idx],
      stringsAsFactors = FALSE
    )
    
    # Apply dynamic config filters using the Rank-Based Backbone strategy
    b_conf <- config$cytoscape_export$baseline
    base_edges_res <- filter_baseline_backbone(base_edges_all, b_conf)
    
    # Process Main Network
    base_edges <- base_edges_res$main
    if (!is.null(base_edges) && nrow(base_edges) > 0) {
      readr::write_csv(base_edges, file.path(cyto_base_dir, paste0(label, "_baseline_network.csv")))
      
      adj_export <- matrix(0, nrow = length(nodes), ncol = length(nodes), dimnames = list(nodes, nodes))
      for(k in 1:nrow(base_edges)) {
        adj_export[base_edges$Source[k], base_edges$Target[k]] <- 1
        adj_export[base_edges$Target[k], base_edges$Source[k]] <- 1
      }
      
      topo_base <- calculate_node_topology(adj_export)
      if (!is.null(topo_base)) {
        topo_base <- annotate_marker_categories(topo_base, config)
        readr::write_csv(topo_base, file.path(cyto_base_dir, paste0(label, "_node_attributes.csv")))
      }
    }
    
    # Process Remaining Network
    base_edges_rem <- base_edges_res$remaining
    if (!is.null(base_edges_rem) && nrow(base_edges_rem) > 0) {
      # Dynamically create the remaining directory only if edges exist
      cyto_base_rem_dir <- file.path(cyto_base_dir, "remaining_baselines")
      if (!dir.exists(cyto_base_rem_dir)) dir.create(cyto_base_rem_dir, recursive = TRUE)
      
      readr::write_csv(base_edges_rem, file.path(cyto_base_rem_dir, paste0(label, "_baseline_network_remaining.csv")))
      
      adj_export_rem <- matrix(0, nrow = length(nodes), ncol = length(nodes), dimnames = list(nodes, nodes))
      for(k in 1:nrow(base_edges_rem)) {
        adj_export_rem[base_edges_rem$Source[k], base_edges_rem$Target[k]] <- 1
        adj_export_rem[base_edges_rem$Target[k], base_edges_rem$Source[k]] <- 1
      }
      
      topo_base_rem <- calculate_node_topology(adj_export_rem)
      if (!is.null(topo_base_rem)) {
        topo_base_rem <- annotate_marker_categories(topo_base_rem, config)
        readr::write_csv(topo_base_rem, file.path(cyto_base_rem_dir, paste0(label, "_baseline_network_remaining_node_attributes.csv")))
      }
    }
  }
}

# 2. DIFFERENTIAL OVERLAYS (Delegated to Workflows)
message("\n--- PHASE 2: SCENARIO EXECUTION ---")
diff_edges_list <- list()
scenario_colors <- c()

for (scen in config$analysis_scenarios) {
  base_ctrl <- baselines[[scen$control_label]]
  base_case <- baselines[[scen$case_label]]
  
  if (is.null(base_ctrl) || is.null(base_case)) {
    message(sprintf("      [Skip] Missing baselines for %s.", scen$id))
    next
  }
  
  spls_drivers <- if (!is.null(payload03$scenarios[[scen$id]]$drivers)) payload03$scenarios[[scen$id]]$drivers else NULL
  
  scen_res <- run_network_scenario_pipeline(
    scenario = scen, base_ctrl = base_ctrl, base_case = base_case, 
    config = config, out_dir = file.path(results_dir, scen$id), spls_drivers = spls_drivers
  )
  
  if (!is.null(scen_res)) {
    payload04$scenarios[[scen$id]] <- scen_res
    if(length(scen_res$sig_edge_ids) > 0) {
      diff_edges_list[[scen$id]] <- scen_res$sig_edge_ids
      tgt_grp <- scen$case_label
      scenario_colors[scen$id] <- if(!is.null(config$colors$groups[[tgt_grp]])) config$colors$groups[[tgt_grp]] else "grey50"
    }
  }
}

# 3. META-ANALYSIS EXPORT
message("\n--- PHASE 3: DIFFERENTIAL OVERLAP ---")
if (length(diff_edges_list) >= 2) {
  meta_dir <- file.path(results_dir, "Meta_Analysis")
  if (!dir.exists(meta_dir)) dir.create(meta_dir, recursive = TRUE)
  
  # --- Categorical Edge Signature Heatmap ---
  # Define the specific scenarios to track chronologically/phenotypically
  target_heat_scenarios <- c("00_HDvsTumor", "01_HDvsEP", "02_HDvsLS", "03_EPvsLS")
  
  message("   [Meta-Analysis] Generating Categorical Edge Signature Heatmap...")
  pdf(file.path(meta_dir, "Edge_Signature_Trajectories.pdf"), width = 10, height = 12)
  tryCatch({
    hm_obj <- viz_plot_categorical_edge_heatmap(
      scenarios_list = payload04$scenarios,
      target_scenarios = target_heat_scenarios,
      config = config,
      title = "Phenotypic Network Rewiring"
    )
    if (!is.null(hm_obj)) ComplexHeatmap::draw(hm_obj, merge_legend = TRUE)
  }, error = function(e) { 
    p_err <- ggplot2::ggplot() + 
      ggplot2::annotate("text", x = 0, y = 0, label = paste("Heatmap Plot Error:", e$message), color = "darkred", size = 5, fontface = "bold") + 
      ggplot2::theme_void()
    print(p_err)
  })
  dev.off()
  
  # --- Legacy: Standard UpSet/Venn Diagram ---
  message("   [Meta-Analysis] Generating UpSet/Venn Overlap...")
  pdf(file.path(meta_dir, "Differential_Networks_Overlap.pdf"), width = 8, height = 8)
  tryCatch({
    viz_plot_differential_overlap(edge_list = diff_edges_list, fill_colors = scenario_colors, title = "Differential Edges Overlap")
  }, error = function(e) { 
    p_err <- ggplot2::ggplot() + 
      ggplot2::annotate("text", x = 0, y = 0, label = paste("Overlap Plot Error:", e$message), color = "darkred", size = 5, fontface = "bold") + 
      ggplot2::theme_void()
    print(p_err)
  })
  dev.off()
  
  payload04$meta_analysis$diff_edges_list <- diff_edges_list
  payload04$meta_analysis$scenario_colors <- scenario_colors
}

payload04_path <- file.path(results_dir, "step04_payload.rds")
saveRDS(payload04, payload04_path)
message(sprintf("   [Output] Step 04 Payload saved: %s", payload04_path))
message("\n=== STEP 4 COMPLETE ===\n")