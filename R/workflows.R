# R/workflows.R
# ==============================================================================
# WORKFLOWS & MACROS MODULE
# Description: Orchestration functions that combine multiple analytical steps 
#              (Filtering, Hypothesis Testing, Multivariate, Viz).
# Dependencies: dplyr, modules_hypothesis, modules_multivariate, modules_viz
# ==============================================================================

library(dplyr)

#' @title Annotate Node Topology with Marker Categories
#' @description Safely maps config-defined marker categories to node attributes.
#'              Handles empty lists, prevents row duplication, and sanitizes types.
#' @param node_df Dataframe containing at least a 'Node' column.
#' @param config Global configuration list.
#' @return Dataframe with an appended 'Marker_Category' column.
annotate_marker_categories <- function(node_df, config) {
  if (is.null(node_df) || !("Node" %in% names(node_df))) return(node_df)
  
  # Protect against pre-existing column conflicts
  if ("Marker_Category" %in% names(node_df)) {
    node_df$Marker_Category <- NULL
  }
  
  if (is.null(config$qc_reporting$marker_categories)) {
    node_df$Marker_Category <- "Other"
    return(node_df)
  }
  
  cat_list <- config$qc_reporting$marker_categories
  
  # Safely construct mapping dataframe bypassing empty lists
  cat_df_list <- lapply(names(cat_list), function(cat_name) {
    nodes_in_cat <- unlist(cat_list[[cat_name]])
    if (length(nodes_in_cat) > 0) {
      data.frame(Node = as.character(nodes_in_cat), 
                 Marker_Category = as.character(cat_name), 
                 stringsAsFactors = FALSE)
    } else {
      NULL
    }
  })
  
  cat_df <- do.call(rbind, Filter(Negate(is.null), cat_df_list))
  
  if (!is.null(cat_df) && nrow(cat_df) > 0) {
    # Resolve conflicting manual assignments in config by keeping the first match
    cat_df <- cat_df[!duplicated(cat_df$Node), ]
    node_df <- dplyr::left_join(node_df, cat_df, by = "Node")
  } else {
    node_df$Marker_Category <- NA_character_
  }
  
  # Impute uncategorized markers
  node_df$Marker_Category[is.na(node_df$Marker_Category)] <- "Other"
  
  return(node_df)
}

#' @title Execute Complete Statistical Pipeline for a Single Scenario
#' @description 
#' Filters data, handles zero-variance, checks dispersion assumptions, 
#' runs PERMANOVA, fits sPLS-DA, and generates plots.
#' 
#' @param scenario List containing scenario definitions from config.
#' @param full_mat Global numeric matrix (Z-scored).
#' @param full_meta Global metadata dataframe.
#' @param config The complete configuration object.
#' @param out_dir Directory path to save scenario-specific PDF plots.
#' @return A list of dataframes: 'dispersion', 'permanova', and 'drivers'.
run_scenario_pipeline <- function(scenario, full_mat, full_meta, config, out_dir) {
  
  # 1. Subsetting Data
  target_groups <- unlist(c(scenario$case_groups, scenario$control_groups))
  
  sub_meta <- full_meta %>% 
    dplyr::filter(Group %in% target_groups) %>% 
    dplyr::mutate(Analysis_Group = ifelse(Group %in% unlist(scenario$case_groups), 
                                          scenario$case_label, 
                                          scenario$control_label)) %>% 
    dplyr::mutate(Analysis_Group = factor(Analysis_Group, 
                                          levels = c(scenario$control_label, scenario$case_label)))
  
  sub_mat <- full_mat[sub_meta$Patient_ID, , drop = FALSE]
  
  if (nrow(sub_mat) == 0) {
    message("      [Skip] No samples found after filtering.")
    return(NULL)
  }
  
  counts <- table(sub_meta$Analysis_Group)
  message(sprintf("      [Info] Samples: %s=%d vs %s=%d", 
                  names(counts)[1], counts[1], names(counts)[2], counts[2]))
  
  # 2. Zero-Variance Filtering (Universal Space Protection)
  var_thresh <- if(!is.null(config$stats$min_variance)) config$stats$min_variance else 1e-6
  vars <- apply(sub_mat, 2, var, na.rm = TRUE)
  low_vars <- vars < var_thresh | is.na(vars)
  
  if (sum(low_vars) > 0) {
    dropped_names <- names(vars)[low_vars]
    message(sprintf("      [Prep] Dropping %d flat features with zero pooled variance: %s", 
                    sum(low_vars), paste(dropped_names, collapse = ", ")))
    sub_mat <- sub_mat[, !low_vars, drop = FALSE]
  }
  
  if (ncol(sub_mat) < 2) {
    message("      [Skip] Less than 2 features remaining after variance filter.")
    return(NULL)
  }
  
  sub_data_input <- cbind(sub_meta, as.data.frame(sub_mat))
  res_list <- list(dispersion = NULL, permanova = NULL, drivers = NULL)
  
  # 3. Beta-Dispersion Check (Strictly Scenario-Specific)
  tryCatch({
    disp_obj <- test_coda_dispersion(
      data_input = sub_data_input, 
      group_col = "Analysis_Group", 
      metadata_cols = colnames(sub_meta),
      n_perm = if(!is.null(config$stats$n_perm)) config$stats$n_perm else 999
    )
    res_list$dispersion <- as.data.frame(disp_obj$anova_table)
  }, error = function(e) message(sprintf("      [Fail] Dispersion test failed: %s", e$message)))
  
  # 4. PERMANOVA
  tryCatch({
    perm_obj <- test_coda_permanova(
      data_input = sub_data_input, 
      group_col = "Analysis_Group", 
      metadata_cols = colnames(sub_meta),
      n_perm = if(!is.null(config$stats$n_perm)) config$stats$n_perm else 999
    )
    res_list$permanova <- as.data.frame(perm_obj)
  }, error = function(e) message(sprintf("      [Fail] PERMANOVA failed: %s", e$message)))
  
  # 5. sPLS-DA & Visualization
  valid_levels <- names(counts)[counts > 0]
  
  # Enforce minimum 3 samples per class to safely execute CV without crashing mixOmics
  if (length(valid_levels) >= 2 && min(counts[valid_levels]) >= 3) {
    tryCatch({
      n_comp <- if(!is.null(config$multivariate$n_comp)) config$multivariate$n_comp else 2
      cv_folds <- if(!is.null(config$multivariate$validation_folds)) config$multivariate$validation_folds else 5
      actual_folds <- min(cv_folds, min(counts[valid_levels]))
      
      # Secondary safety constraint for fold assignment
      if (actual_folds < 3) actual_folds <- 3
      
      spls_res <- run_splsda_model(
        data_z = sub_mat,
        metadata = sub_meta, 
        group_col = "Analysis_Group", 
        n_comp = n_comp, 
        folds = actual_folds, 
        n_repeat = if(!is.null(config$multivariate$n_repeat_cv)) config$multivariate$n_repeat_cv else 10
      )
      
      drivers_df <- extract_plsda_loadings(spls_res)
      
      if (nrow(drivers_df) > 0) {
        # Setup specific palette for plotting
        scen_palette <- c()
        scen_palette[scenario$control_label] <- config$colors$control
        
        # Safe extraction for case color to avoid out of bounds
        case_col <- "firebrick"
        if(!is.null(config$colors$groups[[scenario$case_label]])) {
          case_col <- config$colors$groups[[scenario$case_label]]
        } else if (length(config$colors$cases) > 0) {
          case_col <- config$colors$cases[[1]]
        }
        scen_palette[scenario$case_label] <- case_col
        
        viz_report_plsda(
          pls_res = spls_res, 
          drivers_df = drivers_df, 
          metadata_viz = sub_meta, 
          colors_viz = scen_palette, 
          out_path = file.path(out_dir, paste0(scenario$id, "_Plot_sPLSDA.pdf")),
          group_col = "Analysis_Group"
        )
        
        # Dynamic renaming for Excel readability
        contrast_lbl <- paste0(scenario$case_label, "_VS_", scenario$control_label)
        drivers_df <- drivers_df %>% 
          dplyr::rename_with(.fn = ~ paste0("Weight_", contrast_lbl, "_PC1"), .cols = dplyr::matches("^Comp1_Weight$")) %>%
          dplyr::rename_with(.fn = ~ paste0("Weight_", contrast_lbl, "_PC2"), .cols = dplyr::matches("^Comp2_Weight$"))
        
        res_list$drivers <- drivers_df
      }
    }, error = function(e) message(sprintf("      [Fail] sPLS-DA failed: %s", e$message)))
  } else {
    message("      [Skip] sPLS-DA skipped: Class size too small for reliable CV (<3).")
  }
  
  return(res_list)
}

#' @title Execute Network Pipeline for a Single Scenario
#' @description 
#' Orchestrates differential network computation, topology calculation, 
#' visualization (Hub-Driver, Density, Networks), and Cytoscape export.
#' Now automatically exports remaining networks to dedicated subdirectories 
#' if rank limits are applied.
#' 
#' @param scenario List containing scenario definitions.
#' @param base_ctrl Universal baseline object for control group.
#' @param base_case Universal baseline object for case group.
#' @param config Global configuration list.
#' @param out_dir Directory to save scenario-specific artifacts.
#' @param spls_drivers Dataframe of sPLS-DA drivers from Step 03 payload (optional).
#' @return A list containing edge tables and significant edge IDs.
run_network_scenario_pipeline <- function(scenario, base_ctrl, base_case, config, out_dir, spls_drivers = NULL) {
  
  message(sprintf("\n   [Scenario] Processing Network: %s...", scenario$id))
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  auto_cores <- max(1, parallel::detectCores(logical = TRUE) - 1, na.rm = TRUE)
  target_cores <- if(!is.null(config$stats$n_cores) && config$stats$n_cores != "auto") as.numeric(config$stats$n_cores) else auto_cores
  
  # 1. Differential Overlay
  net_res <- compute_differential_overlay(
    base_ctrl = base_ctrl,
    base_case = base_case,
    n_perm = if(!is.null(config$stats$n_perm)) config$stats$n_perm else 1000,
    seed = if(!is.null(config$stats$seed)) config$stats$seed else 123,
    n_cores = target_cores,
    pvalue_thresh = if(!is.null(config$stats$pvalue_threshold)) config$stats$pvalue_threshold else 0.05,
    min_marginal_cor = if(!is.null(config$stats$min_marginal_cor)) config$stats$min_marginal_cor else 0.10
  )
  
  if (is.null(net_res) || nrow(net_res$edges_table) == 0) {
    message("      [Skip] No differential edges computed.")
    return(NULL)
  }
  
  saveRDS(net_res, file.path(out_dir, paste0(scenario$id, "_Differential_Network.rds")))
  
  # 2. Topology Calculation
  topo_ctrl <- if (sum(net_res$adj_final$ctrl)>0) calculate_node_topology(net_res$adj_final$ctrl) else NULL
  topo_case <- if (sum(net_res$adj_final$case)>0) calculate_node_topology(net_res$adj_final$case) else NULL
  
  # 3. Visualizations
  if (!is.null(spls_drivers) && nrow(spls_drivers) > 0 && (!is.null(topo_ctrl) || !is.null(topo_case))) {
    pdf(file.path(out_dir, paste0(scenario$id, "_Hub_Driver_Report.pdf")), width = 11, height = 8)
    if(!is.null(topo_ctrl)) print(plot_hub_driver_quadrant(integrate_hub_drivers(spls_drivers, topo_ctrl, "Degree"), "Degree", paste("\nNetwork:", scenario$control_label)))
    if(!is.null(topo_case)) print(plot_hub_driver_quadrant(integrate_hub_drivers(spls_drivers, topo_case, "Degree"), "Degree", paste("\nNetwork:", scenario$case_label)))
    if(!is.null(topo_ctrl)) print(plot_hub_driver_quadrant(integrate_hub_drivers(spls_drivers, topo_ctrl, "Betweenness"), "Betweenness", paste("\nNetwork:", scenario$control_label)))
    if(!is.null(topo_case)) print(plot_hub_driver_quadrant(integrate_hub_drivers(spls_drivers, topo_case, "Betweenness"), "Betweenness", paste("\nNetwork:", scenario$case_label)))
    dev.off()
  }
  
  thresh_to_plot <- net_res$applied_threshold
  pdf(file.path(out_dir, paste0(scenario$id, "_Edge_Distribution.pdf")), width = 8, height = 6)
  if(sum(abs(net_res$networks$ctrl) > 0) > 0) {
    print(viz_plot_edge_density(net_res$networks$ctrl, adj_mat = net_res$adj_final$ctrl, threshold = thresh_to_plot, group_label = scenario$control_label))
  }
  if(sum(abs(net_res$networks$case) > 0) > 0) {
    print(viz_plot_edge_density(net_res$networks$case, adj_mat = net_res$adj_final$case, threshold = thresh_to_plot, group_label = scenario$case_label))
  }
  dev.off()
  
  pdf(file.path(out_dir, paste0(scenario$id, "_Plot_Networks.pdf")), width = 12, height = 6)
  if(!is.null(topo_ctrl)) print(plot_network_structure(net_res$adj_final$ctrl, net_res$networks$ctrl, title = paste("Control:", scenario$control_label), min_cor = 0))
  if(!is.null(topo_case)) print(plot_network_structure(net_res$adj_final$case, net_res$networks$case, title = paste("Case:", scenario$case_label), min_cor = 0))
  dev.off()
  
  # 4. CYTOSCAPE EXPORT (Dynamic Filtering via Config)
  cyto_dir <- file.path(out_dir, "cytoscape_export")
  if (!dir.exists(cyto_dir)) dir.create(cyto_dir)
  
  d_conf <- config$cytoscape_export$differential
  diff_edges_res <- filter_differential_export(net_res$edges_table, d_conf)
  
  # Helper to process and export networks dynamically to designated paths
  export_cyto_diff <- function(edges_to_export, file_suffix, target_dir) {
    if (!is.null(edges_to_export) && nrow(edges_to_export) > 0) {
      # Ensure target directory exists before writing
      if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
      
      diff_net <- edges_to_export %>%
        dplyr::select(
          Source = Node1, Target = Node2, Weight = Diff_Score, Significant, P_Value, FDR, Edge_Category,
          dplyr::starts_with("Pcor_"), dplyr::starts_with("Spearman_"),
          dplyr::starts_with("Mech_"), dplyr::starts_with("Is_Stable_"), dplyr::starts_with("StabFreq_")
        ) %>%
        dplyr::rename(Interaction = Edge_Category)
      
      readr::write_csv(diff_net, file.path(target_dir, paste0(scenario$id, file_suffix, ".csv")))
      
      diff_nodes <- unique(c(edges_to_export$Node1, edges_to_export$Node2))
      adj_diff <- matrix(0, nrow = length(diff_nodes), ncol = length(diff_nodes), dimnames = list(diff_nodes, diff_nodes))
      for(k in 1:nrow(edges_to_export)) {
        adj_diff[edges_to_export$Node1[k], edges_to_export$Node2[k]] <- 1
        adj_diff[edges_to_export$Node2[k], edges_to_export$Node1[k]] <- 1
      }
      
      topo_diff <- calculate_node_topology(adj_diff)
      if (!is.null(topo_diff)) {
        all_nodes <- topo_diff %>% dplyr::select(Node, Degree, Betweenness, Closeness, Eigen_Centrality)
        all_nodes <- annotate_marker_categories(all_nodes, config)
        
        if (!is.null(spls_drivers) && nrow(spls_drivers) > 0) {
          drv_summ <- spls_drivers %>% dplyr::group_by(Marker) %>% dplyr::slice_max(order_by = Importance, n = 1, with_ties = FALSE) %>% dplyr::ungroup()
          all_nodes <- dplyr::left_join(all_nodes, drv_summ[, c("Marker", "Importance", "Direction")], by = c("Node" = "Marker"))
        }
        num_cols <- names(all_nodes)[sapply(all_nodes, is.numeric)]
        all_nodes[num_cols] <- lapply(all_nodes[num_cols], function(x) replace(x, is.na(x), 0))
        readr::write_csv(all_nodes, file.path(target_dir, paste0(scenario$id, file_suffix, "_node_attributes.csv")))
      }
    }
  }
  
  # Define sub-directory for remaining edges
  cyto_rem_dir <- file.path(cyto_dir, "remaining_differential")
  
  # Execute helper for both networks targeting correct paths
  export_cyto_diff(diff_edges_res$main, "_diff_network", cyto_dir)
  export_cyto_diff(diff_edges_res$remaining, "_diff_network_remaining", cyto_rem_dir)
  
  # 5. Prepare Return Object (Passed to Meta-Analysis)
  sig_edges <- net_res$edges_table %>% dplyr::filter(Edge_Category != "Weak" & Significant == TRUE) %>% dplyr::arrange(P_Value)
  edge_ids <- if(nrow(sig_edges) > 0) apply(sig_edges, 1, function(r) paste(sort(c(r["Node1"], r["Node2"])), collapse="~")) else character(0)
  
  rewiring_df <- NULL
  if (!is.null(topo_ctrl) && !is.null(topo_case)) {
    rewiring_df <- calculate_jaccard_rewiring(net_res$adj_final$ctrl, net_res$adj_final$case)
  }
  
  return(list(
    edges_table = net_res$edges_table,
    sig_edges = sig_edges,
    sig_edge_ids = unique(edge_ids),
    topo_ctrl = topo_ctrl,
    topo_case = topo_case,
    rewiring = rewiring_df
  ))
}