# R/modules_qc.R
# ==============================================================================
# QUALITY CONTROL MODULE
# Description: Functions for data filtering (Variance, Missingness, Outliers).
# ==============================================================================

#' @title Detect Multivariate Outliers (PCA-based Mahalanobis)
#' @description 
#' Identifies outliers within specific groups using Robust Mahalanobis distance 
#' calculated on the first Principal Components (PCs).
#' 
#' @param mat Numeric matrix (Samples x Markers).
#' @param groups Vector of group labels corresponding to rows of mat.
#' @param conf_level Confidence level for Chi-squared cutoff (default 0.99).
#' @return A logical vector (TRUE = Outlier, FALSE = Keep).
detect_pca_outliers <- function(mat, groups, conf_level = 0.99) {
  
  is_outlier <- rep(FALSE, nrow(mat))
  
  # Ensure groups are factors
  groups <- as.factor(groups)
  
  for (g in levels(groups)) {
    # Get indices for current group
    idx <- which(groups == g)
    
    # Skip if too few samples (need at least 5 for PCA/Covariance stability)
    if (length(idx) < 5) {
      if (length(idx) > 0) {
        message(sprintf("   [QC] Bypassing PCA outlier detection for group '%s' (n=%d < 5).", g, length(idx)))
      }
      next
    } # <--- FIXED: Added missing closing brace here!
    
    sub_mat <- mat[idx, , drop = FALSE]
    
    # 1. Handle NAs: First remove columns that are 100% NA in this subgroup
    # (Median imputation would fail for these, returning NA)
    na_counts <- colSums(is.na(sub_mat))
    valid_cols <- na_counts < nrow(sub_mat)
    sub_mat <- sub_mat[, valid_cols, drop = FALSE]
    
    if (ncol(sub_mat) < 2) next
    
    # 2. Impute remaining sporadic NAs (Median Imputation)
    if (any(is.na(sub_mat))) {
      sub_mat <- apply(sub_mat, 2, function(x) {
        if (all(is.na(x))) return(x) # Should be handled by step 1, but safety first
        x[is.na(x)] <- median(x, na.rm = TRUE)
        return(x)
      })
    }
    
    # 3. Remove zero-variance columns locally (Robust check)
    # Use na.rm = TRUE to prevent NA propagation
    vars <- apply(sub_mat, 2, var, na.rm = TRUE)
    # Check for non-NA and strictly positive variance (use tolerance for float precision)
    keep_vars <- !is.na(vars) & vars > 1e-12
    
    sub_mat <- sub_mat[, keep_vars, drop = FALSE]
    
    if (ncol(sub_mat) < 2) next 
    
    # Run PCA (Scale is important)
    # We catch errors in case of singular matrices
    tryCatch({
      pca_res <- prcomp(sub_mat, scale. = TRUE, center = TRUE)
      
      # Use top PCs explaining ~90% variance or max 5 PCs to avoid noise
      # If N is small, use fewer PCs (N-2)
      max_available_pcs <- ncol(pca_res$x)
      n_pc <- min(ncol(sub_mat), length(idx) - 2, 5) 
      
      # Safety check: ensure n_pc does not exceed available PCs or drop below 1
      if (n_pc > max_available_pcs) n_pc <- max_available_pcs
      if (n_pc < 1) n_pc <- 1
      
      scores <- pca_res$x[, 1:n_pc, drop = FALSE]
      
      # Calculate Robust Mahalanobis Distance
      # We use the standard mahalanobis on PCA scores (orthogonal)
      center_vec <- colMeans(scores)
      cov_mat <- cov(scores)
      
      d2 <- mahalanobis(scores, center_vec, cov_mat)
      
      # Cutoff based on Chi-Squared distribution
      cutoff <- qchisq(conf_level, df = n_pc)
      
      # Flag local outliers
      local_outliers <- d2 > cutoff
      is_outlier[idx[local_outliers]] <- TRUE
      
    }, error = function(e) {
      warning(sprintf("[QC] Outlier detection failed for group '%s': %s", g, e$message))
    })
  }
  
  return(is_outlier)
}

#' @title Row Missingness on a Per-Cohort Panel Basis
#' @description
#' Computes each row's missingness fraction over only those markers that the row's own
#' cohort actually measured. Falls back to the plain union-of-all-markers rate when no
#' panel map is supplied, or for any row whose cohort is not present in the map.
#' @param mat Numeric matrix (Samples x Markers).
#' @param metadata Dataframe aligned to mat, carrying the cohort in a "Group" column.
#' @param panel_map Named list of cohort -> marker names, or NULL.
#' @return Numeric vector of missingness fractions, one per row of mat.
.row_na_own_panel <- function(mat, metadata = NULL, panel_map = NULL) {

  if (is.null(panel_map) || is.null(metadata) || !"Group" %in% colnames(metadata)) {
    return(rowMeans(is.na(mat)))
  }

  cohorts <- as.character(metadata[["Group"]])
  out <- numeric(nrow(mat))

  for (i in seq_len(nrow(mat))) {
    own <- panel_map[[cohorts[i]]]
    cols <- if (is.null(own)) colnames(mat) else intersect(colnames(mat), own)
    # A cohort with no overlap at all would give 0/0; treat it as fully missing.
    out[i] <- if (length(cols) == 0) 1 else mean(is.na(mat[i, cols]))
  }

  return(out)
}

#' @title Run Quality Control Pipeline
#' @description
#' Filters the data matrix based on variance, missingness, and multivariate outliers.
#' Logs dropped items and updates the summary report with group breakdowns.
#' 
#' @param mat_raw Raw numeric matrix (Samples x Markers).
#' @param metadata Dataframe containing Group information (must align with mat_raw).
#' @param qc_config List containing thresholds.
#' @param stratification_col String. Column name for detailed group breakdown (e.g., Original_Source).
#' @param dropped_markers_apriori Dataframe of markers dropped by config.
#' @param dropped_samples_apriori Dataframe of samples dropped by config.
#' @param panel_map Named list of cohort -> marker names present in that cohort's sheet.
#'   When supplied, row missingness is computed over each patient's OWN panel, so a cohort
#'   is not penalised for markers its assay never measured. NULL reproduces the legacy
#'   behaviour (missingness over the union of all panels).
#' @return A list containing the filtered 'data' and the 'report'.
run_qc_pipeline <- function(mat_raw, metadata, qc_config,
                            stratification_col = "Group",
                            dropped_markers_apriori = NULL,
                            dropped_samples_apriori = NULL,
                            panel_map = NULL) {
  
  message("[QC] Running Quality Control...")
  
  # Helper to calculate group counts
  get_group_counts <- function(meta, col_name) {
    if (is.null(col_name) || !col_name %in% colnames(meta)) return(NULL)
    return(table(meta[[col_name]]))
  }
  
  # 0. Initial Statistics
  init_counts <- get_group_counts(metadata, stratification_col)
  
  # Initialize QC Summary Object
  qc_summary <- list(
    n_row_init = nrow(mat_raw),
    n_col_init = ncol(mat_raw),
    breakdown_init = init_counts, 
    
    dropped_markers_apriori = dropped_markers_apriori,
    dropped_samples_apriori = dropped_samples_apriori, 
    
    n_col_zerovar = 0,
    dropped_rows_detail = data.frame(Patient_ID = character(), NA_Percent = numeric(),
                                     Reason = character(), Original_Source = character(),
                                     stringsAsFactors = FALSE),
    dropped_cols_detail = data.frame(Marker = character(), NA_Percent = numeric()),
    dropped_markers_panel = data.frame(Marker = character(), Absent_From = character(),
                                       stringsAsFactors = FALSE),
    n_col_panel_mismatch = 0,
    n_row_dropped = 0,
    n_col_dropped = 0
  )
  
  curr_mat <- mat_raw
  curr_meta <- metadata
  
  # 1. Remove Constant Columns (Zero Variance)
  col_vars <- apply(curr_mat, 2, var, na.rm = TRUE)
  const_cols <- which(col_vars == 0 | is.na(col_vars))
  
  if (length(const_cols) > 0) {
    message(sprintf("   [QC] Dropping %d markers with zero variance.", length(const_cols)))
    qc_summary$n_col_zerovar <- length(const_cols)
    curr_mat <- curr_mat[, -const_cols, drop = FALSE]
  }
  
  # 1b. Drop markers not measured by EVERY cohort (panel mismatch)
  # A marker absent from one cohort's panel sits well below max_na_col_pct once that cohort
  # is a minority of the sample set, so it survives column filtering and is then BPCA-imputed
  # -- fabricating values for an assay that was never run. Dropping it is the honest option.
  if (!is.null(panel_map) && length(panel_map) > 1 &&
      isTRUE(qc_config$require_marker_in_all_panels)) {
    shared_markers <- Reduce(intersect, panel_map)
    panel_mismatch <- setdiff(colnames(curr_mat), shared_markers)

    if (length(panel_mismatch) > 0) {
      message(sprintf("   [QC] Dropping %d markers absent from at least one cohort panel.",
                      length(panel_mismatch)))
      absent_from <- vapply(panel_mismatch, function(mk) {
        paste(names(panel_map)[!vapply(panel_map, function(p) mk %in% p, logical(1))],
              collapse = ", ")
      }, character(1))

      qc_summary$dropped_markers_panel <- data.frame(
        Marker = panel_mismatch,
        Absent_From = unname(absent_from),
        stringsAsFactors = FALSE
      )
      qc_summary$n_col_panel_mismatch <- length(panel_mismatch)
      curr_mat <- curr_mat[, shared_markers[shared_markers %in% colnames(curr_mat)], drop = FALSE]
    }
  }

  # 2. Filter by Missingness (Rows - Patients)
  # Count only SPORADIC missingness: markers that this patient's own cohort panel actually
  # measured. Counting structural gaps here penalises narrow-panel cohorts and, when a new
  # wide-panel cohort is added, can silently wipe out an entire existing group.
  row_na_freq <- .row_na_own_panel(curr_mat, metadata = curr_meta, panel_map = panel_map)
  drop_row_na <- which(row_na_freq > qc_config$max_na_row_pct)

  if (length(drop_row_na) > 0) {
    pids <- rownames(curr_mat)[drop_row_na]
    message(sprintf("   [QC] Dropping %d patients with >%.0f%% missingness (own-panel basis).",
                    length(pids), qc_config$max_na_row_pct * 100))
    
    # Capture Original_Source info for report (Protected against NULL)
    group_info <- if(!is.null(stratification_col) && stratification_col %in% colnames(curr_meta)) {
      as.character(curr_meta[[stratification_col]][drop_row_na])
    } else { rep("N/A", length(pids)) }
    
    qc_summary$dropped_rows_detail <- rbind(qc_summary$dropped_rows_detail, data.frame(
      Patient_ID = pids,
      NA_Percent = round(row_na_freq[drop_row_na] * 100, 2),
      Reason = "High Missingness",
      Original_Source = group_info,
      stringsAsFactors = FALSE
    ))
    
    curr_mat <- curr_mat[-drop_row_na, , drop = FALSE]
    curr_meta <- curr_meta[-drop_row_na, , drop = FALSE]
  }
  
  # 3. Filter by Missingness (Cols - Markers)
  col_na_freq <- colMeans(is.na(curr_mat))
  drop_col_idx <- which(col_na_freq > qc_config$max_na_col_pct)
  
  if (length(drop_col_idx) > 0) {
    mks <- colnames(curr_mat)[drop_col_idx]
    message(sprintf("   [QC] Dropping %d markers with >%.0f%% missingness.", 
                    length(mks), qc_config$max_na_col_pct * 100))
    
    qc_summary$dropped_cols_detail <- data.frame(
      Marker = mks,
      NA_Percent = round(col_na_freq[drop_col_idx] * 100, 2)
    )
    qc_summary$n_col_dropped <- length(drop_col_idx)
    
    curr_mat <- curr_mat[, -drop_col_idx, drop = FALSE]
  }
  
  # 4. Filter Multivariate Outliers (Per Group)
  if (!is.null(qc_config$remove_outliers) && qc_config$remove_outliers) {
    message("   [QC] Checking for multivariate outliers (PCA-based)...")
    
    if ("Group" %in% names(curr_meta)) {
      is_outlier <- detect_pca_outliers(curr_mat, curr_meta$Group, 
                                        conf_level = qc_config$outlier_conf_level)
      
      if (any(is_outlier)) {
        out_pids <- rownames(curr_mat)[is_outlier]
        message(sprintf("   [QC] Dropping %d outliers detected based on group distribution.", length(out_pids)))
        
        # Capture Original_Source info for report (Protected against NULL)
        target_info_col <- if(!is.null(stratification_col) && stratification_col %in% colnames(curr_meta)) stratification_col else "Group"
        group_info <- as.character(curr_meta[[target_info_col]][is_outlier])
        
        qc_summary$dropped_rows_detail <- rbind(qc_summary$dropped_rows_detail, data.frame(
          Patient_ID = out_pids,
          NA_Percent = round(rowMeans(is.na(curr_mat[is_outlier,,drop=FALSE])) * 100, 2),
          Reason = "Outlier",
          Original_Source = group_info,
          stringsAsFactors = FALSE
        ))
        
        curr_mat <- curr_mat[!is_outlier, , drop = FALSE]
        curr_meta <- curr_meta[!is_outlier, , drop = FALSE]
      }
    } else {
      warning("[QC] 'Group' column not found in metadata. Skipping outlier detection.")
    }
  }
  
  # Capture the mapping between Stratification Column (Subgroup) and Group for the report (Protected)
  group_map <- unique(data.frame(
    Subgroup = if(!is.null(stratification_col) && stratification_col %in% colnames(curr_meta)) as.character(curr_meta[[stratification_col]]) else "All",
    Group = if("Group" %in% colnames(curr_meta)) as.character(curr_meta$Group) else "All",
    stringsAsFactors = FALSE
  ))
  
  # Final Stats
  qc_summary$n_row_dropped <- nrow(qc_summary$dropped_rows_detail)
  qc_summary$n_row_final <- nrow(curr_mat)
  qc_summary$n_col_final <- ncol(curr_mat)
  qc_summary$final_markers_names <- colnames(curr_mat)
  
  qc_summary$breakdown_final <- get_group_counts(curr_meta, stratification_col)
  
  qc_summary$group_mapping <- group_map
  
  message(sprintf("   [QC] Final Dimensions: %d Samples x %d Markers", nrow(curr_mat), ncol(curr_mat)))
  
  return(list(
    data = curr_mat,
    metadata = curr_meta, 
    report = qc_summary
  ))
}