# src/01_data_processing.R
# ==============================================================================
# STEP 01: DATA PROCESSING
# Description: Data ingestion + QC + hybrid coda preprocessing
# ==============================================================================

source("R/utils_io.R")     
source("R/modules_coda.R") 
source("R/modules_qc.R")   

message("\n=== PIPELINE STEP 1: INGESTION + QC + HYBRID TRANSFORM ===")

# 1. Load Configuration & Data
# ------------------------------------------------------------------------------
config <- load_config("config/global_params.yml")

raw_data <- load_raw_data(config)
panel_map <- attr(raw_data, "panels")

# 1.1 Cohort Integrity Validation
# ------------------------------------------------------------------------------
# These are hard stops by design. Every one of them guards a failure mode that would
# otherwise let the pipeline finish "successfully" with silently wrong results.
validate_cohort_integrity <- function(raw_data, config) {

  strat_col <- config$stratification$column

  # (a) Every configured cohort actually contributed rows.
  present_cohorts <- unique(as.character(raw_data$Group))
  missing_cohorts <- setdiff(names(config$cohorts), present_cohorts)
  if (length(missing_cohorts) > 0) {
    stop(sprintf("[FATAL] Configured cohort(s) contributed zero samples: %s",
                 paste(missing_cohorts, collapse = ", ")))
  }

  # (b) Patient_ID unique across all sheets. Row-name indexing is used pervasively
  # downstream and returns only the FIRST match for a duplicated ID.
  dup_ids <- raw_data$Patient_ID[duplicated(raw_data$Patient_ID)]
  if (length(dup_ids) > 0) {
    stop(sprintf("[FATAL] %d duplicate Patient_ID(s) across cohort sheets: %s.\nIDs must be unique (prefix them per cohort).",
                 length(unique(dup_ids)), paste(unique(dup_ids), collapse = ", ")))
  }

  # (c) Every group in the data is claimed by the analysis config. A group present in the
  # data but absent from case_groups is excluded from every plot and statistic while STILL
  # shifting the global z-scoring for everyone else.
  declared <- unique(c(unlist(config$control_group), unlist(config$case_groups)))
  observed <- unique(as.character(raw_data[[strat_col]]))
  undeclared <- setdiff(observed, declared)
  if (length(undeclared) > 0) {
    stop(sprintf("[FATAL] Group(s) present in the data but absent from control_group/case_groups: %s.\nThey would be silently excluded from all analyses while still affecting the global scaling.",
                 paste(undeclared, collapse = ", ")))
  }

  # (d) Scenario ids are unique. They become directory names and payload list keys, and
  # collide silently (overwriting each other's outputs).
  scen_ids <- vapply(config$analysis_scenarios, function(s) s$id, character(1))
  if (any(duplicated(scen_ids))) {
    stop(sprintf("[FATAL] Duplicate analysis_scenarios id(s): %s",
                 paste(unique(scen_ids[duplicated(scen_ids)]), collapse = ", ")))
  }

  # (e) Scenario labels must differ within a scenario (duplicate factor levels error out
  # in workflows.R) and every referenced group must exist in the data.
  for (s in config$analysis_scenarios) {
    if (identical(s$case_label, s$control_label)) {
      stop(sprintf("[FATAL] Scenario '%s' has case_label == control_label ('%s').",
                   s$id, s$case_label))
    }
    unknown <- setdiff(unlist(c(s$case_groups, s$control_groups)), observed)
    if (length(unknown) > 0) {
      stop(sprintf("[FATAL] Scenario '%s' references group(s) not present in the data: %s",
                   s$id, paste(unknown, collapse = ", ")))
    }
  }

  message(sprintf("[Check] Cohort integrity OK: %d cohorts, %d samples, %d groups, %d scenarios.",
                  length(present_cohorts), nrow(raw_data), length(observed), length(scen_ids)))
}

validate_cohort_integrity(raw_data, config)

# Snapshot pre-QC counts so we can verify no cohort is decimated by filtering.
counts_pre_qc <- table(as.character(raw_data$Group))

# 2. Setup Initial Matrices
# ------------------------------------------------------------------------------
# Split Metadata and Numeric Matrix
subgroup_col <- config$metadata$subgroup_col
meta_cols <- c("Patient_ID", "Group")
if (!is.null(subgroup_col)) meta_cols <- c(meta_cols, subgroup_col)

# Extract Matrix and Metadata aligned
mat_raw <- as.matrix(raw_data[, setdiff(names(raw_data), meta_cols)])
rownames(mat_raw) <- raw_data$Patient_ID
metadata_raw <- raw_data[, meta_cols]

message(sprintf("[Data] Initial Matrix: %d Samples x %d Markers", nrow(mat_raw), ncol(mat_raw)))

# --- Sample Filtering (Blacklist) ---
target_samples <- rownames(mat_raw)
dropped_samples_apriori <- data.frame()

if (!is.null(config$sample_selection$blacklist) && length(config$sample_selection$blacklist) > 0) {
  blacklist_samples <- config$sample_selection$blacklist
  
  # Identify samples to drop that actually exist in data
  samples_to_drop <- intersect(target_samples, blacklist_samples)
  
  if (length(samples_to_drop) > 0) {
    message(sprintf("[Data] Dropping %d samples defined in blacklist.", length(samples_to_drop)))
    
    # Store details for QC report
    # [FIX] Use match() to align indices and [[]] subsetting to prevent tibble string coercion
    dropped_indices <- match(samples_to_drop, metadata_raw$Patient_ID)
    
    dropped_samples_apriori <- data.frame(
      Patient_ID = samples_to_drop,
      NA_Percent = NA, 
      Reason = "A Priori",
      Original_Source = if (!is.null(subgroup_col) && subgroup_col %in% colnames(metadata_raw)) {
        as.character(metadata_raw[[subgroup_col]][dropped_indices])
      } else {
        "N/A"
      },
      stringsAsFactors = FALSE
    )
    
    # Apply filter
    target_samples <- setdiff(target_samples, samples_to_drop)
    mat_raw <- mat_raw[target_samples, , drop = FALSE]
    metadata_raw <- metadata_raw[match(target_samples, metadata_raw$Patient_ID), , drop = FALSE]
    raw_data <- raw_data %>% filter(Patient_ID %in% target_samples)
  }
}

# --- Marker Filtering (Whitelist/Blacklist) ---
target_markers <- colnames(mat_raw)
if (!is.null(config$marker_selection$whitelist) && length(config$marker_selection$whitelist) > 0) {
  target_markers <- intersect(target_markers, config$marker_selection$whitelist)
}
if (!is.null(config$marker_selection$blacklist) && length(config$marker_selection$blacklist) > 0) {
  target_markers <- setdiff(target_markers, config$marker_selection$blacklist)
}

mat_raw <- mat_raw[, target_markers, drop = FALSE]

# 2.1 Data Assumption Validation
# ------------------------------------------------------------------------------
# Strictly validate data range against configuration to prevent mathematical errors
# in Logit transformation.

target_format <- if (!is.null(config$input_format)) config$input_format else "percentage"
max_val <- max(mat_raw, na.rm = TRUE)

message(sprintf("[Check] Input format configured as: '%s'. Max observed value: %.2f", target_format, max_val))

if (target_format == "percentage") {
  # Expectation: Values 0-100.
  # Warning if max is small (looks like proportion), but allow proceeding (rare populations case).
  if (max_val <= 1.0) {
    warning(sprintf(
      "\n[CAUTION] Config is 'percentage' but max value is %.2f. \nIf these are actually Proportions (0-1), pipeline will divide by 100, crushing variance.\nPlease verify your raw input.", 
      max_val
    ))
  }
} else if (target_format == "proportion") {
  # Expectation: Values 0-1.
  # Critical Error if max > 1.0 (Logit of >1 is NaN).
  if (max_val > 1.0) {
    stop(sprintf(
      "[FATAL] Config is 'proportion' but detected values > 1.0 (Max=%.2f). \nPlease change config to 'percentage' or check raw data.", 
      max_val
    ))
  }
} else {
  stop("[FATAL] Invalid 'input_format' in config. Use 'percentage' or 'proportion'.")
}

# 3. Geometric-Aware Quality Control (Split Strategy)
# ------------------------------------------------------------------------------
message("[QC] Running Split QC Strategy...")

# --- STEP A: Basic QC (Missingness & Zero Variance on RAW DATA) ---
# Temporarily disable outlier detection to focus strictly on cleaning NAs/Zeros first.
# This prevents epsilon-imputation from masking high-missingness samples.
qc_config_basic <- config$qc
qc_config_basic$remove_outliers <- FALSE

message("   [QC-A] Filtering Missingness and Zero Variance on RAW data...")
qc_result <- run_qc_pipeline(
  mat_raw = mat_raw,
  metadata = metadata_raw,
  qc_config = qc_config_basic,
  stratification_col = config$metadata$subgroup_col,
  dropped_markers_apriori = data.frame(),
  dropped_samples_apriori = dropped_samples_apriori,
  panel_map = panel_map
)

# Extract cleaned datasets from Step A
mat_clean_basic  <- qc_result$data
meta_clean_basic <- qc_result$metadata

# --- STEP B: Geometric Outlier Detection (on PROXY DATA) ---
# Now that NAs are handled, we generate the proxy for multivariate checks.
if (config$qc$remove_outliers) {
  message("   [QC-B] Generating geometric proxy for Outlier Detection...")
  
  # 1. Generate Proxy (Fast Mode) on CLEAN data
  proxy_results <- perform_data_transformation(mat_clean_basic, config, mode = "fast")
  mat_proxy <- proxy_results$hybrid_data_z
  
  # 2. Detect Outliers using the Proxy
  message("   [QC-B] Detecting multivariate outliers (PCA-based)...")
  
  # Determine correct grouping column
  group_col <- if(config$metadata$subgroup_col %in% colnames(meta_clean_basic)) config$metadata$subgroup_col else "Group"
  
  is_outlier <- detect_pca_outliers(
    mat = mat_proxy, 
    groups = meta_clean_basic[[group_col]], 
    conf_level = config$qc$outlier_conf_level
  )
  
  # 3. Apply Outlier Filter & Update Report
  if (any(is_outlier)) {
    out_pids <- rownames(mat_clean_basic)[is_outlier]
    message(sprintf("   [QC-B] Dropping %d multivariate outliers.", length(out_pids)))
    
    # Capture info for report
    group_info <- as.character(meta_clean_basic[[group_col]][is_outlier])
    
    # Calculate NA % for these outliers (on the basic clean matrix for transparency)
    na_pcts <- rowMeans(is.na(mat_clean_basic[is_outlier, , drop=FALSE])) * 100
    
    # Append to existing report details
    new_drops <- data.frame(
      Patient_ID = out_pids,
      NA_Percent = round(na_pcts, 2),
      Reason = "Outlier",
      Original_Source = group_info,
      stringsAsFactors = FALSE
    )
    
    qc_result$report$dropped_rows_detail <- rbind(qc_result$report$dropped_rows_detail, new_drops)
    qc_result$report$n_row_dropped <- nrow(qc_result$report$dropped_rows_detail)
    
    # Apply removal
    mat_clean_basic  <- mat_clean_basic[!is_outlier, , drop = FALSE]
    meta_clean_basic <- meta_clean_basic[!is_outlier, , drop = FALSE]
    
    # Update Final Dimensions in Report
    qc_result$report$n_row_final <- nrow(mat_clean_basic)
    
    # Recalculate Final Breakdown
    if (!is.null(qc_result$report$breakdown_final)) {
      qc_result$report$breakdown_final <- table(meta_clean_basic[[config$metadata$subgroup_col]])
    }
  }
}

# --- STEP C: Synchronization ---
# Sync main objects with the final result of the Split QC
mat_raw <- mat_clean_basic
metadata_raw <- meta_clean_basic
raw_data <- raw_data %>% filter(Patient_ID %in% rownames(mat_raw))

message(sprintf("[QC] Final Dimensions: %d Samples x %d Markers", nrow(mat_raw), ncol(mat_raw)))

# --- Post-QC Retention Check ---
# A cohort whose assay panel differs from the others can be silently gutted by missingness
# filtering. Fail loudly rather than analysing a cohort that has quietly disappeared.
counts_post_qc <- table(as.character(metadata_raw$Group))
max_loss <- 0.20
for (coh in names(counts_pre_qc)) {
  n_pre  <- as.integer(counts_pre_qc[[coh]])
  n_post <- if (coh %in% names(counts_post_qc)) as.integer(counts_post_qc[[coh]]) else 0L
  loss   <- 1 - (n_post / n_pre)

  message(sprintf("   [QC] Retention %-16s %3d -> %3d (%.0f%% lost)", coh, n_pre, n_post, loss * 100))

  if (loss > max_loss) {
    stop(sprintf(
      "[FATAL] Cohort '%s' lost %.0f%% of its samples during QC (%d -> %d), above the %.0f%% ceiling.\nThis usually means its assay panel differs from the other cohorts. Inspect QC_Filtering_Report.xlsx before proceeding.",
      coh, loss * 100, n_pre, n_post, max_loss * 100))
  }
}

# Calculate imputation details on the final filtered matrix before saving report
impute_info_df <- data.frame(
  Patient_ID = character(), 
  NA_Percent = numeric(), 
  Original_Source = character(), 
  Imputed_Markers = character(),
  stringsAsFactors = FALSE
)

# Identify samples with NAs
na_counts_per_row <- rowSums(is.na(mat_raw))
rows_with_na <- which(na_counts_per_row > 0)

if (length(rows_with_na) > 0) {
  pids_na <- rownames(mat_raw)[rows_with_na]
  
  # Determine Group/Source column safely
  grp_col <- if(config$metadata$subgroup_col %in% colnames(metadata_raw)) config$metadata$subgroup_col else "Group"
  
  # Construct list of markers per patient
  markers_list <- apply(mat_raw[rows_with_na, , drop=FALSE], 1, function(row_vals) {
    paste(names(row_vals)[is.na(row_vals)], collapse = ", ")
  })
  
  # Match metadata
  meta_indices <- match(pids_na, metadata_raw$Patient_ID)
  sources <- as.character(metadata_raw[[grp_col]][meta_indices])
  
  # Calculate Percentages
  na_pcts <- round((na_counts_per_row[rows_with_na] / ncol(mat_raw)) * 100, 2)
  
  impute_info_df <- data.frame(
    Patient_ID = pids_na,
    NA_Percent = na_pcts,
    Original_Source = sources,
    Imputed_Markers = markers_list,
    stringsAsFactors = FALSE
  )
}

# Inject into report object
qc_result$report$imputed_details <- impute_info_df

# 4. Final Hybrid Transformation (Complete Mode)
# ------------------------------------------------------------------------------
# Now running on strictly filtered data.
# Mode="complete" enables: LOD zero handling, BPCA (imputation), and Logit transformation.

na_rates <- rowMeans(is.na(mat_raw))
impute_indices <- which(na_rates > 0)

if (length(impute_indices) > 0) {
  # Sort indices by missingness (descending) to highlight heavily imputed samples first
  impute_indices <- impute_indices[order(na_rates[impute_indices], decreasing = TRUE)]
  
  # Format string: "PatientID (XX%)"
  log_entries <- sprintf("%s (%.0f%%)", 
                         names(impute_indices), 
                         na_rates[impute_indices] * 100)
  
  message(sprintf("\n[Impute] The following %d samples contain missing values and will undergo BPCA Imputation:", length(log_entries)))
  message(paste(log_entries, collapse = ", "))
  message("") # Add spacing for readability
} else {
  message("\n[Impute] No missing values detected in the filtered dataset. BPCA step will be skipped internally.")
}

message("[Transform] Running Final Data Transformation (Complete Mode) on clean data...")
transform_results <- perform_data_transformation(mat_raw, config, mode = "complete")

# Extract results
mat_hybrid_raw <- transform_results$hybrid_data_raw
mat_hybrid_z   <- transform_results$hybrid_data_z
final_markers  <- transform_results$hybrid_markers
clamped_mask   <- transform_results$clamped_mask

# --- Clamping Details Extraction ---
clamped_info_df <- data.frame(
  Patient_ID = character(), 
  Original_Source = character(), 
  Clamped_Markers = character(),
  stringsAsFactors = FALSE
)

# Identify samples with any clamped markers
if (!is.null(clamped_mask)) {
  rows_with_clamps <- which(rowSums(clamped_mask, na.rm = TRUE) > 0)
  
  if (length(rows_with_clamps) > 0) {
    pids_clamp <- rownames(mat_raw)[rows_with_clamps]
    
    # Determine Group/Source column safely
    grp_col <- if(config$metadata$subgroup_col %in% colnames(metadata_raw)) config$metadata$subgroup_col else "Group"
    meta_indices <- match(pids_clamp, metadata_raw$Patient_ID)
    sources <- as.character(metadata_raw[[grp_col]][meta_indices])
    
    # Construct list of markers per patient
    markers_list <- apply(clamped_mask[rows_with_clamps, , drop=FALSE], 1, function(row_vals) {
      paste(names(row_vals)[row_vals & !is.na(row_vals)], collapse = ", ")
    })
    
    clamped_info_df <- data.frame(
      Patient_ID = pids_clamp,
      Original_Source = sources,
      Clamped_Markers = markers_list,
      stringsAsFactors = FALSE
    )
  }
}

# Inject clamped details into report object
qc_result$report$clamped_details <- clamped_info_df

# Save QC Report (Moved here to include post-transformation data)
out_dir <- file.path(config$output_root, "01_data_processing")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
save_qc_report(qc_result$report, file.path(out_dir, "QC_Filtering_Report.xlsx"), config = config)

# 5. Final Safety Polish
# ------------------------------------------------------------------------------
# Fallback for any residual NAs (e.g., if BPCA didn't converge for a specific point).
if (any(is.na(mat_hybrid_z))) {
  n_na_res <- sum(is.na(mat_hybrid_z))
  message(sprintf("[WARN] Detected %d residual NAs in final matrix. Applying median fallback...", n_na_res))
  
  mat_hybrid_z <- apply(mat_hybrid_z, 2, function(x) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
    return(x)
  })
}

# 6. Save Final Output
# ------------------------------------------------------------------------------
df_hybrid_raw <- cbind(raw_data[, meta_cols], as.data.frame(mat_hybrid_raw))
df_hybrid_z   <- cbind(raw_data[, meta_cols], as.data.frame(mat_hybrid_z))

processed_data <- list(
  metadata        = raw_data[, meta_cols],
  markers         = colnames(mat_raw),      
  raw_matrix      = mat_raw, 
  hybrid_markers  = final_markers,
  hybrid_data_raw = df_hybrid_raw,
  hybrid_data_z   = df_hybrid_z,
  config          = config
)

saveRDS(processed_data, file.path(out_dir, "data_processed.rds"))
message("=== STEP 1 COMPLETE ===\n")