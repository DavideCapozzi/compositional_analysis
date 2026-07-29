# src/06_trustworthy_ai.R
# ==============================================================================
# STEP 06: TRUSTWORTHY / RESPONSIBLE-AI MODULE
# Description: Specification-curve robustness (P1, spine) + conformal calibrated
#              abstention (P2) + explanation-agreement audit (P3), with guarded
#              scaffolds for the clinical-data-gated pillars (P4/P5/P6).
#              Loads data_processed.rds + step03/step04 payloads. Outputs to
#              results_analysis/Trustworthy_AI/ and step06_payload.rds.
# ==============================================================================

source("R/utils_io.R")
source("R/modules_trustworthy.R")

suppressPackageStartupMessages({ library(ggplot2); library(openxlsx) })

message("\n=== PIPELINE STEP 6: TRUSTWORTHY AI ===")

config <- load_config("config/global_params.yml")
tw_cfg <- config$trustworthy
if (is.null(tw_cfg) || !isTRUE(tw_cfg$enabled)) {
  message("[TW] trustworthy module disabled in config; skipping.")
} else {

results_dir <- file.path(config$output_root, "results_analysis")
out_dir <- file.path(results_dir, "Trustworthy_AI")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

DATA <- readRDS(file.path(config$output_root, "01_data_processing", "data_processed.rds"))
payload03 <- readRDS(file.path(results_dir, "step03_payload.rds"))
payload04 <- readRDS(file.path(results_dir, "step04_payload.rds"))

seed <- if (!is.null(config$stats$seed)) config$stats$seed else 2025
set.seed(seed)

# Map scenario id -> case/control group sets from the analysis_scenarios config
scen_map <- setNames(config$analysis_scenarios, vapply(config$analysis_scenarios,
                                                       function(s) s$id, character(1)))

mv <- tw_cfg$multiverse
n_rep <- if (isTRUE(tw_cfg$quick)) mv$n_rep_quick else mv$n_rep_full
grid_cfg <- list(representation = mv$representation,
                 feature_space  = mv$feature_space,
                 edge_threshold = mv$edge_threshold,
                 estimator      = mv$estimator,
                 n_rep = n_rep, k = mv$k, n_perm = mv$n_perm)

message(sprintf("[TW] mode=%s | n_rep=%d | tasks: %s",
                ifelse(isTRUE(tw_cfg$quick), "quick", "full"),
                n_rep, paste(tw_cfg$tasks, collapse = ", ")))

reps <- tw_build_representations(DATA)

# helper: render the two spec-curve panels onto one PDF (gridExtra or 2 pages)
save_spec_pdf <- function(panels, path) {
  pdf(path, width = 9, height = 8)
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(panels$top, panels$bottom, heights = c(1, 1.1))
  } else {
    print(panels$top); print(panels$bottom)
  }
  dev.off()
}

payload06 <- list(specification_curve = list(), conformal = list(),
                  explanation = list(), gated = list())
summary_lines <- c(
  "# Trustworthy AI module — auto-generated summary",
  sprintf("_Generated: %s | mode: %s_", Sys.time(),
          ifelse(isTRUE(tw_cfg$quick), "quick", "full")),
  "",
  "Spine = specification-curve robustness (P1). Supports = conformal abstention",
  "(P2) and explanation-agreement audit (P3). Multiverse axes are leakage-free:",
  "representation x feature-space x edge-threshold x estimator. Group bootstrap",
  "networks are NOT used as per-patient features (would leak); the leakage-free",
  "structure proxy is the pairwise product z_i*z_j.", "")

for (task_id in tw_cfg$tasks) {
  scen <- scen_map[[task_id]]
  if (is.null(scen)) { message(sprintf("[TW] task %s not in analysis_scenarios; skipping.", task_id)); next }
  message(sprintf("\n--- Task %s: %s vs %s ---", task_id, scen$case_label, scen$control_label))

  # Guarded: repeated stratified CV, glmnet and pROC all error on very small or
  # single-class folds. A failing audit task must not abort the module.
  tryCatch({
  task <- tw_make_task(DATA, scen$case_groups, scen$control_groups,
                       config$stratification$column)

  # ---- P1: specification curve ----
  spec <- tw_run_specification_curve(reps, task, grid_cfg, seed = seed)
  m <- spec$meta
  message(sprintf("[P1] n=%d (%d/%d) | AUC %.3f-%.3f (median %.3f) | %d/%d > 0.60 | canonical perm-p=%.4f",
                  m$n, m$n_case, m$n_ctrl, m$auc_min, m$auc_max, m$auc_median,
                  m$n_above_060, m$n_specs, m$canonical_perm_p))
  panels <- tw_plot_spec_curve(spec, sprintf("Specification curve - %s (%s vs %s)",
                                             task_id, scen$case_label, scen$control_label))
  save_spec_pdf(panels, file.path(out_dir, sprintf("Specification_Curve_%s.pdf", task_id)))
  write.xlsx(spec$grid, file.path(out_dir, sprintf("Specification_Curve_%s.xlsx", task_id)))
  payload06$specification_curve[[task_id]] <- spec

  # ---- P2: conformal abstention (canonical z / abundance) ----
  Xz <- reps$z[task$idx, , drop = FALSE]
  cf <- tw_conformal(Xz, task$y, alpha = tw_cfg$conformal$alpha,
                     n_rep = tw_cfg$conformal$n_rep, seed = seed)
  message(sprintf("[P2] conformal @%.0f%%: coverage %.1f%% | informative singletons %.1f%% | abstention %.1f%%",
                  100 * cf$target_cov, 100 * cf$coverage, 100 * cf$singleton, 100 * cf$abstention))
  payload06$conformal[[task_id]] <- cf

  # ---- P3: explanation agreement (4 sources) ----
  mk <- reps$markers
  sources <- list(
    enet       = tw_importance_enet(Xz, task$y, mk),
    rf_perm    = tw_importance_rf(Xz, task$y, mk),
    splsda     = tw_get_splsda_importance(payload03, task_id, mk),
    net_degree = tw_get_network_degree(payload04, task_id, mk))
  agree <- tw_explanation_agreement(sources, top_k = tw_cfg$explanation$top_k)
  message(sprintf("[P3] explanation agreement across %d sources: mean Spearman %.2f | mean top-%d Jaccard %.2f",
                  length(agree$sources), agree$mean_rho, tw_cfg$explanation$top_k, agree$mean_jaccard))
  pdf(file.path(out_dir, sprintf("Explanation_Agreement_%s.pdf", task_id)), width = 7, height = 6)
  print(tw_plot_agreement(agree$rho, sprintf("Explanation rank agreement (Spearman) - %s", task_id), "Spearman"))
  print(tw_plot_agreement(agree$jaccard, sprintf("Explanation top-%d Jaccard - %s", tw_cfg$explanation$top_k, task_id), "Jaccard"))
  dev.off()
  payload06$explanation[[task_id]] <- agree

  # ---- summary block ----
  top_spec <- spec$grid[nrow(spec$grid), ]
  summary_lines <- c(summary_lines,
    sprintf("## %s — %s vs %s (n=%d: %d case / %d control)",
            task_id, scen$case_label, scen$control_label, m$n, m$n_case, m$n_ctrl),
    sprintf("- **P1 robustness:** AUC %.3f-%.3f (median %.3f) across %d specs; %d/%d > 0.60; canonical perm-p = %.4f.",
            m$auc_min, m$auc_max, m$auc_median, m$n_specs, m$n_above_060, m$n_specs, m$canonical_perm_p),
    sprintf("- **P1 best spec:** %s / %s / %s (AUC %.3f) — top of the curve.",
            top_spec$representation, top_spec$feature_space, top_spec$estimator, top_spec$AUC),
    sprintf("- **P2 abstention:** coverage %.1f%% (target %.0f%%); informative singletons %.1f%%; abstention %.1f%%.",
            100 * cf$coverage, 100 * cf$target_cov, 100 * cf$singleton, 100 * cf$abstention),
    sprintf("- **P3 explanation agreement:** mean Spearman %.2f, mean top-%d Jaccard %.2f across {%s}.",
            agree$mean_rho, tw_cfg$explanation$top_k, agree$mean_jaccard, paste(agree$sources, collapse = ", ")),
    "")

  }, error = function(e) {
    message(sprintf("      [WARN] Trustworthy task '%s' failed: %s", task_id, conditionMessage(e)))
    summary_lines <<- c(summary_lines,
      sprintf("## %s — SKIPPED", task_id),
      sprintf("- Audit failed: %s", conditionMessage(e)), "")
  })
}

# ---- conformal & explanation roll-ups to Excel ----
conf_df <- do.call(rbind, lapply(names(payload06$conformal), function(id) {
  c0 <- payload06$conformal[[id]]
  data.frame(task = id, target_coverage = c0$target_cov, coverage = round(c0$coverage, 3),
             singleton = round(c0$singleton, 3), abstention = round(c0$abstention, 3))
}))
if (!is.null(conf_df)) write.xlsx(conf_df, file.path(out_dir, "Conformal_Coverage.xlsx"))

# ---- DATA GATE: P4/P5/P6 ----
gate <- tw_detect_clinical_meta(DATA$metadata, tw_cfg)
payload06$gated$gate <- gate
gate_msgs <- c()
if (isTRUE(tw_cfg$confound$enabled) && !gate$confound)
  gate_msgs <- c(gate_msgs, "- **P4 confound** (batch/biology): PENDING — no acquisition batch/date column in metadata.")
if (isTRUE(tw_cfg$survival$enabled) && !gate$survival)
  gate_msgs <- c(gate_msgs, "- **P5 survival** (calibrated UQ): PENDING — no time-to-event columns in metadata.")
if (isTRUE(tw_cfg$fairness$enabled) && !gate$fairness)
  gate_msgs <- c(gate_msgs, "- **P6 fairness** (subgroup parity): PENDING — no sex/age columns in metadata.")
if (length(gate_msgs)) {
  message("[TW] DATA-GATED pillars pending clinical metadata file:")
  for (gm in gate_msgs) message("     ", sub("^- ", "", gsub("\\*\\*", "", gm)))
  summary_lines <- c(summary_lines, "## Data-gated pillars (awaiting clinical metadata)", gate_msgs, "")
}

# ---- persist ----
saveRDS(payload06, file.path(results_dir, "step06_payload.rds"))
writeLines(summary_lines, file.path(out_dir, "Trustworthy_AI_Summary.md"))
message(sprintf("\n[TW] Done. Outputs in %s", out_dir))

}  # end enabled guard
