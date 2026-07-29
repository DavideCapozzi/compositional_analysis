# R/modules_trustworthy.R
# ==============================================================================
# TRUSTWORTHY / RESPONSIBLE-AI MODULE  (Phase 6)
# Description: Functions supporting the three CORE trust pillars for the
#              ELLIS summer-school module, plus guarded scaffolds for the
#              clinical-data-gated pillars:
#                P1  Specification-curve / multiverse robustness   (SPINE)
#                P2  Conformal calibrated uncertainty / abstention
#                P3  Explanation-agreement audit (4 importance sources)
#                P4  Confound-vs-biology adjudication        (DATA-GATED)
#                P5  Calibrated survival modelling w/ UQ     (DATA-GATED)
#                P6  Fairness audit across subgroups         (DATA-GATED)
#
# Design notes (defensibility):
#  - The multiverse varies LEAKAGE-FREE, reproducible axes only:
#    representation (z / raw-median / raw-z) x feature-space
#    (abundance / +edges / edges-only) x estimator (ridge/enet/lasso/rf),
#    with an unsupervised edge-threshold sub-axis. It does NOT re-bootstrap
#    group networks per cell: group baselines are group-level structures and
#    cannot be per-patient features without leakage. The leakage-free proxy
#    for "network structure" is the pairwise product z_i*z_j (validated in
#    diagnostics/probe_network_features.R).
#  - glmnet standardises internally per fit, so representation scaling does
#    not leak through supervised standardisation; RF is scale-invariant.
#  - Edge pre-selection uses |cor| on the POOLED task matrix (ignores y) ->
#    unsupervised -> leakage-free.
# Dependencies: glmnet, ranger, pROC, ggplot2 (survival/ranger for P5).
# ==============================================================================

suppressPackageStartupMessages({
  library(glmnet)
  library(ranger)
  library(pROC)
})

# --- low-level helpers --------------------------------------------------------

#' AUC of P(case) against a control/case factor.
tw_auc <- function(truth, prob_case) {
  as.numeric(pROC::auc(response = truth, predictor = prob_case,
                       levels = c("control", "case"), direction = "<",
                       quiet = TRUE))
}

#' Stratified k-fold assignment preserving class balance.
tw_strat_folds <- function(y, k) {
  f <- integer(length(y))
  for (lv in levels(y)) {
    i <- which(y == lv); i <- sample(i)
    f[i] <- rep_len(seq_len(k), length(i))
  }
  f
}

#' Median lambda.min over a few CV splits (reused by perm / conformal).
tw_pick_lambda <- function(X, y, alpha = 0.5) {
  ls <- replicate(5, {
    cv <- glmnet::cv.glmnet(X, y, family = "binomial", alpha = alpha, nfolds = 5)
    cv$lambda.min
  })
  stats::median(ls)
}

# --- edge (network-structure proxy) features ---------------------------------

#' All unordered marker pairs as a 2 x K index matrix.
tw_all_pairs <- function(p) utils::combn(p, 2)

#' Unsupervised edge pre-selection: keep pairs whose |Pearson cor| on the
#' pooled matrix X meets `threshold`. Ignores the outcome -> leakage-free.
tw_select_edge_pairs <- function(X, threshold = 0) {
  pr <- tw_all_pairs(ncol(X))
  if (threshold <= 0) return(pr)
  C <- suppressWarnings(stats::cor(X))
  keep <- abs(C[cbind(pr[1, ], pr[2, ])]) >= threshold
  keep[is.na(keep)] <- FALSE
  if (!any(keep)) return(pr[, integer(0), drop = FALSE])
  pr[, keep, drop = FALSE]
}

#' Build pairwise-product (edge-proxy) features for the given pair index matrix.
tw_edge_matrix <- function(X, pairs) {
  if (ncol(pairs) == 0) return(matrix(numeric(0), nrow = nrow(X), ncol = 0))
  M <- X[, pairs[1, ]] * X[, pairs[2, ]]
  M <- matrix(M, nrow = nrow(X))
  colnames(M) <- paste0(colnames(X)[pairs[1, ]], ":", colnames(X)[pairs[2, ]])
  M
}

#' Assemble the design matrix for one multiverse cell.
tw_assemble_features <- function(Xrep, feature_space, edge_threshold) {
  if (feature_space == "abundance") return(Xrep)
  pairs <- tw_select_edge_pairs(Xrep, edge_threshold)
  E <- tw_edge_matrix(Xrep, pairs)
  if (feature_space == "edges_only") return(E)
  cbind(Xrep, E)               # with_edges
}

# --- cross-validated AUC for a given estimator -------------------------------

#' One repeat of stratified k-fold CV AUC. `estimator` in
#' {ridge, enet, lasso, rf}. glmnet uses nested lambda selection per fold.
tw_cv_auc <- function(X, y, estimator, k = 5) {
  f <- tw_strat_folds(y, k); p <- numeric(length(y))
  for (kk in seq_len(k)) {
    tr <- f != kk; te <- !tr
    if (estimator == "rf") {
      rf <- ranger::ranger(x = X[tr, , drop = FALSE], y = y[tr],
                           probability = TRUE, num.trees = 400)
      p[te] <- predict(rf, X[te, , drop = FALSE])$predictions[, "case"]
    } else {
      a <- switch(estimator, ridge = 0, enet = 0.5, lasso = 1)
      cv <- glmnet::cv.glmnet(X[tr, , drop = FALSE], y[tr], family = "binomial",
                              alpha = a, nfolds = 5)
      p[te] <- as.numeric(predict(cv, X[te, , drop = FALSE],
                                  s = "lambda.min", type = "response"))
    }
  }
  tw_auc(y, p)
}

#' Median CV AUC over `n_rep` repeats.
tw_repeated_auc <- function(X, y, estimator, n_rep = 10, k = 5) {
  stats::median(replicate(n_rep, tw_cv_auc(X, y, estimator, k)))
}

#' Permutation-null p-value for an observed AUC (enet, fixed lambda).
tw_perm_pvalue <- function(X, y, obs, n_perm = 200, k = 5) {
  lam <- tw_pick_lambda(X, y)
  null <- replicate(n_perm, {
    yp <- sample(y); f <- tw_strat_folds(yp, k); p <- numeric(length(yp))
    for (kk in seq_len(k)) {
      tr <- f != kk; te <- !tr
      fit <- glmnet::glmnet(X[tr, , drop = FALSE], yp[tr], family = "binomial",
                            alpha = 0.5, lambda = lam)
      p[te] <- as.numeric(predict(fit, X[te, , drop = FALSE], type = "response"))
    }
    tw_auc(yp, p)
  })
  (sum(null >= obs) + 1) / (n_perm + 1)
}

# --- representations & tasks --------------------------------------------------

#' Build the leakage-free preprocessing representations from the processed
#' object. All matrices are row-aligned to `hybrid_data_z`.
#'   z       : logit + BPCA + z-score (pipeline canonical)
#'   raw_med : raw proportions (LOD-adjusted), median-imputed
#'   raw_z   : z-scored raw proportions (transform = none)
tw_build_representations <- function(DATA) {
  mk <- DATA$hybrid_markers
  Z  <- as.matrix(DATA$hybrid_data_z[, mk, drop = FALSE])
  R  <- as.matrix(DATA$hybrid_data_raw[, mk, drop = FALSE])
  for (j in seq_len(ncol(R))) {                 # median imputation for raw path
    na <- is.na(R[, j]); if (any(na)) R[na, j] <- stats::median(R[, j], na.rm = TRUE)
  }
  Rz <- scale(R)
  attr(Rz, "scaled:center") <- NULL; attr(Rz, "scaled:scale") <- NULL
  list(z = Z, raw_med = R, raw_z = as.matrix(Rz), markers = mk)
}

#' Row indices and binary outcome for a case-vs-control task defined by
#' Original_Source group labels.
tw_make_task <- function(DATA, case_groups, control_groups,
                         strat_col = "Original_Source") {
  src <- DATA$hybrid_data_z[[strat_col]]
  idx <- which(src %in% c(case_groups, control_groups))
  y   <- factor(ifelse(src[idx] %in% case_groups, "case", "control"),
                levels = c("control", "case"))
  # Cancer type = the root of the Original_Source label ("NSCLC_LS" -> "NSCLC",
  # "Healthy_Donors" -> "Healthy"). Derived generically so a newly added cohort is not
  # silently bucketed into an existing type.
  ctype <- vapply(strsplit(as.character(src[idx]), "_"), `[`, character(1), 1)

  list(idx = idx, y = y,
       src = src[idx],
       ctype = ctype)
}

# --- P1: specification curve --------------------------------------------------

#' Run the multiverse specification curve for one task.
#' @param reps  output of tw_build_representations()
#' @param task  output of tw_make_task()
#' @param grid_cfg list with vectors: representation, feature_space,
#'                 edge_threshold, estimator; plus n_rep, k.
#' @return list(grid = tidy data.frame with AUC per spec, top = top markers of
#'              the best abundance/enet spec, meta = list)
tw_run_specification_curve <- function(reps, task, grid_cfg, seed = 2025) {
  set.seed(seed)
  y <- task$y; idx <- task$idx
  rep_mats <- list(z = reps$z[idx, , drop = FALSE],
                   raw_med = reps$raw_med[idx, , drop = FALSE],
                   raw_z = reps$raw_z[idx, , drop = FALSE])

  grid <- expand.grid(representation = grid_cfg$representation,
                      feature_space  = grid_cfg$feature_space,
                      edge_threshold = grid_cfg$edge_threshold,
                      estimator      = grid_cfg$estimator,
                      stringsAsFactors = FALSE)
  # edge_threshold is meaningless for abundance-only specs -> collapse to 0
  grid$edge_threshold[grid$feature_space == "abundance"] <- 0
  grid <- unique(grid)

  grid$AUC <- NA_real_; grid$n_feat <- NA_integer_
  for (i in seq_len(nrow(grid))) {
    Xrep <- rep_mats[[grid$representation[i]]]
    X <- tw_assemble_features(Xrep, grid$feature_space[i], grid$edge_threshold[i])
    if (ncol(X) < 2) next
    grid$n_feat[i] <- ncol(X)
    grid$AUC[i] <- tw_repeated_auc(X, y, grid$estimator[i],
                                   n_rep = grid_cfg$n_rep, k = grid_cfg$k)
  }
  grid <- grid[!is.na(grid$AUC), , drop = FALSE]
  grid <- grid[order(grid$AUC), ]
  grid$spec_id <- seq_len(nrow(grid))

  # headline permutation test on the canonical abundance/enet/z spec
  base_X <- rep_mats$z
  obs <- tw_repeated_auc(base_X, y, "enet", n_rep = grid_cfg$n_rep, k = grid_cfg$k)
  pval <- tw_perm_pvalue(base_X, y, obs, n_perm = grid_cfg$n_perm, k = grid_cfg$k)

  list(grid = grid,
       meta = list(
         n = length(y), n_case = sum(y == "case"), n_ctrl = sum(y == "control"),
         auc_min = min(grid$AUC), auc_median = stats::median(grid$AUC),
         auc_max = max(grid$AUC),
         n_above_060 = sum(grid$AUC > 0.60), n_specs = nrow(grid),
         frac_edge_in_top5 = mean(utils::tail(grid$feature_space, 5) != "abundance"),
         canonical_auc = obs, canonical_perm_p = pval))
}

#' Canonical multiverse figure: sorted AUC (top) over the choice grid (bottom).
tw_plot_spec_curve <- function(spec, title) {
  g <- spec$grid
  p_top <- ggplot2::ggplot(g, ggplot2::aes(spec_id, AUC)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 2, colour = "grey60") +
    ggplot2::geom_point(ggplot2::aes(colour = feature_space), size = 1.6) +
    ggplot2::scale_colour_manual(values = c(abundance = "#4682B4",
                                            with_edges = "#CD5C5C",
                                            edges_only = "#9932CC")) +
    ggplot2::labs(title = title,
                  subtitle = sprintf("%d specifications | AUC %.3f-%.3f (median %.3f) | %d/%d > 0.60",
                                     nrow(g), spec$meta$auc_min, spec$meta$auc_max,
                                     spec$meta$auc_median, spec$meta$n_above_060,
                                     spec$meta$n_specs),
                  x = NULL, y = "CV AUC", colour = "Features") +
    ggplot2::theme_bw() + ggplot2::ylim(min(0.4, min(g$AUC)), max(1, max(g$AUC)))

  long <- do.call(rbind, lapply(
    c("representation", "estimator", "feature_space"),
    function(ax) data.frame(spec_id = g$spec_id, axis = ax,
                            level = as.character(g[[ax]]), stringsAsFactors = FALSE)))
  p_bot <- ggplot2::ggplot(long, ggplot2::aes(spec_id, level, fill = level)) +
    ggplot2::geom_tile() +
    ggplot2::facet_grid(axis ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(x = "specification (sorted by AUC)", y = NULL) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "none")

  list(top = p_top, bottom = p_bot)
}

# --- P2: split-conformal calibrated abstention --------------------------------

#' Split-conformal coverage & informativeness for a task. Reports empirical
#' coverage, singleton (informative) rate, and abstention (empty + both) rate.
tw_conformal <- function(X, y, alpha = 0.10, n_rep = 200, seed = 2025) {
  set.seed(seed)
  lam <- tw_pick_lambda(X, y)
  cov <- single <- empty <- both <- numeric(n_rep)
  for (r in seq_len(n_rep)) {
    f <- tw_strat_folds(y, 5)
    tr <- f %in% 1:3; cal <- f == 4; te <- f == 5
    if (length(unique(y[tr])) < 2) { cov[r] <- NA; next }
    fit <- glmnet::glmnet(X[tr, , drop = FALSE], y[tr], family = "binomial",
                          alpha = 0.5, lambda = lam)
    pc <- function(ix) as.numeric(predict(fit, X[ix, , drop = FALSE], type = "response"))
    pcal <- pc(which(cal)); ycal <- y[cal]
    scal <- ifelse(ycal == "case", 1 - pcal, pcal)
    qh <- sort(scal)[ceiling((length(scal) + 1) * (1 - alpha))]
    if (is.na(qh)) qh <- max(scal)
    pte <- pc(which(te)); yte <- y[te]
    in_case <- (1 - pte) <= qh; in_ctrl <- pte <= qh
    set_size <- in_case + in_ctrl
    cov[r]    <- mean(ifelse(yte == "case", in_case, in_ctrl))
    single[r] <- mean(set_size == 1)
    empty[r]  <- mean(set_size == 0)
    both[r]   <- mean(set_size == 2)
  }
  list(alpha = alpha, target_cov = 1 - alpha,
       coverage = mean(cov, na.rm = TRUE),
       singleton = mean(single, na.rm = TRUE),
       empty = mean(empty, na.rm = TRUE),
       both = mean(both, na.rm = TRUE),
       abstention = mean(empty, na.rm = TRUE) + mean(both, na.rm = TRUE))
}

# --- P3: explanation-agreement audit -----------------------------------------

#' Elastic-net |coef| importance over markers (full-data fit).
tw_importance_enet <- function(X, y, markers) {
  cv <- glmnet::cv.glmnet(X, y, family = "binomial", alpha = 0.5, nfolds = 5)
  co <- as.matrix(coef(cv, s = "lambda.min"))
  v <- abs(co[markers, 1]); v[is.na(v)] <- 0; v
}

#' RF permutation importance over markers (full-data fit).
tw_importance_rf <- function(X, y, markers) {
  rf <- ranger::ranger(x = X, y = y, probability = TRUE, num.trees = 500,
                       importance = "permutation")
  v <- rf$variable.importance[markers]; v[is.na(v)] <- 0; v
}

#' sPLS-DA loadings importance from step03 payload.
tw_get_splsda_importance <- function(payload03, scenario_id, markers) {
  sc <- payload03$scenarios[[scenario_id]]
  if (is.null(sc) || is.null(sc$drivers)) return(NULL)
  dr <- sc$drivers
  if (!all(c("Marker", "Importance") %in% names(dr))) return(NULL)
  v <- setNames(rep(0, length(markers)), markers)
  common <- intersect(dr$Marker, markers)
  v[common] <- abs(dr$Importance[match(common, dr$Marker)])
  v
}

#' Network node degree from step04 payload (case topology).
tw_get_network_degree <- function(payload04, scenario_id, markers) {
  sc <- payload04$scenarios[[scenario_id]]
  topo <- if (!is.null(sc$topo_case)) sc$topo_case else sc$topo_ctrl
  if (is.null(topo) || !all(c("Node", "Degree") %in% names(topo))) return(NULL)
  v <- setNames(rep(0, length(markers)), markers)
  common <- intersect(topo$Node, markers)
  v[common] <- topo$Degree[match(common, topo$Node)]
  v
}

#' Assemble the importance sources into an agreement audit. Returns a pairwise
#' Spearman-rank-correlation matrix and a top-k Jaccard matrix across the
#' available sources. Honest either way: low agreement is itself a finding.
tw_explanation_agreement <- function(sources, top_k = 10) {
  sources <- sources[!vapply(sources, is.null, logical(1))]
  nm <- names(sources); m <- length(nm)
  rho <- jac <- matrix(NA_real_, m, m, dimnames = list(nm, nm))
  for (i in seq_len(m)) for (j in seq_len(m)) {
    a <- sources[[i]]; b <- sources[[j]]
    common <- intersect(names(a), names(b))
    rho[i, j] <- suppressWarnings(stats::cor(a[common], b[common], method = "spearman"))
    ta <- names(sort(a, decreasing = TRUE))[seq_len(min(top_k, length(a)))]
    tb <- names(sort(b, decreasing = TRUE))[seq_len(min(top_k, length(b)))]
    jac[i, j] <- length(intersect(ta, tb)) / length(union(ta, tb))
  }
  off <- upper.tri(rho)
  list(rho = rho, jaccard = jac, sources = nm,
       mean_rho = mean(rho[off], na.rm = TRUE),
       mean_jaccard = mean(jac[off], na.rm = TRUE))
}

#' Heatmap of an agreement matrix (Spearman rho or Jaccard).
tw_plot_agreement <- function(mat, title, fill_lab) {
  df <- expand.grid(Source1 = rownames(mat), Source2 = colnames(mat),
                    stringsAsFactors = FALSE)
  df$value <- as.vector(mat)
  ggplot2::ggplot(df, ggplot2::aes(Source1, Source2, fill = value)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", value)), size = 3) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7",
                                  high = "#B2182B", midpoint = 0.5,
                                  limits = c(min(-1, min(mat, na.rm = TRUE)), 1),
                                  name = fill_lab) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

# --- clinical-data gate -------------------------------------------------------

#' Detect whether the metadata carries the optional clinical fields that the
#' DATA-GATED pillars (P4/P5/P6) require. Returns a named logical + the
#' resolved column names so the orchestrator can guard each pillar.
tw_detect_clinical_meta <- function(metadata, cfg = list()) {
  has_col <- function(opts) {
    opts <- opts[!vapply(opts, is.null, logical(1))]
    hit <- intersect(opts, colnames(metadata))
    if (length(hit) > 0) hit[1] else NA_character_
  }
  time_col  <- has_col(c(cfg$survival$time_col, "OS_time", "PFS_time", "time", "survival_time"))
  event_col <- has_col(c(cfg$survival$event_col, "OS_event", "event", "status"))
  batch_col <- has_col(c(cfg$confound$batch_col, "Batch", "batch", "acquisition_date", "Acquisition_Date"))
  sex_col   <- has_col(c("Sex", "sex", "Gender"))
  age_col   <- has_col(c("Age", "age"))
  list(
    survival = !is.na(time_col) && !is.na(event_col),
    confound = !is.na(batch_col),
    fairness = !is.na(sex_col) || !is.na(age_col),
    cols = list(time = time_col, event = event_col, batch = batch_col,
                sex = sex_col, age = age_col))
}

# --- P4/P5/P6 ---
# The data-gated pillars (confound / survival / fairness) are detected by
# tw_detect_clinical_meta() above and reported as PENDING by src/06_trustworthy_ai.R.
# No implementation exists yet: the previous stubs here were never called from anywhere
# and stop()ped if they ever were, so they have been removed rather than left as a
# misleading suggestion that the pillars are wired up. Implement them here when the
# clinical metadata (batch / time-to-event / demographics) actually lands.
