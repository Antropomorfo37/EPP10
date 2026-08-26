# =============================================================================
# bootstrap_stability.R
# -----------------------------------------------------------------------------
# Two-level stability:
#   (1) CLASSIFICATION-STAGE BOOTSTRAP (primary, B=2000, fast)
#       — scores fixed, resample reference cohort to perturb z-score distribution,
#         re-apply classifier. Captures classification-level uncertainty.
#   (2) PIPELINE-STAGE BOOTSTRAP (sensitivity, B=50, parallel, cache resumable)
#       — full re-fit PACE → mFACEs → classification on resampled pseudo-subjects.
#         Captures variance of the mFACEs stage.
# =============================================================================

.libPaths(c("~/.R/library", .libPaths()))

# --- Localiza la raíz del repositorio (Rscript, source() o RStudio) ---------
if (!exists("epp10_path")) {
  .epp10_self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  source(file.path(if (length(.epp10_self)) dirname(sub("^--file=", "", .epp10_self[[1]]))
                   else getwd(), "epp10_paths.R"))
}
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(future); library(future.apply)
})
source(epp10_path("simulate_pseudo_ipd.R"), echo = FALSE)
source(epp10_path("mfaces_dryrun.R"), echo = FALSE)

SEED <- 20260422
REFERENCE <- "no_obese_without_T2DM"
K_CLS <- 3L

# =============================================================================
# PART 1 — CLASSIFICATION-STAGE BOOTSTRAP (B=2000)
# =============================================================================
cat("=== PART 1: Classification-stage bootstrap (B=2000) ===\n\n")

primary <- readRDS(epp10_path("fit_mfaces_primary_results.rds"))
scores_primary   <- primary$retained_primary$mfpca$scores
cohort_vec       <- primary$cohort_vec
K_retained       <- ncol(scores_primary)

# Incretin axis from PRIMARY fit (fixed across bootstrap)
incretin_avail <- intersect(c("GIP_total","GIP_active","GLP1_total","GLP1_active",
                               "PYY_total","PYY_3_36"),
                            names(primary$retained_primary$mfpca$functions))
all_load <- attr(identify_incretin_axis(primary$retained_primary,
                                        incretin = incretin_avail),
                 "all_loadings")
loads_first_K <- all_load[seq_len(min(K_CLS, length(all_load)))]
INC_AXIS <- which.max(loads_first_K)
cat(sprintf("Fixed incretin axis: PC%d  (loading=%.3f)\n\n",
            INC_AXIS, loads_first_K[INC_AXIS]))

classify_stage_boot <- function(scores, cohort, reference, incretin_axis,
                                B, seed, K_cls = 3L) {
  set.seed(seed)
  ref_idx <- which(cohort == reference)
  N <- nrow(scores); n_ref <- length(ref_idx)
  K <- min(K_cls, ncol(scores))

  class_matrix <- matrix(NA_character_, nrow = N, ncol = B)
  for (b in seq_len(B)) {
    # Bootstrap reference indices → new μ_ref, σ_ref
    boot_ref <- sample(ref_idx, n_ref, replace = TRUE)
    mu_b <- colMeans(scores[boot_ref, , drop = FALSE], na.rm = TRUE)
    sd_b <- apply(scores[boot_ref, , drop = FALSE], 2, sd, na.rm = TRUE)
    sd_b[sd_b == 0 | is.na(sd_b)] <- 1
    z_b <- sweep(sweep(scores, 2, mu_b, "-"), 2, sd_b, "/")
    z_b3 <- z_b[, seq_len(K), drop = FALSE]
    class_matrix[, b] <- as.character(classify_by_scores(z_b3,
                                                         incretin_axis = incretin_axis))
  }

  # Modal class per subject + stability
  modal_tbl <- map_dfr(seq_len(N), function(i) {
    tbl <- table(class_matrix[i, ])
    modal <- names(tbl)[which.max(tbl)]
    stab  <- as.numeric(max(tbl)) / B
    tibble(subject_id = i, modal_class = modal, stability = stab)
  })
  attr(modal_tbl, "B") <- B
  attr(modal_tbl, "incretin_axis") <- incretin_axis
  modal_tbl
}

t0 <- Sys.time()
stab_class <- classify_stage_boot(
  scores = scores_primary,
  cohort = cohort_vec,
  reference = REFERENCE,
  incretin_axis = INC_AXIS,
  B = 2000,
  seed = SEED,
  K_cls = K_CLS
)
elapsed_class <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("Classification-stage B=2000 completed in %.1f s\n\n", elapsed_class))

# Aggregate
stab_class$cohort <- cohort_vec
prop_stable <- mean(stab_class$stability >= 0.80, na.rm = TRUE)
cat(sprintf("Proportion of subjects with stability ≥ 0.80 (target): %.1f%%\n",
            100 * prop_stable))

cat("\nStability distribution — summary:\n")
print(summary(stab_class$stability))

cat("\nStability by cohort:\n")
cohort_stab <- stab_class %>% group_by(cohort) %>%
  summarise(n = n(),
            median_stab = median(stability),
            q1_stab = quantile(stability, 0.25),
            q3_stab = quantile(stability, 0.75),
            pct_stable_80 = mean(stability >= 0.80),
            .groups = "drop") %>% arrange(desc(median_stab))
print(cohort_stab)

cat("\nModal class × primary class agreement:\n")
primary_class <- as.character(primary$classification)
agree <- mean(stab_class$modal_class == primary_class)
cat(sprintf("  Modal-class = primary-class: %.1f%%\n", 100 * agree))

# Save
saveRDS(stab_class, epp10_path("stability_classification_stage.rds"))
write_csv(stab_class %>% select(-subject_id),
          epp10_path("stability_classification_stage.csv"))

# =============================================================================
# PART 2 — PIPELINE-STAGE BOOTSTRAP (B=50, full re-fits)
# =============================================================================
cat("\n\n=== PART 2: Pipeline-stage bootstrap (B=50 full re-fits) ===\n")
cat("Parallel via future::multisession(workers=8)\n\n")

CACHE_DIR <- epp10_path("cache_bootstrap_pipeline")
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

B_PIPE <- 50
summary_long <- read_csv(epp10_path("hormones_long_tidy.csv"),
                         show_col_types = FALSE)

# Per-rep seeded function — cache resumable
run_pipeline_rep <- function(b) {
  cache_file <- file.path(CACHE_DIR, sprintf("rep_%04d.rds", b))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  set.seed(SEED + b)
  pipd_full <- simulate_pseudo_ipd(summary_long, M = 1000, rho = 0.5,
                                   cv_mult = 1.0, seed = SEED + b)
  subj_keep <- pipd_full %>% distinct(subject_id, Author, cohort) %>%
    group_by(Author, cohort) %>% slice_sample(n = 50) %>%
    ungroup() %>% pull(subject_id)
  pipd_sub <- pipd_full %>% filter(subject_id %in% subj_keep)

  mfaces <- tryCatch(
    fit_mfaces_joint(pipd_sub,
                     analytes = c("ghrelin_total","ghrelin_acyl","GIP_total",
                                  "GIP_active","GLP1_total","GLP1_active",
                                  "PYY_total","PYY_3_36","glucagon","insulin",
                                  "glucose"),
                     reference_cohort = REFERENCE),
    error = function(e) NULL)
  if (is.null(mfaces)) {
    res <- list(b = b, status = "failed", err = "mfaces error")
    saveRDS(res, cache_file); return(res)
  }
  sens <- retain_by_fve_sensitivity(mfaces, thresholds = 0.90,
                                    n_subjects = nrow(mfaces$scores))
  retained_primary <- sens$by_threshold$fve_0.90
  cohort_v <- tibble(subject_id = rownames(retained_primary$mfpca$scores)) |>
    left_join(distinct(pipd_sub, subject_id, cohort), by = "subject_id") |>
    pull(cohort)
  z_full <- zscore_vs_reference(retained_primary$mfpca$scores, cohort_v,
                                reference = REFERENCE)
  z3 <- z_full[, seq_len(min(K_CLS, ncol(z_full))), drop = FALSE]
  all_load_b <- attr(identify_incretin_axis(retained_primary,
                                            incretin = incretin_avail),
                     "all_loadings")
  loads_first_K_b <- all_load_b[seq_len(min(K_CLS, length(all_load_b)))]
  inc_axis_b <- which.max(loads_first_K_b)
  cls_b <- classify_by_scores(z3, incretin_axis = inc_axis_b)

  # Cohort-level class distribution (ecological, stable across reps)
  dist_by_cohort <- tibble(cohort = cohort_v, cls = as.character(cls_b)) %>%
    count(cohort, cls, name = "n") %>%
    group_by(cohort) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
    select(-n)

  res <- list(b = b, status = "ok",
              K_retained = retained_primary$diagnostics$K_retained,
              incretin_axis = inc_axis_b,
              incretin_loading = loads_first_K_b[inc_axis_b],
              dist_by_cohort = dist_by_cohort)
  saveRDS(res, cache_file)
  res
}

# Determine which reps still need running
cached <- file.exists(file.path(CACHE_DIR,
                                sprintf("rep_%04d.rds", seq_len(B_PIPE))))
to_run <- seq_len(B_PIPE)[!cached]
cat(sprintf("Cached: %d / %d | To run: %d\n", sum(cached), B_PIPE, length(to_run)))

if (length(to_run) > 0) {
  plan(multisession, workers = 8)
  t0 <- Sys.time()
  future_lapply(to_run, run_pipeline_rep, future.seed = TRUE)
  elapsed_pipe <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("Pipeline B=%d completed in %.1f min\n",
              length(to_run), elapsed_pipe / 60))
  plan(sequential)
}

# Aggregate
all_reps <- map(seq_len(B_PIPE),
                ~ readRDS(file.path(CACHE_DIR, sprintf("rep_%04d.rds", .x))))
n_ok <- sum(map_chr(all_reps, "status") == "ok")
cat(sprintf("\nSuccessful reps: %d / %d\n", n_ok, B_PIPE))

# Distribution stability per (cohort, class) across reps
dist_all <- map_dfr(all_reps, function(r) {
  if (r$status != "ok") return(tibble())
  r$dist_by_cohort %>% mutate(b = r$b)
})
dist_summary <- dist_all %>%
  group_by(cohort, cls) %>%
  summarise(median_pct = median(pct),
            q1_pct = quantile(pct, 0.25),
            q3_pct = quantile(pct, 0.75),
            min_pct = min(pct), max_pct = max(pct),
            .groups = "drop") %>%
  arrange(cohort, desc(median_pct))
cat("\nCohort × class prevalence stability (median [IQR] across B=50 reps):\n")
print(dist_summary, n = Inf)

# K retained and incretin loading stability
pipeline_summary <- tibble(
  b = map_dbl(all_reps, ~ .x$b %||% NA_real_),
  K_retained = map_dbl(all_reps, ~ .x$K_retained %||% NA_real_),
  incretin_axis = map_int(all_reps, ~ as.integer(.x$incretin_axis %||% NA)),
  incretin_loading = map_dbl(all_reps, ~ .x$incretin_loading %||% NA_real_)
) %>% filter(!is.na(K_retained))
cat("\nK_retained and incretin-loading stability:\n")
print(summary(pipeline_summary$K_retained))
print(summary(pipeline_summary$incretin_loading))

saveRDS(list(stab_class = stab_class, cohort_stab = cohort_stab,
             pipeline_summary = pipeline_summary,
             dist_summary = dist_summary),
        epp10_path("bootstrap_stability_results.rds"))
cat("\nSaved: bootstrap_stability_results.rds\n")
cat("Cache directory:", CACHE_DIR, "\n")
