# =============================================================================
# bootstrap_B2000_zenodo.R — B=2000 pipeline re-fits for Zenodo archive
# =============================================================================
# Full pipeline re-fits (PACE → mFACEs → classification) for B=2000 pseudo-IPD
# replicates. Cache resumable — every rep written to cache/B2000/rep_NNNN.rds
# immediately, so partial progress survives interruptions.
#
# Expected total compute: ~25 hours with 8 workers, ~200 hours sequential.
# Designed to run in session-persistent background. Re-invoke to resume.
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
B_TARGET <- 2000
CACHE_DIR <- epp10_path("cache_bootstrap_B2000")
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
INCRETIN <- c("GIP_total","GIP_active","GLP1_total","GLP1_active",
              "PYY_total","PYY_3_36")

summary_long <- read_csv(epp10_path("hormones_long_tidy.csv"),
                         show_col_types = FALSE)

run_rep <- function(b) {
  cache_file <- file.path(CACHE_DIR, sprintf("rep_%05d.rds", b))
  if (file.exists(cache_file)) return(b)
  t_start <- Sys.time()
  set.seed(SEED + b)
  pipd_full <- simulate_pseudo_ipd(summary_long, M = 1000, rho = 0.5,
                                   cv_mult = 1.0, seed = SEED + b)
  subj_keep <- pipd_full %>% distinct(subject_id, Author, cohort) %>%
    group_by(Author, cohort) %>% slice_sample(n = 50) %>%
    ungroup() %>% pull(subject_id)
  pipd_sub <- pipd_full %>% filter(subject_id %in% subj_keep)

  res <- tryCatch({
    mfaces <- fit_mfaces_joint(pipd_sub,
        analytes = c("ghrelin_total","ghrelin_acyl","GIP_total","GIP_active",
                     "GLP1_total","GLP1_active","PYY_total","PYY_3_36",
                     "glucagon","insulin","glucose"),
        reference_cohort = REFERENCE)
    sens <- retain_by_fve_sensitivity(mfaces, thresholds = 0.90,
                                      n_subjects = nrow(mfaces$scores))
    retained_primary <- sens$by_threshold$fve_0.90
    cohort_v <- tibble(subject_id = rownames(retained_primary$mfpca$scores)) |>
      left_join(distinct(pipd_sub, subject_id, cohort), by = "subject_id") |>
      pull(cohort)
    z_full <- zscore_vs_reference(retained_primary$mfpca$scores, cohort_v,
                                  reference = REFERENCE)
    z3 <- z_full[, seq_len(min(K_CLS, ncol(z_full))), drop = FALSE]
    all_load <- attr(identify_incretin_axis(retained_primary, incretin = INCRETIN),
                     "all_loadings")
    loads_first_K <- all_load[seq_len(min(K_CLS, length(all_load)))]
    inc_axis <- which.max(loads_first_K)
    cls <- classify_by_scores(z3, incretin_axis = inc_axis)

    list(b = b, status = "ok",
         K_retained = retained_primary$diagnostics$K_retained,
         incretin_axis = inc_axis,
         incretin_loading = loads_first_K[inc_axis],
         per_subject = tibble(pseudo_subject_id = rownames(retained_primary$mfpca$scores),
                              cohort = cohort_v,
                              class = as.character(cls)),
         dist_by_cohort = tibble(cohort = cohort_v,
                                 cls = as.character(cls)) %>%
           count(cohort, cls, name = "n") %>%
           group_by(cohort) %>% mutate(pct = 100 * n / sum(n)) %>%
           ungroup() %>% select(-n),
         elapsed_sec = as.numeric(difftime(Sys.time(), t_start, units = "secs")))
  }, error = function(e) {
    list(b = b, status = "failed", err = conditionMessage(e),
         elapsed_sec = as.numeric(difftime(Sys.time(), t_start, units = "secs")))
  })
  saveRDS(res, cache_file)
  b
}

# Resume: find which reps still need running
cached <- file.exists(file.path(CACHE_DIR,
                                sprintf("rep_%05d.rds", seq_len(B_TARGET))))
to_run <- seq_len(B_TARGET)[!cached]
cat(sprintf("[B2000-Zenodo] Cached: %d / %d | To run: %d\n",
            sum(cached), B_TARGET, length(to_run)))
cat(sprintf("[B2000-Zenodo] Cache dir: %s\n", CACHE_DIR))

if (length(to_run) == 0) {
  cat("[B2000-Zenodo] All reps completed.\n")
} else {
  plan(multisession, workers = 8)
  t0 <- Sys.time()

  # Emit one notification line per 50 completed reps for monitoring
  chunks <- split(to_run, ceiling(seq_along(to_run) / 50))
  for (chunk in chunks) {
    future_lapply(chunk, run_rep, future.seed = TRUE)
    n_now <- sum(file.exists(file.path(CACHE_DIR,
                                        sprintf("rep_%05d.rds",
                                                seq_len(B_TARGET)))))
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("[B2000-Zenodo] progress %d/%d | %.1f min elapsed\n",
                n_now, B_TARGET, el))
  }
  plan(sequential)

  total_elapsed <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
  cat(sprintf("[B2000-Zenodo] Total elapsed: %.2f hours\n", total_elapsed))
}

cat("[B2000-Zenodo] DONE — aggregating results for Zenodo deposit…\n")

# Aggregate
all_reps <- map(seq_len(B_TARGET), function(b) {
  p <- file.path(CACHE_DIR, sprintf("rep_%05d.rds", b))
  if (file.exists(p)) readRDS(p) else NULL
})
n_ok <- sum(map_chr(all_reps, ~ .x$status %||% "missing") == "ok")
n_fail <- sum(map_chr(all_reps, ~ .x$status %||% "missing") == "failed")
cat(sprintf("[B2000-Zenodo] ok=%d | failed=%d | missing=%d / %d\n",
            n_ok, n_fail, B_TARGET - n_ok - n_fail, B_TARGET))

# Pseudo-subject stability across all reps (where pseudo_subject_ids match)
all_per_subject <- map_dfr(all_reps, function(r) {
  if (is.null(r) || r$status != "ok") return(tibble())
  r$per_subject %>% mutate(b = r$b)
})
subj_modal <- all_per_subject %>%
  group_by(pseudo_subject_id, cohort, class) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(pseudo_subject_id, cohort) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  mutate(stability = n / sum(n)) %>% ungroup()
# Note: pseudo_subject_ids regenerate per rep (different seed), so subj-level
# stability is only meaningful per (Author × cohort × rep_idx) keyed draws.
# Use distribution-level stability instead.

# Cohort × class percentage stability across B reps
dist_stability <- map_dfr(all_reps, function(r) {
  if (is.null(r) || r$status != "ok") return(tibble())
  r$dist_by_cohort %>% mutate(b = r$b)
}) %>%
  group_by(cohort, cls) %>%
  summarise(median_pct = median(pct),
            q1_pct = quantile(pct, 0.25),
            q3_pct = quantile(pct, 0.75),
            min_pct = min(pct), max_pct = max(pct),
            n_reps = n(),
            .groups = "drop") %>%
  arrange(cohort, desc(median_pct))

# K and incretin loading stability
pipeline_meta <- tibble(
  b = map_dbl(all_reps, ~ .x$b %||% NA_real_),
  K_retained = map_dbl(all_reps, ~ .x$K_retained %||% NA_real_),
  incretin_axis = map_int(all_reps, ~ as.integer(.x$incretin_axis %||% NA)),
  incretin_loading = map_dbl(all_reps, ~ .x$incretin_loading %||% NA_real_),
  elapsed_sec = map_dbl(all_reps, ~ .x$elapsed_sec %||% NA_real_)
) %>% filter(!is.na(K_retained))

saveRDS(list(dist_stability = dist_stability,
             pipeline_meta = pipeline_meta,
             n_ok = n_ok, n_fail = n_fail, B_TARGET = B_TARGET,
             seed = SEED),
        epp10_path("bootstrap_B2000_results.rds"))
write_csv(dist_stability,
          epp10_path("bootstrap_B2000_dist_stability.csv"))
cat("[B2000-Zenodo] Saved artefacts for Zenodo deposit.\n")
