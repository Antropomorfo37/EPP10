# =============================================================================
# sensitivity_cv.R — CV-multiplier sensitivity for Verification Appendix
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
})
source(epp10_path("simulate_pseudo_ipd.R"), echo = FALSE)
source(epp10_path("mfaces_dryrun.R"), echo = FALSE)

SEED <- 20260422; SUBSAMPLE_N <- 50; RHO <- 0.5
CV_GRID <- c(0.75, 1.25)              # cv_mult=1.0 already in primary
ANALYTES <- c("ghrelin_total","ghrelin_acyl","GIP_total","GIP_active",
              "GLP1_total","GLP1_active","PYY_total","PYY_3_36",
              "glucagon","insulin","glucose")

summary_long <- read_csv(epp10_path("hormones_long_tidy.csv"),
                         show_col_types = FALSE)

run_cv <- function(cv_mult) {
  cat(sprintf("\n=== cv_mult = %.2f ===\n", cv_mult))
  t0 <- Sys.time()
  pipd_full <- simulate_pseudo_ipd(summary_long, M = 1000, rho = RHO,
                                   cv_mult = cv_mult, seed = SEED)
  set.seed(SEED)
  subj_keep <- pipd_full %>% distinct(subject_id, Author, cohort) %>%
    group_by(Author, cohort) %>% slice_sample(n = SUBSAMPLE_N) %>%
    ungroup() %>% pull(subject_id)
  pipd_sub <- pipd_full %>% filter(subject_id %in% subj_keep)

  mfaces <- fit_mfaces_joint(pipd_sub, analytes = ANALYTES,
                             reference_cohort = "no_obese_without_T2DM")
  sens <- retain_by_fve_sensitivity(mfaces, thresholds = c(0.90, 0.95),
                                    n_subjects = nrow(mfaces$scores))
  retained_primary <- sens$by_threshold$fve_0.90
  K <- retained_primary$diagnostics$K_retained
  fve5 <- cumsum(mfaces$values[1:5]) / sum(mfaces$values)

  cohort_vec <- tibble(subject_id = rownames(retained_primary$mfpca$scores)) |>
    left_join(distinct(pipd_sub, subject_id, cohort), by = "subject_id") |>
    pull(cohort)
  z_full <- zscore_vs_reference(retained_primary$mfpca$scores, cohort_vec,
                                reference = "no_obese_without_T2DM")
  z3 <- z_full[, seq_len(min(3, ncol(z_full))), drop = FALSE]
  all_load <- attr(identify_incretin_axis(retained_primary,
                                          incretin = c("GIP_total","GIP_active",
                                                       "GLP1_total","GLP1_active",
                                                       "PYY_total","PYY_3_36")),
                   "all_loadings")
  loads_first_3 <- all_load[seq_len(min(3, length(all_load)))]
  inc_axis <- which.max(loads_first_3)
  cls <- classify_by_scores(z3, incretin_axis = inc_axis)

  prevalence <- tibble(cohort = cohort_vec, cls = cls) |>
    group_by(cohort) |> count(cls, name = "n") |>
    mutate(pct = round(100 * n / sum(n), 1)) |>
    select(-n) |>
    pivot_wider(names_from = cls, values_from = pct, values_fill = 0)

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  cv_mult=%.2f | K@0.90=%d | incretin PC%d loading=%.3f | %.1fs\n",
              cv_mult, K, inc_axis, loads_first_3[inc_axis], elapsed))

  list(cv_mult = cv_mult, K_retained = K, fve_top5 = fve5,
       incretin_axis = inc_axis,
       incretin_loading = loads_first_3[inc_axis],
       prevalence = prevalence,
       cls_vec = cls, cohort_vec = cohort_vec,
       scores_primary = retained_primary$mfpca$scores,
       elapsed_sec = elapsed)
}

results <- map(CV_GRID, run_cv)
names(results) <- sprintf("cv_%.2f", CV_GRID)

# Attach primary cv_mult=1.0
primary <- readRDS(epp10_path("fit_mfaces_primary_results.rds"))
primary_fve5 <- cumsum(primary$mfaces$values[1:5]) / sum(primary$mfaces$values)
results$cv_1.00 <- list(
  cv_mult = 1.0, K_retained = primary$retained_primary$diagnostics$K_retained,
  fve_top5 = primary_fve5,
  incretin_axis = as.integer(primary$incretin_axis),
  incretin_loading = attr(primary$incretin_axis, "loading"),
  prevalence = primary$prevalence_pct,
  cls_vec = primary$classification,
  cohort_vec = primary$cohort_vec,
  scores_primary = primary$retained_primary$mfpca$scores,
  elapsed_sec = NA
)

cat("\n=== AGGREGATE SUMMARY (3 CV values) ===\n")
summary_tbl <- map_dfr(results, function(r) {
  tibble(cv_mult = r$cv_mult, K_retained = r$K_retained,
         FVE_top5 = r$fve_top5[5],
         incretin_PC = r$incretin_axis,
         incretin_loading = r$incretin_loading)
}) %>% arrange(cv_mult)
print(summary_tbl)

cat("\n=== PREVALENCE BY CV ===\n")
for (tag in names(results)) {
  cat(sprintf("\n-- %s --\n", tag))
  print(results[[tag]]$prevalence)
}

saveRDS(results, epp10_path("sensitivity_cv_results.rds"))
cat("\nSaved: sensitivity_cv_results.rds\n")
