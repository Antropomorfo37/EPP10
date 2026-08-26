# =============================================================================
# recover_B2000_results.R — se ejecuta cuando el B=2000 Zenodo complete
# =============================================================================
# Agrega bootstrap_B2000_results.rds y bootstrap_B2000_dist_stability.csv
# Comparar pipeline B=50 vs B=2000 para validar estabilidad del hallazgo
# =============================================================================

.libPaths(c("~/.R/library", .libPaths()))

# --- Localiza la raíz del repositorio (Rscript, source() o RStudio) ---------
if (!exists("epp10_path")) {
  .epp10_self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  source(file.path(if (length(.epp10_self)) dirname(sub("^--file=", "", .epp10_self[[1]]))
                   else getwd(), "epp10_paths.R"))
}
suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(tidyr); library(readr)
})

CACHE_DIR <- epp10_path("cache_bootstrap_B2000")
B_TARGET <- 2000

# 1. Verify all reps cached
cached <- list.files(CACHE_DIR, pattern = "^rep_\\d+\\.rds$")
n_cached <- length(cached)
cat(sprintf("Cache directory: %s\n", CACHE_DIR))
cat(sprintf("Reps cached: %d / %d\n", n_cached, B_TARGET))
if (n_cached < B_TARGET) {
  cat(sprintf("  Missing: %d reps — resume by re-running bootstrap_B2000_zenodo.R\n",
              B_TARGET - n_cached))
}

# 2. Aggregate all reps
cat("\nAggregating reps …\n")
all_reps <- map(seq_len(B_TARGET), function(b) {
  p <- file.path(CACHE_DIR, sprintf("rep_%05d.rds", b))
  if (file.exists(p)) readRDS(p) else NULL
})
n_ok <- sum(map_chr(all_reps, ~ .x$status %||% "missing") == "ok")
n_fail <- sum(map_chr(all_reps, ~ .x$status %||% "missing") == "failed")
cat(sprintf("  ok=%d | failed=%d | missing=%d\n", n_ok, n_fail, B_TARGET-n_ok-n_fail))

# 3. Distribution stability
dist_stability <- map_dfr(all_reps, function(r) {
  if (is.null(r) || r$status != "ok") return(tibble())
  r$dist_by_cohort %>% mutate(b = r$b)
}) %>%
  group_by(cohort, cls) %>%
  summarise(median_pct = median(pct),
            q1_pct = quantile(pct, 0.25),
            q3_pct = quantile(pct, 0.75),
            min_pct = min(pct), max_pct = max(pct),
            n_reps = n(), .groups = "drop") %>%
  arrange(cohort, desc(median_pct))

# 4. K and incretin stability
pipeline_meta <- tibble(
  b = map_dbl(all_reps, ~ .x$b %||% NA_real_),
  K = map_dbl(all_reps, ~ .x$K_retained %||% NA_real_),
  inc_axis = map_int(all_reps, ~ as.integer(.x$incretin_axis %||% NA)),
  inc_loading = map_dbl(all_reps, ~ .x$incretin_loading %||% NA_real_),
  elapsed = map_dbl(all_reps, ~ .x$elapsed_sec %||% NA_real_)
) %>% filter(!is.na(K))

cat(sprintf("\nK retained: median=%d  IQR=[%d, %d]  range=[%d, %d]\n",
            median(pipeline_meta$K),
            quantile(pipeline_meta$K, 0.25),
            quantile(pipeline_meta$K, 0.75),
            min(pipeline_meta$K), max(pipeline_meta$K)))
cat(sprintf("Incretin loading: median=%.3f  IQR=[%.3f, %.3f]\n",
            median(pipeline_meta$inc_loading),
            quantile(pipeline_meta$inc_loading, 0.25, na.rm = TRUE),
            quantile(pipeline_meta$inc_loading, 0.75, na.rm = TRUE)))
cat(sprintf("Total compute: %.1f h (%.1f s/rep × 8 workers)\n",
            sum(pipeline_meta$elapsed, na.rm = TRUE) / 3600,
            median(pipeline_meta$elapsed, na.rm = TRUE)))

# 5. Comparison vs B=50 pipeline
b50 <- readRDS(epp10_path("bootstrap_stability_results.rds"))$dist_summary
cmp <- dist_stability %>%
  inner_join(b50 %>% select(cohort, cls, b50_median = median_pct,
                              b50_q1 = q1_pct, b50_q3 = q3_pct),
             by = c("cohort","cls")) %>%
  mutate(median_drift = median_pct - b50_median,
         iqr_ratio = (q3_pct - q1_pct) / pmax(b50_q3 - b50_q1, 0.1))
cat("\n=== B=2000 vs B=50 drift (medians and IQR ratios) ===\n")
cat("Large |drift| or iqr_ratio > 1.5 indicates B=50 under-estimated uncertainty\n\n")
print(cmp %>% slice_max(abs(median_drift), n = 10), width = Inf)

# 6. Save
write_csv(dist_stability, epp10_path("bootstrap_B2000_dist_stability.csv"))
write_csv(cmp, epp10_path("bootstrap_B2000_vs_B50_comparison.csv"))
saveRDS(list(dist_stability = dist_stability, pipeline_meta = pipeline_meta,
             comparison = cmp, n_ok = n_ok, n_fail = n_fail),
        epp10_path("bootstrap_B2000_results.rds"))

cat("\nArtefactos guardados. Actualizar Apéndice de Verificación S3.6.3 con estos números.\n")
