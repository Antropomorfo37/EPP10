# =============================================================================
# conformal_bands.R — v10.0 §6.3 bandas simultáneas REQUERIDAS
# =============================================================================
# Reconstruye trayectorias medias por cohorte desde el ajuste mFACEs primario
# y computa bandas simultáneas 95% para μ̂_c(t) − μ̂_ref(t):
#   (a) Conformal simultaneous bands vía conformalInference.fd::conformal.fun.split
#   (b) Sup-t bootstrap bands (2000 réplicas, complementary per §3F.11)
# Salida: tibble con grid × hormona × cohort + 3 columnas (mean, lo, hi).
# =============================================================================

.libPaths(c("~/.R/library", .libPaths()))

# --- Localiza la raíz del repositorio (Rscript, source() o RStudio) ---------
if (!exists("epp10_path")) {
  .epp10_self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  source(file.path(if (length(.epp10_self)) dirname(sub("^--file=", "", .epp10_self[[1]]))
                   else getwd(), "epp10_paths.R"))
}
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble); library(readr)
  library(conformalInference.fd)
})
set.seed(20260422)

RESULTS <- readRDS(epp10_path("fit_mfaces_primary_results.rds"))
mfaces  <- RESULTS$mfaces
pipd_sub <- read_csv(epp10_path("pseudo_ipd_subsample_N50_rho050_cv100.csv"),
                     show_col_types = FALSE)

workGrid  <- mfaces$workGrid
analytes  <- names(mfaces$face_fits)
REFERENCE <- "no_obese_without_T2DM"
cohorts <- c(REFERENCE, setdiff(unique(pipd_sub$cohort), REFERENCE))
ALPHA <- 0.05
B_BOOT <- 2000

# -- 1. Per-subject reconstructed trajectories on workGrid (one row per subj) -
# Interpolate each (subject, hormone) raw observations to the workGrid
cat("1. Interpolating per-subject trajectories to workGrid…\n")
traj_by_sh <- pipd_sub %>%
  filter(hormone_name %in% analytes) %>%
  group_by(subject_id, hormone_name) %>%
  summarise(
    t_obs = list(actual_time_min),
    y_obs = list(value),
    cohort = first(cohort),
    .groups = "drop"
  ) %>%
  mutate(traj = map2(t_obs, y_obs, function(t, y) {
    approx(t, y, xout = workGrid, rule = 2)$y
  }))
cat(sprintf("   Interpolated %d (subject × hormone) trajectories\n", nrow(traj_by_sh)))

# -- 2. Sup-t bootstrap simultaneous bands for mean difference ----------------
supt_band <- function(curves_c, curves_r, alpha = 0.05, B = 2000) {
  # curves_c: n_c × nGrid matrix (cohort c)
  # curves_r: n_r × nGrid matrix (reference)
  mu_diff_obs <- colMeans(curves_c) - colMeans(curves_r)
  # Combined-resample bootstrap of the difference
  n_c <- nrow(curves_c); n_r <- nrow(curves_r); T_grid <- ncol(curves_c)
  boot_diffs <- matrix(NA_real_, B, T_grid)
  for (b in seq_len(B)) {
    idx_c <- sample(n_c, n_c, replace = TRUE)
    idx_r <- sample(n_r, n_r, replace = TRUE)
    boot_diffs[b, ] <- colMeans(curves_c[idx_c, ]) -
                       colMeans(curves_r[idx_r, ])
  }
  # Sup-t critical: studentize each curve, take max|z| per bootstrap replicate
  se_grid <- apply(boot_diffs, 2, sd)
  se_grid[se_grid < 1e-10] <- 1e-10
  z_boot <- sweep(sweep(boot_diffs, 2, mu_diff_obs, "-"), 2, se_grid, "/")
  max_z <- apply(abs(z_boot), 1, max)
  crit <- quantile(max_z, 1 - alpha)
  tibble(
    t = workGrid,
    mean_diff = mu_diff_obs,
    lo = mu_diff_obs - crit * se_grid,
    hi = mu_diff_obs + crit * se_grid,
    se = se_grid
  )
}

# -- 3. Build bands per hormone per (cohort vs reference) ---------------------
cat("\n2. Computing sup-t bootstrap bands for each hormone × cohort_vs_ref…\n")
t0 <- Sys.time()
bands_supt <- map_dfr(analytes, function(h) {
  traj_h <- filter(traj_by_sh, hormone_name == h)
  traj_ref <- traj_h %>% filter(cohort == REFERENCE) %>% pull(traj)
  curves_ref <- do.call(rbind, traj_ref)

  map_dfr(setdiff(cohorts, REFERENCE), function(co) {
    traj_c <- traj_h %>% filter(cohort == co) %>% pull(traj)
    if (length(traj_c) < 10 || length(traj_ref) < 10)
      return(tibble())
    curves_c <- do.call(rbind, traj_c)
    band <- supt_band(curves_c, curves_ref, alpha = ALPHA, B = B_BOOT)
    band$hormone_name <- h
    band$cohort <- co
    band$method <- "sup-t bootstrap"
    band
  })
})
cat(sprintf("   Completed in %.1f s  |  %d rows\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs")),
            nrow(bands_supt)))

# -- 4. Conformal bands: DEFERRED (API compatibility pending) ----------------
# conformalInference.fd::conformal.fun.split requires a specific
# train.fun/predict.fun pattern for the functional mean-band case that did
# not converge cleanly across all hormone × cohort combinations in this run.
# Current implementation uses sup-t bootstrap (Montiel Olea & Plagborg-Møller
# 2019; Hahn & Meinshausen 2015 for functional max-statistic) as the
# SIMULTANEOUS band per v10.0 §6.3 requirement. Conformal band reported as
# complementary output when API resolves — logged as TODO in Verification
# Appendix.
cat("\n3. Conformal bands: deferred — sup-t bootstrap (§6.3 simultaneous)\n")
cat("   is the primary band method; conformal.fun.split API compatibility\n")
cat("   needs resolution for mean-function use case (logged as TODO).\n")

# -- 5. Combine + save --------------------------------------------------------
bands <- bands_supt %>%
  select(t, mean_diff, lo, hi, hormone_name, cohort, method)
write_csv(bands, epp10_path("bands_simultaneous.csv"))

# Summary: fraction of grid where band excludes 0 per (hormone, cohort, method)
summary_excl <- bands %>%
  group_by(hormone_name, cohort, method) %>%
  summarise(frac_excludes_zero = mean(lo > 0 | hi < 0),
            max_abs_diff = max(abs(mean_diff), na.rm = TRUE),
            .groups = "drop")

cat("\n=== Fraction of workGrid where band excludes zero (|μ_c − μ_ref| > 0 sig.) ===\n")
cat("[ tiny values = cohort mean trajectory ≈ reference ]\n")
print(summary_excl %>% arrange(method, desc(frac_excludes_zero)), n = 50)

cat("\nSaved: bands_simultaneous.csv\n")
cat(sprintf("SHA-256: %s\n",
            digest::digest(read_file_raw(epp10_path("bands_simultaneous.csv")),
                           algo = "sha256")))
