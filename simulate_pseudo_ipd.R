# =============================================================================
# simulate_pseudo_ipd.R — v10.0 §4.2 pseudo-IPD generator
# =============================================================================
# Generates M synthetic trajectories per (Author, cohort, hormone) from
#   GP(μ_cohort(t), K_ρ(s,t))
# with:
#   • K_ρ(s,t) = ρ^(|s-t| / Δt_ref)   (AR(1) continuous-time, Δt_ref = 30 min)
#   • σ_h(t)   = CV_h × μ_h(t)         (CV from hormone_variability_priors YAML)
# Output schema matches fit_mfaces_joint / fit_all_analytes input requirements.
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
})

# ---- 1. CV priors (frozen per YAML 2026-04-22) -----------------------------
HORMONE_CV <- c(
  ghrelin_total = 0.40, ghrelin_acyl = 0.45,
  GIP_total     = 0.30, GIP_active   = 0.35,
  GLP1_total    = 0.35, GLP1_active  = 0.40,
  PYY_total     = 0.35, PYY_3_36     = 0.40,
  glucagon      = 0.30, insulin      = 0.20, glucose = 0.12
)

# ---- 2. AR(1) kernel -------------------------------------------------------
ar1_kernel <- function(t, rho, delta_ref = 30) {
  d <- as.matrix(dist(t, diag = TRUE, upper = TRUE))
  rho ^ (d / delta_ref)
}

# ---- 3. Multivariate normal sampler (Cholesky) -----------------------------
sample_gp <- function(mu, Sigma, M, seed = NULL,
                      clip_nonneg = TRUE, sigma_floor_frac = 0.01) {
  if (!is.null(seed)) set.seed(seed)
  n_t <- length(mu)
  # Sanitize inputs
  if (any(!is.finite(mu)))     mu[!is.finite(mu)] <- 0
  if (any(!is.finite(Sigma)))  Sigma[!is.finite(Sigma)] <- 0
  # Ridge floor — scaled to matrix content, never below a numerical epsilon
  content_scale <- max(abs(mu), 1, diag(Sigma), na.rm = TRUE)
  diag_floor <- max((sigma_floor_frac * content_scale) ^ 2, 1e-8)
  # Progressive ridge retry (1×, 10×, 100×, 1000×)
  L <- NULL
  for (mult in c(1, 10, 100, 1000)) {
    L <- tryCatch(
      chol(Sigma + diag(diag_floor * mult, n_t)),
      error = function(e) NULL
    )
    if (!is.null(L)) break
  }
  if (is.null(L)) {
    # Final fallback: diagonal-only sampling (treat as independent)
    L <- diag(sqrt(diag(Sigma) + diag_floor), n_t)
  }
  Z <- matrix(rnorm(M * n_t), nrow = n_t, ncol = M)
  Y <- mu + t(L) %*% Z
  if (clip_nonneg) Y[Y < 0] <- 0
  t(Y)
}

# ---- 4. Main function ------------------------------------------------------
# summary_long: tibble from ETL with columns
#   Author, cohort_v10_primary (+other meta), hormone_name,
#   time_min, value, n_subjects
# priors: named numeric vector of CVs, indexed by hormone_name
#
# For each (Author, cohort_v10_primary), we draw M pseudo-subjects,
# where each pseudo-subject has one trajectory per hormone the study reported
# (independent per hormone — no cross-hormone covariance in the master).
#
# Output schema (drop-in for fit_mfaces_joint):
#   subject_id, cohort, hormone_name, timepoint_min, actual_time_min,
#   value, is_censored, pseudo_rho, pseudo_cv_mult,
#   Author, n_subjects_source, ...meta...

simulate_pseudo_ipd <- function(summary_long,
                                priors = HORMONE_CV,
                                M        = 1000,
                                rho      = 0.5,
                                cv_mult  = 1.0,
                                seed     = 20260422,
                                clip_nonneg = TRUE) {
  stopifnot(all(unique(summary_long$hormone_name) %in% names(priors)))

  # Pre-aggregation: collapse duplicate (Author, cohort_v10_primary, hormone, time)
  # entries from source-cohort pooling (e.g., Aukan 2022 Obesity grades I/II/III
  # → cohort_v10_primary = Obesity) using n_subjects-weighted mean.
  summary_long <- summary_long %>%
    dplyr::group_by(Author, cohort_v10_primary, hormone_name, time_min) %>%
    dplyr::summarise(
      value        = stats::weighted.mean(value, w = pmax(n_subjects, 1), na.rm = TRUE),
      n_subjects   = sum(n_subjects, na.rm = TRUE),
      # Take first non-NA of metadata (they should be identical within the group
      # modulo source_cohort, which we keep as concatenation for audit trail)
      source_cohort = paste(sort(unique(source_cohort)), collapse = " | "),
      cohort_v10_sensitivity = first(cohort_v10_sensitivity),
      surgery_status = first(surgery_status),
      weeks_post_surgery = first(weeks_post_surgery),
      weight_loss_modality = first(weight_loss_modality),
      had_t2dm_pre_surgery = first(had_t2dm_pre_surgery),
      .groups = "drop"
    ) %>%
    dplyr::arrange(Author, cohort_v10_primary, hormone_name, time_min)

  groups <- summary_long %>%
    dplyr::distinct(Author, cohort_v10_primary)

  set.seed(seed)
  out <- purrr::pmap_dfr(groups, function(Author, cohort_v10_primary) {
    arm <- dplyr::filter(summary_long,
                         Author == !!Author,
                         cohort_v10_primary == !!cohort_v10_primary)

    # Draw per hormone within this arm; assemble into pseudo-subjects
    per_hormone <- purrr::map_dfr(unique(arm$hormone_name), function(h) {
      rows <- dplyr::filter(arm, hormone_name == h)
      t_obs  <- rows$time_min
      mu_obs <- rows$value
      n_t    <- length(t_obs)
      if (n_t < 1) return(tibble())

      cv_h <- priors[[h]] * cv_mult
      sigma_h <- cv_h * pmax(mu_obs, 0)
      # Build Σ = diag(σ) × C_ρ × diag(σ)
      Crho <- ar1_kernel(t_obs, rho)
      Sigma <- diag(sigma_h, n_t) %*% Crho %*% diag(sigma_h, n_t)

      Y <- sample_gp(mu = mu_obs, Sigma = Sigma, M = M,
                     seed = NULL, clip_nonneg = clip_nonneg)   # M × n_t

      # Pivot to long
      tibble::tibble(
        rep_idx         = rep(seq_len(M), each = n_t),
        timepoint_min   = rep(t_obs, times = M),
        actual_time_min = rep(t_obs, times = M),
        value           = as.vector(t(Y)),
        hormone_name    = h
      )
    })

    if (nrow(per_hormone) == 0) return(tibble())

    # Attach meta: one row's metadata propagates to all hormones of this arm
    meta_row <- arm[1, c("source_cohort","cohort_v10_sensitivity","surgery_status",
                         "weeks_post_surgery","weight_loss_modality",
                         "had_t2dm_pre_surgery","n_subjects")]

    per_hormone %>%
      dplyr::mutate(
        Author            = !!Author,
        cohort            = !!cohort_v10_primary,
        subject_id        = sprintf("%s__%s__%04d",
                                    !!Author, !!cohort_v10_primary, rep_idx),
        is_censored       = FALSE,
        pseudo_rho        = rho,
        pseudo_cv_mult    = cv_mult,
        source_cohort     = meta_row$source_cohort,
        cohort_v10_sensitivity = meta_row$cohort_v10_sensitivity,
        surgery_status    = meta_row$surgery_status,
        weeks_post_surgery = meta_row$weeks_post_surgery,
        weight_loss_modality = meta_row$weight_loss_modality,
        had_t2dm_pre_surgery = meta_row$had_t2dm_pre_surgery,
        n_subjects_source  = meta_row$n_subjects
      ) %>%
      dplyr::select(subject_id, cohort, hormone_name, timepoint_min, actual_time_min,
                    value, is_censored, pseudo_rho, pseudo_cv_mult,
                    Author, n_subjects_source, source_cohort,
                    cohort_v10_sensitivity, surgery_status, weeks_post_surgery,
                    weight_loss_modality, had_t2dm_pre_surgery)
  })

  attr(out, "M")       <- M
  attr(out, "rho")     <- rho
  attr(out, "cv_mult") <- cv_mult
  attr(out, "seed")    <- seed
  attr(out, "priors_used") <- priors * cv_mult
  out
}

# =============================================================================
# Test run on the ETL output (small M for smoke test)
# =============================================================================
if (sys.nframe() == 0L) {
  summary_long <- readr::read_csv(epp10_path("hormones_long_tidy.csv"),
                                  show_col_types = FALSE)
  cat(sprintf("Input summary: %d rows\n", nrow(summary_long)))

  # Smoke test: M=50, rho=0.5, cv_mult=1.0
  pipd <- simulate_pseudo_ipd(summary_long, M = 50, rho = 0.5, cv_mult = 1.0)
  cat(sprintf("\nPseudo-IPD output: %d rows\n", nrow(pipd)))
  cat(sprintf("Unique pseudo-subjects: %d\n", dplyr::n_distinct(pipd$subject_id)))
  cat(sprintf("Cohorts: %d | Hormones: %d\n",
              dplyr::n_distinct(pipd$cohort), dplyr::n_distinct(pipd$hormone_name)))

  cat("\nHead (compact):\n")
  print(dplyr::select(pipd, subject_id, cohort, hormone_name,
                      timepoint_min, value, pseudo_rho) |> head(10))

  cat("\nPseudo-subject counts per cohort:\n")
  print(pipd %>% distinct(subject_id, cohort) %>% count(cohort))

  cat("\nSanity: value statistics by cohort × hormone (first 5 combos):\n")
  print(pipd %>% group_by(cohort, hormone_name) %>%
          summarise(n = n(), mean_val = mean(value),
                    sd_val = sd(value), min_v = min(value),
                    .groups = "drop") %>% head(8))

  # Compare to source means (sanity: pseudo mean should ≈ source mean)
  a1 <- pipd$Author[1]; c1 <- pipd$cohort[1]
  cat(sprintf("\nSource vs pseudo-IPD mean comparison (%s × %s):\n", a1, c1))
  src <- summary_long %>%
    filter(Author == a1, cohort_v10_primary == c1) %>%
    group_by(hormone_name, time_min) %>%
    summarise(source_mean = weighted.mean(value, pmax(n_subjects, 1)),
              .groups = "drop") %>% slice_head(n = 8)
  pi_means <- pipd %>% filter(Author == a1, cohort == c1) %>%
    group_by(hormone_name, timepoint_min) %>%
    summarise(pseudo_mean = mean(value), pseudo_sd = sd(value),
              .groups = "drop") %>%
    rename(time_min = timepoint_min)
  print(left_join(src, pi_means, by = c("hormone_name", "time_min")))

  # Integrity: each subject should have exactly 1 value per (hormone, timepoint)
  dupes <- pipd %>% count(subject_id, hormone_name, timepoint_min) %>% filter(n > 1)
  cat(sprintf("\nDuplicate subject×hormone×time rows (should be 0): %d\n",
              nrow(dupes)))

  cat("\nSHA-256 of generated pseudo-IPD (compact tibble):\n")
  tmpf <- tempfile(fileext = ".csv")
  readr::write_csv(pipd, tmpf)
  cat(digest::digest(readr::read_file_raw(tmpf), algo = "sha256"), "\n")
}
