# =============================================================================
# fit_mfaces_realdata.R — pipeline primario sobre pseudo-IPD desde master CSV
# =============================================================================
# 1. Lee hormones_long_tidy.csv (ETL del master)
# 2. Genera pseudo-IPD M=1000, ρ=0.5, cv_mult=1.0 (archivo Zenodo completo)
# 3. Sub-sample determinístico a N=50 por (Author × cohort_v10_primary) arm
# 4. fit_mfaces_joint sobre el sub-sample
# 5. retain_by_fve sensitivity (0.90 primario, 0.95 secundario)
# 6. Clasificación lean-referenced con eje incretínico argmax
# 7. Emite tabla de prevalencia por cohorte + artefactos para appendix
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
  library(stringr); library(digest)
})

# --- Load helpers (sourcing does NOT run dry-run thanks to sys.nframe guard) -
source(epp10_path("simulate_pseudo_ipd.R"), echo = FALSE)
source(epp10_path("mfaces_dryrun.R"), echo = FALSE)

SEED <- 20260422
SUBSAMPLE_N <- 50
PRIMARY_RHO <- 0.5
PRIMARY_CV_MULT <- 1.0

# ---- Step 1. Read ETL output -----------------------------------------------
cat("=== 1. Read ETL output ===\n")
summary_long <- read_csv(epp10_path("hormones_long_tidy.csv"),
                         show_col_types = FALSE)
cat(sprintf("Summary rows: %d | arms: %d | hormones: %d | cohorts: %d\n",
            nrow(summary_long),
            n_distinct(paste(summary_long$Author, summary_long$cohort_v10_primary)),
            n_distinct(summary_long$hormone_name),
            n_distinct(summary_long$cohort_v10_primary)))

# ---- Step 2. Generate pseudo-IPD primary (M=1000, ρ=0.5, cv_mult=1.0) ------
cat("\n=== 2. Generate pseudo-IPD primary (M=1000) ===\n")
t0 <- Sys.time()
pipd_full <- simulate_pseudo_ipd(summary_long,
                                 M = 1000, rho = PRIMARY_RHO,
                                 cv_mult = PRIMARY_CV_MULT, seed = SEED)
cat(sprintf("Generated in %.1f s | rows=%d | pseudo_subjects=%d\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs")),
            nrow(pipd_full), n_distinct(pipd_full$subject_id)))

archive_path <- epp10_path("pseudo_ipd_primary_M1000_rho050_cv100.csv")
write_csv(pipd_full, archive_path)
cat(sprintf("Archive SHA-256: %s\n",
            digest(read_file_raw(archive_path), algo = "sha256")))

# ---- Step 3. Deterministic sub-sample to N=50 per arm ----------------------
cat("\n=== 3. Sub-sample to N=50 per arm ===\n")
set.seed(SEED)
subj_keep <- pipd_full %>%
  distinct(subject_id, Author, cohort) %>%
  group_by(Author, cohort) %>%
  slice_sample(n = SUBSAMPLE_N, replace = FALSE) %>%
  ungroup() %>%
  pull(subject_id)

pipd_sub <- pipd_full %>% filter(subject_id %in% subj_keep)
cat(sprintf("Sub-sampled: %d rows | %d pseudo-subjects across %d arms\n",
            nrow(pipd_sub), n_distinct(pipd_sub$subject_id),
            n_distinct(paste(pipd_sub$Author, pipd_sub$cohort))))
cat(sprintf("Per-cohort pseudo-subject count:\n"))
print(pipd_sub %>% distinct(subject_id, cohort) %>% count(cohort))

subsample_path <- epp10_path("pseudo_ipd_subsample_N50_rho050_cv100.csv")
write_csv(pipd_sub, subsample_path)

# ---- Step 4. fit_mfaces_joint on the sub-sample ----------------------------
cat("\n=== 4. fit_mfaces_joint on real data sub-sample ===\n")
t0 <- Sys.time()
mfaces <- fit_mfaces_joint(
  long     = pipd_sub,
  analytes = c("ghrelin_total","ghrelin_acyl","GIP_total","GIP_active",
               "GLP1_total","GLP1_active","PYY_total","PYY_3_36",
               "glucagon","insulin","glucose"),
  reference_cohort = "no_obese_without_T2DM"
)
cat(sprintf("mFACEs joint fit: %.1f s\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
print(mfaces)

cat("\nHealth check:\n")
print(mfaces_health_check(mfaces))

# ---- Step 5. Retention sensitivity 0.90/0.95 ------------------------------
cat("\n=== 5. retain_by_fve_sensitivity (0.90 primary, 0.95 sensitivity) ===\n")
sens <- retain_by_fve_sensitivity(mfaces, thresholds = c(0.90, 0.95),
                                  n_subjects = nrow(mfaces$scores))
print(sens)

retained_primary <- sens$by_threshold$fve_0.90

# ---- Step 6. Classification ------------------------------------------------
cat("\n=== 6. Lean-referenced classification ===\n")
cohort_vec <- tibble(subject_id = rownames(retained_primary$mfpca$scores)) |>
  left_join(distinct(pipd_sub, subject_id, cohort), by = "subject_id") |>
  pull(cohort)

z_full <- zscore_vs_reference(retained_primary$mfpca$scores, cohort_vec,
                              reference = "no_obese_without_T2DM")
K_retained <- ncol(z_full)

# PRIMARY: classifier spec (YAML: "Los tres primeros scores ξ1, ξ2, ξ3")
# Applied to first 3 retained PCs per pre-registered rule to avoid multiple-
# testing inflation of Altered when K_retained > 3.
K_CLS <- 3L
z <- z_full[, seq_len(min(K_CLS, K_retained)), drop = FALSE]

incretin_avail <- intersect(c("GIP_total","GIP_active","GLP1_total","GLP1_active",
                               "PYY_total","PYY_3_36"),
                            names(retained_primary$mfpca$functions))
# Incretin axis: argmax restricted to first K_CLS retained PCs (consistent
# with the 3-PC classifier spec)
all_loadings <- attr(identify_incretin_axis(retained_primary, incretin = incretin_avail),
                     "all_loadings")
loadings_first_K <- all_loadings[seq_len(min(K_CLS, length(all_loadings)))]
inc_axis <- which.max(loadings_first_K)
cat(sprintf("Incretin axis (restricted to first %d PCs): PC%d  [loading = %.3f]\n",
            K_CLS, inc_axis, loadings_first_K[inc_axis]))
cat(sprintf("First %d PC incretin loadings: %s\n",
            K_CLS, paste(sprintf("%.3f", loadings_first_K), collapse = ", ")))
cat(sprintf("All %d PC incretin loadings: %s\n",
            length(all_loadings),
            paste(sprintf("%.3f", all_loadings), collapse = ", ")))

cls_primary <- classify_by_scores(z, incretin_axis = inc_axis)

# SENSITIVITY: full K_retained classifier
cls_full <- classify_by_scores(z_full,
                               incretin_axis = which.max(all_loadings))

prevalence_primary <- tibble(cohort = cohort_vec, cls_scores = cls_primary) |>
  group_by(cohort) |>
  count(cls_scores, name = "n") |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  select(-n) |>
  pivot_wider(names_from = cls_scores, values_from = pct, values_fill = 0) |>
  arrange(cohort)
cat(sprintf("\nPRIMARY — Classifier on first %d PCs (per YAML spec):\n", K_CLS))
print(prevalence_primary)

prevalence_full <- tibble(cohort = cohort_vec, cls_scores = cls_full) |>
  group_by(cohort) |>
  count(cls_scores, name = "n") |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  select(-n) |>
  pivot_wider(names_from = cls_scores, values_from = pct, values_fill = 0) |>
  arrange(cohort)
cat(sprintf("\nSENSITIVITY — Classifier on all %d retained PCs:\n", K_retained))
print(prevalence_full)

# Retain for save
cls <- cls_primary
prevalence <- prevalence_primary
prevalence_pct <- prevalence_primary

# ---- Step 7. Save artefacts ------------------------------------------------
cat("\n=== 7. Save artefacts ===\n")
saveRDS(list(mfaces = mfaces, sens = sens,
             retained_primary = retained_primary,
             incretin_axis = inc_axis,
             z_scores = z, classification = cls,
             cohort_vec = cohort_vec,
             prevalence = prevalence,
             prevalence_pct = prevalence_pct),
        epp10_path("fit_mfaces_primary_results.rds"))
cat("Saved: fit_mfaces_primary_results.rds\n")

cat("\n=== PIPELINE COMPLETE ===\n")
