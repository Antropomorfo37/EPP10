# =============================================================================
# fanova_permutation.R — v10.0 §3F.9 formal inference
# =============================================================================
# FANOVA permutation tests on retained mFPC scores:
#   • Univariate F-test per PC (B=5000 permutations)
#   • Multivariate Pillai over ξ₁…ξ_K retained
#   • BH-FDR control α=0.05 stratified by test type
# Contrasts:
#   • Omnibus (all 6 cohorts)
#   • Pairwise each study cohort vs reference no_obese_without_T2DM
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
set.seed(20260422)

RESULTS <- readRDS(epp10_path("fit_mfaces_primary_results.rds"))
scores_primary <- RESULTS$retained_primary$mfpca$scores
cohort_vec     <- RESULTS$cohort_vec
K              <- ncol(scores_primary)
N              <- nrow(scores_primary)
B              <- 5000

cat(sprintf("FANOVA on K=%d retained PCs | N=%d pseudo-subjects | B=%d permutations\n",
            K, N, B))
cat(sprintf("Cohort sizes: "))
cohort_n <- table(cohort_vec); print(cohort_n)

# --- Helper: omnibus F and Pillai on given score matrix × cohort vector -----
fanova_perm_per_pc <- function(xi, cohort, B, seed_base) {
  f_obs <- summary(aov(xi ~ cohort))[[1]]$`F value`[1]
  set.seed(seed_base)
  f_null <- replicate(B, {
    idx <- sample(cohort)
    summary(aov(xi ~ idx))[[1]]$`F value`[1]
  })
  tibble(F_obs = f_obs, p_perm = (1 + sum(f_null >= f_obs, na.rm = TRUE)) / (B + 1))
}

pillai_perm <- function(S, cohort, B, seed_base) {
  observed <- tryCatch(
    summary(manova(S ~ cohort), test = "Pillai")$stats["cohort", "approx F"],
    error = function(e) NA_real_
  )
  if (is.na(observed)) return(tibble(Pillai_F = NA_real_, p_perm = NA_real_))
  set.seed(seed_base)
  null_dist <- replicate(B, {
    idx <- sample(cohort)
    tryCatch(
      summary(manova(S ~ idx), test = "Pillai")$stats["idx", "approx F"],
      error = function(e) NA_real_
    )
  })
  tibble(Pillai_F = observed,
         p_perm = (1 + sum(null_dist >= observed, na.rm = TRUE)) / (B + 1))
}

# =============================================================================
# 1. OMNIBUS — all 6 cohorts
# =============================================================================
cat("\n=== 1. OMNIBUS (all 6 cohorts) ===\n")
t0 <- Sys.time()

# Per-PC F tests
per_pc_omni <- purrr::map_dfr(seq_len(K), function(k) {
  xi <- scores_primary[, k]
  fanova_perm_per_pc(xi, cohort = cohort_vec, B = B,
                     seed_base = 20260422 + k) |>
    dplyr::mutate(pc = paste0("xi", k), contrast = "omnibus")
})

# Multivariate Pillai on retained K
pillai_omni <- pillai_perm(scores_primary, cohort_vec, B = B,
                           seed_base = 20260422 + 1000) |>
  mutate(contrast = "omnibus", test = "Pillai_multivariate")

cat(sprintf("Omnibus computed in %.1f s\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# =============================================================================
# 2. PAIRWISE — each study cohort vs reference
# =============================================================================
cat("\n=== 2. PAIRWISE (each cohort vs no_obese_without_T2DM reference) ===\n")
study_cohorts <- setdiff(unique(cohort_vec), "no_obese_without_T2DM")

pairwise_results <- purrr::map_dfr(study_cohorts, function(co) {
  idx_pair <- cohort_vec %in% c("no_obese_without_T2DM", co)
  S_pair <- scores_primary[idx_pair, , drop = FALSE]
  c_pair <- cohort_vec[idx_pair]

  per_pc <- purrr::map_dfr(seq_len(ncol(S_pair)), function(k) {
    xi <- S_pair[, k]
    fanova_perm_per_pc(xi, cohort = c_pair, B = B,
                       seed_base = 20260422 + 100 + k) |>
      dplyr::mutate(pc = paste0("xi", k),
                    contrast = paste0(co, "_vs_ref"))
  })
  pillai_pair <- pillai_perm(S_pair, c_pair, B = B,
                             seed_base = 20260422 + 200 + nchar(co)) |>
    mutate(contrast = paste0(co, "_vs_ref"), test = "Pillai_multivariate")

  bind_rows(per_pc %>% mutate(test = "F_univariate"),
            pillai_pair %>% rename(F_obs = Pillai_F))
})

# =============================================================================
# 3. Combine + BH-FDR stratified by test type
# =============================================================================
all_results <- bind_rows(
  per_pc_omni %>% mutate(test = "F_univariate"),
  pillai_omni %>% rename(F_obs = Pillai_F),
  pairwise_results
) %>%
  mutate(contrast = factor(contrast,
                           levels = c("omnibus",
                                       paste0(study_cohorts, "_vs_ref"))))

# BH-FDR within each test type (univariate and multivariate separate)
all_results <- all_results %>%
  group_by(test) %>%
  mutate(p_adj_BH = p.adjust(p_perm, method = "BH"),
         significant_FDR_0.05 = p_adj_BH < 0.05) %>%
  ungroup()

# =============================================================================
# 4. Summary tables
# =============================================================================
cat("\n=== 3. OMNIBUS RESULTS ===\n")
omni_tab <- all_results %>% filter(contrast == "omnibus") %>%
  arrange(test, pc) %>%
  select(test, pc, F_obs, p_perm, p_adj_BH, significant_FDR_0.05)
print(omni_tab, n = Inf)

cat("\n=== 4. PAIRWISE — Pillai multivariate per cohort vs reference ===\n")
pillai_pairwise <- all_results %>%
  filter(test == "Pillai_multivariate", contrast != "omnibus") %>%
  select(contrast, F_obs, p_perm, p_adj_BH, significant_FDR_0.05) %>%
  arrange(contrast)
print(pillai_pairwise)

cat("\n=== 5. PAIRWISE — top PC-wise F-tests per contrast (lowest p_adj) ===\n")
per_pc_pairwise <- all_results %>%
  filter(test == "F_univariate", contrast != "omnibus") %>%
  group_by(contrast) %>%
  slice_min(p_adj_BH, n = 5) %>%
  arrange(contrast, p_adj_BH) %>%
  select(contrast, pc, F_obs, p_perm, p_adj_BH, significant_FDR_0.05)
print(per_pc_pairwise, n = Inf)

# Significant-counts summary
cat("\n=== 6. SIGNIFICANCE COUNTS (α_FDR = 0.05) ===\n")
sig_counts <- all_results %>%
  group_by(contrast, test) %>%
  summarise(n_tests = n(),
            n_sig_FDR = sum(significant_FDR_0.05, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(contrast, test)
print(sig_counts, n = Inf)

# =============================================================================
# 5. Save
# =============================================================================
readr::write_csv(all_results, epp10_path("fanova_results.csv"))
saveRDS(list(all_results = all_results, omni_tab = omni_tab,
             pillai_pairwise = pillai_pairwise,
             per_pc_pairwise = per_pc_pairwise,
             sig_counts = sig_counts),
        epp10_path("fanova_results.rds"))
cat("\nSaved: fanova_results.csv, fanova_results.rds\n")
cat("SHA-256 (csv):", digest::digest(
  readr::read_file_raw(epp10_path("fanova_results.csv")),
  algo = "sha256"), "\n")
