## ============================================================
## 19_LOSO_sensitivity.R — Leave-One-Study-Out (LOSO) sensitivity
##
## Para cada uno de los 23 estudios fuente:
##   (a) excluye el estudio del master_long.rds
##   (b) regenera pseudo-IPD multivariada (N=20 por arm)
##   (c) re-ajusta PACE por hormona del broad-panel
##   (d) re-corre Happ-Greven mFPCA broad-panel con pesos Chiou
##   (e) extrae top-4 eigenvalues, FVE acumulado, eigenfunctions
##   (f) empareja Ψ̂_m^(-study) con Ψ̂_m^(full) por matching greedy
##       de inner products ponderados Chiou
##   (g) reporta |corr|, ΔFVE %, status por estudio
##
## Verificación: mean(|corr|) ≥ 0·85 ; max(|ΔFVE|) ≤ 3 pp
##
## H.M. Virgen-Ayala, 2026-05-14
## ============================================================
suppressPackageStartupMessages({
  library(data.table); library(tidyverse); library(ggplot2)
  library(fdapace); library(MFPCA); library(funData)
  library(future.apply); library(cowplot); library(patchwork)
})

WD <- here::here("A4_empirical")
setwd(WD); set.seed(20260514)

## ===== Constantes =====
N_pseudo     <- 20                       # pseudo-subjects per arm
rho_AR1      <- 0.5                      # AR(1) Gaussian-process kernel
times_grid   <- c(0, 15, 30, 60, 90, 120, 150, 180)
joint_broad  <- c("TOTAL GIP","ACTIVE GLP-1","TOTAL GLP-1","Glucagon")
optns_PACE   <- list(dataType="Sparse", kernel="epan",
                     methodMuCovEst="smooth",
                     methodBwMu="GMeanAndGCV",
                     methodBwCov="GMeanAndGCV",
                     methodSelectK="FVE", FVEthreshold=0.95,
                     error=TRUE, methodXi="CE", methodRho="trunc",
                     useBinnedData="OFF", nRegGrid=31, verbose=FALSE)

cohort_FDEP  <- c("LEAN","OB","T2DM","OB+T2DM","PRE-SG","POST-SG",
                  "PRE-RYGBP","POST-RYGBP","POST-CR")

classify_strict <- function(x){
  x <- tolower(trimws(as.character(x)))
  if (grepl("post-rygb|after.*roux", x))                    return("POST-RYGBP")
  if (grepl("pre-rygb|before.*rygb|before.*roux", x))       return("PRE-RYGBP")
  if (grepl("post-sg|after.*sleeve", x))                    return("POST-SG")
  if (grepl("pre-sg|before.*sleeve|before sg", x))          return("PRE-SG")
  if (grepl("post-cr|after.*calor", x))                     return("POST-CR")
  if (grepl("plus.*type 2|with.*t2dm|with type 2", x))      return("OB+T2DM")
  if (grepl("type 2 diabetes|^t2dm", x))                    return("T2DM")
  if (x == "lean" || grepl("no obesity-no t2dm", x))        return("LEAN")
  return("OB")
}

## ===== Carga full-run reference =====
cat("Cargando full-run reference...\n")
dt_full <- readRDS("data/master_long.rds")
dt_full[, stratum := sapply(cohort, classify_strict)]
mfpca_full <- readRDS("data/mFPCA_broad_HappGreven.rds")

Psi_full <- vector("list", 4)
for (m in 1:4){
  Psi_full[[m]] <- lapply(seq_along(joint_broad), function(j){
    mfpca_full$fit$functions[[j]]@X[m, ]
  })
}
lambda_full   <- mfpca_full$fit$values[1:4]
cumFVE_full   <- sum(lambda_full) / sum(mfpca_full$fit$values)
cat(sprintf("Reference: K=%d, cumFVE top-4 = %.3f, eigenvalues = %s\n\n",
            4, cumFVE_full, paste(round(lambda_full,2), collapse=", ")))

## ===== Estudios únicos =====
studies <- sort(unique(dt_full$study))
studies <- studies[nchar(studies) >= 3]
cat(sprintf("Estudios fuente: %d\n", length(studies)))
print(studies)

## ===== Pipeline encapsulado (Layers 2-5) =====
run_FDEP_TP_LOSO <- function(dt_in, N_pseudo = 20, rho_t = 0.5){
  ## Layer 2: pseudo-IPD multivariada
  ## (Replica logic de scripts 14-15-17 sin escribir a disco)
  strata_present <- intersect(cohort_FDEP, unique(dt_in$stratum))
  hormones_present <- intersect(joint_broad, unique(dt_in$hormone))
  if (length(strata_present) < 6 || length(hormones_present) < 3)
    return(list(status = "insufficient_coverage",
                values = NA, cumFVE = NA, Psi = NULL))

  ## Construir mu_h^stratum(t)
  strat_mu_list <- list()
  for (h in joint_broad){
    sub_h <- dt_in[hormone == h & time_min <= 240 & !is.na(value)]
    if (nrow(sub_h) < 4) next
    agg <- sub_h[, .(value = weighted.mean(value, w = N, na.rm = TRUE)),
                 by = .(stratum, time_min)]
    for (st in unique(agg$stratum)){
      arm <- agg[stratum == st][order(time_min)]
      if (nrow(arm) >= 2){
        mu_t <- approx(arm$time_min, arm$value, xout = times_grid, rule = 2)$y
        strat_mu_list[[length(strat_mu_list)+1]] <-
          data.table(stratum = st, hormone = h,
                     time_min = times_grid, mu = mu_t)
      }
    }
  }
  strat_mu_dt <- rbindlist(strat_mu_list, fill = TRUE)
  if (nrow(strat_mu_dt) == 0)
    return(list(status = "no_data", values = NA, cumFVE = NA, Psi = NULL))

  ## Generar pseudo-IPD multivariado
  Sigma_t <- outer(times_grid, times_grid, function(s,t) rho_t^(abs(s-t)/30))
  Lt <- tryCatch(chol(Sigma_t + diag(1e-6, length(times_grid))),
                 error = function(e) NULL)
  if (is.null(Lt))
    return(list(status = "chol_failure", values = NA, cumFVE = NA, Psi = NULL))

  IPD <- list()
  for (st in strata_present){
    h_in_st <- intersect(joint_broad, strat_mu_dt[stratum == st, unique(hormone)])
    if (length(h_in_st) < 3) next
    for (i in 1:N_pseudo){
      subj_id <- paste(st, sprintf("S%03d", i), sep = "_")
      for (h_name in h_in_st){
        mu <- strat_mu_dt[stratum == st & hormone == h_name, mu]
        if (length(mu) != length(times_grid)) next
        sd_h <- 0.15 * (max(mu) - min(mu) + 1)
        perturb <- as.vector(sd_h * (rnorm(length(times_grid)) %*% Lt))
        vals <- pmax(0, mu + perturb)
        IPD[[length(IPD)+1]] <- data.table(
          subject_id = subj_id, stratum = st,
          hormone = h_name, time_min = times_grid, value = vals
        )
      }
    }
  }
  mIPD <- rbindlist(IPD, fill = TRUE)
  if (nrow(mIPD) == 0)
    return(list(status = "empty_IPD", values = NA, cumFVE = NA, Psi = NULL))

  ## Layer 3: PACE por hormona
  PACE_list <- list()
  for (h in joint_broad){
    sub <- mIPD[hormone == h]
    if (nrow(sub) < 50) next
    curves <- split(sub, sub$subject_id)
    Lt_lst <- lapply(curves, function(x) x$time_min)
    Ly_lst <- lapply(curves, function(x) x$value)
    fit <- tryCatch(FPCA(Ly_lst, Lt_lst, optns = optns_PACE),
                    error = function(e) NULL)
    if (is.null(fit)) next
    PACE_list[[h]] <- list(hormone = h, fit = fit, ids = names(curves),
                           cohort = sapply(curves, function(x) x$stratum[1]))
  }
  if (length(PACE_list) < length(joint_broad))
    return(list(status = "PACE_failure", values = NA, cumFVE = NA, Psi = NULL))

  ## Layer 5: Happ-Greven mFPCA broad-panel
  ids_common <- Reduce(intersect, lapply(PACE_list, function(o) o$ids))
  if (length(ids_common) < 30)
    return(list(status = "few_common_subjects", values = NA, cumFVE = NA, Psi = NULL))

  w_chiou <- sapply(joint_broad, function(h) 1/sqrt(sum(PACE_list[[h]]$fit$lambda)))
  uniExp <- lapply(joint_broad, function(h){
    o <- PACE_list[[h]]
    idx <- match(ids_common, o$ids)
    K <- o$fit$selectK
    list(type="given",
         functions = funData(o$fit$workGrid, t(o$fit$phi[, 1:K])),
         scores = o$fit$xiEst[idx, 1:K, drop=FALSE], ortho = TRUE)
  })
  mFDlist <- lapply(joint_broad, function(h){
    o <- PACE_list[[h]]
    idx <- match(ids_common, o$ids)
    K <- o$fit$selectK
    X <- o$fit$xiEst[idx, 1:K, drop=FALSE] %*% t(o$fit$phi[, 1:K, drop=FALSE])
    X <- sweep(X, 2, o$fit$mu, "+")
    funData(argvals = o$fit$workGrid, X = X)
  })
  mFD <- multiFunData(mFDlist)

  mfpca_fit <- tryCatch(
    MFPCA(mFD, M = 12, uniExpansions = uniExp, weights = w_chiou,
          fit = TRUE, bootstrap = FALSE, verbose = FALSE),
    error = function(e) NULL)
  if (is.null(mfpca_fit))
    return(list(status = "MFPCA_failure", values = NA, cumFVE = NA, Psi = NULL))

  ## Extraer top-4 eigenstructure
  values  <- mfpca_fit$values[1:4]
  total   <- sum(mfpca_fit$values)
  cumFVE  <- sum(values) / total
  Psi <- vector("list", 4)
  for (m in 1:4){
    Psi[[m]] <- lapply(seq_along(joint_broad), function(j){
      mfpca_fit$functions[[j]]@X[m, ]
    })
  }
  list(status = "OK", values = values, cumFVE = cumFVE, Psi = Psi,
       n_strata = length(strata_present),
       n_subjects = length(ids_common),
       chiou_weights = w_chiou)
}

## ===== Greedy bipartite matching (4×4) =====
match_eigenfunctions <- function(Psi_loso, Psi_full, w_chiou){
  if (is.null(Psi_loso)) return(list(abs_corr = rep(NA, 4), sigma = NA))
  ## S[m, m'] = Σ_j w_j × <Psi_loso[[m]][[j]], Psi_full[[m']][[j]]>
  S <- matrix(0, 4, 4)
  for (m in 1:4) for (mp in 1:4){
    if (is.null(Psi_loso[[m]]) || is.null(Psi_full[[mp]])) {
      S[m,mp] <- 0; next
    }
    val <- 0
    for (j in seq_along(w_chiou)){
      a <- Psi_loso[[m]][[j]]; b <- Psi_full[[mp]][[j]]
      if (length(a) != length(b)) {
        b <- approx(seq_along(b), b, n = length(a))$y
      }
      val <- val + w_chiou[j] * sum(a * b)
    }
    S[m,mp] <- val
  }
  ## Greedy matching on |S| (4×4 small, exact via brute force permutations)
  perms <- gtools::permutations(4, 4)
  best_score <- -Inf; best_perm <- 1:4
  for (i in 1:nrow(perms)){
    p <- perms[i, ]
    sc <- sum(abs(diag(S[, p])))
    if (sc > best_score){ best_score <- sc; best_perm <- p }
  }
  ## Normalise correlations (S already has weights and absolute inner products)
  ## Compute true |corr| as |S[m,σ(m)]| / sqrt(||a||²·||b||²) (Chiou-weighted)
  abs_corr <- numeric(4)
  for (m in 1:4){
    mp <- best_perm[m]
    norm_a <- 0; norm_b <- 0; ip <- 0
    for (j in seq_along(w_chiou)){
      if (is.null(Psi_loso[[m]]) || is.null(Psi_full[[mp]])) next
      a <- Psi_loso[[m]][[j]]; b <- Psi_full[[mp]][[j]]
      if (length(a) != length(b)) b <- approx(seq_along(b), b, n=length(a))$y
      norm_a <- norm_a + w_chiou[j] * sum(a^2)
      norm_b <- norm_b + w_chiou[j] * sum(b^2)
      ip     <- ip     + w_chiou[j] * sum(a * b)
    }
    abs_corr[m] <- if (norm_a*norm_b > 0) abs(ip) / sqrt(norm_a*norm_b) else NA
  }
  list(abs_corr = abs_corr, sigma = best_perm)
}

## Install gtools if missing
if (!requireNamespace("gtools", quietly = TRUE))
  install.packages("gtools", repos = "https://cloud.r-project.org", quiet = TRUE)

## ===== Outer LOSO loop =====
cat("\n==== LOSO loop: 23 estudios × ~1 min cada uno ====\n")
plan(multisession, workers = 8)
options(future.rng.onMisuse = "ignore")

t0 <- Sys.time()
results_list <- future_lapply(studies, function(study_j){
  cat(sprintf("  [%s] excluyendo...\n", study_j))
  dt_loso <- dt_full[study != study_j]
  res <- tryCatch(run_FDEP_TP_LOSO(dt_loso, N_pseudo = N_pseudo, rho_t = rho_AR1),
                  error = function(e){
                    list(status = paste0("exception:", e$message),
                         values = NA, cumFVE = NA, Psi = NULL)
                  })
  if (res$status == "OK"){
    m <- match_eigenfunctions(res$Psi, Psi_full, res$chiou_weights)
    return(data.table(
      study = study_j,
      n_strata_retained = res$n_strata,
      n_subjects = res$n_subjects,
      FVE_top4_pct = round(res$cumFVE * 100, 2),
      dFVE_pct = round((res$cumFVE - cumFVE_full) * 100, 2),
      abs_corr_Psi1 = round(m$abs_corr[1], 3),
      abs_corr_Psi2 = round(m$abs_corr[2], 3),
      abs_corr_Psi3 = round(m$abs_corr[3], 3),
      abs_corr_Psi4 = round(m$abs_corr[4], 3),
      mean_abs_corr = round(mean(m$abs_corr, na.rm = TRUE), 3),
      status = "OK"))
  } else {
    return(data.table(
      study = study_j, n_strata_retained = NA, n_subjects = NA,
      FVE_top4_pct = NA, dFVE_pct = NA,
      abs_corr_Psi1 = NA, abs_corr_Psi2 = NA,
      abs_corr_Psi3 = NA, abs_corr_Psi4 = NA,
      mean_abs_corr = NA, status = res$status))
  }
}, future.seed = TRUE)
plan(sequential)
t1 <- Sys.time()
cat(sprintf("\nWall time LOSO: %.1f min\n", as.numeric(t1 - t0, units = "mins")))

res_dt <- rbindlist(results_list, fill = TRUE)
fwrite(res_dt, "data/LOSO_sensitivity.csv")
cat("\n==== Resultados LOSO ====\n")
print(res_dt)

## ===== Verification gates =====
ok_runs <- res_dt[status == "OK"]
cat(sprintf("\nRuns OK: %d / %d\n", nrow(ok_runs), nrow(res_dt)))
if (nrow(ok_runs) > 0){
  cat(sprintf("Mean |corr| across all (study × Ψm): %.3f\n",
              mean(c(ok_runs$abs_corr_Psi1, ok_runs$abs_corr_Psi2,
                     ok_runs$abs_corr_Psi3, ok_runs$abs_corr_Psi4),
                   na.rm = TRUE)))
  cat(sprintf("Mean |mean_abs_corr|: %.3f\n",
              mean(ok_runs$mean_abs_corr, na.rm = TRUE)))
  cat(sprintf("Max |ΔFVE|: %.2f pp\n",
              max(abs(ok_runs$dFVE_pct), na.rm = TRUE)))
}

## ===== Figuras =====
if (nrow(ok_runs) > 0){
  long <- melt(ok_runs, id.vars = "study",
               measure.vars = c("abs_corr_Psi1","abs_corr_Psi2",
                                "abs_corr_Psi3","abs_corr_Psi4"),
               variable.name = "eigenfn", value.name = "abs_corr")
  long[, eigenfn := factor(eigenfn,
                            levels = c("abs_corr_Psi1","abs_corr_Psi2",
                                       "abs_corr_Psi3","abs_corr_Psi4"),
                            labels = c("Ψ̂₁","Ψ̂₂","Ψ̂₃","Ψ̂₄"))]
  long[, study := factor(study, levels = studies)]

  p1 <- ggplot(long, aes(study, abs_corr)) +
    geom_col(fill = "#1B9E77", alpha = 0.85, color = "black", linewidth = 0.2) +
    geom_hline(yintercept = 0.85, linetype = "dashed", color = "red") +
    facet_wrap(~ eigenfn, ncol = 2) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = NULL, y = expression(group("|",corr,"|") ~ "vs full-run"),
         title = "DI S3 — LOSO sensitivity: alineamiento de eigenfunctions Ψ̂_m^(-study) vs Ψ̂_m^(full)",
         subtitle = sprintf("Gate: |corr| ≥ 0.85 (línea roja). N_pseudo=%d.", N_pseudo)) +
    theme_cowplot(10) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  ggsave("figures/DI_S3_LOSO_eigenfunctions.png", p1, width = 13, height = 8, dpi = 200)

  p2 <- ggplot(ok_runs, aes(reorder(study, FVE_top4_pct), FVE_top4_pct)) +
    geom_col(fill = "#7570B3", alpha = 0.85, color = "black", linewidth = 0.2) +
    geom_hline(yintercept = cumFVE_full * 100, linetype = "dashed", color = "red") +
    coord_flip() +
    labs(y = "FVE top-4 (%) LOSO", x = NULL,
         title = "DI S3 — Cumulative FVE top-4 con LOSO sensitivity",
         subtitle = sprintf("Línea roja = full-run reference (%.1f%%)", cumFVE_full*100)) +
    theme_cowplot(10)
  ggsave("figures/DI_S3_LOSO_FVE.png", p2, width = 9, height = 8, dpi = 200)
  cat("\n✓ Figuras DI_S3 generadas.\n")
}

saveRDS(list(results = res_dt, ref_cumFVE = cumFVE_full,
             ref_lambda = lambda_full, N_pseudo = N_pseudo),
        "data/LOSO_sensitivity_full.rds")

## Hard gates
if (nrow(ok_runs) > 0){
  mean_corr <- mean(ok_runs$mean_abs_corr, na.rm = TRUE)
  max_dFVE  <- max(abs(ok_runs$dFVE_pct), na.rm = TRUE)
  cat(sprintf("\n==== HARD GATES ====\nmean(|corr|) = %.3f (gate ≥ 0.85)\n", mean_corr))
  cat(sprintf("max(|ΔFVE|) = %.2f pp (gate ≤ 3 pp)\n", max_dFVE))
  if (mean_corr < 0.85)
    warning(sprintf("Gate failed: mean(|corr|) = %.3f < 0.85", mean_corr))
  if (max_dFVE > 3)
    warning(sprintf("Gate failed: max(|ΔFVE|) = %.2f > 3 pp", max_dFVE))
}

cat("\n✓ LOSO sensitivity completo. Outputs: data/LOSO_sensitivity.csv, figures/DI_S3_*.png\n")
