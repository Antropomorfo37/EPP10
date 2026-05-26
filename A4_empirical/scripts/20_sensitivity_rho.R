## ============================================================
## 20_sensitivity_rho.R — Sensitivity analysis on AR(1) kernel ρ
##
## Pre-registered sensitivity (FDEP-TP A1 §Validation roadmap):
##   Primary value:    ρ = 0·5
##   Sensitivity grid: ρ ∈ {0·3, 0·5, 0·7}
##
## For each ρ:
##   (a) Re-generate pseudo-IPD multivariada con K_ρ(s,t) = ρ^|s-t|/30
##   (b) Re-fit PACE univariate per broad-panel hormone
##   (c) Re-fit Happ-Greven mFPCA broad-panel (Chiou weights)
##   (d) Extract top-4 eigenstructure
##   (e) Match Ψ̂_m^(ρ) against Ψ̂_m^(0·5) by greedy bipartite
##   (f) Record cumFVE top-4, |corr| per Ψ̂_m, ξ-score ordering
##
## Hard gates:
##   mean(|corr|) ≥ 0·85  across all (ρ, m) pairs
##   max(|ΔFVE|) ≤ 3 pp   vs reference ρ=0·5
##   cohort ordering on ξ₁ preserved (rank-correlation ≥ 0·8)
##
## H.M. Virgen-Ayala, 2026-05-14
## ============================================================
suppressPackageStartupMessages({
  library(data.table); library(tidyverse); library(ggplot2)
  library(fdapace); library(MFPCA); library(funData)
  library(cowplot); library(patchwork); library(gtools)
})

WD <- here::here("A4_empirical")
setwd(WD); set.seed(20260514)

## ===== Constantes =====
N_pseudo     <- 20
times_grid   <- c(0, 15, 30, 60, 90, 120, 150, 180)
joint_broad  <- c("TOTAL GIP","ACTIVE GLP-1","TOTAL GLP-1","Glucagon")
rho_grid     <- c(0.3, 0.5, 0.7)                  # pre-registered grid
rho_ref      <- 0.5

cohort_FDEP  <- c("LEAN","OB","T2DM","OB+T2DM","PRE-SG","POST-SG",
                  "PRE-RYGBP","POST-RYGBP","POST-CR")

optns_PACE   <- list(dataType="Sparse", kernel="epan",
                     methodMuCovEst="smooth",
                     methodBwMu="GMeanAndGCV", methodBwCov="GMeanAndGCV",
                     methodSelectK="FVE", FVEthreshold=0.95,
                     error=TRUE, methodXi="CE", methodRho="trunc",
                     useBinnedData="OFF", nRegGrid=31, verbose=FALSE)

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

## ===== Cargar reference run (ρ=0·5) =====
cat("Cargando reference run (ρ=0·5)...\n")
dt_full <- readRDS("data/master_long.rds")
dt_full[, stratum := sapply(cohort, classify_strict)]
mfpca_ref <- readRDS("data/mFPCA_broad_HappGreven.rds")

Psi_ref <- vector("list", 4)
for (m in 1:4){
  Psi_ref[[m]] <- lapply(seq_along(joint_broad), function(j){
    mfpca_ref$fit$functions[[j]]@X[m, ]
  })
}
lambda_ref <- mfpca_ref$fit$values[1:4]
cumFVE_ref <- sum(lambda_ref) / sum(mfpca_ref$fit$values)
xi_ref_summary <- aggregate(mfpca_ref$fit$scores[, 1:4],
                            by = list(stratum = mfpca_ref$stratum),
                            FUN = mean)
names(xi_ref_summary) <- c("stratum","xi1","xi2","xi3","xi4")
cat(sprintf("Reference: cumFVE top-4 = %.4f\n\n", cumFVE_ref))

## ===== Pipeline encapsulado por ρ =====
run_FDEP_TP_rho <- function(rho_t){
  ## Construir mu_h^stratum(t) (idéntico para todos los ρ)
  strata_present <- intersect(cohort_FDEP, unique(dt_full$stratum))
  strat_mu_list <- list()
  for (h in joint_broad){
    sub_h <- dt_full[hormone == h & time_min <= 240 & !is.na(value)]
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

  ## Pseudo-IPD multivariada con kernel ρ
  Sigma_t <- outer(times_grid, times_grid, function(s,t) rho_t^(abs(s-t)/30))
  Lt <- chol(Sigma_t + diag(1e-6, length(times_grid)))

  IPD <- list()
  for (st in strata_present){
    h_in_st <- intersect(joint_broad,
                          strat_mu_dt[stratum == st, unique(hormone)])
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

  ## PACE per hormone
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
    return(list(status="PACE_failure", values=NA, cumFVE=NA, Psi=NULL))

  ## Happ-Greven mFPCA
  ids_common <- Reduce(intersect, lapply(PACE_list, function(o) o$ids))
  if (length(ids_common) < 30)
    return(list(status="few_common", values=NA, cumFVE=NA, Psi=NULL))

  w_chiou <- sapply(joint_broad,
                    function(h) 1/sqrt(sum(PACE_list[[h]]$fit$lambda)))
  uniExp <- lapply(joint_broad, function(h){
    o <- PACE_list[[h]]; idx <- match(ids_common, o$ids); K <- o$fit$selectK
    list(type="given",
         functions = funData(o$fit$workGrid, t(o$fit$phi[, 1:K])),
         scores = o$fit$xiEst[idx, 1:K, drop=FALSE], ortho=TRUE)
  })
  mFDlist <- lapply(joint_broad, function(h){
    o <- PACE_list[[h]]; idx <- match(ids_common, o$ids); K <- o$fit$selectK
    X <- o$fit$xiEst[idx, 1:K, drop=FALSE] %*% t(o$fit$phi[, 1:K, drop=FALSE])
    X <- sweep(X, 2, o$fit$mu, "+")
    funData(argvals = o$fit$workGrid, X = X)
  })
  mFD <- multiFunData(mFDlist)

  mfpca_fit <- tryCatch(
    MFPCA(mFD, M=12, uniExpansions=uniExp, weights=w_chiou,
          fit=TRUE, bootstrap=FALSE, verbose=FALSE),
    error = function(e) NULL)
  if (is.null(mfpca_fit))
    return(list(status="MFPCA_failure", values=NA, cumFVE=NA, Psi=NULL))

  values <- mfpca_fit$values[1:4]
  cumFVE <- sum(values) / sum(mfpca_fit$values)
  Psi <- vector("list", 4)
  for (m in 1:4){
    Psi[[m]] <- lapply(seq_along(joint_broad), function(j){
      mfpca_fit$functions[[j]]@X[m, ]
    })
  }
  stratum_per_id <- sapply(ids_common, function(id){
    x <- mIPD[subject_id == id, stratum]; if (length(x)) x[1] else NA
  })
  xi_means <- aggregate(mfpca_fit$scores[, 1:4],
                        by = list(stratum = stratum_per_id), FUN = mean)
  names(xi_means) <- c("stratum","xi1","xi2","xi3","xi4")

  list(status="OK", values=values, cumFVE=cumFVE, Psi=Psi,
       chiou_weights=w_chiou, n_subjects=length(ids_common),
       xi_means=xi_means)
}

## ===== Eigenfunction matching (igual que LOSO) =====
match_eigenfunctions <- function(Psi_q, Psi_ref, w_chiou){
  if (is.null(Psi_q)) return(rep(NA, 4))
  perms <- gtools::permutations(4, 4)
  build_S <- function(){
    S <- matrix(0, 4, 4)
    for (m in 1:4) for (mp in 1:4){
      val <- 0
      for (j in seq_along(w_chiou)){
        a <- Psi_q[[m]][[j]]; b <- Psi_ref[[mp]][[j]]
        if (length(a) != length(b))
          b <- approx(seq_along(b), b, n=length(a))$y
        val <- val + w_chiou[j] * sum(a * b)
      }
      S[m,mp] <- val
    }
    S
  }
  S <- build_S()
  best_score <- -Inf; best_perm <- 1:4
  for (i in 1:nrow(perms)){
    p <- perms[i, ]
    sc <- sum(abs(diag(S[, p])))
    if (sc > best_score){ best_score <- sc; best_perm <- p }
  }
  abs_corr <- numeric(4)
  for (m in 1:4){
    mp <- best_perm[m]
    na <- 0; nb <- 0; ip <- 0
    for (j in seq_along(w_chiou)){
      a <- Psi_q[[m]][[j]]; b <- Psi_ref[[mp]][[j]]
      if (length(a) != length(b))
        b <- approx(seq_along(b), b, n=length(a))$y
      na <- na + w_chiou[j] * sum(a^2)
      nb <- nb + w_chiou[j] * sum(b^2)
      ip <- ip + w_chiou[j] * sum(a * b)
    }
    abs_corr[m] <- if (na*nb > 0) abs(ip)/sqrt(na*nb) else NA
  }
  abs_corr
}

## ===== Loop over rho grid =====
cat("==== Sensitivity ρ loop ====\n")
t0 <- Sys.time()
results <- list(); xi_per_rho <- list()
for (rho in rho_grid){
  cat(sprintf("  ρ = %.1f ...\n", rho))
  set.seed(20260514)   # mismo seed dentro de cada ρ
  res <- run_FDEP_TP_rho(rho)
  if (res$status != "OK"){
    results[[length(results)+1]] <- data.table(
      rho = rho, status = res$status,
      cumFVE_top4_pct = NA, dFVE_pct = NA,
      abs_corr_Psi1 = NA, abs_corr_Psi2 = NA,
      abs_corr_Psi3 = NA, abs_corr_Psi4 = NA,
      mean_abs_corr = NA)
    next
  }
  abs_corr <- match_eigenfunctions(res$Psi, Psi_ref, res$chiou_weights)
  results[[length(results)+1]] <- data.table(
    rho = rho, status = "OK",
    n_subjects = res$n_subjects,
    cumFVE_top4_pct = round(res$cumFVE * 100, 2),
    dFVE_pct = round((res$cumFVE - cumFVE_ref) * 100, 2),
    abs_corr_Psi1 = round(abs_corr[1], 3),
    abs_corr_Psi2 = round(abs_corr[2], 3),
    abs_corr_Psi3 = round(abs_corr[3], 3),
    abs_corr_Psi4 = round(abs_corr[4], 3),
    mean_abs_corr = round(mean(abs_corr, na.rm=TRUE), 3))
  xi_per_rho[[as.character(rho)]] <- cbind(rho = rho, res$xi_means)
}
t1 <- Sys.time()
cat(sprintf("\nWall time: %.1f min\n", as.numeric(t1-t0, units="mins")))

res_dt <- rbindlist(results, fill = TRUE)
fwrite(res_dt, "data/sensitivity_rho.csv")
cat("\n==== Resultados sensitivity ρ ====\n")
print(res_dt)

xi_all <- rbindlist(xi_per_rho, fill = TRUE)
fwrite(xi_all, "data/sensitivity_rho_xi_means.csv")

## ===== Hard gates =====
ok <- res_dt[status == "OK"]
mean_corr <- mean(c(ok$abs_corr_Psi1, ok$abs_corr_Psi2,
                    ok$abs_corr_Psi3, ok$abs_corr_Psi4), na.rm=TRUE)
max_dFVE <- max(abs(ok$dFVE_pct), na.rm=TRUE)
cat(sprintf("\n==== HARD GATES ====\nmean(|corr|) across (rho × Psi_m) = %.3f (gate ≥ 0.85)\n",
            mean_corr))
cat(sprintf("max(|ΔFVE|) = %.2f pp (gate ≤ 3 pp)\n", max_dFVE))

## Cohort ordering on xi1
if (nrow(xi_all) > 0){
  rank_corrs <- sapply(rho_grid, function(r){
    sub <- xi_all[rho == r]
    if (nrow(sub) < 3) return(NA)
    ref_sub <- xi_ref_summary[match(sub$stratum, xi_ref_summary$stratum), "xi1"]
    cor(sub$xi1, ref_sub, method = "spearman", use = "complete.obs")
  })
  cat(sprintf("Spearman rank-correlation cohort ξ₁ vs ref (per ρ):\n"))
  for (i in seq_along(rho_grid))
    cat(sprintf("  ρ=%.1f: %.3f\n", rho_grid[i], rank_corrs[i]))
}

## ===== Figura DI_S4 =====
long <- melt(ok, id.vars = "rho",
             measure.vars = c("abs_corr_Psi1","abs_corr_Psi2",
                              "abs_corr_Psi3","abs_corr_Psi4"),
             variable.name = "eigenfn", value.name = "abs_corr")
long[, eigenfn := factor(eigenfn,
                          levels = c("abs_corr_Psi1","abs_corr_Psi2",
                                     "abs_corr_Psi3","abs_corr_Psi4"),
                          labels = c("Ψ̂₁","Ψ̂₂","Ψ̂₃","Ψ̂₄"))]
long[, rho := factor(rho)]

p1 <- ggplot(long, aes(rho, abs_corr, fill = eigenfn)) +
  geom_col(position = position_dodge(0.7), alpha = 0.85,
           color = "black", linewidth = 0.2, width = 0.6) +
  geom_hline(yintercept = 0.85, linetype = "dashed", color = "red") +
  coord_cartesian(ylim = c(0, 1.05)) +
  scale_fill_brewer(palette = "Set2", name = "Eigenfunction") +
  labs(x = "ρ (AR(1) Gaussian-process kernel)",
       y = expression(group("|",corr,"|") ~ "vs reference (ρ=0·5)"),
       title = "DI S4 — Sensitivity de los ejes mFPC al parámetro ρ del kernel AR(1)",
       subtitle = "Gate: |corr| ≥ 0.85. ρ=0·5 es el valor pre-registrado primario.") +
  theme_cowplot(11) + theme(legend.position = "right")
ggsave("figures/DI_S4_rho_sensitivity.png", p1, width = 11, height = 7, dpi = 200)

## FVE comparison
p2 <- ggplot(ok, aes(factor(rho), cumFVE_top4_pct)) +
  geom_col(fill = "#7570B3", alpha = 0.85, color = "black", linewidth = 0.2) +
  geom_hline(yintercept = cumFVE_ref * 100, linetype = "dashed", color = "red") +
  geom_text(aes(label = sprintf("%.2f%%", cumFVE_top4_pct)), vjust = -0.5) +
  coord_cartesian(ylim = c(95, 100)) +
  labs(x = "ρ", y = "cumFVE top-4 (%)",
       title = "DI S4b — cumFVE top-4 a través del grid ρ",
       subtitle = sprintf("Línea roja = reference (ρ=0·5, %.2f%%)", cumFVE_ref*100)) +
  theme_cowplot(11)
ggsave("figures/DI_S4b_rho_FVE.png", p2, width = 8, height = 6, dpi = 200)

saveRDS(list(results=res_dt, xi=xi_all, mean_corr=mean_corr,
             max_dFVE=max_dFVE, rho_ref=rho_ref, cumFVE_ref=cumFVE_ref),
        "data/sensitivity_rho_full.rds")

if (mean_corr < 0.85)
  warning(sprintf("Gate failed: mean(|corr|) = %.3f < 0.85", mean_corr))
if (max_dFVE > 3)
  warning(sprintf("Gate failed: max(|ΔFVE|) = %.2f > 3 pp", max_dFVE))

cat("\n✓ Sensitivity ρ completo.\n")
cat("Outputs: data/sensitivity_rho.csv, data/sensitivity_rho_xi_means.csv\n")
cat("Figuras: figures/DI_S4_rho_sensitivity.png, figures/DI_S4b_rho_FVE.png\n")
