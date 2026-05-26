## ============================================================
## 15_FDEP_TP_mFPCA_cohort.R — mFPCA con pseudo-IPD multivariados
##
## Generamos pseudo-sujetos donde cada sujeto tiene TODAS las hormonas
## simultáneamente (matching by stratum). Esto permite el operador
## conjunto Happ-Greven directo.
## ============================================================
suppressPackageStartupMessages({
  library(data.table); library(tidyverse); library(ggplot2)
  library(fdapace); library(MFPCA); library(funData)
  library(cowplot); library(patchwork); library(viridis)
})

WD <- here::here("A4_empirical")
setwd(WD); set.seed(20260514)

dt <- readRDS("data/master_long.rds")
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
dt[, stratum := sapply(cohort, classify_strict)]

times_grid <- c(0, 15, 30, 60, 90, 120, 150, 180)
analyte_forms <- c("ACYLATED GHRELIN","TOTAL GHRELIN",
                   "ACTIVE GIP","TOTAL GIP",
                   "ACTIVE GLP-1","TOTAL GLP-1",
                   "PYY3-36","TOTAL PYY",
                   "Glucagon","Insulin","Glucose")
joint_panel <- setdiff(analyte_forms, c("Insulin","Glucose"))

## ===== Construir mu_h^(stratum)(t) por interpolación =====
## Para cada (stratum × hormone), media inter-estudio ponderada por N
strat_mu <- list()
for (h in analyte_forms){
  sub_h <- dt[hormone == h & time_min <= 240 & !is.na(value)]
  if (nrow(sub_h) < 4) next
  agg <- sub_h[, .(value = weighted.mean(value, w = N, na.rm = TRUE),
                   n_sub = sum(N, na.rm = TRUE)),
               by = .(stratum, time_min)]
  ## Interpolar a malla común
  out_strat <- list()
  for (st in unique(agg$stratum)){
    arm <- agg[stratum == st][order(time_min)]
    if (nrow(arm) >= 2){
      mu_t <- approx(arm$time_min, arm$value, xout = times_grid, rule = 2)$y
      out_strat[[st]] <- data.table(stratum = st, hormone = h,
                                    time_min = times_grid, mu = mu_t)
    }
  }
  if (length(out_strat)) strat_mu[[h]] <- rbindlist(out_strat)
}
strat_mu_dt <- rbindlist(strat_mu, fill = TRUE)

## ===== Generar pseudo-IPD multivariado: 1 sujeto = 1 vector de TODAS las hormonas =====
## Para cada (stratum), generamos N pseudo-sujetos. Cada sujeto recibe simultáneamente
## una trayectoria para CADA hormona (basada en mu_h^stratum) con correlación AR(1)
strata_present <- unique(strat_mu_dt$stratum)
hormones_present <- unique(strat_mu_dt$hormone)

generate_multivariate_pseudo <- function(stratum_sel, hormones, times,
                                          N_pseudo = 50, rho_t = 0.5){
  ## Construir matriz mu por hormona
  mu_mat <- matrix(NA, nrow = length(hormones), ncol = length(times),
                   dimnames = list(hormones, times))
  for (h in hormones){
    sub <- strat_mu_dt[stratum == stratum_sel & hormone == h]
    if (nrow(sub)) mu_mat[h, ] <- sub$mu[match(times, sub$time_min)]
  }
  ## Hormonas con datos completos en este stratum
  ok_h <- which(rowSums(is.na(mu_mat)) == 0)
  if (length(ok_h) < 3) return(NULL)
  mu_mat <- mu_mat[ok_h, , drop = FALSE]
  H <- nrow(mu_mat); T_ <- length(times)
  ## Generar perturbaciones AR(1) en tiempo, escalar por SD fisiológica (15%)
  Sigma_t <- outer(times, times, function(s, t) rho_t^(abs(s - t)/30))
  Lt <- chol(Sigma_t + diag(1e-6, T_))
  ## Para cada hormona y cada sujeto, draw independiente
  IPD <- list()
  for (i in 1:N_pseudo){
    subj_id <- paste(stratum_sel, sprintf("S%03d", i), sep = "_")
    for (h_name in rownames(mu_mat)){
      sd_h <- 0.15 * (max(mu_mat[h_name, ]) - min(mu_mat[h_name, ]) + 1)
      Z <- rnorm(T_)
      perturb <- as.vector(sd_h * (Z %*% Lt))
      vals <- pmax(0, mu_mat[h_name, ] + perturb)   # mantener >0
      IPD[[length(IPD) + 1]] <- data.table(
        subject_id = subj_id, stratum = stratum_sel,
        hormone = h_name, time_min = times, value = vals
      )
    }
  }
  rbindlist(IPD)
}

cat("Generando pseudo-IPD multivariado por estrato...\n")
mIPD <- list()
for (st in strata_present){
  res <- generate_multivariate_pseudo(st, analyte_forms, times_grid, N_pseudo = 50)
  if (!is.null(res)){
    mIPD[[st]] <- res
    cat(sprintf("  %s: %d obs (%d sujetos × %d hormonas × %d tiempos)\n",
                st, nrow(res), uniqueN(res$subject_id),
                uniqueN(res$hormone), uniqueN(res$time_min)))
  }
}
mIPD_dt <- rbindlist(mIPD, fill = TRUE)
saveRDS(mIPD_dt, "data/pseudo_IPD_multivariate.rds")
cat(sprintf("Total: %d filas, %d sujetos\n", nrow(mIPD_dt), uniqueN(mIPD_dt$subject_id)))

## ===== PACE univariado por hormona sobre pseudo-IPD multivariado =====
cat("\n[PACE multivariado] por hormona...\n")
optns <- list(dataType="Sparse", kernel="epan", methodMuCovEst="smooth",
              methodBwMu="GMeanAndGCV", methodBwCov="GMeanAndGCV",
              methodSelectK="FVE", FVEthreshold=0.95, error=TRUE,
              methodXi="CE", methodRho="trunc", useBinnedData="OFF",
              nRegGrid=31, verbose=FALSE)
PACE_list <- list()
for (h in analyte_forms){
  sub <- mIPD_dt[hormone == h]
  if (nrow(sub) < 50) next
  curves <- split(sub, sub$subject_id)
  Lt <- lapply(curves, function(x) x$time_min)
  Ly <- lapply(curves, function(x) x$value)
  fit <- tryCatch(FPCA(Ly, Lt, optns = optns), error = function(e) NULL)
  if (is.null(fit)){ cat(h, ": FAIL\n"); next }
  cohort_per <- sapply(curves, function(x) unique(x$stratum)[1])
  PACE_list[[h]] <- list(hormone = h, fit = fit,
                         cohort = unname(cohort_per), ids = names(curves))
  cat(sprintf("  %s: K=%d, FVE=%.2f\n", h, fit$selectK, tail(fit$cumFVE,1)))
}
saveRDS(PACE_list, "data/PACE_multivariate_per_hormone.rds")

## ===== Happ-Greven mFPCA con normalización Chiou =====
joint_hormones <- intersect(joint_panel, names(PACE_list))
cat(sprintf("\n[Happ-Greven mFPCA] sobre %d hormonas joint panel\n", length(joint_hormones)))

ids_common <- Reduce(intersect, lapply(joint_hormones, function(h) PACE_list[[h]]$ids))
cat(sprintf("  Sujetos comunes: %d\n", length(ids_common)))

if (length(ids_common) >= 30 && length(joint_hormones) >= 4){
  ## Pesos Chiou
  w_chiou <- sapply(joint_hormones, function(h){
    1 / sqrt(sum(PACE_list[[h]]$fit$lambda))
  })

  uniExp <- lapply(joint_hormones, function(h){
    o <- PACE_list[[h]]
    idx <- match(ids_common, o$ids)
    list(type = "given",
         functions = funData(o$fit$workGrid, t(o$fit$phi[, 1:o$fit$selectK])),
         scores = o$fit$xiEst[idx, 1:o$fit$selectK, drop = FALSE],
         ortho = TRUE)
  })

  ## Construir funData object multivariado
  mFDlist <- lapply(joint_hormones, function(h){
    o <- PACE_list[[h]]
    idx <- match(ids_common, o$ids)
    K <- o$fit$selectK
    ## X̂ = ξ × φᵀ + μ (reconstrucción truncada)
    X <- o$fit$xiEst[idx, 1:K, drop = FALSE] %*% t(o$fit$phi[, 1:K, drop = FALSE])
    X <- sweep(X, 2, o$fit$mu, "+")
    funData(argvals = o$fit$workGrid, X = X)
  })
  mFD <- multiFunData(mFDlist)

  cat("  Ejecutando MFPCA() (Happ-Greven, 12 componentes)...\n")
  mfpca_fit <- tryCatch(
    MFPCA(mFD, M = 12, uniExpansions = uniExp, weights = w_chiou,
          fit = TRUE, bootstrap = FALSE, verbose = FALSE),
    error = function(e){ message("MFPCA: ", e$message); NULL })

  if (!is.null(mfpca_fit)){
    cumFVE_mFPCA <- cumsum(mfpca_fit$values) / sum(mfpca_fit$values)
    cat(sprintf("  ✓ mFPCA OK. Cumulative FVE top-4: %s\n",
                paste(round(cumFVE_mFPCA[1:4]*100, 1), collapse = ", ")))
    cat(sprintf("  ✓ Per-eigenvalue FVE top-6: %s\n",
                paste(round(mfpca_fit$values[1:6]/sum(mfpca_fit$values)*100, 1),
                      collapse = ", ")))

    stratum_per_id <- sapply(ids_common, function(id){
      x <- mIPD_dt[subject_id == id, stratum]; if (length(x)) x[1] else NA
    })

    mfpca_pkg <- list(fit = mfpca_fit, ids = ids_common,
                      stratum = unname(stratum_per_id),
                      joint_hormones = joint_hormones,
                      chiou_weights = w_chiou,
                      cumFVE = cumFVE_mFPCA)
    saveRDS(mfpca_pkg, "data/mFPCA_HappGreven_final.rds")
    cat("  ✓ Resultados guardados en data/mFPCA_HappGreven_final.rds\n")
  }
}

cat("\n✓ Pipeline mFPCA multivariado completo\n")
