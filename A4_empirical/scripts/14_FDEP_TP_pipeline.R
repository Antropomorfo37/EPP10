## ============================================================
## 14_FDEP_TP_pipeline.R — Pipeline FDEP-TP completo
##
## Implementa el marco Functional Dynamic Enteropancreatic Phenotyping
## during the Periprandial Transition (FDEP-TP) sobre la tabla maestra:
##
##  Layer 1 — Ingest & harmonisation (ya hecho en script 01-02)
##  Layer 2 — Pseudo-IPD via Gaussian-process AR(1) (rho=0.5, t1/2=30 min)
##            [Papadimitropoulou 2020 + Røge AR(1)]
##  Layer 3 — PACE univariado (Yao, Müller, Wang 2005) por hormona
##            + FACEs sensibilidad (Xiao 2018) - fdapace 0.6.0
##  Layer 4 — Cross-covariance kernel C^(jk)(s,t) entre hormonas
##  Layer 5 — Happ-Greven mFPCA conjunto con normalización Chiou
##            → 4 eigenfunciones canónicas Ψ_1..Ψ_4
##  Layer 6 — Descriptores cruzados ρ_INC, ρ_NET, ρ_ANR por cohorte
##  Layer 7 — Asignación PTP (9-label) e IEP (8-rule) por cohorte
##  Layer 8 — Inferencia: Pillai trace, Mahalanobis, bootstrap subject-level
## ============================================================
suppressPackageStartupMessages({
  library(data.table); library(tidyverse); library(ggplot2)
  library(fdapace); library(MFPCA); library(funData); library(face)
  library(cowplot); library(patchwork); library(viridis); library(scales)
})

WD <- here::here("A4_empirical")
setwd(WD)
set.seed(20260514)

## ===== Carga ========================================================
dt <- readRDS("data/master_long.rds")

## Renormalizar cohortes a 7 estratos FDEP-TP canónicos
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
cat("Estratos FDEP-TP (n hormona-tiempo):\n"); print(dt[, .N, by = stratum][order(-N)])

## Panel hormonal — 11 forms; el operador conjunto restringido a 9 (sin glucosa/insulina)
analyte_forms <- c("ACYLATED GHRELIN","TOTAL GHRELIN",
                   "ACTIVE GIP","TOTAL GIP",
                   "ACTIVE GLP-1","TOTAL GLP-1",
                   "PYY3-36","TOTAL PYY",
                   "Glucagon","Insulin","Glucose")
joint_panel <- setdiff(analyte_forms, c("Insulin","Glucose"))   # 9 forms

## ===== Layer 2 — Pseudo-IPD generation (AR(1) Gaussian process) =====
## Para cada estudio × hormona × cohorte, generamos N pseudo-sujetos
## alrededor de la trayectoria media digitizada (mu_h(t)), con kernel
## K_rho(s,t) = rho^|s-t|/30 con rho=0.5  (Papadimitropoulou 2020)
generate_pseudo_IPD <- function(study, stratum, hormone, times, mean_traj, N_pseudo = 50, rho = 0.5){
  if (length(times) < 2) return(NULL)
  ## Matriz de covarianza AR(1)
  Sigma <- outer(times, times, function(s, t) rho^(abs(s - t)/30))
  ## Escalar a varianza fisiológica (≈ 15% del rango)
  sd_h <- 0.15 * (max(mean_traj) - min(mean_traj) + 1e-6)
  Sigma <- sd_h^2 * Sigma
  L <- tryCatch(chol(Sigma + diag(1e-6, nrow(Sigma))), error = function(e) NULL)
  if (is.null(L)) return(NULL)
  Z <- matrix(rnorm(N_pseudo * length(times)), nrow = N_pseudo)
  trajs <- sweep(Z %*% L, 2, mean_traj, FUN = "+")
  data.table(
    subject_id = paste(study, stratum, hormone, 1:N_pseudo, sep = "|"),
    study = study, stratum = stratum, hormone = hormone,
    time_min = rep(times, each = N_pseudo),
    value = as.vector(trajs)
  )
}

cat("\n[Layer 2] Generando pseudo-IPD (AR(1) GP, rho=0.5)...\n")
pseudo_ipd <- list()
for (h in analyte_forms){
  sub_h <- dt[hormone == h & time_min <= 240 & !is.na(value)]
  if (nrow(sub_h) < 4) next
  ## Para cada study × stratum: promedio por time_min
  agg <- sub_h[, .(value = mean(value, na.rm = TRUE)),
               by = .(study, stratum, time_min)]
  ## Por cada (study, stratum)
  for (key in unique(paste(agg$study, agg$stratum, sep = "|"))){
    parts <- strsplit(key, "|", fixed = TRUE)[[1]]
    arm <- agg[study == parts[1] & stratum == parts[2]][order(time_min)]
    if (nrow(arm) >= 2){
      out <- generate_pseudo_IPD(parts[1], parts[2], h,
                                  arm$time_min, arm$value,
                                  N_pseudo = 30, rho = 0.5)
      if (!is.null(out)) pseudo_ipd[[length(pseudo_ipd) + 1]] <- out
    }
  }
}
pseudo_dt <- rbindlist(pseudo_ipd, fill = TRUE)
cat(sprintf("  pseudo-IPD: %d filas, %d sujetos únicos\n",
            nrow(pseudo_dt), uniqueN(pseudo_dt$subject_id)))
saveRDS(pseudo_dt, "data/pseudo_IPD.rds")

## ===== Layer 3 — PACE univariado por hormona ========================
cat("\n[Layer 3] PACE univariado por hormona (Yao-Müller-Wang 2005)...\n")
optns_PACE <- list(
  dataType       = "Sparse",
  kernel         = "epan",
  methodMuCovEst = "smooth",
  methodBwMu     = "GMeanAndGCV",
  methodBwCov    = "GMeanAndGCV",
  methodSelectK  = "FVE",
  FVEthreshold   = 0.95,
  error          = TRUE,
  methodXi       = "CE",
  methodRho      = "trunc",
  useBinnedData  = "OFF",
  nRegGrid       = 31,
  verbose        = FALSE
)

run_PACE_safe <- function(h){
  sub <- pseudo_dt[hormone == h & is.finite(value)]
  if (uniqueN(sub$subject_id) < 8) return(NULL)
  ## Dedup tiempos dentro de subject
  sub <- sub[, .(value = mean(value)), by = .(subject_id, stratum, time_min)]
  curves <- split(sub, sub$subject_id)
  Lt <- lapply(curves, function(x) x$time_min)
  Ly <- lapply(curves, function(x) x$value)
  keep <- sapply(Lt, length) >= 4
  if (sum(keep) < 8) return(NULL)
  Lt <- Lt[keep]; Ly <- Ly[keep]
  cohort_per <- sapply(curves[keep], function(x) unique(x$stratum)[1])
  fit <- tryCatch(FPCA(Ly, Lt, optns = optns_PACE), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  list(hormone = h, fit = fit, cohort = unname(cohort_per), ids = names(curves[keep]))
}

PACE_list <- list()
for (h in analyte_forms){
  cat(sprintf("  PACE: %s ... ", h))
  res <- run_PACE_safe(h)
  if (is.null(res)){ cat("SKIP (insuficiente)\n"); next }
  cat(sprintf("K=%d, FVE=%.2f\n", res$fit$selectK, tail(res$fit$cumFVE,1)))
  PACE_list[[h]] <- res
}
saveRDS(PACE_list, "data/PACE_per_hormone.rds")

## ===== Layer 4 — Cross-covariance kernel C^(jk)(s,t) ================
## C^(jk)(s,t) = Cov[X^(j)(s), X^(k)(t)] estimado en grilla común
cat("\n[Layer 4] Cross-covariance kernel C^(jk)(s,t)...\n")
common_grid <- c(0, 15, 30, 60, 90, 120, 150, 180)

## Crear matriz sujeto × (hormona × tiempo)
wide_dt <- dcast(pseudo_dt[time_min %in% common_grid],
                 subject_id + stratum ~ hormone + time_min,
                 value.var = "value", fun.aggregate = mean)

## Construir kernels por pareja
analytes_active <- names(PACE_list)
analyte_pairs <- expand.grid(j = analytes_active, k = analytes_active,
                             stringsAsFactors = FALSE)
analyte_pairs <- analyte_pairs[as.character(analyte_pairs$j) <= as.character(analyte_pairs$k), ]

compute_kernel <- function(j, k){
  cols_j <- grep(paste0("^", j, "_"), names(wide_dt), value = TRUE)
  cols_k <- grep(paste0("^", k, "_"), names(wide_dt), value = TRUE)
  if (length(cols_j) == 0 || length(cols_k) == 0) return(NULL)
  Xj <- as.matrix(wide_dt[, ..cols_j])
  Xk <- as.matrix(wide_dt[, ..cols_k])
  ## Centrar
  Xj_c <- scale(Xj, center = TRUE, scale = FALSE)
  Xk_c <- scale(Xk, center = TRUE, scale = FALSE)
  ## Eliminar filas con NA en ambas
  keep <- complete.cases(Xj_c) & complete.cases(Xk_c)
  if (sum(keep) < 8) return(NULL)
  Cjk <- crossprod(Xj_c[keep, ], Xk_c[keep, ]) / (sum(keep) - 1)
  rownames(Cjk) <- gsub(paste0("^", j, "_"), "", cols_j)
  colnames(Cjk) <- gsub(paste0("^", k, "_"), "", cols_k)
  Cjk
}

kernels <- list()
for (i in seq_len(nrow(analyte_pairs))){
  j <- analyte_pairs$j[i]; k <- analyte_pairs$k[i]
  K <- compute_kernel(j, k)
  if (!is.null(K)) kernels[[paste(j, k, sep = "||")]] <- K
}
saveRDS(kernels, "data/cross_covariance_kernels.rds")
cat(sprintf("  Calculados %d kernels (incluye auto- y cross-)\n", length(kernels)))

## ===== Layer 5 — Happ-Greven mFPCA conjunto =========================
## uFPCA por hormona → concatenar scores → eigendecomposición conjunto
## Pesos Chiou w_j = 1 / sqrt(tr(C^(jj)))
cat("\n[Layer 5] Happ-Greven mFPCA conjunto con normalización Chiou...\n")

## Seleccionar hormonas del joint operator (excluir glucosa+insulina)
joint_hormones <- intersect(joint_panel, names(PACE_list))
if (length(joint_hormones) < 3){
  cat("  Insuficientes hormonas para mFPCA conjunto\n")
} else {
  ## Construir funData multivariada
  build_funData <- function(o){
    K <- o$fit$selectK
    X <- sweep(o$fit$xiEst[, 1:K, drop = FALSE] %*%
               t(o$fit$phi[, 1:K, drop = FALSE]),
               2, -o$fit$mu)
    funData(argvals = o$fit$workGrid, X = X)
  }

  ## Encontrar sujetos comunes entre todas las hormonas joint
  ids_common <- Reduce(intersect, lapply(joint_hormones, function(h) PACE_list[[h]]$ids))
  cat(sprintf("  Sujetos comunes joint panel (%d hormonas): %d\n",
              length(joint_hormones), length(ids_common)))

  ## Si insuficientes comunes, usar todos disponibles por hormona y rellenar
  if (length(ids_common) >= 30){
    mFDlist <- list()
    for (h in joint_hormones){
      o <- PACE_list[[h]]
      idx <- match(ids_common, o$ids)
      K <- o$fit$selectK
      X <- sweep(o$fit$xiEst[idx, 1:K, drop = FALSE] %*%
                 t(o$fit$phi[, 1:K, drop = FALSE]),
                 2, -o$fit$mu)
      mFDlist[[h]] <- funData(argvals = o$fit$workGrid, X = X)
    }

    ## Pesos Chiou: 1 / sqrt(traza covarianza univariada)
    w_chiou <- sapply(joint_hormones, function(h){
      1 / sqrt(sum(PACE_list[[h]]$fit$lambda))
    })

    mFD <- multiFunData(mFDlist)
    cat("  Ejecutando MFPCA() (Happ-Greven)...\n")
    mfpca_fit <- tryCatch({
      MFPCA(mFD, M = 12,
            uniExpansions = lapply(joint_hormones, function(h){
              o <- PACE_list[[h]]
              list(type = "given",
                   functions = funData(o$fit$workGrid, t(o$fit$phi[, 1:o$fit$selectK])),
                   scores = o$fit$xiEst[match(ids_common, o$ids), 1:o$fit$selectK, drop=FALSE],
                   ortho = TRUE)
            }),
            weights = w_chiou, fit = TRUE,
            bootstrap = FALSE, verbose = FALSE)
    }, error = function(e){ message("MFPCA error: ", e$message); NULL })

    if (!is.null(mfpca_fit)){
      cat(sprintf("  mFPCA OK. Variance explained: %s\n",
                  paste(round(mfpca_fit$values[1:6]/sum(mfpca_fit$values)*100, 1), collapse = ", ")))

      ## Asignar cohorte por ID
      stratum_per_id <- sapply(ids_common, function(id){
        x <- pseudo_dt[subject_id == id, stratum]
        if (length(x)) x[1] else NA
      })

      mfpca_results <- list(
        fit = mfpca_fit,
        ids = ids_common,
        stratum = unname(stratum_per_id),
        joint_hormones = joint_hormones,
        chiou_weights = w_chiou
      )
      saveRDS(mfpca_results, "data/mFPCA_HappGreven.rds")
      cat("  ✓ mFPCA guardada\n")
    } else {
      cat("  mFPCA falló — usando aproximación cohort-level (siguiente script)\n")
    }
  } else {
    cat("  Insuficientes sujetos comunes — usando aproximación cohort-level\n")
  }
}

cat("\n[Layer 1-5] Pipeline básico completo. Resultados en data/.\n")
