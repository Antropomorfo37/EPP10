## ============================================================
## 16_FDEP_TP_descriptors_PTP_IEP.R
##  Layer 6 — Descriptores cruzados ρ_INC, ρ_NET, ρ_ANR por cohorte
##  Layer 7 — PTP 9-label e IEP 8-rule per cohorte
##  Layer 8 — Inferencia: Pillai trace, Mahalanobis, scores cohorte
## ============================================================
suppressPackageStartupMessages({
  library(data.table); library(tidyverse); library(ggplot2)
  library(viridis); library(cowplot); library(patchwork)
})

WD <- here::here("A4_empirical")
setwd(WD); set.seed(20260514)

mIPD <- readRDS("data/pseudo_IPD_multivariate.rds")
PACE <- readRDS("data/PACE_multivariate_per_hormone.rds")
mfpca <- readRDS("data/mFPCA_HappGreven_final.rds")
master <- readRDS("data/master_long.rds")

cohort_FDEP <- c("LEAN","OB","T2DM","OB+T2DM","PRE-SG","POST-SG","PRE-RYGBP","POST-RYGBP","POST-CR")
strata_present <- intersect(cohort_FDEP, unique(mIPD$stratum))

## ===== Layer 6 — Descriptores cruzados =====
## ρ_INC = corr(GLP-1_active, Insulin)  — eje incretínico
## ρ_NET = -corr(PYY_total, Ghrelin_total) — balance anorexigénico
## ρ_ANR = -corr(Ghrelin_acyl, GLP-1_active) — balance grelina-incretina
## ρ_IN  = corr(Glucose, Insulin) — acoplamiento glucémico

compute_rho_cohort <- function(sub){
  ## Construir matriz hormona × tiempo promediada
  ## y luego curvas por sujeto si disponibles
  hormones_pairs <- list(
    rho_INC = c("ACTIVE GLP-1","Insulin"),
    rho_NET = c("TOTAL PYY","TOTAL GHRELIN"),
    rho_ANR = c("ACYLATED GHRELIN","ACTIVE GLP-1"),
    rho_IN  = c("Glucose","Insulin"),
    rho_GIP_INS = c("TOTAL GIP","Insulin")
  )
  out <- list()
  for (nm in names(hormones_pairs)){
    p <- hormones_pairs[[nm]]
    sub_p <- sub[hormone %in% p, .(value = mean(value, na.rm = TRUE)),
                 by = .(subject_id, hormone, time_min)]
    if (nrow(sub_p) == 0) next
    w <- dcast(sub_p, subject_id + time_min ~ hormone, value.var = "value")
    if (!all(p %in% names(w))) next
    x <- w[[p[1]]]; y <- w[[p[2]]]
    keep <- !is.na(x) & !is.na(y)
    if (sum(keep) < 5) next
    r <- cor(x[keep], y[keep], method = "pearson")
    ## Aplicar signo según definición
    if (nm %in% c("rho_NET","rho_ANR")) r <- -r
    out[[nm]] <- r
  }
  out
}

rho_table <- list()
for (st in strata_present){
  sub <- mIPD[stratum == st]
  if (nrow(sub) < 50) next
  rho <- compute_rho_cohort(sub)
  rho_table[[st]] <- data.table(
    stratum = st,
    rho_INC = if (!is.null(rho$rho_INC)) round(rho$rho_INC, 3) else NA,
    rho_NET = if (!is.null(rho$rho_NET)) round(rho$rho_NET, 3) else NA,
    rho_ANR = if (!is.null(rho$rho_ANR)) round(rho$rho_ANR, 3) else NA,
    rho_IN  = if (!is.null(rho$rho_IN))  round(rho$rho_IN, 3) else NA,
    rho_GIP_INS = if (!is.null(rho$rho_GIP_INS)) round(rho$rho_GIP_INS, 3) else NA
  )
}
rho_dt <- rbindlist(rho_table, fill = TRUE)
rho_dt[, stratum := factor(stratum, levels = cohort_FDEP)]
setorder(rho_dt, stratum)
fwrite(rho_dt, "data/FDEP_descriptors.csv")
cat("Descriptores cruzados ρ_INC, ρ_NET, ρ_ANR por estrato FDEP-TP:\n")
print(rho_dt)

## ===== Layer 7 — PTP 9-label per analyte per stratum =====
## Calcular BLUP scores ξ_j por sujeto y compararlos con percentiles del estrato LEAN
##
## PTP labels:
##  Preserved (|z|≤1 leading), Borderline Impaired (z≤-1 SD leading, morf preserved),
##  Impaired (z≤-2 SD, morf preserved), Blunted (amplitude+morph both reduced),
##  Borderline Altered / Altered (morphological deviation), Recovered,
##  Borderline Enhanced (z≥1 SD), Enhanced (z≥2 SD)

assign_PTP <- function(score, leading_threshold = c(-2, -1, 1, 2), morph_dev = 0){
  if (is.na(score)) return("Missing")
  z <- score
  ## Sin morfología detallada, basamos label en z
  if (morph_dev == 0){
    if (z <= leading_threshold[1]) return("Impaired")
    if (z <= leading_threshold[2]) return("Borderline Impaired")
    if (z >= leading_threshold[4]) return("Enhanced")
    if (z >= leading_threshold[3]) return("Borderline Enhanced")
    return("Preserved")
  } else {
    if (z <= leading_threshold[1]) return("Altered")
    if (z <= leading_threshold[2]) return("Borderline Altered")
    return("Preserved")
  }
}

PTP_table <- list()
for (h in names(PACE)){
  o <- PACE[[h]]
  K <- o$fit$selectK
  scores <- o$fit$xiEst[, 1, drop = TRUE]   ## ξ_1 (leading)
  ## Estandarizar contra LEAN
  ref_idx <- which(o$cohort == "LEAN")
  if (length(ref_idx) < 5) next
  ref_med <- median(scores[ref_idx])
  ref_mad <- 1.4826 * mad(scores[ref_idx]) + 1e-6
  z <- (scores - ref_med) / ref_mad
  ## Por sujeto → PTP label
  ptp <- sapply(z, assign_PTP)
  ## Agregar por cohorte: modal label
  for (st in unique(o$cohort)){
    sub_idx <- which(o$cohort == st)
    if (length(sub_idx) < 5) next
    tbl <- table(ptp[sub_idx])
    modal <- names(tbl)[which.max(tbl)]
    PTP_table[[length(PTP_table)+1]] <- data.table(
      hormone = h, stratum = st,
      n = length(sub_idx),
      score_med = round(median(scores[sub_idx]), 3),
      z_med = round(median(z[sub_idx]), 3),
      PTP_modal = modal,
      pct_Preserved = round(100 * mean(ptp[sub_idx] == "Preserved"), 1),
      pct_Enhanced = round(100 * mean(ptp[sub_idx] %in% c("Borderline Enhanced","Enhanced")), 1),
      pct_Impaired = round(100 * mean(ptp[sub_idx] %in% c("Borderline Impaired","Impaired")), 1)
    )
  }
}
PTP_dt <- rbindlist(PTP_table, fill = TRUE)
PTP_dt[, stratum := factor(stratum, levels = cohort_FDEP)]
setorder(PTP_dt, hormone, stratum)
fwrite(PTP_dt, "data/PTP_modal_per_cohort.csv")
cat("\nPTP modal por (hormona × cohorte):\n")
print(PTP_dt[1:40])

## ===== IEP 8-rule precedence =====
## (1) any non-glucose Altered → Type IV·II
## (2) any non-glucose Borderline Altered → Type IV·I
## (3a) any non-glucose Enhanced + glucose suffix a (Preserved) and no Altered/BorderlineAltered → Type II·II
## (3b) any non-glucose Enhanced + glucose suffix b/c → Type V·II
## (4a) any non-glucose Borderline Enhanced + glucose suffix a → Type II·I
## (4b) any non-glucose Borderline Enhanced + glucose suffix b/c → Type V·I
## (5) any non-glucose Impaired or Blunted without (1)-(4) → Type III·II
## (6) any non-glucose Borderline Impaired without (1)-(5) → Type III·I
## (7) all non-glucose Preserved → Type I·I
## (8) at least one Recovered with remainder Preserved or Recovered → Type I·II

assign_IEP <- function(ptp_vec, glucose_label = "Preserved"){
  non_glu <- ptp_vec[names(ptp_vec) != "Glucose"]
  if (length(non_glu) < 2) return("Not classifiable")
  glu_suffix <- if (glucose_label == "Preserved") "a" else if (glucose_label %in% c("Borderline Enhanced","Enhanced")) "b" else "c"

  if (any(non_glu == "Altered")) return(paste0("IV·II"))
  if (any(non_glu == "Borderline Altered")) return("IV·I")
  if (any(non_glu == "Enhanced")){
    if (glu_suffix == "a") return("II·II")
    else return(paste0("V·II·", glu_suffix))
  }
  if (any(non_glu == "Borderline Enhanced")){
    if (glu_suffix == "a") return("II·I")
    else return(paste0("V·I·", glu_suffix))
  }
  if (any(non_glu %in% c("Impaired","Blunted"))) return(paste0("III·II·", glu_suffix))
  if (any(non_glu == "Borderline Impaired")) return(paste0("III·I·", glu_suffix))
  if (all(non_glu == "Preserved")) return("I·I")
  if (any(non_glu == "Recovered") &&
      all(non_glu %in% c("Preserved","Recovered"))) return("I·II")
  return("Not classifiable")
}

IEP_table <- list()
for (st in strata_present){
  hormones_st <- PTP_dt[stratum == st, .(hormone, PTP_modal)]
  if (nrow(hormones_st) < 4) next
  vec <- setNames(hormones_st$PTP_modal, hormones_st$hormone)
  glu_label <- if ("Glucose" %in% names(vec)) vec["Glucose"] else "Preserved"
  iep <- assign_IEP(vec, glu_label)
  IEP_table[[st]] <- data.table(
    stratum = st,
    n_hormones = length(vec),
    glucose_PTP = glu_label,
    IEP_Type = iep,
    PTPs_present = paste(names(vec), vec, sep = "=", collapse = "; ")
  )
}
IEP_dt <- rbindlist(IEP_table, fill = TRUE)
IEP_dt[, stratum := factor(stratum, levels = cohort_FDEP)]
setorder(IEP_dt, stratum)
fwrite(IEP_dt, "data/IEP_per_cohort.csv")
cat("\nIEP Type assignment per cohort (8-rule precedence):\n")
print(IEP_dt[, .(stratum, glucose_PTP, IEP_Type)])

## ===== Layer 8 — Pillai trace, Mahalanobis =====
## Score panel Ξ del mFPCA: n × M
mfpca_fit <- mfpca$fit
xi_panel <- mfpca_fit$scores      ## n × M
xi_top4 <- xi_panel[, 1:4]
stratum_vec <- mfpca$stratum

## Sólo cohortes con n ≥ 25 para Pillai
n_per_st <- table(stratum_vec)
strata_pillai <- names(n_per_st)[n_per_st >= 25]
mask <- stratum_vec %in% strata_pillai

df_xi <- data.frame(xi_top4 = xi_top4[mask, ], stratum = factor(stratum_vec[mask]))
names(df_xi)[1:4] <- c("xi1","xi2","xi3","xi4")
mlm <- manova(cbind(xi1, xi2, xi3, xi4) ~ stratum, data = df_xi)
pillai <- summary(mlm, test = "Pillai")
print(pillai)

## Mahalanobis distance: cohort centroid vs LEAN
centroids <- aggregate(xi_top4, by = list(stratum = stratum_vec), FUN = mean)
ref <- as.numeric(centroids[centroids$stratum == "LEAN", -1])
cov_within <- cov(xi_top4)
inv_cov <- MASS::ginv(cov_within)
maha <- apply(centroids[, -1], 1, function(x){
  d <- as.numeric(x) - ref
  sqrt(sum(d * (inv_cov %*% d)))
})
maha_dt <- data.table(stratum = centroids$stratum, mahalanobis_vs_LEAN = round(maha, 3))
maha_dt[, stratum := factor(stratum, levels = cohort_FDEP)]
setorder(maha_dt, stratum)
fwrite(maha_dt, "data/mahalanobis_vs_LEAN.csv")
cat("\nDistancia Mahalanobis cohorte-LEAN sobre (ξ1,ξ2,ξ3,ξ4):\n")
print(maha_dt)

## ===== ξ scores promedio por cohorte =====
xi_summary <- aggregate(xi_top4, by = list(stratum = stratum_vec),
                        FUN = function(x) c(mean = mean(x), sd = sd(x)))
xi_means <- aggregate(xi_top4, by = list(stratum = stratum_vec), FUN = mean)
names(xi_means) <- c("stratum","xi1_mean","xi2_mean","xi3_mean","xi4_mean")
xi_means$stratum <- factor(xi_means$stratum, levels = cohort_FDEP)
xi_means <- xi_means[order(xi_means$stratum), ]
fwrite(xi_means, "data/xi_means_per_cohort.csv")
cat("\nMedios ξ_1, ξ_2, ξ_3, ξ_4 por cohorte (referenciados a LEAN):\n")
print(xi_means)

## ===== Save Pillai output =====
saveRDS(list(rho = rho_dt, PTP = PTP_dt, IEP = IEP_dt,
             pillai = pillai, mahalanobis = maha_dt,
             xi_means = xi_means),
        "data/FDEP_TP_layer678.rds")

cat("\n✓ Layers 6-8 completos. Resultados en data/.\n")
