## ============================================================
## 18_FDEP_TP_validation.R — Validación empírica FDEP-TP
##
## Calcula:
##  - ρ_INC, ρ_NET, ρ_ANR (descriptores cruzados)
##  - PTP modal por hormona × cohorte
##  - IEP Type assignment por cohorte
##  - Verifica predicciones canónicas del marco A1/A2/A3
##  - Genera figuras DI3, DI4 y tablas T1-T3
## ============================================================
suppressPackageStartupMessages({
  library(data.table); library(tidyverse); library(ggplot2)
  library(viridis); library(cowplot); library(patchwork); library(scales)
})

WD <- here::here("A4_empirical")
setwd(WD); set.seed(20260514)

dt <- readRDS("data/master_long.rds")
mIPD <- readRDS("data/pseudo_IPD_multivariate.rds")
PACE <- readRDS("data/PACE_multivariate_per_hormone.rds")
mfpca <- readRDS("data/mFPCA_broad_HappGreven.rds")
kernels <- readRDS("data/cross_covariance_kernels.rds")

cohort_FDEP <- c("LEAN","OB","T2DM","OB+T2DM","PRE-SG","POST-SG","PRE-RYGBP","POST-RYGBP","POST-CR")
palette <- c(LEAN="#1B9E77", OB="#D95F02", `OB+T2DM`="#7B3294",
             T2DM="#66A61E", `PRE-SG`="#A6761D", `POST-SG`="#E6AB02",
             `PRE-RYGBP`="#999999", `POST-RYGBP`="#E7298A", `POST-CR`="#386CB0")

## ===== Descriptores cruzados ρ_INC, ρ_NET, ρ_ANR =====
compute_rhos <- function(sub){
  pairs <- list(
    rho_INC     = c("ACTIVE GLP-1","Insulin"),
    rho_NET     = c("TOTAL PYY","TOTAL GHRELIN"),
    rho_ANR     = c("ACYLATED GHRELIN","ACTIVE GLP-1"),
    rho_IN      = c("Glucose","Insulin"),
    rho_GIP_INS = c("TOTAL GIP","Insulin"),
    rho_PYY_INS = c("TOTAL PYY","Insulin"),
    rho_GLP_GLU = c("ACTIVE GLP-1","Glucagon")
  )
  out <- list()
  for (nm in names(pairs)){
    p <- pairs[[nm]]
    w <- dcast(sub[hormone %in% p],
               subject_id + time_min ~ hormone, value.var = "value",
               fun.aggregate = mean)
    if (!all(p %in% names(w))) next
    x <- w[[p[1]]]; y <- w[[p[2]]]
    keep <- !is.na(x) & !is.na(y)
    if (sum(keep) < 5) next
    r <- cor(x[keep], y[keep])
    if (nm %in% c("rho_NET","rho_ANR")) r <- -r
    out[[nm]] <- round(r, 3)
  }
  out
}
rho_table <- list()
for (st in cohort_FDEP){
  sub <- mIPD[stratum == st]
  if (nrow(sub) < 50) next
  rho <- compute_rhos(sub)
  if (length(rho) == 0) next
  rho_table[[st]] <- data.table(stratum = st, as.data.table(as.list(rho)))
}
rho_dt <- rbindlist(rho_table, fill = TRUE)
rho_dt[, stratum := factor(stratum, levels = cohort_FDEP)]
setorder(rho_dt, stratum)
fwrite(rho_dt, "data/FDEP_descriptors_final.csv")
cat("Descriptores cruzados ρ por estrato FDEP-TP:\n")
print(rho_dt)

## ===== PTP per hormona-cohorte usando z-score robusto contra LEAN =====
assign_PTP <- function(z, morph_dev = 0){
  if (is.na(z)) return("Missing")
  if (morph_dev == 0){
    if (z <= -2) return("Impaired")
    if (z <= -1) return("Borderline Impaired")
    if (z >= 2)  return("Enhanced")
    if (z >= 1)  return("Borderline Enhanced")
    return("Preserved")
  } else {
    if (z <= -2) return("Altered")
    if (z <= -1) return("Borderline Altered")
    if (z >= 2)  return("Enhanced")
    if (z >= 1)  return("Borderline Enhanced")
    return("Preserved")
  }
}

PTP_table <- list()
for (h in names(PACE)){
  o <- PACE[[h]]
  scores <- o$fit$xiEst[, 1]
  ref_idx <- which(o$cohort == "LEAN")
  if (length(ref_idx) < 5) next
  ref_med <- median(scores[ref_idx])
  ref_mad <- 1.4826 * mad(scores[ref_idx]) + 1e-6
  z_all   <- (scores - ref_med) / ref_mad
  for (st in unique(o$cohort)){
    sub_idx <- which(o$cohort == st)
    if (length(sub_idx) < 5) next
    ptp_st <- sapply(z_all[sub_idx], assign_PTP)
    z_med  <- median(z_all[sub_idx])
    PTP_table[[length(PTP_table)+1]] <- data.table(
      hormone = h, stratum = st, n = length(sub_idx),
      z_med = round(z_med, 3),
      PTP_modal = names(sort(table(ptp_st), decreasing = TRUE))[1],
      pct_Preserved = round(100 * mean(ptp_st == "Preserved"), 1),
      pct_Enhanced  = round(100 * mean(ptp_st %in% c("Enhanced","Borderline Enhanced")), 1),
      pct_Impaired  = round(100 * mean(ptp_st %in% c("Impaired","Borderline Impaired","Blunted")), 1)
    )
  }
}
PTP_dt <- rbindlist(PTP_table, fill = TRUE)
PTP_dt[, stratum := factor(stratum, levels = cohort_FDEP)]
setorder(PTP_dt, hormone, stratum)
fwrite(PTP_dt, "data/PTP_modal_final.csv")
cat("\nPTP modal por hormona × cohorte:\n")
print(PTP_dt)

## ===== IEP Type per cohorte (8-rule precedence) =====
assign_IEP <- function(ptp_vec, glucose_label = "Preserved"){
  non_glu <- ptp_vec[setdiff(names(ptp_vec), "Glucose")]
  if (length(non_glu) < 2) return("Not classifiable")
  glu_suffix <- ifelse(glucose_label == "Preserved", "a",
                ifelse(glucose_label %in% c("Borderline Enhanced","Enhanced"), "b","c"))
  if (any(non_glu == "Altered"))            return("IV·II")
  if (any(non_glu == "Borderline Altered")) return("IV·I")
  if (any(non_glu == "Enhanced")){
    if (glu_suffix == "a") return("II·II")
    else                   return(paste0("V·II·", glu_suffix))
  }
  if (any(non_glu == "Borderline Enhanced")){
    if (glu_suffix == "a") return("II·I")
    else                   return(paste0("V·I·", glu_suffix))
  }
  if (any(non_glu %in% c("Impaired","Blunted")))
    return(paste0("III·II·", glu_suffix))
  if (any(non_glu == "Borderline Impaired"))
    return(paste0("III·I·", glu_suffix))
  if (all(non_glu == "Preserved")) return("I·I")
  if (any(non_glu == "Recovered") && all(non_glu %in% c("Preserved","Recovered")))
    return("I·II")
  return("Not classifiable")
}

IEP_table <- list()
for (st in cohort_FDEP){
  hs <- PTP_dt[stratum == st]
  if (nrow(hs) < 4) next
  vec <- setNames(hs$PTP_modal, hs$hormone)
  glu_label <- if ("Glucose" %in% names(vec)) vec["Glucose"] else "Preserved"
  IEP_table[[st]] <- data.table(
    stratum = st,
    n_hormones = length(vec),
    glucose_PTP = glu_label,
    IEP_Type = assign_IEP(vec, glu_label),
    PTP_summary = paste(names(vec), ":", vec, collapse = "; ")
  )
}
IEP_dt <- rbindlist(IEP_table, fill = TRUE)
IEP_dt[, stratum := factor(stratum, levels = cohort_FDEP)]
setorder(IEP_dt, stratum)
fwrite(IEP_dt, "data/IEP_final.csv")
cat("\nIEP Type assignment (8-rule precedence):\n")
print(IEP_dt[, .(stratum, glucose_PTP, IEP_Type)])

## ===== Verificación de predicciones canónicas A1/A2/A3 =====
## Predicciones FDEP-TP:
##  1. POST-RYGBP modal Type V (43.3 % esperado)
##  2. POST-SG modal Type V (30.7 % esperado)
##  3. OB+T2DM modal Type III (40.0 % esperado)
##  4. POST-CR modal Type III (36.8 % esperado)
##  5. LEAN modal Type I·I (32.2 % esperado)
##  6. ξ₁ post-RYGBP > 0 (PYY-dominant, distal L-cell)
##  7. ξ₄ POST-SG distintivo (ghrelin tone removed)
predictions <- data.table(
  prediction_id = c("POST-RYGBP modal Type V (bariatric)",
                    "POST-SG modal Type V (bariatric)",
                    "OB+T2DM Type III modal",
                    "POST-CR Type III modal",
                    "LEAN Type I·I modal",
                    "ξ₁ POST-RYGBP > LEAN (distal L-cell)",
                    "ξ₄ POST-SG distinctive (ghrelin removed)"),
  framework_value = c("43.3 %","30.7 %","40.0 %","36.8 %","32.2 %","+","strong+")
)

xi_means <- readRDS("data/mFPCA_broad_HappGreven.rds")$fit$scores[, 1:4]
stratum_v <- mfpca$stratum
xi_by <- aggregate(xi_means, by = list(stratum = stratum_v), FUN = mean)

verification <- data.table(
  prediction_id = predictions$prediction_id,
  framework_value = predictions$framework_value,
  empirical_value = c(
    IEP_dt[stratum == "POST-RYGBP", IEP_Type],
    IEP_dt[stratum == "POST-SG", IEP_Type],
    IEP_dt[stratum == "OB+T2DM", IEP_Type],
    IEP_dt[stratum == "POST-CR", IEP_Type],
    IEP_dt[stratum == "LEAN", IEP_Type],
    sprintf("ξ₁ POST-RYGBP=%.2f, LEAN=%.2f",
            xi_by[xi_by$stratum=="POST-RYGBP", 2],
            xi_by[xi_by$stratum=="LEAN", 2]),
    sprintf("ξ₄ POST-SG=%.2f vs LEAN=%.2f",
            xi_by[xi_by$stratum=="POST-SG", 5],
            xi_by[xi_by$stratum=="LEAN", 5])
  )
)
fwrite(verification, "data/FDEP_TP_predictions_verification.csv")
cat("\nVerificación predicciones canónicas FDEP-TP:\n")
print(verification)

## ===== DI3 - Cross-covariance kernels desde pseudo-IPD multivariado =====
times_grid <- c(0,15,30,60,90,120,150,180)
wide_mIPD <- dcast(mIPD[time_min %in% times_grid],
                   subject_id + stratum ~ hormone + time_min,
                   value.var = "value", fun.aggregate = mean)

compute_kernel_v2 <- function(j, k){
  cols_j <- grep(paste0("^", j, "_"), names(wide_mIPD), value = TRUE)
  cols_k <- grep(paste0("^", k, "_"), names(wide_mIPD), value = TRUE)
  if (length(cols_j) == 0 || length(cols_k) == 0) return(NULL)
  Xj <- as.matrix(wide_mIPD[, ..cols_j])
  Xk <- as.matrix(wide_mIPD[, ..cols_k])
  keep <- complete.cases(Xj) & complete.cases(Xk)
  if (sum(keep) < 8) return(NULL)
  Xj_c <- scale(Xj[keep, ], center = TRUE, scale = FALSE)
  Xk_c <- scale(Xk[keep, ], center = TRUE, scale = FALSE)
  Cjk <- crossprod(Xj_c, Xk_c) / (sum(keep) - 1)
  rownames(Cjk) <- gsub(paste0("^", j, "_"), "", cols_j)
  colnames(Cjk) <- gsub(paste0("^", k, "_"), "", cols_k)
  Cjk
}

key_pairs <- list(
  c("ACTIVE GLP-1","Insulin"),
  c("TOTAL GIP","Insulin"),
  c("TOTAL PYY","TOTAL GHRELIN"),
  c("ACYLATED GHRELIN","ACTIVE GLP-1"),
  c("TOTAL PYY","Insulin"),
  c("Glucose","Insulin")
)
heat_list <- list()
for (pr in key_pairs){
  K <- compute_kernel_v2(pr[1], pr[2])
  if (is.null(K)) next
  for (i in 1:nrow(K)) for (j in 1:ncol(K)){
    heat_list[[length(heat_list)+1]] <- data.table(
      time_s = as.numeric(rownames(K)[i]),
      time_t = as.numeric(colnames(K)[j]),
      cov = K[i, j],
      pair = paste(pr, collapse = " ↔ ")
    )
  }
}
heat_dt <- rbindlist(heat_list, fill = TRUE)

p_kernel <- ggplot(heat_dt, aes(time_s, time_t, fill = cov)) +
  geom_raster() +
  geom_text(aes(label = round(cov, 0)), size = 2, color = "black") +
  facet_wrap(~ pair, nrow = 2) +
  scale_fill_distiller(palette = "RdBu", name = expression(C^{(jk)}(s,t))) +
  labs(x = "tiempo s (min, hormona j)", y = "tiempo t (min, hormona k)",
       title = "DI3 — Kernels cross-covariance C^(jk)(s,t) — interacciones entero-pancreáticas",
       subtitle = "Estimados desde la pseudo-IPD multivariada (subjects × hormone × time)") +
  theme_cowplot(11) + theme(legend.position = "right")
ggsave("figures/DI3_cross_covariance_kernels.png", p_kernel, width = 13, height = 8, dpi = 200)

## ===== DI4 - Heatmap descriptores ρ per cohorte =====
rho_long <- melt(rho_dt, id.vars = "stratum",
                 measure.vars = c("rho_INC","rho_NET","rho_ANR","rho_IN","rho_GIP_INS","rho_PYY_INS"),
                 variable.name = "descriptor", value.name = "rho")
rho_long[, descriptor := factor(descriptor,
   levels = c("rho_INC","rho_NET","rho_ANR","rho_IN","rho_GIP_INS","rho_PYY_INS"),
   labels = c("ρ_INC\n(GLP-1↔Ins)","ρ_NET\n(-PYY↔Ghr)","ρ_ANR\n(-Ghr_acyl↔GLP-1)",
              "ρ_IN\n(Glu↔Ins)","ρ_GIP·Ins\n(GIP↔Ins)","ρ_PYY·Ins\n(PYY↔Ins)"))]

p_rho <- ggplot(rho_long, aes(descriptor, stratum, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(rho, 2)), size = 3, color = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = expression(rho),
                       limits = c(-1, 1)) +
  labs(x = NULL, y = NULL,
       title = "DI4 — Descriptores cruzados ρ por estrato FDEP-TP",
       subtitle = "Heatmap de seis acoplamientos funcionales clave por cohorte") +
  theme_cowplot(11) + theme(axis.text.x = element_text(angle = 0))
ggsave("figures/DI4_rho_descriptors.png", p_rho, width = 12, height = 7, dpi = 200)

cat("\n✓ DI3 (cross-covariance kernels) y DI4 (ρ-heatmap) generadas.\n")
cat("\n✓ Validation pipeline FDEP-TP complete.\n")
