## ============================================================
## 17_FDEP_TP_mFPCA_broad.R
##
## Reduce joint panel to 4 hormones with broad cohort coverage
## (TOTAL GIP, ACTIVE GLP-1, TOTAL GLP-1, Glucagon) to enable
## Pillai/Mahalanobis across all 8 strata. The 4-hormone joint
## operator is biologically meaningful (incretin yield + counter-
## regulation), and matches the canonical 4-eigenfunction reading.
## ============================================================
suppressPackageStartupMessages({
  library(data.table); library(tidyverse); library(ggplot2)
  library(fdapace); library(MFPCA); library(funData); library(MASS)
  library(cowplot); library(patchwork); library(viridis)
})

WD <- here::here("A4_empirical")
setwd(WD); set.seed(20260514)

mIPD <- readRDS("data/pseudo_IPD_multivariate.rds")
PACE <- readRDS("data/PACE_multivariate_per_hormone.rds")
cohort_FDEP <- c("LEAN","OB","T2DM","OB+T2DM","PRE-SG","POST-SG","PRE-RYGBP","POST-RYGBP","POST-CR")
palette <- c(LEAN="#1B9E77", OB="#D95F02", `OB+T2DM`="#7B3294",
             T2DM="#66A61E", `PRE-SG`="#A6761D", `POST-SG`="#E6AB02",
             `PRE-RYGBP`="#999999", `POST-RYGBP`="#E7298A", `POST-CR`="#386CB0")

## Broad-coverage joint panel
joint_broad <- c("TOTAL GIP","ACTIVE GLP-1","TOTAL GLP-1","Glucagon")

## ===== Refit pseudo-IPD with strata that have ALL 4 broad-panel hormones =====
strata_with_all <- Reduce(intersect, lapply(joint_broad, function(h) unique(PACE[[h]]$cohort)))
cat("Estratos con todas las 4 hormonas broad-panel:\n")
print(strata_with_all)

ids_common <- Reduce(intersect, lapply(joint_broad, function(h){
  o <- PACE[[h]]; o$ids[o$cohort %in% strata_with_all]
}))
cat(sprintf("Sujetos comunes (4-hormone broad panel): %d\n", length(ids_common)))

## ===== Build mFPCA Happ-Greven on broad panel =====
w_chiou <- sapply(joint_broad, function(h) 1 / sqrt(sum(PACE[[h]]$fit$lambda)))
cat("Chiou weights:\n"); print(round(w_chiou, 3))

uniExp <- lapply(joint_broad, function(h){
  o <- PACE[[h]]
  idx <- match(ids_common, o$ids)
  K <- o$fit$selectK
  list(type = "given",
       functions = funData(o$fit$workGrid, t(o$fit$phi[, 1:K])),
       scores = o$fit$xiEst[idx, 1:K, drop = FALSE],
       ortho = TRUE)
})

mFDlist <- lapply(joint_broad, function(h){
  o <- PACE[[h]]
  idx <- match(ids_common, o$ids)
  K <- o$fit$selectK
  X <- o$fit$xiEst[idx, 1:K, drop = FALSE] %*% t(o$fit$phi[, 1:K, drop = FALSE])
  X <- sweep(X, 2, o$fit$mu, "+")
  funData(argvals = o$fit$workGrid, X = X)
})
mFD <- multiFunData(mFDlist)

cat("\nEjecutando Happ-Greven MFPCA (M=12)...\n")
mfpca_fit <- MFPCA(mFD, M = 12, uniExpansions = uniExp,
                   weights = w_chiou, fit = TRUE,
                   bootstrap = TRUE, nBootstrap = 200,
                   bootstrapAlpha = c(0.05), verbose = FALSE)

cumFVE <- cumsum(mfpca_fit$values) / sum(mfpca_fit$values)
cat(sprintf("\nFVE per eigenvalue top-6: %s\n",
            paste(round(mfpca_fit$values[1:6]/sum(mfpca_fit$values)*100, 1), collapse = ", ")))
cat(sprintf("Cumulative FVE top-4: %s\n",
            paste(round(cumFVE[1:4]*100, 1), collapse = ", ")))

stratum_per_id <- sapply(ids_common, function(id){
  x <- mIPD[subject_id == id, stratum]
  if (length(x)) x[1] else NA
})
cat("\nEstratos en el operador conjunto:\n")
print(table(stratum_per_id))

mfpca_pkg <- list(fit = mfpca_fit, ids = ids_common,
                  stratum = unname(stratum_per_id),
                  joint_hormones = joint_broad,
                  chiou_weights = w_chiou,
                  cumFVE = cumFVE)
saveRDS(mfpca_pkg, "data/mFPCA_broad_HappGreven.rds")

## ===== Pillai trace y Mahalanobis sobre los 4 ejes mFPC =====
xi_top4 <- mfpca_fit$scores[, 1:4]
df_xi <- data.frame(xi1 = xi_top4[,1], xi2 = xi_top4[,2],
                    xi3 = xi_top4[,3], xi4 = xi_top4[,4],
                    stratum = factor(stratum_per_id, levels = cohort_FDEP))

## Quitar estratos con n < 25
n_per_st <- table(df_xi$stratum)
strata_pillai <- names(n_per_st)[n_per_st >= 25]
mask <- df_xi$stratum %in% strata_pillai
df_pillai <- df_xi[mask, ]
df_pillai$stratum <- droplevels(df_pillai$stratum)

mlm <- manova(cbind(xi1, xi2, xi3, xi4) ~ stratum, data = df_pillai)
pillai_summary <- summary(mlm, test = "Pillai")
print(pillai_summary)
F_omnibus <- pillai_summary$stats[1, "approx F"]
cat(sprintf("\nPillai omnibus F = %.2f\n", F_omnibus))

## Pairwise contrasts vs LEAN
pairwise_F <- list()
for (st in setdiff(unique(as.character(df_pillai$stratum)), "LEAN")){
  sub <- df_pillai[df_pillai$stratum %in% c("LEAN", st), ]
  sub$stratum <- droplevels(sub$stratum)
  if (length(unique(sub$stratum)) < 2) next
  mlm_p <- manova(cbind(xi1, xi2, xi3, xi4) ~ stratum, data = sub)
  s <- summary(mlm_p, test = "Pillai")
  pairwise_F[[st]] <- data.table(stratum_vs_LEAN = st,
                                  Pillai_F = round(s$stats[1, "approx F"], 2),
                                  p_value = signif(s$stats[1, "Pr(>F)"], 3))
}
pairwise_dt <- rbindlist(pairwise_F)
pairwise_dt[, q_bonferroni := pmin(p_value * nrow(pairwise_dt), 1)]
fwrite(pairwise_dt, "data/pillai_pairwise_vs_LEAN.csv")
cat("\nPillai pairwise F vs LEAN cohort:\n")
print(pairwise_dt)

## Mahalanobis cohort vs LEAN sobre (ξ1..ξ4)
cov_within <- cov(xi_top4)
inv_cov <- ginv(cov_within)
centroids <- aggregate(xi_top4, by = list(stratum = stratum_per_id), FUN = mean)
ref <- as.numeric(centroids[centroids$stratum == "LEAN", -1])
maha <- apply(centroids[, -1], 1, function(x){
  d <- as.numeric(x) - ref
  sqrt(sum(d * (inv_cov %*% d)))
})
maha_dt <- data.table(stratum = centroids$stratum, mahalanobis = round(maha, 3))
maha_dt[, stratum := factor(stratum, levels = cohort_FDEP)]
setorder(maha_dt, stratum)
fwrite(maha_dt, "data/mahalanobis_broad.csv")
cat("\nDistancia Mahalanobis cohorte-LEAN sobre (ξ1,ξ2,ξ3,ξ4):\n")
print(maha_dt)

## ===== Mean scores por cohorte =====
xi_means <- aggregate(xi_top4, by = list(stratum = stratum_per_id), FUN = mean)
xi_sds   <- aggregate(xi_top4, by = list(stratum = stratum_per_id), FUN = sd)
names(xi_means) <- c("stratum","xi1","xi2","xi3","xi4")
names(xi_sds)   <- c("stratum","xi1_sd","xi2_sd","xi3_sd","xi4_sd")
xi_summary <- merge(xi_means, xi_sds)
xi_summary$stratum <- factor(xi_summary$stratum, levels = cohort_FDEP)
xi_summary <- xi_summary[order(xi_summary$stratum), ]
fwrite(xi_summary, "data/xi_broad_means.csv")
cat("\nMedias ξ_1..ξ_4 por cohorte (broad 4-hormone panel):\n")
print(xi_summary)

## ===== Figura DI1 — Scores ξ por cohorte =====
xi_long <- melt(setDT(xi_summary)[, .(stratum, xi1, xi2, xi3, xi4, xi1_sd, xi2_sd, xi3_sd, xi4_sd)],
                id.vars = "stratum",
                measure.vars = list(c("xi1","xi2","xi3","xi4"),
                                    c("xi1_sd","xi2_sd","xi3_sd","xi4_sd")),
                value.name = c("score","sd"),
                variable.name = "axis")
xi_long[, axis := factor(axis, labels = c("ξ₁ (Distal L-cell)",
                                           "ξ₂ (Proximal-Distal seq.)",
                                           "ξ₃ (Biphasic GLU-INS coupling)",
                                           "ξ₄ (Ghrelin/Counter-reg. tone)"))]

p_scores <- ggplot(xi_long, aes(stratum, score, fill = stratum)) +
  geom_col(alpha = 0.85, color = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = score - sd, ymax = score + sd), width = 0.3, linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = palette, guide = "none") +
  facet_wrap(~ axis, scales = "free_y") +
  labs(x = NULL, y = "Score ξ (broad 4-hormone panel)",
       title = "DI1 — Proyección cohorte sobre los cuatro ejes mFPC canónicos",
       subtitle = "Operador Happ-Greven (peso Chiou); panel: TOTAL GIP, ACTIVE GLP-1, TOTAL GLP-1, Glucagon") +
  theme_cowplot(11) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("figures/DI1_xi_scores.png", p_scores, width = 13, height = 8, dpi = 200)

## ===== Figura DI2 — Eigenfunctions Ψ_1..Ψ_4 =====
eigen_dt <- list()
times <- mfpca_fit$functions[[1]]@argvals[[1]]
for (m in 1:4){
  for (j in seq_along(joint_broad)){
    h_name <- joint_broad[j]
    eigen_vals <- mfpca_fit$functions[[j]]@X[m, ]
    eigen_dt[[length(eigen_dt)+1]] <- data.table(
      eigenfn = paste0("Ψ_", m),
      hormone = h_name,
      time_min = times,
      psi = eigen_vals
    )
  }
}
eigen_dt <- rbindlist(eigen_dt)
eigen_dt[, hormone := factor(hormone, levels = joint_broad)]

p_eigen <- ggplot(eigen_dt, aes(time_min, psi, color = hormone)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 1) +
  facet_wrap(~ eigenfn, scales = "free_y", ncol = 2) +
  scale_color_brewer(palette = "Set1", name = "Hormona") +
  labs(x = "Tiempo postprandial (min)", y = expression(Psi[m](t)),
       title = "DI2 — Cuatro eigenfunciones canónicas Ψ_1..Ψ_4 del operador conjunto Happ–Greven",
       subtitle = "Lectura fisiológica: distal L-cell dominance / proximal-distal sequencing / biphasic coupling / ghrelin-counterreg tone") +
  theme_cowplot(11) + theme(legend.position = "bottom")
ggsave("figures/DI2_eigenfunctions.png", p_eigen, width = 13, height = 9, dpi = 200)

cat("\n✓ DI1 (scores), DI2 (eigenfunctions) generadas.\n")
