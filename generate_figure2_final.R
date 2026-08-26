# =============================================================================
# generate_figure2_final.R — Figura 2 con 4 paneles (joint + IEP + cross-val)
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
  library(ggplot2); library(scales); library(forcats); library(patchwork)
})

theme_jcem <- function() {
  theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          panel.border = element_rect(color="black", size=0.5),
          strip.background = element_rect(fill="grey92", color=NA),
          legend.position = "bottom",
          plot.title = element_text(size=11, face="bold"),
          axis.title = element_text(size=10))
}

cohort_order <- c("no_obese_without_T2DM","Obesity","T2DM",
                  "Obesity_T2DM","SG","RYGBP")

boot   <- readRDS(epp10_path("bootstrap_stability_results.rds"))
fanova <- readRDS(epp10_path("fanova_results.rds"))
iep    <- read_csv(epp10_path("iep_frequency_by_cohort.csv"),
                   show_col_types = FALSE)

# Panel A: joint-mFPC prevalence (B=50 medians)
class_order <- c("Preservado","Impairment_limitrofe","Impaired",
                 "Blunted","Enhanced","Altered")
class_colors <- c(Preservado="#2ECC71", Impairment_limitrofe="#58D68D",
                  Impaired="#F1C40F", Blunted="#3498DB",
                  Enhanced="#E74C3C", Altered="#8E44AD")
prev_df <- boot$dist_summary %>%
  mutate(cohort = factor(cohort, levels = cohort_order),
         cls = factor(cls, levels = class_order))
panel_A <- ggplot(prev_df, aes(cohort, median_pct, fill = cls)) +
  geom_col(position = "stack", color = "white", size = 0.3, width = 0.75) +
  scale_fill_manual(values = class_colors, name = NULL,
                    guide = guide_legend(nrow = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(title = "A — Joint-mFPC 6-class prevalence (B=50 pipeline)",
       x = NULL, y = "% pseudo-subjects") +
  theme_jcem() + theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 8),
                        legend.text = element_text(size = 8))

# Panel B: Pillai F lollipop
pillai_df <- fanova$pillai_pairwise %>%
  mutate(cohort = sub("_vs_ref", "", as.character(contrast))) %>%
  mutate(cohort = factor(cohort, levels = setdiff(cohort_order, "no_obese_without_T2DM"))) %>%
  arrange(desc(F_obs))
panel_B <- ggplot(pillai_df, aes(fct_reorder(cohort, F_obs), F_obs)) +
  geom_segment(aes(xend = cohort, yend = 0), color = "grey50") +
  geom_point(size = 4, color = "steelblue") +
  geom_text(aes(label = sprintf("F=%.0f", F_obs)), hjust = -0.3, size = 3) +
  coord_flip() +
  scale_y_continuous(limits = c(0, max(pillai_df$F_obs) * 1.2),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(title = "B — Pillai F (pairwise vs. no_obese_without_T2DM)",
       x = NULL, y = "Pillai approx. F (FDR<0.05 all)") +
  theme_jcem()

# Panel C: IEP Type I-V stacked bars
iep_long <- iep %>%
  pivot_longer(-cohort, names_to = "type", values_to = "pct") %>%
  mutate(cohort = factor(cohort, levels = cohort_order),
         type = factor(type,
           levels = c("I.I","I.II","II.I","II.II","III.I","III.II",
                      "IV.I","IV.II","V.I","V.II","not_integrable")))
iep_colors <- c("I.I"="#27AE60","I.II"="#2ECC71",
                "II.I"="#F39C12","II.II"="#E67E22",
                "III.I"="#F4D03F","III.II"="#F1C40F",
                "IV.I"="#9B59B6","IV.II"="#6C3483",
                "V.I"="#EC7063","V.II"="#CB4335",
                "not_integrable"="#95A5A6")
panel_C <- ggplot(iep_long %>% filter(!is.na(pct)),
                  aes(cohort, pct, fill = type)) +
  geom_col(position = "stack", color = "white", size = 0.3, width = 0.75) +
  scale_fill_manual(values = iep_colors, name = NULL,
                    guide = guide_legend(nrow = 2)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(title = "C — PTP/IEP Type I–V per cohort (framework v1.0, Apr 2026)",
       x = NULL, y = "% pseudo-subjects") +
  theme_jcem() + theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 8),
                        legend.text = element_text(size = 7))

# Panel D: cross-validation joint Pillai F vs IEP IV.II %
cross_tbl <- pillai_df %>% rename(cohort_study = cohort) %>%
  left_join(iep %>% select(cohort, IV.II) %>% rename(cohort_study = cohort),
            by = "cohort_study") %>%
  filter(!is.na(IV.II))
rho_spear <- round(cor(cross_tbl$F_obs, cross_tbl$IV.II, method = "spearman"), 2)
panel_D <- ggplot(cross_tbl, aes(IV.II, F_obs, label = cohort_study)) +
  geom_point(size = 5, color = "darkorange") +
  ggrepel::geom_text_repel(size = 3.2, box.padding = 0.5, segment.alpha = 0.5) +
  labs(title = sprintf("D — Cross-validation: Pillai F vs IEP Type IV.II (ρ_Spearman = %.2f)",
                        rho_spear),
       x = "% pseudo-subjects in IEP Type IV.II",
       y = "Pillai F (cohort vs. ref)") +
  theme_jcem()

# Compose 2x2
fig2 <- (panel_A | panel_C) / (panel_B | panel_D)
ggsave(epp10_path("figures/Figure2_classification_inference.pdf"),
       fig2, width = 12, height = 9, device = cairo_pdf)
ggsave(epp10_path("figures/Figure2_classification_inference.png"),
       fig2, width = 12, height = 9, dpi = 300)

cat("Figure 2 updated: 4 panels (joint 6-class | IEP Type I-V | Pillai F | cross-val)\n")
