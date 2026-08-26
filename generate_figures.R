# =============================================================================
# generate_figures.R — Figuras 1-3 para el manuscrito medRxiv/JCEM
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
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", size = 0.5),
      strip.background = element_rect(fill = "grey92", color = NA),
      legend.position = "bottom",
      plot.title = element_text(size = 11, face = "bold"),
      axis.title = element_text(size = 10)
    )
}

OUT_DIR <- epp10_path("figures")
dir.create(OUT_DIR, showWarnings = FALSE)

# =============================================================================
# FIGURE 1 — Cohort composition + data-flow bars
# =============================================================================
cat("Figure 1: cohort composition + data flow\n")

map_df <- read_csv(epp10_path("cohort_normalization_map.csv"), show_col_types = FALSE)
long_df <- read_csv(epp10_path("hormones_long_tidy.csv"), show_col_types = FALSE)

cohort_order <- c("no_obese_without_T2DM", "Obesity", "T2DM",
                  "Obesity_T2DM", "SG", "RYGBP")

# Panel A: source labels per canonical cohort
panel_a <- map_df %>%
  mutate(cohort_v10_primary = factor(cohort_v10_primary, levels = cohort_order)) %>%
  count(cohort_v10_primary, name = "n_labels") %>%
  ggplot(aes(cohort_v10_primary, n_labels)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = n_labels), vjust = -0.3, size = 3.2) +
  labs(title = "A — Source labels per canonical cohort",
       x = NULL, y = "Distinct source labels") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_jcem() + theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

# Panel B: observations per cohort × hormone (heatmap)
coverage <- long_df %>%
  mutate(cohort_v10_primary = factor(cohort_v10_primary, levels = cohort_order),
         hormone_name = factor(hormone_name,
           levels = c("glucose","insulin","glucagon","GIP_total","GIP_active",
                      "GLP1_total","GLP1_active","PYY_total","PYY_3_36",
                      "ghrelin_total","ghrelin_acyl"))) %>%
  count(cohort_v10_primary, hormone_name, name = "n_obs")
panel_b <- ggplot(coverage, aes(cohort_v10_primary, hormone_name, fill = n_obs)) +
  geom_tile(color = "white", size = 0.3) +
  geom_text(aes(label = n_obs), size = 2.7) +
  scale_fill_viridis_c(option = "D", trans = "sqrt", name = "n obs") +
  labs(title = "B — Observations per cohort × hormone",
       x = NULL, y = NULL) +
  theme_jcem() + theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

# Panel C: timepoint distribution
panel_c <- long_df %>%
  count(time_min, name = "n") %>%
  ggplot(aes(time_min, n)) +
  geom_col(width = 6, fill = "darkorange", alpha = 0.8) +
  geom_text(aes(label = n), vjust = -0.3, size = 3) +
  labs(title = "C — Temporal sampling density",
       x = "Time since meal (min)", y = "n observations") +
  scale_x_continuous(breaks = c(0, 15, 30, 60, 90, 120, 150, 180)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_jcem()

fig1 <- (panel_a | panel_c) / panel_b + plot_layout(heights = c(1, 1.5))
ggsave(file.path(OUT_DIR, "Figure1_cohort_composition.pdf"), fig1,
       width = 9, height = 7, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Figure1_cohort_composition.png"), fig1,
       width = 9, height = 7, dpi = 300)

# =============================================================================
# FIGURE 2 — Primary classification prevalence + Pillai F
# =============================================================================
cat("Figure 2: classification prevalence + Pillai F\n")

boot <- readRDS(epp10_path("bootstrap_stability_results.rds"))
fanova <- readRDS(epp10_path("fanova_results.rds"))

class_order <- c("Preservado","Impairment_limitrofe","Impaired",
                 "Blunted","Enhanced","Altered")
class_colors <- c(Preservado="#2ECC71", Impairment_limitrofe="#58D68D",
                  Impaired="#F1C40F", Blunted="#3498DB",
                  Enhanced="#E74C3C", Altered="#8E44AD")

# Panel A: stacked bars of primary prevalence (B=50 medians)
prev_df <- boot$dist_summary %>%
  mutate(cohort = factor(cohort, levels = cohort_order),
         cls = factor(cls, levels = class_order))
panel_2a <- ggplot(prev_df, aes(cohort, median_pct, fill = cls)) +
  geom_col(position = "stack", color = "white", size = 0.3, width = 0.75) +
  scale_fill_manual(values = class_colors, name = NULL,
                    guide = guide_legend(nrow = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(title = "A — Classification prevalence per cohort (B=50 medians)",
       x = NULL, y = "% pseudo-subjects in class") +
  theme_jcem() + theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 8))

# Panel B: Pillai F pairwise vs reference (lollipop)
pillai_df <- fanova$pillai_pairwise %>%
  mutate(cohort = sub("_vs_ref", "", as.character(contrast))) %>%
  mutate(cohort = factor(cohort, levels = setdiff(cohort_order, "no_obese_without_T2DM"))) %>%
  arrange(desc(F_obs))
panel_2b <- ggplot(pillai_df, aes(fct_reorder(cohort, F_obs), F_obs)) +
  geom_segment(aes(xend = cohort, yend = 0), color = "grey50") +
  geom_point(size = 4, color = "steelblue") +
  geom_text(aes(label = sprintf("F=%.0f", F_obs)),
            hjust = -0.3, size = 3) +
  coord_flip() +
  scale_y_continuous(limits = c(0, max(pillai_df$F_obs) * 1.2),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(title = "B — Pillai F, pairwise vs. no_obese_without_T2DM",
       x = NULL, y = "Pillai approx. F (FDR < 0.05 for all contrasts)") +
  theme_jcem()

# Panel C: per-PC F tests heatmap (top 10 PCs per contrast)
per_pc_df <- fanova$all_results %>%
  filter(test == "F_univariate", contrast != "omnibus") %>%
  mutate(contrast = sub("_vs_ref", "", as.character(contrast)),
         contrast = factor(contrast, levels = setdiff(cohort_order, "no_obese_without_T2DM")),
         pc_num = as.integer(sub("xi", "", pc)),
         pc = factor(pc, levels = paste0("xi", 1:14)))
panel_2c <- per_pc_df %>%
  filter(pc_num <= 14) %>%
  ggplot(aes(pc, contrast, fill = log10(F_obs))) +
  geom_tile(color = "white", size = 0.3) +
  geom_text(aes(label = round(F_obs)),
            size = 2.5,
            color = ifelse(log10(per_pc_df$F_obs[per_pc_df$pc_num <= 14]) > 1.8, "white", "black")) +
  scale_fill_viridis_c(option = "C", name = "log10(F)") +
  labs(title = "C — F statistic per mFPC component",
       x = "mFPC component", y = NULL) +
  theme_jcem() + theme(axis.text.x = element_text(size = 7))

fig2 <- panel_2a / (panel_2b | panel_2c) + plot_layout(heights = c(1, 1))
ggsave(file.path(OUT_DIR, "Figure2_classification_inference.pdf"), fig2,
       width = 10, height = 8, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Figure2_classification_inference.png"), fig2,
       width = 10, height = 8, dpi = 300)

# =============================================================================
# FIGURE 3 — Trajectory differences with sup-t simultaneous bands
# =============================================================================
cat("Figure 3: trajectory differences + sup-t bands\n")

bands <- read_csv(epp10_path("bands_simultaneous.csv"), show_col_types = FALSE)

# Select 6 hormones for the 2×3 panel: those with strongest cohort-discriminating signal
focus_hormones <- c("ghrelin_total", "GLP1_total", "GIP_total",
                     "PYY_total", "insulin", "glucose")
bands_focus <- bands %>%
  filter(hormone_name %in% focus_hormones) %>%
  mutate(hormone_name = factor(hormone_name, levels = focus_hormones),
         cohort = factor(cohort,
           levels = c("Obesity","T2DM","Obesity_T2DM","SG","RYGBP")))
cohort_colors <- c(Obesity="#5DADE2", T2DM="#F5B041",
                   Obesity_T2DM="#E67E22", SG="#58D68D",
                   RYGBP="#C0392B")

fig3 <- ggplot(bands_focus, aes(t, mean_diff, color = cohort, fill = cohort)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, color = NA) +
  geom_line(size = 0.7) +
  scale_color_manual(values = cohort_colors) +
  scale_fill_manual(values = cohort_colors) +
  facet_wrap(~ hormone_name, scales = "free_y", ncol = 3) +
  labs(title = "μ_cohort(t) − μ_no_obese_without_T2DM(t), 95% simultaneous sup-t bands (B=2000)",
       x = "Time since meal (min)",
       y = "Mean difference from reference",
       color = NULL, fill = NULL) +
  theme_jcem() +
  theme(strip.text = element_text(face = "italic", size = 9),
        legend.position = "top")

ggsave(file.path(OUT_DIR, "Figure3_trajectory_bands.pdf"), fig3,
       width = 11, height = 6.5, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Figure3_trajectory_bands.png"), fig3,
       width = 11, height = 6.5, dpi = 300)

cat("\n=== FIGURES GENERATED ===\n")
cat(sprintf("Output directory: %s\n", OUT_DIR))
for (f in list.files(OUT_DIR, pattern = "Figure")) {
  sz <- round(file.size(file.path(OUT_DIR, f)) / 1024, 1)
  cat(sprintf("  %-50s %s KB\n", f, sz))
}
