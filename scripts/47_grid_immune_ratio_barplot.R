#!/usr/bin/env Rscript
# =============================================================
# Script: 47_grid_immune_ratio_barplot.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-06
#
# Goal:
#   Barplot — proportion de bins immune-enrichis parmi les
#   bins immune-positifs, par condition.
#   Source : grid_enriched_zones_summary.csv
#            (produit par 24_spatial_grid_100um_global.R)
#
# Output : outputs/banksy/inflammatory_niche_step2_grid_100um_global/
#   fig_grid_immune_ratio_barplot.pdf/jpg
# =============================================================

set.seed(1997)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(Cairo)
})

setwd(normalizePath("."))
source("scripts/00_palette.R")

# =============================================================
# 1. Paramètres
# =============================================================

OUT_DIR  <- file.path("outputs", "banksy",
                      "inflammatory_niche_step2_grid_100um_global")
CSV_IN   <- file.path(OUT_DIR, "grid_enriched_zones_summary.csv")

SAMPLE_ORDER  <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")
SAMPLE_LABELS <- c(
  LCMV_1wpi = "LCMV 1 wpi",
  LCMV_3wpi = "LCMV 3 wpi",
  LCMV_6wpi = "LCMV 6 wpi"
)

# Couleurs des conditions depuis GLOBAL_PALETTE
# (Immune (Acod1) = rouge, sinon fallback séquentiel)
COND_COLORS <- c(
  "LCMV 1 wpi"  = "#FDBB84",
  "LCMV 3 wpi"  = "#E34A33",
  "LCMV 6 wpi"  = "#8B0000"
)

REF_LINE <- 0.25   # ligne de référence pointillée

# =============================================================
# 2. Chargement des données
# =============================================================

message("\n=== Chargement : ", CSV_IN)
stopifnot(file.exists(CSV_IN))

df <- read.csv(CSV_IN, stringsAsFactors = FALSE) %>%
  filter(sample %in% SAMPLE_ORDER) %>%
  mutate(
    condition = factor(SAMPLE_LABELS[sample], levels = SAMPLE_LABELS),
    ratio     = ratio_immune_enrich_vs_pos,
    bar_label = paste0(n_immune_enrich_zones, "/", n_immune_pos_zones)
  )

message("  Lignes : ", nrow(df))
print(df %>% select(sample, condition, ratio, bar_label))

# =============================================================
# 3. Figure
# =============================================================

message("\n=== Génération de la figure ===")

p <- ggplot(df, aes(x = condition, y = ratio, fill = condition)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.3) +
  geom_hline(yintercept = REF_LINE,
             linetype = "dashed", color = "grey55", linewidth = 0.6) +
  geom_text(aes(label = bar_label),
            vjust = -0.5, size = 3.5, color = "grey20", fontface = "plain") +
  scale_fill_manual(values = COND_COLORS, guide = "none") +
  scale_y_continuous(
    limits = c(0, max(df$ratio, na.rm = TRUE) * 1.18),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  annotate("text",
           x    = 0.55,
           y    = REF_LINE + 0.012,
           label = "25 %",
           size  = 3, color = "grey45", hjust = 0) +
  labs(
    title    = "Proportion of immune-enriched bins among immune-positive bins",
    subtitle = paste0("Enriched = bins with n_immune > mean(n_immune | n_immune > 0) ",
                      "within sample\n",
                      "Labels: n_enrich / n_pos  |  100 µm grid"),
    x = NULL,
    y = "Ratio (enriched / positive bins)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 8, color = "grey40"),
    axis.text.x      = element_text(size = 11),
    axis.text.y      = element_text(size = 10),
    axis.line        = element_line(linewidth = 0.5, color = "black"),
    axis.ticks       = element_line(linewidth = 0.4, color = "black"),
    plot.margin      = margin(10, 15, 10, 10)
  )

# =============================================================
# 4. Sauvegarde
# =============================================================

stem <- file.path(OUT_DIR, "fig_grid_immune_ratio_barplot")

cairo_pdf(paste0(stem, ".pdf"), width = 6, height = 5)
print(p)
dev.off()

CairoJPEG(paste0(stem, ".jpg"), width = 6 * 150, height = 5 * 150, res = 150)
print(p)
dev.off()

message("  Saved: fig_grid_immune_ratio_barplot.pdf / .jpg")
message("\n=== Done ===")
