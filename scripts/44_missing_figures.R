#!/usr/bin/env Rscript
# ============================================================================
# Script 44 — Régénération des figures manquantes (depuis CSV existants)
#
# Génère 5 figures à partir des CSV déjà calculés — aucun .rds requis.
#
# Figures produites :
#   1. fig_volcano_in_vs_out_lcmv_1wpi.pdf/jpg
#   2. fig_volcano_in_vs_out_lcmv_3wpi.pdf/jpg
#   3. fig_volcano_in_vs_out_lcmv_6wpi.pdf/jpg
#      → outputs/banksy/microglia_dam_niche/
#
#   4. violin_nearest_immune_distance_lam02_res09.pdf/jpg
#      → outputs/banksy/nearest_immune_distance_OK/
#
#   5. fig_grid_ifn_faceted_shared_scale.pdf/jpg
#      → outputs/banksy/inflammatory_niche_step2_grid_100um_global/
#
# Sources CSV :
#   • microglia_dam_niche/DEG_microglia_in_vs_out_lcmv_[1/3/6]wpi.csv
#   • nearest_immune_distance_OK/nearest_immune_distance_per_cell_lam02_res09.csv
#   • inflammatory_niche_step2_grid_100um_global/grid_per_bin.csv
# ============================================================================

set.seed(1997)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(Cairo)
})

base_path <- normalizePath(".")
setwd(base_path)
source("scripts/00_palette.R")

# ------------------------------------------------------------------
# Paramètres communs volcanos
# ------------------------------------------------------------------
FDR_CUTOFF  <- 0.05
FC_CUTOFF   <- 0.25
TOP_N_LABEL <- 15

direction_colors <- c("up" = "#B2182B", "down" = "#2166AC", "ns" = "grey75")

TIMEPOINT_LABELS <- c(
  LCMV_1wpi = "LCMV 1 wpi",
  LCMV_3wpi = "LCMV 3 wpi",
  LCMV_6wpi = "LCMV 6 wpi"
)

out_niche <- file.path("outputs", "banksy", "microglia_dam_niche")

# Helper sauvegarde
save_fig <- function(p, dir, fname, width, height) {
  CairoPDF(file.path(dir, paste0(fname, ".pdf")), width = width, height = height)
  print(p); dev.off()
  CairoJPEG(file.path(dir, paste0(fname, ".jpg")),
            width = width * 150, height = height * 150, res = 150)
  print(p); dev.off()
  message("  Saved: ", fname)
}

# ==================================================================
# FIGURES 1–3 — Volcanos in vs out niche par timepoint
# ==================================================================
message("\n=== Volcanos in vs out niche par timepoint ===\n")

TIMEPOINTS <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")

for (tp in TIMEPOINTS) {
  csv_f <- file.path(out_niche,
                     paste0("DEG_microglia_in_vs_out_", tolower(tp), ".csv"))
  if (!file.exists(csv_f)) {
    warning("CSV manquant, timepoint ignoré : ", csv_f)
    next
  }
  message("  Timepoint : ", tp)

  deg <- read.csv(csv_f, stringsAsFactors = FALSE)

  # Recalcul direction si absent / cohérence
  deg <- deg %>%
    mutate(
      direction = case_when(
        avg_log2FC >  FC_CUTOFF & p_val_adj < FDR_CUTOFF ~ "up",
        avg_log2FC < -FC_CUTOFF & p_val_adj < FDR_CUTOFF ~ "down",
        TRUE ~ "ns"
      ),
      neg_log10_fdr = -log10(p_val_adj + 1e-300)
    )

  n_up   <- sum(deg$direction == "up",   na.rm = TRUE)
  n_down <- sum(deg$direction == "down", na.rm = TRUE)

  lab_genes <- bind_rows(
    deg %>% filter(direction == "up")   %>%
      slice_max(order_by = neg_log10_fdr, n = min(TOP_N_LABEL, n_up)),
    deg %>% filter(direction == "down") %>%
      slice_max(order_by = neg_log10_fdr, n = min(TOP_N_LABEL, n_down))
  )

  tp_label <- TIMEPOINT_LABELS[tp]

  p_vol <- ggplot(deg,
                  aes(x = avg_log2FC, y = neg_log10_fdr, color = direction)) +
    geom_point(size = 1.2, alpha = 0.7, stroke = 0) +
    geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF), linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    geom_text_repel(data = lab_genes,
                    aes(label = gene),
                    size = 2.8, max.overlaps = 25, min.segment.length = 0.2,
                    segment.color = "grey60", segment.size = 0.3) +
    scale_color_manual(
      values = direction_colors,
      labels = c(up   = paste0("Up in niche (", n_up, ")"),
                 down = paste0("Down in niche (", n_down, ")"),
                 ns   = "n.s."),
      name   = NULL
    ) +
    labs(
      title    = paste0("Microglia: In niche vs Out niche — ", tp_label),
      subtitle = sprintf("Wilcoxon | FDR < %.2f | log2FC > %.2f",
                         FDR_CUTOFF, FC_CUTOFF),
      x = "avg log2 FC (In / Out niche)",
      y = "-log10(FDR)"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 12, hjust = 0),
      plot.subtitle = element_text(size = 8, color = "grey40", hjust = 0),
      legend.position = "top"
    )

  slug <- tolower(tp)   # lcmv_1wpi / lcmv_3wpi / lcmv_6wpi
  save_fig(p_vol, out_niche,
           paste0("fig_volcano_in_vs_out_", slug), 7, 6)
}

# ==================================================================
# FIGURE 4 — Violin nearest immune distance
# ==================================================================
message("\n=== Violin nearest immune distance ===\n")

out_dist <- file.path("outputs", "banksy", "nearest_immune_distance_OK")
csv_dist <- file.path(out_dist, "nearest_immune_distance_per_cell_lam02_res09.csv")
stopifnot(file.exists(csv_dist))

dist_df <- read.csv(csv_dist, stringsAsFactors = FALSE) %>%
  mutate(annotation = trimws(annotation))

# Ordre par médiane croissante
anno_order <- dist_df %>%
  group_by(annotation) %>%
  summarise(med = median(nearest_immune_distance_um, na.rm = TRUE), .groups = "drop") %>%
  arrange(med) %>%
  pull(annotation)

plot_df <- dist_df %>%
  mutate(annotation = factor(annotation, levels = anno_order))

color_map <- GLOBAL_PALETTE[levels(plot_df$annotation)]
color_map[is.na(color_map)] <- "#CCCCCC"

p_violin <- ggplot(plot_df,
                   aes(x = annotation, y = nearest_immune_distance_um,
                       fill = annotation)) +
  geom_violin(scale = "width", trim = FALSE, width = 0.8,
              color = "#333333", linewidth = 0.35, alpha = 0.95) +
  stat_summary(fun = median, geom = "point", shape = 95,
               size = 3.2, color = "#333333") +
  scale_fill_manual(values = color_map, drop = FALSE) +
  labs(x = NULL, y = "Dist. from nearest Immune cell (µm)") +
  coord_cartesian(
    ylim = c(0, max(plot_df$nearest_immune_distance_um, na.rm = TRUE) * 1.02)
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 40, hjust = 1, vjust = 1, size = 10),
    axis.text.y     = element_text(size = 10),
    axis.title.y    = element_text(size = 12),
    axis.line       = element_line(linewidth = 0.5, color = "black"),
    axis.ticks      = element_line(linewidth = 0.4, color = "black"),
    plot.margin     = margin(10, 12, 10, 10)
  )

save_fig(p_violin, out_dist,
         "violin_nearest_immune_distance_lam02_res09", 8.2, 4.8)

# ==================================================================
# FIGURE 5 — IFN grid faceted shared scale
# ==================================================================
message("\n=== IFN grid faceted (shared scale) ===\n")

GRID_SIZE_UM <- 100
SAMPLE_ORDER  <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi", "mock_6wpi")
SAMPLE_LABELS <- c(mock_6wpi = "Mock 6 wpi", LCMV_1wpi = "LCMV 1 wpi",
                   LCMV_3wpi = "LCMV 3 wpi", LCMV_6wpi = "LCMV 6 wpi")

out_grid <- file.path("outputs", "banksy",
                      "inflammatory_niche_step2_grid_100um_global")
csv_grid <- file.path(out_grid, "grid_per_bin.csv")
stopifnot(file.exists(csv_grid))

bin_df <- read.csv(csv_grid, stringsAsFactors = FALSE) %>%
  mutate(
    sample       = factor(sample, levels = SAMPLE_ORDER),
    sample_label = SAMPLE_LABELS[as.character(sample)],
    sample_label = factor(sample_label, levels = SAMPLE_LABELS),
    # Coordonnées centrales de la tuile (en µm)
    tile_x       = grid_x * GRID_SIZE_UM + GRID_SIZE_UM / 2,
    tile_y       = grid_y * GRID_SIZE_UM + GRID_SIZE_UM / 2
  )

global_max_ifn <- max(bin_df$n_ifn, na.rm = TRUE)

theme_pub <- function() {
  theme_classic(base_size = 10) +
    theme(
      panel.border     = element_rect(color = "grey50", fill = NA, linewidth = 0.4),
      panel.background = element_rect(fill  = "grey96"),
      axis.line        = element_blank()
    )
}

p_ifn_facet <- ggplot(bin_df,
                      aes(x = tile_x, y = tile_y, fill = n_ifn)) +
  geom_tile(width  = GRID_SIZE_UM, height = GRID_SIZE_UM,
            color  = "grey85", linewidth = 0.07) +
  scale_fill_gradientn(
    name   = "IFN\n(Ifit1)\ncount",
    colors = c("grey96", "#edf8b1", "#7fcdbb", "#1d91c0", "#081d58"),
    limits = c(0, global_max_ifn),
    na.value = "grey96",
    guide  = guide_colorbar(barwidth = 0.7, barheight = 5)
  ) +
  facet_wrap(~ sample_label, nrow = 2) +
  coord_equal() +
  labs(
    title    = "IFN responsive (Ifit1) density — shared scale across samples",
    subtitle = paste0(GRID_SIZE_UM, " µm bins | shared scale: 0 – ",
                      global_max_ifn, " cells/bin"),
    x = "X (µm)", y = "Y (µm)"
  ) +
  theme_pub() +
  theme(
    axis.text        = element_text(size = 6.5),
    strip.text       = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.spacing    = unit(0.6, "cm")
  )

save_fig(p_ifn_facet, out_grid,
         "fig_grid_ifn_faceted_shared_scale", 11, 11)

# ==================================================================
# FIGURE 6 — Immune (Acod1) grid faceted shared scale — en ligne
# ==================================================================
message("\n=== Immune grid faceted (shared scale, en ligne) ===\n")

global_max_immune <- max(bin_df$n_immune, na.rm = TRUE)
message("  Global max n_immune: ", global_max_immune)

p_imm_facet <- ggplot(bin_df,
                      aes(x = tile_x, y = tile_y, fill = n_immune)) +
  geom_tile(width  = GRID_SIZE_UM, height = GRID_SIZE_UM,
            color  = "grey85", linewidth = 0.07) +
  scale_fill_gradientn(
    name   = "Immune\n(Acod1)\ncount",
    colors = c("grey96", "#fee8c8", "#fdbb84", "#e34a33", "#8b0000"),
    limits = c(0, global_max_immune),
    na.value = "white",
    guide  = guide_colorbar(barwidth = 0.7, barheight = 5)
  ) +
  facet_wrap(~ sample_label, nrow = 1) +
  coord_equal() +
  labs(
    title    = "Immune (Acod1) density — shared scale across samples",
    subtitle = paste0(GRID_SIZE_UM, " µm bins | shared scale: 0 – ",
                      global_max_immune, " cells/bin"),
    x = "X (µm)", y = "Y (µm)"
  ) +
  theme_pub() +
  theme(
    axis.text        = element_text(size = 6.5),
    strip.text       = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.spacing    = unit(0.6, "cm"),
    panel.background = element_rect(fill = "white")
  )

save_fig(p_imm_facet, out_grid,
         "fig_grid_immune_faceted_shared_scale", 24, 6)

message("\n=== Script 44 terminé avec succès ===")
message("Fichiers générés dans :")
message("  ", out_niche, "/fig_volcano_in_vs_out_lcmv_[1/3/6]wpi.pdf/jpg")
message("  ", out_dist,  "/violin_nearest_immune_distance_lam02_res09.pdf/jpg")
message("  ", out_grid,  "/fig_grid_ifn_faceted_shared_scale.pdf/jpg")
message("  ", out_grid,  "/fig_grid_immune_faceted_shared_scale.pdf/jpg")
