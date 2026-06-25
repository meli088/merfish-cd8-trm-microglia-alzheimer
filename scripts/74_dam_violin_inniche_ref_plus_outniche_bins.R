#!/usr/bin/env Rscript
# ============================================================================
# Script 74 — DAM upregulated violin plot with IN-niche reference + OUT bins
#
# Goal:
#   For each timepoint (1wpi, 6wpi), plot DAM upregulated module score as violins
#   with x-axis order:
#   In-niche | <50µm | 50-100µm | 100-200µm | 200-300µm | >300µm
#
# Data sources (from script 39 outputs):
#   outputs/banksy/microglia_dam_niche/module_scores_per_cell.csv
#   outputs/banksy/microglia_dam_niche/distance_microglia_to_tcell_per_cell.csv
#
# Output:
#   outputs/banksy/immune_acod1/analysis/figures/
#     fig_dam_violin_inniche_ref_plus_outniche_bins.pdf
#     fig_dam_violin_inniche_ref_plus_outniche_bins.jpg
#
# ============================================================================

set.seed(1997)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(Cairo)
})

source("scripts/00_palette.R")

in_scores_file <- file.path("outputs", "banksy", "microglia_dam_niche", "module_scores_per_cell.csv")
in_dist_file <- file.path("outputs", "banksy", "microglia_dam_niche", "distance_microglia_to_tcell_per_cell.csv")
out_dir <- file.path("outputs", "banksy", "immune_acod1", "analysis", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_scores_file)) stop("Missing file: ", in_scores_file)
if (!file.exists(in_dist_file)) stop("Missing file: ", in_dist_file)

SAMPLE_KEEP <- c("LCMV_1wpi", "LCMV_6wpi")
SAMPLE_LABELS <- c("LCMV_1wpi" = "LCMV 1 wpi", "LCMV_6wpi" = "LCMV 6 wpi")

BIN_BREAKS <- c(0, 50, 100, 200, 300, Inf)
BIN_LABELS <- c("<50µm", "50-100µm", "100-200µm", "200-300µm", ">300µm")
X_LEVELS <- c("In-niche", BIN_LABELS)

message("Reading input tables...")
scores_df <- readr::read_csv(in_scores_file, show_col_types = FALSE) %>%
  mutate(cell_id_orig = sub("^(out|in)_", "", cell_id))

dist_df <- readr::read_csv(in_dist_file, show_col_types = FALSE)

# -------------------------------
# In-niche reference (no distance stratification)
# -------------------------------
in_ref <- scores_df %>%
  filter(niche_status == "In niche", sample %in% SAMPLE_KEEP) %>%
  transmute(
    sample,
    sample_label = SAMPLE_LABELS[sample],
    x_group = "In-niche",
    dam_score = Upregulated_DAM,
    fill_group = "In-niche"
  )

# -------------------------------
# Out-niche bins
# -------------------------------
out_scores <- scores_df %>%
  filter(niche_status == "Out niche", sample %in% SAMPLE_KEEP) %>%
  select(cell_id_orig, sample, dam_score = Upregulated_DAM)

out_dist <- dist_df %>%
  filter(niche_status == "Out niche", sample %in% SAMPLE_KEEP, !is.na(dist_um)) %>%
  select(cell_id, sample, dist_um)

out_join <- inner_join(out_scores, out_dist, by = c("cell_id_orig" = "cell_id", "sample")) %>%
  mutate(
    x_group = cut(
      dist_um,
      breaks = BIN_BREAKS,
      labels = BIN_LABELS,
      right = FALSE,
      include.lowest = TRUE
    ),
    sample_label = SAMPLE_LABELS[sample],
    fill_group = sample_label
  ) %>%
  filter(!is.na(x_group)) %>%
  select(sample, sample_label, x_group, dam_score, fill_group)

plot_df <- bind_rows(in_ref, out_join) %>%
  mutate(
    sample_label = factor(sample_label, levels = c("LCMV 1 wpi", "LCMV 6 wpi")),
    x_group = factor(x_group, levels = X_LEVELS),
    fill_group = factor(fill_group, levels = c("In-niche", "LCMV 1 wpi", "LCMV 6 wpi"))
  )

# Optional sample counts for reference
count_df <- plot_df %>%
  group_by(sample_label, x_group) %>%
  summarise(n = n(), .groups = "drop")

message("Counts by group:")
print(count_df)

fill_colors <- c(
  "In-niche" = "grey65",
  "LCMV 1 wpi" = "#56B4E9",
  "LCMV 6 wpi" = "#F28E2B"
)

p <- ggplot(plot_df, aes(x = x_group, y = dam_score, fill = fill_group)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.65, color = NA) +
  geom_boxplot(width = 0.16, outlier.size = 0.25, outlier.alpha = 0.25,
               color = "grey20", fill = "white") +
  scale_fill_manual(values = fill_colors, name = "Group") +
  facet_wrap(~sample_label, nrow = 1) +
  labs(
    title = "DAM upregulated score: In-niche reference + Out-of-niche distance bins",
    subtitle = "Activated microglia (C1qa) as reference, P2ry12 out-of-niche stratified by distance",
    x = "Distance to nearest T cell / Immune niche",
    y = "DAM upregulated module score (AddModuleScore)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    strip.text = element_text(face = "bold")
  )

pdf_out <- file.path(out_dir, "fig_dam_violin_inniche_ref_plus_outniche_bins.pdf")
jpg_out <- file.path(out_dir, "fig_dam_violin_inniche_ref_plus_outniche_bins.jpg")

CairoPDF(pdf_out, width = 13, height = 5.3)
print(p)
dev.off()

CairoJPEG(jpg_out, width = 13 * 170, height = 5.3 * 170, res = 170)
print(p)
dev.off()

message("Saved:")
message(pdf_out)
message(jpg_out)
