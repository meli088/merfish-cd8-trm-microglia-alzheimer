#!/usr/bin/env Rscript
# ============================================================================
# STEP: Nearest Immune Distance Analysis
# ============================================================================
# Purpose:
#   For each cell, compute the distance to the nearest Immune-annotated cell.
#   Generate violin plot (distance by annotation) and feature plot
#   (UMAP colored by distance gradient).
#
# Inputs:
#  - objects/04_banksy_joint_lam08_after_bloc3.rds (BANKSY object)
#  - ncells_by_sample_lam02_res09_joint_long.csv (annotation mapping)
#  - Composition table with spatial coordinates (reconstructed from object)
#
# Outputs:
#   - nearest_immune_distance_per_cell_lam02_res09.csv
#   - nearest_immune_distance_summary_by_annotation_lam02_res09.csv
#   - violin_nearest_immune_distance_lam02_res09.pdf/jpg
#   - featureplot_nearest_immune_distance_lam02_res09.pdf/jpg
# ============================================================================
library(SingleCellExperiment)
library(SummarizedExperiment)
library(SpatialExperiment)
library(tidyverse)
library(ggplot2)
library(FNN)

base_path <- normalizePath(".")  # Run this script from the project root directory
setwd(base_path)

source("scripts/00_palette.R")

obj_path <- file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds")
csv_annot <- "ncells_by_sample_lam02_res09_joint_long.csv"
out_dir <- file.path("outputs", "banksy", "nearest_immune_distance")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

find_cl_col <- function(se, lam, res) {
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  cn <- colnames(cd)
  lam_pat <- paste0("lam", gsub("\\.", "\\\\.", as.character(lam)))
  res_pat <- paste0("_res", gsub("\\.", "\\\\.", as.character(res)), "$")
  cols <- cn[grepl(lam_pat, cn) & grepl(res_pat, cn)]
  if (length(cols) == 0) stop("Clustering column not found")
  cols[1]
}

cat("\n=== Nearest Immune Distance ===\n")

se <- readRDS(obj_path)
cluster_col <- find_cl_col(se, 0.2, 0.9)

cd <- as.data.frame(SummarizedExperiment::colData(se))
banksy_clusters <- as.numeric(cd[[cluster_col]])
samples <- as.character(cd$sample)
cell_ids <- colnames(se)

spatial_mat <- as.data.frame(SpatialExperiment::spatialCoords(se))
if (!all(c("sdimx", "sdimy") %in% colnames(spatial_mat))) {
  stop("Spatial coordinates must have columns 'sdimx' and 'sdimy'")
}
x_coords <- spatial_mat$sdimx
y_coords <- spatial_mat$sdimy

anno_data <- read.delim(csv_annot, sep = ";", stringsAsFactors = FALSE)
anno_map <- anno_data %>%
  filter(annotation != "" & !is.na(annotation)) %>%
  select(banksy_domain, annotation) %>%
  distinct()

anno_map$cluster_id <- as.numeric(gsub("Domain_", "", anno_map$banksy_domain))
anno_lookup <- setNames(anno_map$annotation, anno_map$cluster_id)

domain_names <- paste0("Domain_", banksy_clusters)
cell_annotations <- anno_lookup[as.character(banksy_clusters)]
cell_annotations[is.na(cell_annotations)] <- "Unknown"

immune_idx <- grep("^Immune", cell_annotations, ignore.case = TRUE)
if (length(immune_idx) == 0) {
  stop("No Immune cells found. Check annotation mapping.")
}

all_coords <- cbind(x_coords, y_coords)
immune_coords <- cbind(x_coords[immune_idx], y_coords[immune_idx])

knn_result <- get.knnx(immune_coords, all_coords, k = 1)
nearest_immune_dist <- knn_result$nn.dist[, 1]

cell_results <- tibble(
  cell_id = cell_ids,
  sample = samples,
  banksy_domain = domain_names,
  annotation = cell_annotations,
  x_coordinate = x_coords,
  y_coordinate = y_coords,
  nearest_immune_distance_um = nearest_immune_dist
)

cell_csv <- file.path(out_dir, "nearest_immune_distance_per_cell_lam02_res09.csv")
write.csv(cell_results, cell_csv, row.names = FALSE)

summary_stats <- cell_results %>%
  group_by(annotation) %>%
  summarise(
    n_cells = n(),
    mean_distance_um = mean(nearest_immune_distance_um, na.rm = TRUE),
    median_distance_um = median(nearest_immune_distance_um, na.rm = TRUE),
    sd_distance_um = sd(nearest_immune_distance_um, na.rm = TRUE),
    min_distance_um = min(nearest_immune_distance_um, na.rm = TRUE),
    max_distance_um = max(nearest_immune_distance_um, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

summary_csv <- file.path(out_dir, "nearest_immune_distance_summary_by_annotation_lam02_res09.csv")
write.csv(summary_stats, summary_csv, row.names = FALSE)

# ------------------------------------------------------------------
# Violin plot styled closer to the example figure
# ------------------------------------------------------------------

anno_order <- summary_stats %>%
  arrange(median_distance_um) %>%
  pull(annotation)

plot_df <- cell_results %>%
  mutate(annotation = factor(annotation, levels = anno_order))

# Couleurs depuis la palette partagée (00_palette.R)
color_map <- GLOBAL_PALETTE[levels(plot_df$annotation)]
missing_annots <- names(color_map)[is.na(color_map)]
if (length(missing_annots) > 0) {
  color_map[missing_annots] <- "#CCCCCC"
}

p_violin <- ggplot(
  plot_df,
  aes(x = annotation, y = nearest_immune_distance_um, fill = annotation)
) +
  geom_violin(
    scale = "width",
    trim = FALSE,
    width = 0.8,
    color = "#333333",
    linewidth = 0.35,
    alpha = 0.95
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 95,
    size = 3.2,
    color = "#333333"
  ) +
  scale_fill_manual(values = color_map, drop = FALSE) +
  labs(
    x = NULL,
    y = "Dist. from nearest Immune cell (µm)"
  ) +
  coord_cartesian(ylim = c(0, max(plot_df$nearest_immune_distance_um, na.rm = TRUE) * 1.02)) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    plot.margin = margin(10, 12, 10, 10)
  )

pdf_violin <- file.path(out_dir, "violin_nearest_immune_distance_lam02_res09.pdf")
jpg_violin <- file.path(out_dir, "violin_nearest_immune_distance_lam02_res09.jpg")
ggsave(pdf_violin, p_violin, width = 8.2, height = 4.8)
ggsave(jpg_violin, p_violin, width = 8.2, height = 4.8, dpi = 300)

# ------------------------------------------------------------------
# Feature plot on UMAP
# ------------------------------------------------------------------
umap_name <- "UMAP_Harmony_lam0.2"
if (!(umap_name %in% reducedDimNames(se))) {
  stop("UMAP_Harmony_lam0.2 not found in reduced dimensions")
}

umap_mat <- reducedDim(se, umap_name)
umap_df <- tibble(
  UMAP_1 = umap_mat[, 1],
  UMAP_2 = umap_mat[, 2],
  nearest_immune_distance_um = nearest_immune_dist,
  annotation = cell_annotations,
  sample = samples
) %>%
  arrange(desc(nearest_immune_distance_um))   # draw far cells first, near cells last

# Optional: cap extreme distances for cleaner color scaling
dist_cap <- quantile(umap_df$nearest_immune_distance_um, 0.99, na.rm = TRUE)
umap_df <- umap_df %>%
  mutate(distance_plot = pmin(nearest_immune_distance_um, dist_cap))

non_immune_umap_df <- umap_df %>%
  filter(!grepl("^Immune", annotation, ignore.case = TRUE))
immune_umap_df <- umap_df %>%
  filter(grepl("^Immune", annotation, ignore.case = TRUE))

p_feature <- ggplot(
  non_immune_umap_df,
  aes(x = UMAP_1, y = UMAP_2, color = distance_plot)
) +
  geom_point(
    size = 0.12,
    alpha = 0.9,
    stroke = 0
  ) +
  geom_point(
    data = immune_umap_df,
    aes(x = UMAP_1, y = UMAP_2),
    color = "#00A651",
    size = 0.16,
    alpha = 0.95,
    stroke = 0,
    inherit.aes = FALSE
  ) +
  scale_color_viridis_c(
    option = "magma",
    begin = 0.03,
    end = 0.88,
    direction = 1,
    name = "Distance to nearest\nImmune cell (µm)",
    limits = c(0, dist_cap),
    oob = scales::squish
  ) +
  coord_equal() +
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(8, 12, 8, 8)
  )

pdf_feature <- file.path(out_dir, "featureplot_nearest_immune_distance_lam02_res09.pdf")
jpg_feature <- file.path(out_dir, "featureplot_nearest_immune_distance_lam02_res09.jpg")
ggsave(pdf_feature, p_feature, width = 6.6, height = 5.8)
ggsave(jpg_feature, p_feature, width = 6.6, height = 5.8, dpi = 400)

p_feature_split <- ggplot(
  non_immune_umap_df,
  aes(x = UMAP_1, y = UMAP_2, color = distance_plot)
) +
  geom_point(
    size = 0.10,
    alpha = 0.9,
    stroke = 0
  ) +
  geom_point(
    data = immune_umap_df,
    aes(x = UMAP_1, y = UMAP_2),
    color = "#00A651",
    size = 0.13,
    alpha = 0.95,
    stroke = 0,
    inherit.aes = FALSE
  ) +
  scale_color_viridis_c(
    option = "magma",
    begin = 0.03,
    end = 0.88,
    direction = 1,
    name = "Distance to nearest\nImmune cell (µm)",
    limits = c(0, dist_cap),
    oob = scales::squish
  ) +
  facet_wrap(~ sample, ncol = 2) +
  coord_equal() +
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    legend.position = "right",
    plot.margin = margin(8, 12, 8, 8)
  )

pdf_feature_split <- file.path(out_dir, "featureplot_nearest_immune_distance_lam02_res09_split.pdf")
jpg_feature_split <- file.path(out_dir, "featureplot_nearest_immune_distance_lam02_res09_split.jpg")
ggsave(pdf_feature_split, p_feature_split, width = 8.5, height = 7.2)
ggsave(jpg_feature_split, p_feature_split, width = 8.5, height = 7.2, dpi = 400)

cat("\nReference Immune cells:", length(immune_idx), "\n")
cat("Coordinates used: sdimx, sdimy\n")
cat("Distances are in microns\n")
cat("Outputs saved in:", out_dir, "\n")