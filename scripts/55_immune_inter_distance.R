#!/usr/bin/env Rscript
# =============================================================
# Script: 55_immune_inter_distance.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-06
#
# Goal:
#   Pour chaque paire de types cellulaires immuns, calculer la distance
#   spatiale (µm) entre chaque cellule source et la cellule cible la
#   plus proche, au sein du même échantillon.
#   Produire :
#     - une heatmap des distances médianes (source × cible, par condition)
#     - des violin plots par paire (source → cible) facettés par condition
#
# Input:
#   objects/08_immune_annotated_lam02_res03.rds
#   SpatialExperiment avec :
#     - spatialCoords(se)  → sdimx / sdimy (µm)
#     - colData(se)$cell_type  (types immuns annotés)
#     - colData(se)$sample     (4 conditions)
#
# Outputs (outputs/banksy/immune_inter_distance/):
#   immune_pairwise_nn_distances.csv   — table par cellule × paire
#   immune_pairwise_nn_summary.csv     — médiane par paire × condition
#   fig_heatmap_median_distance_<sample>.pdf/.jpg
#   fig_violin_pairwise_<from>_<to>.pdf/.jpg
#   fig_heatmap_median_distance_all_samples.pdf/.jpg
# =============================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(FNN)
  library(tidyverse)
  library(ggplot2)
  library(Cairo)
})

base_path <- normalizePath(".")   # Run from project root
setwd(base_path)

source("scripts/00_palette.R")

# =============================================================
# Parameters
# =============================================================

SAMPLE_ORDER <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi", "mock_6wpi")

SAMPLE_LABELS <- c(
  mock_6wpi = "Mock 6 wpi",
  LCMV_1wpi = "LCMV 1 wpi",
  LCMV_3wpi = "LCMV 3 wpi",
  LCMV_6wpi = "LCMV 6 wpi"
)

SAMPLE_COLORS <- c(
  mock_6wpi = "#999999",
  LCMV_1wpi = "#56B4E9",
  LCMV_3wpi = "#E69F00",
  LCMV_6wpi = "#D55E00"
)

# Immune cell types to include (types présents dans 08_immune_annotated)
IMMUNE_TYPES <- c(
  "T cells (Gzmb)",
  "T CD4 (Foxp3)",
  "Mac (Ctss)",
  "Mono (Lyz2)",
  "Microglia (C1qa)",
  "B cells (Cd19)",
  "T cell / Neuron doublet / Cycling 1"
)

out_dir <- file.path("outputs", "banksy", "immune_inter_distance")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================
# Publication theme
# =============================================================

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1,
                                      hjust = 0, margin = margin(b = 4)),
      plot.subtitle    = element_text(size = base_size - 2, color = "grey40",
                                      hjust = 0, lineheight = 1.3,
                                      margin = margin(b = 10)),
      axis.text        = element_text(size = base_size - 1.5),
      axis.title       = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size - 1, face = "bold"),
      legend.text      = element_text(size = base_size - 1.5),
      legend.key.size  = unit(0.45, "cm"),
      panel.grid.major = element_line(color = "grey93", linewidth = 0.3),
      plot.margin      = margin(12, 16, 12, 12)
    )
}

make_slug <- function(x) gsub("[^a-z0-9]+", "_", tolower(trimws(x)))

save_plot <- function(p, path_no_ext, w, h) {
  ggsave(paste0(path_no_ext, ".pdf"), plot = p,
         width = w, height = h, device = cairo_pdf)
  ggsave(paste0(path_no_ext, ".jpg"), plot = p,
         width = w, height = h, dpi = 300)
}

# =============================================================
# 1. Load object
# =============================================================

message("Loading: objects/08_immune_annotated_lam02_res03.rds")
se <- readRDS(file.path("objects", "08_immune_annotated_lam02_res03.rds"))
message("  ", ncol(se), " cells loaded")

cd  <- as.data.frame(colData(se))
xy  <- as.data.frame(SpatialExperiment::spatialCoords(se))

stopifnot(all(c("sdimx", "sdimy") %in% colnames(xy)))
stopifnot("cell_type" %in% colnames(cd))
stopifnot("sample"    %in% colnames(cd))

cell_df <- data.frame(
  cell_id   = colnames(se),
  sample    = factor(as.character(cd$sample), levels = SAMPLE_ORDER),
  cell_type = as.character(cd$cell_type),
  sdimx     = xy$sdimx,
  sdimy     = xy$sdimy,
  stringsAsFactors = FALSE
)

# Keep only recognised immune types
cell_df <- cell_df %>% filter(cell_type %in% IMMUNE_TYPES)
cell_df$cell_type <- factor(cell_df$cell_type, levels = IMMUNE_TYPES)

cat("\nCells per sample (after immune filter):\n")
print(table(cell_df$sample))
cat("\nCells per cell_type:\n")
print(sort(table(cell_df$cell_type), decreasing = TRUE))

# =============================================================
# 2. Compute pairwise NN distances (per sample)
#    For each cell of type FROM, find nearest cell of type TO
#    in the same sample.  When FROM == TO, exclude self (k=2,
#    take 2nd neighbour).
# =============================================================

message("\nComputing pairwise nearest-neighbour distances...")

types_present <- levels(cell_df$cell_type)[
  levels(cell_df$cell_type) %in% unique(cell_df$cell_type)
]

pairs <- expand.grid(from = types_present, to = types_present,
                     stringsAsFactors = FALSE)

results_list <- vector("list", nrow(pairs))

for (i in seq_len(nrow(pairs))) {
  from_type <- pairs$from[i]
  to_type   <- pairs$to[i]
  same_type <- (from_type == to_type)

  rows_list <- vector("list", length(levels(cell_df$sample)))

  for (samp in levels(cell_df$sample)) {

    from_cells <- cell_df %>% filter(sample == samp, cell_type == from_type)
    to_cells   <- cell_df %>% filter(sample == samp, cell_type == to_type)

    if (nrow(from_cells) == 0) next

    # Need at least 2 targets when same type (to exclude self)
    min_to <- if (same_type) 2L else 1L
    if (nrow(to_cells) < min_to) next

    query  <- as.matrix(from_cells[, c("sdimx", "sdimy")])
    target <- as.matrix(to_cells[,   c("sdimx", "sdimy")])

    k_use <- min_to
    knn_res <- FNN::get.knnx(target, query, k = k_use)

    # When same type the 1st neighbour IS the cell itself → take col 2
    dist_col <- if (same_type) 2L else 1L
    dists    <- knn_res$nn.dist[, dist_col]

    rows_list[[match(samp, levels(cell_df$sample))]] <- data.frame(
      cell_id  = from_cells$cell_id,
      sample   = samp,
      from     = from_type,
      to       = to_type,
      dist_um  = dists,
      stringsAsFactors = FALSE
    )
  }

  results_list[[i]] <- bind_rows(rows_list)
}

dist_df <- bind_rows(results_list) %>%
  mutate(
    sample = factor(sample, levels = SAMPLE_ORDER),
    from   = factor(from,   levels = types_present),
    to     = factor(to,     levels = types_present)
  )

cat("\nRows in pairwise distance table:", nrow(dist_df), "\n")

# Save per-cell table
csv_cells <- file.path(out_dir, "immune_pairwise_nn_distances.csv")
write.csv(dist_df, csv_cells, row.names = FALSE)
message("Saved: ", csv_cells)

# =============================================================
# 3. Summary table — median per (from × to × sample)
# =============================================================

summary_df <- dist_df %>%
  group_by(from, to, sample) %>%
  summarise(
    n_cells       = n(),
    median_dist   = median(dist_um, na.rm = TRUE),
    mean_dist     = mean(dist_um,   na.rm = TRUE),
    sd_dist       = sd(dist_um,     na.rm = TRUE),
    q25           = quantile(dist_um, 0.25, na.rm = TRUE),
    q75           = quantile(dist_um, 0.75, na.rm = TRUE),
    .groups       = "drop"
  )

csv_summary <- file.path(out_dir, "immune_pairwise_nn_summary.csv")
write.csv(summary_df, csv_summary, row.names = FALSE)
message("Saved: ", csv_summary)

# =============================================================
# 4. Heatmap — median distance per sample
# =============================================================

message("\nPlotting heatmaps...")

for (samp in levels(dist_df$sample)) {

  hm_data <- summary_df %>%
    filter(sample == samp) %>%
    complete(from, to, fill = list(median_dist = NA, n_cells = 0)) %>%
    mutate(
      from = factor(from, levels = types_present),
      to   = factor(to,   levels = types_present)
    )

  if (nrow(hm_data) == 0) next

  max_dist <- max(hm_data$median_dist, na.rm = TRUE)

  p_hm <- ggplot(hm_data, aes(x = to, y = from, fill = median_dist)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(is.na(median_dist), "",
                                 paste0(round(median_dist), " µm"))),
              size = 2.8, color = "grey10") +
    scale_fill_gradientn(
      name   = "Median NN\ndist. (µm)",
      colors = c("#FFF7EC", "#FDD49E", "#FC8D59", "#D7301F", "#7F0000"),
      na.value = "grey92",
      limits = c(0, max_dist)
    ) +
    scale_x_discrete(position = "top") +
    labs(
      title    = paste0("Pairwise immune NN distances — ", SAMPLE_LABELS[samp]),
      subtitle = "Each cell: distance to nearest cell of target type (same sample, µm)",
      x        = "Target cell type (TO)",
      y        = "Source cell type (FROM)"
    ) +
    coord_fixed() +
    theme_pub(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 40, hjust = 0, size = 8),
      axis.text.y = element_text(size = 8),
      axis.line   = element_blank(),
      axis.ticks  = element_blank()
    )

  fname <- file.path(out_dir, paste0("fig_heatmap_median_distance_", samp))
  n_types <- length(types_present)
  save_plot(p_hm, fname, w = n_types * 1.3 + 2.5, h = n_types * 1.3 + 2.0)
  message("  Saved: fig_heatmap_median_distance_", samp)
}

# =============================================================
# 5. Heatmap — all samples faceted (small multiples)
# =============================================================

hm_all <- summary_df %>%
  complete(from, to, sample, fill = list(median_dist = NA)) %>%
  mutate(
    from   = factor(from,   levels = types_present),
    to     = factor(to,     levels = types_present),
    sample = factor(sample, levels = SAMPLE_ORDER),
    sample_label = SAMPLE_LABELS[as.character(sample)]
  )

max_dist_global <- max(hm_all$median_dist, na.rm = TRUE)

p_hm_all <- ggplot(hm_all, aes(x = to, y = from, fill = median_dist)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradientn(
    name   = "Median NN\ndist. (µm)",
    colors = c("#FFF7EC", "#FDD49E", "#FC8D59", "#D7301F", "#7F0000"),
    na.value = "grey92",
    limits = c(0, max_dist_global)
  ) +
  scale_x_discrete(position = "top") +
  facet_wrap(~ sample_label, ncol = 2) +
  labs(
    title    = "Pairwise immune NN distances — all conditions",
    subtitle = "Median distance (µm) from source to nearest target cell (same sample)",
    x        = "Target cell type (TO)",
    y        = "Source cell type (FROM)"
  ) +
  coord_fixed() +
  theme_pub(base_size = 9) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 0, size = 7),
    axis.text.y  = element_text(size = 7),
    axis.line    = element_blank(),
    axis.ticks   = element_blank(),
    strip.text   = element_text(face = "bold", size = 9),
    strip.background = element_blank()
  )

n_types <- length(types_present)
fname_all <- file.path(out_dir, "fig_heatmap_median_distance_all_samples")
save_plot(p_hm_all, fname_all,
          w = n_types * 1.1 * 2 + 4,
          h = n_types * 1.1 * 2 + 3)
message("Saved: fig_heatmap_median_distance_all_samples")

# =============================================================
# 6. Violin plots — distance distribution by target type,
#    one panel per source type, faceted by sample
# =============================================================

message("\nPlotting violin plots per source type...")

for (from_type in types_present) {

  df_from <- dist_df %>%
    filter(from == from_type) %>%
    mutate(
      sample_label = factor(SAMPLE_LABELS[as.character(sample)],
                            levels = SAMPLE_LABELS[SAMPLE_ORDER]),
      to_label     = as.character(to)
    )

  if (nrow(df_from) == 0) next

  # Order 'to' by overall median distance
  to_order <- df_from %>%
    group_by(to_label) %>%
    summarise(med = median(dist_um, na.rm = TRUE), .groups = "drop") %>%
    arrange(med) %>%
    pull(to_label)

  df_from$to_label <- factor(df_from$to_label, levels = to_order)

  # Colors from GLOBAL_PALETTE for target types
  to_colors <- GLOBAL_PALETTE[to_order]
  to_colors[is.na(to_colors)] <- "#CCCCCC"

  p_vio <- ggplot(df_from,
                  aes(x = to_label, y = dist_um, fill = to_label)) +
    geom_violin(scale = "width", trim = FALSE, width = 0.8,
                color = "#333333", linewidth = 0.3, alpha = 0.9) +
    stat_summary(fun = median, geom = "point",
                 shape = 95, size = 3.5, color = "#333333") +
    scale_fill_manual(values = to_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
    facet_wrap(~ sample_label, ncol = 2, scales = "free_y") +
    labs(
      title    = paste0("Nearest-neighbour distance — from ", from_type),
      subtitle = "Distance (µm) to nearest cell of each target immune type (same sample)",
      x        = "Target cell type",
      y        = "Distance to nearest cell (µm)"
    ) +
    theme_pub(base_size = 10) +
    theme(
      axis.text.x      = element_text(angle = 40, hjust = 1, size = 7.5),
      strip.text       = element_text(face = "bold", size = 9),
      strip.background = element_blank()
    )

  fname_vio <- file.path(out_dir,
                         paste0("fig_violin_from_", make_slug(from_type)))
  save_plot(p_vio, fname_vio,
            w = max(8, length(to_order) * 1.2 + 3),
            h = 7)
  message("  Saved: fig_violin_from_", make_slug(from_type))
}

# =============================================================
# 7. Summary print
# =============================================================

cat("\n=== Top 5 closest immune-cell-type pairs (pooled across samples) ===\n")
top5 <- summary_df %>%
  filter(as.character(from) != as.character(to)) %>%
  group_by(from, to) %>%
  summarise(overall_median = median(median_dist, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(overall_median) %>%
  slice_head(n = 5)
print(top5)

message("\nDone. Outputs in: ", out_dir)
