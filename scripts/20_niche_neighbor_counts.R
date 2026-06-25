#!/usr/bin/env Rscript
# =============================================================
# Script: 20_niche_neighbor_counts.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-05
#
# Goal:
#   Quantify the number of neighbors each cell has within two
#   candidate spatial radii (15 µm and 30 µm), broken down by
#   neighbor annotation category.
#   This compares a strict local scale (15 µm) vs a broader
#   local neighborhood scale (30 µm).
#
# Input:
#   objects/08_immune_annotated_lam02_res03.rds
#     spatialCoords(se) → sdimx / sdimy (microns)
#     colData(se)$cell_type  (17 annotations)
#     colData(se)$sample     (4 conditions)
#
# Method:
#   dbscan::frNN(coords, eps = r) returns, for each cell, the
#   indices of ALL cells within radius r (self excluded).
#   Neighbor counts are then tallied by annotation category.
#
# Radii:
#   15 µm = strict / very local neighborhood
#   30 µm = broader local neighborhood
#
# Outputs (folder: outputs/banksy/inflammatory_niche_step2_neighbor_counts/):
#   neighbor_counts_per_cell.csv          — per-cell × radius table
#   neighbor_counts_summary.csv           — summary by annotation × sample × radius
#   neighbor_counts_immune_focus.csv      — immune-focus subset
#   fig_same_anno_neighbors_global.pdf/jpg
#   fig_same_anno_neighbors_immune.pdf/jpg
#   fig_total_neighbors_comparison.pdf/jpg
# =============================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(tidyverse)
  library(ggplot2)
})

if (!requireNamespace("dbscan", quietly = TRUE)) {
  stop("Package 'dbscan' is required. Install with: install.packages('dbscan')")
}

base_path <- normalizePath(".")  # Run this script from the project root directory
setwd(base_path)

out_dir <- file.path("outputs", "banksy",
                     "inflammatory_niche_step2_neighbor_counts")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE,
                                      showWarnings = FALSE)

# =============================================================
# Parameters
# =============================================================

SAMPLE_ORDER <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi", "mock_6wpi")

sample_labels <- c(
  mock_6wpi = "Mock 6 wpi",
  LCMV_1wpi = "LCMV 1 wpi",
  LCMV_3wpi = "LCMV 3 wpi",
  LCMV_6wpi = "LCMV 6 wpi"
)

sample_colors <- c(
  mock_6wpi = "#999999",
  LCMV_1wpi = "#56B4E9",
  LCMV_3wpi = "#E69F00",
  LCMV_6wpi = "#D55E00"
)

RADII <- c(15L, 30L)   # µm

# Flexible label patterns for neighbor categories
PAT_IMMUNE    <- c("T cells", "Mono", "Mac", "Microglia", "IFN", "Ifit")
PAT_MICROGLIA <- c("Microglia")
PAT_T_CELLS   <- c("T cells")
PAT_MONO      <- c("Mono")
PAT_MAC       <- c("Mac")

PREFERRED_ORDER <- c("T cells", "Mono", "Mac", "Microglia")

# =============================================================
# 1. Load object
# =============================================================

message("Loading: objects/08_immune_annotated_lam02_res03.rds")
se <- readRDS(file.path("objects", "08_immune_annotated_lam02_res03.rds"))
message("  ", ncol(se), " cells | ",
        length(unique(colData(se)$cell_type)), " annotations")

cd  <- as.data.frame(colData(se))
xy  <- as.data.frame(SpatialExperiment::spatialCoords(se))

if (!all(c("sdimx", "sdimy") %in% colnames(xy))) {
  stop("Expected 'sdimx' and 'sdimy' in spatialCoords(se)")
}

cell_df <- data.frame(
  cell_id   = colnames(se),
  sample    = factor(as.character(cd$sample), levels = SAMPLE_ORDER),
  cell_type = as.character(cd$cell_type),
  sdimx     = xy$sdimx,
  sdimy     = xy$sdimy,
  stringsAsFactors = FALSE
)

all_ct <- sort(unique(cell_df$cell_type))

cat("\nCells per sample:\n")
print(table(cell_df$sample))
cat("\nCells per annotation:\n")
print(sort(table(cell_df$cell_type), decreasing = TRUE))

# =============================================================
# 2. Resolve annotation categories from flexible patterns
# =============================================================

resolve_labels <- function(patterns, labels) {
  unique(unlist(lapply(patterns, function(p) {
    grep(p, labels, value = TRUE, ignore.case = TRUE)
  })))
}

immune_cts    <- resolve_labels(PAT_IMMUNE,    all_ct)
microglia_cts <- resolve_labels(PAT_MICROGLIA, all_ct)
t_cells_cts   <- resolve_labels(PAT_T_CELLS,   all_ct)
mono_cts      <- resolve_labels(PAT_MONO,      all_ct)
mac_cts       <- resolve_labels(PAT_MAC,       all_ct)

cat("\nAnnotation categories resolved:\n")
cat("  Immune     :", paste(immune_cts,    collapse = ", "), "\n")
cat("  Microglia  :", paste(microglia_cts, collapse = ", "), "\n")
cat("  T cells    :", paste(t_cells_cts,   collapse = ", "), "\n")
cat("  Mono       :", paste(mono_cts,      collapse = ", "), "\n")
cat("  Mac        :", paste(mac_cts,       collapse = ", "), "\n")

# =============================================================
# 3. Compute neighbor counts per cell × radius
#
#   For each sample separately:
#     dbscan::frNN(coords, eps = r)  → for each cell, indices of
#     all cells within r µm (self excluded by frNN).
#   Tally by category.
# =============================================================

message("\nComputing neighbor counts at radii: ",
        paste(RADII, collapse = ", "), " µm...")
message("  Package: dbscan::frNN | Self excluded by default")

count_rows <- vector("list", length(SAMPLE_ORDER) * length(RADII))
ptr <- 1L

for (samp in levels(cell_df$sample)) {

  idx_s  <- which(cell_df$sample == samp)
  ct_s   <- cell_df$cell_type[idx_s]
  cid_s  <- cell_df$cell_id[idx_s]
  xy_s   <- as.matrix(cell_df[idx_s, c("sdimx", "sdimy")])
  n_s    <- length(idx_s)

  message("  Sample ", samp, " (n = ", n_s, ")")

  for (r in RADII) {

    nn <- dbscan::frNN(xy_s, eps = r, sort = FALSE)
    # nn$id[[i]] = integer vector of indices (1-based, within this sample)
    # that are within radius r of cell i; excludes self

    n_total_r   <- vapply(nn$id, length, integer(1))
    n_same_r    <- vapply(seq_len(n_s), function(i) {
                     sum(ct_s[nn$id[[i]]] == ct_s[i])
                   }, integer(1))
    n_immune_r  <- vapply(nn$id, function(idx) sum(ct_s[idx] %in% immune_cts),    integer(1))
    n_micro_r   <- vapply(nn$id, function(idx) sum(ct_s[idx] %in% microglia_cts), integer(1))
    n_t_cells_r <- vapply(nn$id, function(idx) sum(ct_s[idx] %in% t_cells_cts),   integer(1))
    n_mono_r    <- vapply(nn$id, function(idx) sum(ct_s[idx] %in% mono_cts),       integer(1))
    n_mac_r     <- vapply(nn$id, function(idx) sum(ct_s[idx] %in% mac_cts),        integer(1))

    count_rows[[ptr]] <- data.frame(
      cell_id          = cid_s,
      sample           = samp,
      cell_type        = ct_s,
      radius_um        = r,
      n_total          = n_total_r,
      n_same_anno      = n_same_r,
      n_immune         = n_immune_r,
      n_microglia      = n_micro_r,
      n_t_cells        = n_t_cells_r,
      n_mono           = n_mono_r,
      n_mac            = n_mac_r,
      stringsAsFactors = FALSE
    )
    ptr <- ptr + 1L
  }
}

counts_df <- bind_rows(count_rows) %>%
  mutate(
    sample    = factor(sample,    levels = SAMPLE_ORDER),
    radius_um = factor(radius_um, levels = RADII,
                       labels = paste0(RADII, " µm"))
  )

cat("\n  Total rows:", nrow(counts_df),
    "(", nrow(counts_df) / length(RADII), "cells ×", length(RADII), "radii)\n")

write.csv(counts_df,
          file.path(out_dir, "neighbor_counts_per_cell.csv"),
          row.names = FALSE)
message("Saved: neighbor_counts_per_cell.csv")

# =============================================================
# 4. Summary tables
# =============================================================

message("\nBuilding summary tables...")

# Helper: robust summary of a count column
count_summary <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) {
    return(c(mean = NA_real_, median = NA_real_,
             q75 = NA_real_, q90 = NA_real_))
  }
  c(
    mean   = round(mean(x),                   2),
    median = round(median(x),                 2),
    q75    = round(quantile(x, 0.75), 2),
    q90    = round(quantile(x, 0.90), 2)
  )
}

nn_summary <- counts_df %>%
  group_by(cell_type, sample, radius_um) %>%
  summarise(
    n_cells            = n(),
    mean_total         = round(mean(n_total),     2),
    median_total       = round(median(n_total),   2),
    q75_total          = round(quantile(n_total, 0.75), 2),
    mean_same_anno     = round(mean(n_same_anno), 2),
    median_same_anno   = round(median(n_same_anno), 2),
    q75_same_anno      = round(quantile(n_same_anno, 0.75), 2),
    mean_immune        = round(mean(n_immune),    2),
    median_immune      = round(median(n_immune),  2),
    q75_immune         = round(quantile(n_immune, 0.75), 2),
    mean_microglia     = round(mean(n_microglia), 2),
    median_microglia   = round(median(n_microglia), 2),
    mean_t_cells       = round(mean(n_t_cells),   2),
    median_t_cells     = round(median(n_t_cells), 2),
    mean_mono          = round(mean(n_mono),       2),
    median_mono        = round(median(n_mono),     2),
    mean_mac           = round(mean(n_mac),         2),
    median_mac         = round(median(n_mac),       2),
    .groups = "drop"
  )

write.csv(nn_summary,
          file.path(out_dir, "neighbor_counts_summary.csv"),
          row.names = FALSE)
message("Saved: neighbor_counts_summary.csv")

# Immune-focus subset
immune_focus_cts <- resolve_labels(
  c("T cells", "Mono", "Mac", "Microglia", "IFN", "Ifit"),
  all_ct
)

nn_immune_summary <- nn_summary %>%
  filter(cell_type %in% immune_focus_cts)

write.csv(nn_immune_summary,
          file.path(out_dir, "neighbor_counts_immune_focus.csv"),
          row.names = FALSE)
message("Saved: neighbor_counts_immune_focus.csv")

cat("\nImmune-focus median same-annotation neighbors (15 µm vs 30 µm, pooled conditions):\n")
print(
  nn_immune_summary %>%
    group_by(cell_type, radius_um) %>%
    summarise(
      med_same = round(median(median_same_anno, na.rm = TRUE), 1),
      med_imm  = round(median(median_immune,    na.rm = TRUE), 1),
      .groups  = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = radius_um,
                       values_from = c(med_same, med_imm)),
  n = 20
)

# =============================================================
# 5. Publication theme
# =============================================================

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1,
                                      hjust = 0, margin = margin(b = 4)),
      plot.subtitle    = element_text(size = base_size - 2, color = "grey40",
                                      hjust = 0, lineheight = 1.3,
                                      margin = margin(b = 8)),
      axis.text        = element_text(size = base_size - 1.5),
      axis.title       = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size - 1, face = "bold"),
      legend.text      = element_text(size = base_size - 1.5),
      legend.key.size  = unit(0.4, "cm"),
      panel.grid.major = element_line(color = "grey93", linewidth = 0.3),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      plot.margin      = margin(12, 16, 12, 12)
    )
}

radius_colors <- c("15 µm" = "#2166ac", "30 µm" = "#d6604d")

# Biological cell type ordering
preferred_matched <- unique(unlist(lapply(PREFERRED_ORDER, function(pat) {
  grep(pat, all_ct, value = TRUE, ignore.case = TRUE)
})))
remaining    <- sort(setdiff(all_ct, preferred_matched))
ct_order_bio <- c(preferred_matched, remaining)

# =============================================================
# 6. Plot 1 — Same-annotation neighbor counts, all cell types
#
#   Violin + median point, both radii dodged, conditions pooled.
#   Y-axis: cell type. X-axis: n_same_anno.
#   Separate panels for 15 µm and 30 µm to keep it readable.
# =============================================================

message("\nPlot 1: same-annotation neighbor counts, all cell types...")

plot1_df <- counts_df %>%
  mutate(cell_type = factor(cell_type, levels = rev(ct_order_bio)))

med1_df <- plot1_df %>%
  group_by(cell_type, radius_um) %>%
  summarise(med_val = median(n_same_anno), .groups = "drop")

p1 <- ggplot(plot1_df,
             aes(x = n_same_anno, y = cell_type, fill = radius_um)) +
  geom_violin(color = NA, alpha = 0.75, scale = "width", trim = TRUE,
              position = position_dodge(width = 0.85)) +
  geom_point(data = med1_df,
             aes(x = med_val, y = cell_type, group = radius_um),
             shape = 21, fill = "white", color = "grey20",
             size = 2.0, stroke = 0.7,
             position = position_dodge(width = 0.85)) +
  scale_fill_manual(values = radius_colors, name = "Radius") +
  scale_x_continuous(name = "Same-annotation neighbors within radius",
                     expand = expansion(mult = c(0, 0.04))) +
  labs(
    title    = "Same-annotation neighbor counts — all cell types",
    subtitle = paste0(
      "All conditions pooled | White circle = median\n",
      "Object: 08_immune_annotated_lam02_res03.rds | ",
      ncol(se), " cells | µm (MERFISH)"
    ),
    y = NULL
  ) +
  theme_pub() +
  theme(
    legend.position    = "top",
    axis.text.y        = element_text(size = 8.5),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.4)
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))

n_ct  <- length(ct_order_bio)
h_p1  <- max(6, n_ct * 0.45 + 2.5)

ggsave(file.path(out_dir, "fig_same_anno_neighbors_global.pdf"),
       plot = p1, width = 11, height = h_p1, device = cairo_pdf)
ggsave(file.path(out_dir, "fig_same_anno_neighbors_global.jpg"),
       plot = p1, width = 11, height = h_p1, dpi = 300)
message("  Saved: fig_same_anno_neighbors_global  (11 x ", round(h_p1, 1), " in)")

# =============================================================
# 7. Plot 2 — Immune focus: same-annotation neighbors
#             faceted by condition, 15 vs 30 µm in colour
# =============================================================

message("\nPlot 2: immune focus, faceted by sample...")

immune_pref_order <- ct_order_bio[ct_order_bio %in% immune_focus_cts]

plot2_df <- counts_df %>%
  filter(cell_type %in% immune_focus_cts) %>%
  mutate(
    cell_type = factor(cell_type, levels = rev(immune_pref_order)),
    sample    = factor(sample, levels = SAMPLE_ORDER)
  )

if (nrow(plot2_df) == 0) {
  message("  No immune-focus cells found — skipping Plot 2.")
} else {
  med2_df <- plot2_df %>%
    group_by(cell_type, sample, radius_um) %>%
    summarise(med_val = median(n_same_anno), .groups = "drop")

  p2 <- ggplot(plot2_df,
               aes(x = n_same_anno, y = cell_type, fill = radius_um)) +
    geom_violin(color = NA, alpha = 0.78, scale = "width", trim = TRUE,
                position = position_dodge(width = 0.85)) +
    geom_point(data = med2_df,
               aes(x = med_val, y = cell_type, group = radius_um),
               shape = 21, fill = "white", color = "grey20",
               size = 2.2, stroke = 0.75,
               position = position_dodge(width = 0.85)) +
    scale_fill_manual(values = radius_colors, name = "Radius") +
    scale_x_continuous(name = "Same-annotation neighbors within radius",
                       expand = expansion(mult = c(0, 0.04))) +
    facet_wrap(~ sample, nrow = 1, labeller = as_labeller(sample_labels)) +
    labs(
      title    = "Same-annotation neighbor counts — immune focus populations",
      subtitle = paste0(
        "T cells, Mono, Mac, Microglia (+ IFN if present)\n",
        "White circle = median | µm (MERFISH)"
      ),
      y = NULL
    ) +
    theme_pub() +
    theme(
      legend.position    = "top",
      axis.text.y        = element_text(size = 9),
      axis.text.x        = element_text(size = 8),
      panel.grid.major.x = element_line(color = "grey88", linewidth = 0.4)
    ) +
    guides(fill = guide_legend(override.aes = list(alpha = 1)))

  n_imm <- length(immune_pref_order)
  h_p2  <- max(5, n_imm * 0.6 + 3)

  ggsave(file.path(out_dir, "fig_same_anno_neighbors_immune.pdf"),
         plot = p2, width = 16, height = h_p2, device = cairo_pdf)
  ggsave(file.path(out_dir, "fig_same_anno_neighbors_immune.jpg"),
         plot = p2, width = 16, height = h_p2, dpi = 300)
  message("  Saved: fig_same_anno_neighbors_immune  (16 x ",
          round(h_p2, 1), " in)")
}

# =============================================================
# 8. Plot 3 — Total neighbor counts (all annotations) 15 vs 30
#
#   Point-range (median ± IQR) per annotation × radius,
#   conditions pooled. Clean compact overview.
# =============================================================

message("\nPlot 3: total neighbor count point-range, all annotations...")

ptrange_df <- counts_df %>%
  group_by(cell_type, radius_um) %>%
  summarise(
    med   = median(n_total),
    q25   = quantile(n_total, 0.25),
    q75   = quantile(n_total, 0.75),
    .groups = "drop"
  ) %>%
  mutate(cell_type = factor(cell_type, levels = rev(ct_order_bio)))

p3 <- ggplot(ptrange_df,
             aes(x = med, y = cell_type, color = radius_um,
                 xmin = q25, xmax = q75)) +
  geom_linerange(linewidth = 1.0, alpha = 0.7,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 2.8, shape = 21, fill = "white", stroke = 1.0,
             position = position_dodge(width = 0.6)) +
  scale_color_manual(values = radius_colors, name = "Radius") +
  scale_x_continuous(name = "Total neighbors within radius (median ± IQR)",
                     expand = expansion(mult = c(0, 0.04)),
                     limits = c(0, NA)) +
  labs(
    title    = "Total neighbor counts per cell — all annotations",
    subtitle = paste0(
      "Point = median | Bar = IQR (Q25–Q75) | All conditions pooled\n",
      "All cells within radius r regardless of annotation"
    ),
    y = NULL
  ) +
  theme_pub() +
  theme(
    legend.position    = "top",
    axis.text.y        = element_text(size = 8.5),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.4)
  )

ggsave(file.path(out_dir, "fig_total_neighbors_comparison.pdf"),
       plot = p3, width = 10, height = h_p1, device = cairo_pdf)
ggsave(file.path(out_dir, "fig_total_neighbors_comparison.jpg"),
       plot = p3, width = 10, height = h_p1, dpi = 300)
message("  Saved: fig_total_neighbors_comparison  (10 x ", round(h_p1, 1), " in)")

# =============================================================
# 9. Console summary: 15 vs 30 µm for key immune types
# =============================================================

cat("\n========================================================\n")
cat("Quick comparison: median total vs same-annotation neighbors\n")
cat("========================================================\n")
print(
  nn_summary %>%
    filter(cell_type %in% immune_focus_cts) %>%
    select(cell_type, sample, radius_um, n_cells,
           median_total, median_same_anno,
           median_immune, median_microglia) %>%
    arrange(cell_type, radius_um, sample),
  n = 60
)

message("\n=== Done. Outputs in: ", out_dir, " ===\n")
