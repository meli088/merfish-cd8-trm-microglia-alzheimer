#!/usr/bin/env Rscript
# =============================================================
# Script: 19_niche_threshold_tables.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-05
#
# Goal:
#   Generate quantitative tables to support spatial neighborhood
#   threshold selection before making final niche figures.
#   No downstream enrichment or grid analysis here.
#
# Input:
#   objects/08_immune_annotated_lam02_res03.rds
#     spatialCoords(se)  → sdimx / sdimy (microns, MERFISH)
#     colData(se)$cell_type  (17 annotations)
#     colData(se)$sample     (4 conditions)
#
# Steps:
#   1. Per-cell nearest same-annotation neighbor distance (µm)
#   2. Summary table by annotation × sample (full stats)
#   3. Focused summary: key immune populations only
#   4. Threshold scan: % cells with >= 1 same-annotation neighbor
#      within each candidate radius (10, 15, 20, 25, 30, 40, 50 µm)
#   5. One compact helper plot
#
# Output folder:
#   outputs/banksy/inflammatory_niche_threshold_tables/
#
# Outputs:
#   nn_distances_per_cell.csv            — per-cell NN distances
#   nn_summary_by_annotation_sample.csv  — full stats, all annotations
#   nn_summary_immune_focus.csv          — key immune populations only
#   threshold_scan_long.csv              — % with >= 1 neighbor vs radius
#   fig_threshold_scan.pdf/jpg           — helper plot
# =============================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(FNN)
  library(tidyverse)
  library(ggplot2)
})

base_path <- normalizePath(".")  # Run this script from the project root directory
setwd(base_path)

out_dir <- file.path("outputs", "banksy", "inflammatory_niche_threshold_tables")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE,
                                      showWarnings = FALSE)

# =============================================================
# Parameters
# =============================================================

SAMPLE_ORDER <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi", "mock_6wpi")

CANDIDATE_RADII <- c(10, 15, 20, 25, 30, 40, 50)  # µm

# Focus populations for immune-only summary (flexible matching)
IMMUNE_PATTERNS <- c("T cells", "Mono", "Mac", "Microglia",
                     "IFN", "Ifit", "Interferon")

# =============================================================
# 1. Load object
# =============================================================

message("Loading: objects/08_immune_annotated_lam02_res03.rds")
se <- readRDS(file.path("objects", "08_immune_annotated_lam02_res03.rds"))
message("  ", ncol(se), " cells | ", length(unique(colData(se)$cell_type)),
        " annotations")

cd <- as.data.frame(colData(se))
xy <- as.data.frame(SpatialExperiment::spatialCoords(se))

if (!all(c("sdimx", "sdimy") %in% colnames(xy))) {
  stop("Expected columns 'sdimx' and 'sdimy' in spatialCoords(se)")
}

cat("\nSpatial coordinate ranges:\n")
cat("  sdimx:", round(range(xy$sdimx), 1), "µm\n")
cat("  sdimy:", round(range(xy$sdimy), 1), "µm\n")

cell_df <- data.frame(
  cell_id   = colnames(se),
  sample    = factor(as.character(cd$sample), levels = SAMPLE_ORDER),
  cell_type = as.character(cd$cell_type),
  sdimx     = xy$sdimx,
  sdimy     = xy$sdimy,
  stringsAsFactors = FALSE
)

all_ct <- sort(unique(cell_df$cell_type))

cat("\nCells per annotation:\n")
print(sort(table(cell_df$cell_type), decreasing = TRUE))

# =============================================================
# 2. Per-cell within-annotation nearest-neighbor distance
#
#   Method: FNN::get.knnx(k = 2)
#     k=1 returns self (distance = 0), k=2 returns the nearest OTHER cell.
#   Groups with n == 1 in that sample: nn_dist_um = NA.
#   Units: same as spatialCoords, i.e. microns.
# =============================================================

message("\nComputing within-annotation nearest-neighbor distances (µm)...")
message("  Method: FNN::get.knnx k=2, take 2nd neighbor to exclude self")

nn_rows <- vector("list", length = nrow(cell_df))
ptr <- 1L

for (samp in levels(cell_df$sample)) {
  idx_s <- which(cell_df$sample == samp)
  for (ct in all_ct) {
    idx_ct <- idx_s[cell_df$cell_type[idx_s] == ct]
    n <- length(idx_ct)
    if (n == 0L) next

    dist_val <- if (n == 1L) {
      NA_real_
    } else {
      coords <- as.matrix(cell_df[idx_ct, c("sdimx", "sdimy")])
      FNN::get.knnx(data = coords, query = coords, k = 2L)$nn.dist[, 2L]
    }

    nn_rows[[ptr]] <- data.frame(
      cell_id       = cell_df$cell_id[idx_ct],
      sample        = samp,
      cell_type     = ct,
      sdimx         = cell_df$sdimx[idx_ct],
      sdimy         = cell_df$sdimy[idx_ct],
      nn_dist_um    = if (n == 1L) NA_real_ else round(dist_val, 4),
      n_same_type_in_sample = n,
      stringsAsFactors = FALSE
    )
    ptr <- ptr + 1L
  }
}

nn_df <- bind_rows(nn_rows[!vapply(nn_rows, is.null, logical(1))]) %>%
  mutate(sample = factor(sample, levels = SAMPLE_ORDER))

cat("\n  Total cells:      ", nrow(nn_df), "\n")
cat("  NA (singletons):  ", sum(is.na(nn_df$nn_dist_um)), "\n")
cat("  Valid distances:  ", sum(!is.na(nn_df$nn_dist_um)), "\n")

write.csv(nn_df, file.path(out_dir, "nn_distances_per_cell.csv"),
          row.names = FALSE)
message("Saved: nn_distances_per_cell.csv")

# =============================================================
# 3. Summary table — all annotations × all samples
# =============================================================

message("\nBuilding summary table (all annotations × samples)...")

safe_quantile <- function(x, p) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_real_)
  unname(quantile(x, p))
}

nn_summary <- nn_df %>%
  group_by(cell_type, sample) %>%
  summarise(
    n_cells    = n(),
    n_valid    = sum(!is.na(nn_dist_um)),
    min_um     = if (n_valid > 0) round(min(nn_dist_um,  na.rm = TRUE), 2) else NA_real_,
    q10_um     = round(safe_quantile(nn_dist_um, 0.10), 2),
    q25_um     = round(safe_quantile(nn_dist_um, 0.25), 2),
    median_um  = round(safe_quantile(nn_dist_um, 0.50), 2),
    q75_um     = round(safe_quantile(nn_dist_um, 0.75), 2),
    q90_um     = round(safe_quantile(nn_dist_um, 0.90), 2),
    mean_um    = if (n_valid > 0) round(mean(nn_dist_um,  na.rm = TRUE), 2) else NA_real_,
    max_um     = if (n_valid > 0) round(max(nn_dist_um,   na.rm = TRUE), 2) else NA_real_,
    sd_um      = if (n_valid > 1) round(sd(nn_dist_um,    na.rm = TRUE), 2) else NA_real_,
    .groups = "drop"
  )

cat("\nSummary — median NN distance by annotation (pooled across conditions):\n")
print(
  nn_summary %>%
    group_by(cell_type) %>%
    summarise(
      median_pooled_um = round(median(median_um, na.rm = TRUE), 1),
      min_q10_um       = round(min(q10_um, na.rm = TRUE), 1),
      .groups = "drop"
    ) %>%
    arrange(median_pooled_um),
  n = 30
)

write.csv(nn_summary,
          file.path(out_dir, "nn_summary_by_annotation_sample.csv"),
          row.names = FALSE)
message("Saved: nn_summary_by_annotation_sample.csv")

# =============================================================
# 4. Focused summary — key immune populations only
# =============================================================

message("\nBuilding focused immune summary...")

immune_matched <- unique(unlist(lapply(IMMUNE_PATTERNS, function(pat) {
  grep(pat, all_ct, value = TRUE, ignore.case = TRUE)
})))

cat("\nImmune-focus cell types matched:\n")
print(immune_matched)

nn_immune_summary <- nn_summary %>%
  filter(cell_type %in% immune_matched)

if (nrow(nn_immune_summary) == 0) {
  warning("No immune-focus annotations found. Check IMMUNE_PATTERNS.")
} else {
  cat("\nFocused immune summary:\n")
  print(nn_immune_summary %>% select(cell_type, sample, n_cells, n_valid,
                                      q10_um, q25_um, median_um, q75_um),
        n = 40)
  write.csv(nn_immune_summary,
            file.path(out_dir, "nn_summary_immune_focus.csv"),
            row.names = FALSE)
  message("Saved: nn_summary_immune_focus.csv")
}

# =============================================================
# 5. Threshold scan
#
#   For each cell, and each candidate radius r:
#     has_neighbor_within_r = (nn_dist_um <= r)
#     (NA → FALSE: singleton has no neighbor at any radius)
#
#   Summarize by annotation × sample:
#     pct_with_neighbor = 100 * mean(has_neighbor)
#
#   Output: long-format CSV with columns:
#     cell_type, sample, radius_um, n_cells, n_valid,
#     n_with_neighbor, pct_with_neighbor
# =============================================================

message("\nRunning threshold scan for radii: ",
        paste(CANDIDATE_RADII, collapse = ", "), " µm...")

scan_rows <- vector("list", length(CANDIDATE_RADII))

for (i in seq_along(CANDIDATE_RADII)) {
  r <- CANDIDATE_RADII[i]
  scan_rows[[i]] <- nn_df %>%
    mutate(has_neighbor = !is.na(nn_dist_um) & nn_dist_um <= r) %>%
    group_by(cell_type, sample) %>%
    summarise(
      radius_um        = r,
      n_cells          = n(),
      n_valid          = sum(!is.na(nn_dist_um)),
      n_with_neighbor  = sum(has_neighbor),
      pct_with_neighbor = round(100 * mean(has_neighbor), 2),
      .groups = "drop"
    )
}

scan_df <- bind_rows(scan_rows) %>%
  mutate(
    sample    = factor(sample, levels = SAMPLE_ORDER),
    radius_um = as.integer(radius_um)
  ) %>%
  select(cell_type, sample, radius_um, n_cells, n_valid,
         n_with_neighbor, pct_with_neighbor) %>%
  arrange(cell_type, sample, radius_um)

cat("\nThreshold scan — sample (first rows):\n")
print(head(scan_df, 20))

write.csv(scan_df, file.path(out_dir, "threshold_scan_long.csv"),
          row.names = FALSE)
message("Saved: threshold_scan_long.csv")

# Quick console summary: median pct_with_neighbor across all annotations
cat("\nGlobal median % with >= 1 same-annotation neighbor by radius:\n")
print(
  scan_df %>%
    group_by(radius_um) %>%
    summarise(
      median_pct = round(median(pct_with_neighbor, na.rm = TRUE), 1),
      min_pct    = round(min(pct_with_neighbor,    na.rm = TRUE), 1),
      max_pct    = round(max(pct_with_neighbor,    na.rm = TRUE), 1),
      .groups = "drop"
    )
)

# =============================================================
# 6. Helper plot — threshold scan summary
#
#   Two panels:
#   Left:  median nearest-neighbor distance per annotation (all samples pooled)
#   Right: % cells with >= 1 same-annotation neighbor vs radius,
#          per annotation (lines), focus immune highlighted
# =============================================================

message("\nGenerating helper plot...")

# Panel A: median NN distance bar chart
med_pooled <- nn_summary %>%
  group_by(cell_type) %>%
  summarise(
    med_pooled = round(median(median_um, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  filter(!is.na(med_pooled)) %>%
  arrange(med_pooled) %>%
  mutate(
    cell_type = factor(cell_type, levels = cell_type),
    is_immune = cell_type %in% immune_matched
  )

p_med <- ggplot(med_pooled,
                aes(x = med_pooled, y = cell_type, fill = is_immune)) +
  geom_col(width = 0.7, alpha = 0.85) +
  scale_fill_manual(values = c("TRUE" = "#D55E00", "FALSE" = "steelblue"),
                    labels = c("TRUE" = "Immune focus", "FALSE" = "Other"),
                    name = NULL) +
  scale_x_continuous(
    name = "Median within-annotation\nNN distance (µm, pooled conditions)",
    expand = expansion(mult = c(0, 0.05))
  ) +
  geom_vline(xintercept = CANDIDATE_RADII, color = "grey60",
             linetype = "dashed", linewidth = 0.3) +
  labs(
    title = "A: Median nearest same-annotation neighbor distance",
    y = NULL
  ) +
  theme_classic(base_size = 9.5) +
  theme(
    plot.title      = element_text(face = "bold", size = 10, hjust = 0),
    axis.text.y     = element_text(size = 8),
    legend.position = "top",
    panel.grid.major.x = element_line(color = "grey93", linewidth = 0.3)
  )

# Panel B: threshold scan curves (pooled across conditions)
scan_pooled <- scan_df %>%
  group_by(cell_type, radius_um) %>%
  summarise(
    pct = round(weighted.mean(pct_with_neighbor,
                               w = n_cells, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(is_immune = cell_type %in% immune_matched)

p_scan <- ggplot(scan_pooled,
                 aes(x = radius_um, y = pct, group = cell_type,
                     color = is_immune, linewidth = is_immune,
                     alpha = is_immune)) +
  geom_line() +
  geom_point(size = 1.8) +
  scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "steelblue"),
                     labels = c("TRUE" = "Immune focus", "FALSE" = "Other"),
                     name = NULL) +
  scale_linewidth_manual(values = c("TRUE" = 0.9, "FALSE" = 0.45),
                          guide = "none") +
  scale_alpha_manual(values = c("TRUE" = 1.0, "FALSE" = 0.5),
                     guide = "none") +
  scale_x_continuous(name = "Candidate radius (µm)",
                     breaks = CANDIDATE_RADII) +
  scale_y_continuous(name = "% cells with >= 1 same-annotation neighbor",
                     limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(title = "B: Threshold scan — neighbor coverage vs. radius") +
  theme_classic(base_size = 9.5) +
  theme(
    plot.title      = element_text(face = "bold", size = 10, hjust = 0),
    legend.position = "top",
    panel.grid.major = element_line(color = "grey93", linewidth = 0.3)
  )

# Combine side by side
library(patchwork)
p_combined <- p_med + p_scan +
  plot_layout(widths = c(1, 1.3)) +
  plot_annotation(
    title    = "Spatial neighborhood threshold analysis",
    subtitle = paste0(
      "Object: 08_immune_annotated_lam02_res03.rds | ",
      ncol(se), " cells | Spatial coords in µm (MERFISH)\n",
      "Orange = immune focus populations | Dashed lines in A = candidate radii"
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 8.5, color = "grey40",
                                   lineheight = 1.3)
    )
  )

n_ct  <- nrow(med_pooled)
h_fig <- max(6, n_ct * 0.38 + 2.5)

ggsave(file.path(out_dir, "fig_threshold_scan.pdf"),
       plot = p_combined, width = 14, height = h_fig, device = cairo_pdf)
ggsave(file.path(out_dir, "fig_threshold_scan.jpg"),
       plot = p_combined, width = 14, height = h_fig, dpi = 300)
message("Saved: fig_threshold_scan  (14 x ", round(h_fig, 1), " in)")

message("\n=== Done. Outputs in: ", out_dir, " ===\n")
