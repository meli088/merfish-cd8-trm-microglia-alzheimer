#!/usr/bin/env Rscript
# =============================================================
# Script: 18_inflammatory_niche_step1_distance_same_annotation.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-05
#
# Goal:
#   Characterize spatial proximity structure within each annotated
#   cell type, as a first step toward defining inflammatory niche
#   neighborhoods.
#
# Input:
#   objects/08_immune_annotated_lam02_res03.rds
#   SpatialExperiment with:
#     - spatialCoords(se) → sdimx / sdimy (microns)
#     - colData(se)$cell_type  (17 annotations)
#     - colData(se)$sample     (4 conditions)
#     - colData(se)$Anti.RFP_raw if present (TRM-contact signal)
#
# Steps:
#   1. Plot Anti.RFP_raw expression (spatial + UMAP if available)
#   2. Compute within-annotation nearest-neighbor distance per cell
#      (same sample, exclude self)
#   3. Summarize and plot distance distributions (violin + density)
#   4. Propose candidate neighborhood thresholds from quantiles
#
# Output folder:
#   outputs/banksy/inflammatory_niche_step1_distance_same_annotation/
#
# Outputs:
#   fig_rfp_spatial_<sample>.pdf/jpg   — Anti.RFP_raw spatial map per sample
#   within_annotation_nn_distances.csv  — per-cell table
#   within_annotation_nn_summary.csv    — per annotation × sample summary
#   fig_nn_violin_by_annotation.pdf/jpg — violin plots
#   fig_nn_violin_by_sample.pdf/jpg     — violin plots faceted by sample
#   candidate_thresholds.csv            — quantile-based threshold table
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

out_dir <- file.path("outputs", "banksy",
                     "inflammatory_niche_step1_distance_same_annotation")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE,
                                      showWarnings = FALSE)

# =============================================================
# Parameters
# =============================================================

SAMPLE_ORDER <- c("mock_6wpi", "LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")

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

# Candidate threshold quantiles (used in step 4)
THRESHOLD_QUANTILES <- c(0.10, 0.25, 0.50)

# Annotation color palette — immune populations (bright) + others (muted)
color_map <- c(
  # Immune populations — distinct bright colors
  "T cells (Gzmb)"                      = "#E41A1C",
  "Microglia (C1qa)"                  = "#377EB8",
  "Mac (Ctss)"                        = "#984EA3",
  "Mono (Lyz2)"                       = "#4DAF4A",
  "T cell / Neuron doublet / Cycling 1" = "#F781BF",
  # Glial / neural populations — muted
  "Excitatory neurons (Satb2)"        = "#A6CEE3",
  "Inhibitory neurons (Htr2c)"        = "#B2DF8A",
  "Neurons (Fam107a)"                 = "#CAB2D6",
  "Neurons (Rbfox3)"                  = "#FDBF6F",
  "Glials (Gja1)"                     = "#B15928",
  "Oligo"                             = "#999999",
  "Oligo ? (Nkx2-9)"                  = "#BBBBBB",
  "OPC (Sox10) / Cycling 2"           = "#6A3D9A",
  "Vascular"                          = "#1F78B4",
  "Vascular (Igfbp2)"                 = "#33A02C",
  "Unknown"                           = "#D9D9D9"
)

# =============================================================
# 1. Load object — inspect structure
# =============================================================

message("Loading: objects/08_immune_annotated_lam02_res03.rds")
se <- readRDS(file.path("objects", "08_immune_annotated_lam02_res03.rds"))
message("  ", ncol(se), " cells loaded")

cd <- as.data.frame(colData(se))

cat("\nColData columns:\n")
print(colnames(cd))
cat("\nAssay names:\n")
print(assayNames(se))

# Extract spatial coordinates
xy <- as.data.frame(SpatialExperiment::spatialCoords(se))
if (!all(c("sdimx", "sdimy") %in% colnames(xy))) {
  stop("Spatial coordinates must have columns 'sdimx' and 'sdimy'")
}
cat("\nSpatial coordinate ranges:\n")
cat("  sdimx:", round(min(xy$sdimx), 1), "to", round(max(xy$sdimx), 1), "µm\n")
cat("  sdimy:", round(min(xy$sdimy), 1), "to", round(max(xy$sdimy), 1), "µm\n")

# Check for Anti.RFP_raw in colData
rfp_in_cd <- "Anti.RFP_raw" %in% colnames(cd)
cat("\nAnti.RFP_raw in colData:", rfp_in_cd, "\n")

# Build master cell table (coordinates + annotations)
cell_df <- data.frame(
  cell_id   = colnames(se),
  sample    = factor(as.character(cd$sample), levels = SAMPLE_ORDER),
  cell_type = as.character(cd$cell_type),
  sdimx     = xy$sdimx,
  sdimy     = xy$sdimy,
  stringsAsFactors = FALSE
)

if (rfp_in_cd) {
  cell_df$Anti.RFP_raw <- as.numeric(cd$Anti.RFP_raw)
}

cat("\nCells per sample:\n")
print(table(cell_df$sample))
cat("\nCells per annotation:\n")
print(sort(table(cell_df$cell_type), decreasing = TRUE))

# If Anti.RFP_raw not in colData(se), load it from 03_all_clustered.rds
# (Seurat @meta.data — Anti.RFP_raw was stored there in 00_load_data.R
#  but not transferred to colData when converting to SpatialExperiment)
if (!rfp_in_cd) {
  rfp_obj_path <- file.path("objects", "03_all_clustered.rds")
  if (file.exists(rfp_obj_path)) {
    message("\nAnti.RFP_raw not in colData — loading from 03_all_clustered.rds ...")
    seurat_rfp <- readRDS(rfp_obj_path)
    if ("Anti.RFP_raw" %in% colnames(seurat_rfp@meta.data)) {
      rfp_vals <- seurat_rfp@meta.data[colnames(se), "Anti.RFP_raw"]
      rm(seurat_rfp); gc()
      n_matched <- sum(!is.na(rfp_vals))
      cat("  Cells matched:", n_matched, "/", ncol(se), "\n")
      if (n_matched > 0) {
        cell_df$Anti.RFP_raw <- as.numeric(rfp_vals)
        rfp_in_cd <- TRUE
        cat("  Anti.RFP_raw loaded from 03_all_clustered.rds\n")
      }
    } else {
      rm(seurat_rfp); gc()
      message("  Anti.RFP_raw not found in 03_all_clustered.rds — skipping RFP plots.")
    }
  } else {
    message("  03_all_clustered.rds not found — skipping RFP plots.")
  }
}

# =============================================================
# 2. Publication theme
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

# =============================================================
# 3. Anti.RFP_raw spatial plots (per sample)
# =============================================================

if (rfp_in_cd) {

  message("\nPlotting Anti.RFP_raw spatial maps...")

  rfp_range <- range(cell_df$Anti.RFP_raw, na.rm = TRUE)
  rfp_pos   <- cell_df[!is.na(cell_df$Anti.RFP_raw) &
                          cell_df$Anti.RFP_raw > 0, ]
  cat("  Anti.RFP_raw range:", round(rfp_range[1], 3), "–",
      round(rfp_range[2], 3), "\n")
  cat("  Cells with RFP > 0:", nrow(rfp_pos), "/", nrow(cell_df), "\n")

  for (samp in levels(cell_df$sample)) {

    df_s <- cell_df %>% filter(sample == samp)
    n_rfp_pos <- sum(!is.na(df_s$Anti.RFP_raw) & df_s$Anti.RFP_raw > 0)

    p_rfp <- ggplot(df_s, aes(x = sdimx, y = sdimy,
                               color = Anti.RFP_raw)) +
      geom_point(size = 0.25, alpha = 0.7, shape = 16) +
      scale_color_gradientn(
        name   = "Anti.RFP_raw\n(protein signal)",
        colors = c("grey90", "#ffffb2", "#fd8d3c", "#b10026"),
        na.value = "grey85"
      ) +
      coord_equal() +
      labs(
        title    = paste0("Anti-RFP spatial signal — ", sample_labels[samp]),
        subtitle = paste0(
          "Object: 08_immune_annotated_lam02_res03.rds | Assay: colData$Anti.RFP_raw\n",
          "RFP+ cells (> 0): ", n_rfp_pos, " / ", nrow(df_s),
          " (TRM-contact cells, CRE-lox labeling)"
        ),
        x = "X (µm)", y = "Y (µm)"
      ) +
      theme_pub() +
      theme(
        axis.text     = element_text(size = 7),
        legend.position = "right"
      )

    fname <- file.path(out_dir, paste0("fig_rfp_spatial_", samp))
    ggsave(paste0(fname, ".pdf"), plot = p_rfp,
           width = 7, height = 6, device = cairo_pdf)
    ggsave(paste0(fname, ".jpg"), plot = p_rfp,
           width = 7, height = 6, dpi = 300)
    message("  Saved: fig_rfp_spatial_", samp)
  }

  # All samples faceted
  p_rfp_all <- ggplot(cell_df, aes(x = sdimx, y = sdimy, color = Anti.RFP_raw)) +
    geom_point(size = 0.15, alpha = 0.65, shape = 16) +
    scale_color_gradientn(
      name   = "Anti.RFP_raw",
      colors = c("grey90", "#ffffb2", "#fd8d3c", "#b10026"),
      na.value = "grey85"
    ) +
    facet_wrap(~ sample, nrow = 2, labeller = as_labeller(sample_labels)) +
    coord_equal() +
    labs(
      title    = "Anti-RFP signal across all conditions",
      subtitle = "Object: colData$Anti.RFP_raw | TRM-contact CRE-lox signal",
      x = "X (µm)", y = "Y (µm)"
    ) +
    theme_pub() +
    theme(
      axis.text       = element_text(size = 6.5),
      strip.text      = element_text(face = "bold", size = 9),
      strip.background = element_rect(fill = "grey95", color = NA)
    )

  ggsave(file.path(out_dir, "fig_rfp_spatial_all_samples.pdf"),
         plot = p_rfp_all, width = 14, height = 11, device = cairo_pdf)
  ggsave(file.path(out_dir, "fig_rfp_spatial_all_samples.jpg"),
         plot = p_rfp_all, width = 14, height = 11, dpi = 250)
  message("  Saved: fig_rfp_spatial_all_samples")

} else {
  message("\nAnti.RFP_raw not found in colData — will try 03_all_clustered.rds fallback.")
}

# =============================================================
# 4. Within-annotation nearest-neighbor distances
#
#   For each cell i in annotation A, sample S:
#     Find the nearest cell j in annotation A, sample S, where j != i.
#   
#   Method: FNN::get.knnx(data = coords_A, query = coords_A, k = 2)
#     → second neighbor (k=2) excludes self (nearest at distance 0)
#   Cells in annotations with only 1 cell in that sample → NA.
#
#   Units: same as sdimx/sdimy (microns for MERFISH).
# =============================================================

message("\nComputing within-annotation nearest-neighbor distances...")
message("  Method: FNN::get.knnx, k=2, second distance = nearest non-self")
message("  Unit: microns (MERFISH spatial coordinates)")

samples    <- levels(cell_df$sample)
all_ct     <- unique(cell_df$cell_type)
nn_results <- vector("list", length = nrow(cell_df))

row_ptr <- 1L
for (samp in samples) {
  idx_samp <- which(cell_df$sample == samp)
  for (ct in all_ct) {
    idx_ct <- idx_samp[cell_df$cell_type[idx_samp] == ct]
    n <- length(idx_ct)

    if (n == 0L) next

    if (n == 1L) {
      nn_results[[row_ptr]] <- data.frame(
        cell_id       = cell_df$cell_id[idx_ct],
        sample        = samp,
        cell_type     = ct,
        sdimx         = cell_df$sdimx[idx_ct],
        sdimy         = cell_df$sdimy[idx_ct],
        nn_dist_um    = NA_real_,
        n_same_anno   = 1L,
        stringsAsFactors = FALSE
      )
      row_ptr <- row_ptr + 1L
      next
    }

    # k=2 because the nearest neighbor of itself is at distance 0
    coords_ct <- as.matrix(cell_df[idx_ct, c("sdimx", "sdimy")])
    knn       <- FNN::get.knnx(data = coords_ct, query = coords_ct, k = 2L)
    nn_dist   <- knn$nn.dist[, 2L]  # second neighbor = nearest non-self

    nn_results[[row_ptr]] <- data.frame(
      cell_id       = cell_df$cell_id[idx_ct],
      sample        = samp,
      cell_type     = ct,
      sdimx         = cell_df$sdimx[idx_ct],
      sdimy         = cell_df$sdimy[idx_ct],
      nn_dist_um    = round(nn_dist, 4),
      n_same_anno   = n,
      stringsAsFactors = FALSE
    )
    row_ptr <- row_ptr + 1L
  }
}

nn_df <- bind_rows(nn_results[!sapply(nn_results, is.null)]) %>%
  mutate(sample = factor(sample, levels = SAMPLE_ORDER))

cat("\nPer-cell distance table — head:\n")
print(head(nn_df))
cat("\nNA distances (singleton annotations):",
    sum(is.na(nn_df$nn_dist_um)), "\n")

write.csv(nn_df,
          file.path(out_dir, "within_annotation_nn_distances.csv"),
          row.names = FALSE)
message("Saved: within_annotation_nn_distances.csv")

# =============================================================
# 5. Summary table per annotation × sample
# =============================================================

message("\nComputing summary table...")

nn_summary <- nn_df %>%
  group_by(cell_type, sample) %>%
  summarise(
    n_cells     = n(),
    n_valid     = sum(!is.na(nn_dist_um)),
    mean_nn_um  = round(mean(nn_dist_um, na.rm = TRUE), 2),
    median_nn_um = round(median(nn_dist_um, na.rm = TRUE), 2),
    sd_nn_um    = round(sd(nn_dist_um, na.rm = TRUE), 2),
    q10_nn_um   = round(quantile(nn_dist_um, 0.10, na.rm = TRUE), 2),
    q25_nn_um   = round(quantile(nn_dist_um, 0.25, na.rm = TRUE), 2),
    q75_nn_um   = round(quantile(nn_dist_um, 0.75, na.rm = TRUE), 2),
    max_nn_um   = round(ifelse(all(is.na(nn_dist_um)), NA_real_, max(nn_dist_um, na.rm = TRUE)), 2),
    .groups     = "drop"
  )

cat("\nSummary table (median within-annotation NN distance in µm):\n")
print(
  nn_summary %>%
    select(cell_type, sample, n_valid, median_nn_um) %>%
    tidyr::pivot_wider(names_from = sample,
                       values_from = c(n_valid, median_nn_um),
                       names_glue = "{sample}_{.value}") %>%
    arrange(cell_type),
  n = 30
)

write.csv(nn_summary,
          file.path(out_dir, "within_annotation_nn_summary.csv"),
          row.names = FALSE)
message("Saved: within_annotation_nn_summary.csv")

# =============================================================
# 6. Violin plots — within-annotation NN distance
# =============================================================

message("\nPlotting violin plots...")

# Order annotations by median NN distance (ascending)
ct_median_order <- nn_df %>%
  filter(!is.na(nn_dist_um)) %>%
  group_by(cell_type) %>%
  summarise(med = median(nn_dist_um), .groups = "drop") %>%
  arrange(med) %>%
  pull(cell_type)

# Fill color_map with fallback for any annotation not in palette
missing_cts <- setdiff(all_ct, names(color_map))
if (length(missing_cts) > 0) {
  extra <- setNames(rep("#CCCCCC", length(missing_cts)), missing_cts)
  color_map <- c(color_map, extra)
}

# Exclude singletons for plotting
nn_plot_df <- nn_df %>%
  filter(!is.na(nn_dist_um)) %>%
  mutate(
    cell_type = factor(cell_type, levels = ct_median_order),
    sample    = factor(sample, levels = SAMPLE_ORDER)
  )

# ── Fig A: all conditions pooled ──
p_violin_all <- ggplot(
  nn_plot_df,
  aes(x = cell_type, y = nn_dist_um, fill = cell_type)
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
    y = "Nearest same-annotation neighbor distance (µm)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(angle = 40, hjust = 1, vjust = 1, size = 10),
    axis.text.y  = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.line    = element_line(linewidth = 0.5, color = "black"),
    axis.ticks   = element_line(linewidth = 0.4, color = "black"),
    plot.margin  = margin(10, 12, 10, 10)
  )

ggsave(file.path(out_dir, "fig_nn_violin_by_annotation.pdf"),
       plot = p_violin_all, width = 10, height = 5.5, device = cairo_pdf)
ggsave(file.path(out_dir, "fig_nn_violin_by_annotation.jpg"),
       plot = p_violin_all, width = 10, height = 5.5, dpi = 300)
message("  Saved: fig_nn_violin_by_annotation")

# ── Fig B: by sample (dodged) ──
p_violin_samp <- ggplot(
  nn_plot_df,
  aes(x = cell_type, y = nn_dist_um, fill = sample)
) +
  geom_violin(
    color = "#333333",
    linewidth = 0.3,
    alpha = 0.85,
    scale = "width",
    trim = FALSE,
    position = position_dodge(width = 0.85)
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 95,
    size = 3.2,
    color = "#333333",
    position = position_dodge(width = 0.85)
  ) +
  scale_fill_manual(values = sample_colors, labels = sample_labels,
                    name = "Condition") +
  labs(
    x = NULL,
    y = "Nearest same-annotation neighbor distance (µm)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.x  = element_text(angle = 40, hjust = 1, vjust = 1, size = 10),
    axis.text.y  = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.line    = element_line(linewidth = 0.5, color = "black"),
    axis.ticks   = element_line(linewidth = 0.4, color = "black"),
    plot.margin  = margin(10, 12, 10, 10)
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))

ggsave(file.path(out_dir, "fig_nn_violin_by_sample.pdf"),
       plot = p_violin_samp, width = 13, height = 6, device = cairo_pdf)
ggsave(file.path(out_dir, "fig_nn_violin_by_sample.jpg"),
       plot = p_violin_samp, width = 13, height = 6, dpi = 300)
message("  Saved: fig_nn_violin_by_sample")

# ── Fig C: immune populations only ──
immune_pats <- c("T cells", "Mono", "Mac", "Microglia",
                 "T cell.*Neuron doublet")
immune_cts  <- unique(unlist(lapply(immune_pats, function(p) {
  grep(p, all_ct, value = TRUE, ignore.case = TRUE)
})))

# Sort immune types by ascending median NN distance
immune_order <- ct_median_order[ct_median_order %in% immune_cts]

nn_immune_df <- nn_plot_df %>%
  filter(cell_type %in% immune_cts) %>%
  mutate(cell_type = factor(as.character(cell_type), levels = immune_order))

if (nrow(nn_immune_df) > 0) {
  p_violin_immune <- ggplot(
    nn_immune_df,
    aes(x = cell_type, y = nn_dist_um, fill = cell_type)
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
      y = "Nearest same-annotation neighbor distance (µm)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x  = element_text(angle = 40, hjust = 1, vjust = 1, size = 10),
      axis.text.y  = element_text(size = 10),
      axis.title.y = element_text(size = 12),
      axis.line    = element_line(linewidth = 0.5, color = "black"),
      axis.ticks   = element_line(linewidth = 0.4, color = "black"),
      plot.margin  = margin(10, 12, 10, 10)
    )

  ggsave(file.path(out_dir, "fig_nn_violin_immune_only.pdf"),
         plot = p_violin_immune, width = 8, height = 5, device = cairo_pdf)
  ggsave(file.path(out_dir, "fig_nn_violin_immune_only.jpg"),
         plot = p_violin_immune, width = 8, height = 5, dpi = 300)
  message("  Saved: fig_nn_violin_immune_only")
}

# =============================================================
# 7. Candidate neighborhood thresholds
#
#   Strategy:
#   Compute quantiles of within-annotation NN distances for each
#   annotation, pooled across all conditions. These quantiles serve as
#   candidate radius thresholds r:
#     - r_q10: very tight neighborhood (captures densely packed cells)
#     - r_q25: compact neighborhood
#     - r_q50: median distance (half of same-type cells are within r)
#
#   A threshold r means: two cells of different types are considered
#   "neighbors" if their distance is <= r µm.
#
#   NOTE: Do NOT use these values as final biological decisions.
#   They are a quantitative basis for selecting r in downstream analysis.
# =============================================================

message("\nComputing candidate neighborhood thresholds...")      

thresh_df <- nn_df %>%
  filter(!is.na(nn_dist_um)) %>%
  group_by(cell_type) %>%
  summarise(
    n_cells    = n(),
    q10_um     = round(quantile(nn_dist_um, 0.10), 2),
    q25_um     = round(quantile(nn_dist_um, 0.25), 2),
    q50_um     = round(quantile(nn_dist_um, 0.50), 2),
    mean_um    = round(mean(nn_dist_um), 2),
    sd_um      = round(sd(nn_dist_um), 2),
    .groups    = "drop"
  ) %>%
  arrange(q25_um)

cat("\n=== Candidate neighborhood thresholds (µm) by annotation ===\n")
print(thresh_df, n = 30)

# Global quantiles across all annotations (to get single candidate values)
global_thresh <- nn_df %>%
  filter(!is.na(nn_dist_um)) %>%
  group_by(cell_type) %>%
  summarise(
    q10 = quantile(nn_dist_um, 0.10),
    q25 = quantile(nn_dist_um, 0.25),
    q50 = quantile(nn_dist_um, 0.50),
    .groups = "drop"
  ) %>%
  summarise(
    across(c(q10, q25, q50),
           list(
             min    = ~ round(min(.), 1),
             median = ~ round(median(.), 1),
             max    = ~ round(max(.), 1)
           ))
  )

cat("\nGlobal range of per-annotation quantiles (across all cell types):\n")
print(t(global_thresh))

# Specific thresholds worth inspecting (global medians of per-annotation quantiles)
q10_med <- round(median(thresh_df$q10_um), 1)
q25_med <- round(median(thresh_df$q25_um), 1)
q50_med <- round(median(thresh_df$q50_um), 1)

candidate_global <- data.frame(
  threshold_name = c(
    "median_of_q10_across_annotations",
    "median_of_q25_across_annotations",
    "median_of_q50_across_annotations"
  ),
  value_um = c(q10_med, q25_med, q50_med),
  interpretation = c(
    "Tight: within the 10th-percentile same-type NN distance for a typical annotation",
    "Compact: within the 25th-percentile same-type NN distance for a typical annotation",
    "Median: within the 50th-percentile same-type NN distance for a typical annotation"
  ),
  note = "These are candidate values only. Final threshold requires biological judgment.",
  stringsAsFactors = FALSE
)

cat("\n=== Global candidate thresholds ===\n")
print(candidate_global)

# Save both tables
write.csv(thresh_df,
          file.path(out_dir, "candidate_thresholds_by_annotation.csv"),
          row.names = FALSE)
write.csv(candidate_global,
          file.path(out_dir, "candidate_thresholds_global.csv"),
          row.names = FALSE)
message("Saved: candidate_thresholds_by_annotation.csv")
message("Saved: candidate_thresholds_global.csv")

# Text summary for quick reference
summary_txt <- c(
  "=== Candidate neighborhood thresholds — quick reference ===",
  paste0("Generated: ", Sys.time()),
  paste0("Object: 08_immune_annotated_lam02_res03.rds (", ncol(se),
         " cells, ", length(all_ct), " annotations)"),
  "",
  "Method:",
  "  For each cell, distance to the nearest OTHER cell of the SAME annotation",
  "  in the SAME sample, using Euclidean distance on spatial coords (µm).",
  "  Singleton cells (only one cell of that type in that sample) → NA.",
  "",
  "Candidate global thresholds (median across annotations):",
  paste0("  Tight   (median of q10): ", q10_med, " µm"),
  paste0("  Compact (median of q25): ", q25_med, " µm"),
  paste0("  Median  (median of q50): ", q50_med, " µm"),
  "",
  "See candidate_thresholds_by_annotation.csv for per-annotation values.",
  "The final threshold choice is a biological decision — not made here."
)
writeLines(summary_txt, file.path(out_dir, "threshold_summary.txt"))
message("Saved: threshold_summary.txt")

message("\n=== Done. Outputs in: ", out_dir, " ===\n")
