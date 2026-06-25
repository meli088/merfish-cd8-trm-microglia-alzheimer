#!/usr/bin/env Rscript
# ============================================================================
# Script 56 - Figures 3.3 revised
#
# Section 1 - Bubble plot Ifng on global clusters
# Section 2 - XY IFN-responsive same scale across conditions
# Section 3 - Dotplot canonical markers revised v2
#
# Seed = 1997
# Uses cairo_pdf and source('scripts/00_palette.R')
# ============================================================================

set.seed(1997)

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(Banksy)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

base_path <- normalizePath(".")
setwd(base_path)
source("scripts/00_palette.R")

# ------------------------------------------------------------------
# Global paths and constants
# ------------------------------------------------------------------
PRIMARY_OBJ <- file.path("objects", "04_banksy_joint_lam02_res09.rds")
FALLBACK_OBJ <- file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds")
ANNOT_CSV <- "ncells_by_sample_lam02_res09_joint_long.csv"

OUT_IFN_DIR <- file.path("outputs", "banksy", "ifn_immune_overlay")
OUT_SPATIAL_DIR <- file.path("outputs", "banksy", "spatial_annotations")

if (!dir.exists(OUT_IFN_DIR)) dir.create(OUT_IFN_DIR, recursive = TRUE)
if (!dir.exists(OUT_SPATIAL_DIR)) dir.create(OUT_SPATIAL_DIR, recursive = TRUE)

SAMPLE_ORDER <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi", "mock_6wpi")
SAMPLE_LABELS <- c(
  LCMV_1wpi = "1wpi",
  LCMV_3wpi = "3wpi",
  LCMV_6wpi = "6wpi",
  mock_6wpi = "mock"
)

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
save_fig <- function(p, out_dir, fname, width, height, dpi = 300) {
  pdf_file <- file.path(out_dir, paste0(fname, ".pdf"))
  jpg_file <- file.path(out_dir, paste0(fname, ".jpg"))

  ggsave(pdf_file, p, width = width, height = height, device = cairo_pdf)
  ggsave(jpg_file, p, width = width, height = height, dpi = dpi)

  message("  Saved: ", pdf_file)
  message("  Saved: ", jpg_file)
}

find_cl_col <- function(se, lam, res) {
  all_cols <- clusterNames(se)
  lam_str <- gsub("\\.", "\\\\.", as.character(lam))
  lam_cols <- all_cols[grep(paste0("lam", lam_str), all_cols)]
  if (length(lam_cols) == 0) stop("Cluster column not found for lam=", lam)

  res_vals <- suppressWarnings(as.numeric(sub(".*_res", "", lam_cols)))
  idx <- which(!is.na(res_vals) & abs(res_vals - res) < 1e-8)
  if (length(idx) == 0) stop("Cluster column not found for res=", res)
  lam_cols[idx[1]]
}

pick_assay <- function(se) {
  an <- as.character(SummarizedExperiment::assayNames(se))
  if ("logcounts" %in% an) return("logcounts")
  if ("normcounts" %in% an) return("normcounts")
  if ("counts" %in% an) return("counts")
  an[1]
}

# ------------------------------------------------------------------
# Load object + reconstruct global annotations
# ------------------------------------------------------------------
message("\n=== Loading object and annotations ===")

obj_path <- if (file.exists(PRIMARY_OBJ)) PRIMARY_OBJ else FALLBACK_OBJ
if (!file.exists(obj_path)) {
  stop("Neither primary nor fallback object exists:\n - ", PRIMARY_OBJ, "\n - ", FALLBACK_OBJ)
}
if (obj_path != PRIMARY_OBJ) {
  message("WARNING: primary object missing, using fallback: ", obj_path)
}

if (!file.exists(ANNOT_CSV)) stop("Missing annotation CSV: ", ANNOT_CSV)

se <- readRDS(obj_path)
message("Object loaded: ", obj_path, " | cells=", ncol(se))

cl_col <- find_cl_col(se, lam = 0.2, res = 0.9)
cd <- as.data.frame(SummarizedExperiment::colData(se))

if (!"sample" %in% colnames(cd)) stop("sample column not found in colData")

anno_long <- read.delim(ANNOT_CSV, sep = ";", stringsAsFactors = FALSE) %>%
  dplyr::filter(annotation != "" & !is.na(annotation)) %>%
  dplyr::select(banksy_domain, annotation) %>%
  dplyr::distinct()

anno_map <- setNames(trimws(anno_long$annotation), anno_long$banksy_domain)
domain_labels <- paste0("Domain_", as.character(cd[[cl_col]]))
annotation <- anno_map[domain_labels]
annotation[is.na(annotation)] <- "Non annote"

sample_vec <- as.character(cd$sample)

# Coordinates
xy <- as.data.frame(SpatialExperiment::spatialCoords(se))
if (!all(c("sdimx", "sdimy") %in% colnames(xy))) {
  stop("Spatial coordinates must contain sdimx and sdimy")
}

# ------------------------------------------------------------------
# Section 1 - Ifng bubble plot global clusters
# ------------------------------------------------------------------
message("\n=== Section 1: Ifng bubble plot global clusters ===")

assay_use <- pick_assay(se)
expr_mat <- SummarizedExperiment::assay(se, assay_use)
if (!"Ifng" %in% rownames(expr_mat)) {
  stop("Ifng not found in selected assay rownames (assay=", assay_use, ")")
}

ifng_expr <- as.numeric(expr_mat["Ifng", ])

# nFeature_RNA from colData if present, otherwise compute from counts-like matrix.
if ("nFeature_RNA" %in% colnames(cd)) {
  nfeature <- as.numeric(cd$nFeature_RNA)
} else {
  cnt_assay <- if ("counts" %in% as.character(SummarizedExperiment::assayNames(se))) "counts" else assay_use
  cnt <- SummarizedExperiment::assay(se, cnt_assay)
  nfeature <- Matrix::colSums(cnt > 0)
}

nfeature[!is.finite(nfeature) | nfeature <= 0] <- NA_real_

norm_df <- data.frame(
  sample = sample_vec,
  ifng_expr = ifng_expr,
  nfeature = nfeature,
  annotation = annotation,
  stringsAsFactors = FALSE
)

# Normalize by nFeature_RNA and by sample-specific median depth.
med_nf <- norm_df %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(med_nfeature = median(nfeature, na.rm = TRUE), .groups = "drop")

norm_df <- norm_df %>%
  dplyr::left_join(med_nf, by = "sample") %>%
  dplyr::mutate(
    ifng_norm = ifelse(is.na(nfeature) | nfeature <= 0, NA_real_, ifng_expr / nfeature * med_nfeature),
    ifng_pos = ifng_expr > 0
  )

bubble_df <- norm_df %>%
  dplyr::filter(sample %in% SAMPLE_ORDER) %>%
  dplyr::group_by(annotation, sample) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    pct_expressing = 100 * mean(ifng_pos, na.rm = TRUE),
    median_ifng_pos = ifelse(sum(ifng_pos, na.rm = TRUE) > 0,
                             median(ifng_norm[ifng_pos], na.rm = TRUE),
                             NA_real_),
    .groups = "drop"
  )

ann_order <- bubble_df %>%
  dplyr::group_by(annotation) %>%
  dplyr::summarise(total_pct = sum(pct_expressing, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(total_pct)) %>%
  dplyr::pull(annotation)

bubble_df <- bubble_df %>%
  dplyr::mutate(
    annotation = factor(annotation, levels = rev(ann_order)),
    sample = factor(sample, levels = SAMPLE_ORDER, labels = SAMPLE_LABELS[SAMPLE_ORDER])
  )

bubble_pos <- bubble_df %>% dplyr::filter(pct_expressing > 0)
bubble_zero <- bubble_df %>% dplyr::filter(pct_expressing == 0)

p_bubble <- ggplot() +
  geom_point(data = bubble_zero,
             aes(x = sample, y = annotation),
             shape = 4, size = 1.8, color = "grey75", stroke = 0.7) +
  geom_point(data = bubble_pos,
             aes(x = sample, y = annotation,
                 size = pct_expressing,
                 color = median_ifng_pos),
             alpha = 0.9) +
  scale_size_continuous(name = "% Ifng+", range = c(1.8, 10), breaks = c(0.5, 1, 2, 5, 10)) +
  scale_color_gradientn(
    colors = c("#fee8c8", "#fdbb84", "#e34a33", "#b30000"),
    na.value = "grey80",
    name = "Median Ifng\n(normalized)\namong Ifng+"
  ) +
  labs(
    title = "Ifng bubble plot on global clusters",
    subtitle = "Size = % Ifng+ cells | Color = median Ifng (Ifng+), normalized by sample and nFeature_RNA",
    x = "Condition",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8.5, color = "grey40"),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 8),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

h_bub <- max(5, 0.34 * length(unique(bubble_df$annotation)) + 2)
save_fig(
  p_bubble,
  OUT_IFN_DIR,
  "fig_ifng_bubble_global_clusters",
  width = 9.8,
  height = h_bub
)

# Combined view: merge the global-all-samples and by-condition information
LCMV_ORDER <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")
PLOT_GROUPS <- c("Merged", "1wpi", "3wpi", "6wpi")

bubble_cond <- norm_df %>%
  dplyr::filter(sample %in% LCMV_ORDER) %>%
  dplyr::group_by(annotation, sample) %>%
  dplyr::summarise(
    pct_expressing = 100 * mean(ifng_pos, na.rm = TRUE),
    median_ifng_pos = ifelse(sum(ifng_pos, na.rm = TRUE) > 0,
                             median(ifng_norm[ifng_pos], na.rm = TRUE),
                             NA_real_),
    .groups = "drop"
  ) %>%
  dplyr::mutate(plot_group = dplyr::recode(sample,
                                           "LCMV_1wpi" = "1wpi",
                                           "LCMV_3wpi" = "3wpi",
                                           "LCMV_6wpi" = "6wpi")) %>%
  dplyr::select(annotation, plot_group, pct_expressing, median_ifng_pos)

bubble_merged <- norm_df %>%
  dplyr::filter(sample %in% LCMV_ORDER) %>%
  dplyr::group_by(annotation) %>%
  dplyr::summarise(
    pct_expressing = 100 * mean(ifng_pos, na.rm = TRUE),
    median_ifng_pos = ifelse(sum(ifng_pos, na.rm = TRUE) > 0,
                             median(ifng_norm[ifng_pos], na.rm = TRUE),
                             NA_real_),
    .groups = "drop"
  ) %>%
  dplyr::mutate(plot_group = "Merged") %>%
  dplyr::select(annotation, plot_group, pct_expressing, median_ifng_pos)

bubble_mix <- dplyr::bind_rows(bubble_merged, bubble_cond) %>%
  dplyr::mutate(
    annotation = factor(annotation, levels = rev(ann_order)),
    plot_group = factor(plot_group, levels = PLOT_GROUPS)
  )

bubble_mix_pos <- bubble_mix %>% dplyr::filter(pct_expressing > 0)
bubble_mix_zero <- bubble_mix %>% dplyr::filter(pct_expressing == 0)

p_bubble_mix <- ggplot() +
  geom_point(data = bubble_mix_zero,
             aes(x = plot_group, y = annotation),
             shape = 4, size = 1.8, color = "grey75", stroke = 0.7) +
  geom_point(data = bubble_mix_pos,
             aes(x = plot_group, y = annotation,
                 size = pct_expressing,
                 color = median_ifng_pos),
             alpha = 0.9) +
  scale_size_continuous(name = "% Ifng+", range = c(1.8, 10), breaks = c(0.5, 1, 2, 5, 10)) +
  scale_color_gradientn(
    colors = c("grey95", "#fee8c8", "#fdbb84", "#e34a33", "#b30000"),
    na.value = "grey80",
    name = "Median Ifng\n(normalized)\namong Ifng+"
  ) +
  labs(
    title = "Ifng bubble plot (Merged + conditions)",
    subtitle = "Global clusters | Size = % Ifng+ | Color = median Ifng in Ifng+ cells",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8.5, color = "grey40"),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 8),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

save_fig(
  p_bubble_mix,
  OUT_IFN_DIR,
  "fig_ifng_bubble_global_clusters_merged_plus_conditions",
  width = 10.6,
  height = h_bub
)

# ------------------------------------------------------------------
# Section 2 - XY IFN-responsive same scale with counts
# ------------------------------------------------------------------
message("\n=== Section 2: XY IFN same scale all conditions ===")

xy_df <- data.frame(
  x = xy$sdimx,
  y = xy$sdimy,
  sample = sample_vec,
  annotation = annotation,
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(sample %in% SAMPLE_ORDER)

xy_df$layer <- dplyr::case_when(
  xy_df$annotation == "IFN responsive (Ifit1)" ~ "IFN responsive (Ifit1)",
  TRUE ~ "Other"
)

xy_df$sample <- factor(xy_df$sample, levels = SAMPLE_ORDER, labels = SAMPLE_LABELS[SAMPLE_ORDER])

count_ifn <- xy_df %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(n_ifn = sum(layer == "IFN responsive (Ifit1)"), .groups = "drop")

xy_df <- xy_df %>%
  dplyr::left_join(count_ifn, by = "sample") %>%
  dplyr::mutate(panel_lab = paste0(as.character(sample), "\n", "n = ", n_ifn, " IFN-responsive cells"))

n_ifn_6wpi <- count_ifn %>%
  dplyr::filter(sample == SAMPLE_LABELS["LCMV_6wpi"]) %>%
  dplyr::pull(n_ifn)
if (length(n_ifn_6wpi) == 0) n_ifn_6wpi <- NA_integer_

global_xlim <- range(xy_df$x, na.rm = TRUE)
global_ylim <- range(xy_df$y, na.rm = TRUE)
xpad <- max(diff(global_xlim) * 0.05, 50)
ypad <- max(diff(global_ylim) * 0.05, 50)

p_xy <- ggplot(xy_df, aes(x = x, y = y)) +
  geom_point(data = xy_df %>% dplyr::filter(layer == "Other"),
             color = "#ECECEC", size = 0.25, alpha = 0.35, stroke = 0) +
  geom_point(data = xy_df %>% dplyr::filter(layer == "IFN responsive (Ifit1)"),
             color = "#F28E2B", size = 0.6, alpha = 0.9, stroke = 0) +
  facet_wrap(~ panel_lab, ncol = 2, scales = "fixed") +
  coord_fixed(
    xlim = c(global_xlim[1] - xpad, global_xlim[2] + xpad),
    ylim = c(global_ylim[1] - ypad, global_ylim[2] + ypad),
    expand = FALSE
  ) +
  labs(
    title = "XY IFN-responsive cells (same scale across conditions)",
    subtitle = paste0("Orange = IFN responsive (Ifit1), grey = all other cells | 6wpi appears sparse (n = ", n_ifn_6wpi, ")"),
    x = "X (um)",
    y = "Y (um)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 8.5, color = "grey40", hjust = 0.5),
    strip.text = element_text(face = "bold", size = 9),
    axis.text = element_text(size = 7),
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.35)
  )

save_fig(
  p_xy,
  OUT_IFN_DIR,
  "fig_xy_ifn_same_scale_all_conditions",
  width = 10.8,
  height = 8.2,
  dpi = 450
)

# ------------------------------------------------------------------
# Section 3 - Dotplot canonical markers revised v3
# ------------------------------------------------------------------
message("\n=== Section 3: Dotplot canonical markers revised v3 ===")

markers_v3 <- c(
  "Ptprc", "P2ry12", "Cldn5", "Rbfox3", "Sox10",
  "Aqp4", "Tmem119", "Plp1"
)

markers_use <- markers_v3[markers_v3 %in% rownames(expr_mat)]
missing_markers <- setdiff(markers_v3, markers_use)
if (length(markers_use) == 0) {
  stop("None of requested markers found in object")
}
if (length(missing_markers) > 0) {
  message("Missing markers in object: ", paste(missing_markers, collapse = ", "))
}

dot_df <- do.call(rbind, lapply(markers_use, function(g) {
  vals <- as.numeric(expr_mat[g, ])
  tmp <- data.frame(
    sample = sample_vec,
    annotation = annotation,
    expr = vals,
    marker = g,
    stringsAsFactors = FALSE
  )
  tmp %>%
    dplyr::filter(sample %in% SAMPLE_ORDER) %>%
    dplyr::group_by(sample, annotation, marker) %>%
    dplyr::summarise(
      pct_expr = 100 * mean(expr > 0, na.rm = TRUE),
      mean_expr = mean(expr, na.rm = TRUE),
      .groups = "drop"
    )
}))

ann_levels <- order_annotations(unique(dot_df$annotation), extended = TRUE)
ann_levels <- ann_levels[ann_levels %in% unique(dot_df$annotation)]

dot_df <- dot_df %>%
  dplyr::mutate(
    sample = factor(sample, levels = SAMPLE_ORDER, labels = SAMPLE_LABELS[SAMPLE_ORDER]),
    annotation = factor(annotation, levels = ann_levels),
    marker = factor(marker, levels = markers_use)
  )

cap99 <- quantile(dot_df$mean_expr, 0.99, na.rm = TRUE)
dot_df$mean_expr_cap <- pmin(dot_df$mean_expr, cap99)

p_dot <- ggplot(dot_df, aes(x = marker, y = annotation)) +
  geom_point(aes(size = pct_expr, color = mean_expr_cap), alpha = 0.9) +
  facet_wrap(~ sample, ncol = 2) +
  scale_size_continuous(name = "% expressing", range = c(0.5, 6)) +
  scale_color_gradient(
    low = "#ffe082",
    high = "#d32f2f",
    name = "Mean expr\n(capped q99)"
  ) +
  labs(
    title = "Canonical markers revised v3",
    subtitle = "Conditions ordered: 1wpi, 3wpi, 6wpi, mock",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8.5, color = "grey40"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank()
  )

w_dot <- max(10, length(markers_use) * 0.55 + 5)
h_dot <- max(6, length(unique(dot_df$annotation)) * 0.26 + 3)

save_fig(
  p_dot,
  OUT_SPATIAL_DIR,
  "dotplot_canonical_markers_revised_v3",
  width = w_dot,
  height = h_dot
)

# Merged version (all samples pooled into a single panel)
dot_df_merged <- do.call(rbind, lapply(markers_use, function(g) {
  vals <- as.numeric(expr_mat[g, ])
  tmp <- data.frame(
    annotation = annotation,
    expr = vals,
    marker = g,
    stringsAsFactors = FALSE
  )
  tmp %>%
    dplyr::group_by(annotation, marker) %>%
    dplyr::summarise(
      pct_expr = 100 * mean(expr > 0, na.rm = TRUE),
      mean_expr = mean(expr, na.rm = TRUE),
      .groups = "drop"
    )
}))

dot_df_merged <- dot_df_merged %>%
  dplyr::mutate(
    annotation = factor(annotation, levels = ann_levels),
    marker = factor(marker, levels = markers_use)
  )

cap99_merged <- quantile(dot_df_merged$mean_expr, 0.99, na.rm = TRUE)
dot_df_merged$mean_expr_cap <- pmin(dot_df_merged$mean_expr, cap99_merged)

p_dot_merged <- ggplot(dot_df_merged, aes(x = marker, y = annotation)) +
  geom_point(aes(size = pct_expr, color = mean_expr_cap), alpha = 0.9) +
  scale_size_continuous(name = "% expressing", range = c(0.5, 6)) +
  scale_color_gradient(
    low = "#ffe082",
    high = "#d32f2f",
    name = "Mean expr\n(capped q99)"
  ) +
  labs(
    title = "Canonical markers revised v3 (all samples merged)",
    subtitle = "All conditions pooled: 1wpi + 3wpi + 6wpi + mock",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8.5, color = "grey40"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 7),
    panel.grid.minor = element_blank()
  )

save_fig(
  p_dot_merged,
  OUT_SPATIAL_DIR,
  "dotplot_canonical_markers_revised_v3_merged",
  width = w_dot,
  height = h_dot
)

message("\n=== Script 56 completed successfully ===\n")
