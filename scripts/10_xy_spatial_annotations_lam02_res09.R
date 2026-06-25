#!/usr/bin/env Rscript
# ============================================================================
# STEP: XY Spatial Maps by Annotation (lambda = 0.2, res = 0.9)
# ============================================================================

library(SingleCellExperiment)
library(SummarizedExperiment)
library(SpatialExperiment)
library(tidyverse)
library(ggplot2)
library(grid)
library(patchwork)

base_path <- normalizePath(".")  # Run this script from the project root directory
setwd(base_path)

# Palette et ordre des annotations (fichier partagé)
source("scripts/00_palette.R")

obj_path <- file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds")
csv_annot <- "ncells_by_sample_lam02_res09_joint_long.csv"
out_dir <- file.path("outputs", "banksy", "spatial_annotations")
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

safe_slug <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

cat("\n=== XY Spatial Maps by Annotation ===\n")
cat("Using object:", obj_path, "\n")

# ------------------------------------------------------------------
# Load object
# ------------------------------------------------------------------
se <- readRDS(obj_path)
cluster_col <- find_cl_col(se, 0.2, 0.9)

cd <- as.data.frame(SummarizedExperiment::colData(se))
banksy_clusters <- as.numeric(cd[[cluster_col]])
samples <- as.character(cd$sample)

spatial_mat <- as.data.frame(SpatialExperiment::spatialCoords(se))
if (!all(c("sdimx", "sdimy") %in% colnames(spatial_mat))) {
  stop("Spatial coordinates must have columns 'sdimx' and 'sdimy'")
}

# ------------------------------------------------------------------
# Load annotation mapping
# ------------------------------------------------------------------
anno_data <- read.delim(csv_annot, sep = ";", stringsAsFactors = FALSE)

anno_map <- anno_data %>%
  filter(annotation != "" & !is.na(annotation)) %>%
  select(banksy_domain, annotation) %>%
  distinct()

anno_map$cluster_id <- as.numeric(gsub("Domain_", "", anno_map$banksy_domain))
anno_lookup <- setNames(anno_map$annotation, anno_map$cluster_id)

cell_annotations <- anno_lookup[as.character(banksy_clusters)]
cell_annotations[is.na(cell_annotations)] <- "Non annote"

# ------------------------------------------------------------------
# Build plotting data
# ------------------------------------------------------------------
plot_df <- tibble(
  x = spatial_mat$sdimx,
  y = spatial_mat$sdimy,
  sample = samples,
  annotation = cell_annotations
)

sample_order <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi", "mock_6wpi")
samples_present <- unique(plot_df$sample)
sample_levels <- c(
  sample_order[sample_order %in% samples_present],
  setdiff(samples_present, sample_order)
)
plot_df$sample <- factor(plot_df$sample, levels = sample_levels)

all_annotations <- order_annotations(unique(plot_df$annotation))

# Shared XY limits across sample panels (same spatial scale within each figure)
global_x_range <- range(plot_df$x, na.rm = TRUE)
global_y_range <- range(plot_df$y, na.rm = TRUE)
global_x_pad <- max(diff(global_x_range) * 0.08, 50)
global_y_pad <- max(diff(global_y_range) * 0.08, 50)
global_xlim <- c(global_x_range[1] - global_x_pad, global_x_range[2] + global_x_pad)
global_ylim <- c(global_y_range[1] - global_y_pad, global_y_range[2] + global_y_pad)

cat("Cells:", nrow(plot_df), "\n")
cat("Annotations:", length(all_annotations), "\n")
cat("Using shared XY limits across sample panels with fixed micron aspect ratio.\n")

# Save annotation counts
annot_counts <- plot_df %>%
  count(annotation, sample, name = "n_cells")

write.csv(
  annot_counts,
  file.path(out_dir, "xy_spatial_annotation_counts_lam02_res09.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Compute convex hull (unified gray layer) per sample for tissue silhouette
# ------------------------------------------------------------------
hull_df <- do.call(rbind, lapply(levels(plot_df$sample), function(s) {
  d <- plot_df %>% filter(sample == s)
  if (nrow(d) < 3) return(NULL)
  h <- d[chull(d$x, d$y), c("x", "y")]
  h$sample <- s
  h
}))
if (!is.null(hull_df) && nrow(hull_df) > 0) {
  hull_df$sample <- factor(hull_df$sample, levels = levels(plot_df$sample))
} else {
  hull_df <- tibble(x = numeric(), y = numeric(), sample = factor(levels(plot_df$sample)))
}

# ------------------------------------------------------------------
# Plot one figure per annotation
# ------------------------------------------------------------------
for (a in all_annotations) {
  message("Plotting annotation: ", a)

  df_sel <- plot_df %>% filter(annotation == a)
  slug <- safe_slug(a)

  n_sel <- nrow(df_sel)
  message("  selected cells: ", n_sel)

  sample_plots <- lapply(levels(plot_df$sample), function(sample_name) {
    df_sample <- plot_df %>% filter(sample == sample_name)
    df_other <- df_sample %>% filter(annotation != a)
    df_sel_sample <- df_sample %>% filter(annotation == a)

    if (nrow(df_sample) == 0) {
      return(NULL)
    }

    hull_sample <- hull_df %>% filter(sample == sample_name)

    ggplot() +
      (if (nrow(hull_sample) > 0) geom_polygon(
        data = hull_sample,
        aes(x = x, y = y),
        fill = "#f8f8f8",
        color = NA
      ) else NULL) +
      geom_point(
        data = df_sel_sample,
        aes(x = x, y = y),
        color = "#D7191C",
        size = 0.20,
        alpha = 0.98,
        stroke = 0
      ) +
      coord_fixed(
        xlim = global_xlim,
        ylim = global_ylim,
        expand = FALSE
      ) +
      labs(
        title = sample_name,
        x = "X coordinate (um)",
        y = "Y coordinate (um)"
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.ticks = element_line(linewidth = 0.4, color = "black"),
        plot.margin = margin(6, 8, 6, 8)
      )
  })

  sample_plots <- sample_plots[!vapply(sample_plots, is.null, logical(1))]
  if (length(sample_plots) == 0) {
    next
  }

  p <- wrap_plots(sample_plots, ncol = 2) +
    plot_annotation(title = paste0("Spatial distribution of ", a, " by sample")) &
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
    )

  pdf_file <- file.path(
    out_dir,
    paste0("xy_spatial_", slug, "_by_sample_lam02_res09.pdf")
  )
  jpg_file <- file.path(
    out_dir,
    paste0("xy_spatial_", slug, "_by_sample_lam02_res09.jpg")
  )

  ggsave(pdf_file, p, width = 9.2, height = 8.0, device = cairo_pdf)
  ggsave(jpg_file, p, width = 9.2, height = 8.0, dpi = 450)
}

cat("\nSaved outputs to:", out_dir, "\n")
cat("Per-annotation plots generated:", length(all_annotations), "\n")
cat("Each annotation figure contains 4 sample panels with fixed micron aspect ratio.\n")

# ------------------------------------------------------------------
# Dot plot: canonical cell-type markers by annotation
# ------------------------------------------------------------------
cat("\n=== Dot plot: canonical markers by annotation ===\n")

canonical_markers <- c(
  "Egr1", "Gad2", "Satb2", "Chat", "Pdgfra", "Bcas1", "Laptm5",
  "Gja1", "Pecam1", "Vtn", "Dcn", "Myh11", "Vcan", "Ptprc",
  "Cd8a", "Cd3e", "Tmem119", "P2ry12", "Mog", "Olig2", "Gfap", "Stat1",
  "Pip1", "Slc17a8", "Rxfp1", "Neurod1", "Nefm", "Htr2c", "Dkk3",
  "Crhbp", "Arhgap36", "Adora2a", "Acod1", "Ifit1", "Col1a1", "Fgfr3", "Plp1", "Enpp2"
)
# "Plp1" "Enpp2"

# Pick the best available normalized assay
avail_assays <- SummarizedExperiment::assayNames(se)
assay_use <- if ("logcounts" %in% avail_assays) {
  "logcounts"
} else if ("normcounts" %in% avail_assays) {
  "normcounts"
} else {
  avail_assays[1]
}
cat("Using assay:", assay_use, "\n")

expr_mat <- SummarizedExperiment::assay(se, assay_use)

# Keep only markers present in the object
markers_use <- canonical_markers[canonical_markers %in% rownames(expr_mat)]
missing_markers <- setdiff(canonical_markers, markers_use)
if (length(missing_markers) > 0) {
  cat("Markers not found in object:", paste(missing_markers, collapse = ", "), "\n")
}
cat("Markers found:", length(markers_use), "/", length(canonical_markers), "\n")

# Build per-annotation summary (mean expression + fraction expressed)
annotations_vec <- cell_annotations  # already computed above

dot_data <- do.call(rbind, lapply(markers_use, function(gene) {
  expr <- as.numeric(expr_mat[gene, ])
  do.call(rbind, lapply(unique(annotations_vec), function(ann) {
    idx <- annotations_vec == ann
    e <- expr[idx]
    data.frame(
      gene       = gene,
      annotation = ann,
      mean_expr  = mean(e, na.rm = TRUE),
      pct_expr   = mean(e > 0, na.rm = TRUE) * 100,
      stringsAsFactors = FALSE
    )
  }))
}))

# Order genes as supplied, annotations in ANNOTATION_ORDER (Unknown last)
dot_data$gene <- factor(dot_data$gene, levels = markers_use)
ann_levels <- order_annotations(unique(dot_data$annotation), extended = TRUE)
ann_levels <- ann_levels[ann_levels %in% unique(dot_data$annotation)]
dot_data$annotation <- factor(dot_data$annotation, levels = ann_levels)

# Cap mean expression at the 99th percentile to prevent outlier markers
# from dominating the color scale
cap_val <- quantile(dot_data$mean_expr, 0.99, na.rm = TRUE)
dot_data$mean_expr_capped <- pmin(dot_data$mean_expr, cap_val)
cat(sprintf("Color scale capped at 99th percentile: %.3f\n", cap_val))

p_dot <- ggplot(dot_data, aes(x = gene, y = annotation)) +
  geom_point(aes(size = pct_expr, color = mean_expr_capped)) +
  scale_color_gradient(low = "#ffee01", high = "#ff001e",
                       name = "Mean expr.\n(capped 99p)") +
  scale_size_continuous(range = c(0.5, 6), name = "% cells\nexpressing") +
  labs(
    title = "Canonical cell-type markers by annotation",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y  = element_text(size = 9),
    plot.title   = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
  )

dot_pdf <- file.path(out_dir, "dotplot_canonical_markers_by_annotation_lam02_res09.pdf")
dot_jpg <- file.path(out_dir, "dotplot_canonical_markers_by_annotation_lam02_res09.jpg")

n_ann  <- length(levels(dot_data$annotation))
n_gene <- length(markers_use)
fig_w  <- max(8, n_gene * 0.42 + 3)
fig_h  <- max(4, n_ann  * 0.38 + 2)

ggsave(dot_pdf, p_dot, width = fig_w, height = fig_h, device = cairo_pdf)
ggsave(dot_jpg, p_dot, width = fig_w, height = fig_h, dpi = 300)
cat("Dot plot saved to:", dot_pdf, "\n")

cat("Done.\n")