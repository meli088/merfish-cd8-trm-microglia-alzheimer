#!/usr/bin/env Rscript
# =============================================================
# Script: 22_rfp_plots_readable.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-05
#
# Goal:
#   Generate interpretable, publication-style Anti.RFP_raw plots
#   from objects/03_all_clustered.rds (Seurat object), where the
#   protein signal is stored in @meta.data$Anti.RFP_raw.
#
#   Problem addressed:
#   The raw color scale is dominated by a few extreme outliers,
#   making the distribution unreadable on UMAP and spatial plots.
#   Solution: clip at quantile (q99 / q99.5) for display and
#   also provide log1p-transformed version.
#
# Outputs (folder: outputs/banksy/rfp_plots_readable/):
#   rfp_summary_stats.csv                — global distribution stats
#   rfp_summary_by_sample.csv            — per-sample stats
#   fig_rfp_umap_clipped.pdf/jpg         — UMAP, color clipped at q99
#   fig_rfp_umap_log1p.pdf/jpg           — UMAP, log1p-transformed
#   fig_rfp_umap_by_sample.pdf/jpg       — UMAP split by sample (q99 clip)
#   fig_rfp_spatial_by_sample.pdf/jpg    — spatial (one panel per sample)
#   fig_rfp_spatial_log1p_by_sample.pdf/jpg
#   fig_rfp_distribution.pdf/jpg         — violin + density per sample
# =============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

base_path <- normalizePath(".")  # Run this script from the project root directory
setwd(base_path)

out_dir <- file.path("outputs", "banksy", "rfp_plots_readable")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE,
                                      showWarnings = FALSE)

# =============================================================
# Parameters
# =============================================================

SAMPLE_ORDER  <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi", "mock_6wpi")
sample_labels <- c(mock_6wpi  = "Mock 6 wpi",
                   LCMV_1wpi  = "LCMV 1 wpi",
                   LCMV_3wpi  = "LCMV 3 wpi",
                   LCMV_6wpi  = "LCMV 6 wpi")
sample_colors <- c(mock_6wpi  = "#999999",
                   LCMV_1wpi  = "#56B4E9",
                   LCMV_3wpi  = "#E69F00",
                   LCMV_6wpi  = "#D55E00")

CLIP_QUANTILE     <- 0.99    # display clipping for colour scale
LOG1P_LABEL       <- "log1p(Anti.RFP_raw)"
CLIPPED_LABEL     <- paste0("Anti.RFP_raw\n(clipped at q",
                             CLIP_QUANTILE * 100, ")")
PT_SIZE_UMAP      <- 0.4
PT_SIZE_SPATIAL   <- 0.35
PT_ALPHA          <- 0.75

# =============================================================
# Theme helpers
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
      plot.margin      = margin(12, 16, 12, 12),
      panel.grid.major = element_line(color = "grey93", linewidth = 0.3)
    )
}

theme_spatial <- function(base_size = 10) {
  theme_void(base_size = base_size) +
    theme(
      plot.title        = element_text(face = "bold", size = base_size + 1,
                                       hjust = 0.5, margin = margin(b = 4)),
      plot.subtitle     = element_text(size = base_size - 2, color = "grey40",
                                       hjust = 0.5, lineheight = 1.3,
                                       margin = margin(b = 6)),
      strip.text        = element_text(face = "bold", size = base_size),
      strip.background  = element_rect(fill = "grey95", color = NA),
      legend.title      = element_text(size = base_size - 1, face = "bold"),
      legend.text       = element_text(size = base_size - 1.5),
      plot.margin       = margin(8, 12, 8, 8)
    )
}

rfp_gradient_clip <- function(clip_val) {
  scale_color_gradientn(
    name   = CLIPPED_LABEL,
    colors = c("grey92", "#ffffb2", "#fd8d3c", "#e31a1c", "#67000d"),
    limits = c(0, clip_val),
    oob    = scales::squish,
    na.value = "grey85",
    guide  = guide_colorbar(barwidth = 0.8, barheight = 5,
                            title.position = "top")
  )
}

rfp_gradient_log1p <- scale_color_gradientn(
  name   = LOG1P_LABEL,
  colors = c("grey92", "#ffffb2", "#fd8d3c", "#e31a1c", "#67000d"),
  na.value = "grey85",
  guide  = guide_colorbar(barwidth = 0.8, barheight = 5,
                          title.position = "top")
)

save_fig <- function(plot, name, width, height,
                     dpi_jpg = 250, cairo = TRUE) {
  pdf_path <- file.path(out_dir, paste0(name, ".pdf"))
  jpg_path <- file.path(out_dir, paste0(name, ".jpg"))
  if (cairo) {
    ggsave(pdf_path, plot = plot, width = width, height = height,
           device = cairo_pdf)
  } else {
    ggsave(pdf_path, plot = plot, width = width, height = height)
  }
  ggsave(jpg_path, plot = plot, width = width, height = height,
         dpi = dpi_jpg)
  message("  Saved: ", basename(pdf_path))
}

# =============================================================
# 1. Load the Seurat object
# =============================================================

obj_path <- file.path("objects", "03_all_clustered.rds")
stopifnot("Object file not found" = file.exists(obj_path))

message("Loading ", obj_path, " ...")
obj <- readRDS(obj_path)
message("  Loaded: ", ncol(obj), " cells, ", nrow(obj), " features")
message("  Class : ", class(obj)[1])

# Confirm Anti.RFP_raw is present
if (!"Anti.RFP_raw" %in% colnames(obj@meta.data)) {
  stop("Anti.RFP_raw not found in @meta.data. ",
       "Run 21_rfp_signal_check.R first to identify the correct object.")
}
message("  Anti.RFP_raw column confirmed in @meta.data.")

# =============================================================
# 2. Build working data frame
# =============================================================

md <- obj@meta.data
md$cell_id <- rownames(md)

if (!("sample" %in% colnames(md))) {
  stop("No 'sample' column in meta.data.")
}
md$sample <- factor(md$sample, levels = SAMPLE_ORDER)

rfp_raw <- as.numeric(md$Anti.RFP_raw)

# UMAP coordinates
has_umap <- "umap" %in% tolower(names(obj@reductions))
if (has_umap) {
  umap_name <- names(obj@reductions)[tolower(names(obj@reductions)) == "umap"][1]
  umap_emb  <- as.data.frame(obj@reductions[[umap_name]]@cell.embeddings)
  colnames(umap_emb)[1:2] <- c("UMAP_1", "UMAP_2")
  md$UMAP_1 <- umap_emb[rownames(md), "UMAP_1"]
  md$UMAP_2 <- umap_emb[rownames(md), "UMAP_2"]
  message("  UMAP reduction found: '", umap_name, "'")
} else {
  message("  No UMAP reduction found — UMAP plots will be skipped.")
}

# Spatial coordinates — loop over per-sample image slots
message("  Extracting spatial coordinates per sample...")
coord_list <- lapply(SAMPLE_ORDER, function(samp) {
  tryCatch({
    obj_sub  <- subset(obj, subset = sample == samp)
    fov_name <- Images(obj_sub)[1]
    if (is.null(fov_name) || is.na(fov_name)) return(NULL)
    coords   <- GetTissueCoordinates(obj_sub, image = fov_name)
    if (!all(c("x", "y", "cell") %in% colnames(coords))) return(NULL)
    rownames(coords) <- coords$cell
    data.frame(
      cell_id = rownames(obj_sub@meta.data),
      sdimx   = coords[rownames(obj_sub@meta.data), "x"],
      sdimy   = coords[rownames(obj_sub@meta.data), "y"],
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message("    Coord extraction failed for ", samp, ": ", conditionMessage(e))
    NULL
  })
})
names(coord_list) <- SAMPLE_ORDER

coord_df <- bind_rows(coord_list[!sapply(coord_list, is.null)])
has_spatial <- nrow(coord_df) > 0
if (has_spatial) {
  md <- left_join(md, coord_df, by = "cell_id")
  message("  Spatial coordinates merged: ",
          sum(!is.na(md$sdimx)), " / ", nrow(md), " cells")
} else {
  message("  Spatial coordinates not available.")
}

# =============================================================
# 3. Distribution statistics
# =============================================================

message("\n--- Anti.RFP_raw distribution (global) ---")

quants <- quantile(rfp_raw, probs = c(0.25, 0.5, 0.75, 0.90,
                                       0.95, 0.99, 0.995),
                   na.rm = TRUE)

stats_global <- data.frame(
  stat  = c("n_cells", "n_nonzero", "pct_above_zero",
            "min", "mean", "median", "max",
            "q25", "q50", "q75", "q90", "q95", "q99", "q99.5"),
  value = c(
    length(rfp_raw),
    sum(rfp_raw > 0, na.rm = TRUE),
    round(100 * mean(rfp_raw > 0, na.rm = TRUE), 2),
    min(rfp_raw, na.rm = TRUE),
    round(mean(rfp_raw, na.rm = TRUE), 4),
    round(median(rfp_raw, na.rm = TRUE), 4),
    max(rfp_raw, na.rm = TRUE),
    round(quants["25%"], 4),
    round(quants["50%"], 4),
    round(quants["75%"], 4),
    round(quants["90%"], 4),
    round(quants["95%"], 4),
    round(quants["99%"], 4),
    round(quants["99.5%"], 4)
  )
)

print(stats_global, row.names = FALSE)
write.csv(stats_global, file.path(out_dir, "rfp_summary_stats.csv"),
          row.names = FALSE)
message("  Saved: rfp_summary_stats.csv")

# Determine clip value for display
clip_val <- as.numeric(quants["99%"])
message("  Display clip value (q", CLIP_QUANTILE * 100, "): ", clip_val)
pct_clipped <- round(100 * mean(rfp_raw > clip_val, na.rm = TRUE), 2)
message("  Percentage of cells clipped: ", pct_clipped, "%")

q95_global <- as.numeric(quants["95%"])

# Per-sample stats
stats_by_sample <- md %>%
  group_by(sample) %>%
  summarise(
    n_cells          = n(),
    n_above_zero     = sum(Anti.RFP_raw > 0, na.rm = TRUE),
    pct_above_zero   = round(100 * mean(Anti.RFP_raw > 0, na.rm = TRUE), 2),
    n_above_q95_glob = sum(Anti.RFP_raw > q95_global, na.rm = TRUE),
    pct_above_q95    = round(100 * mean(Anti.RFP_raw > q95_global,
                                         na.rm = TRUE), 2),
    mean_rfp         = round(mean(Anti.RFP_raw, na.rm = TRUE), 4),
    median_rfp       = round(median(Anti.RFP_raw, na.rm = TRUE), 4),
    max_rfp          = max(Anti.RFP_raw, na.rm = TRUE),
    q90              = round(quantile(Anti.RFP_raw, 0.90, na.rm = TRUE), 4),
    q95              = round(quantile(Anti.RFP_raw, 0.95, na.rm = TRUE), 4),
    q99              = round(quantile(Anti.RFP_raw, 0.99, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  mutate(sample_label = sample_labels[as.character(sample)])

print(as.data.frame(stats_by_sample), row.names = FALSE)
write.csv(stats_by_sample,
          file.path(out_dir, "rfp_summary_by_sample.csv"),
          row.names = FALSE)
message("  Saved: rfp_summary_by_sample.csv")

# =============================================================
# 4. Distribution figure — violin + jitter per sample
# =============================================================

message("\n--- Fig: Distribution per sample ---")

dist_df <- md %>%
  select(sample, Anti.RFP_raw) %>%
  mutate(rfp_log1p = log1p(Anti.RFP_raw))

# Violin — raw (clipped display)
p_violin_raw <- ggplot(
    dist_df,
    aes(x = sample, y = pmin(Anti.RFP_raw, clip_val),
        fill = sample, color = sample)
  ) +
  geom_violin(alpha = 0.55, linewidth = 0.4, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.8,
               linewidth = 0.4, color = "grey30",
               fill = "white") +
  scale_fill_manual(values  = sample_colors, guide = "none") +
  scale_color_manual(values = sample_colors, guide = "none") +
  scale_x_discrete(labels = sample_labels) +
  labs(
    title    = "Anti.RFP_raw distribution by sample",
    subtitle = paste0("Color scale clipped at q", CLIP_QUANTILE * 100,
                      " = ", round(clip_val, 2),
                      " (", pct_clipped, "% of cells above)"),
    x = NULL,
    y = paste0("Anti.RFP_raw (clipped at q", CLIP_QUANTILE * 100, ")")
  ) +
  theme_pub()

# Violin — log1p
p_violin_log <- ggplot(
    dist_df,
    aes(x = sample, y = rfp_log1p, fill = sample, color = sample)
  ) +
  geom_violin(alpha = 0.55, linewidth = 0.4, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.8,
               linewidth = 0.4, color = "grey30",
               fill = "white") +
  scale_fill_manual(values  = sample_colors, guide = "none") +
  scale_color_manual(values = sample_colors, guide = "none") +
  scale_x_discrete(labels = sample_labels) +
  labs(
    title    = "log1p(Anti.RFP_raw) distribution by sample",
    subtitle = "log1p transformation to stabilize variance",
    x = NULL,
    y = LOG1P_LABEL
  ) +
  theme_pub()

p_dist <- p_violin_raw / p_violin_log +
  plot_annotation(
    title    = "Anti.RFP_raw — Distribution overview",
    subtitle = paste0("objects/03_all_clustered.rds | ",
                      nrow(md), " cells"),
    theme    = theme(plot.title = element_text(face = "bold", size = 13),
                     plot.subtitle = element_text(size = 9, color = "grey40"))
  )

save_fig(p_dist, "fig_rfp_distribution", width = 9, height = 9)

# =============================================================
# 5. UMAP plots
# =============================================================

if (has_umap) {
  message("\n--- Fig: UMAP plots ---")

  umap_df <- md %>%
    filter(!is.na(UMAP_1), !is.na(UMAP_2)) %>%
    arrange(Anti.RFP_raw)   # plot low values first so high values visible

  # 5a. UMAP clipped at q99
  p_umap_clip <- ggplot(umap_df,
      aes(x = UMAP_1, y = UMAP_2,
          color = pmin(Anti.RFP_raw, clip_val))) +
    geom_point(size = PT_SIZE_UMAP, alpha = PT_ALPHA, shape = 16) +
    rfp_gradient_clip(clip_val) +
    labs(
      title    = "Anti.RFP_raw — UMAP (clipped colour scale)",
      subtitle = paste0("Display clipped at q", CLIP_QUANTILE * 100,
                        " = ", round(clip_val, 2),
                        " | raw values preserved in CSV"),
      x = "UMAP 1", y = "UMAP 2"
    ) +
    theme_pub() +
    theme(axis.text = element_blank(), axis.ticks = element_blank())

  save_fig(p_umap_clip, "fig_rfp_umap_clipped", width = 7, height = 6)

  # 5b. UMAP log1p
  p_umap_log <- ggplot(umap_df,
      aes(x = UMAP_1, y = UMAP_2, color = log1p(Anti.RFP_raw))) +
    geom_point(size = PT_SIZE_UMAP, alpha = PT_ALPHA, shape = 16) +
    rfp_gradient_log1p +
    labs(
      title    = "Anti.RFP_raw — UMAP (log1p scale)",
      subtitle = "Colour = log1p(Anti.RFP_raw); no clipping applied",
      x = "UMAP 1", y = "UMAP 2"
    ) +
    theme_pub() +
    theme(axis.text = element_blank(), axis.ticks = element_blank())

  save_fig(p_umap_log, "fig_rfp_umap_log1p", width = 7, height = 6)

  # 5c. UMAP split by sample (clipped)
  p_umap_split <- ggplot(umap_df,
      aes(x = UMAP_1, y = UMAP_2,
          color = pmin(Anti.RFP_raw, clip_val))) +
    geom_point(size = PT_SIZE_UMAP * 0.85, alpha = PT_ALPHA, shape = 16) +
    rfp_gradient_clip(clip_val) +
    facet_wrap(~ sample, nrow = 2,
               labeller = as_labeller(sample_labels)) +
    labs(
      title    = "Anti.RFP_raw — UMAP by sample (clipped colour scale)",
      subtitle = paste0("Display clipped at q", CLIP_QUANTILE * 100,
                        " = ", round(clip_val, 2)),
      x = "UMAP 1", y = "UMAP 2"
    ) +
    theme_pub() +
    theme(
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      strip.text       = element_text(face = "bold", size = 9),
      strip.background = element_rect(fill = "grey95", color = NA)
    )

  save_fig(p_umap_split, "fig_rfp_umap_by_sample",
           width = 11, height = 9)
}

# =============================================================
# 6. Spatial plots
# =============================================================

if (has_spatial) {
  message("\n--- Fig: Spatial plots ---")

  spatial_df <- md %>%
    filter(!is.na(sdimx), !is.na(sdimy)) %>%
    arrange(Anti.RFP_raw)  # low values first

  # 6a. Spatial — clipped, faceted by sample
  p_spatial_clip <- ggplot(spatial_df,
      aes(x = sdimx, y = sdimy,
          color = pmin(Anti.RFP_raw, clip_val))) +
    geom_point(size = PT_SIZE_SPATIAL, alpha = PT_ALPHA, shape = 16) +
    rfp_gradient_clip(clip_val) +
    facet_wrap(~ sample, nrow = 2,
               labeller = as_labeller(sample_labels),
               scales = "free") +
    labs(
      title    = "Anti.RFP_raw — Spatial distribution by sample",
      subtitle = paste0("Display clipped at q", CLIP_QUANTILE * 100,
                        " = ", round(clip_val, 2),
                        " | Coordinates in µm"),
      x = "X (µm)", y = "Y (µm)"
    ) +
    theme_spatial() +
    theme(
      axis.title = element_text(size = 8, color = "grey40"),
      axis.text  = element_text(size = 6, color = "grey60")
    )

  save_fig(p_spatial_clip, "fig_rfp_spatial_by_sample",
           width = 13, height = 10)

  # 6b. Spatial — log1p, faceted by sample
  p_spatial_log <- ggplot(spatial_df,
      aes(x = sdimx, y = sdimy, color = log1p(Anti.RFP_raw))) +
    geom_point(size = PT_SIZE_SPATIAL, alpha = PT_ALPHA, shape = 16) +
    rfp_gradient_log1p +
    facet_wrap(~ sample, nrow = 2,
               labeller = as_labeller(sample_labels),
               scales = "free") +
    labs(
      title    = "Anti.RFP_raw — Spatial distribution by sample (log1p scale)",
      subtitle = "Colour = log1p(Anti.RFP_raw); coordinates in µm",
      x = "X (µm)", y = "Y (µm)"
    ) +
    theme_spatial() +
    theme(
      axis.title = element_text(size = 8, color = "grey40"),
      axis.text  = element_text(size = 6, color = "grey60")
    )

  save_fig(p_spatial_log, "fig_rfp_spatial_log1p_by_sample",
           width = 13, height = 10)

  # 6c. Per-sample spatial panels with individual scales (better readability)
  message("  Generating per-sample spatial panels with individual scales...")

  panel_plots <- lapply(SAMPLE_ORDER, function(samp) {
    df_s <- spatial_df %>% filter(sample == samp)
    if (nrow(df_s) == 0) return(NULL)

    clip_s <- quantile(df_s$Anti.RFP_raw, CLIP_QUANTILE, na.rm = TRUE)
    n_pos  <- sum(df_s$Anti.RFP_raw > 0, na.rm = TRUE)
    pct_pos <- round(100 * n_pos / nrow(df_s), 1)

    ggplot(df_s,
        aes(x = sdimx, y = sdimy,
            color = pmin(Anti.RFP_raw, clip_s))) +
      geom_point(size = PT_SIZE_SPATIAL * 0.9, alpha = PT_ALPHA,
                 shape = 16) +
      scale_color_gradientn(
        name   = CLIPPED_LABEL,
        colors = c("grey92", "#ffffb2", "#fd8d3c", "#e31a1c", "#67000d"),
        limits = c(0, clip_s),
        oob    = scales::squish,
        na.value = "grey85",
        guide  = guide_colorbar(barwidth = 0.6, barheight = 4)
      ) +
      coord_equal() +
      labs(
        title    = sample_labels[samp],
        subtitle = paste0(nrow(df_s), " cells | ",
                          n_pos, " (", pct_pos, "%) RFP > 0\n",
                          "Clip: q", CLIP_QUANTILE * 100,
                          " = ", round(clip_s, 1))
      ) +
      theme_spatial(base_size = 9) +
      theme(plot.title = element_text(size = 10, face = "bold",
                                      color = sample_colors[samp]))
  })

  panel_plots <- panel_plots[!sapply(panel_plots, is.null)]

  if (length(panel_plots) > 0) {
    p_panels <- wrap_plots(panel_plots, nrow = 2) +
      plot_annotation(
        title    = "Anti.RFP_raw — Spatial (per-sample scale, clipped at q99)",
        subtitle = "Each panel uses its own q99 clip for maximum local contrast",
        theme    = theme(
          plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9, color = "grey40")
        )
      )
    save_fig(p_panels, "fig_rfp_spatial_panels_indiv_scale",
             width = 14, height = 11)
  }
}

# =============================================================
# Done
# =============================================================

message("\n============================================================")
message("All outputs saved in: ", out_dir)
message("Files produced:")
message("  rfp_summary_stats.csv")
message("  rfp_summary_by_sample.csv")
message("  fig_rfp_distribution.pdf/jpg")
if (has_umap) {
  message("  fig_rfp_umap_clipped.pdf/jpg")
  message("  fig_rfp_umap_log1p.pdf/jpg")
  message("  fig_rfp_umap_by_sample.pdf/jpg")
}
if (has_spatial) {
  message("  fig_rfp_spatial_by_sample.pdf/jpg")
  message("  fig_rfp_spatial_log1p_by_sample.pdf/jpg")
  message("  fig_rfp_spatial_panels_indiv_scale.pdf/jpg")
}
message("============================================================\n")
