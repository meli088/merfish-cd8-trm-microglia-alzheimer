#!/usr/bin/env Rscript

set.seed(1997)

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Cairo)
})

base_path <- normalizePath(".")
setwd(base_path)

source("scripts/00_palette.R")

OBJ_FILE <- file.path("objects", "08_immune_annotated_lam02_res03.rds")
OUT_DIR <- file.path("outputs", "banksy", "immune_acod1", "analysis", "figures")
OUT_BASE <- file.path(OUT_DIR, "xy_zoom_niche_highlighted_immune_subclusters")

TIME_ORDER <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")
FOCUS_TYPE <- "Cycling CD8 T cells (Gzma)"
KEY_TYPES <- c(
  "Cycling CD8 T cells (Gzma)",
  "Regulatory T cells (Foxp3)",
  "Activated microglia (C1qa)"
)

# Zoom half-window in micrometers (user requested ~400-500 um)
HALF_WINDOW_UM <- 1000
SCALE_BAR_UM <- 100

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(OBJ_FILE)) {
  stop("Missing object: ", OBJ_FILE)
}

resolve_annotation_color <- function(label) {
  if (label %in% names(IMMUNE_PALETTE)) {
    return(unname(IMMUNE_PALETTE[[label]]))
  }
  if (label %in% names(GLOBAL_PALETTE)) {
    return(unname(GLOBAL_PALETTE[[label]]))
  }
  "#6F6F6F"
}

message("Loading: ", OBJ_FILE)
se <- readRDS(OBJ_FILE)

cd <- as.data.frame(SummarizedExperiment::colData(se))
sc <- as.data.frame(SpatialExperiment::spatialCoords(se))

if (!all(c("sample", "cell_type") %in% colnames(cd))) {
  stop("colData must contain 'sample' and 'cell_type'")
}
if (!all(c("sdimx", "sdimy") %in% colnames(sc))) {
  stop("spatialCoords must contain 'sdimx' and 'sdimy'")
}

df <- dplyr::bind_cols(
  cd %>% transmute(
    sample = as.character(sample),
    cell_type = as.character(cell_type)
  ),
  sc %>% transmute(x = as.numeric(sdimx), y = as.numeric(sdimy))
) %>%
  filter(
    sample %in% TIME_ORDER,
    !is.na(cell_type), cell_type != "",
    !is.na(x), !is.na(y)
  )

if (nrow(df) == 0) {
  stop("No cells left after filtering to requested samples")
}

centroids <- df %>%
  filter(cell_type == FOCUS_TYPE) %>%
  group_by(sample) %>%
  summarise(cx = mean(x), cy = mean(y), n = n(), .groups = "drop")

missing_tp <- setdiff(TIME_ORDER, centroids$sample)
if (length(missing_tp) > 0) {
  stop(
    "Missing Cycling CD8 cells in sample(s): ",
    paste(missing_tp, collapse = ", ")
  )
}

key_colors <- setNames(vapply(KEY_TYPES, resolve_annotation_color, character(1)), KEY_TYPES)

make_zoom_plot <- function(tp) {
  ctr <- centroids %>% filter(sample == tp)

  xlim_tp <- c(ctr$cx - HALF_WINDOW_UM, ctr$cx + HALF_WINDOW_UM)
  ylim_tp <- c(ctr$cy - HALF_WINDOW_UM, ctr$cy + HALF_WINDOW_UM)

  # Light-gray underlay to suggest tissue/cortex contours in the zoom window.
  d_bg_contour <- df %>%
    filter(sample == tp, x >= xlim_tp[1], x <= xlim_tp[2], y >= ylim_tp[1], y <= ylim_tp[2])

  d_tp <- d_bg_contour %>%
    filter(sample == tp, x >= xlim_tp[1], x <= xlim_tp[2], y >= ylim_tp[1], y <= ylim_tp[2]) %>%
    mutate(
      display = dplyr::case_when(
        cell_type == "Cycling CD8 T cells (Gzma)" ~ "Cycling CD8 T cells (Gzma)",
        cell_type == "Regulatory T cells (Foxp3)" ~ "Regulatory T cells (Foxp3)",
        cell_type == "Activated microglia (C1qa)" ~ "Activated microglia (C1qa)",
        TRUE ~ "Other"
      ),
      display = factor(display, levels = c("Other", KEY_TYPES)),
      pt_size = dplyr::case_when(
        display %in% KEY_TYPES ~ 0.82,
        TRUE ~ 0.58
      )
    )

  d_key <- d_tp %>% filter(display != "Other")

  sb_margin <- 40
  x1 <- xlim_tp[2] - sb_margin
  x0 <- x1 - SCALE_BAR_UM
  y0 <- ylim_tp[1] + sb_margin

  ggplot() +
    geom_point(
      data = d_bg_contour,
      aes(x = x, y = y),
      color = "#E6E6E6",
      size = 0.62,
      alpha = 0.70,
      stroke = 0
    ) +
    geom_point(
      data = d_tp,
      aes(x = x, y = y),
      color = "#8E8E8E",
      size = 0.56,
      alpha = 0.95,
      stroke = 0
    ) +
    geom_point(
      data = d_key,
      aes(x = x, y = y, color = display, size = pt_size),
      alpha = 0.9,
      stroke = 0
    ) +
    geom_segment(
      aes(x = x0, xend = x1, y = y0, yend = y0),
      inherit.aes = FALSE,
      linewidth = 0.9,
      color = "black"
    ) +
    annotate(
      "text",
      x = (x0 + x1) / 2,
      y = y0 + 22,
      label = paste0(SCALE_BAR_UM, " um"),
      size = 3,
      vjust = 0,
      color = "black"
    ) +
    scale_color_manual(values = key_colors, breaks = KEY_TYPES, drop = FALSE, name = "Annotation") +
    scale_size_identity(guide = "none") +
    coord_fixed(xlim = xlim_tp, ylim = ylim_tp, expand = FALSE) +
    labs(
      title = tp,
      x = "X (um)",
      y = "Y (um)"
    ) +
    theme_classic(base_size = 10) +
    theme(
      panel.background = element_rect(fill = "#F5F5F5", color = NA),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      plot.subtitle = element_text(size = 8.5, hjust = 0.5, color = "grey40"),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm"),
      panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.4)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
}

plots <- lapply(TIME_ORDER, make_zoom_plot)
p_final <- wrap_plots(plots, nrow = 1, guides = "collect") +
  plot_annotation(
    title = "XY zoom niche highlighted (immune subclusters)",
    subtitle = sprintf("Fenetre +/- %d um autour du centroide Cycling CD8 (CD8/T CD4/microglies colories)", HALF_WINDOW_UM)
  ) &
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12))

cairo_pdf(paste0(OUT_BASE, ".pdf"), width = 14, height = 5.4)
print(p_final)
dev.off()

CairoJPEG(
  paste0(OUT_BASE, ".jpg"),
  width = round(14 * 180),
  height = round(5.4 * 180),
  res = 180,
  quality = 95
)
print(p_final)
dev.off()

message("Saved: ", paste0(OUT_BASE, ".pdf"))
message("Saved: ", paste0(OUT_BASE, ".jpg"))
