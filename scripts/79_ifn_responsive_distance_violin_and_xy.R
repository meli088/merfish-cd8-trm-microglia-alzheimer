#!/usr/bin/env Rscript
# ============================================================================
# Script 79 — IFN responsive cluster: 2 figures
#
# FIGURE 1: Violin plot — distance from nearest Immune (Acod1) cell
#           by global annotation (all cell types, all samples merged)
#           → same style as existing fig from script 09
#           → annotations ordered by median distance (closest left)
#
# FIGURE 2: XY spatial maps — Immune (Acod1) in red + IFN responsive (Ifit1)
#           in blue, 4 conditions (Mock, 1wpi, 3wpi, 6wpi)
#           → centered coordinates, shared identical axis limits
#
# Input:  objects/04_banksy_joint_lam08_after_bloc3.rds
#         ncells_by_sample_lam02_res09_joint_long.csv
#
# Outputs: outputs/banksy/ifn_immune_overlay/
#   fig_ifn_nearest_immune_distance_violin.pdf/.jpg
#   fig_ifn_xy_immune_acod1_ifit1_by_condition.pdf/.jpg
# ============================================================================

set.seed(1997)

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(Banksy)
  library(FNN)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(Cairo)
})

source("scripts/00_palette.R")

obj_path  <- file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds")
csv_annot <- "ncells_by_sample_lam02_res09_joint_long.csv"
out_dir   <- file.path("outputs", "banksy", "ifn_immune_overlay")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dist_dir  <- file.path("outputs", "banksy", "nearest_immune_distance")
dir.create(dist_dir, recursive = TRUE, showWarnings = FALSE)

SAMPLE_ORDER <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi", "mock_6wpi")
SAMPLE_LABELS <- c(
  mock_6wpi = "Mock 6 wpi",
  LCMV_1wpi = "LCMV 1 wpi",
  LCMV_3wpi = "LCMV 3 wpi",
  LCMV_6wpi = "LCMV 6 wpi"
)

save_fig <- function(p, fname, width, height) {
  CairoPDF(file.path(out_dir, paste0(fname, ".pdf")), width = width, height = height)
  print(p); dev.off()
  CairoJPEG(file.path(out_dir, paste0(fname, ".jpg")),
            width = width * 170, height = height * 170, res = 170)
  print(p); dev.off()
  message("Saved: ", fname)
}

find_cl_col <- function(se, lam, res) {
  all_cols <- Banksy::clusterNames(se)
  lam_str  <- gsub("\\.", "\\\\.", as.character(lam))
  lam_cols <- all_cols[grep(paste0("lam", lam_str), all_cols)]
  res_vals <- suppressWarnings(as.numeric(sub(".*_res", "", lam_cols)))
  idx      <- which(!is.na(res_vals) & abs(res_vals - res) < 1e-8)
  lam_cols[idx[1]]
}

# ==================================================================
# LOAD
# ==================================================================
message("Loading global object...")
se <- readRDS(obj_path)
message("  Cells: ", ncol(se))

cl_col <- find_cl_col(se, 0.2, 0.9)
message("  Cluster column: ", cl_col)

cd      <- as.data.frame(colData(se))
samples <- as.character(cd$sample)
coords  <- as.data.frame(spatialCoords(se))
colnames(coords) <- c("x", "y")

# Reconstruct annotations
annot_raw <- readr::read_delim(csv_annot, delim = ";",
                               locale = locale(decimal_mark = "."),
                               show_col_types = FALSE, trim_ws = TRUE) %>%
  dplyr::select(-matches("^Unnamed")) %>%
  dplyr::filter(!is.na(annotation), annotation != "") %>%
  dplyr::mutate(banksy_domain = as.character(banksy_domain),
                annotation = trimws(annotation)) %>%
  dplyr::distinct(banksy_domain, annotation)

anno_lookup    <- setNames(annot_raw$annotation, annot_raw$banksy_domain)
domain_labels  <- paste0("Domain_", as.character(cd[[cl_col]]))
cell_annot     <- ifelse(!is.na(anno_lookup[domain_labels]) & anno_lookup[domain_labels] != "",
                         anno_lookup[domain_labels], "Non annote")

# Keep global labels and only normalize ependymal naming
cell_annot_global <- dplyr::recode(
  cell_annot,
  "Choroid plexus Epi (Enpp2)" = "Ependymal (Enpp2)",
  "Prolif neural/glial (Ccdc153)" = "Ependymal (Ccdc153)"
)

# ==================================================================
# LOAD immune object for reference masks in distance / XY overlays
# ==================================================================
message("Loading immune subtype object (08)...")
se_immune <- readRDS("objects/08_immune_annotated_lam02_res03.rds")
immune_ct <- setNames(
  as.character(SummarizedExperiment::colData(se_immune)$cell_type),
  colnames(se_immune)
)
message("  Immune subtypes: ", length(unique(immune_ct)))
immune_mask_global <- colnames(se) %in% colnames(se_immune)
message("  Immune reference cells (obj 08 IDs): ", sum(immune_mask_global))

# ==================================================================
# SECTION 1 — Nearest Immune distance violin
# ==================================================================
message("\n=== SECTION 1: Nearest Immune distance violin ===")

dist_csv <- file.path(dist_dir, "nearest_immune_distance_per_cell_globallabels_lam02_res09.csv")

if (file.exists(dist_csv)) {
  message("  Loading cached global-label distance CSV...")
  dist_df <- readr::read_csv(dist_csv, show_col_types = FALSE)
} else {
  message("  Computing distances...")
  immune_mask <- immune_mask_global
  message("  Immune reference cells (from obj 08): ", sum(immune_mask))

  all_coords_mat    <- as.matrix(coords[, c("x", "y")])
  immune_coords_mat <- all_coords_mat[immune_mask, , drop = FALSE]

  knn_res <- FNN::get.knnx(data = immune_coords_mat, query = all_coords_mat, k = 1)
  dists   <- as.numeric(knn_res$nn.dist[, 1])

  dist_df <- tibble(
    cell_id    = colnames(se),
    sample     = samples,
    annotation = cell_annot_global,
    x          = coords$x,
    y          = coords$y,
    nearest_immune_distance_um = dists
  )

  readr::write_csv(dist_df, dist_csv)
  message("  Saved: ", dist_csv)
}

# Apply label harmonization even when reading cached CSV
dist_df <- dist_df %>%
  dplyr::mutate(annotation = dplyr::recode(
    annotation,
    "Ependymal (Enpp2)" = "Choroid plexus Epi (Enpp2)",
    "Prolif neural/glial (Ccdc153)" = "Ependymal (Ccdc153)"
  ))

# Aggregate per annotation: median order (closest first)
anno_order <- dist_df %>%
  dplyr::filter(!is.na(nearest_immune_distance_um)) %>%
  dplyr::group_by(annotation) %>%
  dplyr::summarise(med = median(nearest_immune_distance_um, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(med) %>%
  dplyr::pull(annotation)

# Build colour map from shared palette, adding recent labels not in GLOBAL_PALETTE
EXTRA_COLORS <- c(
  "Choroid plexus Epi (Enpp2)"               = "#FFD700",
  "Fibroblasts/VSMC (Col1a1)"                = "#9467BD",
  "Ependymal (Ccdc153)"                      = "#BCBD22",
  "Neurons (Slc17a8)"                        = "#E7CB94",
  "Neurons (Crhbp)"                          = "#DBDB8D",
  "Neurons (Dkk3)"                           = "#9EDAE5"
)

FULL_PALETTE <- c(GLOBAL_PALETTE, IMMUNE_PALETTE, EXTRA_COLORS)

color_map <- setNames(
  rep("#CCCCCC", length(anno_order)),
  anno_order
)
for (nm in names(FULL_PALETTE)) {
  if (nm %in% anno_order) color_map[nm] <- FULL_PALETTE[nm]
}

plot_dist <- dist_df %>%
  dplyr::filter(!is.na(nearest_immune_distance_um)) %>%
  dplyr::mutate(annotation = factor(annotation, levels = anno_order))

p_violin <- ggplot(plot_dist, aes(x = annotation, y = nearest_immune_distance_um, fill = annotation)) +
  geom_violin(scale = "width", trim = FALSE, width = 0.8,
              color = "#333333", linewidth = 0.35, alpha = 0.95) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 3.2, color = "#333333") +
  scale_fill_manual(values = color_map, drop = FALSE) +
  labs(
    x = NULL,
    y = "Dist. from nearest Immune cell (\u00b5m)"
  ) +
  coord_cartesian(ylim = c(0, quantile(plot_dist$nearest_immune_distance_um, 0.995, na.rm = TRUE))) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    plot.margin = margin(10, 12, 10, 10)
  )

save_fig(p_violin, "fig_ifn_nearest_immune_distance_violin", width = 14, height = 5.2)

# ==================================================================
# SECTION 2 — XY spatial: Immune (Acod1) + IFN responsive (Ifit1)
# ==================================================================
message("\n=== SECTION 2: XY spatial maps ===")

IFN_LABEL    <- "IFN responsive (Ifit1)"
IMMUNE_LABEL <- "Immune (Acod1)"

xy_df <- tibble(
  cell_id     = colnames(se),
  x          = coords$x,
  y          = coords$y,
  sample     = samples,
  annotation = cell_annot_global,
  is_immune  = immune_mask_global
) %>%
  dplyr::filter((annotation == IFN_LABEL) | is_immune,
                sample %in% SAMPLE_ORDER) %>%
  dplyr::mutate(
    sample_label = factor(SAMPLE_LABELS[sample], levels = SAMPLE_LABELS[SAMPLE_ORDER]),
    color_group  = ifelse(annotation == IFN_LABEL, IFN_LABEL, IMMUNE_LABEL)
  )

# Build full cell df for centering (use all cells, not just highlighted)
all_df <- tibble(
  x      = coords$x,
  y      = coords$y,
  sample = samples
) %>%
  dplyr::filter(sample %in% SAMPLE_ORDER)

# Center each sample on its own median, then compute shared limits
all_df <- all_df %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(
    x_c = x - median(x, na.rm = TRUE),
    y_c = y - median(y, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

# Join centering offsets into xy_df
center_offsets <- all_df %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    x_med = median(x, na.rm = TRUE),
    y_med = median(y, na.rm = TRUE),
    .groups = "drop"
  )

xy_df <- xy_df %>%
  dplyr::left_join(center_offsets, by = "sample") %>%
  dplyr::mutate(
    x_c = x - x_med,
    y_c = y - y_med
  )

# Compute shared span from the 1st-99th quantile of all cells across all conditions
EXTRA_MARGIN <- 800

sample_spans <- sapply(SAMPLE_ORDER, function(s) {
  d <- all_df %>% dplyr::filter(sample == s)
  if (nrow(d) < 10) return(NA_real_)
  xq <- quantile(d$x_c, probs = c(0.005, 0.995), na.rm = TRUE)
  yq <- quantile(d$y_c, probs = c(0.005, 0.995), na.rm = TRUE)
  max(diff(xq), diff(yq))
})

target_half <- (max(sample_spans, na.rm = TRUE) / 2) + EXTRA_MARGIN
global_xlim <- c(-target_half, target_half)
global_ylim <- c(-target_half, target_half)

message(sprintf("  Shared axis limits: [%.0f, %.0f]", global_xlim[1], global_xlim[2]))

# Scale bar
scalebar_um <- 500
scalebar_x1 <- global_xlim[2] - 0.05 * diff(global_xlim)
scalebar_x0 <- scalebar_x1 - scalebar_um
scalebar_y  <- global_ylim[1] + 0.07 * diff(global_ylim)

scalebar_df <- tibble(
  sample_label = factor(SAMPLE_LABELS[SAMPLE_ORDER], levels = SAMPLE_LABELS[SAMPLE_ORDER]),
  x0 = scalebar_x0, x1 = scalebar_x1,
  y = scalebar_y, y_text = scalebar_y - 0.025 * diff(global_ylim)
)

HIGHLIGHT_PALETTE <- c(
  "IFN responsive (Ifit1)" = "#56B4E9",
  "Immune (Acod1)"          = "#D62728"
)

# For each panel: grey background = ALL cells in condition
panels <- lapply(SAMPLE_ORDER, function(smp) {
  slabel <- SAMPLE_LABELS[[smp]]

  df_bg <- all_df %>%
    dplyr::filter(sample == smp)
  df_fg <- xy_df %>%
    dplyr::filter(sample == smp) %>%
    dplyr::mutate(color_group = factor(color_group, levels = c("Immune (Acod1)", "IFN responsive (Ifit1)")))

  # Draw IFN first, then Immune on top
  df_ifn    <- df_fg %>% dplyr::filter(color_group == "IFN responsive (Ifit1)")
  df_immune <- df_fg %>% dplyr::filter(color_group == "Immune (Acod1)")

  n_total <- nrow(df_bg)
  pt_bg   <- if (n_total > 80000) 0.3 else 0.5

  ggplot() +
    geom_point(data = df_bg, aes(x = x_c, y = y_c), color = "grey88",
               size = pt_bg, alpha = 0.6, stroke = 0) +
    geom_point(data = df_ifn, aes(x = x_c, y = y_c), color = "#56B4E9",
               size = 0.7, alpha = 0.85, stroke = 0) +
    geom_point(data = df_immune, aes(x = x_c, y = y_c), color = "#D62728",
               size = 0.8, alpha = 0.9, stroke = 0) +
    geom_segment(
      data = scalebar_df %>% dplyr::filter(sample_label == slabel),
      aes(x = x0, xend = x1, y = y, yend = y),
      color = "black", linewidth = 0.6
    ) +
    annotate("text", x = (scalebar_x0 + scalebar_x1) / 2,
             y = scalebar_y - 0.025 * diff(global_ylim),
             label = paste0(scalebar_um, " \u00b5m"), size = 2.8, color = "black") +
    coord_fixed(ratio = 1, xlim = global_xlim, ylim = global_ylim) +
    labs(title = slabel, x = "X coordinate (\u00b5m)", y = "Y coordinate (\u00b5m)") +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 8),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white")
    )
})
names(panels) <- SAMPLE_ORDER

# Order: Mock | LCMV 1wpi | LCMV 3wpi | LCMV 6wpi
# Arrange 2x2: top-left=Mock, top-right=1wpi, bottom-left=3wpi, bottom-right=6wpi
p_final <- wrap_plots(
  panels[["mock_6wpi"]],
  panels[["LCMV_1wpi"]],
  panels[["LCMV_3wpi"]],
  panels[["LCMV_6wpi"]],
  ncol = 2
) +
  plot_annotation(
    caption = "\u25cf Immune (Acod1)      \u25cf IFN responsive (Ifit1)",
    theme = theme(plot.caption = element_text(size = 10, hjust = 0.5, color = "black"))
  )

save_fig(p_final, "fig_ifn_xy_immune_acod1_ifit1_by_condition", width = 11, height = 12)

message("\n=== Script 79 completed ===")
