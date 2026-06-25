#!/usr/bin/env Rscript
# ============================================================================
# Script 41 — xy_spatial_all_annotations_combined (revised)
#
# Problèmes corrigés vs version originale :
#   - Points trop petits → taille augmentée
#   - Légende illisible (trop de types) → regroupement des neurones +
#     affichage uniquement des grandes populations
#   - Immune (Acod1) mis au premier plan en rouge vif, taille +50%
#
# Layout : facet 2×2 (un panneau par échantillon)
# Populations affichées dans la légende :
#   Immune (Acod1), IFN responsive (Ifit1), Microglia (P2ry12),
#   Astrocytes (Fgfr3), Astrocytes (Gfap), Oligodendrocytes (Plp1),
#   OPC (Pdgfra), Vascular (Cldn5), Neurons (tous sous-types),
#   Autres (toutes autres annotations en gris)
#
# Output : outputs/banksy/spatial_annotations/
#   xy_spatial_all_annotations_combined_lam02_res09.pdf/jpg
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

obj_path  <- file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds")
csv_annot <- "ncells_by_sample_lam02_res09_joint_long.csv"
out_dir   <- file.path("outputs", "banksy", "spatial_annotations")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
find_cl_col <- function(se, lam, res) {
  all_cols <- clusterNames(se)
  lam_str  <- gsub("\\.", "\\\\.", as.character(lam))
  lam_cols <- all_cols[grep(paste0("lam", lam_str), all_cols)]
  if (length(lam_cols) == 0) stop("Cluster column not found for lam=", lam)
  res_vals <- suppressWarnings(as.numeric(sub(".*_res", "", lam_cols)))
  idx      <- which(!is.na(res_vals) & abs(res_vals - res) < 1e-8)
  if (length(idx) == 0) stop("Cluster column not found for res=", res)
  lam_cols[idx[1]]
}

# ------------------------------------------------------------------
# Chargement de l'objet + reconstruction des annotations
# ------------------------------------------------------------------
cat("\n=== Chargement de l'objet ===\n")
se <- readRDS(obj_path)
cat("Cells:", ncol(se), "\n")

cl_col <- find_cl_col(se, 0.2, 0.9)
cat("Cluster column:", cl_col, "\n")

cd              <- as.data.frame(SummarizedExperiment::colData(se))
banksy_clusters <- as.numeric(cd[[cl_col]])
samples_vec     <- as.character(cd$sample)

spatial_mat <- as.data.frame(SpatialExperiment::spatialCoords(se))
colnames(spatial_mat) <- c("x", "y")

anno_data <- read.delim(csv_annot, sep = ";", stringsAsFactors = FALSE)
anno_map  <- anno_data %>%
  dplyr::filter(annotation != "" & !is.na(annotation)) %>%
  dplyr::select(banksy_domain, annotation) %>%
  dplyr::distinct()
anno_map$cluster_id <- as.numeric(gsub("Domain_", "", anno_map$banksy_domain))
anno_lookup         <- setNames(trimws(anno_map$annotation), anno_map$cluster_id)

cell_annotations <- anno_lookup[as.character(banksy_clusters)]
cell_annotations[is.na(cell_annotations)] <- "Non annote"

sample_order <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi", "mock_6wpi")
SAMPLE_LABELS <- c(
  mock_6wpi = "Mock 6 wpi",
  LCMV_1wpi = "LCMV 1 wpi",
  LCMV_3wpi = "LCMV 3 wpi",
  LCMV_6wpi = "LCMV 6 wpi"
)

SCALE_BAR_UM <- 500
ZOOM_MIN_SPAN_UM <- 5000
EXTRA_MARGIN_UM <- 1000

plot_df <- data.frame(
  x          = spatial_mat$x,
  y          = spatial_mat$y,
  sample     = samples_vec,
  annotation = cell_annotations,
  stringsAsFactors = FALSE
)

samples_present <- unique(plot_df$sample)
sample_levels   <- c(
  sample_order[sample_order %in% samples_present],
  setdiff(samples_present, sample_order)
)
plot_df$sample_label <- factor(
  SAMPLE_LABELS[plot_df$sample],
  levels = SAMPLE_LABELS[sample_levels]
)

# Recentrer chaque échantillon et appliquer une fenêtre de zoom partagée
# (mêmes échelles + ratio spatial conservé)
plot_df <- plot_df %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(
    x_plot = x - median(x, na.rm = TRUE),
    y_plot = y - median(y, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

sample_spans <- sapply(sample_levels, function(smp) {
  d <- plot_df %>% dplyr::filter(sample == smp)
  if (nrow(d) < 10) return(NA_real_)
  xq <- quantile(d$x_plot, probs = c(0.01, 0.99), na.rm = TRUE)
  yq <- quantile(d$y_plot, probs = c(0.01, 0.99), na.rm = TRUE)
  max(diff(xq), diff(yq))
})
target_span <- max(ZOOM_MIN_SPAN_UM, as.numeric(stats::quantile(sample_spans, 0.90, na.rm = TRUE)) * 1.05)
half_span <- target_span / 2
global_xlim <- c(-half_span - EXTRA_MARGIN_UM, half_span + EXTRA_MARGIN_UM)
global_ylim <- c(-half_span - EXTRA_MARGIN_UM, half_span + EXTRA_MARGIN_UM)

scalebar_x1 <- global_xlim[2] - 0.08 * diff(global_xlim)
scalebar_x0 <- scalebar_x1 - SCALE_BAR_UM
scalebar_y  <- global_ylim[1] + 0.10 * diff(global_ylim)
scalebar_y_text <- scalebar_y - 0.03 * diff(global_ylim)

scalebar_df <- data.frame(
  sample = sample_levels,
  x0 = scalebar_x0,
  x1 = scalebar_x1,
  y = scalebar_y,
  y_text = scalebar_y_text,
  label = paste0(SCALE_BAR_UM, " µm"),
  stringsAsFactors = FALSE
)

cat("Annotations dans l'objet:\n")
print(sort(table(plot_df$annotation), decreasing = TRUE))

# ------------------------------------------------------------------
# Regroupement des annotations
# ------------------------------------------------------------------
# Patterns de sous-types neuronaux à regrouper
NEURON_PATTERNS <- c(
  "Neurons \\(", "Neurons\\(",           # e.g. "Neurons (Nefm)"
  "Excitatory neurons", "Inhibitory neurons",
  "Neurons \\(Fam107a\\)", "Neurons \\(Rbfox3\\)"
)

# Grandes populations à conserver dans la légende avec leur couleur d'origine
KEEP_LABELS <- c(
  "Immune (Acod1)",
  "IFN responsive (Ifit1)",
  "Microglia (P2ry12)",
  "Astrocytes (Fgfr3)",
  "Astrocytes (Gfap)",
  "Oligodendrocytes (Plp1)",
  "OPC (Pdgfra)",
  "Vascular (Cldn5)"
)

is_neuron <- grepl(paste(NEURON_PATTERNS, collapse = "|"), plot_df$annotation)

plot_df$display_annotation <- dplyr::case_when(
  plot_df$annotation %in% KEEP_LABELS ~ plot_df$annotation,
  is_neuron                            ~ "Neurons",
  TRUE                                 ~ "Autres"
)

# ------------------------------------------------------------------
# Palette pour les display_annotations
# ------------------------------------------------------------------
DISPLAY_LEVELS <- c(
  "Immune (Acod1)",
  "IFN responsive (Ifit1)",
  "Microglia (P2ry12)",
  "Astrocytes (Fgfr3)",
  "Astrocytes (Gfap)",
  "Oligodendrocytes (Plp1)",
  "OPC (Pdgfra)",
  "Vascular (Cldn5)",
  "Neurons",
  "Autres"
)

DISPLAY_PALETTE <- c(
  GLOBAL_PALETTE[KEEP_LABELS],            # couleurs d'origine
  "Neurons" = "#AEC7E8",                  # bleu clair unifié
  "Autres"  = "#D9D9D9"                   # gris neutre
)
# Supprimer les éventuels NA (labels absents de GLOBAL_PALETTE)
DISPLAY_PALETTE <- DISPLAY_PALETTE[!is.na(DISPLAY_PALETTE)]

plot_df$display_annotation <- factor(plot_df$display_annotation,
                                      levels = DISPLAY_LEVELS)

cat("\nDistribution des display_annotations:\n")
print(sort(table(plot_df$display_annotation), decreasing = TRUE))

# ------------------------------------------------------------------
# Taille des points selon la population
# ------------------------------------------------------------------
# Adapter selon le nombre total de cellules (toutes coupes mergées)
total_cells <- nrow(plot_df)
pt_base     <- if (total_cells > 150000) 0.35 else if (total_cells > 80000) 0.55 else 0.8
pt_immune   <- pt_base * 1.6   # Immune (Acod1) légèrement plus grands

# ------------------------------------------------------------------
# Construction des panneaux par échantillon
# ------------------------------------------------------------------
# On construit séparément les couches "fond" et "Immune (Acod1)" pour
# garantir que Immune soit toujours au premier plan dans chaque panneau.

panels <- lapply(sample_levels, function(smp) {
  df_smp <- plot_df %>% dplyr::filter(sample == smp)
  if (nrow(df_smp) == 0) return(NULL)

  df_bg     <- df_smp %>% dplyr::filter(display_annotation != "Immune (Acod1)")
  df_immune <- df_smp %>% dplyr::filter(display_annotation == "Immune (Acod1)")

  p <- ggplot() +
    # 1) Arrière-plan : toutes les populations sauf Immune (Acod1)
    geom_point(
      data  = df_bg,
      aes(x = x_plot, y = y_plot, color = display_annotation),
      size  = pt_base, alpha = 0.7, stroke = 0
    ) +
    # 2) Premier plan : Immune (Acod1) en rouge vif
    geom_point(
      data  = df_immune,
      aes(x = x_plot, y = y_plot, color = display_annotation),
      size  = pt_immune, alpha = 0.95, stroke = 0
    ) +
    geom_segment(
      data = scalebar_df %>% dplyr::filter(sample == smp),
      aes(x = x0, xend = x1, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 1.8,
      color = "white"
    ) +
    geom_segment(
      data = scalebar_df %>% dplyr::filter(sample == smp),
      aes(x = x0, xend = x1, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.8,
      color = "black"
    ) +
    geom_text(
      data = scalebar_df %>% dplyr::filter(sample == smp),
      aes(x = (x0 + x1) / 2, y = y_text, label = label),
      inherit.aes = FALSE,
      size = 2.8,
      color = "black",
      vjust = 1
    ) +
    scale_color_manual(
      values = DISPLAY_PALETTE,
      drop   = FALSE,
      name   = "Annotation"
    ) +
    coord_fixed(
      xlim   = global_xlim,
      ylim   = global_ylim,
      expand = FALSE
    ) +
    labs(
      title = SAMPLE_LABELS[smp],
      x = NULL, y = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title       = element_text(face = "bold", size = 10, hjust = 0.5),
      axis.text        = element_text(size = 6),
      axis.ticks       = element_line(linewidth = 0.3),
      axis.line        = element_line(linewidth = 0.3),
      panel.border     = element_rect(color = "grey70", fill = NA, linewidth = 0.3),
      legend.position  = "none"   # légende commune via patchwork
    )

  p
})

panels <- panels[!sapply(panels, is.null)]

# ------------------------------------------------------------------
# Figure combinée 2×2 avec légende commune
# ------------------------------------------------------------------
# On récupère la légende depuis un panneau avec légende visible
p_legend_source <- ggplot(plot_df,
                           aes(x = x_plot, y = y_plot, color = display_annotation)) +
  geom_point(size = 2) +
  scale_color_manual(
    values = DISPLAY_PALETTE,
    drop   = FALSE,
    name   = "Annotation"
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 3.5, alpha = 1),
    ncol = 1
  )) +
  theme_void() +
  theme(
    legend.title     = element_text(face = "bold", size = 9),
    legend.text      = element_text(size = 8),
    legend.key.size  = unit(0.45, "cm"),
    legend.spacing.y = unit(0.15, "cm")
  )

library(grid)
extract_legend <- function(p) {
  g    <- ggplotGrob(p)
  leg  <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
  leg
}
legend_grob <- extract_legend(p_legend_source)

# Assembler les 4 panneaux en une seule ligne (ncol = 4)
# Chaque panneau a coord_fixed() → même hauteur physique garantie par
# plot_layout(heights = 1) ; la largeur de chaque panneau est identique.
p_grid <- wrap_plots(panels, ncol = 4) +
  plot_annotation(
    title    = "Spatial distribution of cell types by sample",
    subtitle = sprintf("BANKSY lam02 res09 | centered per sample | shared window ~%.0f µm | scale bar %d µm | n=%s",
                       target_span, SCALE_BAR_UM, format(total_cells, big.mark = ",")),
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 8, colour = "grey40", hjust = 0.5)
    )
  ) &
  # Forcer chaque panneau à la même hauteur (aspect géré par coord_fixed)
  theme(plot.margin = margin(4, 4, 4, 4))

# Combiner grille + légende à droite
p_final <- wrap_elements(p_grid) +
  wrap_elements(legend_grob) +
  plot_layout(widths = c(8, 1))

# ------------------------------------------------------------------
# Export : 4 panneaux côte à côte → largeur 24 pouces, hauteur 6
# ------------------------------------------------------------------
pdf_out <- file.path(out_dir, "xy_spatial_all_annotations_combined_lam02_res09.pdf")
jpg_out <- file.path(out_dir, "xy_spatial_all_annotations_combined_lam02_res09.jpg")

cairo_pdf(pdf_out, width = 24, height = 6)
print(p_final)
dev.off()
cat("Saved:", pdf_out, "\n")

jpeg(jpg_out, width = 24 * 150, height = 6 * 150, res = 150, quality = 95)
print(p_final)
dev.off()
cat("Saved:", jpg_out, "\n")

cat("\n=== Script 41 terminé avec succès ===\n")
