#!/usr/bin/env Rscript
# ============================================================================
# Script 38 — Spatial figures révisées
#
# SECTION 1 : XY zoom autour de la niche Immune (Acod1) par échantillon LCMV
#             Output : outputs/banksy/spatial_annotations/
#                      xy_zoom_immune_niche_[sample].pdf/.jpg
#
# SECTION 2 : Feature plots Map2 et Ptprc côte à côte (tous échantillons)
#             Output : outputs/banksy/spatial_annotations/
#                      xy_map2_ptprc_all_samples.pdf/.jpg
#
# SECTION 3 : Zoom microglies (P2ry12) proches des T cells / Immune (Acod1)
#             à 1wpi (distance < 50 µm), fenêtre ±500 µm autour du centroïde
#             Output : outputs/banksy/spatial_annotations/
#                      xy_zoom_microglia_close_tcell_1wpi.pdf/.jpg
# ============================================================================

set.seed(1997)

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(FNN)
})

base_path <- normalizePath(".")
setwd(base_path)

source("scripts/00_palette.R")

obj_path  <- file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds")
csv_annot <- "ncells_by_sample_lam02_res09_joint_long.csv"
out_dir   <- file.path("outputs", "banksy", "spatial_annotations")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------------------------------------------------------
# Helper : trouve la colonne de clustering (lam0.2, res0.9)
# ------------------------------------------------------------------
find_cl_col <- function(se, lam, res) {
  cd      <- as.data.frame(SummarizedExperiment::colData(se))
  cn      <- colnames(cd)
  lam_pat <- paste0("lam", gsub("\\.", "\\\\.", as.character(lam)))
  res_pat <- paste0("_res", gsub("\\.", "\\\\.", as.character(res)), "$")
  cols    <- cn[grepl(lam_pat, cn) & grepl(res_pat, cn)]
  if (length(cols) == 0) stop("Clustering column not found for lam=", lam, " res=", res)
  cols[1]
}

# ==================================================================
# Chargement de l'objet et reconstruction des annotations
# ==================================================================
cat("\n=== Chargement de l'objet ===\n")
se <- readRDS(obj_path)
cat("Object loaded:", obj_path, "\n")

cluster_col     <- find_cl_col(se, 0.2, 0.9)
cd              <- as.data.frame(SummarizedExperiment::colData(se))
banksy_clusters <- as.numeric(cd[[cluster_col]])
samples_vec     <- as.character(cd$sample)

spatial_mat <- as.data.frame(SpatialExperiment::spatialCoords(se))
if (!all(c("sdimx", "sdimy") %in% colnames(spatial_mat))) {
  stop("Coordonnées spatiales introuvables : colonnes 'sdimx' et 'sdimy' attendues")
}

# Reconstruction annotation Domain -> label
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
lcmv_samples <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")

# Data frame global
plot_df <- data.frame(
  x          = spatial_mat$sdimx,
  y          = spatial_mat$sdimy,
  sample     = samples_vec,
  annotation = cell_annotations,
  stringsAsFactors = FALSE
)

cat("Cellules totales :", nrow(plot_df), "\n")
cat("Annotations présentes :\n")
print(sort(table(plot_df$annotation)))

# ==================================================================
# SECTION 1 — XY zoom sur la niche Immune (Acod1)
# ==================================================================
cat("\n=== SECTION 1 : XY zoom niche Immune (Acod1) ===\n")

HALF_WINDOW_UM <- 2000  # demi-fenêtre centrée sur le centroïde des cellules immunes

# Ordre et palette pour les annotations présentes
ann_levels_all <- order_annotations(unique(plot_df$annotation), extended = TRUE)

# Couleurs : GLOBAL_PALETTE + gris pour absentes
build_palette <- function(annots) {
  pal  <- GLOBAL_PALETTE[names(GLOBAL_PALETTE) %in% annots]
  miss <- setdiff(annots, names(pal))
  if (length(miss) > 0) {
    extra <- setNames(
      colorRampPalette(c("#CCCCCC", "#888888"))(length(miss)),
      miss
    )
    pal <- c(pal, extra)
  }
  pal
}

for (smp in lcmv_samples) {
  cat("  Traitement échantillon :", smp, "\n")

  df_smp <- plot_df %>% dplyr::filter(sample == smp)

  # Cellules Immune (Acod1) dans cet échantillon
  immune_cells <- df_smp %>% dplyr::filter(annotation == "Immune (Acod1)")

  if (nrow(immune_cells) == 0) {
    cat("    Aucune cellule Immune (Acod1) dans", smp, "— échantillon ignoré.\n")
    next
  }

  # Centroïde des cellules immunes + fenêtre fixe
  cx   <- mean(immune_cells$x)
  cy   <- mean(immune_cells$y)
  xmin <- cx - HALF_WINDOW_UM
  xmax <- cx + HALF_WINDOW_UM
  ymin <- cy - HALF_WINDOW_UM
  ymax <- cy + HALF_WINDOW_UM

  zoom_scale_bar_um <- if (HALF_WINDOW_UM <= 1200) 200 else 500
  sb_xmin <- xmax - zoom_scale_bar_um - 100
  sb_xmax <- xmax - 100
  sb_y    <- ymin + 100

  # Toutes les cellules dans cette zone
  df_zoom <- df_smp %>%
    dplyr::filter(x >= xmin, x <= xmax, y >= ymin, y <= ymax)

  cat("    Cellules dans la zone :", nrow(df_zoom),
      "(dont Immune (Acod1) :", sum(df_zoom$annotation == "Immune (Acod1)"), ")\n")

  # Ordre et palette locaux
  local_anns <- order_annotations(unique(df_zoom$annotation), extended = TRUE)
  local_anns <- local_anns[local_anns %in% unique(df_zoom$annotation)]
  df_zoom$annotation <- factor(df_zoom$annotation, levels = local_anns)

  pal_local <- build_palette(local_anns)

  # Taille des points adaptée à la densité
  pt_size <- if (nrow(df_zoom) > 20000) 0.2 else if (nrow(df_zoom) > 5000) 0.4 else 0.6

  # Séparer les cellules immunes pour les afficher au premier plan
  df_bg     <- df_zoom %>% dplyr::filter(annotation != "Immune (Acod1)")
  df_immune <- df_zoom %>% dplyr::filter(annotation == "Immune (Acod1)")

  p_zoom <- ggplot() +
    # Fond : toutes les autres cellules
    geom_point(data = df_bg,
               aes(x = x, y = y, color = annotation),
               size = pt_size, alpha = 0.6, stroke = 0) +
    # Premier plan : cellules Immune (Acod1) en rouge, légèrement plus grandes
    geom_point(data = df_immune,
               aes(x = x, y = y),
               color = "#D62728", size = pt_size * 1.8, alpha = 0.9, stroke = 0) +
    # Entrée manuelle pour la légende
    geom_point(data = df_immune,
               aes(x = x, y = y, color = annotation),
               size = pt_size * 1.8, alpha = 0.9, stroke = 0) +
    scale_color_manual(values = pal_local, name = "Annotation",
                       drop = FALSE) +
    coord_fixed() +
    annotate("segment",
             x = sb_xmin, xend = sb_xmax,
             y = sb_y, yend = sb_y,
             linewidth = 1.2, color = "black") +
    annotate("text",
             x = (sb_xmin + sb_xmax) / 2,
             y = sb_y + 80,
             label = paste0(zoom_scale_bar_um, " µm"),
             size = 3.5,
             fontface = "bold",
             color = "black") +
    labs(
      title    = paste0("Zoom niche Immune (Acod1) — ", smp),
      subtitle = sprintf("Fenêtre ±%dµm autour du centroïde | %d cellules",
                         HALF_WINDOW_UM, nrow(df_zoom)),
      x = "X (µm)", y = "Y (µm)"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle    = element_text(size = 9, hjust = 0.5, color = "grey40"),
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 8),
      legend.key.size  = unit(0.4, "cm"),
      axis.text        = element_text(size = 8),
      panel.border     = element_rect(color = "grey60", fill = NA, linewidth = 0.4)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

  # Fenêtre carrée => aspect ratio 1:1
  fig_w <- 7
  fig_h <- 7

  slug_smp <- gsub("[^a-z0-9]+", "_", tolower(smp))
  pdf_out  <- file.path(out_dir, paste0("xy_zoom_immune_niche_", slug_smp, ".pdf"))
  jpg_out  <- file.path(out_dir, paste0("xy_zoom_immune_niche_", slug_smp, ".jpg"))

  cairo_pdf(pdf_out, width = fig_w, height = fig_h)
  print(p_zoom)
  dev.off()
  cat("    Saved:", pdf_out, "\n")

  jpeg(jpg_out, width = round(fig_w * 150), height = round(fig_h * 150),
       res = 150, quality = 95)
  print(p_zoom)
  dev.off()
  cat("    Saved:", jpg_out, "\n")
}

# ==================================================================
# SECTION 2 — Feature plots Map2 et Ptprc (tous échantillons)
# ==================================================================
cat("\n=== SECTION 2 : Feature plots Map2 et Ptprc ===\n")

# Sélection de l'assay normalisé
avail_assays <- SummarizedExperiment::assayNames(se)
assay_use <- if ("logcounts" %in% avail_assays) {
  "logcounts"
} else if ("normcounts" %in% avail_assays) {
  "normcounts"
} else {
  avail_assays[1]
}
cat("Assay utilisé :", assay_use, "\n")

expr_mat <- SummarizedExperiment::assay(se, assay_use)

genes_plot <- c("Map2", "Ptprc")
missing_g  <- setdiff(genes_plot, rownames(expr_mat))
if (length(missing_g) > 0) {
  cat("Gènes absents de l'objet :", paste(missing_g, collapse = ", "), "\n")
}
genes_ok <- intersect(genes_plot, rownames(expr_mat))

# Clipping à q99 par gène
q99 <- function(x) quantile(x, 0.99, na.rm = TRUE)

# Shared XY limits across sample facets (same spatial scale within this panel)
gx <- range(spatial_mat$sdimx, na.rm = TRUE)
gy <- range(spatial_mat$sdimy, na.rm = TRUE)
xpad <- max(diff(gx) * 0.03, 30)
ypad <- max(diff(gy) * 0.03, 30)
xlim_all <- c(gx[1] - xpad, gx[2] + xpad)
ylim_all <- c(gy[1] - ypad, gy[2] + ypad)

# Ordre des échantillons pour les facettes
samples_present <- unique(samples_vec)
sample_levels   <- c(
  sample_order[sample_order %in% samples_present],
  setdiff(samples_present, sample_order)
)

feature_plots <- lapply(genes_ok, function(gene) {
  expr   <- as.numeric(expr_mat[gene, ])
  cap    <- q99(expr)
  expr_c <- pmin(expr, cap)

  df_gene <- data.frame(
    x      = spatial_mat$sdimx,
    y      = spatial_mat$sdimy,
    expr   = expr_c,
    sample = factor(samples_vec, levels = sample_levels),
    stringsAsFactors = FALSE
  )

  # Ordonner par expression croissante pour afficher les cellules exprimantes
  # au premier plan
  df_gene <- df_gene[order(df_gene$expr), ]

  pt_size <- if (nrow(df_gene) > 50000) 0.15 else if (nrow(df_gene) > 20000) 0.25 else 0.4

  ggplot(df_gene, aes(x = x, y = y, color = expr)) +
    geom_point(size = pt_size, stroke = 0, alpha = 0.8) +
    scale_color_gradient(
      low  = "grey90",
      high = "#CC0000",
      name = sprintf("%s\n(log, q99\nclip: %.2f)", gene, cap),
      limits = c(0, cap),
      oob   = scales::squish
    ) +
    facet_wrap(~sample, nrow = 1, scales = "fixed") +
    coord_fixed(xlim = xlim_all, ylim = ylim_all, expand = FALSE) +
    labs(
      title = gene,
      x     = "X (µm)",
      y     = "Y (µm)"
    ) +
    theme_classic(base_size = 10) +
    theme(
      aspect.ratio     = 1,
      plot.title       = element_text(face = "bold.italic", size = 13, hjust = 0.5),
      strip.text       = element_text(face = "bold", size = 9),
      strip.background = element_blank(),
      axis.text        = element_text(size = 6),
      legend.title     = element_text(size = 8),
      legend.text      = element_text(size = 7),
      panel.border     = element_rect(color = "grey60", fill = NA, linewidth = 0.4)
    )
})

if (length(feature_plots) == 2) {
  p_combined <- feature_plots[[1]] / feature_plots[[2]] +
    plot_annotation(
      title    = "Spatial expression: Map2 and Ptprc per sample",
      theme    = theme(plot.title = element_text(face = "bold", size = 13,
                                                  hjust = 0.5))
    )
} else if (length(feature_plots) == 1) {
  p_combined <- feature_plots[[1]]
} else {
  stop("Aucun des gènes Map2 / Ptprc n'est présent dans l'objet.")
}

mp_pdf <- file.path(out_dir, "xy_map2_ptprc_all_samples.pdf")
mp_jpg <- file.path(out_dir, "xy_map2_ptprc_all_samples.jpg")

cairo_pdf(mp_pdf, width = 18, height = 10)
print(p_combined)
dev.off()
cat("Saved:", mp_pdf, "\n")

jpeg(mp_jpg, width = 18 * 150, height = 10 * 150, res = 150, quality = 95)
print(p_combined)
dev.off()
cat("Saved:", mp_jpg, "\n")

# ==================================================================
# SECTION 3 — Zoom microglies proches des T cells / Immune à 1wpi
# ==================================================================
cat("\n=== SECTION 3 : Zoom microglies proches des T cells à 1wpi ===\n")

df_1wpi <- plot_df %>% dplyr::filter(sample == "LCMV_1wpi")

if (nrow(df_1wpi) == 0) {
  cat("Aucune cellule LCMV_1wpi trouvée; section 3 ignorée.\n")
} else {
  mg_df <- df_1wpi %>% dplyr::filter(annotation == "Microglia (P2ry12)")
  tc_df <- df_1wpi %>%
    dplyr::filter(grepl("T cells|Immune \\(Acod1\\)", annotation))

  cat("Microglia (P2ry12) à 1wpi:", nrow(mg_df), "\n")
  cat("T cells/Immune (Acod1) à 1wpi:", nrow(tc_df), "\n")

  if (nrow(mg_df) == 0 || nrow(tc_df) == 0) {
    cat("Population manquante pour le calcul de distance; section 3 ignorée.\n")
  } else {
    mg_xy <- as.matrix(mg_df[, c("x", "y")])
    tc_xy <- as.matrix(tc_df[, c("x", "y")])
    knn <- FNN::get.knnx(data = tc_xy, query = mg_xy, k = 1)
    mg_df$dist_to_tcell_um <- as.numeric(knn$nn.dist[, 1])

    mg_close <- mg_df %>% dplyr::filter(dist_to_tcell_um < 50)
    cat("Microglia proches (<50µm):", nrow(mg_close), "\n")

    if (nrow(mg_close) == 0) {
      cat("Aucune microglie proche des T cells (<50µm); section 3 ignorée.\n")
    } else {
      cx <- mean(mg_close$x, na.rm = TRUE)
      cy <- mean(mg_close$y, na.rm = TRUE)
      half_window <- 500

      xlim_zoom <- c(cx - half_window, cx + half_window)
      ylim_zoom <- c(cy - half_window, cy + half_window)

      zoom_scale_bar_um <- if (half_window <= 1200) 200 else 500
      sb_xmin <- xlim_zoom[2] - zoom_scale_bar_um - 100
      sb_xmax <- xlim_zoom[2] - 100
      sb_y    <- ylim_zoom[1] + 100

      df_zoom <- df_1wpi %>%
        dplyr::filter(
          x >= xlim_zoom[1], x <= xlim_zoom[2],
          y >= ylim_zoom[1], y <= ylim_zoom[2]
        )

      ann_levels <- order_annotations(unique(df_zoom$annotation), extended = TRUE)
      ann_levels <- ann_levels[ann_levels %in% unique(df_zoom$annotation)]
      df_zoom$annotation <- factor(df_zoom$annotation, levels = ann_levels)

      pal_zoom <- GLOBAL_PALETTE[names(GLOBAL_PALETTE) %in% ann_levels]
      miss <- setdiff(ann_levels, names(pal_zoom))
      if (length(miss) > 0) {
        pal_zoom <- c(
          pal_zoom,
          setNames(colorRampPalette(c("#CCCCCC", "#888888"))(length(miss)), miss)
        )
      }

      pt_size <- if (nrow(df_zoom) > 20000) 0.2 else if (nrow(df_zoom) > 5000) 0.35 else 0.5

      p_zoom_tcell <- ggplot(df_zoom, aes(x = x, y = y, color = annotation)) +
        geom_point(size = pt_size, alpha = 0.85, stroke = 0) +
        scale_color_manual(values = pal_zoom, drop = FALSE, name = "Annotation") +
        coord_fixed(xlim = xlim_zoom, ylim = ylim_zoom, expand = FALSE) +
        annotate("segment",
                 x = sb_xmin, xend = sb_xmax,
                 y = sb_y, yend = sb_y,
                 linewidth = 1.2, color = "black") +
        annotate("text",
                 x = (sb_xmin + sb_xmax) / 2,
                 y = sb_y + 80,
                 label = paste0(zoom_scale_bar_um, " µm"),
                 size = 3.5,
                 fontface = "bold",
                 color = "black") +
        labs(
          title = "Zoom 1wpi: microglia proches des T cells / Immune (Acod1)",
          subtitle = sprintf("Microglia <50µm de T cells: n=%d | Fenêtre ±%dµm autour du centroïde",
                             nrow(mg_close), half_window),
          x = "X (µm)",
          y = "Y (µm)"
        ) +
        theme_classic(base_size = 11) +
        theme(
          plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.4, "cm"),
          panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.4)
        ) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

      zoom_pdf <- file.path(out_dir, "xy_zoom_microglia_close_tcell_1wpi.pdf")
      zoom_jpg <- file.path(out_dir, "xy_zoom_microglia_close_tcell_1wpi.jpg")

      cairo_pdf(zoom_pdf, width = 7, height = 7)
      print(p_zoom_tcell)
      dev.off()
      cat("Saved:", zoom_pdf, "\n")

      jpeg(zoom_jpg, width = 7 * 150, height = 7 * 150, res = 150, quality = 95)
      print(p_zoom_tcell)
      dev.off()
      cat("Saved:", zoom_jpg, "\n")
    }
  }
}

cat("\n=== Script 38 terminé avec succès ===\n")
