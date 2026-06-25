#!/usr/bin/env Rscript
# ============================================================================
# Script 37 — Dotplot canonical markers (revised) +
#             Composition Immune vs Non-Immune overtime
#
# SECTION 1 : Dotplot avec marqueurs canoniques révisés incluant marqueurs
#             macrophages (Ms4a7, Tgfbi, Clec12a, Sall1, Sparc)
#             Output : outputs/banksy/spatial_annotations/
#                      dotplot_canonical_markers_revised.pdf/.jpg
#
# SECTION 2 : Stacked barplot Immune vs Non-Immune, 2 barres par condition,
#             chaque barre colorée par sous-population
#             Output : outputs/banksy/umap_annotated/
#                      composition_immune_vs_nonimmune_overtime.pdf/.jpg
# ============================================================================

set.seed(1997)

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(Cairo)
})

base_path <- normalizePath(".")
setwd(base_path)

source("scripts/00_palette.R")

RENAME_ANNOT <- c("Prolif neural/glial (Ccdc153)" = "Ependymal (Ccdc153)")
GLOBAL_PALETTE_LOCAL <- c(
  GLOBAL_PALETTE,
  "Ependymal (Ccdc153)" = GLOBAL_PALETTE[["Prolif neural/glial (Ccdc153)"]]
)

obj_path  <- file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds")
csv_annot <- "ncells_by_sample_lam02_res09_joint_long.csv"
out_dir1  <- file.path("outputs", "banksy", "spatial_annotations")
out_dir2  <- file.path("outputs", "banksy", "umap_annotated")

for (d in c(out_dir1, out_dir2)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ------------------------------------------------------------------
# Helper : trouve la colonne de clustering (lam0.2, res0.9)
# ------------------------------------------------------------------
find_cl_col <- function(se, lam, res) {
  cd  <- as.data.frame(SummarizedExperiment::colData(se))
  cn  <- colnames(cd)
  lam_pat <- paste0("lam", gsub("\\.", "\\\\.", as.character(lam)))
  res_pat <- paste0("_res", gsub("\\.", "\\\\.", as.character(res)), "$")
  cols <- cn[grepl(lam_pat, cn) & grepl(res_pat, cn)]
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
samples         <- as.character(cd$sample)

# Reconstruction de la correspondance Domain -> annotation depuis le CSV
anno_data <- read.delim(csv_annot, sep = ";", stringsAsFactors = FALSE)
anno_map  <- anno_data %>%
  filter(annotation != "" & !is.na(annotation)) %>%
  dplyr::select(banksy_domain, annotation) %>%
  dplyr::distinct()

anno_map$cluster_id <- as.numeric(gsub("Domain_", "", anno_map$banksy_domain))
anno_lookup         <- setNames(trimws(anno_map$annotation), anno_map$cluster_id)

cell_annotations <- anno_lookup[as.character(banksy_clusters)]
cell_annotations[is.na(cell_annotations)] <- "Non annote"
cell_annotations <- dplyr::recode(cell_annotations, !!!RENAME_ANNOT)

cat("Cellules chargées :", length(cell_annotations), "\n")
cat("Annotations présentes :\n")
print(sort(table(cell_annotations)))

# ==================================================================
# SECTION 1 — Dotplot marqueurs canoniques révisé
# ==================================================================
cat("\n=== SECTION 1 : Dotplot canonical markers (revised) ===\n")

markers <- c(
  "Ptprc", "Map2", "Aqp4", "Plp1", "Sox10",
  "P2ry12", "Tmem119", "Cldn5", "Rbfox3"
)

# Sélection du meilleur assay normalisé disponible
avail_assays <- SummarizedExperiment::assayNames(se)
assay_use <- if ("logcounts" %in% avail_assays) {
  "logcounts"
} else if ("normcounts" %in% avail_assays) {
  "normcounts"
} else {
  avail_assays[1]
}
cat("Assay utilisé :", assay_use, "\n")

expr_mat    <- SummarizedExperiment::assay(se, assay_use)
markers_use <- markers[markers %in% rownames(expr_mat)]
missing_m   <- setdiff(markers, markers_use)
if (length(missing_m) > 0) {
  cat("Marqueurs absents de l'objet :", paste(missing_m, collapse = ", "), "\n")
}
cat("Marqueurs trouvés :", length(markers_use), "/", length(markers), "\n")

# Calcul : expression moyenne + fraction de cellules exprimant > 0
# pour chaque (gène × annotation)
dot_data <- do.call(rbind, lapply(markers_use, function(gene) {
  expr <- as.numeric(expr_mat[gene, ])
  do.call(rbind, lapply(unique(cell_annotations), function(ann) {
    idx <- cell_annotations == ann
    e   <- expr[idx]
    data.frame(
      gene       = gene,
      annotation = ann,
      mean_expr  = mean(e, na.rm = TRUE),
      pct_expr   = mean(e > 0, na.rm = TRUE) * 100,
      stringsAsFactors = FALSE
    )
  }))
}))

# Ordres des axes
dot_data$gene <- factor(dot_data$gene, levels = markers_use)
ann_levels    <- order_annotations(unique(dot_data$annotation), extended = TRUE)
ann_levels    <- ann_levels[ann_levels %in% unique(dot_data$annotation)]
dot_data$annotation <- factor(dot_data$annotation, levels = ann_levels)

# Plafonnement de l'expression au 99e percentile
cap_val <- quantile(dot_data$mean_expr, 0.99, na.rm = TRUE)
dot_data$mean_expr_capped <- pmin(dot_data$mean_expr, cap_val)
cat(sprintf("Échelle couleur plafonnée au 99e percentile : %.3f\n", cap_val))

p_dot <- ggplot(dot_data, aes(x = gene, y = annotation)) +
  geom_point(aes(size = pct_expr, color = mean_expr_capped)) +
  scale_color_gradient(
    low  = "#ffee01",
    high = "#ff001e",
    name = "Mean expr.\n(capped 99p)"
  ) +
  scale_size_continuous(range = c(0.5, 6), name = "% cells\nexpressing") +
  labs(
    title = "Canonical cell-type markers by global annotation",
    x     = NULL,
    y     = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y      = element_text(size = 9),
    plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
  )

n_ann  <- length(levels(dot_data$annotation))
n_gene <- length(markers_use)
fig_w  <- max(9, n_gene * 0.45 + 3)
fig_h  <- max(5, n_ann  * 0.40 + 2)

dot_pdf <- file.path(out_dir1, "dotplot_canonical_markers_revised_v2.pdf")
dot_jpg <- file.path(out_dir1, "dotplot_canonical_markers_revised_v2.jpg")

CairoPDF(dot_pdf, width = fig_w, height = fig_h)
print(p_dot)
dev.off()
cat("Saved:", dot_pdf, "\n")

CairoJPEG(dot_jpg, width = round(fig_w * 150), height = round(fig_h * 150),
    res = 150, quality = 95)
print(p_dot)
dev.off()
cat("Saved:", dot_jpg, "\n")

# ==================================================================
# SECTION 2 — Composition Immune vs Non-Immune overtime
# ==================================================================
cat("\n=== SECTION 2 : Composition Immune vs Non-Immune overtime ===\n")

immune_annots <- c("Immune (Acod1)", "IFN responsive (Ifit1)", "Microglia (P2ry12)")

sample_order  <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi", "mock_6wpi")
samples_present <- unique(samples)
sample_levels <- c(
  sample_order[sample_order %in% samples_present],
  setdiff(samples_present, sample_order)
)

comp_df <- data.frame(
  sample     = factor(samples,          levels = sample_levels),
  annotation = trimws(cell_annotations),
  stringsAsFactors = FALSE
)
comp_df$category <- ifelse(
  comp_df$annotation %in% immune_annots,
  "Immune",
  "Non-Immune"
)
comp_df$category <- factor(comp_df$category, levels = c("Immune", "Non-Immune"))

# Proportion de chaque sous-population au sein de chaque (sample × catégorie)
prop_df <- comp_df %>%
  dplyr::count(sample, category, annotation) %>%
  dplyr::group_by(sample, category) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()

# Ordre des annotations selon ANNOTATION_ORDER
ann_ord <- order_annotations(unique(prop_df$annotation), extended = TRUE)
ann_ord <- ann_ord[ann_ord %in% unique(prop_df$annotation)]
prop_df$annotation <- factor(prop_df$annotation, levels = ann_ord)

# Palette : GLOBAL_PALETTE + couleurs grises pour les absentes
pal_use     <- GLOBAL_PALETTE_LOCAL[names(GLOBAL_PALETTE_LOCAL) %in% levels(prop_df$annotation)]
missing_col <- setdiff(levels(prop_df$annotation), names(pal_use))
if (length(missing_col) > 0) {
  extra_cols <- setNames(
    colorRampPalette(c("#CCCCCC", "#777777"))(length(missing_col)),
    missing_col
  )
  pal_use <- c(pal_use, extra_cols)
}

p_comp <- ggplot(prop_df, aes(x = sample, y = prop, fill = annotation)) +
  geom_bar(stat = "identity", position = "stack", width = 0.75) +
  facet_wrap(~category, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = pal_use, name = "Annotation") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  labs(
    title = "Cell-type composition: Immune vs Non-Immune over time",
    x     = "Condition",
    y     = "Proportion within category"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y      = element_text(size = 9),
    plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.text       = element_text(face = "bold", size = 11),
    strip.background = element_blank(),
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    legend.key.size  = unit(0.4, "cm")
  )

comp_pdf <- file.path(out_dir2, "composition_immune_vs_nonimmune_overtime.pdf")
comp_jpg <- file.path(out_dir2, "composition_immune_vs_nonimmune_overtime.jpg")

cairo_pdf(comp_pdf, width = 13, height = 7)
print(p_comp)
dev.off()
cat("Saved:", comp_pdf, "\n")

jpeg(comp_jpg, width = 13 * 150, height = 7 * 150, res = 150, quality = 95)
print(p_comp)
dev.off()
cat("Saved:", comp_jpg, "\n")

# ==================================================================
# SECTION 3 — Figure supplémentaire : XY spatial LCMV_1wpi + barplot
#             composition overtime (grandes populations)
# ==================================================================
cat("\n=== SECTION 3 : Figure supplémentaire composition + XY ===\n")

library(SpatialExperiment)
library(patchwork)

# --- Coordonnées spatiales (depuis l'objet déjà chargé) ---
spatial_mat3 <- as.data.frame(SpatialExperiment::spatialCoords(se))
colnames(spatial_mat3) <- c("x", "y")

base_df3 <- data.frame(
  x          = spatial_mat3$x,
  y          = spatial_mat3$y,
  sample     = samples,
  annotation = trimws(cell_annotations),
  stringsAsFactors = FALSE
)

# --- Regroupement en grandes populations (comme script 41) ---
NEURON_PAT3 <- c("Neurons \\(", "Neurons\\(",
                  "Excitatory neurons", "Inhibitory neurons")
KEEP3 <- c("Immune (Acod1)", "IFN responsive (Ifit1)", "Microglia (P2ry12)",
            "Astrocytes (Fgfr3)", "Astrocytes (Gfap)", "Oligodendrocytes (Plp1)",
            "OPC (Pdgfra)", "Vascular (Cldn5)")

base_df3$group <- dplyr::case_when(
  base_df3$annotation %in% KEEP3                                        ~ base_df3$annotation,
  grepl(paste(NEURON_PAT3, collapse = "|"), base_df3$annotation)         ~ "Neurons",
  TRUE                                                                   ~ "Autres"
)

GROUP_LEVELS3 <- c(KEEP3, "Neurons", "Autres")
GROUP_PAL3    <- c(
  GLOBAL_PALETTE[KEEP3],
  "Neurons" = "#AEC7E8",
  "Autres"  = "#D9D9D9"
)
GROUP_PAL3 <- GROUP_PAL3[!is.na(GROUP_PAL3)]
base_df3$group <- factor(base_df3$group, levels = GROUP_LEVELS3)

# ---- Partie gauche : XY tous les 4 échantillons (facetté) ----
SAMPLE_LABELS3 <- c(
  mock_6wpi = "Mock 6wpi",
  LCMV_1wpi = "LCMV 1wpi",
  LCMV_3wpi = "LCMV 3wpi",
  LCMV_6wpi = "LCMV 6wpi"
)

df_all4 <- base_df3 %>%
  dplyr::filter(sample %in% sample_order) %>%
  dplyr::mutate(
    sample_label = factor(SAMPLE_LABELS3[sample],
                          levels = SAMPLE_LABELS3[sample_order])
  )

n_all4  <- nrow(df_all4)
pt3     <- if (n_all4 > 150000) 0.15 else if (n_all4 > 80000) 0.25 else 0.4

df_all4_bg  <- df_all4 %>% dplyr::filter(group != "Immune (Acod1)")
df_all4_imm <- df_all4 %>% dplyr::filter(group == "Immune (Acod1)")

p_xy <- ggplot() +
  geom_point(data = df_all4_bg,
             aes(x = x, y = y, color = group),
             size = pt3, alpha = 0.7, stroke = 0) +
  geom_point(data = df_all4_imm,
             aes(x = x, y = y, color = group),
             size = pt3 * 1.6, alpha = 0.95, stroke = 0) +
  scale_color_manual(values = GROUP_PAL3, drop = FALSE, name = "Annotation") +
  facet_wrap(~sample_label, nrow = 1, scales = "free") +
  labs(
    title = "Spatial distribution — all samples",
    x = "X (µm)", y = "Y (µm)"
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 3, alpha = 1), ncol = 1
  )) +
  theme_classic(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold", size = 11, hjust = 0.5),
    strip.text       = element_text(face = "bold", size = 9),
    strip.background = element_blank(),
    axis.text        = element_text(size = 6),
    aspect.ratio     = 1,
    legend.title     = element_text(face = "bold", size = 9),
    legend.text      = element_text(size = 8),
    legend.key.size  = unit(0.4, "cm"),
    panel.border     = element_rect(color = "grey60", fill = NA, linewidth = 0.4)
  )

# ---- Partie droite : barplot composition par échantillon ----
comp3 <- base_df3 %>%
  dplyr::count(sample, group) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(sample %in% sample_order)

comp3$sample_label <- factor(
  SAMPLE_LABELS3[comp3$sample],
  levels = SAMPLE_LABELS3[sample_order[sample_order %in% comp3$sample]]
)
comp3$group <- factor(comp3$group, levels = GROUP_LEVELS3)

p_bar <- ggplot(comp3, aes(x = sample_label, y = prop, fill = group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = GROUP_PAL3, drop = FALSE, guide = "none") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  labs(
    title = "Cell type composition over time",
    x = NULL, y = "Proportion"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
    axis.text.x      = element_text(angle = 35, hjust = 1, size = 9),
    axis.text.y      = element_text(size = 9),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
  )

# ---- Assemblage patchwork : XY facetté (3) + barplot (1) ----
p_supp <- p_xy + p_bar +
  plot_layout(widths = c(3, 1)) +
  plot_annotation(
    title = "Spatial distribution and cell-type composition — LCMV infection",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5)
    )
  )

supp_pdf <- file.path(out_dir2, "fig_supp_composition_xy_overtime.pdf")
supp_jpg <- file.path(out_dir2, "fig_supp_composition_xy_overtime.jpg")

cairo_pdf(supp_pdf, width = 18, height = 5)
print(p_supp)
dev.off()
cat("Saved:", supp_pdf, "\n")

jpeg(supp_jpg, width = 18 * 150, height = 5 * 150, res = 150, quality = 95)
print(p_supp)
dev.off()
cat("Saved:", supp_jpg, "\n")

cat("\n=== Script 37 terminé avec succès ===\n")
