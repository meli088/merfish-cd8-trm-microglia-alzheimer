#!/usr/bin/env Rscript
# ============================================================================
# Script 42 — Assemblage Fig 2 (multi-panneaux PDF)
#
# Fig 2 — Early life infection induces the emergence and persistence
#          of inflammatory niches
#
# Panels :
#   A  umap_grid_res03
#   B  heatmap_annotated_lam02_res09_top5
#   C  dotplot_canonical_markers_res03
#   D  composition_immune_vs_nonimmune_overtime
#   E  composition_stacked_res03
#   F  xy_zoom_immune_niche_lcmv_1wpi
#   G  fig_volcano_in_vs_out_all_timepoints_merged
#   H  fig_module_score_merged_violin_upregulated_dam
#   I  fig_volcano_in_vs_out_lcmv_1wpi + fig_volcano_in_vs_out_lcmv_6wpi
#   J  fig_module_score_merged_violin_downregulated_dam
#   K  fig_distance_microglia_to_tcell_overtime + fig_dam_score_by_distance_bins
#
# Output : outputs/fig2_inflammatory_niches.pdf
# ============================================================================

set.seed(1997)

suppressPackageStartupMessages({
  library(magick)
  library(grid)
  library(ggplot2)
  library(patchwork)
  library(Cairo)
})

base_path <- normalizePath(".")
setwd(base_path)

out_pdf <- file.path("outputs", "fig2_inflammatory_niches.pdf")
if (!dir.exists("outputs")) dir.create("outputs")

# ------------------------------------------------------------------
# Chemins des figures sources
# ------------------------------------------------------------------
root <- file.path("outputs", "banksy")

paths <- list(
  A  = file.path(root, "immune_acod1", "analysis", "figures",
                 "umap_grid_res03_OK.jpg"),
  B  = file.path(root, "umap_annotated_OK",
                 "heatmap_annotated_lam02_res09_top5_NOT_USED.jpg"),
  C  = file.path(root, "immune_acod1", "analysis", "figures",
                 "dotplot_canonical_markers_res03_OK.jpg"),
  D  = file.path(root, "umap_annotated",
                 "composition_immune_vs_nonimmune_overtime.jpg"),
  E  = file.path(root, "ifn", "analysis", "figures",
                 "composition_stacked_res03.jpg"),
  F1 = file.path(root, "spatial_annotations",
                 "xy_zoom_immune_niche_lcmv_1wpi.jpg"),
  F2 = file.path(root, "spatial_annotations",
                 "xy_zoom_immune_niche_lcmv_3wpi.jpg"),
  F3 = file.path(root, "spatial_annotations",
                 "xy_zoom_immune_niche_lcmv_6wpi.jpg"),
  G  = file.path(root, "microglia_dam_niche",
                 "fig_volcano_in_vs_out_all_timepoints_merged.jpg"),
  H  = file.path(root, "microglia_dam_niche",
                 "fig_module_score_merged_violin_upregulated_dam.jpg"),
  I1 = file.path(root, "microglia_dam_niche",
                 "fig_volcano_in_vs_out_lcmv_1wpi.jpg"),
  I2 = file.path(root, "microglia_dam_niche",
                 "fig_volcano_in_vs_out_lcmv_6wpi.jpg"),
  J  = file.path(root, "microglia_dam_niche",
                 "fig_module_score_merged_violin_downregulated_dam.jpg"),
  K1 = file.path(root, "microglia_dam_niche",
                 "fig_distance_microglia_to_tcell_overtime.jpg"),
  K2 = file.path(root, "microglia_dam_niche",
                 "fig_dam_score_by_distance_bins.jpg"),
  # PAGE 8 : DEG barplots par domaine (T cells / Mac / Microglia)
  L1 = file.path(root, "immune_acod1", "domain_DEG_evolution", "res03", "heatmaps",
                 "barplot_domain_1_t_cells__gzmb__vs_1wpi.jpg"),
  L2 = file.path(root, "immune_acod1", "domain_DEG_evolution", "res03", "heatmaps",
                 "barplot_domain_1_t_cells__gzmb__vs_3wpi.jpg"),
  L3 = file.path(root, "immune_acod1", "domain_DEG_evolution", "res03", "heatmaps",
                 "barplot_domain_3_mac__ctss__vs_1wpi.jpg"),
  L4 = file.path(root, "immune_acod1", "domain_DEG_evolution", "res03", "heatmaps",
                 "barplot_domain_3_mac__ctss__vs_3wpi.jpg"),
  L5 = file.path(root, "immune_acod1", "domain_DEG_evolution", "res03", "heatmaps",
                 "barplot_domain_8_microglia__c1qa__vs_1wpi.jpg"),
  L6 = file.path(root, "immune_acod1", "domain_DEG_evolution", "res03", "heatmaps",
                 "barplot_domain_8_microglia__c1qa__vs_3wpi.jpg"),
  # PAGE 9 : Volcanos par type cellulaire (6wpi vs 1wpi)
  M1 = file.path(root, "immune_niche_volcano_by_celltype",
                 "fig_volcano_t_cells_gzmb__lcmv_6wpi_vs_lcmv_1wpi_OK.jpg"),
  M2 = file.path(root, "immune_niche_volcano_by_celltype",
                 "fig_volcano_mac_ctss__lcmv_6wpi_vs_lcmv_1wpi_OK.jpg"),
  M3 = file.path(root, "immune_niche_volcano_by_celltype",
                 "fig_volcano_microglia_c1qa__lcmv_6wpi_vs_lcmv_1wpi.jpg"),
  # PAGE 10 : Volcanos P2ry12 microglie conditions
  N1 = file.path(root, "microglia_conditions_deg",
                 "fig_volcano_microglia_lcmv_1wpi_vs_mock_6wpi_OK.jpg"),
  N2 = file.path(root, "microglia_conditions_deg",
                 "fig_volcano_microglia_lcmv_6wpi_vs_mock_6wpi_OK.jpg"),
  N3 = file.path(root, "microglia_conditions_deg",
                 "fig_volcano_microglia_lcmv_6wpi_vs_lcmv_1wpi_OK.jpg"),
  # PAGE 11 : Volcanos in vs out niche par timepoint
  O1 = file.path(root, "microglia_dam_niche",
                 "fig_volcano_in_vs_out_lcmv_1wpi.jpg"),
  O2 = file.path(root, "microglia_dam_niche",
                 "fig_volcano_in_vs_out_lcmv_3wpi.jpg"),
  O3 = file.path(root, "microglia_dam_niche",
                 "fig_volcano_in_vs_out_lcmv_6wpi.jpg")
)

# Vérification
missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing) > 0) {
  warning("Fichiers manquants : ", paste(missing, collapse = ", "))
}
message("Figures trouvées : ", sum(file.exists(unlist(paths))), " / ", length(paths))

# ------------------------------------------------------------------
# Helper : lire une image JPG → grob ggplot (rasterGrob)
# ------------------------------------------------------------------
img_to_gg <- function(path, label = NULL) {
  if (!file.exists(path)) {
    # Placeholder gris si fichier absent
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste0("MISSING\n", basename(path)),
                 size = 4, color = "red", hjust = 0.5) +
        theme_void() +
        theme(panel.background = element_rect(fill = "grey95", color = "grey70"))
    )
  }
  img <- magick::image_read(path)
  rg  <- grid::rasterGrob(as.raster(img), interpolate = TRUE)
  p   <- ggplot() +
    annotation_custom(rg, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    theme_void()
  if (!is.null(label)) {
    p <- p + labs(tag = label) +
      theme(
        plot.tag          = element_text(face = "bold", size = 14,
                                         hjust = 0, vjust = 1),
        plot.tag.position = c(0.01, 0.99)
      )
  }
  p
}

# ------------------------------------------------------------------
# Charger toutes les images
# ------------------------------------------------------------------
pA  <- img_to_gg(paths$A,  "A")
pB  <- img_to_gg(paths$B,  "B")
pC  <- img_to_gg(paths$C,  "C")
pD  <- img_to_gg(paths$D,  "D")
pE  <- img_to_gg(paths$E,  "E")
pF1 <- img_to_gg(paths$F1, "F")
pF2 <- img_to_gg(paths$F2)
pF3 <- img_to_gg(paths$F3)
pG  <- img_to_gg(paths$G,  "G")
pH  <- img_to_gg(paths$H,  "H")
pI1 <- img_to_gg(paths$I1, "I")
pI2 <- img_to_gg(paths$I2)
pJ  <- img_to_gg(paths$J,  "J")
pK1 <- img_to_gg(paths$K1, "K")
pK2 <- img_to_gg(paths$K2)

# PAGE 8 — DEG barplots
pL1 <- img_to_gg(paths$L1, "L")
pL2 <- img_to_gg(paths$L2)
pL3 <- img_to_gg(paths$L3)
pL4 <- img_to_gg(paths$L4)
pL5 <- img_to_gg(paths$L5)
pL6 <- img_to_gg(paths$L6)

# PAGE 9 — Volcanos cell types
pM1 <- img_to_gg(paths$M1, "M")
pM2 <- img_to_gg(paths$M2)
pM3 <- img_to_gg(paths$M3)

# PAGE 10 — Volcanos P2ry12 conditions
pN1 <- img_to_gg(paths$N1, "N")
pN2 <- img_to_gg(paths$N2)
pN3 <- img_to_gg(paths$N3)

# PAGE 11 — Volcanos in vs out
pO1 <- img_to_gg(paths$O1, "O")
pO2 <- img_to_gg(paths$O2)
pO3 <- img_to_gg(paths$O3)

# ------------------------------------------------------------------
# Assemblage en pages thématiques
# ------------------------------------------------------------------
# PAGE 1 : Panels A–C (UMAP / heatmap / dotplot)
page1 <- (pA | pB | pC) +
  plot_annotation(
    title = "Fig 2 — Panels A–C",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

# PAGE 2 : Panels D–E (composition)
page2 <- (pD / pE) +
  plot_annotation(
    title = "Fig 2 — Panels D–E",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

# PAGE 3 : Panel F (zooms spatiaux 1wpi / 3wpi / 6wpi)
page3 <- (pF1 | pF2 | pF3) +
  plot_annotation(
    title = "Fig 2 — Panel F",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

# PAGE 4 : Panels G–H (volcano mergé + module score UP)
page4 <- (pG | pH) +
  plot_annotation(
    title = "Fig 2 — Panels G–H",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

# PAGE 5 : Panel I (volcanos 1wpi + 6wpi)
page5 <- (pI1 | pI2) +
  plot_annotation(
    title = "Fig 2 — Panel I",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

# PAGE 6 : Panel J (module score DOWN)
page6 <- pJ +
  plot_annotation(
    title = "Fig 2 — Panel J",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

# PAGE 7 : Panel K (distance + DAM bins)
page7 <- (pK1 | pK2) +
  plot_annotation(
    title = "Fig 2 — Panel K",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

# PAGE 8 : DEG barplots (T cells / Mac / Microglia — 3 lignes × 2 comparaisons)
page8 <- ((pL1 | pL2) / (pL3 | pL4) / (pL5 | pL6)) +
  plot_annotation(
    title = "Fig 2 — Panel L : DEG évolution par domaine",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

# PAGE 9 : Volcanos par type cellulaire (6wpi vs 1wpi)
page9 <- (pM1 | pM2 | pM3) +
  plot_annotation(
    title = "Fig 2 — Panel M : Volcanos T cells / Mac / Microglia C1qa (6wpi vs 1wpi)",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

# PAGE 10 : Volcanos P2ry12 conditions
page10 <- (pN1 | pN2 | pN3) +
  plot_annotation(
    title = "Fig 2 — Panel N : Volcanos P2ry12 (1wpi vs mock | 6wpi vs mock | 6wpi vs 1wpi)",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

# PAGE 11 : Volcanos in vs out niche
page11 <- (pO1 | pO2 | pO3) +
  plot_annotation(
    title = "Fig 2 — Panel O : Volcanos in vs out niche (1wpi | 3wpi | 6wpi)",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

# ------------------------------------------------------------------
# Export multi-page PDF
# ------------------------------------------------------------------
message("\nExport PDF : ", out_pdf)

cairo_pdf(out_pdf, width = 16, height = 9, onefile = TRUE)
print(page1)
print(page2)
print(page3)
print(page4)
print(page5)
print(page6)
print(page7)
print(page8)
print(page9)
print(page10)
print(page11)
dev.off()

message("Saved: ", out_pdf)
message("\n=== Script 42 terminé avec succès ===\n")
