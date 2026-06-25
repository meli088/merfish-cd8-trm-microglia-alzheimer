#!/usr/bin/env Rscript
# ============================================================================
# Script 43 — Assemblage Fig 3 (multi-panneaux PDF)
#
# Fig 3 — Immune niches are surrounded by IFN-responsive cells
#
# Panels :
#   J  fig1_xy_ifn_immune_overlay_OK
#   K  violin_nearest_immune_distance_lam02_res09_OK
#   L  umap_grid_res03  (IFN sub-clustering — ifn/analysis/figures/)
#   M  xy_map2_ptprc_all_samples
#   N  fig_grid_ifn_faceted_shared_scale_OK
#   O  fig_ifng_dotplot_all_celltypes
#   P  fig1_ifng_bubble_OK
#   Q  fig_ifng_tcell_overtime
#
# Output : outputs/fig3_ifn_responsive_niches.pdf
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

out_pdf <- file.path("outputs", "fig3_ifn_responsive_niches.pdf")
if (!dir.exists("outputs")) dir.create("outputs")

# ------------------------------------------------------------------
# Chemins des figures sources
# ------------------------------------------------------------------
root <- file.path("outputs", "banksy")

paths <- list(
  J = file.path(root, "ifn_immune_overlay",
                "fig1_xy_ifn_immune_overlay_OK.jpg"),
  K = file.path(root, "nearest_immune_distance_OK",
                "violin_nearest_immune_distance_lam02_res09_OK.jpg"),
  L = file.path(root, "ifn", "analysis", "figures",
                "umap_grid_res03.jpg"),
  M = file.path(root, "spatial_annotations",
                "xy_map2_ptprc_all_samples.jpg"),
  N = file.path(root, "inflammatory_niche_step2_grid_100um_global",
                "fig_grid_ifn_faceted_shared_scale_OK.jpg"),
  O = file.path(root, "ifn_immune_overlay",
                "fig_ifng_dotplot_all_celltypes.jpg"),
  P = file.path(root, "ifn_immune_overlay",
                "fig1_ifng_bubble_OK.jpg"),
  Q = file.path(root, "ifn_immune_overlay",
                "fig_ifng_tcell_overtime.jpg")
)

# Vérification
missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing) > 0) warning("Fichiers manquants : ", paste(missing, collapse = ", "))
message("Figures trouvées : ", sum(file.exists(unlist(paths))), " / ", length(paths))

# ------------------------------------------------------------------
# Helper : image JPG → ggplot avec tag de panel
# ------------------------------------------------------------------
img_to_gg <- function(path, label = NULL) {
  if (!file.exists(path)) {
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
pJ <- img_to_gg(paths$J, "J")
pK <- img_to_gg(paths$K, "K")
pL <- img_to_gg(paths$L, "L")
pM <- img_to_gg(paths$M, "M")
pN <- img_to_gg(paths$N, "N")
pO <- img_to_gg(paths$O, "O")
pP <- img_to_gg(paths$P, "P")
pQ <- img_to_gg(paths$Q, "Q")

# ------------------------------------------------------------------
# Assemblage en pages thématiques
# ------------------------------------------------------------------
# PAGE 1 : Panel J — XY overlay IFN/Immune
page1 <- pJ +
  plot_annotation(
    title = "Fig 3 — Panel J : XY IFN-Immune overlay",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

# PAGE 2 : Panel K — Violin nearest immune distance
page2 <- pK +
  plot_annotation(
    title = "Fig 3 — Panel K : Nearest immune distance",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

# PAGE 3 : Panel L — UMAP IFN sub-clustering res03
page3 <- pL +
  plot_annotation(
    title = "Fig 3 — Panel L : UMAP IFN sub-clustering (res03)",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

# PAGE 4 : Panel M — Map2 / Ptprc XY (tous échantillons)
page4 <- pM +
  plot_annotation(
    title = "Fig 3 — Panel M : Map2 & Ptprc spatial",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

# PAGE 5 : Panel N — IFN grid faceted
page5 <- pN +
  plot_annotation(
    title = "Fig 3 — Panel N : IFN grid (faceted, shared scale)",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

# PAGE 6 : Panels O–P — Ifng dotplot + bubble plot
page6 <- (pO | pP) +
  plot_annotation(
    title = "Fig 3 — Panels O–P",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

# PAGE 7 : Panel Q — Ifng T cells overtime
page7 <- pQ +
  plot_annotation(
    title = "Fig 3 — Panel Q : Ifng+ T cells overtime",
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
dev.off()

message("Saved: ", out_pdf)
message("\n=== Script 43 terminé avec succès ===\n")
