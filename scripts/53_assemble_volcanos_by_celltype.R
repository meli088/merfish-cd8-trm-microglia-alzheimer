#!/usr/bin/env Rscript
# =============================================================
# Script: 53_assemble_volcanos_by_celltype.R
# Goal:   PDF de tous les volcanos DEG par cell type (6wpi vs 1wpi)
# Output: outputs/banksy/panels/supp/fig_supp_volcanos_by_celltype.pdf
# =============================================================

suppressPackageStartupMessages({ library(magick) })
setwd(normalizePath("."))

VOL_DIR  <- file.path("outputs", "banksy", "immune_niche_volcano_by_celltype")
OUT_DIR  <- file.path("outputs", "banksy", "panels", "supp")
OUT_FILE <- file.path(OUT_DIR, "fig_supp_volcanos_by_celltype.pdf")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

PAGE_W_PX   <- 2480L
PAGE_H_PX   <- 1754L
HEADER_H_PX <- 80L
FONT_LBL    <- 24L
FONT_SRC    <- 16L

# Toutes les figures volcano 6wpi vs 1wpi, dans l'ordre biologique
CONTRAST <- "lcmv_6wpi_vs_lcmv_1wpi"

PAIRS <- list(
  list(
    label = "DEG Summary — nombre de DEGs par type cellulaire",
    left  = list(path = file.path(VOL_DIR, "fig_ndeg_summary_by_celltype.jpg"),
                 name = "fig_ndeg_summary_by_celltype.jpg"),
    right = NULL
  ),
  list(
    label = "T cells (Gzmb)  |  T CD4 (Foxp3) — LCMV 6wpi vs 1wpi",
    left  = list(path = file.path(VOL_DIR, "fig_volcano_t_cells_gzmb__lcmv_6wpi_vs_lcmv_1wpi.jpg"),
                 name = "T cells (Gzmb)"),
    right = list(path = file.path(VOL_DIR, "fig_volcano_t_cd4_foxp3__lcmv_6wpi_vs_lcmv_1wpi.jpg"),
                 name = "T CD4 (Foxp3)")
  ),
  list(
    label = "Microglia (C1qa)  |  Mac (Ctss) — LCMV 6wpi vs 1wpi",
    left  = list(path = file.path(VOL_DIR, "fig_volcano_microglia_c1qa__lcmv_6wpi_vs_lcmv_1wpi.jpg"),
                 name = "Microglia (C1qa)"),
    right = list(path = file.path(VOL_DIR, "fig_volcano_mac_ctss__lcmv_6wpi_vs_lcmv_1wpi.jpg"),
                 name = "Mac (Ctss)")
  ),
  list(
    label = "Glials (Gja1)  |  Vascular — LCMV 6wpi vs 1wpi",
    left  = list(path = file.path(VOL_DIR, "fig_volcano_glials_gja1__lcmv_6wpi_vs_lcmv_1wpi.jpg"),
                 name = "Glials (Gja1)"),
    right = list(path = file.path(VOL_DIR, "fig_volcano_vascular_lcmv_6wpi_vs_lcmv_1wpi.jpg"),
                 name = "Vascular")
  ),
  list(
    label = "Vascular (Igfbp2)  |  Oligo — LCMV 6wpi vs 1wpi",
    left  = list(path = file.path(VOL_DIR, "fig_volcano_vascular_igfbp2__lcmv_6wpi_vs_lcmv_1wpi.jpg"),
                 name = "Vascular (Igfbp2)"),
    right = list(path = file.path(VOL_DIR, "fig_volcano_oligo_lcmv_6wpi_vs_lcmv_1wpi.jpg"),
                 name = "Oligo")
  ),
  list(
    label = "Excitatory neurons (Satb2)  |  Inhibitory neurons (Htr2c) — LCMV 6wpi vs 1wpi",
    left  = list(path = file.path(VOL_DIR, "fig_volcano_excitatory_neurons_satb2__lcmv_6wpi_vs_lcmv_1wpi.jpg"),
                 name = "Excitatory neurons (Satb2)"),
    right = list(path = file.path(VOL_DIR, "fig_volcano_inhibitory_neurons_htr2c__lcmv_6wpi_vs_lcmv_1wpi.jpg"),
                 name = "Inhibitory neurons (Htr2c)")
  ),
  list(
    label = "Neurons (Fam107a) — LCMV 6wpi vs 1wpi",
    left  = list(path = file.path(VOL_DIR, "fig_volcano_neurons_fam107a__lcmv_6wpi_vs_lcmv_1wpi.jpg"),
                 name = "Neurons (Fam107a)"),
    right = NULL
  )
)

# --- helpers --------------------------------------------------
make_header <- function(label) {
  image_blank(PAGE_W_PX, HEADER_H_PX, color = "#F2F2F2") %>%
    image_annotate(label,
                   gravity = "NorthWest", location = "+20+18",
                   size = FONT_LBL, font = "sans",
                   weight = 700, color = "#1A1A1A")
}

load_fit <- function(path, w_px, h_px) {
  if (!file.exists(path)) {
    warning("MISSING: ", path)
    return(
      image_blank(w_px, h_px, color = "white") %>%
        image_annotate(paste0("Missing:\n", basename(path)),
                       gravity = "Center", size = 26, color = "red3")
    )
  }
  img <- image_read(path)
  img <- image_resize(img, paste0(w_px, "x", h_px, ">"))
  image_extent(img, paste0(w_px, "x", h_px), gravity = "Center", color = "white")
}

# --- assemblage -----------------------------------------------
message("\n=== Assemblage : ", OUT_FILE)
img_h_px <- PAGE_H_PX - HEADER_H_PX
pages <- list()

for (i in seq_along(PAIRS)) {
  pg <- PAIRS[[i]]
  message(sprintf("  Page %d/%d : %s", i, length(PAIRS), pg$label))

  header <- make_header(pg$label)

  if (is.null(pg$right)) {
    # Page centrée (1 seule figure)
    body <- load_fit(pg$left$path, PAGE_W_PX, img_h_px)
  } else {
    half_w    <- (PAGE_W_PX - 4L) %/% 2L
    left_img  <- load_fit(pg$left$path,  half_w, img_h_px)
    right_img <- load_fit(pg$right$path, half_w, img_h_px)
    sep       <- image_blank(4L, img_h_px, color = "#CCCCCC")
    body      <- image_append(c(left_img, sep, right_img), stack = FALSE)
  }
  pages[[i]] <- image_append(c(header, body), stack = TRUE)
}

all_pages <- do.call(c, pages)
image_write(all_pages, path = OUT_FILE, format = "pdf", density = 200, quality = 85)
message("\n=== Sauvegarde : ", normalizePath(OUT_FILE))
message("    Pages : ", length(PAIRS))
