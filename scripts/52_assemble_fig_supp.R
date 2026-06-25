#!/usr/bin/env Rscript
# =============================================================
# Script: 52_assemble_fig_supp.R
# Goal:   Assembler les figures supplémentaires :
#         Supp 1 — RFP antibody signal intensity par cluster
#         Supp 2 — Volcano DEGs T cells (Gzmb) 6wpi vs 1wpi
#         Supp 3 — Volcano DEGs Mac (Ctss) 6wpi vs 1wpi
#         Supp 4 — Volcano DEGs Microglia (C1qa) 6wpi vs 1wpi
# Output: outputs/banksy/panels/supp/fig_supp.pdf
# =============================================================

suppressPackageStartupMessages({ library(magick) })
setwd(normalizePath("."))

ROOT    <- file.path("outputs", "banksy")
OUT_DIR <- file.path(ROOT, "panels", "supp")
OUT_FILE <- file.path(OUT_DIR, "fig_supp.pdf")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

PAGE_W_PX   <- 2480L
PAGE_H_PX   <- 1754L
HEADER_H_PX <- 120L   # légèrement plus grand pour le texte de légende
FONT_LBL    <- 27L
FONT_SRC    <- 17L
FONT_CAP    <- 19L

PAGES <- list(

  # --- Supp 1 : RFP signal ---
  list(
    type    = "side_by_side",
    label   = "Figure Supp 1 — RFP antibody signal intensity across annotated clusters",
    caption = paste0(
      "RFP antibody signal intensity was quantified across all annotated clusters ",
      "as a proxy for TRM contact history under the CRE-lox labeling system."
    ),
    paths = c(
      file.path(ROOT, "rfp_plots_readable", "fig_rfp_fc_all_clusters_barplot.jpg"),
      file.path(ROOT, "rfp_plots_readable", "fig_rfp_fc_all_clusters_heatmap.jpg")
    ),
    names = c("fig_rfp_fc_all_clusters_barplot.jpg",
              "fig_rfp_fc_all_clusters_heatmap.jpg")
  ),

  # --- Supp 2 : Volcano T cells (Gzmb) ---
  list(
    type    = "single",
    label   = "Figure Supp 2 — DEGs T cells (Gzmb) : LCMV 6wpi vs LCMV 1wpi",
    caption = "Full DEG tables for all sub-populations are provided in Supplementary Table 1.",
    path    = file.path(ROOT, "immune_niche_volcano_by_celltype",
                        "fig_volcano_t_cells_gzmb__lcmv_6wpi_vs_lcmv_1wpi.jpg")
  ),

  # --- Supp 3 : Volcano Mac (Ctss) ---
  list(
    type    = "single",
    label   = "Figure Supp 3 — DEGs Mac (Ctss) : LCMV 6wpi vs LCMV 1wpi",
    caption = "Full DEG tables for all sub-populations are provided in Supplementary Table 1.",
    path    = file.path(ROOT, "immune_niche_volcano_by_celltype",
                        "fig_volcano_mac_ctss__lcmv_6wpi_vs_lcmv_1wpi.jpg")
  ),

  # --- Supp 4 : Volcano Microglia (C1qa) ---
  list(
    type    = "single",
    label   = "Figure Supp 4 — DEGs Microglia (C1qa) : LCMV 6wpi vs LCMV 1wpi",
    caption = "Full DEG tables for all sub-populations are provided in Supplementary Table 1.",
    path    = file.path(ROOT, "immune_niche_volcano_by_celltype",
                        "fig_volcano_microglia_c1qa__lcmv_6wpi_vs_lcmv_1wpi.jpg")
  )
)

# --- helpers --------------------------------------------------
make_header <- function(label, caption = NULL) {
  img <- image_blank(PAGE_W_PX, HEADER_H_PX, color = "#F2F2F2") %>%
    image_annotate(label,
                   gravity = "NorthWest", location = "+20+12",
                   size = FONT_LBL, font = "sans",
                   weight = 700, color = "#1A1A1A")
  if (!is.null(caption) && nchar(caption) > 0) {
    img <- img %>%
      image_annotate(caption,
                     gravity = "NorthWest",
                     location = paste0("+20+", FONT_LBL + 20L),
                     size = FONT_CAP, font = "sans",
                     weight = 400, color = "#444444")
  }
  img
}

load_fit <- function(path, w_px, h_px) {
  if (!file.exists(path)) {
    warning("MISSING: ", path)
    return(
      image_blank(w_px, h_px, color = "white") %>%
        image_annotate(paste0("File not found:\n", basename(path)),
                       gravity = "Center", size = 28, color = "red3")
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

for (i in seq_along(PAGES)) {
  pg <- PAGES[[i]]
  message(sprintf("  Page %d/%d : %s", i, length(PAGES), pg$label))
  caption <- if (!is.null(pg$caption)) pg$caption else NULL

  if (pg$type == "single") {
    header     <- make_header(pg$label, caption)
    body       <- load_fit(pg$path, PAGE_W_PX, img_h_px)
    pages[[i]] <- image_append(c(header, body), stack = TRUE)

  } else if (pg$type == "side_by_side") {
    header    <- make_header(pg$label, caption)
    half_w    <- (PAGE_W_PX - 4L) %/% 2L
    left_img  <- load_fit(pg$paths[[1]], half_w, img_h_px)
    right_img <- load_fit(pg$paths[[2]], half_w, img_h_px)
    sep       <- image_blank(4L, img_h_px, color = "#CCCCCC")
    body      <- image_append(c(left_img, sep, right_img), stack = FALSE)
    pages[[i]] <- image_append(c(header, body), stack = TRUE)
  }
}

all_pages <- do.call(c, pages)
image_write(all_pages, path = OUT_FILE, format = "pdf", density = 200, quality = 85)
message("\n=== Sauvegarde : ", normalizePath(OUT_FILE))
message("    Pages : ", length(PAGES))
