#!/usr/bin/env Rscript
# =============================================================
# Script: 54_assemble_fig3_panels.R
# Goal:   Assembler les panneaux de la Fig 3 en PDF multi-pages
# Output: outputs/banksy/panels/fig3_panel_draft.pdf
# =============================================================

suppressPackageStartupMessages({ library(magick) })
setwd(normalizePath("."))

ROOT     <- file.path("outputs", "banksy")
OUT_DIR  <- file.path(ROOT, "panels")
OUT_FILE <- file.path(OUT_DIR, "fig3_panel_draft.pdf")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

PAGE_W_PX   <- 2480L
PAGE_H_PX   <- 1754L
HEADER_H_PX <- 92L
FONT_LBL    <- 27L
FONT_SRC    <- 17L

PAGES <- list(
  list(
    type  = "single",
    label = "Fig 3 — Panel J : XY overlay IFN + cellules immunes",
    path  = file.path(ROOT, "ifn_immune_overlay",
                      "fig1_xy_ifn_immune_overlay_OK.jpg")
  ),
  list(
    type  = "single",
    label = "Fig 3 — Panel K : Violin distance aux cellules immunes (lam02 res09)",
    path  = file.path(ROOT, "nearest_immune_distance_OK",
                      "violin_nearest_immune_distance_lam02_res09.jpg")
  ),
  list(
    type  = "single",
    label = "Fig 3 — Panel L : UMAP sous-clustering IFN — tous timepoints fusionnés (res03)",
    path  = file.path(ROOT, "ifn", "analysis", "figures",
                      "umap_merged_res03.jpg")
  ),
  list(
    type  = "single",
    label = "Fig 3 — Panel M : XY Map2 + Ptprc (tous échantillons)",
    path  = file.path(ROOT, "spatial_annotations",
                      "xy_map2_ptprc_all_samples.jpg")
  ),
  list(
    type   = "side_by_side",
    label  = "Fig 3 — Panel N : Grid IFN faceted (shared scale)  |  Immune ratio barplot",
    paths  = c(
      file.path(ROOT, "inflammatory_niche_step2_grid_100um_global",
                "fig_grid_ifn_faceted_shared_scale.jpg"),
      file.path(ROOT, "inflammatory_niche_step2_grid_100um_global",
                "fig_grid_immune_ratio_barplot.jpg")
    ),
    names  = c("fig_grid_ifn_faceted_shared_scale.jpg",
               "fig_grid_immune_ratio_barplot.jpg")
  ),
  list(
    type  = "single",
    label = "Fig 3 — Panel O : IFNg dotplot tous types cellulaires",
    path  = file.path(ROOT, "ifn_immune_overlay",
                      "fig_ifng_dotplot_all_celltypes.jpg")
  ),
  list(
    type  = "single",
    label = "Fig 3 — Panel P : IFNg bubble plot",
    path  = file.path(ROOT, "ifn_immune_overlay",
                      "fig1_ifng_bubble_OK.jpg")
  ),
  list(
    type  = "single",
    label = "Fig 3 — Panel Q : IFNg T cells overtime",
    path  = file.path(ROOT, "ifn_immune_overlay",
                      "fig_ifng_tcell_overtime.jpg")
  )
)

# --- helpers --------------------------------------------------
make_header <- function(label, src_txt) {
  image_blank(PAGE_W_PX, HEADER_H_PX, color = "#F2F2F2") %>%
    image_annotate(label,
                   gravity = "NorthWest", location = "+20+12",
                   size = FONT_LBL, font = "sans",
                   weight = 700, color = "#1A1A1A") %>%
    image_annotate(paste0("Source: ", src_txt),
                   gravity = "NorthWest",
                   location = paste0("+20+", FONT_LBL + 16L),
                   size = FONT_SRC, font = "sans",
                   weight = 400, color = "#666666")
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
message("\n=== Assemblage HD : ", OUT_FILE)
img_h_px <- PAGE_H_PX - HEADER_H_PX
pages <- list()

for (i in seq_along(PAGES)) {
  pg <- PAGES[[i]]
  message(sprintf("  Page %d/%d : %s", i, length(PAGES), pg$label))

  if (pg$type == "single") {
    header     <- make_header(pg$label, basename(pg$path))
    body       <- load_fit(pg$path, PAGE_W_PX, img_h_px)
    pages[[i]] <- image_append(c(header, body), stack = TRUE)

  } else if (pg$type == "side_by_side") {
    header    <- make_header(pg$label, paste(pg$names, collapse = "  |  "))
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
message("\n=== Sauvegardé : ", normalizePath(OUT_FILE))
message("    Pages : ", length(PAGES))
