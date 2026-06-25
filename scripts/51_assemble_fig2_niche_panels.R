#!/usr/bin/env Rscript
# =============================================================
# Script: 51_assemble_fig2_niche_panels.R
# Goal:   PDF de revue — figures niche immune + microglie
#         (composition, spatial, volcano, module score, distance)
# Output: outputs/banksy/panels/fig2_niche_panels.pdf
# =============================================================

suppressPackageStartupMessages({ library(magick) })
setwd(normalizePath("."))

ROOT     <- file.path("outputs", "banksy")
OUT_DIR  <- file.path(ROOT, "panels")
OUT_FILE <- file.path(OUT_DIR, "fig2_niche_panels.pdf")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

PAGE_W_PX   <- 2480L
PAGE_H_PX   <- 1754L
HEADER_H_PX <- 92L
FONT_LBL    <- 27L
FONT_SRC    <- 17L

PAGES <- list(
  list(
    type  = "single",
    label = "Composition immune vs non-immune par timepoint (global)",
    path  = file.path(ROOT, "umap_annotated",
                      "composition_immune_vs_nonimmune_overtime.jpg")
  ),
  list(
    type  = "single",
    label = "Composition sous-clusters immune par timepoint (res03)",
    path  = file.path(ROOT, "immune_acod1", "analysis", "figures",
                      "composition_stacked_res03.jpg")
  ),
  list(
    type   = "triple",
    label  = "XY zoom niche immune — 1wpi | 3wpi | 6wpi",
    paths  = c(
      file.path(ROOT, "spatial_annotations", "xy_zoom_immune_niche_lcmv_1wpi.jpg"),
      file.path(ROOT, "spatial_annotations", "xy_zoom_immune_niche_lcmv_3wpi.jpg"),
      file.path(ROOT, "spatial_annotations", "xy_zoom_immune_niche_lcmv_6wpi.jpg")
    ),
    names  = c("xy_zoom_immune_niche_lcmv_1wpi.jpg",
               "xy_zoom_immune_niche_lcmv_3wpi.jpg",
               "xy_zoom_immune_niche_lcmv_6wpi.jpg")
  ),
  list(
    type  = "single",
    label = "Volcano microglie in vs out niche (tous timepoints)",
    path  = file.path(ROOT, "microglia_dam_niche",
                      "fig_volcano_in_vs_out_all_timepoints_merged.jpg")
  ),
  list(
    type  = "single",
    label = "Module score DAM up-regulated — violin par timepoint",
    path  = file.path(ROOT, "microglia_dam_niche",
                      "fig_module_score_merged_violin_upregulated_dam.jpg")
  ),
  list(
    type   = "side_by_side",
    label  = "Distance microglie -> T cells + DAM score par bins de distance",
    paths  = c(
      file.path(ROOT, "microglia_dam_niche",
                "fig_distance_microglia_to_tcell_overtime.jpg"),
      file.path(ROOT, "microglia_dam_niche",
                "fig_dam_score_by_distance_bins.jpg")
    ),
    names  = c("fig_distance_microglia_to_tcell_overtime.jpg",
               "fig_dam_score_by_distance_bins.jpg")
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
message("\n=== Assemblage : ", OUT_FILE)
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

  } else if (pg$type == "triple") {
    header  <- make_header(pg$label, paste(pg$names, collapse = "  |  "))
    third_w <- (PAGE_W_PX - 8L) %/% 3L
    img1 <- load_fit(pg$paths[[1]], third_w, img_h_px)
    img2 <- load_fit(pg$paths[[2]], third_w, img_h_px)
    img3 <- load_fit(pg$paths[[3]], third_w, img_h_px)
    sep  <- image_blank(4L, img_h_px, color = "#CCCCCC")
    body <- image_append(c(img1, sep, img2, sep, img3), stack = FALSE)
    pages[[i]] <- image_append(c(header, body), stack = TRUE)
  }
}

all_pages <- do.call(c, pages)
image_write(all_pages, path = OUT_FILE, format = "pdf", density = 200, quality = 85)
message("\n=== Sauvegarde : ", normalizePath(OUT_FILE))
message("    Pages : ", length(PAGES))
