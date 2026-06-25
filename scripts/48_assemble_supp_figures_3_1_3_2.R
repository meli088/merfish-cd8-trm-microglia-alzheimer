#!/usr/bin/env Rscript
# =============================================================
# Script: 48_assemble_supp_figures_3_1_3_2.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-06
#
# Goal:
#   Assembler les figures supplémentaires des parties 3.1 et 3.2
#   en un seul PDF multi-pages (une figure par page).
#   Chaque page : label "Supp Fig Sx" + nom du fichier source
#   en en-tête, figure en dessous.
#
# Output : outputs/banksy/panels/supp_figures_3_1_3_2_draft.pdf
# =============================================================

set.seed(1997)

suppressPackageStartupMessages({
  library(magick)
})

setwd(normalizePath("."))

OUT_DIR  <- file.path("outputs", "banksy", "panels")
OUT_FILE <- file.path(OUT_DIR, "supp_figures_3_1_3_2_draft.pdf")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

ROOT <- file.path("outputs", "banksy")

# =============================================================
# Définition des pages (label, chemin)
# =============================================================

PAGES <- list(
  # --- PARTIE 3.1 ---
  list(
    label = "Supp Fig S1 — RFP log2FC all clusters (barplot)",
    path  = file.path(ROOT, "rfp_plots_readable",
                      "fig_rfp_fc_all_clusters_barplot.jpg")
  ),
  # --- PARTIE 3.2 ---
  list(
    label = "Supp Fig S2 — nDEG summary by domain (immune subclusters res03)",
    path  = file.path(ROOT, "immune_acod1", "domain_DEG_evolution",
                      "res03", "heatmaps", "fig_ndeg_summary_by_domain.jpg")
  ),
  list(
    label = "Supp Fig S3 — Volcano T cells (Gzmb) : LCMV 6wpi vs 1wpi",
    path  = file.path(ROOT, "immune_niche_volcano_by_celltype",
                      "fig_volcano_t_cells_gzmb__lcmv_6wpi_vs_lcmv_1wpi_OK.jpg")
  ),
  list(
    label = "Supp Fig S4 — Volcano Mac (Ctss) : LCMV 6wpi vs 1wpi",
    path  = file.path(ROOT, "immune_niche_volcano_by_celltype",
                      "fig_volcano_mac_ctss__lcmv_6wpi_vs_lcmv_1wpi_OK.jpg")
  ),
  list(
    label = "Supp Fig S5 — Volcano Microglia (C1qa) : LCMV 6wpi vs 1wpi",
    path  = file.path(ROOT, "immune_niche_volcano_by_celltype",
                      "fig_volcano_microglia_c1qa__lcmv_6wpi_vs_lcmv_1wpi.jpg")
  ),
  list(
    label = "Supp Fig S6 — Volcano Microglia in vs out niche : LCMV 1wpi",
    path  = file.path(ROOT, "microglia_dam_niche",
                      "fig_volcano_in_vs_out_lcmv_1wpi.jpg")
  ),
  list(
    label = "Supp Fig S7 — Volcano Microglia in vs out niche : LCMV 3wpi",
    path  = file.path(ROOT, "microglia_dam_niche",
                      "fig_volcano_in_vs_out_lcmv_3wpi.jpg")
  ),
  list(
    label = "Supp Fig S8 — Volcano Microglia in vs out niche : LCMV 6wpi",
    path  = file.path(ROOT, "microglia_dam_niche",
                      "fig_volcano_in_vs_out_lcmv_6wpi.jpg")
  ),
  list(
    label = "Supp Fig S9 — Module score violin (down-regulated DAM)",
    path  = file.path(ROOT, "microglia_dam_niche",
                      "fig_module_score_merged_violin_downregulated_dam.jpg")
  ),
  list(
    label = "Supp Fig S10 — Volcano Microglia : LCMV 1wpi vs Mock 6wpi",
    path  = file.path(ROOT, "microglia_conditions_deg",
                      "fig_volcano_microglia_lcmv_1wpi_vs_mock_6wpi.jpg")
  ),
  list(
    label = "Supp Fig S11 — Volcano Microglia : LCMV 6wpi vs Mock 6wpi",
    path  = file.path(ROOT, "microglia_conditions_deg",
                      "fig_volcano_microglia_lcmv_6wpi_vs_mock_6wpi.jpg")
  ),
  list(
    label = "Supp Fig S12 — Volcano Microglia : LCMV 6wpi vs 1wpi",
    path  = file.path(ROOT, "microglia_conditions_deg",
                      "fig_volcano_microglia_lcmv_6wpi_vs_lcmv_1wpi.jpg")
  )
)

# =============================================================
# Assemblage PDF via magick — résolution native des images
# =============================================================
# Dimensions page cible en pixels (300 dpi, A4 landscape)
PAGE_W_PX <- 3508L   # 11.69 * 300
PAGE_H_PX <- 2480L   # 8.27  * 300
HEADER_H_PX <- 130L  # hauteur du bandeau en-tête en pixels
FONT_SIZE_LBL <- 38L
FONT_SIZE_SRC <- 24L

message("\n=== Assemblage HD : ", OUT_FILE)

pages <- list()

for (i in seq_along(PAGES)) {
  pg    <- PAGES[[i]]
  label <- pg$label
  path  <- pg$path

  message(sprintf("  Page %2d / %d : %s", i, length(PAGES), label))

  # ── Bandeau en-tête ──────────────────────────────────────────
  header <- image_blank(PAGE_W_PX, HEADER_H_PX, color = "#F2F2F2") %>%
    image_annotate(label,
                   gravity   = "NorthWest",
                   location  = "+20+12",
                   size      = FONT_SIZE_LBL,
                   font      = "sans",
                   weight    = 700,
                   color     = "#1A1A1A") %>%
    image_annotate(paste0("Source: ", basename(path)),
                   gravity   = "NorthWest",
                   location  = paste0("+20+", FONT_SIZE_LBL + 16L),
                   size      = FONT_SIZE_SRC,
                   font      = "sans",
                   weight    = 400,
                   color     = "#666666")

  # ── Image source ─────────────────────────────────────────────
  img_h_px <- PAGE_H_PX - HEADER_H_PX

  if (!file.exists(path)) {
    warning("  MISSING: ", path)
    body <- image_blank(PAGE_W_PX, img_h_px, color = "white") %>%
      image_annotate(paste0("File not found:\n", path),
                     gravity = "Center", size = 28,
                     color = "red3")
  } else {
    raw_img <- image_read(path)
    # Redimensionner pour tenir dans la zone image en conservant le ratio
    body <- image_resize(raw_img,
                         paste0(PAGE_W_PX, "x", img_h_px, ">"))
    # Coller sur fond blanc si plus petit que la zone
    body <- image_extent(body,
                         paste0(PAGE_W_PX, "x", img_h_px),
                         gravity = "Center",
                         color   = "white")
  }

  # ── Assemblage vertical ──────────────────────────────────────
  page <- image_append(c(header, body), stack = TRUE)
  pages[[i]] <- page
}

# ── Écriture PDF multi-pages ──────────────────────────────────
all_pages <- do.call(c, pages)
image_write(all_pages, path = OUT_FILE, format = "pdf",
            density = 300, quality = 100)

message("\n=== Sauvegardé : ", normalizePath(OUT_FILE))
message("    Pages : ", length(PAGES))
