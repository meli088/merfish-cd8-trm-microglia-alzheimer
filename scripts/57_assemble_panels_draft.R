#!/usr/bin/env Rscript

set.seed(1997)

suppressPackageStartupMessages({
  library(magick)
})

base_path <- normalizePath(".")
setwd(base_path)

banksy_root <- file.path("outputs", "banksy")
panels_dir <- file.path(banksy_root, "panels")
if (!dir.exists(panels_dir)) dir.create(panels_dir, recursive = TRUE)

PAGE_W <- 2480
PAGE_H <- 3508
MARGIN_X <- 70
HEADER_Y_LABEL <- 90
HEADER_Y_SOURCE <- 150
BODY_TOP <- 220
BODY_H <- PAGE_H - BODY_TOP - 70
BODY_W <- PAGE_W - 2 * MARGIN_X

build_not_found_tile <- function(msg, w, h) {
  tile <- image_blank(width = w, height = h, color = "white")
  tile <- image_annotate(tile,
                         text = msg,
                         gravity = "center",
                         size = 48,
                         color = "#444444")
  tile
}

fit_image_no_upscale <- function(img, slot_w, slot_h) {
  # Keep aspect ratio; shrink only if larger than slot, never enlarge.
  img <- image_resize(img, geometry = paste0(slot_w, "x", slot_h, ">"))
  info <- image_info(img)
  canvas <- image_blank(width = slot_w, height = slot_h, color = "white")
  ox <- floor((slot_w - info$width[1]) / 2)
  oy <- floor((slot_h - info$height[1]) / 2)
  image_composite(canvas, img, offset = paste0("+", ox, "+", oy))
}

best_match_from_candidates <- function(candidates, requested) {
  if (length(candidates) == 1) return(candidates[1])

  req_dir <- dirname(gsub("\\\\", "/", requested))
  req_tokens <- unique(unlist(strsplit(tolower(req_dir), "/", fixed = TRUE)))
  req_tokens <- req_tokens[nzchar(req_tokens)]

  cand_df <- data.frame(path = candidates, stringsAsFactors = FALSE)
  cand_df$path_low <- tolower(gsub("\\\\", "/", cand_df$path))
  cand_df$score <- vapply(cand_df$path_low, function(p) {
    sum(vapply(req_tokens, function(tok) grepl(tok, p, fixed = TRUE), logical(1)))
  }, numeric(1))
  cand_df$len <- nchar(cand_df$path)
  cand_df <- cand_df[order(-cand_df$score, cand_df$len), ]
  cand_df$path[1]
}

resolve_path <- function(requested_path) {
  req <- gsub("\\\\", "/", requested_path)
  if (file.exists(req)) return(req)

  rel <- req
  rel <- sub("^\\./", "", rel)
  rel <- sub("^outputs/banksy/", "", rel)

  candidate_rel <- file.path(banksy_root, rel)
  if (file.exists(candidate_rel)) return(candidate_rel)

  target_name <- basename(req)
  all_files <- list.files(banksy_root, recursive = TRUE, full.names = TRUE)
  matches <- all_files[tolower(basename(all_files)) == tolower(target_name)]

  if (length(matches) == 0) return(NA_character_)
  best_match_from_candidates(matches, requested_path)
}

compose_single_page <- function(panel_label, source_text, image_path = NULL, placeholder_text = NULL) {
  page <- image_blank(width = PAGE_W, height = PAGE_H, color = "white")
  page <- image_annotate(page,
                         text = panel_label,
                         gravity = "north",
                         location = paste0("+0+", HEADER_Y_LABEL),
                         size = 52,
                         color = "#111111")
  page <- image_annotate(page,
                         text = source_text,
                         gravity = "north",
                         location = paste0("+0+", HEADER_Y_SOURCE),
                         size = 30,
                         color = "#333333")

  if (!is.null(placeholder_text)) {
    tile <- build_not_found_tile(placeholder_text, BODY_W, BODY_H)
    page <- image_composite(page, tile, offset = paste0("+", MARGIN_X, "+", BODY_TOP))
    return(page)
  }

  if (is.null(image_path) || !isTRUE(file.exists(image_path))) {
    nf_msg <- paste0("FILE NOT FOUND: ", source_text)
    tile <- build_not_found_tile(nf_msg, BODY_W, BODY_H)
    page <- image_composite(page, tile, offset = paste0("+", MARGIN_X, "+", BODY_TOP))
    return(page)
  }

  img <- image_read(image_path)
  tile <- fit_image_no_upscale(img, BODY_W, BODY_H)
  page <- image_composite(page, tile, offset = paste0("+", MARGIN_X, "+", BODY_TOP))
  page
}

compose_double_page <- function(panel_label, source_left, source_right, path_left, path_right) {
  source_line <- paste0(source_left, " | ", source_right)
  page <- image_blank(width = PAGE_W, height = PAGE_H, color = "white")
  page <- image_annotate(page,
                         text = panel_label,
                         gravity = "north",
                         location = paste0("+0+", HEADER_Y_LABEL),
                         size = 52,
                         color = "#111111")
  page <- image_annotate(page,
                         text = source_line,
                         gravity = "north",
                         location = paste0("+0+", HEADER_Y_SOURCE),
                         size = 24,
                         color = "#333333")

  gap <- 40
  slot_w <- floor((BODY_W - gap) / 2)
  slot_h <- BODY_H
  left_x <- MARGIN_X
  right_x <- MARGIN_X + slot_w + gap

  get_tile <- function(path, source_text) {
    if (is.na(path) || !file.exists(path)) {
      return(build_not_found_tile(paste0("FILE NOT FOUND: ", source_text), slot_w, slot_h))
    }
    img <- image_read(path)
    fit_image_no_upscale(img, slot_w, slot_h)
  }

  left_tile <- get_tile(path_left, source_left)
  right_tile <- get_tile(path_right, source_right)

  page <- image_composite(page, left_tile, offset = paste0("+", left_x, "+", BODY_TOP))
  page <- image_composite(page, right_tile, offset = paste0("+", right_x, "+", BODY_TOP))
  page
}

assemble_pdf <- function(output_pdf, pages_spec) {
  pages <- list()

  for (pg in pages_spec) {
    if (pg$type == "placeholder") {
      pages[[length(pages) + 1]] <- compose_single_page(
        panel_label = pg$label,
        source_text = pg$source,
        placeholder_text = pg$placeholder
      )
      next
    }

    if (pg$type == "double") {
      p_left <- resolve_path(pg$left)
      p_right <- resolve_path(pg$right)
      pages[[length(pages) + 1]] <- compose_double_page(
        panel_label = pg$label,
        source_left = pg$left,
        source_right = pg$right,
        path_left = p_left,
        path_right = p_right
      )
      next
    }

    resolved <- resolve_path(pg$source)
    pages[[length(pages) + 1]] <- compose_single_page(
      panel_label = pg$label,
      source_text = pg$source,
      image_path = if (is.na(resolved)) NULL else resolved
    )
  }

  all_pages <- do.call(c, pages)
  image_write(all_pages, path = output_pdf, format = "pdf", density = "300")
  message("Created: ", output_pdf)
}

fig1_pages <- list(
  list(type = "placeholder", label = "Fig 1 - Panel A", source = "[SCHEMA EXPERIMENTAL]", placeholder = "To be added manually"),
  list(type = "image", label = "Fig 1 - Panel B", source = "outputs/banksy/umap_annotated/UMAP_annotated_from_ncellscsv_lam02_res09.jpg"),
  list(type = "image", label = "Fig 1 - Panel C", source = "outputs/banksy/umap_annotated/heatmap_annotated_lam02_res09_top2_OK.jpg"),
  list(type = "image", label = "Fig 1 - Panel D", source = "outputs/banksy/spatial_annotations/dotplot_canonical_markers_revised_v2.jpg"),
  list(type = "image", label = "Fig 1 - Panel E", source = "outputs/banksy/spatial_annotations/xy_spatial_all_annotations_combined_lam02_res09.jpg"),
  list(type = "image", label = "Fig 1 - Panel F", source = "outputs/banksy/umap_annotated/enrichment_vs_mock_top10_OK.jpg"),
  list(type = "image", label = "Fig 1 - Panel G", source = "outputs/banksy/inflammatory_niche_step2_grid_100um_global/fig_grid_immune_ratio_barplot.jpg"),
  list(type = "image", label = "Fig 1 - Panel H", source = "outputs/banksy/rfp_plots_readable/fig_rfp_fc_summary_5groups.jpg")
)

fig2_pages <- list(
  list(type = "image", label = "Fig 2 - Panel A", source = "outputs/banksy/immune_acod1/analysis/figures/umap_merged_res03.jpg"),
  list(type = "image", label = "Fig 2 - Panel B", source = "outputs/banksy/immune_acod1/analysis/figures/heatmap_immune_deg_top2_res03.jpg"),
  list(type = "image", label = "Fig 2 - Panel C", source = "outputs/banksy/immune_acod1/analysis/figures/dotplot_canonical_markers_res03_OK_v2.jpg"),
  list(type = "image", label = "Fig 2 - Panel D", source = "outputs/banksy/spatial_annotations/composition_immune_vs_nonimmune_overtime.jpg"),
  list(type = "image", label = "Fig 2 - Panel E", source = "outputs/banksy/immune_acod1/analysis/figures/composition_stacked_res03_OK.jpg"),
  list(type = "placeholder", label = "Fig 2 - Panel F", source = "[PLACEHOLDER]", placeholder = "XY zoom niche highlighted - to be remade"),
  list(type = "image", label = "Fig 2 - Panel G", source = "outputs/banksy/microglia_dam_niche/fig_module_score_merged_violin_upregulated_dam.jpg"),
    list(type = "double", label = "Fig 2 - Panel H",
      left = "outputs/banksy/microglia_dam_niche/fig_dam_score_violin_in_vs_out_niche.jpg",
      right = "outputs/banksy/microglia_dam_niche/fig_dam_score_barplot_in_vs_out_niche.jpg"),
  list(type = "double", label = "Fig 2 - Panel I",
       left = "outputs/banksy/microglia_dam_niche/fig_distance_microglia_to_tcell_overtime.jpg",
       right = "outputs/banksy/microglia_dam_niche/fig_dam_score_by_distance_bins.jpg")
)

fig3_pages <- list(
  list(type = "image", label = "Fig 3 - Panel J", source = "outputs/banksy/ifn_immune_overlay/fig1_xy_ifn_immune_overlay_OK.jpg"),
  list(type = "image", label = "Fig 3 - Panel K", source = "outputs/banksy/nearest_immune_distance_OK/violin_nearest_immune_distance_lam02_res09.jpg"),
  list(type = "image", label = "Fig 3 - Panel L", source = "outputs/banksy/immune_acod1/analysis/figures/umap_merged_res03.jpg"),
  list(type = "image", label = "Fig 3 - Panel M", source = "outputs/banksy/spatial_annotations/xy_map2_ptprc_all_samples.jpg"),
  list(type = "placeholder", label = "Fig 3 - Panel N", source = "[PLACEHOLDER]", placeholder = "Bubble Ifng global clusters - to be remade"),
  list(type = "image", label = "Fig 3 - Panel O", source = "outputs/banksy/ifn_immune_overlay/fig_ifng_dotplot_all_celltypes.jpg"),
  list(type = "image", label = "Fig 3 - Panel P", source = "outputs/banksy/ifn_immune_overlay/fig_ifng_tcell_overtime.jpg")
)

supp_pages <- list(
  list(type = "image", label = "Supp S1", source = "outputs/banksy/rfp_plots_readable/fig_rfp_fc_all_clusters_barplot.jpg"),
  list(type = "image", label = "Supp S2", source = "outputs/banksy/immune_niche_volcano_by_celltype/fig_ndeg_summary_by_celltype.jpg"),
  list(type = "image", label = "Supp S3", source = "outputs/banksy/immune_niche_volcano_by_celltype/fig_volcano_t_cells_gzmb__lcmv_6wpi_vs_lcmv_1wpi.jpg"),
  list(type = "image", label = "Supp S4", source = "outputs/banksy/immune_niche_volcano_by_celltype/fig_volcano_mac_ctss__lcmv_6wpi_vs_lcmv_1wpi.jpg"),
  list(type = "image", label = "Supp S5", source = "outputs/banksy/immune_niche_volcano_by_celltype/fig_volcano_microglia_c1qa__lcmv_6wpi_vs_lcmv_1wpi.jpg"),
  list(type = "image", label = "Supp S6", source = "outputs/banksy/microglia_dam_niche/fig_volcano_in_vs_out_lcmv_1wpi.jpg"),
  list(type = "image", label = "Supp S7", source = "outputs/banksy/microglia_dam_niche/fig_volcano_in_vs_out_lcmv_3wpi.jpg"),
  list(type = "image", label = "Supp S8", source = "outputs/banksy/microglia_dam_niche/fig_volcano_in_vs_out_lcmv_6wpi.jpg"),
  list(type = "image", label = "Supp S9", source = "outputs/banksy/microglia_dam_niche/fig_module_score_merged_violin_downregulated_dam.jpg"),
  list(type = "image", label = "Supp S10", source = "outputs/banksy/microglia_conditions_deg/fig_volcano_microglia_lcmv_1wpi_vs_mock_6wpi.jpg"),
  list(type = "image", label = "Supp S11", source = "outputs/banksy/microglia_conditions_deg/fig_volcano_microglia_lcmv_6wpi_vs_mock_6wpi.jpg"),
  list(type = "image", label = "Supp S12", source = "outputs/banksy/microglia_conditions_deg/fig_volcano_microglia_lcmv_6wpi_vs_lcmv_1wpi.jpg")
)

assemble_pdf(file.path(panels_dir, "fig1_panel_draft.pdf"), fig1_pages)
assemble_pdf(file.path(panels_dir, "fig2_panel_draft.pdf"), fig2_pages)
assemble_pdf(file.path(panels_dir, "fig3_panel_draft.pdf"), fig3_pages)
assemble_pdf(file.path(panels_dir, "supp_figures_all_draft.pdf"), supp_pages)

message("Done.")
