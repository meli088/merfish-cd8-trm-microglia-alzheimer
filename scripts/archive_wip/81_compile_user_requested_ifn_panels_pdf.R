#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jpeg)
  library(grid)
})

base_path <- normalizePath(".")
setwd(base_path)

img_paths <- c(
  "outputs/banksy/ifn_immune_overlay/fig_ifn_nearest_immune_distance_violin.jpg",
  "outputs/banksy/umap_annotated/UMAP_annotated_from_ncellscsv_lam02_res09.jpg",
  "outputs/banksy/ifn_immune_overlay/fig_xy_ifn_same_scale_all_conditions.jpg",
  "outputs/banksy/ifn_immune_overlay/fig_ifng_dotplot_all_celltypes_timepoints.jpg"
)

missing <- img_paths[!file.exists(img_paths)]
if (length(missing) > 0) {
  stop("Missing input figure(s):\n", paste(missing, collapse = "\n"))
}

out_pdf <- "outputs/banksy/ifn_immune_overlay/fig_ifn_requested_4panels_combined.pdf"

imgs <- lapply(img_paths, jpeg::readJPEG)

cairo_pdf(filename = out_pdf, width = 16, height = 11)
grid::grid.newpage()
lay <- grid::grid.layout(nrow = 2, ncol = 2)
grid::pushViewport(grid::viewport(layout = lay))

panel_pos <- list(c(1, 1), c(1, 2), c(2, 1), c(2, 2))
for (i in seq_along(imgs)) {
  rr <- panel_pos[[i]][1]
  cc <- panel_pos[[i]][2]
  grid::pushViewport(grid::viewport(layout.pos.row = rr, layout.pos.col = cc))
  grid::grid.raster(imgs[[i]], interpolate = TRUE)
  grid::popViewport()
}

dev.off()
message("Saved: ", out_pdf)
