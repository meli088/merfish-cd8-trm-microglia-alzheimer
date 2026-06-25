#!/usr/bin/env Rscript

set.seed(1997)

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(tidyverse)
})

base_path <- normalizePath(".")
setwd(base_path)

source("scripts/00_palette.R")

safe_label <- "immune_acod1"
obj07 <- file.path("objects", "07_immune_banksy_lam02_immune_acod1.rds")
if (!file.exists(obj07)) obj07 <- file.path("objects", "07_immune_banksy_lam02_after_bloc1.rds")
obj08 <- file.path("objects", "08_immune_annotated_lam02_res03.rds")

out_dir <- file.path("outputs", "banksy", safe_label, "analysis", "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

find_cl_col <- function(se_obj, lam, res) {
  cd <- SummarizedExperiment::colData(se_obj)
  cn <- colnames(cd)
  res_txt <- sprintf("%.1f", res)
  lam_pat <- paste0("lam", gsub("\\.", "\\\\.", as.character(lam)))
  res_pat <- paste0("_res", gsub("\\.", "\\\\.", res_txt), "$")
  cols <- cn[grepl(lam_pat, cn) & grepl(res_pat, cn)]
  if (length(cols) == 0) return(NULL)
  cols[1]
}

build_ct_colors <- function(ct_levels) {
  palette20 <- c(
    "#ff0004", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
    "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
    "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
    "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A"
  )
  cols <- setNames(rep(palette20, length.out = length(ct_levels)), ct_levels)
  for (ct in ct_levels) {
    if (ct %in% names(GLOBAL_PALETTE)) cols[[ct]] <- GLOBAL_PALETTE[[ct]]
  }
  cols
}

message("Loading: ", obj07)
se <- readRDS(obj07)

# Same behavior as script 14: focus on LCMV time course only.
se <- se[, as.character(SummarizedExperiment::colData(se)$sample) != "mock_6wpi"]

assay_name <- if ("counts" %in% assayNames(se)) "counts" else assayNames(se)[1]
so <- as.Seurat(se, counts = assay_name, data = NULL)
default_assay <- if ("RNA" %in% SeuratObject::Assays(so)) "RNA" else SeuratObject::Assays(so)[1]
DefaultAssay(so) <- default_assay

needs_norm <- tryCatch({
  dm <- GetAssayData(so, assay = default_assay, layer = "data")
  nrow(dm) == 0 || ncol(dm) == 0
}, error = function(e) TRUE)
if (needs_norm) {
  so <- NormalizeData(so, assay = default_assay, verbose = FALSE)
}

if (!"umap" %in% names(so@reductions)) {
  if ("UMAP" %in% names(so@reductions)) {
    so[["umap"]] <- so[["UMAP"]]
  } else {
    avail <- names(so@reductions)
    harm_red <- avail[grepl("^[Hh]armony", avail)][1]
    pca_red <- avail[grepl("^PCA", avail)][1]
    umap_src <- if (!is.na(harm_red)) harm_red else if (!is.na(pca_red)) pca_red else avail[1]
    n_dims <- min(40L, ncol(so@reductions[[umap_src]]))
    so <- RunUMAP(so, reduction = umap_src, dims = seq_len(n_dims), verbose = FALSE)
  }
}

umap_mat <- as.data.frame(so@reductions[["umap"]]@cell.embeddings)
colnames(umap_mat) <- c("UMAP1", "UMAP2")
umap_mat$sample <- so@meta.data$sample

# Override labels from object 08 with latest annotations.
cell_type_override <- NULL
if (file.exists(obj08)) {
  se08 <- readRDS(obj08)
  if ("cell_type" %in% colnames(colData(se08))) {
    common_cells <- intersect(colnames(so), colnames(se08))
    if (length(common_cells) > 0) {
      cell_type_override <- setNames(
        as.character(colData(se08)[common_cells, "cell_type"]),
        common_cells
      )
    }
  }
}

res <- 0.3
res_id <- gsub("\\.", "", sprintf("%.1f", res))
cl_col <- find_cl_col(se, 0.2, res)
if (is.null(cl_col)) stop("No clust column found for res=0.3")

raw_levels <- paste0("Domain_", sort(unique(as.integer(colData(se)[[cl_col]]))))
domain_labels <- paste0("Domain_", as.character(colData(se)[[cl_col]]))

cell_disp <- domain_labels
if (!is.null(cell_type_override)) {
  so_cells <- colnames(so)
  ovr_idx <- which(so_cells %in% names(cell_type_override))
  if (length(ovr_idx) > 0) {
    cell_disp[ovr_idx] <- cell_type_override[so_cells[ovr_idx]]
  }
}

bio_lvls <- unique(cell_disp[!grepl("^Domain_", cell_disp)])
dom_lvls <- unique(cell_disp[grepl("^Domain_", cell_disp)])
disp_order <- c(order_annotations(bio_lvls), dom_lvls)

udf <- umap_mat %>%
  mutate(
    cluster = factor(cell_disp, levels = disp_order),
    sample = factor(sample, levels = c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi"))
  )

dom_colors <- build_ct_colors(levels(udf$cluster))

p_grid <- ggplot(udf, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(size = 0.8, alpha = 0.7, stroke = 0) +
  scale_color_manual(values = dom_colors, drop = FALSE) +
  facet_wrap(~sample, ncol = 2) +
  labs(
    title = "Immune (Acod1) — BANKSY sub-clusters (res=0.3)",
    color = "Cell type"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 11)
  )

p_merged <- ggplot(udf, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(size = 0.8, alpha = 0.7, stroke = 0) +
  scale_color_manual(values = dom_colors, drop = FALSE) +
  labs(
    title = "Immune (Acod1) — BANKSY sub-clusters (res=0.3) — all timepoints",
    color = "Cell type"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
  theme_classic(base_size = 10) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 11)
  )

ggsave(file.path(out_dir, paste0("umap_grid_res", res_id, ".pdf")), p_grid, width = 10, height = 8)
ggsave(file.path(out_dir, paste0("umap_grid_res", res_id, ".jpg")), p_grid, width = 10, height = 8, dpi = 300)
ggsave(file.path(out_dir, paste0("umap_merged_res", res_id, ".pdf")), p_merged, width = 7, height = 6)
ggsave(file.path(out_dir, paste0("umap_merged_res", res_id, ".jpg")), p_merged, width = 7, height = 6, dpi = 300)

message("Saved: ", file.path(out_dir, paste0("umap_grid_res", res_id, ".pdf/.jpg")))
message("Saved: ", file.path(out_dir, paste0("umap_merged_res", res_id, ".pdf/.jpg")))
