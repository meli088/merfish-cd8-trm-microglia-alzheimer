#!/usr/bin/env Rscript
# =============================================================
# Script: 50_immune_deg_heatmap.R
# Project: LCMV MERFISH â€” TRM-Microglia niche analysis
# Author: MÃ©lina Farshchi
# Date: 2026-06
#
# Goal:
#   Heatmap DEG top2 par cluster immune (res=0.3), mÃªme format que
#   heatmap_annotated_lam02_res09_top2 :
#   AverageExpression â†’ z-score â†’ pheatmap avec gaps_row par cluster.
#
# Output: outputs/banksy/immune_acod1/analysis/figures/
#           heatmap_immune_deg_top2_res03.jpg / .pdf
# =============================================================

set.seed(1997)

suppressPackageStartupMessages({
  library(Seurat)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(Banksy)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(Cairo)
})

base_path <- normalizePath(".")
setwd(base_path)
`%in%` <- base::`%in%`
Assays  <- SeuratObject::Assays
source("scripts/00_palette.R")

# -------------------------------------------------------
# ParamÃ¨tres
# -------------------------------------------------------
TOP_N      <- 2
OUT_DIR    <- file.path("outputs", "banksy", "immune_acod1", "analysis", "figures")
ANNOT_FILE <- file.path("outputs", "banksy", "immune_acod1",
                         "annotation_map_immune_res03.csv")
OBJ_FILE   <- file.path("objects", "08_immune_annotated_lam02_res03.rds")

# -------------------------------------------------------
# Helper heatmap (identique au script 05c)
# -------------------------------------------------------
save_heatmap_pair <- function(avg_mat, title, out_base, gaps_row = NULL,
                              labels_row = NULL, w = 12, h = 10,
                              label_scale = 1) {
  if (is.null(avg_mat) || nrow(avg_mat) == 0 || ncol(avg_mat) == 0) {
    message("Heatmap ignorÃ©e: matrice vide"); return(invisible(NULL))
  }
  avg_scaled <- t(scale(t(avg_mat)))
  avg_scaled <- pmin(pmax(avg_scaled, -2), 2)
  avg_scaled[is.nan(avg_scaled)] <- 0

  n_genes  <- nrow(avg_scaled)
  n_groups <- ncol(avg_scaled)
  w_auto   <- max(w, 7 + 0.6 * n_groups)
  h_auto   <- max(h, 5 + 0.2 * n_genes)
  fs_row   <- max(5, min(9, 260 / max(25, n_genes))) * label_scale
  fs_col   <- max(7, min(10, 130 / max(8, n_groups))) * label_scale

  labels_col_wrapped <- vapply(
    colnames(avg_scaled),
    function(x) paste(strwrap(x, width = 14), collapse = "\n"),
    character(1)
  )

  render_heatmap <- function() {
    pheatmap(
      avg_scaled,
      color          = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      cluster_rows   = FALSE,
      cluster_cols   = FALSE,
      show_rownames  = TRUE,
      show_colnames  = TRUE,
      labels_row     = labels_row,
      labels_col     = labels_col_wrapped,
      fontsize_row   = fs_row,
      fontsize_col   = fs_col,
      angle_col      = 45,
      border_color   = NA,
      main           = title,
      gaps_row       = gaps_row,
      legend         = FALSE
    )
  }

  out_pdf <- paste0(out_base, ".pdf")
  out_jpg <- paste0(out_base, ".jpg")

  CairoPDF(out_pdf, width = w_auto, height = h_auto)
  render_heatmap()
  dev.off()

  CairoJPEG(out_jpg, width = round(w_auto * 150), height = round(h_auto * 150),
            res = 150, quality = 95)
  render_heatmap()
  dev.off()

  message("Saved: ", basename(out_pdf))
  message("Saved: ", basename(out_jpg))
}

# -------------------------------------------------------
# Charger l'objet 08 + convertir en Seurat
# -------------------------------------------------------
message("Loading: ", OBJ_FILE)
se <- readRDS(OBJ_FILE)
message("  ", ncol(se), " cells")

assay_nm <- if ("counts" %in% assayNames(se)) "counts" else assayNames(se)[1]
cnt      <- methods::as(assay(se, assay_nm), "dgCMatrix")
rownames(cnt) <- rownames(se)
colnames(cnt) <- colnames(se)

so <- suppressWarnings(CreateSeuratObject(counts = cnt, project = "immune"))
rm(cnt)
DefaultAssay(so) <- "RNA"
so <- NormalizeData(so, verbose = FALSE)
so <- JoinLayers(so)

# TransfÃ©rer cell_type depuis colData
so$cell_type <- as.character(colData(se)$cell_type)

# -------------------------------------------------------
# Annotation map : mapping Domain_X -> nom biologique
# (mÃªme filtrage que script 14)
# -------------------------------------------------------
annot_map <- read.table(ANNOT_FILE, sep = ";", header = TRUE,
                        stringsAsFactors = FALSE, strip.white = TRUE)
annot_lookup <- annot_map %>%
  filter(!is.na(annotation), annotation != "",
         !grepl("Contamination|TRM_CD8|Microglia_|Astrocyte|Oligodendrocyte",
                annotation)) %>%
  distinct(cluster, .keep_all = TRUE)
message("  Annotation map: ", nrow(annot_lookup), " clusters")

# -------------------------------------------------------
# Exclure les types non annotés / contaminants
# (on garde tout ce qui n'est PAS exclu, y compris T CD4 (Foxp3))
# -------------------------------------------------------
exclude_pattern <- "Contamination|TRM_CD8|Microglia_|Non annote|^$"
cell_types_all  <- unique(so$cell_type)
cell_types_all  <- cell_types_all[!is.na(cell_types_all) & cell_types_all != ""]
kept_types <- cell_types_all[!grepl(exclude_pattern, cell_types_all)]
so_sub <- subset(so, subset = cell_type %in% kept_types)
so_sub$cell_type <- factor(so_sub$cell_type)
Idents(so_sub) <- "cell_type"
message("  Clusters retenus: ", nlevels(so_sub$cell_type), " â€” ", ncol(so_sub), " cells")
print(table(so_sub$cell_type))

# -------------------------------------------------------
# Dotplot canonical markers (all conditions merged)
# -------------------------------------------------------
canonical_markers <- c(
  "Cd8a", "Cd3e", "Foxp3", "Ccr7",
  "H2-Aa", "H2-Ab1", "Ms4a7",
  "Lyz2", "C1qa", "P2ry12", "Hexb", "Tmem119",
  "Cd19", "Gzma", "Gzmb", "Ptprc"
)

expr_mat <- GetAssayData(so_sub, assay = "RNA", layer = "data")
markers_use <- canonical_markers[canonical_markers %in% rownames(expr_mat)]
message("  Canonical markers found: ", length(markers_use), " / ", length(canonical_markers))

domain_levels_plot <- order_annotations(levels(so_sub$cell_type), extended = TRUE)
domain_levels_plot <- domain_levels_plot[domain_levels_plot %in% levels(so_sub$cell_type)]
domain_vec <- factor(as.character(so_sub$cell_type), levels = domain_levels_plot)

dot_data <- do.call(rbind, lapply(markers_use, function(gene) {
  expr <- as.numeric(expr_mat[gene, ])
  do.call(rbind, lapply(domain_levels_plot, function(dom) {
    e <- expr[domain_vec == dom]
    data.frame(
      gene      = gene,
      domain    = dom,
      mean_expr = if (length(e) > 0) mean(e, na.rm = TRUE) else 0,
      pct_expr  = if (length(e) > 0) mean(e > 0, na.rm = TRUE) * 100 else 0,
      stringsAsFactors = FALSE
    )
  }))
}))

dot_data$gene   <- factor(dot_data$gene, levels = markers_use)
dot_data$domain <- factor(dot_data$domain, levels = domain_levels_plot)
cap_val <- quantile(dot_data$mean_expr, 0.99, na.rm = TRUE)
dot_data$mean_expr_capped <- pmin(dot_data$mean_expr, cap_val)

p_dot <- ggplot(dot_data, aes(x = gene, y = domain)) +
  geom_point(aes(size = pct_expr, color = mean_expr_capped)) +
  scale_color_gradient(low = "#ffee01", high = "#ff001e",
                       name = "Mean expr.\n(capped 99p)") +
  scale_size_continuous(range = c(0.3, 6), name = "% cells\nexpressing") +
  labs(
    title = "Immune niche canonical markers by cell type (res=0.3)",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y      = element_text(size = 9),
    plot.title       = element_text(face = "bold", size = 11, hjust = 0.5),
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
  )

fig_w_dot <- max(8, length(markers_use) * 0.40 + 3)
fig_h_dot <- max(5, length(domain_levels_plot) * 0.40 + 2)
dot_base <- file.path(OUT_DIR, "dotplot_canonical_markers_res03_OK_v2")

CairoPDF(paste0(dot_base, ".pdf"), width = fig_w_dot, height = fig_h_dot)
print(p_dot)
dev.off()

CairoJPEG(paste0(dot_base, ".jpg"),
          width = round(fig_w_dot * 150), height = round(fig_h_dot * 150),
          res = 150, quality = 95)
print(p_dot)
dev.off()

message("Saved: ", basename(paste0(dot_base, ".pdf")))
message("Saved: ", basename(paste0(dot_base, ".jpg")))

# -------------------------------------------------------
# FindAllMarkers
# -------------------------------------------------------
message("\nFindAllMarkers...")
markers <- FindAllMarkers(
  so_sub,
  only.pos        = TRUE,
  min.pct         = 0.10,
  logfc.threshold = 0.25,
  verbose         = FALSE
)
message("  ", nrow(markers), " DEGs")

# -------------------------------------------------------
# AverageExpression + heatmap top2
# -------------------------------------------------------
avg <- AverageExpression(
  so_sub,
  assays   = "RNA",
  layer    = "data",
  group.by = "cell_type"
)$RNA

# Ordre colonnes selon palette
cluster_order <- order_annotations(colnames(avg), extended = TRUE)
cluster_order <- cluster_order[cluster_order %in% colnames(avg)]
avg <- avg[, cluster_order, drop = FALSE]

for (top_n in c(5, 2)) {
  topN <- markers %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    filter(gene %in% rownames(avg),
           cluster %in% cluster_order) %>%
    mutate(cluster = factor(cluster, levels = cluster_order)) %>%
    arrange(cluster, desc(avg_log2FC), gene)

  if (nrow(topN) == 0) {
    message("Pas de gÃ¨nes pour top", top_n, " â€” skip"); next
  }

  heat_mat    <- avg[topN$gene, , drop = FALSE]
  row_labels  <- make.unique(as.character(topN$gene), sep = "__")
  rownames(heat_mat) <- row_labels

  block_sizes <- as.integer(table(topN$cluster)[levels(topN$cluster)])
  block_sizes <- block_sizes[!is.na(block_sizes)]
  gaps_row    <- if (length(block_sizes) > 1) cumsum(block_sizes)[-length(block_sizes)] else NULL

  out_base <- file.path(OUT_DIR,
    paste0("heatmap_immune_deg_top", top_n, "_res03"))

  save_heatmap_pair(
    heat_mat,
    title      = paste0("DEG immune sub-clusters top", top_n, " (res=0.3)"),
    out_base   = out_base,
    gaps_row   = gaps_row,
    labels_row = as.character(topN$gene),
    w = 12,
    h = max(10, 0.45 * nrow(heat_mat)),
    label_scale = ifelse(top_n == 2, 1.3, 1)
  )
}

message("\n=== Done ===")
