#!/usr/bin/env Rscript
# =============================================================
# Script: 14_immune_analysis.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Description:
#   For each BANKSY resolution (0.1 – 0.5):
#     1. FindAllMarkers (Wilcoxon) → full DEG CSV + top-10 summary CSV
#        (top-10 used for manual domain annotation)
#     2. UMAP 2×2 grid — one panel per sample (mock, 1wpi, 3wpi, 6wpi),
#        domains colored by number, same palette across all panels,
#        no text labels on points
#
# Input:  objects/07_immune_banksy_lam02_<label>.rds
# Output: outputs/banksy/<label>/analysis/
#           DEGs/DEGs_all_res0X.csv
#           DEGs/DEGs_top10_per_domain_res0X.csv
#           figures/umap_grid_res0X.jpg  /  .pdf
# =============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(tidyverse)
  library(patchwork)
})

base_path <- normalizePath(".")  # Run this script from the project root directory
setwd(base_path)

# Palette et ordre des annotations (fichier partagé)
source("scripts/00_palette.R")

# -------------------------------------------------------
# Parameters
# -------------------------------------------------------
IMMUNE_LABEL  <- "Immune (Acod1)"
FOLDER_NAME   <- NULL   # optional: override folder/file slug with --folder microglia
TOP_N_DEG     <- 5   # top N markers per domain for dotplot / summary CSV
N_HEAT_GENES  <- 5    # top N per domain for heatmap rows
DOTPLOT_ONLY  <- FALSE  # set TRUE via --dotplot-only to skip DEGs/UMAP/composition

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  for (i in seq_along(args)) {
    a <- args[i]
    if (a %in% c("--label", "-l") && (i + 1) <= length(args)) {
      IMMUNE_LABEL <- args[i + 1]
    }
    if (a == "--folder" && (i + 1) <= length(args)) {
      FOLDER_NAME <- args[i + 1]
    }
    if (a == "--dotplot-only") {
      DOTPLOT_ONLY <- TRUE
    }
  }
}

slugify_label <- function(label) {
  label <- tolower(trimws(label))
  label <- gsub("[^a-z0-9]+", "_", label)
  gsub("^_+|_+$", "", label)
}

safe_label <- if (!is.null(FOLDER_NAME)) tolower(trimws(FOLDER_NAME)) else slugify_label(IMMUNE_LABEL)
label_root <- file.path("outputs", "banksy", safe_label)
out_dir    <- file.path(label_root, "analysis")
fig_dir    <- file.path(out_dir, "figures")
deg_dir    <- file.path(out_dir, "DEGs")

for (d in c(out_dir, fig_dir, deg_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# -------------------------------------------------------
# Inputs
# -------------------------------------------------------
banksy_candidates <- c(
  file.path("objects", paste0("07_immune_banksy_lam02_", safe_label, ".rds")),
  file.path("objects", "07_immune_banksy_lam02_after_bloc1.rds")
)
banksy_obj <- banksy_candidates[file.exists(banksy_candidates)][1]
if (is.na(banksy_obj)) stop("BANKSY object not found for label: ", IMMUNE_LABEL)

# -------------------------------------------------------
# Helpers
# -------------------------------------------------------
safe_res_name <- function(x) gsub("\\.", "", sprintf("%.1f", x))

find_cl_col <- function(se_obj, lam, res) {
  cd      <- SummarizedExperiment::colData(se_obj)
  cn      <- colnames(cd)
  res_txt <- sprintf("%.1f", res)
  lam_pat <- paste0("lam", gsub("\\.", "\\\\.", as.character(lam)))
  res_pat <- paste0("_res", gsub("\\.", "\\\\.", res_txt), "$")
  cols    <- cn[grepl(lam_pat, cn) & grepl(res_pat, cn)]
  if (length(cols) == 0) return(NULL)
  cols[1]
}

palette20 <- c(
  "#ff0004", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
  "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
  "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A"
)
# Si un mapping annotation->domain est fourni, on construit une palette
# biologiquement cohérente depuis GLOBAL_PALETTE.
build_dom_colors <- function(domain_levels, annot_lookup = NULL) {
  cols <- setNames(
    rep(palette20, length.out = length(domain_levels)),
    domain_levels
  )
  if (!is.null(annot_lookup)) {
    for (dom in domain_levels) {
      bio <- annot_lookup[[dom]]
      if (!is.null(bio) && bio %in% names(GLOBAL_PALETTE)) {
        cols[[dom]] <- GLOBAL_PALETTE[[bio]]
      }
    }
  }
  cols
}

# Palette pour labels biologiques directs (cell_type names)
build_ct_colors <- function(ct_levels) {
  cols <- setNames(rep(palette20, length.out = length(ct_levels)), ct_levels)
  for (ct in ct_levels) {
    if (ct %in% names(GLOBAL_PALETTE)) cols[[ct]] <- GLOBAL_PALETTE[[ct]]
  }
  cols
}
sample_order <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")

# -------------------------------------------------------
# Load BANKSY object
# -------------------------------------------------------
message("Loading: ", banksy_obj)
se         <- readRDS(banksy_obj)
# Exclude mock_6wpi — analysis focuses on LCMV time course
se         <- se[, as.character(SummarizedExperiment::colData(se)$sample) != "mock_6wpi"]
sample_vec <- as.character(SummarizedExperiment::colData(se)$sample)
message("  ", ncol(se), " cells | ", length(unique(sample_vec)), " samples (mock_6wpi excluded)")

# Charger le mapping Domain_X -> nom biologique (produit par 17_annotate_immune_object.R)
annot_csv_path <- file.path("outputs", "banksy", safe_label, "annotation_map_immune_res03.csv")
annot_lookup   <- NULL
if (file.exists(annot_csv_path)) {
  annot_df     <- read.csv(annot_csv_path, sep = ";", stringsAsFactors = FALSE, check.names = FALSE)
  annot_lookup <- setNames(trimws(annot_df$annotation), trimws(annot_df$cluster))
  annot_lookup <- annot_lookup[nchar(annot_lookup) > 0]
  message("  Annotation map loaded: ", length(annot_lookup), " entries from ", basename(annot_csv_path))
} else {
  message("  No annotation map found — figures will show Domain_X labels")
}

# -------------------------------------------------------
# Convert to Seurat
# -------------------------------------------------------
assay_name   <- if ("counts" %in% assayNames(se)) "counts" else assayNames(se)[1]
so           <- as.Seurat(se, counts = assay_name, data = NULL)
default_assay <- if ("RNA" %in% SeuratObject::Assays(so)) "RNA" else SeuratObject::Assays(so)[1]
DefaultAssay(so) <- default_assay

# Normalize if needed
needs_norm <- tryCatch({
  dm <- GetAssayData(so, assay = default_assay, layer = "data")
  nrow(dm) == 0 || ncol(dm) == 0
}, error = function(e) TRUE)
if (needs_norm) {
  so <- NormalizeData(so, assay = default_assay, verbose = FALSE)
  message("  NormalizeData done")
}

# UMAP — re-use UMAP from script 13 if stored as a reducedDim,
# otherwise compute from the Harmony embedding.
if (!"umap" %in% names(so@reductions)) {
  if ("UMAP" %in% names(so@reductions)) {
    so[["umap"]] <- so[["UMAP"]]
  } else {
    avail     <- names(so@reductions)
    harm_red  <- avail[grepl("^[Hh]armony", avail)][1]
    pca_red   <- avail[grepl("^PCA",        avail)][1]
    umap_src  <- if (!is.na(harm_red)) harm_red else if (!is.na(pca_red)) pca_red else avail[1]
    n_dims    <- min(40L, ncol(so@reductions[[umap_src]]))
    so <- RunUMAP(so, reduction = umap_src, dims = seq_len(n_dims), verbose = FALSE)
    message("  RunUMAP on '", umap_src, "' (", n_dims, " dims)")
  }
}

# Pre-extract UMAP coordinates once (shared across all resolutions)
umap_mat <- as.data.frame(so@reductions[["umap"]]@cell.embeddings)
colnames(umap_mat) <- c("UMAP1", "UMAP2")
umap_mat$sample <- so@meta.data$sample

RES_SEQ <- seq(0.1, 0.5, by = 0.1)

# Canonical cell-type markers for dot plot
canonical_markers <- c(
 "Cd3e", "Cd4", "Cd8a", "Cd8b", "Trbc1", "Klrg1",  # T/NK
 "Cd3g", "Itgam", "Itgax", "Ptprc", "Aqp4", "Map2", "Ly6c1", "Hexb", "Tmem119", "P2ry12", "Ccr7", "H2-Aa", "H2-Ab1" 
)

# -------------------------------------------------------
# Per-resolution loop
# -------------------------------------------------------
for (res in RES_SEQ) {
  cl_col <- find_cl_col(se, 0.2, res)
  if (is.null(cl_col)) {
    message("  res=", res, " — no cluster column, skipping")
    next
  }

  res_id      <- safe_res_name(res)
  cluster_var <- paste0("banksy_res", res_id)
  message("\n=== res=", res, " (", cl_col, ") ===")

  # Labels techniques Domain_X (conservés pour DEGs)
  domain_labels <- paste0("Domain_", as.character(colData(se)[[cl_col]]))
  raw_levels    <- paste0("Domain_", sort(unique(as.integer(colData(se)[[cl_col]]))))
  so@meta.data[[cluster_var]] <- factor(domain_labels, levels = raw_levels)
  so <- SetIdent(so, value = cluster_var)
  domain_levels <- raw_levels

  # Labels biologiques pour les figures (Domain_X -> cell_type si mapping dispo)
  display_var <- paste0(cluster_var, "_display")
  if (!is.null(annot_lookup)) {
    mapped     <- annot_lookup[raw_levels]
    disp_lvls  <- ifelse(!is.na(mapped) & nchar(mapped) > 0, mapped, raw_levels)
    bio_lvls   <- unique(disp_lvls[!grepl("^Domain_", disp_lvls)])
    dom_lvls   <- unique(disp_lvls[ grepl("^Domain_", disp_lvls)])
    disp_order <- c(order_annotations(bio_lvls), dom_lvls)
    cell_disp  <- annot_lookup[domain_labels]
    cell_disp  <- ifelse(!is.na(cell_disp) & nchar(cell_disp) > 0, cell_disp, domain_labels)
  } else {
    disp_order <- raw_levels
    cell_disp  <- domain_labels
  }
  so@meta.data[[display_var]] <- factor(cell_disp, levels = disp_order)
  dom_colors <- build_ct_colors(disp_order)

  # ── 1. DEGs ───────────────────────────────────────────────
  if (DOTPLOT_ONLY) {
    message("  --dotplot-only: skipping DEGs / UMAP / composition")
  } else {
  message("  FindAllMarkers...")
  markers <- tryCatch(
    FindAllMarkers(
      so,
      only.pos        = FALSE,
      min.pct         = 0.10,
      logfc.threshold = 0.05,
      return.thresh   = 1,
      test.use        = "wilcox",
      verbose         = FALSE
    ),
    error = function(e) {
      message("  FindAllMarkers ERROR: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(markers) || nrow(markers) == 0) {
    message("  No markers found — skipping")
    next
  }

  markers <- markers %>% arrange(cluster, desc(avg_log2FC))

  # Full table
  tryCatch(
    write.csv(markers, file.path(deg_dir, paste0("DEGs_all_res", res_id, ".csv")), row.names = FALSE),
    error = function(e) message("  WARNING: could not write DEGs_all_res", res_id, ".csv (", conditionMessage(e), ") — file may be open in Excel")
  )
  message("  DEGs: ", nrow(markers), " rows")

  # Top-N summary per domain (for manual annotation)
  top_n <- markers %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = TOP_N_DEG, with_ties = FALSE) %>%
    ungroup()
  tryCatch(
    write.csv(top_n, file.path(deg_dir, paste0("DEGs_top", TOP_N_DEG, "_per_domain_res", res_id, ".csv")), row.names = FALSE),
    error = function(e) message("  WARNING: could not write DEGs_top summary (", conditionMessage(e), ")")
  )
  message("  Top-", TOP_N_DEG, " summary: ",
          n_distinct(top_n$cluster), " domains × ", TOP_N_DEG, " genes")

  # ── 2. UMAP 2×2 grid ──────────────────────────────────────
  present <- sample_order[sample_order %in% unique(so@meta.data$sample)]

  udf <- umap_mat %>%
    mutate(
      cluster = so@meta.data[[display_var]],
      sample  = factor(sample, levels = present)
    )

  p_grid <- ggplot(udf, aes(UMAP1, UMAP2, color = cluster)) +
    geom_point(size = 0.8, alpha = 0.7, stroke = 0) +
    scale_color_manual(values = dom_colors, drop = FALSE) +
    facet_wrap(~sample, ncol = 2) +
    labs(
      title = paste0(IMMUNE_LABEL, " — BANKSY sub-clusters (res=", res, ")"),
      color = "Cell type"
    ) +
    guides(color = guide_legend(
      override.aes = list(size = 3, alpha = 1),
      ncol = 1
    )) +
    theme_classic(base_size = 10) +
    theme(
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = 10),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      legend.position  = "right",
      plot.title       = element_text(face = "bold", size = 11)
    )

  ggsave(
    file.path(fig_dir, paste0("umap_grid_res", res_id, ".pdf")),
    p_grid, width = 10, height = 8
  )
  ggsave(
    file.path(fig_dir, paste0("umap_grid_res", res_id, ".jpg")),
    p_grid, width = 10, height = 8, dpi = 300
  )
  message("  UMAP 2×2 grid saved")

  # ── 3. Composition — 100% stacked bar (all samples) ─────────
  comp <- so@meta.data %>%
    as.data.frame() %>%
    transmute(
      sample  = as.character(sample),
      cluster = .data[[display_var]]
    ) %>%
    count(sample, cluster, name = "n_cells") %>%
    tidyr::complete(sample, cluster, fill = list(n_cells = 0L)) %>%
    group_by(sample) %>%
    mutate(pct = round(100 * n_cells / sum(n_cells), 2)) %>%
    ungroup() %>%
    mutate(sample = factor(sample, levels = present))

  # Save global composition CSV
  tryCatch(
    write.csv(comp, file.path(deg_dir, paste0("composition_res", res_id, ".csv")), row.names = FALSE),
    error = function(e) message("  WARNING: could not write composition CSV (", conditionMessage(e), ")")
  )

  # 100% stacked bar — all samples in one figure
  p_comp <- ggplot(comp, aes(x = sample, y = pct, fill = cluster)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    scale_fill_manual(values = dom_colors) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      x     = NULL,
      y     = "% of cells",
      fill  = "Cell type",
      title = paste0(IMMUNE_LABEL, " — Cell type composition (res=", res, ")")
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(
    file.path(fig_dir, paste0("composition_stacked_res", res_id, ".pdf")),
    p_comp, width = 7, height = 5
  )
  ggsave(
    file.path(fig_dir, paste0("composition_stacked_res", res_id, ".jpg")),
    p_comp, width = 7, height = 5, dpi = 300
  )
  message("  Composition stacked bar saved")

  # ── 4. Domain composition over time (line plot) ───────────
  p_time <- comp %>%
    ggplot(aes(x = sample, y = pct, color = cluster, group = cluster)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = dom_colors) +
    labs(
      x     = NULL,
      y     = "% of cells",
      color = "Cell type",
      title = paste0(IMMUNE_LABEL, " — Cell type composition over time (res=", res, ")")
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(
    file.path(fig_dir, paste0("composition_overtime_res", res_id, ".pdf")),
    p_time, width = 8, height = 5
  )
  ggsave(
    file.path(fig_dir, paste0("composition_overtime_res", res_id, ".jpg")),
    p_time, width = 8, height = 5, dpi = 300
  )
  message("  Composition over time saved")
  } # end if (!DOTPLOT_ONLY)

  # ── 5. Canonical markers dot plot ─────────────────────────
  expr_mat    <- GetAssayData(so, assay = default_assay, layer = "data")
  markers_use <- canonical_markers[canonical_markers %in% rownames(expr_mat)]
  message("  Canonical markers found: ", length(markers_use), " / ", length(canonical_markers))

  if (length(markers_use) > 0) {
    domain_vec        <- as.character(so@meta.data[[display_var]])
    domain_levels_plot <- levels(so@meta.data[[display_var]])

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

    dot_data$gene   <- factor(dot_data$gene,   levels = markers_use)
    dot_data$domain <- factor(dot_data$domain, levels = domain_levels_plot)

    cap_val <- quantile(dot_data$mean_expr, 0.99, na.rm = TRUE)
    dot_data$mean_expr_capped <- pmin(dot_data$mean_expr, cap_val)

    p_dot <- ggplot(dot_data, aes(x = gene, y = domain)) +
      geom_point(aes(size = pct_expr, color = mean_expr_capped)) +
      scale_color_gradient(low = "#ffee01", high = "#ff001e",
                           name = "Mean expr.\n(capped 99p)") +
      scale_size_continuous(range = c(0.3, 6), name = "% cells\nexpressing") +
      labs(
        title = paste0(IMMUNE_LABEL, " — Canonical markers by cell type (res=", res, ")"),
        x = NULL, y = NULL
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

    n_dom  <- length(domain_levels)
    n_gene <- length(markers_use)
    fig_w  <- max(8, n_gene * 0.40 + 3)
    fig_h  <- max(4, n_dom  * 0.40 + 2)

    ggsave(
      file.path(fig_dir, paste0("dotplot_canonical_markers_res", res_id, ".pdf")),
      p_dot, width = fig_w, height = fig_h
    )
    ggsave(
      file.path(fig_dir, paste0("dotplot_canonical_markers_res", res_id, ".jpg")),
      p_dot, width = fig_w, height = fig_h, dpi = 300
    )
    message("  Canonical markers dot plot saved")
  }
}

message("\nDone. Outputs in: ", out_dir)
