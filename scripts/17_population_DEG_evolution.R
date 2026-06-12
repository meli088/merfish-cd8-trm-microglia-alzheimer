#!/usr/bin/env Rscript
# =============================================================
# Script: 17_population_DEG_evolution.R
# Description:
#   DEGs 6wpi vs mock/1wpi/3wpi sur TOUTE la population
#   (pas par domaine BANKSY).
#   Produit :
#     - CSV DEGs par comparaison
#     - CSV summary (log2FC x 3 comparaisons)
#     - Heatmap z-score (union top N gènes)
#
# Usage:
#   Rscript scripts/17_population_DEG_evolution.R --folder microglia
#   Rscript scripts/17_population_DEG_evolution.R --folder immune_acod1 --top 15
# =============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(tidyverse)
  library(pheatmap)
})

base_path <- normalizePath(".")  # Run this script from the project root directory
setwd(base_path)

# -------------------------------------------------------
# Default parameters
# -------------------------------------------------------
IMMUNE_LABEL <- "Immune (Acod1)"
FOLDER_NAME  <- NULL
TOP_N_UNION  <- 10
FC_THRESH    <- 0.05
MIN_PCT      <- 0.10

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  for (i in seq_along(args)) {
    a <- args[i]
    if (a %in% c("--label", "-l") && (i + 1) <= length(args)) IMMUNE_LABEL <- args[i + 1]
    if (a == "--folder"           && (i + 1) <= length(args)) FOLDER_NAME  <- args[i + 1]
    if (a == "--top"              && (i + 1) <= length(args)) TOP_N_UNION  <- as.integer(args[i + 1])
  }
}

slugify   <- function(x) gsub("^_+|_+$", "", gsub("[^a-z0-9]+", "_", tolower(trimws(x))))
safe_label <- if (!is.null(FOLDER_NAME)) tolower(trimws(FOLDER_NAME)) else slugify(IMMUNE_LABEL)
out_dir   <- file.path("outputs", "banksy", safe_label, "population_DEG_evolution")
csv_dir   <- file.path(out_dir, "DEG_csv")
fig_dir   <- file.path(out_dir, "heatmaps")
for (d in c(csv_dir, fig_dir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

message("=== 17_population_DEG_evolution.R | ", safe_label, " ===")

# -------------------------------------------------------
# Load BANKSY object
# -------------------------------------------------------
cands <- c(
  file.path("objects", paste0("07_immune_banksy_lam02_", safe_label, ".rds")),
  file.path("objects", "07_immune_banksy_lam02_after_bloc1.rds")
)
banksy_obj <- cands[file.exists(cands)][1]
if (is.na(banksy_obj)) stop("BANKSY object not found for: ", safe_label)
se <- readRDS(banksy_obj)
message("  ", ncol(se), " cells loaded from: ", banksy_obj)

assay_name    <- if ("counts" %in% assayNames(se)) "counts" else assayNames(se)[1]
so            <- as.Seurat(se, counts = assay_name, data = NULL)
default_assay <- if ("RNA" %in% SeuratObject::Assays(so)) "RNA" else SeuratObject::Assays(so)[1]
DefaultAssay(so) <- default_assay

needs_norm <- tryCatch({
  dm <- GetAssayData(so, assay = default_assay, layer = "data")
  nrow(dm) == 0 || ncol(dm) == 0
}, error = function(e) TRUE)
if (needs_norm) {
  so <- NormalizeData(so, assay = default_assay, verbose = FALSE)
  message("  NormalizeData done")
}

# -------------------------------------------------------
# Sample info
# -------------------------------------------------------
sample_vec    <- as.character(so@meta.data$sample)
avail_samples <- unique(sample_vec)
if (!"LCMV_6wpi" %in% avail_samples) stop("LCMV_6wpi not found in object")

sample_order <- c("mock_6wpi", "LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")
present      <- sample_order[sample_order %in% avail_samples]
comparisons  <- setdiff(present, "LCMV_6wpi")
comp_tags    <- gsub("LCMV_", "", comparisons)
comp_tags    <- gsub("_", " ", comp_tags)
comp_slugs   <- gsub("[^a-z0-9]", "_", tolower(comp_tags))

message("  Samples: ", paste(present, collapse = ", "))
message("  Comparisons: 6wpi vs ", paste(comparisons, collapse = ", "))

n_6wpi <- sum(sample_vec == "LCMV_6wpi")
message("  LCMV_6wpi cells: ", n_6wpi)
if (n_6wpi < 1) stop("No LCMV_6wpi cells")

# Set ident = sample
so <- SetIdent(so, value = "sample")

# -------------------------------------------------------
# FindMarkers per comparison
# -------------------------------------------------------
all_comp_deg <- list()

for (i in seq_along(comparisons)) {
  comp      <- comparisons[i]
  comp_slug <- comp_slugs[i]
  comp_tag  <- comp_tags[i]

  n_comp <- sum(sample_vec == comp)
  if (n_comp < 10) {
    message("  Skip vs ", comp, ": only ", n_comp, " cells")
    next
  }

  mk <- tryCatch(
    FindMarkers(so, ident.1 = "LCMV_6wpi", ident.2 = comp,
                only.pos = FALSE, min.pct = MIN_PCT,
                logfc.threshold = FC_THRESH, return.thresh = 1,
                test.use = "wilcox", verbose = FALSE),
    error = function(e) { message("  ERROR vs ", comp, ": ", conditionMessage(e)); NULL }
  )
  if (is.null(mk) || nrow(mk) == 0) {
    message("  No markers vs ", comp)
    next
  }

  mk$gene       <- rownames(mk)
  mk$direction  <- ifelse(mk$avg_log2FC > 0, "up_6wpi", "down_6wpi")
  mk$comparison <- paste0("6wpi_vs_", comp_slug)
  mk <- mk %>%
    filter(p_val_adj < 0.05) %>%
    arrange(desc(avg_log2FC)) %>%
    select(gene, avg_log2FC, direction, pct.1, pct.2, p_val, p_val_adj, comparison)

  all_comp_deg[[comp]] <- mk
  n_up   <- sum(mk$avg_log2FC > 0)
  n_down <- sum(mk$avg_log2FC < 0)
  message("  vs ", comp_tag, ": ", nrow(mk), " DEGs (", n_up, " up / ", n_down, " down)")

  fname <- file.path(csv_dir, paste0("DEG_", safe_label, "_6wpi_vs_", comp_slug, ".csv"))
  tryCatch(
    write.csv(mk, fname, row.names = FALSE),
    error = function(e) message("  WARNING write: ", conditionMessage(e))
  )
}

if (length(all_comp_deg) == 0) stop("No significant DEGs for any comparison")

# -------------------------------------------------------
# Heatmap: union top N genes per comparison
# -------------------------------------------------------
all_deg_df <- do.call(rbind, all_comp_deg)

top_genes_union <- all_deg_df %>%
  group_by(comparison) %>%
  slice_max(abs(avg_log2FC), n = TOP_N_UNION, with_ties = FALSE) %>%
  ungroup() %>%
  pull(gene) %>%
  unique()

message("  Union top genes for heatmap: ", length(top_genes_union))
if (length(top_genes_union) < 2) stop("Too few genes for heatmap")

# AverageExpression per sample
avg_exp <- tryCatch(
  AverageExpression(so, assays = default_assay, features = top_genes_union,
                    group.by = "sample", verbose = FALSE)[[default_assay]],
  error = function(e) { stop("AverageExpression ERROR: ", conditionMessage(e)) }
)

# Fix underscore -> dash replacement by Seurat
colnames(avg_exp) <- gsub("-", "_", colnames(avg_exp))

# All 4 timepoints in order; fill absent ones with a zero column
cols_use <- present
for (col in cols_use[!cols_use %in% colnames(avg_exp)]) {
  zero_col <- matrix(0, nrow = nrow(avg_exp), ncol = 1,
                     dimnames = list(rownames(avg_exp), col))
  avg_exp <- cbind(avg_exp, zero_col)
}
genes_ok <- top_genes_union[top_genes_union %in% rownames(avg_exp)]
avg_mat  <- as.matrix(avg_exp[genes_ok, cols_use, drop = FALSE])

if (nrow(avg_mat) < 2 || ncol(avg_mat) < 2) stop("Matrix too small")

# Rename columns for readability
col_rename <- c(mock_6wpi = "Mock", LCMV_1wpi = "1 wpi", LCMV_3wpi = "3 wpi", LCMV_6wpi = "6 wpi")
colnames(avg_mat) <- col_rename[colnames(avg_mat)]

# Z-score by gene (row)
mat_z <- t(scale(t(avg_mat)))
mat_z[is.nan(mat_z)] <- 0

# Row direction (kept for summary CSV only)
row_dir <- setNames(sapply(rownames(mat_z), function(g) {
  n_up   <- sum(sapply(all_comp_deg, function(df) g %in% df$gene[df$direction == "up_6wpi"]))
  n_down <- sum(sapply(all_comp_deg, function(df) g %in% df$gene[df$direction == "down_6wpi"]))
  if (n_up > n_down) "up_6wpi"
  else if (n_down > n_up) "down_6wpi"
  else {
    fc_vals <- sapply(all_comp_deg, function(df) {
      idx <- which(df$gene == g)
      if (length(idx) == 0) NA_real_ else df$avg_log2FC[idx[1]]
    })
    fc_vals <- fc_vals[!is.na(fc_vals)]
    if (length(fc_vals) == 0 || fc_vals[which.max(abs(fc_vals))] < 0) "down_6wpi" else "up_6wpi"
  }
}), rownames(mat_z))

# Order: up block then down block
gene_order <- c(
  names(row_dir)[row_dir == "up_6wpi"],
  names(row_dir)[row_dir == "down_6wpi"]
)
mat_z_ord <- mat_z[gene_order, , drop = FALSE]

# -------------------------------------------------------
# Summary CSV
# -------------------------------------------------------
summary_wide <- tibble(gene = gene_order, direction_majority = row_dir[gene_order])
for (j in seq_along(comparisons)) {
  cj <- comparisons[j]; sj <- comp_slugs[j]
  if (cj %in% names(all_comp_deg)) {
    df_j <- all_comp_deg[[cj]] %>% select(gene, avg_log2FC, p_val_adj)
    colnames(df_j)[2:3] <- c(paste0("log2FC_vs_", sj), paste0("padj_vs_", sj))
    summary_wide <- left_join(summary_wide, df_j, by = "gene")
  } else {
    summary_wide[[paste0("log2FC_vs_", sj)]] <- NA_real_
    summary_wide[[paste0("padj_vs_", sj)]]   <- NA_real_
  }
}
padj_cols <- grep("^padj_", colnames(summary_wide), value = TRUE)
summary_wide$sig_in_comparisons <- apply(
  summary_wide[, padj_cols, drop = FALSE], 1,
  function(x) {
    idx <- which(!is.na(x) & x < 0.05)
    if (length(idx) == 0) "none" else paste(comp_slugs[idx], collapse = "|")
  }
)
write.csv(summary_wide,
          file.path(csv_dir, paste0("summary_genes_", safe_label, "_population.csv")),
          row.names = FALSE)
message("  Summary table: ", nrow(summary_wide), " genes saved")

# -------------------------------------------------------
# Heatmap
# -------------------------------------------------------
p_heat <- pheatmap::pheatmap(
  mat_z_ord,
  color          = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  cluster_rows   = TRUE,
  cluster_cols   = FALSE,
  show_colnames  = TRUE,
  show_rownames  = TRUE,
  fontsize_row   = 8,
  fontsize_col   = 10,
  border_color   = "grey90",
  main           = paste0(safe_label, " — gene expression dynamics over time"),
  silent         = TRUE
)

fig_h      <- max(5, nrow(mat_z_ord) * 0.35 + 3)
fname_base <- file.path(fig_dir, paste0("heatmap_", safe_label, "_population"))
ggsave(paste0(fname_base, ".pdf"), p_heat$gtable, width = 7, height = fig_h)
ggsave(paste0(fname_base, ".jpg"), p_heat$gtable, width = 7, height = fig_h, dpi = 300)
message("  Heatmap saved: ", basename(fname_base))

message("\nDone. Outputs in: ", out_dir)
