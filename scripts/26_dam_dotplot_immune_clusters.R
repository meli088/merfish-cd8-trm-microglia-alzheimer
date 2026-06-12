#!/usr/bin/env Rscript
# =============================================================
# Script: 26_dam_dotplot_immune_clusters.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-05
#
# Goal:
#   Compact DAM signature dotplot across immune niche clusters
#   (colData$cell_type from object 08), using curated signature
#   gene lists from script 25.
#
# Inputs:
#   objects/08_immune_annotated_lam02_res03.rds
#   outputs/banksy/dam_signature_curation/Upregulated_DAM.csv
#   outputs/banksy/dam_signature_curation/Downregulated_DAM.csv
#   outputs/banksy/dam_signature_curation/Microglia_signature_union.csv
#
# Outputs (outputs/banksy/dam_dotplot_immune_clusters/):
#   filtered_Upregulated_DAM.csv
#   filtered_Downregulated_DAM.csv
#   filtered_Microglia_signature_union.csv
#   selected_genes_summary.csv
#   fig_dam_dotplot_combined.pdf / .jpg      [main figure]
#   fig_dam_dotplot_microglia_only.pdf / .jpg [only if combined is too crowded]
#
# Notes:
#   - Expression values in 'counts' assay are log1p-like (already normalised)
#   - Dot size  = % cells expressing (> 0) in each cluster
#   - Dot color = z-score of average expression across clusters (per gene)
#   - Genes matched to panel via case-harmonised matching (toTitleCase/tolower)
#     because signature files use human ALLCAPS while panel uses mouse Titlecase
#   - Genes ranked by mean expression in full object; top 10 kept per sig
#   - No DEG, no heatmaps, no between-condition analysis here
# =============================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(tidyverse)
  library(ggplot2)
  library(Cairo)
})

base_path <- normalizePath(".")  # Run this script from the project root directory
setwd(base_path)

out_dir <- file.path("outputs", "banksy", "dam_dotplot_immune_clusters")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

TOP_N         <- 10L    # max genes per signature block
PCT_THRESHOLD <-  0     # expression threshold (> 0 = any detection)

# =============================================================
# 1. Load annotated immune niche object
# =============================================================

message("Loading: objects/08_immune_annotated_lam02_res03.rds")
se <- readRDS(file.path("objects", "08_immune_annotated_lam02_res03.rds"))
message("  ", ncol(se), " cells | ", nrow(se), " genes")

# Identify expression assay
assay_name <- if ("counts" %in% assayNames(se)) "counts" else assayNames(se)[1]
message("  Using assay: ", assay_name)

expr_mat <- assay(se, assay_name)   # genes × cells; log1p-like

cell_type_vec <- as.character(colData(se)$cell_type)
panel_genes   <- rownames(se)

cat("\ncell_type labels (cells per cluster):\n")
ct_tbl <- sort(table(cell_type_vec), decreasing = TRUE)
print(ct_tbl)

# Preferred cell-type order — pattern-based (grepl) so labels like
# "Microglia (C1qa)" or "T CD8 (Gzmb)" are correctly placed.
# Patterns are matched in priority order; first match wins.
CELL_PATTERN_ORDER <- c(
  "Microglia", "Homeostatic.Mg", "Homeostatic.micro",
  "DAM",
  "Reactive.Mg", "Reactive.micro",
  "IFN.responsive.Mg", "IFN.responsive.micro", "IFN",
  "Mono",
  "Mac",
  "DC",
  "NK", "ILC",
  "T cells", "Treg", "T cell",
  "B cell", "Plasma"
)

build_ct_order <- function(labels) {
  assigned <- rep(NA_integer_, length(labels))
  for (i in seq_along(CELL_PATTERN_ORDER)) {
    pat  <- CELL_PATTERN_ORDER[i]
    hits <- grepl(pat, labels, ignore.case = TRUE, perl = TRUE)
    assigned[hits & is.na(assigned)] <- i   # first pattern wins
  }
  # Labels matched by a pattern: order by pattern priority, then label alpha
  matched  <- labels[!is.na(assigned)]
  matched  <- matched[order(assigned[!is.na(assigned)], matched)]
  # Remaining labels: alphabetical
  rest     <- sort(labels[is.na(assigned)])
  c(matched, rest)
}

ct_levels <- build_ct_order(unique(cell_type_vec))
cat("\nCluster order for y-axis:\n")
cat(paste(ct_levels, collapse = "\n"), "\n\n")

# =============================================================
# 2. Read curated signature files
# =============================================================

sig_dir <- file.path("outputs", "banksy", "dam_signature_curation")

read_sig <- function(fname) {
  f <- file.path(sig_dir, fname)
  stopifnot("Signature file not found" = file.exists(f))
  read.csv(f, stringsAsFactors = FALSE)$gene
}

sig_raw <- list(
  "Upregulated DAM"           = read_sig("Upregulated_DAM.csv"),
  "Downregulated DAM"         = read_sig("Downregulated_DAM.csv"),
  "Microglia signature union" = read_sig("Microglia_signature_union.csv")
)

for (nm in names(sig_raw)) {
  message("  ", nm, ": ", length(sig_raw[[nm]]), " genes (raw)")
}

# =============================================================
# 3. Filter each signature to genes present in the MERFISH panel
#    Uses case-harmonised matching: signature genes may be in human
#    ALLCAPS (e.g. TREM2) while the panel uses mouse Titlecase (Trem2).
#    Strategy: toTitleCase(tolower()) both sides, then look up original
#    panel gene names — so downstream outputs always use panel names.
# =============================================================

# Build Titlecase index of panel genes
panel_title    <- tools::toTitleCase(tolower(panel_genes))
panel_name_map <- setNames(panel_genes, panel_title)   # Titlecase → original

harm_intersect <- function(sig_genes) {
  sig_title <- tools::toTitleCase(tolower(sig_genes))
  hits      <- sig_title[sig_title %in% panel_title]
  matched   <- panel_name_map[hits]
  matched   <- matched[!is.na(matched)]
  unique(as.character(matched))   # original panel gene names, no duplicates
}

sig_filtered <- lapply(sig_raw, harm_intersect)

for (nm in names(sig_filtered)) {
  message("  ", nm, ": ", length(sig_filtered[[nm]]),
          " genes after case-harmonised panel filter (from ",
          length(sig_raw[[nm]]), " raw)")
}

# Save filtered CSVs
safe_name <- function(nm) gsub("[^A-Za-z0-9]", "_", nm)
filtered_csv_names <- c(
  "Upregulated DAM"           = "filtered_Upregulated_DAM.csv",
  "Downregulated DAM"         = "filtered_Downregulated_DAM.csv",
  "Microglia signature union" = "filtered_Microglia_signature_union.csv"
)

for (nm in names(sig_filtered)) {
  write.csv(
    data.frame(gene = sig_filtered[[nm]]),
    file.path(out_dir, filtered_csv_names[nm]),
    row.names = FALSE
  )
  message("  Saved: ", filtered_csv_names[nm])
}

# =============================================================
# 4. Rank genes by between-cluster variance of avg expression → top 10
#    Rationale: selects genes most discriminative across clusters,
#    rather than simply the most globally abundant ones.
# =============================================================

# Pre-compute per-gene × per-cluster average expression
cluster_avg <- function(genes_vec) {
  g_ok <- genes_vec[genes_vec %in% rownames(expr_mat)]
  if (length(g_ok) == 0) return(matrix(nrow = 0, ncol = 0))
  sub <- expr_mat[g_ok, , drop = FALSE]
  do.call(cbind, lapply(ct_levels, function(ct) {
    idx <- which(cell_type_vec == ct)
    if (length(idx) == 0) rep(0, length(g_ok))
    else rowMeans(sub[, idx, drop = FALSE])
  }))
}

# All filtered genes across all sigs (union) — compute avg once
all_filt_genes <- unique(unlist(sig_filtered))
avg_mat        <- cluster_avg(all_filt_genes)   # genes × clusters
rownames(avg_mat) <- all_filt_genes
gene_mean_expr    <- rowMeans(expr_mat[all_filt_genes, , drop = FALSE])

# Between-cluster variance of per-cluster average expression
bcv <- apply(avg_mat, 1, var)   # one value per gene

# Build summary table with all filtered genes
summary_rows <- lapply(names(sig_filtered), function(nm) {
  g_filt <- sig_filtered[[nm]]
  if (length(g_filt) == 0) return(NULL)
  g_filt <- g_filt[g_filt %in% rownames(avg_mat)]
  df <- data.frame(
    signature              = nm,
    gene                   = g_filt,
    mean_expression        = round(gene_mean_expr[g_filt], 5),
    between_cluster_var    = round(bcv[g_filt], 7),
    rank_in_sig            = rank(-bcv[g_filt], ties.method = "first"),
    stringsAsFactors = FALSE
  )
  df$selected_for_plot <- df$rank_in_sig <= TOP_N
  df
}) %>% bind_rows()

write.csv(summary_rows,
          file.path(out_dir, "selected_genes_summary.csv"),
          row.names = FALSE)
message("\nSaved: selected_genes_summary.csv")

cat("\nGenes selected per signature (ranked by between-cluster variance):\n")
print(summary_rows %>%
        filter(selected_for_plot) %>%
        select(signature, gene, mean_expression, between_cluster_var, rank_in_sig) %>%
        arrange(signature, rank_in_sig))

# Build per-signature selected lists (highest between-cluster variance first)
sig_selected <- lapply(names(sig_filtered), function(nm) {
  summary_rows %>%
    filter(signature == nm, selected_for_plot) %>%
    arrange(rank_in_sig) %>%
    pull(gene)
}) %>% setNames(names(sig_filtered))

# =============================================================
# 5. Compute dot statistics per gene × cell_type
# =============================================================
# dot_size  = % cells with expr > 0 in the cluster
# dot_color = z-score of mean expression across clusters (per gene)

compute_dot_stats <- function(genes, expr_matrix, cluster_vec, cluster_levels) {

  if (length(genes) == 0) return(NULL)

  genes <- genes[genes %in% rownames(expr_matrix)]
  if (length(genes) == 0) return(NULL)

  sub_mat <- expr_matrix[genes, , drop = FALSE]

  rows <- lapply(genes, function(g) {
    g_expr <- sub_mat[g, ]
    lapply(cluster_levels, function(ct) {
      idx <- which(cluster_vec == ct)
      if (length(idx) == 0) return(NULL)
      vals <- g_expr[idx]
      data.frame(
        gene        = g,
        cell_type   = ct,
        avg_expr    = mean(vals),
        pct_express = mean(vals > PCT_THRESHOLD) * 100,
        n_cells     = length(idx),
        stringsAsFactors = FALSE
      )
    }) %>% bind_rows()
  }) %>% bind_rows()

  # Z-score per gene across clusters (for colour scale)
  rows <- rows %>%
    group_by(gene) %>%
    mutate(
      z_score = {
        mn  <- mean(avg_expr)
        sdd <- sd(avg_expr)
        if (!is.na(sdd) && sdd > 0) (avg_expr - mn) / sdd else rep(0, n())
      }
    ) %>%
    ungroup()

  rows
}

# =============================================================
# 6. Main figure: combined dotplot (all 3 signature blocks)
# =============================================================

# Flatten all selected genes (in block order)
# A gene may appear in more than one signature after case-harmonised matching;
# keep first occurrence only (first-block assignment) to avoid duplicate factor levels.
block_order <- names(sig_selected)

sig_block_map_raw <- unlist(lapply(names(sig_selected), function(nm) {
  setNames(rep(nm, length(sig_selected[[nm]])), sig_selected[[nm]])
}))

all_genes_ordered_raw <- names(sig_block_map_raw)          # preserves block order
dedup_idx             <- !duplicated(all_genes_ordered_raw) # first occurrence wins
all_genes_ordered     <- all_genes_ordered_raw[dedup_idx]
sig_block_map         <- sig_block_map_raw[dedup_idx]

n_genes_total <- length(all_genes_ordered)
message("\nTotal genes for combined dotplot (after dedup): ", n_genes_total)

dot_df <- compute_dot_stats(all_genes_ordered, expr_mat,
                            cell_type_vec, ct_levels)

dot_df$signature <- sig_block_map[dot_df$gene]
dot_df$gene      <- factor(dot_df$gene, levels = all_genes_ordered)
dot_df$cell_type <- factor(dot_df$cell_type, levels = ct_levels)
dot_df$signature <- factor(dot_df$signature, levels = block_order)

# Colour scale limits (symmetric z-score, capped)
z_cap   <- 2.5
z_range <- c(-z_cap, z_cap)

# Decide figure size
n_ct    <- length(ct_levels)   # y-axis clusters
n_genes <- n_genes_total       # x-axis genes

fig_w <- max(8, n_genes * 0.55 + 3)
fig_h <- max(5, n_ct   * 0.40 + 2.5)

# Publication-style combined dotplot with facet strips per signature block
p_combined <- ggplot(dot_df,
                     aes(x = gene, y = cell_type,
                         size = pct_express, colour = z_score)) +
  geom_point() +
  facet_grid(. ~ signature,
             scales = "free_x",
             space  = "free_x") +
  scale_colour_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    limits   = z_range,
    oob      = scales::squish,
    name     = "Scaled\nexpression\n(z-score)"
  ) +
  scale_size_continuous(
    range  = c(0.3, 6),
    limits = c(0, 100),
    name   = "% cells\nexpressing"
  ) +
  labs(
    title = "DAM / Microglia signature — immune niche clusters",
    x     = NULL,
    y     = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 8,
                                     face = "italic"),
    axis.text.y       = element_text(size = 9),
    strip.text        = element_text(size = 8.5, face = "bold"),
    strip.background  = element_rect(fill = "#f0f0f0", colour = "grey60"),
    panel.grid.major  = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(face = "bold", size = 11, hjust = 0),
    legend.position   = "right",
    legend.key.size   = unit(0.9, "lines"),
    legend.title      = element_text(size = 8.5),
    legend.text       = element_text(size = 8)
  )

pdf_path <- file.path(out_dir, "fig_dam_dotplot_combined.pdf")
jpg_path <- file.path(out_dir, "fig_dam_dotplot_combined.jpg")

CairoPDF(pdf_path, width = fig_w, height = fig_h)
print(p_combined)
dev.off()
message("Saved: fig_dam_dotplot_combined.pdf  (", round(fig_w, 1),
        " × ", round(fig_h, 1), " in)")

CairoJPEG(jpg_path, width = fig_w * 150, height = fig_h * 150, res = 150)
print(p_combined)
dev.off()
message("Saved: fig_dam_dotplot_combined.jpg")

# =============================================================
# 7. Optional secondary figure — Microglia signature only
#    Generated only when the combined figure has >= 25 genes
#    (combined figure likely too wide to read comfortably)
# =============================================================

MG_ONLY_THRESHOLD <- 25   # generate secondary figure if combined has this many genes

if (n_genes_total >= MG_ONLY_THRESHOLD) {
  message("\nCombined figure has ", n_genes_total,
          " genes — also generating Microglia signature-only figure.")

  mg_genes <- sig_selected[["Microglia signature union"]]

  if (length(mg_genes) > 0) {
    dot_mg <- compute_dot_stats(mg_genes, expr_mat, cell_type_vec, ct_levels)
    dot_mg$gene      <- factor(dot_mg$gene, levels = mg_genes)
    dot_mg$cell_type <- factor(dot_mg$cell_type, levels = ct_levels)

    fig_w2 <- max(6, length(mg_genes) * 0.65 + 3)
    fig_h2 <- fig_h   # same height as combined

    p_mg <- ggplot(dot_mg,
                   aes(x = gene, y = cell_type,
                       size = pct_express, colour = z_score)) +
      geom_point() +
      scale_colour_gradient2(
        low      = "#2166AC",
        mid      = "white",
        high     = "#B2182B",
        midpoint = 0,
        limits   = z_range,
        oob      = scales::squish,
        name     = "Scaled\nexpression\n(z-score)"
      ) +
      scale_size_continuous(
        range  = c(0.3, 6),
        limits = c(0, 100),
        name   = "% cells\nexpressing"
      ) +
      labs(
        title = "Microglia signature (union) — immune niche clusters",
        x     = NULL,
        y     = NULL
      ) +
      theme_bw(base_size = 10) +
      theme(
        axis.text.x      = element_text(angle = 45, hjust = 1, size = 8,
                                        face = "italic"),
        axis.text.y      = element_text(size = 9),
        panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold", size = 11, hjust = 0),
        legend.position  = "right",
        legend.key.size  = unit(0.9, "lines"),
        legend.title     = element_text(size = 8.5),
        legend.text      = element_text(size = 8)
      )

    pdf2 <- file.path(out_dir, "fig_dam_dotplot_microglia_only.pdf")
    jpg2 <- file.path(out_dir, "fig_dam_dotplot_microglia_only.jpg")

    CairoPDF(pdf2, width = fig_w2, height = fig_h2)
    print(p_mg)
    dev.off()
    message("Saved: fig_dam_dotplot_microglia_only.pdf  (",
            round(fig_w2, 1), " × ", round(fig_h2, 1), " in)")

    CairoJPEG(jpg2, width = fig_w2 * 150, height = fig_h2 * 150, res = 150)
    print(p_mg)
    dev.off()
    message("Saved: fig_dam_dotplot_microglia_only.jpg")

  } else {
    message("  No Microglia signature genes passed panel filter — skipping secondary figure.")
  }

} else {
  message("\nCombined figure has only ", n_genes_total,
          " genes — secondary figure not needed.")
}

# =============================================================
# Done
# =============================================================

message("\n=== Done. Outputs in: ", out_dir, " ===\n")
cat("Files written:\n")
cat("  filtered_Upregulated_DAM.csv\n")
cat("  filtered_Downregulated_DAM.csv\n")
cat("  filtered_Microglia_signature_union.csv\n")
cat("  selected_genes_summary.csv\n")
cat("  fig_dam_dotplot_combined.pdf / .jpg\n")
if (n_genes_total >= MG_ONLY_THRESHOLD)
  cat("  fig_dam_dotplot_microglia_only.pdf / .jpg\n")
