#!/usr/bin/env Rscript
# =============================================================
# Script: 28_dam_dotplot_microglia_conditions.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-05
#
# Goal:
#   Evaluate how DAM programs vary across conditions within the full
#   Microglia (P2ry12) population.
#   Dotplot: y = condition, x = selected DAM genes (per signature block).
#
# Strategy:
#   1. Load global BANKSY object (04_banksy_joint_lam08_after_bloc3.rds)
#   2. Reconstruct broad annotations from ncells_by_sample_lam02_res09_joint_long.csv
#   3. Subset to "Microglia (P2ry12)" cells
#   4. Group by condition (colData$sample)
#   5. Match curated DAM signatures to panel via case-harmonised matching
#   6. Rank genes by between-condition variance of avg expression; top 10 per sig
#   7. Compact dotplot across conditions
#
# Inputs:
#   objects/04_banksy_joint_lam08_after_bloc3.rds
#   ncells_by_sample_lam02_res09_joint_long.csv
#   outputs/banksy/dam_signature_curation/Upregulated_DAM.csv
#   outputs/banksy/dam_signature_curation/Downregulated_DAM.csv
#   outputs/banksy/dam_signature_curation/Microglia_signature_union.csv
#
# Outputs (outputs/banksy/dam_dotplot_microglia_conditions/):
#   filtered_Upregulated_DAM.csv
#   filtered_Downregulated_DAM.csv
#   filtered_Microglia_signature_union.csv
#   selected_genes_summary.csv
#   fig_dam_microglia_conditions.pdf / .jpg      [main figure]
#   fig_dam_upregulated_only.pdf / .jpg          [only if combined too crowded]
# =============================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(Banksy)
  library(tidyverse)
  library(ggplot2)
  library(Cairo)
})

base_path <- normalizePath(".")  # Run this script from the project root directory
setwd(base_path)

out_dir <- file.path("outputs", "banksy", "dam_dotplot_microglia_conditions")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================
# Parameters
# =============================================================

MICROGLIA_LABEL <- "Microglia (P2ry12)"
LAM             <- 0.2
RES_TARGET      <- 0.9
TOP_N           <- 10L
PCT_THRESHOLD   <- 0       # > 0 = any detection

SAMPLE_ORDER  <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi", "mock_6wpi")
SAMPLE_LABELS <- c(
  mock_6wpi = "Mock 6 wpi",
  LCMV_1wpi = "LCMV 1 wpi",
  LCMV_3wpi = "LCMV 3 wpi",
  LCMV_6wpi = "LCMV 6 wpi"
)

# =============================================================
# 1. Load global BANKSY object
# =============================================================

obj_file <- file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds")
stopifnot("Object file not found" = file.exists(obj_file))

message("Loading: ", obj_file)
se <- readRDS(obj_file)
message("  ", ncol(se), " cells | class: ", class(se)[1])

# =============================================================
# 2. Find BANKSY cluster column (lambda=0.2, res=0.9)
# =============================================================

find_cl_col <- function(se, lam, res) {
  all_cols <- clusterNames(se)
  lam_str  <- gsub("\\.", "\\\\.", as.character(lam))
  lam_cols <- all_cols[grep(paste0("lam", lam_str), all_cols)]
  if (length(lam_cols) == 0) return(NULL)
  res_vals <- suppressWarnings(as.numeric(sub(".*_res", "", lam_cols)))
  idx      <- which(!is.na(res_vals) & abs(res_vals - res) < 1e-8)
  if (length(idx) == 0) return(NULL)
  lam_cols[idx[1]]
}

cl_col <- find_cl_col(se, LAM, RES_TARGET)
if (is.null(cl_col)) {
  stop("BANKSY cluster column not found for lambda=", LAM,
       " res=", RES_TARGET,
       "\nAvailable: ", paste(clusterNames(se), collapse = ", "))
}
message("  Cluster column: ", cl_col)

# =============================================================
# 3. Reconstruct annotation mapping from CSV
# =============================================================

csv_path <- "ncells_by_sample_lam02_res09_joint_long.csv"
stopifnot("Annotation CSV not found" = file.exists(csv_path))

annot_long <- read_delim(
  csv_path,
  delim           = ";",
  locale          = locale(decimal_mark = "."),
  show_col_types  = FALSE,
  trim_ws         = TRUE
) %>%
  select(-matches("^Unnamed")) %>%
  mutate(
    banksy_domain = as.character(banksy_domain),
    annotation    = trimws(as.character(annotation))
  )

annotation_map <- annot_long %>%
  filter(!is.na(annotation), annotation != "") %>%
  distinct(banksy_domain, annotation)

message("  Annotation mappings: ", nrow(annotation_map))

# =============================================================
# 4. Assign annotation to each cell
# =============================================================

domain_labels <- paste0("Domain_", as.character(colData(se)[[cl_col]]))
anno_lookup   <- setNames(annotation_map$annotation,
                          annotation_map$banksy_domain)

annotation <- ifelse(
  !is.na(anno_lookup[domain_labels]) & anno_lookup[domain_labels] != "",
  anno_lookup[domain_labels],
  "Non annote"
)

cat("\nAnnotation distribution (all cells):\n")
print(sort(table(annotation), decreasing = TRUE))

if (!MICROGLIA_LABEL %in% annotation) {
  stop("Label '", MICROGLIA_LABEL, "' not found in annotations.\n",
       "Present: ", paste(sort(unique(annotation)), collapse = ", "))
}

# =============================================================
# 5. Subset to Microglia (P2ry12) cells
# =============================================================

mg_idx <- which(annotation == MICROGLIA_LABEL)
se_mg  <- se[, mg_idx]
message("\n  Microglia (P2ry12) cells: ", ncol(se_mg))

sample_vec <- as.character(colData(se_mg)$sample)
cat("\nCells per condition:\n")
print(table(sample_vec))

# Identify expression assay
assay_name <- if ("counts" %in% assayNames(se_mg)) "counts" else assayNames(se_mg)[1]
message("  Using assay: ", assay_name)

expr_mat  <- assay(se_mg, assay_name)    # genes × cells; log1p-like
panel_genes <- rownames(se_mg)

# Conditions present (keep SAMPLE_ORDER for those actually present)
cond_levels <- SAMPLE_ORDER[SAMPLE_ORDER %in% unique(sample_vec)]
cat("\nCondition order for y-axis: ", paste(cond_levels, collapse = " | "), "\n\n")

# =============================================================
# 6. Read curated signature files
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
# 7. Case-harmonised panel matching
#    Signature files use human ALLCAPS; panel uses mouse Titlecase.
#    Match via toTitleCase(tolower()); preserve original panel names.
# =============================================================

panel_title    <- tools::toTitleCase(tolower(panel_genes))
panel_name_map <- setNames(panel_genes, panel_title)   # Titlecase → original

harm_intersect <- function(sig_genes) {
  sig_title <- tools::toTitleCase(tolower(sig_genes))
  hits      <- sig_title[sig_title %in% panel_title]
  matched   <- panel_name_map[hits]
  unique(as.character(matched[!is.na(matched)]))
}

sig_filtered <- lapply(sig_raw, harm_intersect)

filtered_csv_names <- c(
  "Upregulated DAM"           = "filtered_Upregulated_DAM.csv",
  "Downregulated DAM"         = "filtered_Downregulated_DAM.csv",
  "Microglia signature union" = "filtered_Microglia_signature_union.csv"
)

for (nm in names(sig_filtered)) {
  message("  ", nm, ": ", length(sig_filtered[[nm]]),
          " genes after harmonised filter (from ", length(sig_raw[[nm]]), " raw)")
  write.csv(
    data.frame(gene = sig_filtered[[nm]]),
    file.path(out_dir, filtered_csv_names[nm]),
    row.names = FALSE
  )
}
message("  Saved: filtered_*.csv")

# =============================================================
# 8. Rank genes by between-condition variance of avg expression
#    → top 10 per signature for plotting
# =============================================================

# Per-gene × per-condition average expression
all_filt_genes <- unique(unlist(sig_filtered))
all_filt_genes <- all_filt_genes[all_filt_genes %in% rownames(expr_mat)]

cond_avg_mat <- do.call(cbind, lapply(cond_levels, function(cond) {
  idx  <- which(sample_vec == cond)
  if (length(idx) == 0) return(rep(NA_real_, length(all_filt_genes)))
  rowMeans(expr_mat[all_filt_genes, idx, drop = FALSE])
}))
rownames(cond_avg_mat) <- all_filt_genes
colnames(cond_avg_mat) <- cond_levels

gene_mean_expr <- rowMeans(expr_mat[all_filt_genes, , drop = FALSE])

# Between-condition variance of per-condition average expression
bcv <- apply(cond_avg_mat, 1, var, na.rm = TRUE)

# Build summary table
summary_rows <- lapply(names(sig_filtered), function(nm) {
  g_filt <- sig_filtered[[nm]]
  g_filt <- g_filt[g_filt %in% rownames(cond_avg_mat)]
  if (length(g_filt) == 0) return(NULL)
  df <- data.frame(
    signature              = nm,
    gene                   = g_filt,
    mean_expression        = round(gene_mean_expr[g_filt], 5),
    between_condition_var  = round(bcv[g_filt], 7),
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

cat("\nGenes selected per signature (ranked by between-condition variance):\n")
print(summary_rows %>%
        filter(selected_for_plot) %>%
        select(signature, gene, mean_expression, between_condition_var, rank_in_sig) %>%
        arrange(signature, rank_in_sig))

sig_selected <- lapply(names(sig_filtered), function(nm) {
  summary_rows %>%
    filter(signature == nm, selected_for_plot) %>%
    arrange(rank_in_sig) %>%
    pull(gene)
}) %>% setNames(names(sig_filtered))

# =============================================================
# 9. Compute dot statistics per gene × condition
# =============================================================

compute_dot_stats <- function(genes) {
  genes <- genes[genes %in% rownames(expr_mat)]
  if (length(genes) == 0) return(NULL)

  sub_mat <- expr_mat[genes, , drop = FALSE]

  rows <- lapply(genes, function(g) {
    g_expr <- sub_mat[g, ]
    lapply(cond_levels, function(cond) {
      idx  <- which(sample_vec == cond)
      if (length(idx) == 0) return(NULL)
      vals <- g_expr[idx]
      data.frame(
        gene        = g,
        condition   = cond,
        avg_expr    = mean(vals),
        pct_express = mean(vals > PCT_THRESHOLD) * 100,
        n_cells     = length(idx),
        stringsAsFactors = FALSE
      )
    }) %>% bind_rows()
  }) %>% bind_rows()

  # Z-score per gene across conditions
  rows %>%
    group_by(gene) %>%
    mutate(
      z_score = {
        mn  <- mean(avg_expr)
        sdd <- sd(avg_expr)
        if (!is.na(sdd) && sdd > 0) (avg_expr - mn) / sdd else rep(0, n())
      }
    ) %>%
    ungroup()
}

# =============================================================
# 10. Main figure: combined dotplot (all 3 signature blocks)
# =============================================================

block_order <- names(sig_selected)

# Build gene order with deduplication (first-block assignment wins)
sig_block_map_raw <- unlist(lapply(names(sig_selected), function(nm) {
  setNames(rep(nm, length(sig_selected[[nm]])), sig_selected[[nm]])
}))
all_genes_raw   <- names(sig_block_map_raw)
dedup_idx       <- !duplicated(all_genes_raw)
all_genes       <- all_genes_raw[dedup_idx]
sig_block_map   <- sig_block_map_raw[dedup_idx]

n_genes_total <- length(all_genes)
message("\nTotal genes for combined dotplot (after dedup): ", n_genes_total)

dot_df <- compute_dot_stats(all_genes)

dot_df$signature <- sig_block_map[dot_df$gene]
dot_df$gene      <- factor(dot_df$gene,      levels = all_genes)
dot_df$condition <- factor(dot_df$condition, levels = rev(cond_levels))
                                              # rev so mock is at bottom, LCMV 6 at top
dot_df$signature <- factor(dot_df$signature, levels = block_order)

# Condition display labels on y-axis
dot_df$condition_label <- factor(
  SAMPLE_LABELS[as.character(dot_df$condition)],
  levels = SAMPLE_LABELS[rev(cond_levels)]
)

z_cap   <- 2.5
z_range <- c(-z_cap, z_cap)

n_cond  <- length(cond_levels)   # y-axis rows
n_genes <- n_genes_total          # x-axis columns

fig_w <- max(8, n_genes * 0.55 + 3)
fig_h <- max(3.5, n_cond * 0.55 + 2)

p_combined <- ggplot(dot_df,
                     aes(x = gene, y = condition_label,
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
    range  = c(0.5, 7),
    limits = c(0, 100),
    name   = "% cells\nexpressing"
  ) +
  labs(
    title    = "DAM / Microglia signature — Microglia (P2ry12) across conditions",
    subtitle = paste0("n = ", ncol(se_mg), " microglia cells  |  ",
                      "top ", TOP_N, " genes per signature block",
                      " (ranked by between-condition variance)"),
    x = NULL,
    y = NULL
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
    plot.subtitle     = element_text(size = 8, colour = "grey40", hjust = 0),
    legend.position   = "right",
    legend.key.size   = unit(0.9, "lines"),
    legend.title      = element_text(size = 8.5),
    legend.text       = element_text(size = 8)
  )

CairoPDF(file.path(out_dir, "fig_dam_microglia_conditions.pdf"),
         width = fig_w, height = fig_h)
print(p_combined)
dev.off()
message("Saved: fig_dam_microglia_conditions.pdf  (",
        round(fig_w, 1), " × ", round(fig_h, 1), " in)")

CairoJPEG(file.path(out_dir, "fig_dam_microglia_conditions.jpg"),
           width = fig_w * 150, height = fig_h * 150, res = 150)
print(p_combined)
dev.off()
message("Saved: fig_dam_microglia_conditions.jpg")

# =============================================================
# 11. Optional secondary figure — Upregulated DAM only
#     Generated when combined figure has >= 25 genes (too wide)
# =============================================================

UP_ONLY_THRESHOLD <- 25

if (n_genes_total >= UP_ONLY_THRESHOLD) {
  message("\nCombined figure has ", n_genes_total,
          " genes — also generating Upregulated DAM-only figure.")

  up_genes <- sig_selected[["Upregulated DAM"]]

  if (length(up_genes) > 0) {
    dot_up <- compute_dot_stats(up_genes)
    dot_up$gene            <- factor(dot_up$gene, levels = up_genes)
    dot_up$condition       <- factor(dot_up$condition, levels = rev(cond_levels))
    dot_up$condition_label <- factor(
      SAMPLE_LABELS[as.character(dot_up$condition)],
      levels = SAMPLE_LABELS[rev(cond_levels)]
    )

    fig_w2 <- max(6, length(up_genes) * 0.65 + 3)
    fig_h2 <- fig_h

    p_up <- ggplot(dot_up,
                   aes(x = gene, y = condition_label,
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
        range  = c(0.5, 7),
        limits = c(0, 100),
        name   = "% cells\nexpressing"
      ) +
      labs(
        title    = "Upregulated DAM — Microglia (P2ry12) across conditions",
        subtitle = paste0("n = ", ncol(se_mg), " microglia cells  |  ",
                          "top ", TOP_N, " genes by between-condition variance"),
        x = NULL,
        y = NULL
      ) +
      theme_bw(base_size = 10) +
      theme(
        axis.text.x      = element_text(angle = 45, hjust = 1, size = 8,
                                        face = "italic"),
        axis.text.y      = element_text(size = 9),
        panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold", size = 11, hjust = 0),
        plot.subtitle    = element_text(size = 8, colour = "grey40", hjust = 0),
        legend.position  = "right",
        legend.key.size  = unit(0.9, "lines"),
        legend.title     = element_text(size = 8.5),
        legend.text      = element_text(size = 8)
      )

    CairoPDF(file.path(out_dir, "fig_dam_upregulated_only.pdf"),
             width = fig_w2, height = fig_h2)
    print(p_up)
    dev.off()

    CairoJPEG(file.path(out_dir, "fig_dam_upregulated_only.jpg"),
              width = fig_w2 * 150, height = fig_h2 * 150, res = 150)
    print(p_up)
    dev.off()
    message("Saved: fig_dam_upregulated_only.pdf/jpg  (",
            round(fig_w2, 1), " × ", round(fig_h2, 1), " in)")

  } else {
    message("  No Upregulated DAM genes passed panel filter — skipping secondary figure.")
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
cat("  fig_dam_microglia_conditions.pdf / .jpg\n")
if (n_genes_total >= UP_ONLY_THRESHOLD)
  cat("  fig_dam_upregulated_only.pdf / .jpg\n")
