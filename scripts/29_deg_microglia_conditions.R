#!/usr/bin/env Rscript
# =============================================================
# Script: 29_deg_microglia_conditions.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-05
#
# Goal:
#   DEG analysis within the global Microglia (P2ry12) population,
#   comparing each LCMV time-point to mock.
#   Produces one DEG CSV + one volcano plot per contrast,
#   plus a compact summary table.
#
# Contrasts:
#   LCMV_1wpi vs mock_6wpi
#   LCMV_3wpi vs mock_6wpi
#   LCMV_6wpi vs mock_6wpi
#
# Inputs:
#   objects/04_banksy_joint_lam08_after_bloc3.rds
#   ncells_by_sample_lam02_res09_joint_long.csv
#
# Outputs (outputs/banksy/deg_microglia_conditions/):
#   DEG_LCMV_1wpi_vs_mock.csv
#   DEG_LCMV_3wpi_vs_mock.csv
#   DEG_LCMV_6wpi_vs_mock.csv
#   DEG_summary.csv
#   fig_volcano_LCMV_1wpi_vs_mock.pdf / .jpg
#   fig_volcano_LCMV_3wpi_vs_mock.pdf / .jpg
#   fig_volcano_LCMV_6wpi_vs_mock.pdf / .jpg
# =============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(Banksy)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(Cairo)
})

base_path <- normalizePath(".")  # Run this script from the project root directory
setwd(base_path)

out_dir <- file.path("outputs", "banksy", "deg_microglia_conditions")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================
# Parameters
# =============================================================

MICROGLIA_LABEL <- "Microglia (P2ry12)"
LAM             <- 0.2
RES_TARGET      <- 0.9
REF_CONDITION   <- "mock_6wpi"
SAMPLE_ORDER    <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi", "mock_6wpi")
COMPARISONS     <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")  # each vs REF_CONDITION

# FindMarkers settings
MIN_PCT      <- 0.05    # low threshold — MERFISH panels are small
FC_THRESH    <- 0.1     # log2FC threshold for gene screening
TEST_USE     <- "wilcox"

# Volcano display
FDR_CUTOFF   <- 0.05
FC_CUTOFF    <- 0.25    # |log2FC| threshold for "significant" labelling
N_LABEL      <- 15      # max gene labels per volcano

SAMPLE_LABELS <- c(
  mock_6wpi = "Mock 6 wpi",
  LCMV_1wpi = "LCMV 1 wpi",
  LCMV_3wpi = "LCMV 3 wpi",
  LCMV_6wpi = "LCMV 6 wpi"
)

COND_COLORS <- c(
  LCMV_1wpi = "#56B4E9",
  LCMV_3wpi = "#E69F00",
  LCMV_6wpi = "#D55E00"
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

# =============================================================
# 4. Assign annotation and subset to Microglia (P2ry12)
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
  stop("Label '", MICROGLIA_LABEL, "' not found.\n",
       "Present: ", paste(sort(unique(annotation)), collapse = ", "))
}

mg_idx <- which(annotation == MICROGLIA_LABEL)
se_mg  <- se[, mg_idx]
message("\n  Microglia (P2ry12) cells: ", ncol(se_mg))

sample_vec <- as.character(colData(se_mg)$sample)
cat("\nCells per condition:\n")
print(table(sample_vec))

# Check reference is present
if (!REF_CONDITION %in% sample_vec) {
  stop("Reference condition '", REF_CONDITION, "' not found in microglia cells.")
}

# Keep only comparisons with enough cells
valid_comps <- COMPARISONS[sapply(COMPARISONS, function(c) sum(sample_vec == c) >= 10)]
if (length(valid_comps) == 0) stop("No comparison has >= 10 microglia cells.")
message("  Valid comparisons: ",
        paste(paste0(valid_comps, " vs ", REF_CONDITION), collapse = ", "))

# =============================================================
# 5. Convert to Seurat for FindMarkers
# =============================================================

assay_name    <- if ("counts" %in% assayNames(se_mg)) "counts" else assayNames(se_mg)[1]
# suppressWarnings: as.Seurat() renames reducedDim keys (PC→PC_, UMAP_*→...)
# to comply with Seurat key rules — cosmetic only, does not affect FindMarkers
so            <- suppressWarnings(as.Seurat(se_mg, counts = assay_name, data = NULL))
default_assay <- if ("RNA" %in% SeuratObject::Assays(so)) "RNA" else SeuratObject::Assays(so)[1]
DefaultAssay(so) <- default_assay

# Normalize if data layer is empty
needs_norm <- tryCatch({
  dm <- GetAssayData(so, assay = default_assay, layer = "data")
  nrow(dm) == 0 || ncol(dm) == 0
}, error = function(e) TRUE)
if (needs_norm) {
  so <- NormalizeData(so, assay = default_assay, verbose = FALSE)
  message("  NormalizeData done")
}

# Set identity = sample
so <- SetIdent(so, value = "sample")

# =============================================================
# 6. DEG per contrast via FindMarkers (Wilcoxon)
# =============================================================

all_deg  <- list()
slug_map <- setNames(
  gsub("[^a-z0-9]", "_", tolower(paste0(valid_comps, "_vs_mock"))),
  valid_comps
)

for (comp in valid_comps) {
  slug <- slug_map[[comp]]
  n1   <- sum(sample_vec == comp)
  n2   <- sum(sample_vec == REF_CONDITION)
  message("\n--- ", comp, " vs ", REF_CONDITION,
          "  (n=", n1, " vs n=", n2, ") ---")

  mk <- tryCatch(
    FindMarkers(
      so,
      ident.1          = comp,
      ident.2          = REF_CONDITION,
      only.pos         = FALSE,
      min.pct          = MIN_PCT,
      logfc.threshold  = FC_THRESH,
      return.thresh    = 1,
      test.use         = TEST_USE,
      verbose          = FALSE
    ),
    error = function(e) {
      message("  ERROR: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(mk) || nrow(mk) == 0) {
    message("  No results returned — skipping.")
    next
  }

  # Add columns
  mk$gene       <- rownames(mk)
  mk$contrast   <- paste0(comp, "_vs_", REF_CONDITION)
  mk$direction  <- case_when(
    mk$p_val_adj <= FDR_CUTOFF & mk$avg_log2FC >  FC_CUTOFF ~ "up",
    mk$p_val_adj <= FDR_CUTOFF & mk$avg_log2FC < -FC_CUTOFF ~ "down",
    TRUE ~ "ns"
  )
  mk <- mk %>%
    arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
    select(gene, avg_log2FC, pct.1, pct.2,
           p_val, p_val_adj, direction, contrast)

  all_deg[[comp]] <- mk

  n_up   <- sum(mk$direction == "up")
  n_down <- sum(mk$direction == "down")
  message("  ", nrow(mk), " genes tested | ",
          n_up, " up | ", n_down, " down (FDR<=", FDR_CUTOFF,
          ", |log2FC|>", FC_CUTOFF, ")")

  csv_path_out <- file.path(out_dir, paste0("DEG_", slug, ".csv"))
  write.csv(mk, csv_path_out, row.names = FALSE)
  message("  Saved: DEG_", slug, ".csv")
}

if (length(all_deg) == 0) stop("No DEG results for any contrast.")

# =============================================================
# 7. DEG summary table
# =============================================================

summary_df <- lapply(names(all_deg), function(comp) {
  mk   <- all_deg[[comp]]
  slug <- slug_map[[comp]]
  data.frame(
    contrast                     = paste0(comp, "_vs_", REF_CONDITION),
    n_genes_tested               = nrow(mk),
    n_up_FDR_0.05                = sum(mk$p_val_adj <= 0.05 & mk$avg_log2FC > 0),
    n_down_FDR_0.05              = sum(mk$p_val_adj <= 0.05 & mk$avg_log2FC < 0),
    n_up_FDR_0.05_logFC_0.25    = sum(mk$p_val_adj <= 0.05 & mk$avg_log2FC >  FC_CUTOFF),
    n_down_FDR_0.05_logFC_0.25  = sum(mk$p_val_adj <= 0.05 & mk$avg_log2FC < -FC_CUTOFF),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

cat("\n=== DEG Summary ===\n")
print(summary_df, row.names = FALSE)
write.csv(summary_df, file.path(out_dir, "DEG_summary.csv"), row.names = FALSE)
message("Saved: DEG_summary.csv")

# =============================================================
# 7b. Detection context — QC table + interpretation note + boxplot
# =============================================================

# Per-cell: n genes detected (count > 0) and total counts
cnt_mg           <- assay(se_mg, assay_name)
ngenes_per_cell  <- colSums(cnt_mg > 0)
ncounts_per_cell <- colSums(cnt_mg)

detect_per_cond <- tibble(
  sample   = factor(sample_vec, levels = SAMPLE_ORDER),
  n_genes  = as.numeric(ngenes_per_cell),
  n_counts = as.numeric(ncounts_per_cell)
) %>%
  group_by(sample) %>%
  summarise(
    n_cells           = n(),
    mean_genes_det    = round(mean(n_genes),   1),
    median_genes_det  = median(n_genes),
    mean_counts_det   = round(mean(n_counts),  1),
    median_counts_det = median(n_counts),
    .groups = "drop"
  ) %>%
  mutate(condition_label = SAMPLE_LABELS[as.character(sample)]) %>%
  select(sample, condition_label, n_cells,
         mean_genes_det, median_genes_det,
         mean_counts_det, median_counts_det)

# Merge DEG counts (ref condition gets NAs)
deg_per_cond <- summary_df %>%
  mutate(sample = sub("_vs_.*", "", contrast)) %>%
  select(sample, n_genes_tested,
         n_up_FDR_0.05, n_down_FDR_0.05,
         n_up_FDR_0.05_logFC_0.25, n_down_FDR_0.05_logFC_0.25)

context_df <- detect_per_cond %>%
  left_join(deg_per_cond, by = "sample") %>%
  arrange(factor(as.character(sample), levels = SAMPLE_ORDER))

cat("\n=== Detection context ===\n")
print(as.data.frame(context_df), row.names = FALSE, na.print = "-")

write.csv(context_df, file.path(out_dir, "DEG_detection_context.csv"), row.names = FALSE)
message("Saved: DEG_detection_context.csv")

# --- Interpretation note ---
note_lines <- c(
  "DEG Interpretation Note — Microglia (P2ry12)",
  "============================================",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "1. Detection differences across conditions",
  "   The number of genes detected per microglia cell varies across conditions",
  "   (see DEG_detection_context.csv, columns mean_genes_det / median_genes_det).",
  "   This reflects both biological changes (transcriptional activation) and",
  "   technical variation inherent to spatial transcriptomics (MERFISH panel).",
  "",
  "2. Asymmetry in DEG results",
  "   A strong asymmetry (many upregulated genes, few or no downregulated genes)",
  "   can partly reflect higher detection / transcriptome complexity in the",
  "   LCMV condition relative to mock, rather than purely biological upregulation.",
  "   Genes that appear 'upregulated' may include transcripts that are simply",
  "   more frequently captured in conditions with higher overall detection.",
  "",
  "3. LCMV_3wpi vs mock — special caution",
  "   If LCMV_3wpi shows substantially higher mean/median genes detected than",
  "   mock_6wpi, upregulated DEGs from that contrast should be interpreted with",
  "   particular caution: the detection advantage could inflate the apparent",
  "   number of upregulated genes.",
  "",
  "4. Recommended interpretation approach",
  "   - Focus on genes with large |log2FC| (> 0.5) in addition to FDR significance.",
  "   - Cross-reference with the DAM signature dotplot (script 28) which uses",
  "     z-scored average expression (less sensitive to raw detection rates).",
  "   - Treat the direction 'up in LCMV' as 'higher expression OR higher detection'",
  "     unless independently validated.",
  "",
  "Reference: DEG_detection_context.csv, DEG_summary.csv"
)
writeLines(note_lines, file.path(out_dir, "DEG_interpretation_note.txt"))
message("Saved: DEG_interpretation_note.txt")

# --- Detection boxplot (violin + boxplot) ---
detect_cells_df <- tibble(
  sample  = factor(sample_vec, levels = SAMPLE_ORDER),
  n_genes = as.numeric(ngenes_per_cell)
) %>%
  mutate(
    condition_label = factor(
      SAMPLE_LABELS[as.character(sample)],
      levels = SAMPLE_LABELS[SAMPLE_ORDER]
    )
  )

cond_palette <- c(
  "Mock 6 wpi"  = "grey60",
  "LCMV 1 wpi"  = "#56B4E9",
  "LCMV 3 wpi"  = "#E69F00",
  "LCMV 6 wpi"  = "#D55E00"
)

p_det <- ggplot(detect_cells_df,
                aes(x = condition_label, y = n_genes, fill = condition_label)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.45, colour = NA) +
  geom_boxplot(width = 0.18, outlier.size = 0.5, outlier.alpha = 0.4,
               colour = "grey20", fill = "white") +
  scale_fill_manual(values = cond_palette, guide = "none") +
  labs(
    title = "Microglia (P2ry12) \u2014 detected genes per cell",
    x     = NULL,
    y     = "N genes detected (count > 0)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title         = element_text(face = "bold", size = 10, hjust = 0),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank()
  )

CairoPDF(file.path(out_dir, "fig_detection_boxplot_microglia.pdf"), width = 5, height = 4)
print(p_det)
dev.off()
CairoJPEG(file.path(out_dir, "fig_detection_boxplot_microglia.jpg"),
          width = 5 * 150, height = 4 * 150, res = 150)
print(p_det)
dev.off()
message("Saved: fig_detection_boxplot_microglia.pdf / .jpg")

# =============================================================
# 8. Volcano plots — one per contrast
# =============================================================

direction_colors <- c(
  "up"   = "#B2182B",
  "down" = "#2166AC",
  "ns"   = "grey75"
)

make_volcano <- function(mk, comp) {
  slug  <- slug_map[[comp]]
  label <- SAMPLE_LABELS[comp]
  ref_l <- SAMPLE_LABELS[REF_CONDITION]

  # -log10(p_val_adj), floor at machine minimum to avoid -Inf
  mk$neg_log10_fdr <- -log10(pmax(mk$p_val_adj, .Machine$double.eps))

  # Cap extreme values for display (keeps axis readable)
  x_cap <- max(3.5, quantile(abs(mk$avg_log2FC), 0.999, na.rm = TRUE) * 1.05)
  y_cap <- max(5,   quantile(mk$neg_log10_fdr,   0.999, na.rm = TRUE) * 1.05)

  mk$x_plot <- pmin(pmax(mk$avg_log2FC,   -x_cap), x_cap)
  mk$y_plot <- pmin(mk$neg_log10_fdr, y_cap)

  # Genes to label: most significant among up/down, prefer highest |log2FC|
  sig_mk <- mk %>% filter(direction != "ns")
  if (nrow(sig_mk) > 0) {
    lab_genes <- sig_mk %>%
      arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
      slice_head(n = N_LABEL) %>%
      pull(gene)
  } else {
    # Fallback: label the most significant genes regardless of FC
    lab_genes <- mk %>%
      arrange(p_val_adj) %>%
      slice_head(n = min(N_LABEL, 10)) %>%
      pull(gene)
  }

  mk$label <- ifelse(mk$gene %in% lab_genes, mk$gene, NA_character_)

  n_up   <- sum(mk$direction == "up",   na.rm = TRUE)
  n_down <- sum(mk$direction == "down", na.rm = TRUE)

  p <- ggplot(mk, aes(x = x_plot, y = y_plot,
                      colour = direction, label = label)) +
    geom_point(size = 1.2, alpha = 0.7, stroke = 0) +
    geom_text_repel(
      data        = mk %>% filter(!is.na(label)),
      aes(label   = label),
      size        = 2.6,
      fontface    = "italic",
      max.overlaps = 20,
      box.padding  = 0.3,
      point.padding = 0.2,
      segment.size  = 0.3,
      segment.color = "grey50",
      show.legend   = FALSE
    ) +
    geom_hline(yintercept = -log10(FDR_CUTOFF),
               linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF),
               linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    scale_colour_manual(
      values = direction_colors,
      labels = c(up   = paste0("Up in ", label, " (n=", n_up, ")"),
                 down = paste0("Down in ", label, " (n=", n_down, ")"),
                 ns   = "Not significant"),
      name   = NULL
    ) +
    scale_x_continuous(
      limits = c(-x_cap, x_cap),
      expand = expansion(mult = 0.02)
    ) +
    scale_y_continuous(
      limits = c(0, y_cap),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      title    = paste0("Microglia (P2ry12): ", label, " vs ", ref_l),
      subtitle = paste0(nrow(mk), " genes tested  |  FDR \u2264 ", FDR_CUTOFF,
                        "  |  |log\u2082FC| > ", FC_CUTOFF),
      x        = paste0("log\u2082FC (", label, " / ", ref_l, ")"),
      y        = expression(-log[10](FDR))
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title      = element_text(face = "bold", size = 11, hjust = 0),
      plot.subtitle   = element_text(size = 8, colour = "grey40", hjust = 0),
      legend.position = "bottom",
      legend.text     = element_text(size = 8),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey93", linewidth = 0.3)
    )

  # Save
  fname <- paste0("fig_volcano_", slug)
  CairoPDF(file.path(out_dir, paste0(fname, ".pdf")), width = 6, height = 5.5)
  print(p)
  dev.off()
  CairoJPEG(file.path(out_dir, paste0(fname, ".jpg")),
            width = 6 * 150, height = 5.5 * 150, res = 150)
  print(p)
  dev.off()
  message("  Saved: ", fname, ".pdf / .jpg")
}

message("\nGenerating volcano plots...")
for (comp in names(all_deg)) {
  make_volcano(all_deg[[comp]], comp)
}

# =============================================================
# Done
# =============================================================

message("\n=== Done. Outputs in: ", out_dir, " ===\n")
cat("Files written:\n")
for (comp in names(all_deg)) {
  slug <- slug_map[[comp]]
  cat("  DEG_", slug, ".csv\n", sep = "")
  cat("  fig_volcano_", slug, ".pdf / .jpg\n", sep = "")
}
cat("  DEG_summary.csv\n")
cat("  DEG_detection_context.csv\n")
cat("  DEG_interpretation_note.txt\n")
cat("  fig_detection_boxplot_microglia.pdf / .jpg\n")
