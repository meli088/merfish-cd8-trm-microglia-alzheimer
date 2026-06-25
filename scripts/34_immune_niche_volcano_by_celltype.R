#!/usr/bin/env Rscript
# =============================================================
# Script: 34_immune_niche_volcano_by_celltype.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-05
#
# Goal:
#   For each annotated cell type in the immune sub-cluster object
#   (08_immune_annotated_lam02_res03.rds), compare gene expression
#   between temporal conditions and produce one volcano plot per
#   cell type × contrast.
#
# Contrast (for each cell type):
#   LCMV_6wpi vs LCMV_1wpi
#
# Inputs:
#   objects/08_immune_annotated_lam02_res03.rds
#
# Outputs (outputs/banksy/immune_niche_volcano_by_celltype/):
#   [Sec 2] DEG_[celltype_slug]_[contrast_slug].csv
#           fig_volcano_[celltype_slug]_[contrast_slug].pdf/.jpg
#   [Sec 3] fig_ndeg_summary_by_celltype.pdf/.jpg
#           ndeg_summary_table.csv
# =============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(Cairo)
})

set.seed(1997)

base_path <- normalizePath(".")   # Run from project root
setwd(base_path)

source("scripts/00_palette.R")

# =============================================================
# Global parameters
# =============================================================

SAMPLE_ORDER <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi", "mock_6wpi")
FDR_CUTOFF   <- 0.05
FC_CUTOFF    <- 0.25
TOP_N_LABEL  <- 15
MIN_PCT      <- 0.05
FC_THRESH    <- 0.1
SEED         <- 1997
MIN_CELLS    <- 10     # minimum cells per group to attempt FindMarkers

out_dir <- "outputs/banksy/immune_niche_volcano_by_celltype"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

SAMPLE_LABELS <- c(
  mock_6wpi = "Mock 6 wpi",
  LCMV_1wpi = "LCMV 1 wpi",
  LCMV_3wpi = "LCMV 3 wpi",
  LCMV_6wpi = "LCMV 6 wpi"
)

COND_PALETTE <- c(
  "Mock 6 wpi"  = "grey60",
  "LCMV 1 wpi"  = "#56B4E9",
  "LCMV 3 wpi"  = "#E69F00",
  "LCMV 6 wpi"  = "#D55E00"
)

direction_colors <- c(
  "up"   = "#B2182B",
  "down" = "#2166AC",
  "ns"   = "grey75"
)

CONTRASTS <- list(
  list(id1 = "LCMV_6wpi", id2 = "LCMV_1wpi")
)

# =============================================================
# Helpers
# =============================================================

make_slug <- function(x) {
  gsub("[^a-z0-9]+", "_", tolower(trimws(x)))
}

contrast_slug <- function(id1, id2) {
  paste0(make_slug(id1), "_vs_", make_slug(id2))
}

save_fig <- function(p, fname, width = 7, height = 5.5) {
  CairoPDF(file.path(out_dir, paste0(fname, ".pdf")), width = width, height = height)
  print(p)
  dev.off()
  CairoJPEG(file.path(out_dir, paste0(fname, ".jpg")),
            width = width * 150, height = height * 150, res = 150)
  print(p)
  dev.off()
  message("  Saved: ", fname, ".pdf / .jpg")
}

make_volcano <- function(mk, id1, id2, ct_label) {
  slug <- contrast_slug(id1, id2)
  lab1 <- SAMPLE_LABELS[id1]
  lab2 <- SAMPLE_LABELS[id2]

  mk$neg_log10_fdr <- -log10(pmax(mk$p_val_adj, .Machine$double.eps))
  x_cap <- max(3.5, quantile(abs(mk$avg_log2FC), 0.999, na.rm = TRUE) * 1.05)
  y_cap <- max(5,   quantile(mk$neg_log10_fdr,   0.999, na.rm = TRUE) * 1.05)
  mk$x_plot <- pmin(pmax(mk$avg_log2FC, -x_cap), x_cap)
  mk$y_plot <- pmin(mk$neg_log10_fdr, y_cap)

  sig_mk <- mk %>% filter(direction != "ns")
  if (nrow(sig_mk) > 0) {
    lab_genes <- sig_mk %>%
      arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
      slice_head(n = TOP_N_LABEL) %>%
      pull(gene)
  } else {
    lab_genes <- mk %>%
      arrange(p_val_adj) %>%
      slice_head(n = min(TOP_N_LABEL, 10)) %>%
      pull(gene)
  }
  mk$label <- ifelse(mk$gene %in% lab_genes, mk$gene, NA_character_)

  n_up   <- sum(mk$direction == "up",   na.rm = TRUE)
  n_down <- sum(mk$direction == "down", na.rm = TRUE)

  # Cell type colour from palette (fallback to grey)
  ct_colour <- if (ct_label %in% names(GLOBAL_PALETTE)) {
    GLOBAL_PALETTE[[ct_label]]
  } else {
    "grey40"
  }

  p <- ggplot(mk, aes(x = x_plot, y = y_plot, colour = direction)) +
    geom_point(size = 1.2, alpha = 0.7, stroke = 0) +
    geom_text_repel(
      data          = mk %>% filter(!is.na(label)),
      aes(label     = label),
      size          = 2.6,
      fontface      = "italic",
      max.overlaps  = 20,
      box.padding   = 0.3,
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
      labels = c(
        up   = paste0("Up in ", lab1, " (n=", n_up,   ")"),
        down = paste0("Down in ", lab1, " (n=", n_down, ")"),
        ns   = "Not significant"
      ),
      name = NULL
    ) +
    scale_x_continuous(limits = c(-x_cap, x_cap),
                       expand = expansion(mult = 0.02)) +
    scale_y_continuous(limits = c(0, y_cap),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(
      title    = paste0(ct_label, ": ", lab1, " vs ", lab2),
      subtitle = paste0(nrow(mk), " genes tested  |  FDR \u2264 ", FDR_CUTOFF,
                        "  |  |log2FC| > ", FC_CUTOFF),
      x        = paste0("log2FC (", lab1, " / ", lab2, ")"),
      y        = expression(-log[10](FDR))
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title       = element_text(face = "bold", size = 11, hjust = 0,
                                      colour = ct_colour),
      plot.subtitle    = element_text(size = 8, colour = "grey40", hjust = 0),
      legend.position  = "bottom",
      legend.text      = element_text(size = 8),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey93", linewidth = 0.3)
    )

  p
}

# =============================================================
# SECTION 1 — Load annotated immune object
# =============================================================

message("\n=== SECTION 1: Loading annotated immune object ===\n")

obj_file <- file.path("objects", "08_immune_annotated_lam02_res03.rds")
stopifnot("Object file not found" = file.exists(obj_file))

message("Loading: ", obj_file)
se <- readRDS(obj_file)
message("  ", ncol(se), " cells | class: ", class(se)[1])

# --- Confirm cell_type column ---
if (!"cell_type" %in% colnames(colData(se))) {
  stop("'cell_type' column not found in colData. Run script 17 first.")
}
if (!"sample" %in% colnames(colData(se))) {
  stop("'sample' column not found in colData.")
}

cell_types_all <- as.character(colData(se)$cell_type)
samples_all    <- as.character(colData(se)$sample)

cat("\nCells per cell type × condition:\n")
print(table(cell_types_all, samples_all))

# --- Convert to Seurat and normalize ---
assay_name <- if ("counts" %in% assayNames(se)) "counts" else assayNames(se)[1]
so <- suppressWarnings(as.Seurat(se, counts = assay_name, data = NULL))
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

# Add cell_type and sample to metadata (from original colData)
so$cell_type <- cell_types_all
so$sample    <- samples_all

cell_types <- sort(unique(cell_types_all))
message("\nCell types to process (", length(cell_types), "):")
for (ct in cell_types) message("  - ", ct)

# =============================================================
# SECTION 2 — DEG per cell type × contrast
# =============================================================

message("\n=== SECTION 2: DEG per cell type × contrast ===\n")

# Accumulate nDEG counts for Section 3 summary
ndeg_records <- list()

for (ct in cell_types) {

  ct_slug <- make_slug(ct)
  cells_ct <- which(so$cell_type == ct)
  so_ct    <- so[, cells_ct]

  message("\n--- Cell type: ", ct, " (", length(cells_ct), " cells) ---")

  # Set idents to sample for FindMarkers
  so_ct <- SetIdent(so_ct, value = "sample")

  for (ctr in CONTRASTS) {
    id1  <- ctr$id1
    id2  <- ctr$id2
    slug <- contrast_slug(id1, id2)

    n1 <- sum(so_ct$sample == id1, na.rm = TRUE)
    n2 <- sum(so_ct$sample == id2, na.rm = TRUE)

    message("  [", slug, "]  n=", n1, " vs n=", n2)

    # --- Skip if insufficient cells ---
    if (n1 < MIN_CELLS || n2 < MIN_CELLS) {
      message("    Skipping: insufficient cells (need >= ", MIN_CELLS,
              " per group, got ", n1, " vs ", n2, ")")
      ndeg_records <- c(ndeg_records, list(data.frame(
        cell_type = ct,
        contrast  = slug,
        n_up      = NA_integer_,
        n_down    = NA_integer_,
        n_total   = NA_integer_,
        skipped   = TRUE
      )))
      next
    }

    # --- FindMarkers ---
    mk <- tryCatch(
      FindMarkers(
        so_ct,
        ident.1         = id1,
        ident.2         = id2,
        only.pos        = FALSE,
        min.pct         = MIN_PCT,
        logfc.threshold = FC_THRESH,
        return.thresh   = 1,
        test.use        = "wilcox",
        verbose         = FALSE
      ),
      error = function(e) {
        message("    ERROR in FindMarkers: ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(mk) || nrow(mk) == 0) {
      message("    No results returned — skipping.")
      ndeg_records <- c(ndeg_records, list(data.frame(
        cell_type = ct,
        contrast  = slug,
        n_up      = 0L,
        n_down    = 0L,
        n_total   = 0L,
        skipped   = FALSE
      )))
      next
    }

    # --- Annotate results ---
    mk$gene      <- rownames(mk)
    mk$contrast  <- paste0(id1, "_vs_", id2)
    mk$cell_type <- ct
    mk$direction <- case_when(
      mk$p_val_adj <= FDR_CUTOFF & mk$avg_log2FC >  FC_CUTOFF ~ "up",
      mk$p_val_adj <= FDR_CUTOFF & mk$avg_log2FC < -FC_CUTOFF ~ "down",
      TRUE ~ "ns"
    )
    mk <- mk %>%
      arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
      select(gene, cell_type, contrast,
             avg_log2FC, pct.1, pct.2, p_val, p_val_adj, direction)

    n_up   <- sum(mk$direction == "up")
    n_down <- sum(mk$direction == "down")
    message("    ", nrow(mk), " genes tested | ", n_up, " up | ", n_down, " down",
            " (FDR<=", FDR_CUTOFF, ", |log2FC|>", FC_CUTOFF, ")")

    # --- Save CSV ---
    csv_name <- paste0("DEG_", ct_slug, "_", slug, ".csv")
    write.csv(mk, file.path(out_dir, csv_name), row.names = FALSE)
    message("    Saved: ", csv_name)

    # --- Volcano plot ---
    p <- make_volcano(mk, id1, id2, ct)
    fig_name <- paste0("fig_volcano_", ct_slug, "_", slug)
    save_fig(p, fig_name)

    # --- Accumulate nDEG for summary ---
    ndeg_records <- c(ndeg_records, list(data.frame(
      cell_type = ct,
      contrast  = slug,
      n_up      = n_up,
      n_down    = n_down,
      n_total   = n_up + n_down,
      skipped   = FALSE
    )))
  }
}

# =============================================================
# SECTION 3 — Summary figure: nDEG by cell type (6wpi vs 1wpi)
# =============================================================

message("\n=== SECTION 3: Summary nDEG figure ===\n")

if (length(ndeg_records) == 0) {
  message("No nDEG records to plot — skipping summary.")
} else {

  summary_df <- bind_rows(ndeg_records) %>%
    mutate(
      n_total_plot = ifelse(skipped, 0L, n_total),
      n_up_plot = ifelse(skipped | is.na(n_up), 0L, n_up),
      n_down_plot = ifelse(skipped | is.na(n_down), 0L, n_down),
      contrast_label = factor("6wpi vs 1wpi", levels = "6wpi vs 1wpi")
    )

  # Order cell types by total nDEG for 6wpi vs 1wpi
  ct_order <- summary_df %>%
    group_by(cell_type) %>%
    summarise(total = sum(n_total_plot, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total)) %>%
    pull(cell_type)

  summary_df <- summary_df %>%
    mutate(cell_type = factor(cell_type, levels = rev(ct_order)))

  # Save summary table
  write.csv(summary_df,
            file.path(out_dir, "ndeg_summary_table.csv"),
            row.names = FALSE)
  message("  Saved: ndeg_summary_table.csv")

  # --- Barplot ---
  summary_long <- summary_df %>%
    select(cell_type, n_up_plot, n_down_plot, n_total_plot) %>%
    pivot_longer(
      cols = c(n_up_plot, n_down_plot),
      names_to = "direction",
      values_to = "n_deg"
    ) %>%
    mutate(
      direction = factor(
        direction,
        levels = c("n_down_plot", "n_up_plot"),
        labels = c("Down in 6wpi", "Up in 6wpi")
      )
    )

  p_summary <- ggplot(
    summary_long,
    aes(x = cell_type, y = n_deg, fill = direction)
  ) +
    geom_col(width = 0.7, colour = "grey30", linewidth = 0.25) +
    geom_text(
      data = summary_df,
      aes(x = cell_type,
          y = n_total_plot,
          label = ifelse(
            n_total_plot > 0,
            paste0(n_up_plot, " up / ", n_down_plot, " down"),
            ""
          )),
      hjust    = -0.15,
      size     = 2.8,
      colour   = "grey20",
      inherit.aes = FALSE
    ) +
    coord_flip() +
    scale_fill_manual(
      values = c("Down in 6wpi" = "#2166AC", "Up in 6wpi" = "#B2182B"),
      name = NULL
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
      title    = "Number of significant DEGs per cell type",
      subtitle = paste0("FDR \u2264 ", FDR_CUTOFF, "  |  |log2FC| > ", FC_CUTOFF,
                        "  |  Contrast: LCMV 6wpi vs LCMV 1wpi"),
      x        = NULL,
      y        = "N significant DEGs"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title       = element_text(face = "bold", size = 11, hjust = 0),
      plot.subtitle    = element_text(size = 8, colour = "grey40", hjust = 0),
      legend.position  = "right",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
      axis.text.y      = element_text(size = 9)
    )

  n_ct_shown <- length(unique(summary_df$cell_type[!summary_df$skipped]))
  fig_height <- max(4.5, 0.45 * n_ct_shown + 2)

  save_fig(p_summary, "fig_ndeg_summary_by_celltype",
           width = 8, height = fig_height)

  # --- Print summary to console ---
  cat("\n=== nDEG Summary ===\n")
  print(as.data.frame(
    summary_df %>%
      select(cell_type, contrast_label, n_up, n_down, n_total_plot, skipped) %>%
      arrange(cell_type, contrast_label)
  ))
}

message("\n=== Script 34 complete ===")
message("Outputs in: ", out_dir)
