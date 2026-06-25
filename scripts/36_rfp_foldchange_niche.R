#!/usr/bin/env Rscript
# =============================================================
# Script: 36_rfp_foldchange_niche.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-06
#
# Goal:
#   Test whether RFP is enriched in the immune niche and IFN-responsive
#   cells using an intra-sample fold change (cancels slide batch effects).
#
#   Starting point: objects/03_rfp_analysis.rds
#   (built by 34_rebuild_protein_metadata.R — protein columns correctly joined
#    via EntityID + CellType already transferred from SPE; no back-fill needed)
#
#   Niche groups : NICHE_GROUPS (default: "Immune (Acod1)", "IFN responsive (Ifit1)")
#   "Hors niche"  : all other annotated cells (Non annote excluded everywhere)
#
#   Three RFP measures:
#     raw      : Anti.RFP_raw
#     highpass : Anti.RFP_high_pass
#     vol_norm : Anti.RFP_raw / volume  (guards against cell-size confounding)
#
# ─────────────────────────────────────────────────────────────
# Analysis 1 — Fold change of intensity (intra-sample)
#   log2FC_median = log2(median_in / median_out)
#   log2FC_mean   = log2(mean_in   / mean_out)
#   Columns: sample, group, measure,
#            median_in, median_out, mean_in, mean_out,
#            log2FC_median, log2FC_mean, n_in, n_out
#
# Analysis 2 — Fold enrichment of positivity (self-referenced threshold)
#   threshold per sample = q99(RFP in ALL annotated cells of THAT sample)
#   fold_enrich = pct_pos_in / pct_pos_out
#   Columns: sample, group, measure,
#            threshold_q99, pct_pos_in, pct_pos_out, fold_enrich, n_in, n_out
#
# ─────────────────────────────────────────────────────────────
# Outputs: outputs/rfp_foldchange_niche/
#   rfp_fc_intensity.csv          — Analysis 1 results
#   rfp_fc_positivity.csv         — Analysis 2 results
#   fig_fc_intensity.pdf/jpg      — barplot log2FC intensity
#   fig_fc_positivity.pdf/jpg     — barplot fold enrichment positivity
# =============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(readr)
})

setwd(normalizePath("."))   # Run from project root

# =============================================================
# 1. Parameters
# =============================================================

# Built by 34_rebuild_protein_metadata.R: correct protein columns (EntityID
# join) + CellType from SPE. Use this, NOT 03_all_clustered.rds which still
# carries the corrupted (row-order) protein assignment.
OBJ_SEURAT   <- file.path("objects", "03_rfp_analysis.rds")
OUT_DIR      <- file.path("outputs", "rfp_foldchange_niche")

SAMPLE_ORDER <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi", "mock_6wpi")

# Cell types defining the niche of interest
NICHE_GROUPS <- c("Immune (Acod1)", "IFN responsive (Ifit1)")

# Label used to mark unannotated cells (excluded from all analyses)
UNANNOTATED  <- "Non annote"

# Self-referenced positivity threshold (per sample)
THRESHOLD_Q  <- 0.99

# Minimum cells per group×sample to compute statistics (warn if below)
MIN_CELLS    <- 5

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE,
                                      showWarnings = FALSE)

# =============================================================
# 2. Helper: save ggplot as PDF + JPG
# =============================================================

save_plot <- function(p, stem, width = 12, height = 8) {
  ggsave(paste0(stem, ".pdf"), p, width = width, height = height,
         device = "pdf")
  ggsave(paste0(stem, ".jpg"), p, width = width, height = height,
         dpi = 150, device = "jpeg")
  message("  Saved: ", basename(stem), ".pdf / .jpg")
}

# =============================================================
# 3. Load Seurat object
# =============================================================

message("\n=== Loading Seurat object ===")
message("  Path: ", OBJ_SEURAT)
obj <- readRDS(OBJ_SEURAT)
message("  Cells  : ", ncol(obj))
message("  Samples: ", paste(sort(unique(obj$sample)), collapse = ", "))

# Sanity check: protein columns and CellType must already be present
required_cols <- c("Anti.RFP_raw", "Anti.RFP_high_pass", "volume", "CellType")
missing_cols  <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing column(s) in object: ", paste(missing_cols, collapse = ", "),
       "\n  -> Run scripts/34_rebuild_protein_metadata.R first.")
}
message("  Required columns present: ",
        paste(required_cols, collapse = ", "))
message("  CellType table:")
print(sort(table(obj$CellType), decreasing = TRUE))

# =============================================================
# 4. Verify niche group labels
# =============================================================

message("\n=== Verifying niche group labels ===")
all_celltypes <- unique(obj$CellType)
missing_groups <- setdiff(NICHE_GROUPS, all_celltypes)
if (length(missing_groups) > 0) {
  message("  WARNING: The following NICHE_GROUPS were not found in CellType:")
  for (g in missing_groups) message("    - '", g, "'")
  message("  Available CellType values:")
  for (ct in sort(all_celltypes)) message("    '", ct, "'")
  stop("Fix NICHE_GROUPS to match exact CellType labels above.")
}
message("  All niche groups found: ", paste(NICHE_GROUPS, collapse = " | "))

# =============================================================
# 5. Build vol_norm column (Anti.RFP_raw / volume)
# =============================================================

message("\n=== Building volume-normalised RFP column ===")
raw_vals  <- as.numeric(obj@meta.data$Anti.RFP_raw)
vol_vals  <- as.numeric(obj@meta.data$volume)

n_zero_vol <- sum(!is.na(vol_vals) & vol_vals <= 0)
if (n_zero_vol > 0) {
  message("  WARNING: ", n_zero_vol, " cells with volume <= 0 -> vol_norm = NA")
}
obj$RFP_raw_over_volume <- ifelse(!is.na(vol_vals) & vol_vals > 0,
                                   raw_vals / vol_vals,
                                   NA_real_)
message("  RFP_raw_over_volume: NAs = ",
        sum(is.na(obj$RFP_raw_over_volume)),
        " / ", ncol(obj))

# =============================================================
# 6. Define RFP measures
# =============================================================

MEASURES <- list(
  raw      = "Anti.RFP_raw",
  highpass = "Anti.RFP_high_pass",
  vol_norm = "RFP_raw_over_volume"
)

# =============================================================
# 7. Set up working data frame (annotated cells only)
# =============================================================

message("\n=== Building analysis data frame ===")

obj$sample <- factor(obj$sample,
                     levels = intersect(SAMPLE_ORDER, unique(as.character(obj$sample))))

ana_df <- as.data.frame(obj@meta.data) %>%
  filter(CellType != UNANNOTATED) %>%
  mutate(sample = factor(sample, levels = levels(obj$sample)))

message("  Annotated cells: ", nrow(ana_df),
        " / ", ncol(obj), " total")
message("  Samples in data: ",
        paste(levels(droplevels(ana_df$sample)), collapse = ", "))

# Quick n per group × sample
message("\n  Cells per CellType × sample:")
ct_tbl <- ana_df %>%
  count(sample, CellType) %>%
  tidyr::pivot_wider(names_from = sample, values_from = n, values_fill = 0)
print(as.data.frame(ct_tbl))

# =============================================================
# 8. Analysis 1 — Fold change of intensity (intra-sample)
# =============================================================

message("\n=== Analysis 1: Fold change of intensity ===")

fc_rows <- list()

for (samp in levels(ana_df$sample)) {
  df_s <- ana_df %>% filter(sample == samp)

  for (grp in NICHE_GROUPS) {
    cells_in  <- df_s %>% filter(CellType == grp)
    cells_out <- df_s %>% filter(CellType != grp)

    n_in  <- nrow(cells_in)
    n_out <- nrow(cells_out)

    if (n_in < MIN_CELLS) {
      message("  WARNING: n_in = ", n_in, " (< ", MIN_CELLS, ") for ",
              grp, " in ", samp, " -- skipping intensity FC")
      next
    }
    if (n_out < MIN_CELLS) {
      message("  WARNING: n_out = ", n_out, " (< ", MIN_CELLS, ") for ",
              "hors niche in ", samp, " -- skipping intensity FC")
      next
    }

    for (mname in names(MEASURES)) {
      col <- MEASURES[[mname]]
      v_in  <- cells_in[[col]][!is.na(cells_in[[col]])]
      v_out <- cells_out[[col]][!is.na(cells_out[[col]])]

      med_in  <- median(v_in)
      med_out <- median(v_out)
      mea_in  <- mean(v_in)
      mea_out <- mean(v_out)

      log2fc_med <- if (med_out  > 0) log2(med_in  / med_out)  else NA_real_
      log2fc_mea <- if (mea_out  > 0) log2(mea_in  / mea_out)  else NA_real_

      fc_rows[[length(fc_rows) + 1]] <- data.frame(
        sample        = samp,
        group         = grp,
        measure       = mname,
        median_in     = med_in,
        median_out    = med_out,
        mean_in       = mea_in,
        mean_out      = mea_out,
        log2FC_median = log2fc_med,
        log2FC_mean   = log2fc_mea,
        n_in          = n_in,
        n_out         = n_out,
        stringsAsFactors = FALSE
      )
    }
  }
}

fc_df <- bind_rows(fc_rows)
fc_df$sample  <- factor(fc_df$sample,  levels = levels(ana_df$sample))
fc_df$group   <- factor(fc_df$group,   levels = NICHE_GROUPS)
fc_df$measure <- factor(fc_df$measure, levels = names(MEASURES))

message("  Rows: ", nrow(fc_df))
print(fc_df %>% select(sample, group, measure, log2FC_median, log2FC_mean, n_in))

csv_fc <- file.path(OUT_DIR, "rfp_fc_intensity.csv")
write_csv(fc_df, csv_fc)
message("  Saved: rfp_fc_intensity.csv")

# =============================================================
# 9. Analysis 2 — Fold enrichment of positivity (self-referenced q99)
# =============================================================

message("\n=== Analysis 2: Fold enrichment of positivity ===")

# Per-sample, per-measure q99 threshold (computed on ALL annotated cells)
thresholds <- list()
for (samp in levels(ana_df$sample)) {
  df_s <- ana_df %>% filter(sample == samp)
  for (mname in c("raw", "highpass")) {
    col <- MEASURES[[mname]]
    vals <- df_s[[col]][!is.na(df_s[[col]])]
    thresholds[[paste(samp, mname, sep = "__")]] <- quantile(vals, probs = THRESHOLD_Q)
  }
}

enrich_rows <- list()

for (samp in levels(ana_df$sample)) {
  df_s <- ana_df %>% filter(sample == samp)

  for (grp in NICHE_GROUPS) {
    cells_in  <- df_s %>% filter(CellType == grp)
    cells_out <- df_s %>% filter(CellType != grp)

    n_in  <- nrow(cells_in)
    n_out <- nrow(cells_out)

    if (n_in < MIN_CELLS) {
      message("  WARNING: n_in = ", n_in, " (< ", MIN_CELLS, ") for ",
              grp, " in ", samp, " -- skipping positivity enrichment")
      next
    }
    if (n_out < MIN_CELLS) {
      message("  WARNING: n_out = ", n_out, " (< ", MIN_CELLS, ") for ",
              "hors niche in ", samp, " -- skipping positivity enrichment")
      next
    }

    for (mname in c("raw", "highpass")) {
      col  <- MEASURES[[mname]]
      thr  <- thresholds[[paste(samp, mname, sep = "__")]]

      v_in  <- cells_in[[col]]
      v_out <- cells_out[[col]]

      pct_in  <- 100 * sum(!is.na(v_in)  & v_in  > thr) / n_in
      pct_out <- 100 * sum(!is.na(v_out) & v_out > thr) / n_out

      fe <- if (pct_out > 0) pct_in / pct_out else NA_real_

      enrich_rows[[length(enrich_rows) + 1]] <- data.frame(
        sample        = samp,
        group         = grp,
        measure       = mname,
        threshold_q99 = as.numeric(thr),
        pct_pos_in    = round(pct_in,  3),
        pct_pos_out   = round(pct_out, 3),
        fold_enrich   = round(fe, 4),
        n_in          = n_in,
        n_out         = n_out,
        stringsAsFactors = FALSE
      )
    }
  }
}

enrich_df <- bind_rows(enrich_rows)
enrich_df$sample  <- factor(enrich_df$sample,  levels = levels(ana_df$sample))
enrich_df$group   <- factor(enrich_df$group,   levels = NICHE_GROUPS)
enrich_df$measure <- factor(enrich_df$measure, levels = c("raw", "highpass"))

message("  Rows: ", nrow(enrich_df))
print(enrich_df %>% select(sample, group, measure, pct_pos_in, pct_pos_out,
                            fold_enrich, n_in))

csv_enrich <- file.path(OUT_DIR, "rfp_fc_positivity.csv")
write_csv(enrich_df, csv_enrich)
message("  Saved: rfp_fc_positivity.csv")

# =============================================================
# 10. Figure 1 — Barplot log2FC intensity
#     x = sample, facet_grid(group ~ .), colour = measure
#     Horizontal line at y = 0. Two sub-panels: median / mean.
# =============================================================

message("\n=== Figure 1: log2FC intensity barplot ===")

MEASURE_COLORS <- c(
  raw      = "#4393C3",   # blue
  highpass = "#D6604D",   # red-orange
  vol_norm = "#74C476"    # green
)

dodge <- position_dodge(width = 0.7)

make_fc_panel <- function(df, ycol, ylab) {
  ggplot(df, aes(x = sample, y = .data[[ycol]],
                 fill = measure, group = measure)) +
    geom_col(position = dodge, width = 0.65, color = "white", linewidth = 0.2) +
    geom_hline(yintercept = 0, linewidth = 0.5, color = "grey30") +
    facet_grid(group ~ ., labeller = label_wrap_gen(width = 22)) +
    scale_fill_manual(values = MEASURE_COLORS,
                      labels = c(raw = "Anti.RFP_raw",
                                 highpass = "Anti.RFP_high_pass",
                                 vol_norm = "raw / volume")) +
    labs(x = NULL, y = ylab, fill = "Mesure RFP") +
    theme_classic(base_size = 10) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      axis.text.x      = element_text(angle = 35, hjust = 1, size = 9),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(face = "bold", size = 9),
      legend.position  = "bottom"
    )
}

p_fc_med <- make_fc_panel(fc_df, "log2FC_median",
                           ylab = "log2FC intensité (médiane in / out)")
p_fc_mea <- make_fc_panel(fc_df, "log2FC_mean",
                           ylab = "log2FC intensité (moyenne in / out)")

p_fc <- (p_fc_med + ggtitle("Médiane")) +
        (p_fc_mea + ggtitle("Moyenne")) +
  plot_annotation(
    title    = "Fold change d'intensité RFP — niche vs. hors niche (intra-échantillon)",
    subtitle = paste0("log2(médiane_in / médiane_out) et log2(moyenne_in / moyenne_out)  |  ",
                      "hors niche = toutes cellules annotées hors groupe"),
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.title      = element_text(face = "bold", size = 11),
      plot.subtitle   = element_text(size = 9, color = "grey40")
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

save_plot(p_fc,
          file.path(OUT_DIR, "fig_fc_intensity"),
          width = 14, height = 5 * length(NICHE_GROUPS))

# =============================================================
# 11. Console summary
# =============================================================

message("\n", strrep("=", 60))
message("SUMMARY — 36_rfp_foldchange_niche")
message(strrep("=", 60))
message("  Annotated cells (excl. '", UNANNOTATED, "'): ", nrow(ana_df))
message("  Samples: ", paste(levels(ana_df$sample), collapse = ", "))
message("  Niche groups: ", paste(NICHE_GROUPS, collapse = " | "))
message("")
message("  Analysis 1 (FC intensity) rows: ", nrow(fc_df))
message("  Analysis 2 (FC positivity) rows: ", nrow(enrich_df))
message("")
message("  CSVs:")
message("    ", normalizePath(csv_fc,     mustWork = FALSE))
message("    ", normalizePath(csv_enrich, mustWork = FALSE))
message("")
message("  Figures:")
message("    fig_fc_intensity.pdf/jpg")
message("")
message("  Output folder: ", normalizePath(OUT_DIR))
message(strrep("=", 60))
