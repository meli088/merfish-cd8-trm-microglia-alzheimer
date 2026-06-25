#!/usr/bin/env Rscript
# =============================================================
# Script: 35_rfp_by_celltype.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-06
#
# Goal:
#   RFP signal analysis by cell type (Zoé-style):
#     - Threshold  : q99 of mock (raw) / q99 of mock (highpass)
#     - RFP_positive     : Anti.RFP_raw        > p99_mock_raw  -> "RFP_pos"/"RFP_neg"
#     - RFP_positive_high: Anti.RFP_high_pass  > p99_mock_high -> "RFP_pos"/"RFP_neg"
#
#   Figures produced:
#     fig_vln_rfp_raw_by_celltype.pdf/jpg   - VlnPlot raw, group by CellType
#     fig_vln_rfp_raw_by_sample.pdf/jpg     - VlnPlot raw, group by sample
#     fig_vln_rfp_high_by_celltype.pdf/jpg  - VlnPlot highpass, group by CellType
#     fig_vln_rfp_high_by_sample.pdf/jpg    - VlnPlot highpass, group by sample
#     fig_dim_rfp_positive.pdf/jpg          - DimPlot RFP_positive (powderblue/tomato)
#     fig_dim_rfp_positive_high.pdf/jpg     - DimPlot RFP_positive_high (powderblue/darkorange)
#     fig_feature_rfp_raw.pdf/jpg           - FeaturePlot Anti.RFP_raw
#     fig_feature_rfp_high.pdf/jpg          - FeaturePlot Anti.RFP_high_pass
#     rfp_positive_by_celltype_sample.csv   - RFP+ % by CellType x sample
#
# NOTE - Anti.RFP_high_pass:
#   Exists in source cell_metadata.csv files (100% non-NA across all 4 samples)
#   but was not transferred in 00_load_data.R. This script back-fills it
#   in-memory from the source CSVs (section 5b). Join key: obj$cell = EntityID.
#
# Outputs: outputs/rfp_by_celltype/
# =============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(readr)
})

setwd(normalizePath("."))   # Run from project root

# =============================================================
# 1. Parameters
# =============================================================

# Built by 34_rebuild_protein_metadata.R from 02_all_normalized.rds:
# correct EntityID join for protein channels + CellType from SPE
OBJ_SEURAT   <- file.path("objects", "03_rfp_analysis.rds")
OUT_DIR      <- file.path("outputs", "rfp_by_celltype")

SAMPLE_ORDER <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi", "mock_6wpi")
MOCK_PATTERN <- "mock"   # case-insensitive grep on $sample
THRESHOLD_Q  <- 0.99     # quantile applied directly to raw mock values

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
message("  Cells: ", ncol(obj))
message("  Samples: ", paste(sort(unique(obj$sample)), collapse = ", "))

# =============================================================
# 4. Diagnostic: metadata columns + sample labels
# =============================================================

message("\n=== Metadata columns ===")
cat(paste(colnames(obj@meta.data), collapse = "\n"), "\n\n")

message("=== Sample column values ===")
if ("sample" %in% colnames(obj@meta.data)) {
  print(sort(table(obj$sample), decreasing = TRUE))
} else {
  stop("'sample' column not found in metadata.")
}

message("=== DAPI/PolyT columns in metadata ===")
print(grep("DAPI|PolyT", colnames(obj@meta.data), value = TRUE))

# =============================================================
# 5. Auto-detect RFP columns
# =============================================================

message("\n=== Detecting RFP columns ===")
meta_cols <- colnames(obj@meta.data)

rfp_raw_candidates  <- meta_cols[grepl("rfp", meta_cols, ignore.case = TRUE) &
                                   !grepl("high|pass", meta_cols, ignore.case = TRUE)]
rfp_high_candidates <- meta_cols[grepl("rfp", meta_cols, ignore.case = TRUE) &
                                    grepl("high|pass", meta_cols, ignore.case = TRUE)]

if (length(rfp_raw_candidates) == 0) {
  stop("No RFP raw column found in metadata.")
}
rfp_raw_col  <- rfp_raw_candidates[1]
rfp_high_col <- if (length(rfp_high_candidates) > 0) rfp_high_candidates[1] else NULL

message("  RFP raw  : ", rfp_raw_col)
message("  RFP high : ", if (is.null(rfp_high_col)) "NOT FOUND" else rfp_high_col)

if (is.null(rfp_high_col)) {
  stop("Anti.RFP_high_pass not found in object. ",
       "Run scripts/34_rebuild_protein_metadata.R first.")
}

if (!"CellType" %in% colnames(obj@meta.data)) {
  stop("'CellType' not found in object. ",
       "Run scripts/34_rebuild_protein_metadata.R first.")
}
message("  CellType table:")
print(sort(table(obj$CellType), decreasing = TRUE))

# =============================================================
# 6. Guardrails: PolyT/DAPI checks + safe divisions (avoid Inf/NaN)
# =============================================================

message("\n=== PolyT/DAPI guardrails before division ===")

required_norm_cols <- c("PolyT_high_pass", "DAPI_high_pass")
missing_norm_cols <- setdiff(required_norm_cols, colnames(obj@meta.data))
if (length(missing_norm_cols) > 0) {
  stop("Missing required metadata column(s): ",
       paste(missing_norm_cols, collapse = ", "))
}

polyT <- as.numeric(obj@meta.data$PolyT_high_pass)
dapi  <- as.numeric(obj@meta.data$DAPI_high_pass)
rfp_raw_vals  <- as.numeric(obj@meta.data[[rfp_raw_col]])
rfp_high_vals <- as.numeric(obj@meta.data[[rfp_high_col]])

message("  PolyT_high_pass: NAs = ", sum(is.na(polyT)),
  " | <=0 = ", sum(!is.na(polyT) & polyT <= 0))
message("  DAPI_high_pass : NAs = ", sum(is.na(dapi)),
  " | <=0 = ", sum(!is.na(dapi) & dapi <= 0))

# Safe division prevents Inf/NaN when denominator is missing or non-positive.
obj$RFP_raw_over_PolyT <- ifelse(!is.na(polyT) & polyT > 0,
         rfp_raw_vals / polyT,
         NA_real_)
obj$RFP_high_over_PolyT <- ifelse(!is.na(polyT) & polyT > 0,
          rfp_high_vals / polyT,
          NA_real_)
obj$RFP_raw_over_DAPI <- ifelse(!is.na(dapi) & dapi > 0,
        rfp_raw_vals / dapi,
        NA_real_)
obj$RFP_high_over_DAPI <- ifelse(!is.na(dapi) & dapi > 0,
         rfp_high_vals / dapi,
         NA_real_)

message("  Safe ratio columns added: ",
  paste(c("RFP_raw_over_PolyT", "RFP_high_over_PolyT",
    "RFP_raw_over_DAPI", "RFP_high_over_DAPI"), collapse = ", "))

# =============================================================
# 7. Compute q99 thresholds in mock (Zoe approach -- no trim)
# =============================================================

message("\n=== Computing q99 thresholds in mock ===")

all_samples  <- unique(obj$sample)
mock_samples <- all_samples[grepl(MOCK_PATTERN, all_samples, ignore.case = TRUE)]
if (length(mock_samples) == 0) {
  stop("No mock sample found with pattern '", MOCK_PATTERN, "'.")
}
message("  Mock sample(s): ", paste(mock_samples, collapse = ", "))

mock_cells    <- which(obj$sample %in% mock_samples)

rfp_raw_mock  <- obj@meta.data[mock_cells, rfp_raw_col]
p99_mock_raw  <- quantile(rfp_raw_mock[!is.na(rfp_raw_mock)], probs = THRESHOLD_Q)
message("  p99 mock (raw)  : ",
        format(round(p99_mock_raw, 0), big.mark = ",", scientific = FALSE))

rfp_high_mock <- obj@meta.data[mock_cells, rfp_high_col]
p99_mock_high <- quantile(rfp_high_mock[!is.na(rfp_high_mock)], probs = THRESHOLD_Q)
message("  p99 mock (high) : ",
        format(round(p99_mock_high, 0), big.mark = ",", scientific = FALSE))

# =============================================================
# 7b. Diagnostic: RFP value distributions (all samples, mock highlighted)
#
#   Purpose: verify that the object contains meaningful RFP signal and
#   that the q99 mock threshold falls in a sensible place.
#   Two panels: Anti.RFP_raw | Anti.RFP_high_pass
#   Each panel: density curves per sample (mock = solid/red, LCMV = dashed)
#               + vertical red dashed line at q99 mock threshold
#               + x-axis clipped at q99.9 to avoid extreme outlier stretch
#   Output: fig_diag_rfp_distributions.pdf/jpg
# =============================================================

message("\n=== Diagnostic: RFP distributions ===")

diag_df <- as.data.frame(obj@meta.data[, c("sample", rfp_raw_col, rfp_high_col),
                                         drop = FALSE])
diag_df$sample <- factor(diag_df$sample,
                          levels = intersect(SAMPLE_ORDER, unique(diag_df$sample)))

SAMPLE_COLORS <- c(mock_6wpi  = "#E41A1C",   # red   — control
                   LCMV_1wpi  = "#56B4E9",   # blue
                   LCMV_3wpi  = "#E69F00",   # orange
                   LCMV_6wpi  = "#D55E00")   # vermillion
SAMPLE_LT     <- c(mock_6wpi  = "solid",
                   LCMV_1wpi  = "dashed",
                   LCMV_3wpi  = "dashed",
                   LCMV_6wpi  = "dashed")

make_diag_density <- function(df, col, threshold, col_label) {
  vals  <- df[[col]][!is.na(df[[col]])]
  x_max <- quantile(vals, probs = 0.999)   # clip extreme right tail
  x_min <- 0

  ggplot(df[!is.na(df[[col]]), ],
         aes(x = .data[[col]], color = sample, linetype = sample)) +
    geom_density(linewidth = 0.8, adjust = 1.5) +
    geom_vline(xintercept = threshold, color = "red",
               linetype = "dashed", linewidth = 0.9) +
    annotate("text",
             x = threshold, y = Inf,
             label = paste0("q99 mock\n= ",
                            format(round(threshold, 0), big.mark = ",")),
             hjust = -0.05, vjust = 1.3, color = "red", size = 3) +
    coord_cartesian(xlim = c(x_min, x_max)) +
    scale_x_continuous(
      labels = scales::label_number(scale = 1e-6, suffix = "M")
    ) +
    scale_color_manual(values    = SAMPLE_COLORS,
                       breaks    = names(SAMPLE_COLORS),
                       na.value  = "grey60") +
    scale_linetype_manual(values   = SAMPLE_LT,
                          breaks   = names(SAMPLE_LT)) +
    labs(x = paste0(col_label, " intensity"),
         y = "Density",
         color    = NULL,
         linetype = NULL,
         title    = col_label) +
    theme_classic(base_size = 10) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      legend.position  = "bottom",
      legend.key.width = unit(1.5, "cm"),
      plot.title       = element_text(face = "bold")
    )
}

p_diag_raw  <- make_diag_density(diag_df, rfp_raw_col,
                                  p99_mock_raw,  rfp_raw_col)
p_diag_high <- make_diag_density(diag_df, rfp_high_col,
                                  p99_mock_high, rfp_high_col)

p_diag <- p_diag_raw + p_diag_high +
  plot_annotation(
    title    = "RFP signal distributions — all samples",
    subtitle = paste0(
      "Density clipped at q99.9  |  Red dashed = q",
      THRESHOLD_Q * 100, " mock threshold  |  ",
      "n_cells = ", nrow(diag_df)
    ),
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.title      = element_text(face = "bold", size = 12),
      plot.subtitle   = element_text(size = 9, color = "grey40")
    )
  )

save_plot(p_diag,
          file.path(OUT_DIR, "fig_diag_rfp_distributions"),
          width = 14, height = 6)

# =============================================================
# 8. Classify RFP_positive + RFP_positive_high
# =============================================================

message("\n=== Classifying RFP+ cells ===")

obj$RFP_positive <- factor(
  ifelse(obj@meta.data[[rfp_raw_col]] > p99_mock_raw, "RFP_pos", "RFP_neg"),
  levels = c("RFP_neg", "RFP_pos")
)

obj$RFP_positive_high <- factor(
  ifelse(obj@meta.data[[rfp_high_col]] > p99_mock_high, "RFP_pos", "RFP_neg"),
  levels = c("RFP_neg", "RFP_pos")
)

n_total    <- ncol(obj)
n_pos      <- sum(obj$RFP_positive      == "RFP_pos", na.rm = TRUE)
n_pos_high <- sum(obj$RFP_positive_high == "RFP_pos", na.rm = TRUE)

message("  RFP_pos (raw)  : ", n_pos,      " / ", n_total,
        " (", round(100 * n_pos      / n_total, 2), "%)")
message("  RFP_pos (high) : ", n_pos_high, " / ", n_total,
        " (", round(100 * n_pos_high / n_total, 2), "%)")

# =============================================================
# 9. Violin plots (ggplot, CellType on x, faceted by sample)
# =============================================================

message("\n=== Violin plots ===")

sample_levels <- intersect(SAMPLE_ORDER, unique(as.character(obj$sample)))
obj$sample    <- factor(obj$sample, levels = sample_levels)

# Exclude "Non annote" for cleaner violins
vln_df <- as.data.frame(obj@meta.data) %>%
  filter(CellType != "Non annote") %>%
  mutate(sample = factor(sample, levels = sample_levels))

make_vln <- function(df, ycol, ylab, threshold, thr_label) {
  ggplot(df, aes(x = CellType, y = .data[[ycol]], fill = CellType)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.3) +
    geom_hline(yintercept = threshold, linetype = "dashed",
               color = "red", linewidth = 0.7) +
    annotate("text", x = Inf, y = threshold, label = thr_label,
             hjust = 1.05, vjust = -0.4, color = "red", size = 3) +
    facet_wrap(~ sample, ncol = 2) +
    scale_y_continuous(
      labels = scales::label_number(scale = 1e-6, suffix = "M")
    ) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 10) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(face = "bold"),
      legend.position  = "none"
    )
}

p_vln_raw <- make_vln(vln_df, rfp_raw_col,
                       ylab      = paste0(rfp_raw_col, " intensity"),
                       threshold = p99_mock_raw,
                       thr_label = paste0("q99 mock = ",
                                          format(round(p99_mock_raw, 0),
                                                 big.mark = ",")))
save_plot(p_vln_raw,
          file.path(OUT_DIR, "fig_vln_rfp_raw_by_celltype"),
          width = 14, height = 10)

p_vln_high <- make_vln(vln_df, rfp_high_col,
                        ylab      = paste0(rfp_high_col, " intensity"),
                        threshold = p99_mock_high,
                        thr_label = paste0("q99 mock = ",
                                           format(round(p99_mock_high, 0),
                                                  big.mark = ",")))
save_plot(p_vln_high,
          file.path(OUT_DIR, "fig_vln_rfp_high_by_celltype"),
          width = 14, height = 10)

# Quick overview by sample (Seurat VlnPlot)
p_vln_raw_sample <- VlnPlot(obj, features = rfp_raw_col,
                              group.by = "sample", pt.size = 0) +
  geom_hline(yintercept = p99_mock_raw, linetype = "dashed",
             color = "red", linewidth = 0.7) +
  labs(title = paste0(rfp_raw_col, "  |  red dashed = q99 mock"), x = NULL) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        legend.position = "none")
save_plot(p_vln_raw_sample,
          file.path(OUT_DIR, "fig_vln_rfp_raw_by_sample"),
          width = 10, height = 6)

p_vln_high_sample <- VlnPlot(obj, features = rfp_high_col,
                               group.by = "sample", pt.size = 0) +
  geom_hline(yintercept = p99_mock_high, linetype = "dashed",
             color = "red", linewidth = 0.7) +
  labs(title = paste0(rfp_high_col, "  |  red dashed = q99 mock"), x = NULL) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        legend.position = "none")
save_plot(p_vln_high_sample,
          file.path(OUT_DIR, "fig_vln_rfp_high_by_sample"),
          width = 10, height = 6)

# =============================================================
# 10. Summary CSV: RFP+ count and % by CellType x sample
# =============================================================

message("\n=== Summary CSV ===")

summary_df <- obj@meta.data %>%
  as.data.frame() %>%
  filter(CellType != "Non annote") %>%
  group_by(CellType, sample) %>%
  summarise(
    n_cells          = n(),
    n_rfp_pos_raw    = sum(RFP_positive      == "RFP_pos", na.rm = TRUE),
    pct_rfp_pos_raw  = round(100 * n_rfp_pos_raw  / n_cells, 2),
    n_rfp_pos_high   = sum(RFP_positive_high == "RFP_pos", na.rm = TRUE),
    pct_rfp_pos_high = round(100 * n_rfp_pos_high / n_cells, 2),
    .groups = "drop"
  ) %>%
  arrange(CellType, sample)

csv_path <- file.path(OUT_DIR, "rfp_positive_by_celltype_sample.csv")
write_csv(summary_df, csv_path)
message("  Saved: rfp_positive_by_celltype_sample.csv (", nrow(summary_df), " rows)")

# =============================================================
# 13. Console summary
# =============================================================

message("\n", strrep("=", 60))
message("SUMMARY")
message(strrep("=", 60))
message("  Cells total        : ", n_total)
message("  Mock sample(s)     : ", paste(mock_samples, collapse = ", "))
message("  q99 mock (raw)     : ",
        format(round(p99_mock_raw, 0), big.mark = ",", scientific = FALSE))
message("  q99 mock (high)    : ",
        format(round(p99_mock_high, 0), big.mark = ",", scientific = FALSE))
message("  RFP_pos (raw)      : ", n_pos,
        " (", round(100 * n_pos / n_total, 2), "% of all cells)")
message("  RFP_pos (highpass) : ", n_pos_high,
        " (", round(100 * n_pos_high / n_total, 2), "% of all cells)")
message("")
message("  Output folder: ", normalizePath(OUT_DIR))
message(strrep("=", 60))

# =============================================================
# 14. Spatial map x/y — RFP_positive (back-fill coords from CSVs)
#
#   center_x / center_y from source cell_metadata.csv joined by EntityID.
#   One panel per sample, free scales (slide2 vs slide4 coordinate systems differ).
#   RFP_neg plotted first (grey), RFP_pos on top (tomato) for visibility.
# =============================================================

message("\n=== Spatial map x/y ===")

SAMPLE_CSV_XY <- list(
  mock_6wpi = file.path("data", "slide4", "region_R2", "cell_metadata.csv"),
  LCMV_1wpi = file.path("data", "slide4", "region_R3", "cell_metadata.csv"),
  LCMV_3wpi = file.path("data", "slide2", "region_R1", "cell_metadata.csv"),
  LCMV_6wpi = file.path("data", "slide2", "region_R2", "cell_metadata.csv")
)

x_vec <- setNames(rep(NA_real_, ncol(obj)), colnames(obj))
y_vec <- setNames(rep(NA_real_, ncol(obj)), colnames(obj))

for (sname in names(SAMPLE_CSV_XY)) {
  csv_path <- SAMPLE_CSV_XY[[sname]]
  if (!file.exists(csv_path)) {
    message("  WARNING: CSV not found for ", sname, " -- skipping coords")
    next
  }
  obj_idx      <- which(obj$sample == sname)
  if (length(obj_idx) == 0) next
  obj_barcodes <- obj@meta.data$cell[obj_idx]

  csv_df <- read.csv(csv_path, stringsAsFactors = FALSE,
                     colClasses = c(EntityID = "character"))
  m <- match(obj_barcodes, csv_df$EntityID)
  x_vec[obj_idx] <- csv_df$center_x[m]
  y_vec[obj_idx] <- csv_df$center_y[m]
  message("  ", sname, ": ", sum(!is.na(m)), "/", length(obj_idx), " coords matched")
}

xy_df <- data.frame(
  x        = x_vec,
  y        = y_vec,
  rfp_pos  = obj$RFP_positive,
  sample   = factor(obj$sample, levels = sample_levels),
  stringsAsFactors = FALSE
)
xy_df <- xy_df[!is.na(xy_df$x) & !is.na(xy_df$y), ]

df_neg <- xy_df[xy_df$rfp_pos == "RFP_neg", ]
df_pos <- xy_df[xy_df$rfp_pos == "RFP_pos", ]

p_xy <- ggplot() +
  geom_point(data = df_neg, aes(x = x, y = y),
             color = "grey80", size = 0.15, alpha = 0.4) +
  geom_point(data = df_pos, aes(x = x, y = y),
             color = "tomato", size = 0.4, alpha = 0.8) +
  facet_wrap(~ sample, scales = "free", ncol = 2) +
  labs(title = "Cellules RFP+ dans l'espace tissulaire",
       subtitle = paste0("Rouge = RFP_pos (q99 mock raw = ",
                         format(round(p99_mock_raw, 0), big.mark = ","), ")  |  ",
                         "Gris = RFP_neg"),
       x = "x (µm)", y = "y (µm)") +
  theme_void(base_size = 10) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold", size = 10),
    plot.title       = element_text(face = "bold", size = 11, hjust = 0),
    plot.subtitle    = element_text(size = 9, color = "grey40", hjust = 0),
    aspect.ratio     = 1
  )

save_plot(p_xy,
          file.path(OUT_DIR, "fig_spatial_rfp_positive"),
          width = 14, height = 12)
