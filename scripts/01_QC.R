# =============================================================
# Script: 01_QC.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-04
# Description: Quality control of MERFISH data (slide4).
#              QC is performed per sample BEFORE merging.
#              Key metrics: transcripts/cell, genes/cell,
#              cell volume, blank gene detection, spatial
#              distribution of QC metrics.
# Input:  objects/00_mock_6wpi.rds
#         objects/00_LCMV_1wpi.rds
# Output: objects/01_mock_6wpi_qc.rds
#         objects/01_LCMV_1wpi_qc.rds
#         objects/01_slide4_merged_qc.rds
#         outputs/QC/  (figures + CSV summary)
# =============================================================


# -------------------------------------------------------
# 1. Libraries
# -------------------------------------------------------

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(patchwork)

OBJ_DIR <- "objects"
QC_DIR  <- "outputs/QC"
dir.create(QC_DIR, showWarnings = FALSE, recursive = TRUE)


# -------------------------------------------------------
# 2. Load objects
# -------------------------------------------------------

message("Loading objects...")
mock      <- readRDS(file.path(OBJ_DIR, "00_mock_6wpi.rds"))
lcmv_1wpi <- readRDS(file.path(OBJ_DIR, "00_LCMV_1wpi.rds"))

samples <- list(
  mock_6wpi = mock,
  LCMV_1wpi = lcmv_1wpi
)

# --- Early guard: check core QC columns exist ---
for (sample_name in names(samples)) {
  obj <- samples[[sample_name]]
  required_cols <- c("nCount_Vizgen", "nFeature_Vizgen")
  missing <- required_cols[!required_cols %in% colnames(obj@meta.data)]
  if (length(missing) > 0) {
    stop(sample_name, ": missing required columns: ", paste(missing, collapse = ", "),
         "\nCheck that LoadVizgen() ran correctly in 00_load_data.R")
  }
}
message("Core QC columns present in all samples.")


# -------------------------------------------------------
# 3. Visualize QC metrics — safe log10(x + 1) plots
# Note: blank gene QC is not applicable here — the 54 blank barcodes
# used by Vizgen for internal QC are not exported in cell_by_gene.csv
# and therefore not present in the Seurat object. Blank statistics
# are available in the Vizgen HTML report.
# -------------------------------------------------------

message("\n--- Visualizing QC metrics ---")

# Helper: extract spatial coordinates from FOV (Vizgen objects store coords
# in the FOV slot, not in metadata)
get_spatial_coords <- function(obj) {
  tryCatch({
    coords <- GetTissueCoordinates(obj, image = Images(obj)[1])
    if (!all(c("x", "y", "cell") %in% colnames(coords))) return(NULL)
    coords_indexed           <- coords
    rownames(coords_indexed) <- coords_indexed$cell
    # Align to object metadata order
    data.frame(
      x = coords_indexed[rownames(obj@meta.data), "x"],
      y = coords_indexed[rownames(obj@meta.data), "y"],
      row.names = rownames(obj@meta.data)
    )
  }, error = function(e) NULL)
}

# Helper: plot stain signal safely
plot_stain <- function(md, col, color, subtitle) {
  if (!col %in% colnames(md)) {
    return(ggplot() + labs(title = paste(col, "— not found")) + theme_void())
  }
  vals <- md[[col]]
  vals_safe <- vals[is.finite(vals) & vals > 0]
  if (length(vals_safe) == 0) {
    return(ggplot() + labs(title = paste(col, "— all zero or NA")) + theme_void())
  }
  ggplot(data.frame(x = log10(vals_safe + 1)), aes(x = x)) +
    geom_histogram(bins = 100, fill = color, color = "white", linewidth = 0.1) +
    labs(title = col, subtitle = subtitle,
         x = "log10(intensity + 1)", y = "N cells") +
    theme_classic()
}

for (sample_name in names(samples)) {
  obj <- samples[[sample_name]]
  md  <- obj@meta.data

  # 4a. Transcript count
  p_counts <- ggplot(md, aes(x = log10(nCount_Vizgen + 1))) +
    geom_histogram(bins = 100, fill = "#2E75B6", color = "white", linewidth = 0.1) +
    geom_vline(xintercept = log10(median(md$nCount_Vizgen) + 1),
               color = "red", linetype = "dashed") +
    labs(title = paste(sample_name, "— Transcripts per cell"),
         subtitle = paste0("Median: ", median(md$nCount_Vizgen),
                           " | N cells: ", nrow(md)),
         x = "log10(nCount_Vizgen + 1)", y = "N cells") +
    theme_classic()

  # 4b. Gene count
  p_genes <- ggplot(md, aes(x = nFeature_Vizgen)) +
    geom_histogram(bins = 100, fill = "#70AD47", color = "white", linewidth = 0.1) +
    geom_vline(xintercept = median(md$nFeature_Vizgen),
               color = "red", linetype = "dashed") +
    labs(title = paste(sample_name, "— Genes per cell"),
         subtitle = paste0("Median: ", median(md$nFeature_Vizgen)),
         x = "nFeature_Vizgen", y = "N cells") +
    theme_classic()

  # 4c. Cell volume
  if ("volume" %in% colnames(md)) {
    vol_vals <- md$volume[is.finite(md$volume) & md$volume > 0]
    p_vol <- ggplot(data.frame(v = vol_vals), aes(x = log10(v + 1))) +
      geom_histogram(bins = 100, fill = "#FFC000", color = "white", linewidth = 0.1) +
      geom_vline(xintercept = log10(median(vol_vals) + 1),
                 color = "red", linetype = "dashed") +
      labs(title = paste(sample_name, "— Cell volume"),
           subtitle = paste0("Median: ", round(median(vol_vals), 0), " µm³"),
           x = "log10(volume + 1)", y = "N cells") +
      theme_classic()
  } else {
    p_vol <- ggplot() + labs(title = "Volume not available") + theme_void()
  }

  # 4d. Counts vs genes scatter
  p_scatter <- ggplot(md, aes(x = log10(nCount_Vizgen + 1), y = nFeature_Vizgen)) +
    geom_point(size = 0.1, alpha = 0.3, color = "#404040") +
    labs(title = paste(sample_name, "— Counts vs Genes"),
         x = "log10(nCount_Vizgen + 1)", y = "nFeature_Vizgen") +
    theme_classic()

  # 4e. Stain signals
  p_rfp  <- plot_stain(md, "Anti.RFP_raw",  "#E74C3C", "TRM-contact cells")
  p_iba1 <- plot_stain(md, "Anti.IBA1_raw", "#9B59B6", "Microglia")
  p_cd8  <- plot_stain(md, "Anti.CD8_raw",  "#2980B9", "CD8+ T cells (TRM)")

  # 4f. Spatial distribution — coordinates from FOV slot
  spatial_coords <- get_spatial_coords(obj)
  if (!is.null(spatial_coords) && sum(!is.na(spatial_coords$x)) > 0) {
    md_spatial        <- md
    md_spatial$x      <- spatial_coords[rownames(md), "x"]
    md_spatial$y      <- spatial_coords[rownames(md), "y"]
    md_spatial        <- md_spatial[!is.na(md_spatial$x), ]

    p_spatial <- ggplot(md_spatial, aes(x = x, y = y,
                                         color = log10(nCount_Vizgen + 1))) +
      geom_point(size = 0.05, alpha = 0.5) +
      scale_color_viridis_c(option = "magma", name = "log10(counts+1)") +
      labs(title = paste(sample_name, "— Spatial transcript density"),
           x = "X (µm)", y = "Y (µm)") +
      theme_classic() +
      coord_equal()
  } else {
    message(sample_name, ": spatial coordinates not available — skipping spatial plot")
    p_spatial <- ggplot() + labs(title = "Centroids not available") + theme_void()
  }

  # Combine and save
  qc_panel <- (p_counts | p_genes | p_vol) /
              (p_scatter | plot_spacer() | plot_spacer()) /
              (p_rfp | p_iba1 | p_cd8) /
              p_spatial

  ggsave(
    filename = file.path(QC_DIR, paste0("QC_", sample_name, "_before_filtering.pdf")),
    plot     = qc_panel,
    width    = 18, height = 22
  )
  message("QC plot saved for: ", sample_name)
}


# -------------------------------------------------------
# 5. Summary statistics — saved as CSV
# -------------------------------------------------------

message("\n--- QC Summary ---")

summary_rows <- lapply(names(samples), function(sample_name) {
  md <- samples[[sample_name]]@meta.data
  data.frame(
    sample             = sample_name,
    n_cells            = nrow(md),
    median_counts      = median(md$nCount_Vizgen),
    median_genes       = median(md$nFeature_Vizgen),
    min_counts         = min(md$nCount_Vizgen),
    max_counts         = max(md$nCount_Vizgen),
    pct5_counts        = quantile(md$nCount_Vizgen, 0.05),
    pct95_counts       = quantile(md$nCount_Vizgen, 0.95),
    pct5_genes         = quantile(md$nFeature_Vizgen, 0.05),
    pct95_genes        = quantile(md$nFeature_Vizgen, 0.95),
    median_blank_ratio = if ("blank_ratio" %in% colnames(md))
                           round(median(md$blank_ratio), 5) else NA,
    median_volume      = if ("volume" %in% colnames(md))
                           round(median(md$volume, na.rm = TRUE), 1) else NA
  )
})

summary_df <- do.call(rbind, summary_rows)
print(summary_df)
write.csv(summary_df,
          file = file.path(QC_DIR, "QC_summary_before_filtering.csv"),
          row.names = FALSE)
message("Summary CSV saved.")


# -------------------------------------------------------
# 6. Apply filtering thresholds
# -------------------------------------------------------
# IMPORTANT: inspect QC plots BEFORE running this section.
# Values below are starting points — adjust based on your distributions.
#
# MERFISH-specific notes:
# - Median ~16-21 transcripts/cell here — much lower than scRNAseq
#   Do NOT use scRNAseq thresholds (200+ minimum) here
# - Volume filter removes segmentation artifacts

message("\n--- Applying QC filters ---")

qc_thresholds <- list(
  nCount_min      = 5,      # debris / empty segmentations
  nCount_max      = 500,    # potential doublets / merged cells
  nFeature_min    = 5,      # very low quality cells
  volume_min      = 50,     # µm³ — segmentation noise
  volume_max      = 10000,  # µm³ — merged cells
  blank_ratio_max = 0.1     # > 10% blank = poor detection quality
)

message("Thresholds applied:")
print(as.data.frame(qc_thresholds))

for (sample_name in names(samples)) {
  obj      <- samples[[sample_name]]
  md       <- obj@meta.data
  n_before <- ncol(obj)

  keep <- rep(TRUE, nrow(md))

  keep <- keep & (md$nCount_Vizgen   >= qc_thresholds$nCount_min)
  keep <- keep & (md$nCount_Vizgen   <= qc_thresholds$nCount_max)
  keep <- keep & (md$nFeature_Vizgen >= qc_thresholds$nFeature_min)

  if ("volume" %in% colnames(md)) {
    keep <- keep & !is.na(md$volume) &
            (md$volume >= qc_thresholds$volume_min) &
            (md$volume <= qc_thresholds$volume_max)
  }

  if ("blank_ratio" %in% colnames(md)) {
    keep <- keep & (md$blank_ratio <= qc_thresholds$blank_ratio_max)
  }

  # Safety: NA -> FALSE
  keep[is.na(keep)] <- FALSE

  obj_filtered <- obj[, keep]
  n_after      <- ncol(obj_filtered)

  message(sample_name, ": ", n_before, " -> ", n_after,
          " cells (", round((1 - n_after / n_before) * 100, 1), "% removed)")

  samples[[sample_name]] <- obj_filtered
}


# -------------------------------------------------------
# 7. Before / After comparison plots
# -------------------------------------------------------

for (sample_name in names(samples)) {
  obj_before <- if (sample_name == "mock_6wpi") mock else lcmv_1wpi
  obj_after  <- samples[[sample_name]]

  df_both <- rbind(
    data.frame(counts = obj_before@meta.data$nCount_Vizgen, status = "Before"),
    data.frame(counts = obj_after@meta.data$nCount_Vizgen,  status = "After")
  )
  df_both$status <- factor(df_both$status, levels = c("Before", "After"))

  p_compare <- ggplot(df_both, aes(x = log10(counts + 1), fill = status)) +
    geom_histogram(bins = 80, alpha = 0.6, position = "identity",
                   color = "white", linewidth = 0.1) +
    scale_fill_manual(values = c("Before" = "#95A5A6", "After" = "#2E75B6")) +
    labs(title    = paste(sample_name, "— Before vs After filtering"),
         subtitle = paste0("Before: ", nrow(obj_before@meta.data),
                           " | After: ", nrow(obj_after@meta.data), " cells"),
         x = "log10(nCount_Vizgen + 1)", y = "N cells") +
    theme_classic()

  ggsave(
    filename = file.path(QC_DIR, paste0("QC_", sample_name, "_before_vs_after.pdf")),
    plot     = p_compare,
    width    = 7, height = 5
  )
}
message("Before/after comparison plots saved.")


# -------------------------------------------------------
# 8. Post-filtering summary
# -------------------------------------------------------

post_rows <- lapply(names(samples), function(s) {
  md <- samples[[s]]@meta.data
  data.frame(
    sample        = s,
    n_cells_after = nrow(md),
    median_counts = median(md$nCount_Vizgen),
    median_genes  = median(md$nFeature_Vizgen)
  )
})
post_df <- do.call(rbind, post_rows)
write.csv(post_df,
          file = file.path(QC_DIR, "QC_summary_after_filtering.csv"),
          row.names = FALSE)
message("Post-filtering summary saved.")


# -------------------------------------------------------
# 9. Merge and save
# -------------------------------------------------------

message("\nMerging filtered objects...")

slide4_merged_qc <- merge(
  samples$mock_6wpi,
  y            = samples$LCMV_1wpi,
  add.cell.ids = c("mock_6wpi", "LCMV_1wpi"),
  project      = "LCMV_MERFISH_slide4"
)

message("Merged: ", ncol(slide4_merged_qc), " total cells")
print(table(slide4_merged_qc@meta.data$sample))

saveRDS(samples$mock_6wpi,   file.path(OBJ_DIR, "01_mock_6wpi_qc.rds"))
saveRDS(samples$LCMV_1wpi,   file.path(OBJ_DIR, "01_LCMV_1wpi_qc.rds"))
saveRDS(slide4_merged_qc,    file.path(OBJ_DIR, "01_slide4_merged_qc.rds"))

message("\nDone. Objects saved to: ", OBJ_DIR)
message("QC figures and CSVs saved to: ", QC_DIR)
message("\nIMPORTANT: review QC plots before proceeding.")
message("Adjust thresholds in section 6 if needed and re-run.")
message("\nNext step: 02_normalization.R")