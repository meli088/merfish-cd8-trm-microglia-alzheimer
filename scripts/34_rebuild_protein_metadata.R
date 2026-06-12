#!/usr/bin/env Rscript
# =============================================================
# Script: 34_rebuild_protein_metadata.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-06
#
# Purpose:
#   The protein metadata columns (Anti.RFP_raw, Anti.IBA1_raw, Anti.CD8_raw,
#   DAPI_high_pass, Anti.RFP_high_pass, PolyT_high_pass) were incorrectly
#   assigned in 00_load_data.R: the CSV was joined on "cell_id" which does
#   not exist in the Vizgen output, so the if() guard was never triggered
#   and values were assigned in raw CSV row order (= wrong cell order).
#   Every protein value was assigned to the wrong cell.
#
#   This script:
#     1. Loads objects/02_all_normalized.rds  (last clean Seurat in pipeline;
#        UMAP + Harmony + SCT are correct — gene-based, unaffected by bug)
#     2. Drops all 6 corrupted protein columns
#     3. Re-imports them from source CSVs using the correct EntityID join
#        (same approach as 35_rfp_by_celltype.R section 5b, confirmed 100% match)
#     4. Transfers CellType from 05_joint_annotated_lam02_res09.rds (SPE)
#     5. Saves objects/03_rfp_analysis.rds  (NEW — does NOT overwrite anything)
#
#   35_rfp_by_celltype.R then loads 03_rfp_analysis.rds directly.
#   No back-fill needed at runtime.
#
# Input:
#   objects/02_all_normalized.rds
#   objects/05_joint_annotated_lam02_res09.rds
#   data/slide4/region_R2/cell_metadata.csv  (mock_6wpi)
#   data/slide4/region_R3/cell_metadata.csv  (LCMV_1wpi)
#   data/slide2/region_R1/cell_metadata.csv  (LCMV_3wpi)
#   data/slide2/region_R2/cell_metadata.csv  (LCMV_6wpi)
#
# Output:
#   objects/03_rfp_analysis.rds
# =============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
})

setwd(normalizePath("."))

OBJ_IN  <- file.path("objects", "02_all_normalized.rds")
OBJ_SPE <- file.path("objects", "05_joint_annotated_lam02_res09.rds")
OBJ_OUT <- file.path("objects", "03_rfp_analysis.rds")

SAMPLE_CSV <- list(
  mock_6wpi = file.path("data", "slide4", "region_R2", "cell_metadata.csv"),
  LCMV_1wpi = file.path("data", "slide4", "region_R3", "cell_metadata.csv"),
  LCMV_3wpi = file.path("data", "slide2", "region_R1", "cell_metadata.csv"),
  LCMV_6wpi = file.path("data", "slide2", "region_R2", "cell_metadata.csv")
)

# Protein columns to drop (corrupted) and re-import
PROTEIN_COLS <- c("Anti.RFP_raw", "Anti.IBA1_raw", "Anti.CD8_raw",
                  "DAPI_high_pass", "Anti.RFP_high_pass", "PolyT_high_pass")

# CSV column name (after read.csv dot-substitution) -> R metadata name
# read.csv() converts dashes to dots: "Anti-RFP_raw" -> "Anti.RFP_raw"
CSV_TO_META <- c(
  "Anti.RFP_raw"        = "Anti.RFP_raw",
  "Anti.IBA1_raw"       = "Anti.IBA1_raw",
  "Anti.CD8_raw"        = "Anti.CD8_raw",
  "DAPI_high_pass"      = "DAPI_high_pass",
  "Anti.RFP_high_pass"  = "Anti.RFP_high_pass",
  "PolyT_high_pass"     = "PolyT_high_pass"
)

# =============================================================
# 1. Load base Seurat object
# =============================================================

message("\n=== Loading 02_all_normalized.rds ===")
obj <- readRDS(OBJ_IN)
message("  Cells: ", ncol(obj))
message("  Reductions: ", paste(names(obj@reductions), collapse = ", "))
message("  Samples: ", paste(sort(unique(obj$sample)), collapse = ", "))

if (!"cell" %in% colnames(obj@meta.data)) {
  stop("'cell' column not found — needed as EntityID join key.")
}

# =============================================================
# 2. Drop corrupted protein columns
# =============================================================

message("\n=== Dropping corrupted protein columns ===")
cols_present <- intersect(PROTEIN_COLS, colnames(obj@meta.data))
message("  Found: ", paste(cols_present, collapse = ", "))
for (col in cols_present) {
  obj@meta.data[[col]] <- NULL
}
message("  Dropped ", length(cols_present), " columns.")

# =============================================================
# 3. Re-import protein values from source CSVs (correct EntityID join)
# =============================================================

message("\n=== Re-importing protein metadata from source CSVs ===")

missing_csv <- names(SAMPLE_CSV)[!file.exists(unlist(SAMPLE_CSV))]
if (length(missing_csv) > 0) {
  stop("Source CSV(s) not found for: ", paste(missing_csv, collapse = ", "))
}

# Initialise vectors (one per protein column)
new_vals <- lapply(PROTEIN_COLS, function(col) {
  setNames(rep(NA_real_, ncol(obj)), colnames(obj))
})
names(new_vals) <- PROTEIN_COLS

for (sname in names(SAMPLE_CSV)) {
  csv_path     <- SAMPLE_CSV[[sname]]
  obj_idx      <- which(obj$sample == sname)

  if (length(obj_idx) == 0) {
    message("  ", sname, ": 0 cells in object -- skipping")
    next
  }
  obj_barcodes <- obj@meta.data$cell[obj_idx]

  csv_df <- read.csv(csv_path, stringsAsFactors = FALSE,
                     colClasses = c(EntityID = "character"))

  # Match cells
  m           <- match(obj_barcodes, csv_df$EntityID)
  n_matched   <- sum(!is.na(m))
  pct_matched <- round(100 * n_matched / length(obj_idx), 2)
  message("  ", sname, ": ", n_matched, "/", length(obj_idx),
          " (", pct_matched, "%) cells matched via EntityID")

  if (pct_matched < 95) {
    stop("Match rate too low (", pct_matched, "%) for ", sname,
         ". Check that obj$cell contains raw EntityID strings.")
  }

  # Assign each protein column
  for (csv_col in names(CSV_TO_META)) {
    meta_col <- CSV_TO_META[[csv_col]]
    if (!csv_col %in% colnames(csv_df)) {
      message("    WARNING: '", csv_col, "' not found in CSV for ", sname, " -- skipping")
      next
    }
    new_vals[[meta_col]][obj_idx] <- csv_df[[csv_col]][m]
  }
}

# Add back to object
for (col in PROTEIN_COLS) {
  n_na <- sum(is.na(new_vals[[col]]))
  if (n_na > 0) message("  WARNING: ", n_na, " NAs in ", col)
  obj[[col]] <- new_vals[[col]]
}
message("  Protein columns re-imported: ", paste(PROTEIN_COLS, collapse = ", "))

# Quick sanity check: per-sample summary of Anti.RFP_raw
message("\n  Anti.RFP_raw per-sample summary (should differ between mock and LCMV):")
rfp_summary <- tapply(obj$Anti.RFP_raw, obj$sample, function(x) {
  round(c(median = median(x, na.rm = TRUE), mean = mean(x, na.rm = TRUE),
          q99 = quantile(x, 0.99, na.rm = TRUE)), 0)
})
print(do.call(rbind, rfp_summary))

# =============================================================
# 4. Transfer CellType from SPE
# =============================================================

message("\n=== Transferring CellType from SPE ===")
spe <- readRDS(OBJ_SPE)
message("  SPE cells: ", ncol(spe))

cd <- as.data.frame(colData(spe))
if (!"annotation" %in% colnames(cd)) {
  stop("'annotation' column not found in SPE colData.")
}

anno_vec  <- setNames(as.character(cd$annotation), colnames(spe))
matched   <- intersect(colnames(obj), names(anno_vec))
unmatched <- ncol(obj) - length(matched)

obj$CellType <- NA_character_
obj$CellType[matched] <- anno_vec[matched]
obj$CellType[is.na(obj$CellType)] <- "Non annote"

message("  Matched  : ", length(matched), " / ", ncol(obj), " cells")
message("  Unmatched: ", unmatched, " -> 'Non annote'")
message("  CellType table:")
print(sort(table(obj$CellType), decreasing = TRUE))

rm(spe, cd, anno_vec); gc(verbose = FALSE)

# =============================================================
# 5. Save
# =============================================================

message("\n=== Saving ===")
message("  Output: ", OBJ_OUT)
saveRDS(obj, OBJ_OUT)
message("  Done. Object size: ",
        round(file.size(OBJ_OUT) / 1024^3, 2), " GB")
message("\n  -> Load in 35_rfp_by_celltype.R with:")
message("     OBJ_SEURAT <- \"objects/03_rfp_analysis.rds\"")