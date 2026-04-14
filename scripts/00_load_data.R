# =============================================================
# Script: 00_load_data.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-04
# Description: Load raw MERFISH data from Slide4 (most recent run),
#              correct sample labels (R2 = Mock_6wpi, R3 = LCMV_1wpi),
#              add auxiliary stain metadata (RFP, IBA1, CD8),
#              merge into a single Seurat object and save.
# Input:  Raw Vizgen files per region:
#           - cell_by_gene.csv
#           - cell_metadata.csv
#           - cell_boundaries.parquet
# Output: objects/00_slide4_merged.rds
# Notes:  - Requires custom Seurat fork for Vizgen compatibility
#           (alikhuseynov/seurat, ref = 'vizgen_seurat5')
#         - Slide4 region assignment:
#             R2 = Mock_6wpi  (74,129 cells, median 16 transcripts)
#             R3 = LCMV_1wpi  (81,095 cells, median 21 transcripts)
# =============================================================


# -------------------------------------------------------
# 0. Install custom Seurat if not already installed
# -------------------------------------------------------

# Set to TRUE only if you intentionally want to reinstall Seurat from GitHub.

INSTALL_SEURAT <- FALSE   # ← set to TRUE only on first run

if (INSTALL_SEURAT) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  if (requireNamespace("Seurat", quietly = TRUE)) remove.packages("Seurat")
  remotes::install_github("alikhuseynov/seurat", ref = "vizgen_seurat5")
  message("Seurat custom installed. Restart R before continuing.")
  stop("Restart R now, then re-run the script with INSTALL_SEURAT = FALSE.")
}
# install.packages(c("sfarrow", "tidyverse", "sf"), type = "binary")
# sapply(c("sfarrow", "tidyverse", "sf"), requireNamespace, quietly = TRUE)

# -------------------------------------------------------
# 1. Libraries
# -------------------------------------------------------

library(Seurat)
library(SeuratObject)
library(data.table)
library(arrow)        # for reading .parquet files
library(dplyr)
library(ggplot2)
library(future)
library(BiocParallel)
library(progressr)
library(magrittr)
library(sf)
# Allow larger global variables for parallel processing
options(future.globals.maxSize = 8000 * 1024^2)


# -------------------------------------------------------
# 2. Paths
# -------------------------------------------------------

# Root directory where Globus data was downloaded
DATA_ROOT <- "data/slide4"

# Output directory for saved Seurat objects
OBJ_DIR <- "objects"
dir.create(OBJ_DIR, showWarnings = FALSE)

# MERFISH parameters
mol.type    <- "microns"
coord.space <- "micron"
z.stack     <- 3L   # z-plane used for segmentation (Zoé's choice)


# -------------------------------------------------------
# 3. Load Region R2 — Mock_6wpi
#    (CORRECTED: R2 is Mock, not LCMV_1wpi as in Zoé's original code)
# -------------------------------------------------------

message("Loading R2 — Mock_6wpi...")

mock <- LoadVizgen(
  data.dir              = file.path(DATA_ROOT, "region_R2"),
  fov                   = "mock",
  assay                 = "Vizgen",
  metadata              = c("volume", "fov"),
  mol.type              = mol.type,
  type                  = c("segmentations", "centroids"),
  z                     = z.stack,
  add.zIndex            = TRUE,
  update.object         = TRUE,
  use.BiocParallel      = TRUE,
  workers.MulticoreParam = 10,
  verbose               = TRUE
)

# Add sample-level metadata
mock@meta.data$section   <- "slide4_R2"
mock@meta.data$sample    <- "mock_6wpi"
mock@meta.data$condition <- "mock"
mock@meta.data$timepoint <- "6wpi"
mock@meta.data$slide     <- "slide4"
mock@meta.data$cell      <- rownames(mock@meta.data)

# Add auxiliary stain intensities from cell_metadata.csv
cell_meta_mock <- read.csv(
  file.path(DATA_ROOT, "region_R2", "cell_metadata.csv"),
  stringsAsFactors = FALSE
)

# Align row order to Seurat object
if ("cell_id" %in% colnames(cell_meta_mock)) {
  cell_meta_mock <- cell_meta_mock[
    match(mock@meta.data$cell, cell_meta_mock$cell_id), 
  ]
}

# Add stain signals — key for downstream analysis:
# Anti-RFP  : marks cells in direct contact with TRM (CRE-lox labeling system)
# Anti-IBA1 : marks all microglia
# Anti-CD8  : marks CD8+ T cells (TRM)
mock[["Anti.RFP_raw"]]  <- cell_meta_mock$Anti.RFP_raw
mock[["Anti.IBA1_raw"]] <- cell_meta_mock$Anti.IBA1_raw
mock[["Anti.CD8_raw"]]  <- cell_meta_mock$Anti.CD8_raw
mock[["DAPI_high_pass"]]  <- cell_meta_mock$DAPI_high_pass
mock[["PolyT_high_pass"]] <- cell_meta_mock$PolyT_high_pass

message("Mock_6wpi loaded: ", ncol(mock), " cells")


# -------------------------------------------------------
# 4. Load Region R3 — LCMV_1wpi
# -------------------------------------------------------

message("Loading R3 — LCMV_1wpi...")

lcmv_1wpi <- LoadVizgen(
  data.dir              = file.path(DATA_ROOT, "region_R3"),
  fov                   = "lcmv_1wpi",
  assay                 = "Vizgen",
  metadata              = c("volume", "fov"),
  mol.type              = mol.type,
  type                  = c("segmentations", "centroids"),
  z                     = z.stack,
  add.zIndex            = TRUE,
  update.object         = TRUE,
  use.BiocParallel      = TRUE,
  workers.MulticoreParam = 10,
  verbose               = TRUE
)

# Add sample-level metadata
lcmv_1wpi@meta.data$section   <- "slide4_R3"
lcmv_1wpi@meta.data$sample    <- "LCMV_1wpi"
lcmv_1wpi@meta.data$condition <- "LCMV"
lcmv_1wpi@meta.data$timepoint <- "1wpi"
lcmv_1wpi@meta.data$slide     <- "slide4"
lcmv_1wpi@meta.data$cell      <- rownames(lcmv_1wpi@meta.data)

# Add auxiliary stain intensities
cell_meta_lcmv1 <- read.csv(
  file.path(DATA_ROOT, "region_R3", "cell_metadata.csv"),
  stringsAsFactors = FALSE
)

if ("cell_id" %in% colnames(cell_meta_lcmv1)) {
  cell_meta_lcmv1 <- cell_meta_lcmv1[
    match(lcmv_1wpi@meta.data$cell, cell_meta_lcmv1$cell_id), 
  ]
}

lcmv_1wpi[["Anti.RFP_raw"]]  <- cell_meta_lcmv1$Anti.RFP_raw
lcmv_1wpi[["Anti.IBA1_raw"]] <- cell_meta_lcmv1$Anti.IBA1_raw
lcmv_1wpi[["Anti.CD8_raw"]]  <- cell_meta_lcmv1$Anti.CD8_raw
lcmv_1wpi[["DAPI_high_pass"]]  <- cell_meta_lcmv1$DAPI_high_pass
lcmv_1wpi[["PolyT_high_pass"]] <- cell_meta_lcmv1$PolyT_high_pass

message("LCMV_1wpi loaded: ", ncol(lcmv_1wpi), " cells")


# -------------------------------------------------------
# 5. Sanity checks before merging
# -------------------------------------------------------

message("\n--- Sanity checks ---")

# Check cell counts match report
message("Expected from Vizgen report:")
message("  Mock_6wpi  : 74,129 cells")
message("  LCMV_1wpi  : 81,095 cells")
message("Loaded:")
message("  Mock_6wpi  : ", ncol(mock), " cells")
message("  LCMV_1wpi  : ", ncol(lcmv_1wpi), " cells")

# Check that gene panels are identical between regions
genes_mock    <- rownames(mock)
genes_lcmv1   <- rownames(lcmv_1wpi)
if (!identical(sort(genes_mock), sort(genes_lcmv1))) {
  warning("Gene panels differ between regions! Check input files.")
} else {
  message("Gene panels match between regions: ", length(genes_mock), " genes")
}

# Check metadata columns are complete (no NAs in stain signals)
for (col in c("Anti.RFP_raw", "Anti.IBA1_raw", "Anti.CD8_raw")) {
  n_na_mock  <- sum(is.na(mock@meta.data[[col]]))
  n_na_lcmv1 <- sum(is.na(lcmv_1wpi@meta.data[[col]]))
  message(col, " — NAs in mock: ", n_na_mock, 
          " | NAs in lcmv_1wpi: ", n_na_lcmv1)
}

# Check sample labels are correct
message("\nSample labels check:")
print(table(mock@meta.data$sample))
print(table(lcmv_1wpi@meta.data$sample))


# -------------------------------------------------------
# 6. Merge into a single Seurat object
# -------------------------------------------------------

message("\nMerging objects...")

slide4_merged <- merge(
  mock,
  y          = lcmv_1wpi,
  add.cell.ids = c("mock_6wpi", "LCMV_1wpi"),
  project    = "LCMV_MERFISH_slide4"
)

# Verify merge
message("Merged object: ", ncol(slide4_merged), " total cells")
message("Sample distribution:")
print(table(slide4_merged@meta.data$sample))
print(table(slide4_merged@meta.data$condition))


# -------------------------------------------------------
# 7. Save
# -------------------------------------------------------

saveRDS(mock,          file = file.path(OBJ_DIR, "00_mock_6wpi.rds"))
saveRDS(lcmv_1wpi,    file = file.path(OBJ_DIR, "00_LCMV_1wpi.rds"))
saveRDS(slide4_merged, file = file.path(OBJ_DIR, "00_slide4_merged.rds"))

message("\nDone. Objects saved to: ", OBJ_DIR)
