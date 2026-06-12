# =============================================================
# Script: 00_load_data.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-04
# Description: Load raw MERFISH data from Slide4 and Slide2,
#              correct sample labels, add stain metadata,
#              merge into a single Seurat object and save.
# Slide 4: R2 = Mock_6wpi  |  R3 = LCMV_1wpi
# Slide 2: R1 = LCMV_3wpi  |  R2 = LCMV_6wpi
# =============================================================


# -------------------------------------------------------
# 0. Install custom Seurat if not already installed
# -------------------------------------------------------
INSTALL_SEURAT <- FALSE
if (INSTALL_SEURAT) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  if (requireNamespace("Seurat", quietly = TRUE)) remove.packages("Seurat")
  remotes::install_github("alikhuseynov/seurat", ref = "vizgen_seurat5")
  message("Seurat custom installed. Restart R before continuing.")
  stop("Restart R now, then re-run the script with INSTALL_SEURAT = FALSE.")
}

"LoadVizgen" %in% getNamespaceExports("Seurat")
# TRUE

# -------------------------------------------------------
# 1. Libraries
# -------------------------------------------------------

library("hdf5r")
packageVersion("Seurat")
library(Seurat)
library(SeuratObject)
library(data.table)
library(arrow)
library(dplyr)
library(ggplot2)
library(future)
library(BiocParallel)
library(progressr)
library(magrittr)
library(sf)
options(future.globals.maxSize = 8000 * 1024^2)


# -------------------------------------------------------
# 2. Paths
# -------------------------------------------------------
DATA_SLIDE4 <- "data/slide4"
DATA_SLIDE2 <- "data/slide2"

OBJ_DIR <- "objects"
dir.create(OBJ_DIR, showWarnings = FALSE)

mol.type    <- "microns"
coord.space <- "micron"
z.stack     <- 3L
WORKERS     <- 10L


# -------------------------------------------------------
# 3. Load Region R2 — Mock_6wpi (Slide 4)
# -------------------------------------------------------
message("Loading R2 — Mock_6wpi...")

mock <- LoadVizgen(
  data.dir              = file.path(DATA_SLIDE4, "region_R2"),
  fov                   = "mock",
  assay                 = "Vizgen",
  metadata              = c("volume", "fov"),
  mol.type              = mol.type,
  type                  = c("segmentations", "centroids"),
  z                     = z.stack,
  add.zIndex            = TRUE,
  update.object         = TRUE,
  use.BiocParallel      = TRUE,
  workers.MulticoreParam = WORKERS,
  verbose               = TRUE
)

mock@meta.data$section   <- "slide4_R2"
mock@meta.data$sample    <- "mock_6wpi"
mock@meta.data$condition <- "mock"
mock@meta.data$timepoint <- "6wpi"
mock@meta.data$slide     <- "slide4"
mock@meta.data$cell      <- rownames(mock@meta.data)

cell_meta_mock <- read.csv(
  file.path(DATA_SLIDE4, "region_R2", "cell_metadata.csv"),
  stringsAsFactors = FALSE,
  colClasses = c(EntityID = "character")
)
cell_meta_mock <- cell_meta_mock[match(mock@meta.data$cell, cell_meta_mock$EntityID), ]
mock[["Anti.RFP_raw"]]    <- cell_meta_mock$Anti.RFP_raw
mock[["Anti.IBA1_raw"]]   <- cell_meta_mock$Anti.IBA1_raw
mock[["Anti.CD8_raw"]]    <- cell_meta_mock$Anti.CD8_raw
mock[["DAPI_high_pass"]]      <- cell_meta_mock$DAPI_high_pass
mock[["Anti.RFP_high_pass"]] <- cell_meta_mock$Anti.RFP_high_pass
mock[["PolyT_high_pass"]]     <- cell_meta_mock$PolyT_high_pass

message("Mock_6wpi loaded: ", ncol(mock), " cells")


# -------------------------------------------------------
# 4. Load Region R3 — LCMV_1wpi (Slide 4)
# -------------------------------------------------------
message("Loading R3 — LCMV_1wpi...")

lcmv_1wpi <- LoadVizgen(
  data.dir              = file.path(DATA_SLIDE4, "region_R3"),
  fov                   = "lcmv_1wpi",
  assay                 = "Vizgen",
  metadata              = c("volume", "fov"),
  mol.type              = mol.type,
  type                  = c("segmentations", "centroids"),
  z                     = z.stack,
  add.zIndex            = TRUE,
  update.object         = TRUE,
  use.BiocParallel      = TRUE,
  workers.MulticoreParam = WORKERS,
  verbose               = TRUE
)

lcmv_1wpi@meta.data$section   <- "slide4_R3"
lcmv_1wpi@meta.data$sample    <- "LCMV_1wpi"
lcmv_1wpi@meta.data$condition <- "LCMV"
lcmv_1wpi@meta.data$timepoint <- "1wpi"
lcmv_1wpi@meta.data$slide     <- "slide4"
lcmv_1wpi@meta.data$cell      <- rownames(lcmv_1wpi@meta.data)

cell_meta_lcmv1 <- read.csv(
  file.path(DATA_SLIDE4, "region_R3", "cell_metadata.csv"),
  stringsAsFactors = FALSE,
  colClasses = c(EntityID = "character")
)
cell_meta_lcmv1 <- cell_meta_lcmv1[match(lcmv_1wpi@meta.data$cell, cell_meta_lcmv1$EntityID), ]
lcmv_1wpi[["Anti.RFP_raw"]]    <- cell_meta_lcmv1$Anti.RFP_raw
lcmv_1wpi[["Anti.IBA1_raw"]]   <- cell_meta_lcmv1$Anti.IBA1_raw
lcmv_1wpi[["Anti.CD8_raw"]]    <- cell_meta_lcmv1$Anti.CD8_raw
lcmv_1wpi[["DAPI_high_pass"]]      <- cell_meta_lcmv1$DAPI_high_pass
lcmv_1wpi[["Anti.RFP_high_pass"]] <- cell_meta_lcmv1$Anti.RFP_high_pass
lcmv_1wpi[["PolyT_high_pass"]]     <- cell_meta_lcmv1$PolyT_high_pass

message("LCMV_1wpi loaded: ", ncol(lcmv_1wpi), " cells")


# -------------------------------------------------------
# 5. Load Region R1 — LCMV_3wpi (Slide 2)
# -------------------------------------------------------
message("Loading Slide2 R1 — LCMV_3wpi...")

lcmv_3wpi <- LoadVizgen(
  data.dir              = file.path(DATA_SLIDE2, "region_R1"),
  fov                   = "lcmv_3wpi",
  assay                 = "Vizgen",
  metadata              = c("volume", "fov"),
  mol.type              = mol.type,
  type                  = c("segmentations", "centroids"),
  z                     = z.stack,
  add.zIndex            = TRUE,
  update.object         = TRUE,
  use.BiocParallel      = TRUE,
  workers.MulticoreParam = WORKERS,
  verbose               = TRUE
)

lcmv_3wpi@meta.data$section   <- "slide2_R1"
lcmv_3wpi@meta.data$sample    <- "LCMV_3wpi"
lcmv_3wpi@meta.data$condition <- "LCMV"
lcmv_3wpi@meta.data$timepoint <- "3wpi"
lcmv_3wpi@meta.data$slide     <- "slide2"
lcmv_3wpi@meta.data$cell      <- rownames(lcmv_3wpi@meta.data)

cell_meta_lcmv3 <- read.csv(
  file.path(DATA_SLIDE2, "region_R1", "cell_metadata.csv"),
  stringsAsFactors = FALSE,
  colClasses = c(EntityID = "character")
)
cell_meta_lcmv3 <- cell_meta_lcmv3[match(lcmv_3wpi@meta.data$cell, cell_meta_lcmv3$EntityID), ]
lcmv_3wpi[["Anti.RFP_raw"]]    <- cell_meta_lcmv3$Anti.RFP_raw
lcmv_3wpi[["Anti.IBA1_raw"]]   <- cell_meta_lcmv3$Anti.IBA1_raw
lcmv_3wpi[["Anti.CD8_raw"]]    <- cell_meta_lcmv3$Anti.CD8_raw
lcmv_3wpi[["DAPI_high_pass"]]      <- cell_meta_lcmv3$DAPI_high_pass
lcmv_3wpi[["Anti.RFP_high_pass"]] <- cell_meta_lcmv3$Anti.RFP_high_pass
lcmv_3wpi[["PolyT_high_pass"]]     <- cell_meta_lcmv3$PolyT_high_pass

message("LCMV_3wpi loaded: ", ncol(lcmv_3wpi), " cells")


# -------------------------------------------------------
# 6. Load Region R2 — LCMV_6wpi (Slide 2)
# -------------------------------------------------------
message("Loading Slide2 R2 — LCMV_6wpi...")

lcmv_6wpi <- LoadVizgen(
  data.dir              = file.path(DATA_SLIDE2, "region_R2"),
  fov                   = "lcmv_6wpi",
  assay                 = "Vizgen",
  metadata              = c("volume", "fov"),
  mol.type              = mol.type,
  type                  = c("segmentations", "centroids"),
  z                     = z.stack,
  add.zIndex            = TRUE,
  update.object         = TRUE,
  use.BiocParallel      = TRUE,
  workers.MulticoreParam = WORKERS,
  verbose               = TRUE
)

lcmv_6wpi@meta.data$section   <- "slide2_R2"
lcmv_6wpi@meta.data$sample    <- "LCMV_6wpi"
lcmv_6wpi@meta.data$condition <- "LCMV"
lcmv_6wpi@meta.data$timepoint <- "6wpi"
lcmv_6wpi@meta.data$slide     <- "slide2"
lcmv_6wpi@meta.data$cell      <- rownames(lcmv_6wpi@meta.data)

cell_meta_lcmv6 <- read.csv(
  file.path(DATA_SLIDE2, "region_R2", "cell_metadata.csv"),
  stringsAsFactors = FALSE,
  colClasses = c(EntityID = "character")
)
cell_meta_lcmv6 <- cell_meta_lcmv6[match(lcmv_6wpi@meta.data$cell, cell_meta_lcmv6$EntityID), ]
lcmv_6wpi[["Anti.RFP_raw"]]    <- cell_meta_lcmv6$Anti.RFP_raw
lcmv_6wpi[["Anti.IBA1_raw"]]   <- cell_meta_lcmv6$Anti.IBA1_raw
lcmv_6wpi[["Anti.CD8_raw"]]    <- cell_meta_lcmv6$Anti.CD8_raw
lcmv_6wpi[["DAPI_high_pass"]]      <- cell_meta_lcmv6$DAPI_high_pass
lcmv_6wpi[["Anti.RFP_high_pass"]] <- cell_meta_lcmv6$Anti.RFP_high_pass
lcmv_6wpi[["PolyT_high_pass"]]     <- cell_meta_lcmv6$PolyT_high_pass

message("LCMV_6wpi loaded: ", ncol(lcmv_6wpi), " cells")


# -------------------------------------------------------
# 7. Sanity checks
# -------------------------------------------------------
message("\n--- Sanity checks ---")

message("Loaded:")
message("  mock_6wpi  : ", ncol(mock),      " cells")
message("  LCMV_1wpi  : ", ncol(lcmv_1wpi), " cells")
message("  LCMV_3wpi  : ", ncol(lcmv_3wpi), " cells")
message("  LCMV_6wpi  : ", ncol(lcmv_6wpi), " cells")

# Panels identiques entre les 4 ?
genes_mock  <- rownames(mock)
genes_lcmv1 <- rownames(lcmv_1wpi)
genes_lcmv3 <- rownames(lcmv_3wpi)
genes_lcmv6 <- rownames(lcmv_6wpi)

if (any(c(
  !identical(sort(genes_mock), sort(genes_lcmv1)),
  !identical(sort(genes_mock), sort(genes_lcmv3)),
  !identical(sort(genes_mock), sort(genes_lcmv6))
))) {
  warning("Gene panels differ between samples! Check input files.")
} else {
  message("Gene panels match across all 4 samples: ", length(genes_mock), " genes")
}

# NAs dans les stains
for (col in c("Anti.RFP_raw", "Anti.IBA1_raw", "Anti.CD8_raw")) {
  for (obj_name in c("mock", "lcmv_1wpi", "lcmv_3wpi", "lcmv_6wpi")) {
    obj_tmp <- get(obj_name)
    n_na <- sum(is.na(obj_tmp@meta.data[[col]]))
    if (n_na > 0) message(" NAs dans ", obj_name, " — ", col, " : ", n_na)
  }
}

message("\nSample labels check:")
print(table(mock@meta.data$sample))
print(table(lcmv_1wpi@meta.data$sample))
print(table(lcmv_3wpi@meta.data$sample))
print(table(lcmv_6wpi@meta.data$sample))


# -------------------------------------------------------
# 8. Merge — 4 échantillons
# -------------------------------------------------------
message("\nMerging all objects...")

all_merged <- merge(
  mock,
  y            = list(lcmv_1wpi, lcmv_3wpi, lcmv_6wpi),
  add.cell.ids = c("mock_6wpi", "LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi"),
  project      = "LCMV_MERFISH_all"
)

message("Merged object: ", ncol(all_merged), " total cells")
message("Sample distribution:")
print(table(all_merged@meta.data$sample))
print(table(all_merged@meta.data$condition))
print(table(all_merged@meta.data$slide))


# -------------------------------------------------------
# 9. Save
# -------------------------------------------------------
saveRDS(mock,       file = file.path(OBJ_DIR, "00_mock_6wpi.rds"))
saveRDS(lcmv_1wpi,  file = file.path(OBJ_DIR, "00_LCMV_1wpi.rds"))
saveRDS(lcmv_3wpi,  file = file.path(OBJ_DIR, "00_LCMV_3wpi.rds"))
saveRDS(lcmv_6wpi,  file = file.path(OBJ_DIR, "00_LCMV_6wpi.rds"))
saveRDS(all_merged, file = file.path(OBJ_DIR, "00_all_merged.rds"))

message("\nDone. Objects saved to: ", OBJ_DIR)
message("  00_mock_6wpi.rds")
message("  00_LCMV_1wpi.rds")
message("  00_LCMV_3wpi.rds")
message("  00_LCMV_6wpi.rds")
message("  00_all_merged.rds = objet principal pour la suite")