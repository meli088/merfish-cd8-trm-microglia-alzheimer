# =============================================================
# Script: 04a_banksy_compute.R
# Project: LCMV MERFISH - TRM-Microglia niche analysis
# Author: Melina Farshchi
# Date: 2026-05
# Description: BANKSY neighborhood computation, PCA, and Harmony
#              correction for joint all-cells analysis.
#
#   1. Loads 02_all_normalized.rds (Seurat, SCTransform)
#   2. Converts each sample to SpatialExperiment
#   3. computeBanksy() per sample (preserves spatial neighborhoods)
#   4. Merges all samples with cbind()
#   5. For lambda 0.2 and 0.8:
#        runBanksyPCA() -> RunHarmony(sample + slide)
#   6. Saves 04_banksy_joint.rds
#
# Output reducedDimNames(se_joint):
#   "PCA_M1_lam0.2", "Harmony_lam0.2"
#   "PCA_M1_lam0.8", "Harmony_lam0.8"
#
# Next step: 04_banksy.R (04b_banksy_clustering.R)
# =============================================================

# -------------------------------------------------------
# 1. Libraries
# -------------------------------------------------------
library(Banksy)
library(SpatialExperiment)
library(SingleCellExperiment)
library(Seurat)
library(SeuratObject)
library(harmony)
library(scater)
library(ggplot2)
library(dplyr)

# -------------------------------------------------------
# 2. Parameters
# -------------------------------------------------------
OBJ_DIR    <- "objects"
BANKSY_DIR <- "outputs/banksy"
dir.create(BANKSY_DIR, showWarnings = FALSE, recursive = TRUE)

SEED         <- 1997
SAMPLE_NAMES <- c("mock_6wpi", "LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")
K_GEOM       <- 30
LAMBDAS      <- c(0.2, 0.8)
N_PCS        <- 30

# -------------------------------------------------------
# 3. Helpers
# -------------------------------------------------------
find_pca_dim <- function(se, lam) {
  dims    <- reducedDimNames(se)
  lam_chr <- as.character(lam)
  lam_nd  <- gsub("\\.", "", lam_chr)
  patterns <- c(
    paste0("^PCA_M1_lam", gsub("\\.", "\\\\.", lam_chr), "$"),
    paste0("^PCA_M1_lam", lam_nd, "$"),
    paste0("^PCA_M0_lam", gsub("\\.", "\\\\.", lam_chr), "$"),
    paste0("^PCA.*lam",   gsub("\\.", "\\\\.", lam_chr), "$"),
    paste0("^PCA.*lam",   lam_nd, "$")
  )
  for (pat in patterns) {
    hit <- dims[grepl(pat, dims)]
    if (length(hit) > 0) return(hit[1])
  }
  pca_hits <- dims[grepl("^PCA", dims)]
  if (length(pca_hits) > 0) return(tail(pca_hits, 1))
  NULL
}

run_gc <- function(step = "") {
  invisible(gc(verbose = FALSE))
  if (nzchar(step)) message("  GC: ", step)
}

# -------------------------------------------------------
# 4. Load normalized Seurat object
# -------------------------------------------------------
message("Loading objects/02_all_normalized.rds ...")
obj <- readRDS(file.path(OBJ_DIR, "02_all_normalized.rds"))
message("Seurat object: ", ncol(obj), " cells | ", nrow(obj), " genes")
print(table(obj@meta.data$sample))
SeuratObject::DefaultAssay(obj) <- "Vizgen"

# Join layers so counts are in a single matrix (required before NormalizeData).
# NOTE: JoinLayers() on a v5 assay that has only a counts layer will automatically
# populate the data layer with log1p(counts). We overwrite it explicitly below.
obj <- JoinLayers(obj)
message("Layers after JoinLayers: ", paste(Layers(obj), collapse = ", "))

# -- Normalisation : log1p(counts) ------------------------------------------
# JoinLayers() sur un assay Seurat v5 qui n'a qu'un layer "counts" peuple
# automatiquement le layer "data" avec log1p(counts) - sans normalisation
# par taille de librairie. C'est ce qui s'est passe dans les runs originaux
# (valeurs observees : log1p(1)=0.693, log1p(2)=1.099, etc.).
# On le fait ici explicitement pour etre transparent sur ce qui a ete utilise.
counts_mat <- GetAssayData(obj, assay = "Vizgen", layer = "counts")
LayerData(obj, assay = "Vizgen", layer = "data") <- log1p(counts_mat)
rm(counts_mat)
message("Normalisation : log1p(counts) stocke dans Vizgen 'data' layer.")

# -------------------------------------------------------
# 5. Build per-sample SpatialExperiment + computeBanksy
# -------------------------------------------------------
message("\n=== computeBanksy per sample (k_geom=", K_GEOM, ") ===")

se_list <- lapply(SAMPLE_NAMES, function(sname) {
  message("  ", sname, " ...")
  obj_s <- subset(obj, subset = sample == sname)

  # RC-normalized counts from Vizgen data layer (non-negative, per BANKSY requirements)
  gcm <- GetAssayData(obj_s, assay = "Vizgen", layer = "data")

  # Spatial coordinates
  fov_name <- Images(obj_s)[1]
  coords   <- GetTissueCoordinates(obj_s, image = fov_name)
  rownames(coords) <- coords$cell

  locs <- as.matrix(data.frame(
    sdimx = coords[rownames(obj_s@meta.data), "x"],
    sdimy = coords[rownames(obj_s@meta.data), "y"],
    row.names = rownames(obj_s@meta.data)
  ))

  valid <- !is.na(locs[, 1]) & !is.na(locs[, 2])
  if (sum(!valid) > 0) {
    message("    Removing ", sum(!valid), " cells with NA coordinates")
    gcm   <- gcm[, valid]
    locs  <- locs[valid, ]
    obj_s <- obj_s[, valid]
  }

  keep_cols <- intersect(c("sample", "condition", "slide"),
                         colnames(obj_s@meta.data))
  mdata <- obj_s@meta.data[, keep_cols, drop = FALSE]

  se_s <- SpatialExperiment(
    assay         = list(normcounts = as.matrix(gcm)),
    spatialCoords = locs,
    colData       = mdata,
    sample_id     = sname
  )

  Banksy::computeBanksy(
    se_s,
    compute_agf = TRUE,
    k_geom      = K_GEOM,
    assay_name  = "normcounts"
  )
})

run_gc("after computeBanksy")
message("cbind() all samples...")
se_joint <- do.call(cbind, se_list)
rm(se_list, obj); run_gc("after cbind")
message("Joint object: ", ncol(se_joint), " cells")
message("assays: ", paste(assayNames(se_joint), collapse = ", "))

# -------------------------------------------------------
# 6. runBanksyPCA + Harmony for each lambda
# -------------------------------------------------------
message("\n=== PCA + Harmony ===")

for (lam in LAMBDAS) {
  harmony_name <- paste0("Harmony_lam", lam)
  message("\n--- Lambda = ", lam, " ---")

  message("  runBanksyPCA ...")
  se_joint <- Banksy::runBanksyPCA(
    se_joint,
    lambda = lam,
    npcs   = N_PCS,
    seed   = SEED
  )

  pca_name <- find_pca_dim(se_joint, lam)
  if (is.null(pca_name)) {
    stop(
      "No PCA reducedDim found after runBanksyPCA for lambda=", lam,
      ". Available: ", paste(reducedDimNames(se_joint), collapse = ", ")
    )
  }
  message("  PCA: ", pca_name,
          " (", nrow(reducedDim(se_joint, pca_name)),
          " x ", ncol(reducedDim(se_joint, pca_name)), ")")

  message("  RunHarmony (sample + slide) ...")
  pca_mat <- reducedDim(se_joint, pca_name)

  batch_vars <- intersect(c("sample", "slide"),
                          names(colData(se_joint)))
  if (length(batch_vars) == 0) {
    stop("No batch variables (sample/slide) found in colData for Harmony.")
  }
  message("  Correcting for: ", paste(batch_vars, collapse = " + "))

  harmony_emb <- harmony::RunHarmony(
    data_mat  = pca_mat,
    meta_data = as.data.frame(colData(se_joint)),
    vars_use  = batch_vars,
    verbose   = FALSE
  )
  reducedDim(se_joint, harmony_name) <- harmony_emb
  message("  Stored: ", harmony_name)

  # UMAP on Harmony embedding
  umap_name <- paste0("UMAP_lam", lam)
  message("  runUMAP (", harmony_name, ") ...")
  se_joint <- scater::runUMAP(
    se_joint,
    dimred     = harmony_name,
    n_dimred   = N_PCS,
    name       = umap_name,
    BPPARAM    = BiocParallel::SerialParam(RNGseed = SEED)
  )
  message("  Stored: ", umap_name)

  # QC plot - UMAP coloured by sample
  umap_coords <- as.data.frame(reducedDim(se_joint, umap_name))
  colnames(umap_coords) <- c("UMAP1", "UMAP2")
  umap_coords$sample <- as.character(colData(se_joint)$sample)

  p_umap <- ggplot(umap_coords[sample(nrow(umap_coords)), ],
                   aes(UMAP1, UMAP2, color = sample)) +
    geom_point(size = 0.2, alpha = 0.4) +
    scale_color_brewer(palette = "Set1", name = "Sample") +
    labs(
      title    = paste0("UMAP - all cells, lambda = ", lam),
      subtitle = paste0("Harmony correction (sample + slide) | N_PCS = ", N_PCS)
    ) +
    theme_classic() +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

  plot_path_pdf <- file.path(BANKSY_DIR, paste0("umap_lam", gsub("\\.", "", as.character(lam)), "_by_sample.pdf"))
  plot_path_jpg <- file.path(BANKSY_DIR, paste0("umap_lam", gsub("\\.", "", as.character(lam)), "_by_sample.jpg"))
  ggsave(plot_path_pdf, p_umap, width = 7, height = 6, device = cairo_pdf)
  ggsave(plot_path_jpg, p_umap, width = 7, height = 6, dpi = 300)
  message("  UMAP plot saved: ", plot_path_jpg)

  run_gc(paste0("after lambda=", lam))
}

# -------------------------------------------------------
# 7. Save
# -------------------------------------------------------
message("\nReducedDims: ", paste(reducedDimNames(se_joint), collapse = ", "))
message("Cluster columns: ", length(clusterNames(se_joint)))

out_path <- file.path(OBJ_DIR, "04_banksy_joint.rds")
saveRDS(se_joint, out_path)
message("Saved: ", out_path)
message("\nDone. se_joint has ", ncol(se_joint), " cells and ",
        length(reducedDimNames(se_joint)), " reducedDims.")
message("Next step: 04_banksy.R (04b_banksy_clustering.R)")
