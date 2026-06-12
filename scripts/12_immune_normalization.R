#!/usr/bin/env Rscript
# =============================================================
# Script: 12_immune_normalization.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Description: Normalize the joint immune subset.
#              Strictly mirrors 02_normalization.R:
#              - Vizgen assay (not RNA)
#              - JoinLayers + split by sample
#              - SCTransform per layer on Vizgen assay
#              - PCA on all panel genes
#              - IntegrateLayers(HarmonyIntegration)
#              - UMAP
#              Then rebuilds a SpatialExperiment for BANKSY (13).
#
# Input:  objects/05_immune_subset_lam02_res09_<label>.rds
# Output: objects/06_immune_normalized_lam02_<label>.rds      (Seurat)
#         objects/06_immune_normalized_lam02_<label>_spe.rds  (SpatialExperiment)
#         outputs/banksy/<label>/normalization/
# =============================================================

library(Seurat)
library(SeuratObject)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(SpatialExperiment)
library(ggplot2)
library(dplyr)
library(patchwork)
library(harmony)

base_path <- normalizePath(".")  # Run this script from the project root directory
setwd(base_path)

# -------------------------------------------------------
# Parameters
# -------------------------------------------------------
IMMUNE_LABEL <- "Immune (Acod1)"  # default; override with --label
FOLDER_NAME  <- NULL   # optional: override folder/file slug with --folder microglia

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  for (i in seq_along(args)) {
    a <- args[i]
    if (a %in% c("--label", "-l") && (i + 1) <= length(args)) {
      IMMUNE_LABEL <- args[i + 1]
    }
    if (a == "--folder" && (i + 1) <= length(args)) {
      FOLDER_NAME <- args[i + 1]
    }
  }
}

slugify_label <- function(label) {
  label <- tolower(trimws(label))
  label <- gsub("[^a-z0-9]+", "_", label)
  gsub("^_+|_+$", "", label)
}

safe_label <- if (!is.null(FOLDER_NAME)) tolower(trimws(FOLDER_NAME)) else slugify_label(IMMUNE_LABEL)
label_root <- file.path("outputs", "banksy", safe_label)
out_dir    <- file.path(label_root, "normalization")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

N_PCS <- 40   # adjust after inspecting elbow_plot.pdf

options(future.globals.maxSize = 2000 * 1024^2)

# -------------------------------------------------------
# 1. Load immune subset (SpatialExperiment)
# -------------------------------------------------------
in_candidates <- c(
  file.path("objects", paste0("05_immune_subset_lam02_res09_", safe_label, ".rds")),
  file.path("objects", "05_immune_subset_lam02_res09.rds")
)
in_file <- in_candidates[file.exists(in_candidates)][1]
if (is.na(in_file)) stop("Input immune subset not found for label: ", IMMUNE_LABEL)

spe <- readRDS(in_file)

ss <- as.data.frame(SummarizedExperiment::colData(spe))
uniq_samples <- sort(unique(as.character(ss$sample)))
message("Immune subset: ", ncol(spe), " cells across ", length(uniq_samples), " samples")
message("Samples: ", paste(uniq_samples, collapse = ", "))

coord_df <- as.data.frame(SpatialExperiment::spatialCoords(spe))

# -------------------------------------------------------
# 2. Build Seurat object with Vizgen assay (mirrors 02)
# -------------------------------------------------------
if (!"counts" %in% assayNames(spe)) {
  stop("No raw counts assay found in immune subset — refusing to proceed")
}

so <- CreateSeuratObject(
  counts    = assay(spe, "counts"),
  assay     = "Vizgen",
  meta.data = as.data.frame(SummarizedExperiment::colData(spe))
)
DefaultAssay(so) <- "Vizgen"

message("Seurat object: ", ncol(so), " cells | ", nrow(so), " genes")
print(table(so@meta.data$sample))

# -------------------------------------------------------
# 3. SCTransform on Vizgen assay (whole-object)
#    The immune subset may contain very small samples
#    (e.g. mock: 58 cells) which break per-layer SCTransform.
#    We use whole-object SCTransform on the Vizgen assay —
#    same assay name as 02_normalization.R.
# -------------------------------------------------------
so <- JoinLayers(so)
message("\nRunning SCTransform (Vizgen assay, whole-object)...")
so <- SCTransform(
  so,
  assay      = "Vizgen",
  clip.range = c(-10, 10),
  verbose    = FALSE
)

if (!("SCT" %in% SeuratObject::Assays(so))) {
  stop("SCTransform did not create the SCT assay.")
}
message("SCTransform done. Variable features: ", length(VariableFeatures(so)))

# -------------------------------------------------------
# 4. PCA on all panel genes (mirrors 02)
# -------------------------------------------------------
message("\nRunning PCA on all panel genes...")
all_genes <- rownames(so)
so <- RunPCA(so, features = all_genes, npcs = 50, verbose = FALSE)

p_elbow <- ElbowPlot(so, ndims = 50) +
  labs(title    = "Elbow plot — PCA (immune subset)",
       subtitle = paste0("Set N_PCS in script 12 — currently ", N_PCS)) +
  geom_vline(xintercept = N_PCS, color = "red", linetype = "dashed", alpha = 0.6) +
  annotate("text", x = N_PCS + 1, y = 2, label = paste0("N_PCS = ", N_PCS),
           color = "red", size = 3, hjust = 0, alpha = 0.8) +
  theme_classic()

ggsave(file.path(out_dir, "elbow_plot.pdf"), p_elbow, width = 7, height = 5)
message("Elbow plot saved — inspect before adjusting N_PCS!")

# -------------------------------------------------------
# 5. Harmony batch correction (mirrors 02 variables)
#    RunHarmony is used directly because IntegrateLayers
#    requires split layers which are incompatible with
#    whole-object SCTransform on small-sample subsets.
# -------------------------------------------------------
batch_vars <- intersect(c("sample", "slide"), colnames(so@meta.data))
if (length(batch_vars) == 0) stop("No batch variables found for Harmony.")

message("\nRunning Harmony (", N_PCS, " PCs, correcting for: ",
        paste(batch_vars, collapse = " + "), ")...")

so <- RunHarmony(
  object         = so,
  group.by.vars  = batch_vars,
  reduction.use  = "pca",
  verbose        = FALSE
)
message("Harmony done.")

# -------------------------------------------------------
# 6. Neighbor graph + UMAP (mirrors 02)
# -------------------------------------------------------
message("\nBuilding neighbor graph and UMAP...")
so <- FindNeighbors(so, reduction = "harmony", dims = 1:N_PCS, verbose = FALSE)
so <- RunUMAP(so,     reduction = "harmony", dims = 1:N_PCS, verbose = FALSE)
message("UMAP done.")

# -------------------------------------------------------S
# 7. Save Seurat object
# -------------------------------------------------------
obj_out <- file.path("objects", paste0("06_immune_normalized_lam02_", safe_label, ".rds"))
saveRDS(so, obj_out)
message("Saved Seurat object: ", obj_out)

# -------------------------------------------------------
# 8. Convert to SpatialExperiment for downstream BANKSY (13)
# -------------------------------------------------------
message("\nConverting to SpatialExperiment for BANKSY...")
# For conversion we use the Vizgen (raw counts) assay only.
# Temporarily drop SCT to avoid JoinLayers issues with SCTAssay.
so_for_spe <- so
DefaultAssay(so_for_spe) <- "Vizgen"
so_for_spe[["SCT"]] <- NULL

so_for_spe <- JoinLayers(so_for_spe)

tryCatch(
  so_for_spe[["Vizgen"]] <- as(so_for_spe[["Vizgen"]], "Assay"),
  error = function(e) message("Note: Vizgen assay v5 cast skipped")
)

sce <- as.SingleCellExperiment(so_for_spe)

if ("pca"     %in% reducedDimNames(sce)) reducedDim(sce, "PCA")     <- reducedDim(sce, "pca")
if ("harmony" %in% reducedDimNames(sce)) reducedDim(sce, "HARMONY") <- reducedDim(sce, "harmony")
if ("umap"    %in% reducedDimNames(sce)) reducedDim(sce, "UMAP")    <- reducedDim(sce, "umap")

spe_out <- SpatialExperiment(
  assays        = list(counts = assay(sce, "counts")),
  rowData       = rowData(sce),
  colData       = colData(sce),
  spatialCoords = as.matrix(coord_df[, c("sdimx", "sdimy")])
)
reducedDims(spe_out) <- reducedDims(sce)

spe_file <- file.path("objects", paste0("06_immune_normalized_lam02_", safe_label, "_spe.rds"))
saveRDS(spe_out, spe_file)
message("Saved SpatialExperiment: ", spe_file)

# -------------------------------------------------------
# 10. QC plots
# -------------------------------------------------------
p_sample <- DimPlot(
  so, reduction = "umap", group.by = "sample", pt.size = 0.1, raster = FALSE
) + labs(title    = paste0("UMAP — ", IMMUNE_LABEL, " by sample"),
         subtitle = paste0("Harmony integration | N_PCS = ", N_PCS)) +
  theme_minimal()

pdf(file.path(out_dir, "immune_normalization_QC.pdf"), width = 10, height = 8)
print(p_sample)
print(p_elbow)
dev.off()

ggsave(file.path(out_dir, "umap_by_sample.jpg"), p_sample, width = 7, height = 6, dpi = 300)
ggsave(file.path(out_dir, "pca_elbow.jpg"),      p_elbow,  width = 6, height = 4, dpi = 300)

message("QC plots saved: ", out_dir)
message("\nImmune normalization complete.")
message("N_PCS used: ", N_PCS, " — adjust at top of script if needed after inspecting elbow plot.")
