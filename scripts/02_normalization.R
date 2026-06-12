# =============================================================
# Script: 02_normalization.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-04
# Description: Normalization and dimensionality reduction.
#              Follows Feng et al. (Immunity 2025) pipeline:
#              - SCTransform per layer (one layer = one sample)
#              - PCA on all 496 genes
#              - Harmony batch correction via IntegrateLayers
#              - UMAP to verify biology vs batch separation
# Input:  objects/01_all_merged_qc.rds
# Output: objects/02_all_normalized.rds
#         outputs/normalization/
# =============================================================


# -------------------------------------------------------
# 1. Libraries
# -------------------------------------------------------
# remotes::install_github("satijalab/sctransform", ref = "main")
# packageVersion("sctransform")

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(patchwork)
library(harmony)

# Fail fast on known Seurat/SeuratObject API mismatch
if (packageVersion("Seurat") < "5.4.0" || packageVersion("SeuratObject") < "5.4.0") {
  stop(
    "Incompatible package versions for SCTransform on layered assays. ",
    "Please use Seurat >= 5.4.0 and SeuratObject >= 5.4.0. ",
    "Current: Seurat=", as.character(packageVersion("Seurat")),
    ", SeuratObject=", as.character(packageVersion("SeuratObject"))
  )
}

OBJ_DIR  <- "objects"
NORM_DIR <- "outputs/normalization"
dir.create(NORM_DIR, showWarnings = FALSE, recursive = TRUE)

options(future.globals.maxSize = 2000 * 1024^2)

FEATURE_COLS <- c("#FFFFCC", "#FED976", "#FD8D3C", "#E31A1C", "#800026")
CONT_COLS    <- c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#084594")


# -------------------------------------------------------
# 2. Load QC-filtered merged object
# -------------------------------------------------------
message("Loading QC-filtered merged object...")
obj <- readRDS(file.path(OBJ_DIR, "01_all_merged_qc.rds"))

message("Object: ", ncol(obj), " cells | ", nrow(obj), " genes")
print(table(obj@meta.data$sample))
print(table(obj@meta.data$slide))


# -------------------------------------------------------
# 3. Split into layers (one layer per sample)
# -------------------------------------------------------
message("\nSplitting into layers by sample...")
obj <- JoinLayers(obj)
obj[["Vizgen"]] <- split(obj[["Vizgen"]], f = obj$sample)

message("Layers created:")
print(Layers(obj))


# -------------------------------------------------------
# 4. Normalization — SCTransform per layer
# -------------------------------------------------------
message("\nRunning SCTransform per layer...")

obj <- SCTransform(
  obj,
  assay      = "Vizgen",
  clip.range = c(-10, 10),
  verbose    = FALSE
)

if (!("SCT" %in% Assays(obj))) {
  stop("SCTransform did not create the SCT assay. Check package versions and input assay layers.")
}

message("SCTransform done.")
message("SCT assay: ", nrow(obj[["SCT"]]), " variable features identified")


# -------------------------------------------------------
# 4.5. Inspect highly variable genes
# -------------------------------------------------------
hvg <- VariableFeatures(obj)
message("Top 20 HVGs: ", paste(head(hvg, 20), collapse = ", "))


# -------------------------------------------------------
# 5. PCA on ALL panel genes
# -------------------------------------------------------
message("\nRunning PCA on all panel genes...")

all_genes <- rownames(obj)

obj <- RunPCA(
  obj,
  features = all_genes,
  npcs     = 50,
  verbose  = FALSE
)

p_elbow <- ElbowPlot(obj, ndims = 50) +
  labs(title    = "Elbow plot — PCA",
       subtitle = "Set N_PCS in section 6 based on where the curve flattens") +
  geom_vline(xintercept = 30, color = "red", linetype = "dashed", alpha = 0.5) +
  annotate("text", x = 31, y = 2, label = "Feng: 30 PCs",
           color = "red", size = 3, hjust = 0, alpha = 0.7) +
  theme_classic()

ggsave(file.path(NORM_DIR, "elbow_plot.pdf"), p_elbow, width = 7, height = 5)
message("Elbow plot saved — inspect before continuing!")

# -------------------------------------------------------
# STOP HERE — inspect elbow_plot.pdf before running section 6
# -------------------------------------------------------


# -------------------------------------------------------
# 6. Harmony batch correction
# -------------------------------------------------------
# N_PCS chosen from elbow plot.
# Harmony corrects for both sample AND slide batch effects.
# vars_use = c("sample", "slide") handles:
#   - sample-level differences (Mock vs LCMV timepoints)
#   - slide-level differences (slide4 vs slide2 technical batch)

N_PCS <- 19   # <-- adjust after inspecting elbow plot

message("\nRunning Harmony batch correction (", N_PCS, " PCs)...")
message("Correcting for: sample + slide")

obj <- IntegrateLayers(
  object         = obj,
  method         = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction  = "harmony",
  group.by.vars  = c("sample", "slide"),
  verbose        = FALSE
)

message("Harmony done.")


# -------------------------------------------------------
# 7. Neighborhood graph and UMAP
# -------------------------------------------------------
message("\nBuilding neighbor graph and UMAP on Harmony embeddings...")

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:N_PCS, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:N_PCS, verbose = FALSE)

message("UMAP done.")


# -------------------------------------------------------
# 8. Visualizations
# -------------------------------------------------------
message("\nGenerating visualizations...")

p_sample <- DimPlot(
  obj, reduction = "umap", group.by = "sample",
  pt.size = 0.05, raster = FALSE
) + labs(
  title    = "UMAP — colored by sample",
  subtitle = paste0("Good: samples overlap | N_PCS = ", N_PCS)
)

p_condition <- DimPlot(
  obj, reduction = "umap", group.by = "condition",
  pt.size = 0.05, raster = FALSE
) + labs(title = "UMAP — LCMV vs Mock")

# Timepoint — nouvelle variable utile avec 4 échantillons
p_timepoint <- DimPlot(
  obj, reduction = "umap", group.by = "timepoint",
  pt.size = 0.05, raster = FALSE
) + labs(title = "UMAP — Timepoint (1/3/6wpi + Mock)")

# Slide — vérifier que les 2 slides se mélangent bien
p_slide <- DimPlot(
  obj, reduction = "umap", group.by = "slide",
  pt.size = 0.05, raster = FALSE
) + labs(title = "UMAP — Slide (slide2 vs slide4)")

p_counts <- FeaturePlot(
  obj, features = "nCount_Vizgen",
  reduction = "umap", pt.size = 0.05,
  min.cutoff = 0, max.cutoff = "q95", raster = FALSE
) + scale_color_gradientn(colors = CONT_COLS) +
  labs(title = "UMAP — nCount_Vizgen")

p_genes <- FeaturePlot(
  obj, features = "nFeature_Vizgen",
  reduction = "umap", pt.size = 0.05,
  min.cutoff = 0, max.cutoff = "q95", raster = FALSE
) + scale_color_gradientn(colors = CONT_COLS) +
  labs(title = "UMAP — nFeature_Vizgen")

p_iba1 <- FeaturePlot(
  obj, features = "Anti.IBA1_raw",
  reduction = "umap", pt.size = 0.05,
  min.cutoff = "q10", max.cutoff = "q95", raster = FALSE
) + scale_color_gradientn(colors = FEATURE_COLS) +
  labs(title = "UMAP — Anti-IBA1 (microglia)")

p_rfp <- FeaturePlot(
  obj, features = "Anti.RFP_raw",
  reduction = "umap", pt.size = 0.05,
  min.cutoff = "q10", max.cutoff = "q95", raster = FALSE
) + scale_color_gradientn(colors = FEATURE_COLS) +
  labs(title = "UMAP — Anti-RFP (TRM-contact)")

p_panel <- (p_sample    | p_condition) /
           (p_timepoint | p_slide) /
           (p_counts    | p_genes) /
           (p_iba1      | p_rfp)

ggsave(file.path(NORM_DIR, paste0("UMAP_N", N_PCS, "PCs.pdf")),
       p_panel, width = 14, height = 24)

message("UMAP plots saved.")


# -------------------------------------------------------
# 9. Save
# -------------------------------------------------------
saveRDS(obj, file.path(OBJ_DIR, "02_all_normalized.rds"))

message("\nDone. Object saved: objects/02_all_normalized.rds")
message("N_PCS used: ", N_PCS)
message("Figures in: ", NORM_DIR)
