# =============================================================
# Script: 02_normalization.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-04
# Description: Normalization and dimensionality reduction of the
#              QC-filtered merged object (114,588 cells).
#              Follows Feng et al. (Immunity 2025) pipeline:
#              - SCTransform per layer (one layer = one sample)
#              - PCA on all 496 genes
#              - Harmony batch correction via IntegrateLayers
#              - UMAP to verify biology vs batch separation
# Input:  objects/01_slide4_merged_qc.rds
# Output: objects/02_slide4_normalized.rds
#         outputs/normalization/  (figures)
#
# WORKFLOW:
#   Run sections 1-5 first (up to elbow plot)
#   → inspect elbow plot to choose N_PCS
#   → set N_PCS in section 6
#   → run sections 6-9
# =============================================================


# -------------------------------------------------------
# 1. Libraries
# -------------------------------------------------------

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(patchwork)
library(harmony)

OBJ_DIR  <- "objects"
NORM_DIR <- "outputs/normalization"
dir.create(NORM_DIR, showWarnings = FALSE, recursive = TRUE)

options(future.globals.maxSize = 2000 * 1024^2)

# Color palettes
# FEATURE_COLS: yellow (low) -> orange -> dark red (high)
# CONT_COLS:    white (low) -> dark blue (high) for metadata
FEATURE_COLS <- c("#FFFFCC", "#FED976", "#FD8D3C", "#E31A1C", "#800026")
CONT_COLS    <- c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#084594")


# -------------------------------------------------------
# 2. Load QC-filtered merged object
# -------------------------------------------------------

message("Loading QC-filtered merged object...")
obj <- readRDS(file.path(OBJ_DIR, "01_slide4_merged_qc.rds"))

message("Object: ", ncol(obj), " cells | ", nrow(obj), " genes")
print(table(obj@meta.data$sample))


# -------------------------------------------------------
# 3. Split into layers (one layer per sample)
# -------------------------------------------------------
# Feng et al.: "the data from different samples were segmented
# into separate layers. Subsequently, SCTransform was applied
# to each layer with clip.range = c(-10, 10)"

message("\nSplitting into layers by sample...")
obj <- JoinLayers(obj)
obj[["Vizgen"]] <- split(obj[["Vizgen"]], f = obj$sample)

message("Layers created:")
print(Layers(obj))


# -------------------------------------------------------
# 4. Normalization — SCTransform per layer
# -------------------------------------------------------
# SCTransform in one step:
#   (a) normalizes counts (corrects for sequencing depth)
#   (b) identifies highly variable genes (HVGs)
#   (c) scales the data (Pearson residuals)
# clip.range = c(-10, 10) limits extreme residuals

message("\nRunning SCTransform per layer...")

obj <- SCTransform(
  obj,
  assay      = "Vizgen",
  clip.range = c(-10, 10),
  verbose    = FALSE
)

message("SCTransform done.")
message("SCT assay: ", nrow(obj[["SCT"]]), " variable features identified")


# -------------------------------------------------------
# 4.5. Inspect highly variable genes (HVGs)
# -------------------------------------------------------
# SCTransform identified these HVGs internally.
# We visualize them to understand which genes drive variation,
# even though we will use ALL 496 genes for PCA (Feng et al.).


hvg <- VariableFeatures(obj)
message("Top 20 HVGs: ", paste(head(hvg, 20), collapse = ", "))


# -------------------------------------------------------
# 5. PCA on ALL panel genes
# -------------------------------------------------------
# Feng et al.: "PCA was performed using all 500 genes"
# We use all 496 genes — no HVG selection needed since the
# MERFISH panel was pre-selected for biological relevance.
# 50 PCs computed; N_PCS chosen after inspecting elbow plot.

message("\nRunning PCA on all panel genes...")

all_genes <- rownames(obj)

obj <- RunPCA(
  obj,
  features = all_genes,
  npcs     = 50,
  verbose  = FALSE
)

# Elbow plot — INSPECT THIS before setting N_PCS below
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
# Choose N_PCS = where the curve clearly flattens
# -------------------------------------------------------


# -------------------------------------------------------
# 6. Harmony batch correction
# -------------------------------------------------------
# ADJUST N_PCS based on elbow plot before running this section
# Current choice: 19 (curve flattens around PC 10-19 on our data)
# Feng et al. used 30 (different dataset, more heterogeneous)

N_PCS <- 19   # <-- adjust here after inspecting elbow plot

message("\nRunning Harmony batch correction (", N_PCS, " PCs)...")

obj <- IntegrateLayers(
  object         = obj,
  method         = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction  = "harmony",
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
# 8. Visualizations — assess batch correction
# -------------------------------------------------------
# Key checks:
# - sample UMAP: Mock and LCMV should overlap (not separate)
# - nCount UMAP: should be uniform (no gradient = good normalization)
# - IBA1/RFP: protein signals for reference

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

p_counts <- FeaturePlot(
  obj, features = "nCount_Vizgen",
  reduction = "umap", pt.size = 0.05,
  min.cutoff = 0, max.cutoff = "q95", raster = FALSE
) + scale_color_gradientn(colors = CONT_COLS) +
  labs(title = "UMAP — nCount_Vizgen (should be uniform)")

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
  labs(title = "UMAP — Anti-IBA1 (microglia protein signal)")

p_rfp <- FeaturePlot(
  obj, features = "Anti.RFP_raw",
  reduction = "umap", pt.size = 0.05,
  min.cutoff = "q10", max.cutoff = "q95", raster = FALSE
) + scale_color_gradientn(colors = FEATURE_COLS) +
  labs(title = "UMAP — Anti-RFP (TRM-contact protein signal)")

p_panel <- (p_sample | p_condition) /
           (p_counts  | p_genes) /
           (p_iba1    | p_rfp)

ggsave(file.path(NORM_DIR, paste0("UMAP_N", N_PCS, "PCs.pdf")),
       p_panel, width = 14, height = 18)

message("UMAP plots saved (filename includes N_PCS for comparison).")


# -------------------------------------------------------
# 9. Save
# -------------------------------------------------------

saveRDS(obj, file.path(OBJ_DIR, "02_slide4_normalized.rds"))

message("\nDone. Object saved: objects/02_slide4_normalized.rds")
message("N_PCS used: ", N_PCS)
message("Figures in: ", NORM_DIR)
