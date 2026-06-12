# =============================================================
# Script: 03_clustering_seurat.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-04
# =============================================================

package.version("Seurat")
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(Banksy)
library(SpatialExperiment)
library(scuttle)
library(harmony)

OBJ_DIR     <- "objects"
CLUSTER_DIR <- "outputs/clustering"
dir.create(CLUSTER_DIR, showWarnings = FALSE, recursive = TRUE)

options(future.globals.maxSize = 2000 * 1024^2)

FEATURE_COLS <- c("#FFFFCC", "#FED976", "#FD8D3C", "#E31A1C", "#800026")
CONT_COLS    <- c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#084594")


# -------------------------------------------------------
# 2. Load normalized object
# -------------------------------------------------------
message("Loading normalized object...")
obj <- readRDS(file.path(OBJ_DIR, "02_all_normalized.rds"))

message("Object: ", ncol(obj), " cells")
print(table(obj@meta.data$sample))
print(table(obj@meta.data$slide))

SeuratObject::DefaultAssay(obj) <- "Vizgen"


# =============================================================
# PART A — SEURAT CLUSTERING
# =============================================================

N_PCS <- 19

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:N_PCS, verbose = FALSE)

obj <- FindClusters(obj, graph.name = "SCT_snn",
                    resolution = 0.3, algorithm = 1, verbose = FALSE)
obj$seurat_clusters_res0.3 <- obj$seurat_clusters

obj <- FindClusters(obj, graph.name = "SCT_snn",
                    resolution = 0.5, algorithm = 1, verbose = FALSE)
obj$seurat_clusters_res0.5 <- obj$seurat_clusters

# ── RESOLUTION ACTIVE : 0.5 ──────────────────────────────────
Idents(obj) <- "seurat_clusters_res0.5"

message("Clusters found at res 0.3: ", length(unique(obj$seurat_clusters_res0.3)))
message("Clusters found at res 0.5: ", length(unique(obj$seurat_clusters_res0.5)))

# ── Clustree — stabilité des clusters selon la résolution ──
if (!requireNamespace("clustree", quietly = TRUE)) install.packages("clustree")
library(clustree)

# Calculer plusieurs résolutions pour le clustree
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)) {
  obj <- FindClusters(obj, graph.name = "SCT_snn",
                      resolution = res, algorithm = 1, verbose = FALSE)
}

p_clustree <- clustree(obj, prefix = "SCT_snn_res.")

ggsave(file.path(CLUSTER_DIR, "clustree_resolutions.pdf"),
       p_clustree, width = 14, height = 18)
message("Clustree saved.")

# -------------------------------------------------------
# 4. UMAP
# -------------------------------------------------------
p_cl03 <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters_res0.3",
                   label = TRUE, label.size = 4, pt.size = 0.05, raster = FALSE) + NoLegend() +
  labs(title = "Seurat clusters — resolution 0.3")

p_cl05 <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters_res0.5",
                   label = TRUE, label.size = 3, pt.size = 0.05, raster = FALSE) + NoLegend() +
  labs(title = "Seurat clusters — resolution 0.5")

p_sample    <- DimPlot(obj, reduction = "umap", group.by = "sample",
                        pt.size = 0.05, raster = FALSE) + labs(title = "Sample")
p_timepoint <- DimPlot(obj, reduction = "umap", group.by = "timepoint",
                        pt.size = 0.05, raster = FALSE) + labs(title = "Timepoint")

ggsave(file.path(CLUSTER_DIR, "UMAP_seurat_clusters.pdf"),
       (p_cl03 | p_cl05) / (p_sample | p_timepoint), width = 14, height = 12)
message("Cluster UMAP saved.")


# -------------------------------------------------------
# 5. Module scores

# Joindre les layers Vizgen avant les module scores
SeuratObject::DefaultAssay(obj) <- "Vizgen"
obj <- NormalizeData(obj, verbose = FALSE)
obj <- JoinLayers(obj, assay = "Vizgen")

# Vérifier
Layers(obj)
# doit afficher "counts" et "data" — plus de "counts.mock_6wpi" etc.
# -------------------------------------------------------
panel_genes <- rownames(obj)

homeo_genes    <- c("Cx3cr1", "P2ry12", "P2ry13", "Tmem119", "Hexb", "Csf1r", "Sall1", "Siglech")
dam_genes      <- c("Trem2", "Apoe", "Lpl", "Cst7", "Cd9", "Tyrobp", "Ctsd", "Itgax",
                    "Clec7a", "Axl", "Gpnmb", "Lgals3", "Spp1")
trm_cyto_genes <- c("Itgae", "Gzmb", "Gzmk", "Prf1", "Nkg7")
exhaust_genes  <- c("Pdcd1", "Lag3", "Havcr2", "Tox", "Tigit")
# ── Astrocyte pan + réactifs fusionnés ───────────────────────
astro_genes    <- c("Aldh1l1", "Aqp4", "S100b", "Gfap", "Vim", "Lcn2", "Serpina3n")
oligo_genes    <- c("Mog", "Mbp", "Plp1", "Sox10", "Olig2", "Mal", "Ermn", "Opalin")
neuron_genes   <- c("Snap25", "Rbfox3", "Syt1", "Syt2", "Map2", "Nefh", "Nefm")

filter_genes <- function(genes, panel) {
  found <- genes[genes %in% panel]
  message("  ", length(found), "/", length(genes), " : ", paste(found, collapse = ", "))
  found
}

message("Filtering gene sets by the 496-gene panel :")
message("Homeostatic:")  ; homeo_genes     <- filter_genes(homeo_genes, panel_genes)
message("DAM:")          ; dam_genes       <- filter_genes(dam_genes, panel_genes)
message("TRM_Cytotoxic:"); trm_cyto_genes  <- filter_genes(trm_cyto_genes, panel_genes)
message("Exhaustion:")   ; exhaust_genes   <- filter_genes(exhaust_genes, panel_genes)
message("Astrocyte:")    ; astro_genes     <- filter_genes(astro_genes, panel_genes)
message("Oligo:")        ; oligo_genes     <- filter_genes(oligo_genes, panel_genes)
message("Neurones:")     ; neuron_genes    <- filter_genes(neuron_genes, panel_genes)

score_list  <- list(homeo_genes, dam_genes, trm_cyto_genes,
                    exhaust_genes, astro_genes, oligo_genes, neuron_genes)
score_names <- c("Homeostatic", "DAM", "TRM_Cytotoxic",
                 "Exhaustion", "Astrocyte", "Oligodendrocyte", "Neuron")

SeuratObject::DefaultAssay(obj) <- "Vizgen"
for (i in seq_along(score_list)) {
  if (length(score_list[[i]]) > 0) {
    obj <- AddModuleScore(obj, features = list(score_list[[i]]),
                          name = paste0(score_names[i], "_score"), ctrl = 20)
    message("Score calculated: ", score_names[i])
  }
}

score_cols <- paste0(score_names, "_score1")
score_cols <- score_cols[score_cols %in% colnames(obj@meta.data)]

score_plots <- lapply(score_cols, function(s) {
  FeaturePlot(obj, features = s, reduction = "umap",
              pt.size = 0.03, min.cutoff = "q10", max.cutoff = "q95", raster = FALSE) +
    scale_color_gradientn(colors = FEATURE_COLS) +
    labs(title = gsub("_score1", "", s))
})

p_scores <- wrap_plots(score_plots, ncol = 3)
ggsave(file.path(CLUSTER_DIR, "UMAP_module_scores.pdf"),
       p_scores, width = 18, height = ceiling(length(score_cols)/3) * 5)
message("Module score UMAPs saved.")

p_vln_dam <- VlnPlot(obj, features = "DAM_score1",
                      group.by = "seurat_clusters_res0.5", pt.size = 0) +
  labs(title = "DAM score by cluster") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")

p_vln_homeo <- VlnPlot(obj, features = "Homeostatic_score1",
                         group.by = "seurat_clusters_res0.5", pt.size = 0) +
  labs(title = "Homeostatic score by cluster")

ggsave(file.path(CLUSTER_DIR, "violin_DAM_homeostatic_by_cluster.pdf"),
       p_vln_dam / p_vln_homeo, width = 14, height = 10)
message("Violin plots saved.")


# -------------------------------------------------------
# 6. Cluster composition
# -------------------------------------------------------
comp_df <- obj@meta.data %>%
  group_by(seurat_clusters_res0.5, timepoint) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters_res0.5) %>%
  mutate(pct = n / sum(n) * 100)

p_comp <- ggplot(comp_df, aes(x = seurat_clusters_res0.5, y = pct, fill = timepoint)) +
  geom_col() +
  scale_fill_manual(values = c("6wpi" = "#95A5A6", "1wpi" = "#E74C3C",
                                "3wpi" = "#E67E22", "6wpi_mock" = "#BDC3C7")) +
  labs(title = "Cluster composition — par timepoint",
       subtitle = "Verify Harmony did not erase biological signal",
       x = "Cluster", y = "% of cells") +
  theme_classic()

ggsave(file.path(CLUSTER_DIR, "cluster_composition_by_timepoint.pdf"),
       p_comp, width = 14, height = 5)
message("Cluster composition plot saved.")


# -------------------------------------------------------
# 7. FindAllMarkers
# -------------------------------------------------------
message("\nFinding cluster markers...")

DefaultAssay(obj) <- "Vizgen"
obj <- NormalizeData(obj, verbose = FALSE)
obj <- SeuratObject::JoinLayers(obj, assay = "Vizgen")

Idents(obj) <- "seurat_clusters_res0.5"

all_markers <- FindAllMarkers(obj, assay = "Vizgen", only.pos = TRUE,
                               min.pct = 0.1, logfc.threshold = 0.25,
                               test.use = "wilcox", verbose = FALSE)

message("Markers found: ", nrow(all_markers))
write.csv(all_markers, file = file.path(CLUSTER_DIR, "cluster_markers_res0.5.csv"),
          row.names = FALSE)

if (nrow(all_markers) > 0) {
  top5 <- all_markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5) %>% ungroup()
  message("Top 5 markers per cluster:")
  print(top5 %>% select(cluster, gene, avg_log2FC, p_val_adj))
} else {
  top5 <- data.frame()
}


# -------------------------------------------------------
# 8. Heatmaps style Feng
# -------------------------------------------------------
heatmap_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

make_pheatmap <- function(seurat_obj, group_col, assay = "Vizgen",
                           marker_genes, title, filename) {
  avg <- AverageExpression(seurat_obj, assays = assay, layer = "data",
                            group.by = group_col)[[assay]]
  genes_present <- marker_genes[marker_genes %in% rownames(avg)]
  if (length(genes_present) == 0) { message("No genes for: ", title); return(invisible(NULL)) }
  avg_sub    <- avg[genes_present, , drop = FALSE]
  avg_scaled <- t(scale(t(avg_sub)))
  avg_scaled <- pmin(pmax(avg_scaled, -2), 2)
  pheatmap(avg_scaled, color = heatmap_colors, cluster_rows = TRUE, cluster_cols = TRUE,
           show_rownames = TRUE, show_colnames = TRUE, fontsize_row = 9, fontsize_col = 10,
           border_color = NA, main = title, filename = filename, width = 8, height = 6)
  message("Saved: ", filename)
}

SeuratObject::DefaultAssay(obj) <- "Vizgen"
obj <- NormalizeData(obj, verbose = FALSE)

top2_seurat <- all_markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 2) %>%
  ungroup() %>% pull(gene) %>% unique()

obj@meta.data$cluster_label <- paste0("Cluster_", obj@meta.data$seurat_clusters_res0.5)

make_pheatmap(obj, "cluster_label", marker_genes = top2_seurat,
              title = "DEGs par cluster — tous échantillons",
              filename = file.path(CLUSTER_DIR, "heatmap_seurat_merged.pdf"))

for (sname in c("mock_6wpi", "LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")) {
  obj_sub <- subset(obj, subset = sample == sname)
  obj_sub@meta.data$cluster_label <- paste0("Cluster_", obj_sub@meta.data$seurat_clusters_res0.5)
  make_pheatmap(obj_sub, "cluster_label", marker_genes = top2_seurat,
                title = paste("DEGs par cluster —", sname),
                filename = file.path(CLUSTER_DIR, paste0("heatmap_seurat_", gsub("_", "", sname), ".pdf")))
}
message("Heatmaps Seurat sauvegardées.")


# -------------------------------------------------------
# 9. Dotplot
# -------------------------------------------------------
canonical_markers <- c(
  "Cx3cr1", "P2ry12", "Tmem119", "Hexb",
  "Trem2", "Apoe", "Lpl", "Cst7",
  "Cd8a", "Cd3e", "Gzmb",
  "Gfap", "Aqp4", "S100b",
  "Mog", "Plp1", "Olig2",
  "Snap25", "Rbfox3",
  "Pecam1", "Cldn5"
)
canonical_markers <- canonical_markers[canonical_markers %in% panel_genes]

p_dot <- DotPlot(obj, features = canonical_markers,
                  group.by = "seurat_clusters_res0.5") + RotatedAxis() +
  scale_color_gradientn(colors = FEATURE_COLS) +
  labs(title = "Canonical markers by cluster — use to annotate",
       subtitle = "Dot size = % expressing | Color = mean expression")

ggsave(file.path(CLUSTER_DIR, "dotplot_canonical_markers.pdf"),
       p_dot, width = 16, height = 7)
message("Dotplot saved.")


# -------------------------------------------------------
# 9.5. Sub-clustering — cluster immunitaire
# -------------------------------------------------------
# !! Vérifier le numéro du cluster immunitaire dans le dotplot !!
# C'est le cluster qui exprime P2ry12 + Gzmb ensemble

IMMUNE_CLUSTER <- "7"  # cluster immunitaire/myéloïde — P2ry12, Tmem119, Hexb, Csf1r

message("\n--- SUB-CLUSTERING CLUSTER ", IMMUNE_CLUSTER, " ---")

obj_c2 <- subset(obj, subset = seurat_clusters_res0.5 == IMMUNE_CLUSTER)
message("Cluster ", IMMUNE_CLUSTER, ": ", ncol(obj_c2), " cells")

SeuratObject::DefaultAssay(obj_c2) <- "Vizgen"
obj_c2 <- FindNeighbors(obj_c2, reduction = "harmony", dims = 1:N_PCS, verbose = FALSE)
obj_c2 <- FindClusters(obj_c2, graph.name = "SCT_snn",
                        resolution = 0.3, verbose = FALSE)
obj_c2 <- RunUMAP(obj_c2, reduction = "harmony", dims = 1:N_PCS, verbose = FALSE)

message("Sub-clusters found: ", length(unique(Idents(obj_c2))))

p1 <- DimPlot(obj_c2, label = TRUE, pt.size = 0.5) +
  labs(title = paste0("Sous-clusters du cluster ", IMMUNE_CLUSTER))
p2 <- FeaturePlot(obj_c2, features = c("Hexb", "Gzmb", "Cd8a", "Havcr2"),
                   ncol = 2, min.cutoff = 0, max.cutoff = "q95") &
  scale_color_gradientn(colors = FEATURE_COLS)

ggsave(file.path(CLUSTER_DIR, "subcluster2_UMAP.pdf"), p1, width = 7, height = 6)
ggsave(file.path(CLUSTER_DIR, "subcluster2_markers.pdf"), p2, width = 10, height = 8)

SeuratObject::DefaultAssay(obj_c2) <- "Vizgen"
obj_c2 <- NormalizeData(obj_c2, verbose = FALSE)
obj_c2 <- JoinLayers(obj_c2, assay = "Vizgen")

markers_c2 <- FindAllMarkers(obj_c2, assay = "Vizgen", only.pos = TRUE,
                               min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
write.csv(markers_c2, file.path(CLUSTER_DIR, "subcluster2_markers.csv"), row.names = FALSE)
message("Top 3 markers per sub-cluster:")
print(markers_c2 %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 3) %>%
        select(cluster, gene, avg_log2FC, p_val_adj) %>% print(n = 40))


# -------------------------------------------------------
# 10. Manual annotation
# -------------------------------------------------------
# !! À mettre à jour après inspection des marqueurs du sous-clustering !!
# Les labels ci-dessous sont basés sur les résultats précédents (run slide4+slide2)
# Vérifier avec le dotplot et les marqueurs de la résolution 0.5

message("\n--- MANUAL ANNOTATION ---")

# Labels pour les sous-clusters du cluster immunitaire
# Basés sur les marqueurs observés au run précédent :

subcluster_labels <- c(
  "0" = "Contamination_neuron",
  "1" = "Microglia_homeostatic",
  "2" = "Contamination_oligo",
  "3" = "Contamination_neuron",
  "4" = "TRM_CD8",
  "5" = "Contamination_vascular",
  "6" = "Contamination_astrocyte",
  "7" = "Contamination_neuron",
  "8" = "Contamination_OPC"
)

obj_c2@meta.data$cell_type_c2 <-
  subcluster_labels[as.character(obj_c2@meta.data$seurat_clusters)]

# Annotation clusters principaux (résolution 0.5)
# !! À ajuster après inspection du dotplot !!
cell_type_map <- c(
  "0"  = "Neuron",
  "1"  = "Astrocyte",
  "2"  = "Vascular",          # endothelial + pericyte
  "3"  = "Oligodendrocyte",
  "4"  = "Neuron",
  "5"  = "Neuron",
  "6"  = "Neuron",
  "7"  = "Immune_myeloid",    # à sous-clustériser
  "8"  = "Neuron",
  "9"  = "Neuron",
  "10" = "Neuron",
  "11" = "Microglia_reactive",
  "12" = "Neuron",
  "13" = "OPC",
  "14" = "Astrocyte",
  "15" = "Neuron",
  "16" = "Neuron",
  "17" = "Oligodendrocyte"    # NFOL / immature
)

obj@meta.data$cell_type <-
  cell_type_map[as.character(obj@meta.data$seurat_clusters_res0.5)]

cells_in_obj  <- rownames(obj@meta.data)[obj@meta.data$seurat_clusters_res0.5 == IMMUNE_CLUSTER]
cells_in_obj2 <- rownames(obj_c2@meta.data)
matched       <- intersect(cells_in_obj, cells_in_obj2)

message("Cluster ", IMMUNE_CLUSTER, " cells: ", length(cells_in_obj),
        " | Matched: ", length(matched))
obj@meta.data[matched, "cell_type"] <- obj_c2@meta.data[matched, "cell_type_c2"]

unmatched <- setdiff(cells_in_obj, matched)
if (length(unmatched) > 0) {
  message("Warning: ", length(unmatched), " unmatched — labelled 'Unknown'")
  obj@meta.data[unmatched, "cell_type"] <- "Unknown"
}

message("\nCell type distribution:")
print(table(obj@meta.data$cell_type, useNA = "ifany"))
Idents(obj) <- "cell_type"

ct_colors <- c(
  "Neuron"                    = "#E41A1C",
  "Astrocyte"                 = "#4DAF4A",
  "Oligodendrocyte"           = "#377EB8",
  "OPC"                       = "#66C2A5",
  "Vascular"                  = "#FC8D62",
  "Immune_myeloid"            = "#A6761D",
  "Microglia_homeostatic"     = "#1B9E77",
  "Microglia_reactive"        = "#D95F02",
  "TRM_CD8"                   = "#7570B3",
  "Contamination_neuron"      = "#CCCCCC",
  "Contamination_oligo"       = "#DDDDDD",
  "Contamination_OPC"         = "#EEEEEE",
  "Contamination_vascular"    = "#BBBBBB",
  "Contamination_astrocyte"   = "#AAAAAA",
  "Unknown"                   = "#000000"
)

present_types <- names(ct_colors)[names(ct_colors) %in% unique(obj@meta.data$cell_type)]

p_annot <- DimPlot(obj, reduction = "umap", group.by = "cell_type",
                    cols = ct_colors[present_types], label = TRUE, label.size = 3,
                    pt.size = 0.05, raster = FALSE) + NoLegend() +
  labs(title = "Annotation finale — types cellulaires",
       subtitle = paste0("Slide2 + Slide4 | ", ncol(obj), " cellules"))

p_annot_legend <- DimPlot(obj, reduction = "umap", group.by = "cell_type",
                            cols = ct_colors[present_types], pt.size = 0.03, raster = FALSE) +
  labs(title = "Annotation finale — avec légende")

ggsave(file.path(CLUSTER_DIR, "UMAP_annotated_final.pdf"), p_annot, width = 10, height = 8)
ggsave(file.path(CLUSTER_DIR, "UMAP_annotated_final_legend.pdf"), p_annot_legend, width = 14, height = 8)
message("Annotation UMAPs saved.")

saveRDS(obj,    file.path(OBJ_DIR, "03_all_clustered.rds"))
saveRDS(obj_c2, file.path(OBJ_DIR, "03_cluster2_subclustered.rds"))
message("Objects saved.")


# =============================================================
# PART B — BANKSY SPATIAL CLUSTERING
# =============================================================
message("\n\n--- BANKSY SPATIAL CLUSTERING ---")

run_banksy_for_sample <- function(seurat_obj, sample_name, obj_dir, cluster_dir) {
  message("\nRunning BANKSY for: ", sample_name)
  obj_s <- subset(seurat_obj, subset = sample == sample_name)
  gcm   <- GetAssayData(obj_s, assay = "Vizgen", layer = "data")

  fov_name       <- Images(obj_s)[1]
  coords         <- GetTissueCoordinates(obj_s, image = fov_name)
  coords_indexed <- coords
  rownames(coords_indexed) <- coords_indexed$cell

  locs <- as.matrix(data.frame(
    sdimx = coords_indexed[rownames(obj_s@meta.data), "x"],
    sdimy = coords_indexed[rownames(obj_s@meta.data), "y"],
    row.names = rownames(obj_s@meta.data)
  ))

  valid <- !is.na(locs[, 1]) & !is.na(locs[, 2])
  if (sum(!valid) > 0) {
    message("  Removing ", sum(!valid), " cells with NA coordinates")
    gcm   <- gcm[, valid]
    locs  <- locs[valid, ]
    obj_s <- obj_s[, valid]
  }

  mdata <- obj_s@meta.data[, c("sample", "condition", "seurat_clusters_res0.5")]
  se <- SpatialExperiment(assay = list(counts = as.matrix(gcm)),
                           spatialCoords = locs, colData = mdata, sample_id = sample_name)
  aname <- "normcounts"
  assay(se, aname) <- as.matrix(gcm)

  k_geom <- 30 ; compute_agf <- TRUE ; lambda <- 0.8 ; npcs <- 19 ; SEED <- 1997

  message("  Computing BANKSY neighborhood matrix...")
  se <- computeBanksy(se, assay_name = aname, compute_agf = compute_agf, k_geom = k_geom)

  message("  Running BANKSY PCA...")
  se <- runBanksyPCA(se, use_agf = TRUE, lambda = lambda, npcs = npcs, seed = SEED)

  message("  Skipping Harmony (single sample) — using BANKSY PCA directly...")
  reducedDim(se, "Harmony_BANKSY") <- reducedDim(se, paste0("PCA_M1_lam", lambda))

  message("  Running BANKSY UMAP...")
  se <- runBanksyUMAP(se, use_agf = TRUE, lambda = lambda, npcs = npcs)
  se <- runBanksyUMAP(se, dimred = "Harmony_BANKSY")

  message("  Clustering...")
  se <- clusterBanksy(se, dimred = "Harmony_BANKSY", resolution = c(0.1, 0.2, 0.3), seed = SEED)

  cluster_col <- clusterNames(se)[grep("res0.2", clusterNames(se))]
  if (length(cluster_col) == 0) cluster_col <- clusterNames(se)[2]

  coords_df         <- as.data.frame(spatialCoords(se))
  coords_df$cluster <- colData(se)[[cluster_col]]

  p_spatial <- ggplot(coords_df, aes(x = sdimx, y = sdimy, color = cluster)) +
    geom_point(size = 0.3, alpha = 0.7) + coord_equal() +
    labs(title = paste("BANKSY domains —", sample_name),
         subtitle = "lambda=0.8 | k_geom=30 | resolution=0.2",
         x = "X (µm)", y = "Y (µm)") +
    theme_classic() + guides(color = guide_legend(override.aes = list(size = 3)))

  ggsave(file.path(cluster_dir, paste0("BANKSY_spatial_", sample_name, ".pdf")),
         p_spatial, width = 8, height = 7)

  banksy_clusters <- colData(se)[[cluster_col]]
  names(banksy_clusters) <- colnames(se)
  return(list(se = se, clusters = banksy_clusters))
}

banksy_mock  <- run_banksy_for_sample(obj, "mock_6wpi",  OBJ_DIR, CLUSTER_DIR)
banksy_lcmv1 <- run_banksy_for_sample(obj, "LCMV_1wpi",  OBJ_DIR, CLUSTER_DIR)
banksy_lcmv3 <- run_banksy_for_sample(obj, "LCMV_3wpi",  OBJ_DIR, CLUSTER_DIR)
banksy_lcmv6 <- run_banksy_for_sample(obj, "LCMV_6wpi",  OBJ_DIR, CLUSTER_DIR)

# Clustree BANKSY (comme Seurat) pour vérifier la stabilité des domaines
make_banksy_clustree <- function(se_obj, sample_name, cluster_dir) {
  cns <- clusterNames(se_obj)
  cns_res <- cns[grepl("res", cns)]

  if (length(cns_res) < 2) {
    message("BANKSY clustree skipped for ", sample_name, " (need >=2 resolution columns)")
    return(NULL)
  }

  res_vals <- sub(".*res([0-9]+\\.?[0-9]*).*$", "\\1", cns_res)
  keep <- grepl("^[0-9]+\\.?[0-9]*$", res_vals)
  cns_res <- cns_res[keep]
  res_vals <- as.numeric(res_vals[keep])

  if (length(cns_res) < 2) {
    message("BANKSY clustree skipped for ", sample_name, " (resolution parsing failed)")
    return(NULL)
  }

  ord <- order(res_vals)
  cns_res <- cns_res[ord]
  res_vals <- res_vals[ord]

  clust_df <- as.data.frame(colData(se_obj)[, cns_res, drop = FALSE])
  new_names <- paste0("banksy_res", res_vals)
  colnames(clust_df) <- make.unique(new_names)
  clust_df[] <- lapply(clust_df, as.factor)

  p_tree <- clustree(clust_df, prefix = "banksy_res") +
    labs(title = paste0("BANKSY clustree — ", sample_name))

  out_file <- file.path(cluster_dir, paste0("BANKSY_clustree_", sample_name, ".pdf"))
  ggsave(out_file, p_tree, width = 12, height = 14)
  message("BANKSY clustree saved: ", out_file)
  p_tree
}

p_tree_mock  <- make_banksy_clustree(banksy_mock$se,  "mock_6wpi", CLUSTER_DIR)
p_tree_lcmv1 <- make_banksy_clustree(banksy_lcmv1$se, "LCMV_1wpi", CLUSTER_DIR)
p_tree_lcmv3 <- make_banksy_clustree(banksy_lcmv3$se, "LCMV_3wpi", CLUSTER_DIR)
p_tree_lcmv6 <- make_banksy_clustree(banksy_lcmv6$se, "LCMV_6wpi", CLUSTER_DIR)

trees <- list(p_tree_mock, p_tree_lcmv1, p_tree_lcmv3, p_tree_lcmv6)
trees <- trees[!vapply(trees, is.null, logical(1))]
if (length(trees) > 0) {
  p_banksy_clustree <- wrap_plots(trees, ncol = 2)
  ggsave(file.path(CLUSTER_DIR, "BANKSY_clustree_all_samples.pdf"),
         p_banksy_clustree, width = 18, height = 14)
  message("BANKSY combined clustree saved.")
}

obj@meta.data$banksy_domain <- NA_character_
obj@meta.data[names(banksy_mock$clusters),  "banksy_domain"] <- paste0("mock_",  banksy_mock$clusters)
obj@meta.data[names(banksy_lcmv1$clusters), "banksy_domain"] <- paste0("lcmv1_", banksy_lcmv1$clusters)
obj@meta.data[names(banksy_lcmv3$clusters), "banksy_domain"] <- paste0("lcmv3_", banksy_lcmv3$clusters)
obj@meta.data[names(banksy_lcmv6$clusters), "banksy_domain"] <- paste0("lcmv6_", banksy_lcmv6$clusters)

saveRDS(banksy_mock$se,  file.path(OBJ_DIR, "03_banksy_mock_6wpi.rds"))
saveRDS(banksy_lcmv1$se, file.path(OBJ_DIR, "03_banksy_LCMV_1wpi.rds"))
saveRDS(banksy_lcmv3$se, file.path(OBJ_DIR, "03_banksy_LCMV_3wpi.rds"))
saveRDS(banksy_lcmv6$se, file.path(OBJ_DIR, "03_banksy_LCMV_6wpi.rds"))
message("BANKSY objects saved.")


# -------------------------------------------------------
# DEG entre domaines BANKSY
# -------------------------------------------------------
run_banksy_deg <- function(seurat_obj, sample_name, prefix, cluster_dir) {
  obj_sub <- subset(seurat_obj, subset = sample == sample_name)
  obj_sub <- subset(obj_sub, subset = !is.na(banksy_domain))
  obj_sub@meta.data$banksy_short <- paste0("Domain_",
    gsub(paste0(prefix, "_"), "", obj_sub@meta.data$banksy_domain))
  Idents(obj_sub) <- "banksy_short"

  SeuratObject::DefaultAssay(obj_sub) <- "Vizgen"
  obj_sub <- NormalizeData(obj_sub, verbose = FALSE)
  obj_sub <- JoinLayers(obj_sub, assay = "Vizgen")

  markers <- FindAllMarkers(obj_sub, assay = "Vizgen", only.pos = TRUE,
                             min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
  write.csv(markers, file.path(cluster_dir, paste0("banksy_domain_markers_", prefix, ".csv")),
            row.names = FALSE)
  message("BANKSY markers saved for: ", sample_name)
  return(list(obj = obj_sub, markers = markers))
}

res_mock  <- run_banksy_deg(obj, "mock_6wpi",  "mock",  CLUSTER_DIR)
res_lcmv1 <- run_banksy_deg(obj, "LCMV_1wpi",  "lcmv1", CLUSTER_DIR)
res_lcmv3 <- run_banksy_deg(obj, "LCMV_3wpi",  "lcmv3", CLUSTER_DIR)
res_lcmv6 <- run_banksy_deg(obj, "LCMV_6wpi",  "lcmv6", CLUSTER_DIR)

for (res in list(
  list(res = res_mock,  name = "mock_6wpi",  prefix = "mock"),
  list(res = res_lcmv1, name = "LCMV_1wpi",  prefix = "lcmv1"),
  list(res = res_lcmv3, name = "LCMV_3wpi",  prefix = "lcmv3"),
  list(res = res_lcmv6, name = "LCMV_6wpi",  prefix = "lcmv6")
)) {
  markers <- res$res$markers
  if (nrow(markers) > 0) {
    top2 <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 2) %>%
      pull(gene) %>% unique()
    make_pheatmap(res$res$obj, "banksy_short", marker_genes = top2,
                  title = paste("DEGs par domaine BANKSY —", res$name),
                  filename = file.path(CLUSTER_DIR, paste0("heatmap_banksy_", res$prefix, ".pdf")))
  }
}

# Comparaison résolutions BANKSY — mock
se_mock <- banksy_mock$se
plots <- lapply(c("res0.1", "res0.2", "res0.3"), function(res_str) {
  col <- clusterNames(se_mock)[grep(res_str, clusterNames(se_mock))]
  if (length(col) == 0) return(NULL)
  coords_df <- as.data.frame(spatialCoords(se_mock))
  coords_df$cluster <- colData(se_mock)[[col]]
  n_clusters <- length(unique(coords_df$cluster))
  ggplot(coords_df, aes(x = sdimx, y = sdimy, color = cluster)) +
    geom_point(size = 0.2, alpha = 0.6) + coord_equal() +
    labs(title = paste0("BANKSY Mock — ", res_str, " (", n_clusters, " domaines)")) +
    theme_classic()
})
p_compare <- wrap_plots(plots, ncol = 3)
ggsave(file.path(CLUSTER_DIR, "BANKSY_resolution_comparison_mock.pdf"),
       p_compare, width = 21, height = 7)
message("BANKSY resolution comparison saved.")


# =============================================================
# SAVE FINAL OBJECT
# =============================================================
saveRDS(obj, file.path(OBJ_DIR, "03_all_clustered.rds"))

message("\nDone.")
message("  03_all_clustered.rds")
message("  03_cluster2_subclustered.rds")
message("  03_banksy_mock_6wpi.rds / LCMV_1wpi / LCMV_3wpi / LCMV_6wpi")