# =============================================================
# Script: 05c_umap_annotated_from_ncells_csv_res09.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-04
#
# Description:
#   Génère un UMAP BANKSY annoté à partir du CSV long
#   ncells_by_sample_lam02_res09_joint_long.csv.
#
#   Le script :
#     1) repart de se_joint
#     2) récupère la colonne BANKSY pour lambda=0.2, res=0.9
#     3) lit le CSV long par sample
#     4) reconstruit un mapping unique banksy_domain -> annotation
#     5) annote le UMAP
#     6) sauve :
#          - un UMAP global annoté
#          - un UMAP splitté par sample
#
# Input:
#   objects/04_banksy_joint_lam08_after_bloc3.rds
#   ncells_by_sample_lam02_res09_joint_long.csv
#
# Output:
#   outputs/banksy/umap_annotated/
#     UMAP_annotated_from_ncellscsv_lam02_res09.jpg
#     UMAP_annotated_from_ncellscsv_lam02_res09.pdf
#     UMAP_annotated_from_ncellscsv_lam02_res09_split.jpg
#     UMAP_annotated_from_ncellscsv_lam02_res09_split.pdf
#     annotation_map_lam02_res09.csv
# =============================================================

# -------------------------------------------------------
# 1. Libraries
# -------------------------------------------------------
library(Banksy)
library(SpatialExperiment)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(dplyr)
library(readr)
library(pheatmap)

# -------------------------------------------------------
# 2. Parameters
# -------------------------------------------------------
OBJ_DIR   <- "objects"
OUT_DIR   <- file.path("outputs", "banksy", "umap_annotated")
CSV_PATH  <- "ncells_by_sample_lam02_res09_joint_long.csv"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

LAM          <- 0.2
RES_TARGET   <- 0.9
SEED         <- 1997
SAMPLE_ORDER <- c("mock_6wpi", "LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")

# Palette et ordre des annotations (fichier partagé)
source("scripts/00_palette.R")

# Ordre des colonnes pour la heatmap DEG — aligné sur ANNOTATION_ORDER
CLUSTER_ORDER <- ANNOTATION_ORDER

# -------------------------------------------------------
# 3. Helpers
# -------------------------------------------------------
find_cl_col <- function(se, lam, res) {
  all_cols <- clusterNames(se)
  lam_cols <- all_cols[grep(
    paste0("lam", gsub("\\.", "\\\\.", as.character(lam))),
    all_cols
  )]
  if (length(lam_cols) == 0) return(NULL)

  res_vals <- suppressWarnings(as.numeric(sub(".*_res", "", lam_cols)))
  idx <- which(!is.na(res_vals) & abs(res_vals - res) < 1e-8)
  if (length(idx) == 0) return(NULL)

  lam_cols[idx[1]]
}

save_plot_pdf_jpg <- function(plot_obj, base_path, width, height) {
  ggsave(
    filename = paste0(base_path, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height
  )
  ggsave(
    filename = paste0(base_path, ".jpg"),
    plot = plot_obj,
    device = "jpeg",
    width = width,
    height = height,
    dpi = 300,
    quality = 95
  )
  message("Saved: ", basename(base_path), ".pdf / .jpg")
}

save_heatmap_jpg <- function(avg_mat, title, out_jpg, gaps_row = NULL, labels_row = NULL, w = 12, h = 10) {
  if (is.null(avg_mat) || nrow(avg_mat) == 0 || ncol(avg_mat) == 0) {
    message("Heatmap ignorée: matrice vide")
    return(invisible(NULL))
  }

  avg_scaled <- t(scale(t(avg_mat)))
  avg_scaled <- pmin(pmax(avg_scaled, -2), 2)
  avg_scaled[is.nan(avg_scaled)] <- 0

  n_genes <- nrow(avg_scaled)
  n_groups <- ncol(avg_scaled)

  # Auto-ajustement pour éviter le chevauchement des labels.
  w_auto <- max(w, 7 + 0.6 * n_groups)
  h_auto <- max(h, 5 + 0.2 * n_genes)
  fs_row <- max(5, min(9, 260 / max(25, n_genes)))
  fs_col <- max(7, min(10, 130 / max(8, n_groups)))

  labels_col_wrapped <- vapply(
    colnames(avg_scaled),
    function(x) paste(strwrap(x, width = 14), collapse = "\n"),
    character(1)
  )

  pheatmap(
    avg_scaled,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    labels_row = labels_row,
    labels_col = labels_col_wrapped,
    fontsize_row = fs_row,
    fontsize_col = fs_col,
    angle_col = 45,
    border_color = NA,
    main = title,
    gaps_row = gaps_row,
    legend = FALSE,
    filename = out_jpg,
    width = w_auto,
    height = h_auto
  )

  message("Saved: ", basename(out_jpg))
}

# -------------------------------------------------------
# 4. Load se_joint
# -------------------------------------------------------
message("Loading se_joint...")
se_joint <- readRDS(file.path(OBJ_DIR, "04_banksy_joint_lam08_after_bloc3.rds"))
message("  se_joint: ", ncol(se_joint), " cells")
message("  reducedDims: ", paste(reducedDimNames(se_joint), collapse = ", "))

# -------------------------------------------------------
# 5. Ensure UMAP exists
# -------------------------------------------------------
harmony_name <- paste0("Harmony_lam", LAM)
umap_name    <- paste0("UMAP_Harmony_lam", LAM)

if (!harmony_name %in% reducedDimNames(se_joint)) {
  stop("ReducedDim introuvable : ", harmony_name)
}

if (!umap_name %in% reducedDimNames(se_joint)) {
  message("UMAP absent, calcul via runBanksyUMAP...")
  se_joint <- runBanksyUMAP(se_joint, dimred = harmony_name, seed = SEED)

  new_umaps <- reducedDimNames(se_joint)[grep("UMAP", reducedDimNames(se_joint))]
  latest_umap <- tail(new_umaps, 1)

  if (latest_umap != umap_name) {
    reducedDim(se_joint, umap_name) <- reducedDim(se_joint, latest_umap)
  }

  saveRDS(se_joint, file.path(OBJ_DIR, "04_banksy_joint_lam08_after_bloc3.rds"))
  message("  UMAP calculé et objet resauvegardé")
} else {
  message("UMAP déjà présent : ", umap_name)
}

# -------------------------------------------------------
# 6. Find BANKSY column for res=0.9
# -------------------------------------------------------
col <- find_cl_col(se_joint, LAM, RES_TARGET)
if (is.null(col)) {
  stop("Colonne BANKSY introuvable pour lambda=", LAM, " et res=", RES_TARGET)
}
message("Colonne BANKSY utilisée : ", col)

# -------------------------------------------------------
# 7. Read annotation CSV long
# -------------------------------------------------------
message("Loading annotation CSV...")
annot_long <- read_delim(
  CSV_PATH,
  delim = ";",
  locale = locale(decimal_mark = "."),
  show_col_types = FALSE,
  trim_ws = TRUE
) %>%
  select(-matches("^Unnamed")) %>%
  mutate(
    banksy_domain = as.character(banksy_domain),
    annotation = trimws(as.character(annotation)),
    sample = as.character(sample)
  )

message("  n lignes CSV : ", nrow(annot_long))

required_cols <- c("sample", "banksy_domain", "annotation")
missing_cols <- setdiff(required_cols, colnames(annot_long))
if (length(missing_cols) > 0) {
  stop("Colonnes manquantes dans le CSV : ", paste(missing_cols, collapse = ", "))
}

# -------------------------------------------------------
# 8. Build unique mapping Domain -> annotation
# -------------------------------------------------------
annotation_map <- annot_long %>%
  filter(!is.na(annotation), annotation != "") %>%
  distinct(banksy_domain, annotation)

# Vérifier s'il y a des domaines annotés plusieurs fois différemment
duplicates <- annotation_map %>%
  count(banksy_domain) %>%
  filter(n > 1)

if (nrow(duplicates) > 0) {
  message("\nATTENTION : certains domaines ont plusieurs annotations :")
  print(duplicates)

  print(
    annotation_map %>%
      semi_join(duplicates, by = "banksy_domain") %>%
      arrange(banksy_domain)
  )

  stop("Mapping ambigu dans le CSV. Corriger les annotations avant de continuer.")
}

write_csv(
  annotation_map,
  file.path(OUT_DIR, "annotation_map_lam02_res09.csv")
)
message("Saved: annotation_map_lam02_res09.csv")

# -------------------------------------------------------
# 9. Build UMAP dataframe
# -------------------------------------------------------
umap_coords <- as.data.frame(reducedDim(se_joint, umap_name))
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_coords$cell_id   <- colnames(se_joint)
umap_coords$sample    <- factor(as.character(se_joint$sample), levels = SAMPLE_ORDER)

umap_coords$banksy_domain <- paste0(
  "Domain_",
  as.character(colData(se_joint)[[col]])
)

umap_coords <- umap_coords %>%
  left_join(annotation_map, by = "banksy_domain")

umap_coords$annotation[is.na(umap_coords$annotation) | umap_coords$annotation == ""] <- "Non annote"

# Ordonner les niveaux selon ANNOTATION_ORDER (inconnus en fin)
annotation_levels <- order_annotations(unique(as.character(umap_coords$annotation)))
umap_coords$annotation <- factor(umap_coords$annotation, levels = annotation_levels)

message("\nRésumé annotations UMAP :")
print(table(umap_coords$annotation, useNA = "ifany"))

# -------------------------------------------------------
# 10. Plot global UMAP
# -------------------------------------------------------
n_types <- n_distinct(umap_coords$annotation)
n_domains <- n_distinct(umap_coords$banksy_domain)

p_global <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = annotation)) +
  geom_point(size = 0.05, alpha = 0.5) +
  scale_color_manual(values = GLOBAL_PALETTE, na.value = "grey70", drop = FALSE, na.translate = FALSE) +
  labs(
    title = paste0("BANKSY UMAP annoté — lambda=", LAM, " | res=", RES_TARGET),
    subtitle = paste0(
      n_types, " types cellulaires | ",
      n_domains, " domaines | ",
      nrow(umap_coords), " cellules"
    ),
    color = "Type cellulaire"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

save_plot_pdf_jpg(
  p_global,
  file.path(OUT_DIR, "UMAP_annotated_from_ncellscsv_lam02_res09"),
  width = 10,
  height = 8
)

# -------------------------------------------------------
# 11. Plot split by sample
# -------------------------------------------------------
p_split <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = annotation)) +
  geom_point(size = 0.05, alpha = 0.5) +
  scale_color_manual(values = GLOBAL_PALETTE, na.value = "grey70", drop = FALSE, na.translate = FALSE) +
  facet_wrap(~ sample, ncol = 2) +
  labs(
    title = paste0("BANKSY UMAP annoté par échantillon — lambda=", LAM, " | res=", RES_TARGET),
    subtitle = paste0(n_types, " types cellulaires"),
    color = "Type cellulaire"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

save_plot_pdf_jpg(
  p_split,
  file.path(OUT_DIR, "UMAP_annotated_from_ncellscsv_lam02_res09_split"),
  width = 12,
  height = 10
)

# -------------------------------------------------------
# 12. DEG + heatmap from annotated labels
# -------------------------------------------------------
message("\n=== DEG + HEATMAP ===")

message("Loading Seurat object...")
obj <- readRDS(file.path(OBJ_DIR, "03_all_clustered.rds"))
message("  obj: ", ncol(obj), " cells")

DefaultAssay(obj) <- "Vizgen"
obj <- NormalizeData(obj, verbose = FALSE)
obj <- JoinLayers(obj, assay = "Vizgen")

# Transfer annotation labels from UMAP table to Seurat cells.
annotation_by_cell <- setNames(as.character(umap_coords$annotation), umap_coords$cell_id)
matched <- intersect(colnames(obj), names(annotation_by_cell))
obj@meta.data[matched, "banksy_annotation"] <- annotation_by_cell[matched]
obj@meta.data[is.na(obj@meta.data$banksy_annotation), "banksy_annotation"] <- "Non annote"

obj_sub <- subset(obj, subset = banksy_annotation != "Non annote")
n_groups <- length(unique(obj_sub@meta.data$banksy_annotation))
if (n_groups < 2) {
  stop("Pas assez de groupes annotés pour DEG (n<2).")
}
obj_sub@meta.data$banksy_annotation <- factor(obj_sub@meta.data$banksy_annotation)
Idents(obj_sub) <- "banksy_annotation"

deg_csv <- file.path(OUT_DIR, "DEG_annotated_lam02_res09.csv")
heat_jpg <- file.path(OUT_DIR, "heatmap_annotated_lam02_res09.jpg")

message("Running FindAllMarkers...")
markers <- FindAllMarkers(
  obj_sub,
  assay = "Vizgen",
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25,
  verbose = FALSE
)
write.csv(markers, deg_csv, row.names = FALSE)
message("Saved: ", basename(deg_csv), " (", nrow(markers), " lignes)")

top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(cluster, desc(avg_log2FC), gene)

cluster_order <- unique(as.character(top5$cluster))
if (!is.null(CLUSTER_ORDER)) {
  cluster_order <- c(intersect(CLUSTER_ORDER, cluster_order), setdiff(cluster_order, CLUSTER_ORDER))
}

avg <- AverageExpression(
  obj_sub,
  assays = "Vizgen",
  layer = "data",
  group.by = "banksy_annotation"
)$Vizgen

cluster_order <- cluster_order[cluster_order %in% colnames(avg)]
avg <- avg[, cluster_order, drop = FALSE]

genes_ok <- top5$gene[top5$gene %in% rownames(avg)]
if (length(genes_ok) > 0) {
  for (top_n in c(5, 2)) {
    topN <- markers %>%
      group_by(cluster) %>%
      slice_max(avg_log2FC, n = top_n, with_ties = FALSE) %>%
      ungroup() %>%
      filter(gene %in% genes_ok, cluster %in% cluster_order) %>%
      mutate(cluster = factor(cluster, levels = cluster_order)) %>%
      arrange(cluster, desc(avg_log2FC), gene)

    heat_mat <- avg[topN$gene, , drop = FALSE]
    row_labels <- make.unique(as.character(topN$gene), sep = "__")
    rownames(heat_mat) <- row_labels

    block_sizes <- as.integer(table(topN$cluster))
    gaps_row <- if (length(block_sizes) > 1) cumsum(block_sizes)[-length(block_sizes)] else NULL

    save_heatmap_jpg(
      heat_mat,
      title = paste0("DEG annotated top", top_n, " (lambda=0.2, res=0.9)"),
      out_jpg = file.path(OUT_DIR, paste0("heatmap_annotated_lam02_res09_top", top_n, ".jpg")),
      gaps_row = gaps_row,
      labels_row = as.character(topN$gene),
      w = 12,
      h = max(10, 0.45 * nrow(heat_mat))
    )
  }
} else {
  message("Heatmap ignorée: aucun gène top5 présent dans AverageExpression")
}

# -------------------------------------------------------
# 13. Done
# -------------------------------------------------------
message("\n=== Done ===")
message("Outputs : ", OUT_DIR)

