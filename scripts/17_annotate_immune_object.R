#!/usr/bin/env Rscript
# =============================================================
# Script: 17_annotate_immune_object.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Description:
#   Ajoute les annotations manuelles dans le colData de l'objet
#   sous-clustering immune.
#
#   Input:  objects/07_immune_banksy_lam02_immune_acod1.rds
#           outputs/banksy/immune_acod1/annotation_map_immune_res03.csv
#             (fichier manuel, format semicolon: cluster;annotation)
#   Output: objects/08_immune_annotated_lam02_res03.rds
#
#   NOTE: Ne jamais regénérer annotation_map_immune_res03.csv depuis ce
#         script. Ce fichier est la source de vérité des annotations
#         manuelles — l'éditer directement si besoin.
# =============================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(dplyr)
})

base_path <- normalizePath(".")  # Run this script from the project root directory
setwd(base_path)

# Palette et ordre des annotations (fichier partagé)
source("scripts/00_palette.R")

# -------------------------------------------------------
# Charger l'objet
# -------------------------------------------------------
message("Loading: objects/07_immune_banksy_lam02_immune_acod1.rds")
se <- readRDS(file.path("objects", "07_immune_banksy_lam02_immune_acod1.rds"))
message("  ", ncol(se), " cells")

cd <- as.data.frame(colData(se))
cat("Colonnes disponibles:\n")
print(colnames(cd))

# -------------------------------------------------------
# Trouver la colonne de clustering res=0.3
# -------------------------------------------------------
cl_cols <- colnames(cd)[grepl("res0\\.?3$", colnames(cd))]
cat("\nColonnes res0.3 trouvées:\n")
print(cl_cols)

# Choisir celle avec ~17 clusters
cl_col <- NULL
for (col in cl_cols) {
  n <- length(unique(cd[[col]]))
  cat(col, "→", n, "clusters\n")
  if (n >= 15 && n <= 20) cl_col <- col
}

if (is.null(cl_col)) {
  stop("Aucune colonne avec ~17 clusters trouvée à res=0.3. Vérifier manuellement.")
}
message("\nColonne sélectionnée: ", cl_col)
domain_labels <- paste0("Domain_", as.character(cd[[cl_col]]))

# -------------------------------------------------------
# Charger les annotations manuelles (source de vérité)
# -------------------------------------------------------
annot_file <- file.path("outputs", "banksy", "immune_acod1", "annotation_map_immune_res03.csv")
if (!file.exists(annot_file)) {
  stop("Fichier d'annotations introuvable: ", annot_file,
       "\nCe fichier doit être maintenu manuellement (Domain_X;nom biologique).")
}
mapping_df <- read.table(annot_file, sep = ";", header = TRUE,
                         stringsAsFactors = FALSE, strip.white = TRUE)
stopifnot(all(c("cluster", "annotation") %in% colnames(mapping_df)))
anno_lookup <- setNames(mapping_df$annotation, mapping_df$cluster)
message("Annotations manuelles chargées (", nrow(mapping_df), " entrées): ", annot_file)

# -------------------------------------------------------
# Ajouter cell_type et banksy_domain_immune dans colData
# -------------------------------------------------------
cell_type <- ifelse(
  !is.na(anno_lookup[domain_labels]) & anno_lookup[domain_labels] != "",
  anno_lookup[domain_labels],
  domain_labels
)

se$cell_type     <- cell_type
se$banksy_domain_immune <- domain_labels

cat("\nDistribution cell_type:\n")
print(sort(table(se$cell_type), decreasing = TRUE))

# -------------------------------------------------------
# Sauvegarder
# -------------------------------------------------------
out_file <- file.path("objects", "08_immune_annotated_lam02_res03.rds")
saveRDS(se, out_file)
message("\nSaved: ", out_file)
message("=== Done ===\n")
