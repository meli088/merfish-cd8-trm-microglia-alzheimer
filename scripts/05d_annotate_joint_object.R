#!/usr/bin/env Rscript
# =============================================================
# Script: 05d_annotate_joint_object.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
#
# Description:
#   Ajoute les annotations globales (lambda=0.2, res=0.9) dans le
#   colData de l'objet BANKSY joint, puis sauvegarde un objet annoté.
#
#   La source de vérité des annotations est le CSV long :
#     ncells_by_sample_lam02_res09_joint_long.csv
#   Ne jamais régénérer ce CSV depuis ce script — l'éditer directement.
#
#   Colonnes ajoutées dans colData :
#     banksy_domain   — "Domain_XX" (cluster brut)
#     annotation      — label biologique issu du CSV
#                       (= banksy_domain si domaine non annoté)
#
# Input:
#   objects/04_banksy_joint_lam08_after_bloc3.rds
#   ncells_by_sample_lam02_res09_joint_long.csv
#
# Output:
#   objects/05_joint_annotated_lam02_res09.rds
# =============================================================

suppressPackageStartupMessages({
  library(Banksy)
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(dplyr)
  library(readr)
})

setwd(normalizePath("."))   # Run from project root

# Palette et ordre des annotations (fichier partagé)
source("scripts/00_palette.R")

# -------------------------------------------------------
# Paramètres
# -------------------------------------------------------
OBJ_IN   <- file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds")
CSV_PATH <- "ncells_by_sample_lam02_res09_joint_long.csv"
OBJ_OUT  <- file.path("objects", "05_joint_annotated_lam02_res09.rds")
LAM      <- 0.2
RES      <- 0.9

# -------------------------------------------------------
# Helper : trouver la colonne de clustering BANKSY
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

# -------------------------------------------------------
# Charger l'objet
# -------------------------------------------------------
message("Loading: ", OBJ_IN)
se <- readRDS(OBJ_IN)
message("  ", ncol(se), " cells")

# -------------------------------------------------------
# Trouver la colonne de clustering lambda=0.2, res=0.9
# -------------------------------------------------------
cl_col <- find_cl_col(se, LAM, RES)
if (is.null(cl_col)) {
  stop("Colonne BANKSY introuvable pour lambda=", LAM, " et res=", RES,
       "\nColonnes disponibles : ", paste(clusterNames(se), collapse = ", "))
}
message("Colonne BANKSY utilisée : ", cl_col)

domain_labels <- paste0("Domain_", as.character(colData(se)[[cl_col]]))
message("  Domaines distincts : ", length(unique(domain_labels)))

# -------------------------------------------------------
# Charger les annotations depuis le CSV long
# -------------------------------------------------------
if (!file.exists(CSV_PATH)) {
  stop("CSV introuvable : ", CSV_PATH)
}

annot_long <- read_delim(
  CSV_PATH,
  delim = ";",
  locale = locale(decimal_mark = "."),
  show_col_types = FALSE,
  trim_ws = TRUE
) %>%
  select(-matches("^Unnamed")) %>%
  mutate(
    banksy_domain = trimws(as.character(banksy_domain)),
    annotation    = trimws(as.character(annotation))
  )

message("CSV chargé (", nrow(annot_long), " lignes)")

# Construire le mapping unique domaine -> annotation
# (prendre la première annotation non-vide par domaine)
annotation_map <- annot_long %>%
  filter(!is.na(annotation), annotation != "", annotation != "NA") %>%
  distinct(banksy_domain, annotation)

# Vérifier les domaines annotés de façon contradictoire
duplicates <- annotation_map %>% count(banksy_domain) %>% filter(n > 1)
if (nrow(duplicates) > 0) {
  message("\nATTENTION : domaines avec plusieurs annotations différentes :")
  print(annotation_map %>% semi_join(duplicates, by = "banksy_domain") %>% arrange(banksy_domain))
  stop("Mapping ambigu dans le CSV — corriger avant de continuer.")
}

anno_lookup <- setNames(annotation_map$annotation, annotation_map$banksy_domain)
message("Mapping chargé : ", length(anno_lookup), " domaines annotés")

# Vérifier la couverture
domains_in_obj   <- unique(domain_labels)
domains_in_csv   <- names(anno_lookup)
domains_missing  <- setdiff(domains_in_obj, domains_in_csv)
if (length(domains_missing) > 0) {
  message("Domaines non annotés dans le CSV (garderont leur label brut) : ",
          paste(sort(domains_missing), collapse = ", "))
}

# -------------------------------------------------------
# Ajouter les colonnes dans colData
# -------------------------------------------------------
annotation <- ifelse(
  !is.na(anno_lookup[domain_labels]) & anno_lookup[domain_labels] != "",
  anno_lookup[domain_labels],
  domain_labels
)

se$banksy_domain <- domain_labels
se$annotation    <- annotation

cat("\nDistribution annotation :\n")
print(sort(table(se$annotation), decreasing = TRUE))

# -------------------------------------------------------
# Sauvegarder
# -------------------------------------------------------
saveRDS(se, OBJ_OUT)
message("\nSaved: ", OBJ_OUT)
message("=== Done ===\n")
