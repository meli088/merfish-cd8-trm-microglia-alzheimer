#!/usr/bin/env Rscript
# =============================================================
# Script: 13_immune_banksy.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Description: BANKSY clustering on the joint immune subset.
#              Strictly mirrors old_scripts/04c_banksy_no_neurons.R:
#              - computeBanksy per sample (preserves spatial neighborhoods)
#              - cbind() all samples
#              - runBanksyPCA() on merged object
#              - RunHarmony() directly on PCA matrix (data_mat API)
#              - clusterBanksy() for res 0.1 – 0.5
#              - Clustree + per-resolution composition CSVs
#
# Input:  objects/06_immune_normalized_lam02_<label>_spe.rds
# Output: objects/07_immune_banksy_lam02_<label>.rds
#         outputs/banksy/<label>/banksy/
# =============================================================

library(Banksy)
library(SpatialExperiment)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(clustree)
library(ggplot2)
library(patchwork)
library(harmony)
library(dplyr)
library(tidyr)

base_path <- normalizePath(".")  # Run this script from the project root directory
setwd(base_path)

# -------------------------------------------------------
# Parameters (mirrors 04c_banksy_no_neurons.R)
# -------------------------------------------------------
IMMUNE_LABEL <- "Immune (Acod1)"  # default; override with --label
FOLDER_NAME  <- NULL   # optional: override folder/file slug with --folder microglia
SEED     <- 1997
K_GEOM   <- 30
LAM      <- 0.2
ALL_RES  <- seq(0.1, 0.5, by = 0.1)

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
    if (a == "--k_geom" && (i + 1) <= length(args)) {
      K_GEOM <- as.integer(args[i + 1])
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
BANKSY_DIR <- file.path(label_root, "banksy")
if (!dir.exists(BANKSY_DIR)) dir.create(BANKSY_DIR, recursive = TRUE, showWarnings = FALSE)

# (SEED, K_GEOM, LAM, ALL_RES defaults set above before CLI parsing)

# -------------------------------------------------------
# Helpers (mirrors 04c_banksy_no_neurons.R)
# -------------------------------------------------------
find_cl_col <- function(se, lam, res) {
  all_cols <- clusterNames(se)
  lam_cols <- all_cols[grep(
    paste0("lam", gsub("\\.", "\\\\.", as.character(lam))), all_cols
  )]
  if (length(lam_cols) == 0) return(NULL)
  res_vals <- suppressWarnings(as.numeric(sub(".*_res", "", lam_cols)))
  idx <- which(!is.na(res_vals) & abs(res_vals - res) < 1e-8)
  if (length(idx) == 0) return(NULL)
  lam_cols[idx[1]]
}

find_pca_dim <- function(se, lam) {
  dims     <- reducedDimNames(se)
  lam_chr  <- as.character(lam)
  lam_nd   <- gsub("\\.", "", lam_chr)
  patterns <- c(
    paste0("^PCA_M1_lam", gsub("\\.", "\\\\.", lam_chr), "$"),
    paste0("^PCA_M1_lam", lam_nd, "$"),
    paste0("^PCA_M0_lam", gsub("\\.", "\\\\.", lam_chr), "$"),
    paste0("^PCA_M0_lam", lam_nd, "$"),
    paste0("^PCA.*lam",  gsub("\\.", "\\\\.", lam_chr), "$"),
    paste0("^PCA.*lam",  lam_nd, "$")
  )
  for (pat in patterns) {
    hit <- dims[grepl(pat, dims)]
    if (length(hit) > 0) return(hit[1])
  }
  pca_hits <- dims[grepl("^PCA", dims)]
  if (length(pca_hits) > 0) return(tail(pca_hits, 1))
  NULL
}

already_clustered <- function(se, lam, res) {
  !is.null(find_cl_col(se, lam, res))
}

run_gc <- function(step = "") {
  invisible(gc(verbose = FALSE))
  if (nzchar(step)) message("  GC: ", step)
}

save_fig <- function(p, basename, width = 10, height = 8) {
  ggsave(filename = file.path(BANKSY_DIR, paste0(basename, ".jpg")),
         plot = p, device = "jpeg", width = width, height = height,
         dpi = 300, quality = 95)
  message("  Saved: ", basename, ".jpg")
}

# -------------------------------------------------------
# 1. Load SPE from script 12
# -------------------------------------------------------
in_candidates <- c(
  file.path("objects", paste0("06_immune_normalized_lam02_", safe_label, "_spe.rds")),
  file.path("objects", "06_immune_normalized_lam02_sce.rds")
)
in_file <- in_candidates[file.exists(in_candidates)][1]
if (is.na(in_file)) stop("Normalized SPE not found for label: ", IMMUNE_LABEL)

message("Loading immune normalized SPE: ", in_file)
se_all <- readRDS(in_file)

sample_vec  <- as.character(SummarizedExperiment::colData(se_all)$sample)
uniq_samples <- sort(unique(sample_vec))
message("Cells: ", ncol(se_all), " | Genes: ", nrow(se_all))
message("Joint immune subset across ", length(uniq_samples), " samples: ",
        paste(uniq_samples, collapse = ", "))

# -------------------------------------------------------
# 2. Clean inherited BANKSY assays / reducedDims / clusters
#    (object comes from script 12 which may carry PCA/HARMONY/UMAP
#     from Seurat — those were computed across all samples jointly
#     and are not per-sample spatial neighborhoods)
# -------------------------------------------------------
message("\nCleaning inherited reductions and BANKSY assays...")

for (d in reducedDimNames(se_all)) reducedDim(se_all, d) <- NULL
for (cl in clusterNames(se_all))   colData(se_all)[[cl]] <- NULL
for (a in c("H0", "H1")) {
  if (a %in% assayNames(se_all)) {
    assay(se_all, a) <- NULL
    message("  Removed assay: ", a)
  }
}

assay_name <- if ("counts" %in% assayNames(se_all)) "counts" else assayNames(se_all)[1]
message("Using assay: ", assay_name)
message("assays after clean: ", paste(assayNames(se_all), collapse = ", "))

# -------------------------------------------------------
# 3. computeBanksy per sample (mirrors 04c exactly)
#    Each tissue section has its own x/y coordinate space.
#    Computing neighbors across sections is biologically wrong.
# -------------------------------------------------------
lam_str      <- sub("\\.", "", as.character(LAM))
harmony_name <- paste0("Harmony_lam", LAM)

message("\n--- Lambda = ", LAM, " ---")
message("computeBanksy per sample (k_geom=", K_GEOM, ")...")

se_list <- lapply(uniq_samples, function(sname) {
  idx  <- se_all$sample == sname
  se_s <- se_all[, idx]
  message("  ", sname, ": ", sum(idx), " cells")
  Banksy::computeBanksy(
    se_s,
    compute_agf = TRUE,
    k_geom      = K_GEOM,
    assay_name  = assay_name
  )
})

se_merged <- do.call(cbind, se_list)
rm(se_list); run_gc("after computeBanksy")
message("Merged: ", ncol(se_merged), " cells")

# -------------------------------------------------------
# 4. runBanksyPCA on merged object (mirrors 04c)
# -------------------------------------------------------
message("\nrunBanksyPCA (lambda=", LAM, ")...")
se_merged <- Banksy::runBanksyPCA(
  se_merged,
  lambda = LAM,
  npcs   = 40,
  seed   = SEED
)

pca_name <- find_pca_dim(se_merged, LAM)
if (is.null(pca_name)) {
  stop("No PCA reducedDim found after runBanksyPCA. Available: ",
       paste(reducedDimNames(se_merged), collapse = ", "))
}
message("BANKSY PCA: ", pca_name, " (",
        nrow(reducedDim(se_merged, pca_name)), " x ",
        ncol(reducedDim(se_merged, pca_name)), ")")

# -------------------------------------------------------
# 5. RunHarmony directly on PCA matrix (mirrors 04c)
#    Using data_mat API to avoid SCE reducedDim dispatch issues.
# -------------------------------------------------------
message("\nRunHarmony on BANKSY PCA (sample only — slide confounded with timepoint)...")

pca_mat <- reducedDim(se_merged, pca_name)
cd_df   <- as.data.frame(SummarizedExperiment::colData(se_merged))

batch_vars <- "sample"
message("Correcting for: ", paste(batch_vars, collapse = " + "))


harmony_emb <- harmony::RunHarmony(
  data_mat  = pca_mat,
  meta_data = cd_df,
  vars_use  = batch_vars,
  verbose   = FALSE
)
reducedDim(se_merged, harmony_name) <- harmony_emb
message("Harmony done: ", harmony_name)

# -------------------------------------------------------
# 6. clusterBanksy for res 0.1 – 0.5 (mirrors 04c logic)
# -------------------------------------------------------
message("\nClustering resolutions: ", paste(ALL_RES, collapse = ", "), "...")

for (res in ALL_RES) {
  if (already_clustered(se_merged, LAM, res)) {
    message("  res=", res, " — already done, skip")
    next
  }
  run_gc(paste0("before res=", res))

  se_merged <- tryCatch(
    Banksy::clusterBanksy(
      se_merged,
      dimred     = harmony_name,
      resolution = res,
      seed       = SEED
    ),
    error = function(e) {
      message("  res=", res, " — ERROR: ", conditionMessage(e))
      se_merged
    }
  )
  run_gc(paste0("after res=", res))
  message("  res=", res, " — OK")
}

message("Clustering columns: ", length(clusterNames(se_merged)))

# -------------------------------------------------------
# 7. Save final object
# -------------------------------------------------------
final_out <- file.path("objects", paste0("07_immune_banksy_lam02_", safe_label, ".rds"))
saveRDS(se_merged, final_out)
message("\nSaved: ", final_out)

# -------------------------------------------------------
# 8. Clustree (mirrors 04c)
# -------------------------------------------------------
all_cl_cols <- clusterNames(se_merged)
lam_cols    <- all_cl_cols[grepl(
  paste0("lam", gsub("\\.", "\\\\.", as.character(LAM))), all_cl_cols
)]

if (length(lam_cols) >= 2) {
  message("\n=== CLUSTREE ===")
  res_vals <- suppressWarnings(as.numeric(sub(".*_res", "", lam_cols)))
  ord      <- order(res_vals)
  lam_cols <- lam_cols[ord]
  res_vals <- res_vals[ord]

  cl_df <- as.data.frame(
    lapply(lam_cols, function(col) as.character(colData(se_merged)[[col]]))
  )
  colnames(cl_df) <- paste0("res", res_vals)

  p_ctree <- clustree(cl_df, prefix = "res")
  save_fig(p_ctree, paste0("clustree_BANKSY_lam", lam_str), width = 14, height = 22)
}

# -------------------------------------------------------
# 9. Per-resolution composition tables
# -------------------------------------------------------
cd <- as.data.frame(SummarizedExperiment::colData(se_merged))

for (r in ALL_RES) {
  cl_col <- find_cl_col(se_merged, LAM, r)
  if (is.null(cl_col)) next

  # Build full sample × cluster grid so samples with 0 cells in a
  # cluster still appear (with n_cells=0, pct=0) instead of being dropped.
  out_tab <- cd %>%
    transmute(
      sample  = as.character(sample),
      cluster = paste0("Domain_", as.character(.data[[cl_col]]))
    ) %>%
    count(sample, cluster, name = "n_cells") %>%
    tidyr::complete(sample, cluster, fill = list(n_cells = 0L)) %>%
    group_by(sample) %>%
    mutate(pct = round(100 * n_cells / sum(n_cells), 2)) %>%
    ungroup() %>%
    arrange(sample, cluster)

  write.csv(
    out_tab,
    file.path(BANKSY_DIR, paste0("composition_res", gsub("\\.", "", sprintf("%.1f", r)), ".csv")),
    row.names = FALSE
  )
}

message("\nBANKSY clustering complete.")
message("Outputs in: ", BANKSY_DIR)