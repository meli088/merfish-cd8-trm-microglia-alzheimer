#!/usr/bin/env Rscript
# =============================================================
# Script: 30_microglia_qc_table.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-05
#
# Goal:
#   Compact QC/summary table for the Microglia (P2ry12) population:
#     - n cells per condition
#     - mean / median n genes detected (count > 0) per cell per condition
#     - n genes tested per DEG contrast (from DEG_summary.csv if available)
#
# Input:
#   objects/04_banksy_joint_lam08_after_bloc3.rds
#   ncells_by_sample_lam02_res09_joint_long.csv
#   outputs/banksy/deg_microglia_conditions/DEG_summary.csv  (optional)
#
# Output:
#   outputs/banksy/deg_microglia_conditions/microglia_qc_table.csv
# =============================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(Banksy)
  library(tidyverse)
  library(Matrix)
})

base_path <- normalizePath(".")  # Run this script from the project root directory
setwd(base_path)

out_dir <- file.path("outputs", "banksy", "deg_microglia_conditions")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================
# Parameters (must match script 28 / 29)
# =============================================================

MICROGLIA_LABEL <- "Microglia (P2ry12)"
LAM             <- 0.2
RES_TARGET      <- 0.9
REF_CONDITION   <- "mock_6wpi"
SAMPLE_ORDER    <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi", "mock_6wpi")

SAMPLE_LABELS <- c(
  mock_6wpi = "Mock 6 wpi",
  LCMV_1wpi = "LCMV 1 wpi",
  LCMV_3wpi = "LCMV 3 wpi",
  LCMV_6wpi = "LCMV 6 wpi"
)

# =============================================================
# 1. Load object
# =============================================================

obj_file <- file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds")
stopifnot("Object file not found" = file.exists(obj_file))

message("Loading: ", obj_file)
se <- readRDS(obj_file)
message("  ", ncol(se), " cells | ", nrow(se), " genes")

# =============================================================
# 2. Find cluster column (lambda=0.2, res=0.9)
# =============================================================

find_cl_col <- function(se, lam, res) {
  all_cols <- clusterNames(se)
  lam_str  <- gsub("\\.", "\\\\.", as.character(lam))
  lam_cols <- all_cols[grep(paste0("lam", lam_str), all_cols)]
  if (length(lam_cols) == 0) return(NULL)
  res_vals <- suppressWarnings(as.numeric(sub(".*_res", "", lam_cols)))
  idx      <- which(!is.na(res_vals) & abs(res_vals - res) < 1e-8)
  if (length(idx) == 0) return(NULL)
  lam_cols[idx[1]]
}

cl_col <- find_cl_col(se, LAM, RES_TARGET)
if (is.null(cl_col)) {
  stop("Cluster column not found for lambda=", LAM, " res=", RES_TARGET,
       "\nAvailable: ", paste(clusterNames(se), collapse = ", "))
}
message("  Cluster column: ", cl_col)

# =============================================================
# 3. Reconstruct annotation
# =============================================================

csv_path <- "ncells_by_sample_lam02_res09_joint_long.csv"
stopifnot("Annotation CSV not found" = file.exists(csv_path))

annot_long <- read_delim(
  csv_path,
  delim          = ";",
  locale         = locale(decimal_mark = "."),
  show_col_types = FALSE,
  trim_ws        = TRUE
) %>%
  select(-matches("^Unnamed")) %>%
  mutate(
    banksy_domain = as.character(banksy_domain),
    annotation    = trimws(as.character(annotation))
  )

annotation_map <- annot_long %>%
  filter(!is.na(annotation), annotation != "") %>%
  distinct(banksy_domain, annotation)

domain_labels <- paste0("Domain_", as.character(colData(se)[[cl_col]]))
anno_lookup   <- setNames(annotation_map$annotation, annotation_map$banksy_domain)
annotation    <- ifelse(
  !is.na(anno_lookup[domain_labels]) & anno_lookup[domain_labels] != "",
  anno_lookup[domain_labels],
  "Non annote"
)

if (!MICROGLIA_LABEL %in% annotation) {
  stop("'", MICROGLIA_LABEL, "' not found. Present: ",
       paste(sort(unique(annotation)), collapse = ", "))
}

# =============================================================
# 4. Subset to Microglia
# =============================================================

mg_idx <- which(annotation == MICROGLIA_LABEL)
se_mg  <- se[, mg_idx]
message("  Microglia cells: ", ncol(se_mg))

sample_vec <- factor(
  as.character(colData(se_mg)$sample),
  levels = SAMPLE_ORDER
)

# =============================================================
# 5. Compute nGenes per cell (genes with count > 0)
# =============================================================

assay_name <- if ("counts" %in% assayNames(se_mg)) "counts" else assayNames(se_mg)[1]
cnt        <- assay(se_mg, assay_name)   # genes × cells

# Works for both dense and sparse matrices
if (inherits(cnt, "sparseMatrix")) {
  ngenes_per_cell <- diff(cnt@p)          # column-wise nnz for dgCMatrix
} else {
  ngenes_per_cell <- colSums(cnt > 0)
}

# =============================================================
# 6. QC table per condition
# =============================================================

qc_df <- tibble(
  sample    = as.character(sample_vec),
  n_genes   = ngenes_per_cell
) %>%
  group_by(sample) %>%
  summarise(
    n_cells          = n(),
    mean_genes_det   = round(mean(n_genes, na.rm = TRUE), 1),
    median_genes_det = median(n_genes, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sample        = factor(sample, levels = SAMPLE_ORDER),
    condition_label = SAMPLE_LABELS[as.character(sample)]
  ) %>%
  arrange(sample) %>%
  select(sample, condition_label, n_cells, mean_genes_det, median_genes_det)

# =============================================================
# 7. Merge with DEG genes-tested per contrast (if available)
# =============================================================

deg_summary_path <- file.path(out_dir, "DEG_summary.csv")

if (file.exists(deg_summary_path)) {
  deg_sum <- read.csv(deg_summary_path, stringsAsFactors = FALSE) %>%
    mutate(
      # extract the LCMV condition from "LCMV_Xwpi_vs_mock_6wpi"
      sample = sub("_vs_.*", "", contrast)
    ) %>%
    select(sample, n_genes_tested,
           n_up_FDR_0.05, n_down_FDR_0.05,
           n_up_FDR_0.05_logFC_0.25, n_down_FDR_0.05_logFC_0.25)

  qc_df <- qc_df %>%
    left_join(deg_sum, by = "sample")

  message("  DEG_summary.csv merged — n_genes_tested column added.")
} else {
  message("  DEG_summary.csv not found — run script 29 first to add DEG counts.")
}

# =============================================================
# 8. Print and save
# =============================================================

cat("\n=== Microglia (P2ry12) — QC / DEG summary ===\n")
print(as.data.frame(qc_df), row.names = FALSE, na.print = "-")

out_csv <- file.path(out_dir, "microglia_qc_table.csv")
write.csv(qc_df, out_csv, row.names = FALSE)
message("\nSaved: ", out_csv)
