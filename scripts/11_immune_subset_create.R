#!/usr/bin/env Rscript

library(SingleCellExperiment)
library(SummarizedExperiment)
library(SpatialExperiment)
library(tidyverse)

base_path <- normalizePath(".")  # Run this script from the project root directory
setwd(base_path)

OBJ_PATH <- file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds")
CSV_ANNOT <- "ncells_by_sample_lam02_res09_joint_long.csv"

LAM <- 0.2
RES <- 0.9
# Default label; can be overridden with --label "Immune (Acod1)"
IMMUNE_LABEL <- "Immune (Acod1)"
FOLDER_NAME  <- NULL   # optional: override folder/file slug with --folder microglia

# Allow CLI override: Rscript scripts/11_immune_subset_create.R --label "Immune (Acod1)" --folder immune
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
subset_dir <- file.path(label_root, "subset")
if (!dir.exists(subset_dir)) dir.create(subset_dir, recursive = TRUE, showWarnings = FALSE)

OUT_FILE <- file.path("objects", paste0("05_immune_subset_lam02_res09_", safe_label, ".rds"))
OUT_SUM <- file.path(subset_dir, "immune_subset_summary.csv")

find_cl_col <- function(se, lam, res) {
  all_cols <- colnames(as.data.frame(SummarizedExperiment::colData(se)))
  lam_pat <- paste0("lam", gsub("\\.", "\\\\.", as.character(lam)))
  res_pat <- paste0("res", gsub("\\.", "\\\\.", as.character(res)))
  candidates <- all_cols[grepl(lam_pat, all_cols) & grepl(res_pat, all_cols)]
  if (length(candidates) == 0) {
    stop("No BANKSY clustering column found for lam=", lam, " and res=", res)
  }
  candidates[1]
}

cat("\n=== Create immune subset from BANKSY joint object ===\n")

se <- readRDS(OBJ_PATH)
cd <- as.data.frame(SummarizedExperiment::colData(se))

cluster_col <- find_cl_col(se, LAM, RES)
cat("Using clustering column:", cluster_col, "\n")

# Verify joint nature: report unique samples covered by this joint object
unique_samples <- sort(unique(as.character(cd$sample)))
cat("Joint object contains", length(unique_samples), "samples:\n")
print(unique_samples)

anno_data <- read.delim(CSV_ANNOT, sep = ";", stringsAsFactors = FALSE)

if (!all(c("banksy_domain", "annotation") %in% colnames(anno_data))) {
  stop("CSV must contain columns: banksy_domain and annotation")
}

anno_map <- anno_data %>%
  filter(annotation != "" & !is.na(annotation)) %>%
  mutate(annotation = trimws(annotation)) %>%
  select(banksy_domain, annotation) %>%
  distinct()

dup_check <- anno_map %>%
  count(banksy_domain) %>%
  filter(n > 1)

if (nrow(dup_check) > 0) {
  stop("Ambiguous annotation mapping: some banksy_domain values have multiple annotations")
}

anno_map$cluster_id <- as.numeric(gsub("Domain_", "", anno_map$banksy_domain))
anno_lookup <- setNames(anno_map$annotation, anno_map$cluster_id)


banksy_clusters <- as.numeric(cd[[cluster_col]])
cell_annotations <- anno_lookup[as.character(banksy_clusters)]
cell_annotations <- trimws(cell_annotations)
cell_annotations[is.na(cell_annotations)] <- "Unknown"

available <- sort(unique(cell_annotations))
print(available)

# match label case-insensitively after trimming
norm_annos <- tolower(trimws(cell_annotations))
norm_label <- tolower(trimws(IMMUNE_LABEL))
if (!(norm_label %in% norm_annos)) {
  # helpful message with available annotations
  stop("IMMUNE_LABEL not found in annotations: ", IMMUNE_LABEL, "\nAvailable: ", paste(available, collapse = ", "))
}

keep_idx <- which(norm_annos == norm_label)
if (length(keep_idx) == 0) {
  stop("No cells found for immune label: ", IMMUNE_LABEL)
}

se_immune <- se[, keep_idx]
colData(se_immune)$annotation <- cell_annotations[keep_idx]
colData(se_immune)$banksy_domain <- paste0("Domain_", banksy_clusters[keep_idx])

saveRDS(se_immune, OUT_FILE)

# write per-sample summary CSV
summary_tbl <- as.data.frame(table(sample = as.character(colData(se_immune)$sample)), stringsAsFactors = FALSE)
colnames(summary_tbl) <- c("sample", "n_cells")
summary_tbl <- summary_tbl %>% arrange(match(sample, c("mock_6wpi", "LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")))
if (!dir.exists(dirname(OUT_SUM))) dir.create(dirname(OUT_SUM), recursive = TRUE, showWarnings = FALSE)
write.csv(summary_tbl, OUT_SUM, row.names = FALSE)

cat("Immune label used:", IMMUNE_LABEL, "\n")
cat("Cells kept:", ncol(se_immune), "\n")
cat("Per-sample counts:\n")
print(summary_tbl)
cat("Saved immune subset to:", OUT_FILE, "\n")
cat("Saved per-sample summary to:", OUT_SUM, "\n")