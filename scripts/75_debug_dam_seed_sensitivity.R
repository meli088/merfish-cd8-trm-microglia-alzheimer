#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(readr)
})

se <- readRDS("objects/08_immune_annotated_lam02_res03.rds")
cell_type <- as.character(SummarizedExperiment::colData(se)$cell_type)
sample_vec <- as.character(SummarizedExperiment::colData(se)$sample)

use_label <- if ("Activated microglia (C1qa)" %in% unique(cell_type)) {
  "Activated microglia (C1qa)"
} else {
  "Microglia (C1qa)"
}

keep <- which(cell_type == use_label & sample_vec %in% c("LCMV_1wpi", "LCMV_6wpi"))
se_sub <- se[, keep]
assay_in <- if ("counts" %in% SummarizedExperiment::assayNames(se_sub)) "counts" else SummarizedExperiment::assayNames(se_sub)[1]

up_raw <- readr::read_csv("outputs/banksy/dam_signature_curation/Upregulated_DAM.csv", show_col_types = FALSE)$gene
down_raw <- readr::read_csv("outputs/banksy/dam_signature_curation/Downregulated_DAM.csv", show_col_types = FALSE)$gene

harm_intersect <- function(sig_genes, panel_genes) {
  panel_title <- tools::toTitleCase(tolower(panel_genes))
  panel_map <- setNames(panel_genes, panel_title)
  sig_title <- tools::toTitleCase(tolower(sig_genes))
  hits <- sig_title[sig_title %in% panel_title]
  matched <- panel_map[hits]
  unique(as.character(matched[!is.na(matched)]))
}

run_once <- function(seed_val, with_seed) {
  set.seed(seed_val)
  so <- suppressWarnings(as.Seurat(se_sub, counts = assay_in, data = NULL))
  da <- if ("RNA" %in% SeuratObject::Assays(so)) "RNA" else SeuratObject::Assays(so)[1]
  DefaultAssay(so) <- da
  so <- NormalizeData(so, assay = da, verbose = FALSE)
  so$sample <- as.character(SummarizedExperiment::colData(se_sub)$sample)

  up_genes <- harm_intersect(up_raw, rownames(so))
  down_genes <- harm_intersect(down_raw, rownames(so))

  if (with_seed) {
    so <- AddModuleScore(
      so,
      features = list(Upregulated_DAM = up_genes, Downregulated_DAM = down_genes),
      name = "",
      ctrl = 5,
      seed = 1997,
      verbose = FALSE
    )
  } else {
    so <- AddModuleScore(
      so,
      features = list(Upregulated_DAM = up_genes, Downregulated_DAM = down_genes),
      name = "",
      ctrl = 5,
      verbose = FALSE
    )
  }

  up_1 <- so@meta.data[so$sample == "LCMV_1wpi", "1"]
  up_6 <- so@meta.data[so$sample == "LCMV_6wpi", "1"]
  dn_1 <- so@meta.data[so$sample == "LCMV_1wpi", "2"]
  dn_6 <- so@meta.data[so$sample == "LCMV_6wpi", "2"]

  c(
    p_up = wilcox.test(up_6, up_1, exact = FALSE)$p.value,
    p_down = wilcox.test(dn_6, dn_1, exact = FALSE)$p.value,
    med_up_1 = median(up_1),
    med_up_6 = median(up_6),
    med_dn_1 = median(dn_1),
    med_dn_6 = median(dn_6)
  )
}

cat("=== Non-locked behavior (no seed in AddModuleScore) ===\n")
for (s in c(1, 2, 3, 42, 1997, 7777)) {
  r <- run_once(s, with_seed = FALSE)
  cat(sprintf("init_seed=%d | p_up=%.6g | p_down=%.6g | med_dn(1->6)=%.3f->%.3f\n",
              s, r[["p_up"]], r[["p_down"]], r[["med_dn_1"]], r[["med_dn_6"]]))
}

cat("\n=== Locked behavior (seed fixed inside AddModuleScore) ===\n")
for (s in c(1, 2, 3, 42, 1997, 7777)) {
  r <- run_once(s, with_seed = TRUE)
  cat(sprintf("init_seed=%d | p_up=%.6g | p_down=%.6g | med_dn(1->6)=%.3f->%.3f\n",
              s, r[["p_up"]], r[["p_down"]], r[["med_dn_1"]], r[["med_dn_6"]]))
}
