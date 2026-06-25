#!/usr/bin/env Rscript

set.seed(1997)

suppressPackageStartupMessages({
  library(Seurat)
  library(SpatialExperiment)
  library(SummarizedExperiment)
})

OBJ_FILE <- file.path("objects", "08_immune_annotated_lam02_res03.rds")
UP_SIG_FILE <- file.path("outputs", "banksy", "dam_signature_curation", "Upregulated_DAM.csv")
DOWN_SIG_FILE <- file.path("outputs", "banksy", "dam_signature_curation", "Downregulated_DAM.csv")
OUT_FILE <- file.path("outputs", "banksy", "immune_acod1", "analysis", "figures", "dam_module_score_stats_6wpi_vs_1wpi.csv")

TIME_REF <- "LCMV_1wpi"
TIME_CASE <- "LCMV_6wpi"

read_sig <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE)
  unique(as.character(x$gene))
}

match_signature_genes <- function(sig_genes, panel_genes) {
  panel_title <- tools::toTitleCase(tolower(panel_genes))
  map_panel <- setNames(panel_genes, panel_title)
  sig_title <- tools::toTitleCase(tolower(sig_genes))
  hits <- sig_title[sig_title %in% panel_title]
  out <- unique(as.character(map_panel[hits]))
  out[!is.na(out)]
}

calc_one <- function(se, pop_label, sig_genes, sig_label) {
  cd <- as.data.frame(colData(se))
  idx <- which(as.character(cd$cell_type) == pop_label & as.character(cd$sample) %in% c(TIME_REF, TIME_CASE))
  se_sub <- se[, idx]

  assay_use <- if ("counts" %in% assayNames(se_sub)) "counts" else assayNames(se_sub)[1]
  so <- suppressWarnings(as.Seurat(se_sub, counts = assay_use, data = NULL))
  da <- if ("RNA" %in% SeuratObject::Assays(so)) "RNA" else SeuratObject::Assays(so)[1]
  DefaultAssay(so) <- da
  so <- NormalizeData(so, assay = da, verbose = FALSE)
  so$sample <- as.character(colData(se_sub)$sample)

  sig_hits <- match_signature_genes(sig_genes, rownames(so))
  ctrl_n <- max(1L, min(5L, floor((nrow(so) - length(sig_hits)) / 2)))

  score_name <- paste0(gsub("[^A-Za-z0-9]", "_", sig_label), "_score")
  so <- AddModuleScore(so, features = list(sig_hits), name = score_name, ctrl = ctrl_n, seed = 1997)
  score_col <- paste0(score_name, "1")

  d1 <- so[[score_col, drop = TRUE]][so$sample == TIME_REF]
  d2 <- so[[score_col, drop = TRUE]][so$sample == TIME_CASE]

  pval <- tryCatch(wilcox.test(d2, d1, exact = FALSE)$p.value, error = function(e) NA_real_)

  data.frame(
    population = pop_label,
    signature = sig_label,
    n_1wpi = length(d1),
    n_6wpi = length(d2),
    median_1wpi = median(d1, na.rm = TRUE),
    median_6wpi = median(d2, na.rm = TRUE),
    delta_median_6wpi_minus_1wpi = median(d2, na.rm = TRUE) - median(d1, na.rm = TRUE),
    p_value_wilcoxon = pval,
    stringsAsFactors = FALSE
  )
}

se <- readRDS(OBJ_FILE)
up <- read_sig(UP_SIG_FILE)
down <- read_sig(DOWN_SIG_FILE)

res <- rbind(
  calc_one(se, "Activated microglia (C1qa)", up, "Upregulated DAM"),
  calc_one(se, "Activated microglia (C1qa)", down, "Downregulated DAM"),
  calc_one(se, "Antigen-presenting myeloid cells (Cd74)", up, "Upregulated DAM"),
  calc_one(se, "Antigen-presenting myeloid cells (Cd74)", down, "Downregulated DAM")
)

write.csv(res, OUT_FILE, row.names = FALSE)
cat("Saved:", OUT_FILE, "\n")
print(res)
