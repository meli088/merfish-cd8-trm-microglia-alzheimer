#!/usr/bin/env Rscript

set.seed(1997)

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(dplyr)
  library(ggplot2)
  library(Cairo)
})

base_path <- normalizePath(".")
setwd(base_path)

OBJ_IMMUNE <- file.path("objects", "08_immune_annotated_lam02_res03.rds")
SIG_FILE <- file.path("outputs", "banksy", "dam_signature_curation", "Upregulated_DAM.csv")
OUT_DIR <- file.path("outputs", "banksy", "microglia_dam_niche")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TARGET_LABEL <- "Activated microglia (C1qa)"
TP_REF <- "LCMV_1wpi"
TP_CASE <- "LCMV_6wpi"

save_fig <- function(p, out_base, width = 6.2, height = 4.8) {
  CairoPDF(paste0(out_base, ".pdf"), width = width, height = height)
  print(p)
  dev.off()
  CairoJPEG(paste0(out_base, ".jpg"), width = width * 150, height = height * 150, res = 150)
  print(p)
  dev.off()
}

read_sig <- function(path) {
  stopifnot(file.exists(path))
  x <- read.csv(path, stringsAsFactors = FALSE)
  if (!"gene" %in% colnames(x)) stop("Missing column 'gene' in signature file")
  unique(as.character(x$gene))
}

match_signature_genes <- function(sig_genes, panel_genes) {
  panel_title <- tools::toTitleCase(tolower(panel_genes))
  map_panel <- setNames(panel_genes, panel_title)
  sig_title <- tools::toTitleCase(tolower(sig_genes))
  sig_hits <- sig_title[sig_title %in% panel_title]
  out <- unique(as.character(map_panel[sig_hits]))
  out[!is.na(out)]
}

message("Loading: ", OBJ_IMMUNE)
se <- readRDS(OBJ_IMMUNE)
cd <- as.data.frame(colData(se))
stopifnot("cell_type missing" = "cell_type" %in% colnames(cd))
stopifnot("sample missing" = "sample" %in% colnames(cd))

idx <- which(as.character(cd$cell_type) == TARGET_LABEL & as.character(cd$sample) %in% c(TP_REF, TP_CASE))
if (length(idx) == 0) stop("No cells found for target label/timepoints")
se_sub <- se[, idx]

assay_use <- if ("counts" %in% assayNames(se_sub)) "counts" else assayNames(se_sub)[1]
so <- suppressWarnings(as.Seurat(se_sub, counts = assay_use, data = NULL))
da <- if ("RNA" %in% SeuratObject::Assays(so)) "RNA" else SeuratObject::Assays(so)[1]
DefaultAssay(so) <- da
so <- NormalizeData(so, assay = da, verbose = FALSE)

so$sample <- as.character(colData(se_sub)$sample)

sig_raw <- read_sig(SIG_FILE)
sig_hits <- match_signature_genes(sig_raw, rownames(so))
if (length(sig_hits) == 0) stop("No signature genes matched panel")

ctrl_n <- max(1L, min(5L, floor((nrow(so) - length(sig_hits)) / 2)))
so <- AddModuleScore(
  so,
  features = list(sig_hits),
  name = "Upregulated_DAM",
  ctrl = ctrl_n,
  seed = 1997
)

plot_df <- data.frame(
  sample = as.character(so$sample),
  dam_score = so$Upregulated_DAM1,
  stringsAsFactors = FALSE
) %>%
  filter(sample %in% c(TP_REF, TP_CASE)) %>%
  mutate(sample = factor(sample, levels = c(TP_REF, TP_CASE), labels = c("LCMV 1 wpi", "LCMV 6 wpi")))

pval <- tryCatch(
  wilcox.test(
    plot_df$dam_score[plot_df$sample == "LCMV 6 wpi"],
    plot_df$dam_score[plot_df$sample == "LCMV 1 wpi"],
    exact = FALSE
  )$p.value,
  error = function(e) NA_real_
)

p <- ggplot(plot_df, aes(x = sample, y = dam_score, fill = sample)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.58, colour = NA) +
  geom_boxplot(width = 0.18, outlier.size = 0.35, outlier.alpha = 0.35, colour = "grey20", fill = "white") +
  scale_fill_manual(values = c("LCMV 1 wpi" = "#56B4E9", "LCMV 6 wpi" = "#F28E2B"), guide = "none") +
  labs(
    title = "DAM signature score in Activated microglia (C1qa): 6 wpi vs 1 wpi",
    subtitle = paste0("Wilcoxon p = ", ifelse(is.na(pval), "NA", formatC(pval, digits = 3, format = "e"))),
    x = NULL,
    y = "DAM module score (AddModuleScore)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, colour = "grey40"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

out_base <- file.path(OUT_DIR, "fig_dam_signature_activated_microglia_c1qa_6wpi_vs_1wpi_violin")
save_fig(p, out_base)

write.csv(
  plot_df,
  file.path(OUT_DIR, "dam_signature_activated_microglia_c1qa_6wpi_vs_1wpi_scores.csv"),
  row.names = FALSE
)

stats_df <- data.frame(
  comparison = "LCMV_6wpi_vs_LCMV_1wpi",
  p_value_wilcox = pval,
  n_1wpi = sum(plot_df$sample == "LCMV 1 wpi"),
  n_6wpi = sum(plot_df$sample == "LCMV 6 wpi"),
  stringsAsFactors = FALSE
)
write.csv(
  stats_df,
  file.path(OUT_DIR, "dam_signature_activated_microglia_c1qa_6wpi_vs_1wpi_stats.csv"),
  row.names = FALSE
)

message("Saved DAM outputs in: ", OUT_DIR)
