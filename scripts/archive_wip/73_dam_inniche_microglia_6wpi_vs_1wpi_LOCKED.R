#!/usr/bin/env Rscript
# ============================================================================
# Script 73 (LOCKED) — In-niche microglia (C1qa) DAM score: 6wpi vs 1wpi
#
# This script is the single authoritative run for manuscript values.
# Fixed choices:
#   - population: in-niche C1qa only (object 08)
#   - timepoints: LCMV_1wpi vs LCMV_6wpi
#   - signatures: Upregulated_DAM and Downregulated_DAM
#   - method: NormalizeData + AddModuleScore(ctrl=5, seed=1997)
#   - test: Wilcoxon (exact=FALSE)
#   - seed: 1997
#
# Outputs:
#   outputs/banksy/immune_acod1/analysis/figures/
#     fig_dam_score_inniche_microglia_6wpi_vs_1wpi_LOCKED.pdf/.jpg
#     dam_score_inniche_microglia_6wpi_vs_1wpi_LOCKED_stats.csv
# ============================================================================

set.seed(1997)

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

source("scripts/00_palette.R")

OUT_DIR <- "outputs/banksy/immune_acod1/analysis/figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TIMEPOINTS <- c("LCMV_1wpi", "LCMV_6wpi")
TIME_LABELS <- c("LCMV_1wpi" = "LCMV 1 wpi", "LCMV_6wpi" = "LCMV 6 wpi")
TIME_COLORS <- c("LCMV 1 wpi" = "#56B4E9", "LCMV 6 wpi" = "#F28E2B")

pval_stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "n.s.",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "n.s."
  )
}

harm_intersect <- function(sig_genes, panel_genes) {
  panel_title <- tools::toTitleCase(tolower(panel_genes))
  panel_map <- setNames(panel_genes, panel_title)
  sig_title <- tools::toTitleCase(tolower(sig_genes))
  hits <- sig_title[sig_title %in% panel_title]
  matched <- panel_map[hits]
  unique(as.character(matched[!is.na(matched)]))
}

message("Loading in-niche object...")
se <- readRDS("objects/08_immune_annotated_lam02_res03.rds")
cell_type <- as.character(SummarizedExperiment::colData(se)$cell_type)
sample_vec <- as.character(SummarizedExperiment::colData(se)$sample)

# Locked label priority
if ("Activated microglia (C1qa)" %in% unique(cell_type)) {
  use_label <- "Activated microglia (C1qa)"
} else if ("Microglia (C1qa)" %in% unique(cell_type)) {
  use_label <- "Microglia (C1qa)"
} else {
  stop("No in-niche C1qa label found in object 08")
}

keep <- which(cell_type == use_label & sample_vec %in% TIMEPOINTS)
if (length(keep) == 0) stop("No cells for selected label and timepoints")

se_sub <- se[, keep]
assay_in <- if ("counts" %in% SummarizedExperiment::assayNames(se_sub)) "counts" else SummarizedExperiment::assayNames(se_sub)[1]
so <- suppressWarnings(as.Seurat(se_sub, counts = assay_in, data = NULL))
dassay <- if ("RNA" %in% SeuratObject::Assays(so)) "RNA" else SeuratObject::Assays(so)[1]
DefaultAssay(so) <- dassay
so <- NormalizeData(so, assay = dassay, verbose = FALSE)

so$sample <- as.character(SummarizedExperiment::colData(se_sub)$sample)

message("Loading DAM signatures...")
up_raw <- readr::read_csv("outputs/banksy/dam_signature_curation/Upregulated_DAM.csv", show_col_types = FALSE)$gene
down_raw <- readr::read_csv("outputs/banksy/dam_signature_curation/Downregulated_DAM.csv", show_col_types = FALSE)$gene

up_genes <- harm_intersect(up_raw, rownames(so))
down_genes <- harm_intersect(down_raw, rownames(so))

if (length(up_genes) == 0 || length(down_genes) == 0) {
  stop("Signature mapping failed")
}

message(sprintf("Using label: %s", use_label))
message(sprintf("Cells: %d", ncol(so)))
message(sprintf("Up genes mapped: %d", length(up_genes)))
message(sprintf("Down genes mapped: %d", length(down_genes)))

so <- AddModuleScore(
  so,
  features = list(Upregulated_DAM = up_genes, Downregulated_DAM = down_genes),
  name = "",
  ctrl = 5,
  seed = 1997,
  verbose = FALSE
)

score_df <- tibble(
  sample = as.character(so$sample),
  Upregulated_DAM = so@meta.data[, "1"],
  Downregulated_DAM = so@meta.data[, "2"]
) %>%
  mutate(timepoint = factor(TIME_LABELS[sample], levels = c("LCMV 1 wpi", "LCMV 6 wpi")))

p_up <- wilcox.test(
  score_df$Upregulated_DAM[score_df$timepoint == "LCMV 6 wpi"],
  score_df$Upregulated_DAM[score_df$timepoint == "LCMV 1 wpi"],
  exact = FALSE
)$p.value

p_down <- wilcox.test(
  score_df$Downregulated_DAM[score_df$timepoint == "LCMV 6 wpi"],
  score_df$Downregulated_DAM[score_df$timepoint == "LCMV 1 wpi"],
  exact = FALSE
)$p.value

make_panel <- function(df, col_name, title_txt, pval) {
  y_max <- max(df[[col_name]], na.rm = TRUE)
  y_min <- min(df[[col_name]], na.rm = TRUE)
  y_txt <- y_max + 0.08 * (y_max - y_min)

  ggplot(df, aes(x = timepoint, y = .data[[col_name]], fill = timepoint)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.6, color = NA) +
    geom_boxplot(width = 0.2, outlier.size = 0.3, outlier.alpha = 0.25,
                 color = "grey20", fill = "white") +
    annotate("text", x = 1.5, y = y_txt, label = pval_stars(pval), size = 5, fontface = "bold") +
    scale_fill_manual(values = TIME_COLORS, guide = "none") +
    labs(
      title = title_txt,
      subtitle = sprintf("Wilcoxon p=%s", formatC(pval, digits = 3, format = "e")),
      x = NULL,
      y = "DAM module score"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8, color = "grey40"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

p1 <- make_panel(score_df, "Upregulated_DAM", "Upregulated DAM", p_up)
p2 <- make_panel(score_df, "Downregulated_DAM", "Downregulated DAM", p_down)

n_by_time <- table(score_df$timepoint)

p_final <- p1 + p2 +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "In-niche microglia (C1qa): DAM score 6wpi vs 1wpi [LOCKED]",
    subtitle = sprintf("n(1wpi)=%d, n(6wpi)=%d", n_by_time[["LCMV 1 wpi"]], n_by_time[["LCMV 6 wpi"]]),
    theme = theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 9, color = "grey35", hjust = 0.5)
    )
  )

pdf_out <- file.path(OUT_DIR, "fig_dam_score_inniche_microglia_6wpi_vs_1wpi_LOCKED.pdf")
jpg_out <- file.path(OUT_DIR, "fig_dam_score_inniche_microglia_6wpi_vs_1wpi_LOCKED.jpg")

ggsave(pdf_out, p_final, width = 10.5, height = 4.8, device = cairo_pdf)
ggsave(jpg_out, p_final, width = 10.5, height = 4.8, dpi = 300)

stats <- tibble(
  comparison = c("Upregulated DAM", "Downregulated DAM"),
  population = use_label,
  n_1wpi = c(n_by_time[["LCMV 1 wpi"]], n_by_time[["LCMV 1 wpi"]]),
  n_6wpi = c(n_by_time[["LCMV 6 wpi"]], n_by_time[["LCMV 6 wpi"]]),
  genes_mapped = c(length(up_genes), length(down_genes)),
  median_1wpi = c(
    median(score_df$Upregulated_DAM[score_df$timepoint == "LCMV 1 wpi"], na.rm = TRUE),
    median(score_df$Downregulated_DAM[score_df$timepoint == "LCMV 1 wpi"], na.rm = TRUE)
  ),
  median_6wpi = c(
    median(score_df$Upregulated_DAM[score_df$timepoint == "LCMV 6 wpi"], na.rm = TRUE),
    median(score_df$Downregulated_DAM[score_df$timepoint == "LCMV 6 wpi"], na.rm = TRUE)
  ),
  delta_median_6wpi_minus_1wpi = median_6wpi - median_1wpi,
  p_value_wilcoxon = c(p_up, p_down),
  stars = c(pval_stars(p_up), pval_stars(p_down)),
  seed = 1997,
  ctrl = 5
)

stats_out <- file.path(OUT_DIR, "dam_score_inniche_microglia_6wpi_vs_1wpi_LOCKED_stats.csv")
readr::write_csv(stats, stats_out)

message("Saved figure:")
message(pdf_out)
message(jpg_out)
message("Saved stats:")
message(stats_out)
print(stats)
