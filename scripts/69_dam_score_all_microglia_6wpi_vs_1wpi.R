#!/usr/bin/env Rscript
# ============================================================================
# Script 69 — DAM module score for all microglia (IN + OUT niche merged)
# Comparison: LCMV 6wpi vs LCMV 1wpi
#
# Inputs:
#   objects/08_immune_annotated_lam02_res03.rds     (Activated microglia C1qa)
#   objects/04_banksy_joint_lam08_after_bloc3.rds   (Microglia P2ry12)
#   ncells_by_sample_lam02_res09_joint_long.csv      (BANKSY annotations)
#   outputs/banksy/dam_signature_curation/Upregulated_DAM.csv
#   outputs/banksy/dam_signature_curation/Downregulated_DAM.csv
#
# Output:
#   outputs/banksy/immune_acod1/analysis/figures/
#     fig_dam_score_all_microglia_6wpi_vs_1wpi.pdf
#     fig_dam_score_all_microglia_6wpi_vs_1wpi.jpg
#
# ============================================================================

set.seed(1997)

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(Banksy)
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
    is.na(p)  ~ "n.s.",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "n.s."
  )
}

find_cl_col <- function(se, lam, res) {
  all_cols <- Banksy::clusterNames(se)
  lam_str  <- gsub("\\.", "\\\\.", as.character(lam))
  lam_cols <- all_cols[grep(paste0("lam", lam_str), all_cols)]
  if (length(lam_cols) == 0) return(NULL)
  res_vals <- suppressWarnings(as.numeric(sub(".*_res", "", lam_cols)))
  idx      <- which(!is.na(res_vals) & abs(res_vals - res) < 1e-8)
  if (length(idx) == 0) return(NULL)
  lam_cols[idx[1]]
}

message("Loading objects...")
se_immune <- readRDS("objects/08_immune_annotated_lam02_res03.rds")
se_global <- readRDS("objects/04_banksy_joint_lam08_after_bloc3.rds")

message("Reconstructing BANKSY annotations for global object...")
cl_col <- find_cl_col(se_global, lam = 0.2, res = 0.9)
if (is.null(cl_col)) {
  stop("BANKSY cluster column not found for lambda=0.2 res=0.9")
}

annot_long <- readr::read_delim(
  "ncells_by_sample_lam02_res09_joint_long.csv",
  delim = ";",
  locale = locale(decimal_mark = "."),
  show_col_types = FALSE,
  trim_ws = TRUE
) %>%
  dplyr::select(-matches("^Unnamed")) %>%
  dplyr::mutate(
    banksy_domain = as.character(banksy_domain),
    annotation = trimws(as.character(annotation))
  )

annotation_map <- annot_long %>%
  dplyr::filter(!is.na(annotation), annotation != "") %>%
  dplyr::distinct(banksy_domain, annotation)

domain_labels <- paste0("Domain_", as.character(colData(se_global)[[cl_col]]))
anno_lookup <- setNames(annotation_map$annotation, annotation_map$banksy_domain)
annotation_global <- ifelse(
  !is.na(anno_lookup[domain_labels]) & anno_lookup[domain_labels] != "",
  anno_lookup[domain_labels],
  "Non annote"
)

message("Selecting microglia populations and timepoints...")
cell_type_niche <- as.character(colData(se_immune)$cell_type)
in_candidates <- c("Activated microglia (C1qa)", "Microglia (C1qa)")
in_label <- in_candidates[in_candidates %in% unique(cell_type_niche)][1]
if (is.na(in_label)) {
  stop("In-niche C1qa label not found in object 08")
}

in_keep <- which(cell_type_niche == in_label & as.character(colData(se_immune)$sample) %in% TIMEPOINTS)
out_keep <- which(annotation_global == "Microglia (P2ry12)" & as.character(colData(se_global)$sample) %in% TIMEPOINTS)

if (length(in_keep) == 0 || length(out_keep) == 0) {
  stop("No cells found for requested microglia groups/timepoints")
}

se_in <- se_immune[, in_keep]
se_out <- se_global[, out_keep]

message(sprintf("In-niche cells kept: %d", ncol(se_in)))
message(sprintf("Out-niche cells kept: %d", ncol(se_out)))

message("Merging IN + OUT into All microglia...")
assay_in <- if ("counts" %in% assayNames(se_in)) "counts" else assayNames(se_in)[1]
assay_out <- if ("counts" %in% assayNames(se_out)) "counts" else assayNames(se_out)[1]

cnt_in <- assay(se_in, assay_in)
cnt_out <- assay(se_out, assay_out)
common_genes <- intersect(rownames(cnt_in), rownames(cnt_out))
if (length(common_genes) == 0) {
  stop("No common genes found between IN and OUT microglia objects")
}

cnt_all <- cbind(cnt_in[common_genes, , drop = FALSE], cnt_out[common_genes, , drop = FALSE])
meta_all <- dplyr::bind_rows(
  tibble(cell = colnames(cnt_in), sample = as.character(colData(se_in)$sample), niche = "In niche"),
  tibble(cell = colnames(cnt_out), sample = as.character(colData(se_out)$sample), niche = "Out niche")
)

so_all <- CreateSeuratObject(counts = cnt_all, project = "all_microglia")
so_all$sample <- meta_all$sample[match(colnames(so_all), meta_all$cell)]
so_all$niche_status <- meta_all$niche[match(colnames(so_all), meta_all$cell)]

so_all <- NormalizeData(so_all, verbose = FALSE)

message("Loading and harmonizing DAM signatures...")
up_raw <- readr::read_csv("outputs/banksy/dam_signature_curation/Upregulated_DAM.csv", show_col_types = FALSE)$gene
down_raw <- readr::read_csv("outputs/banksy/dam_signature_curation/Downregulated_DAM.csv", show_col_types = FALSE)$gene

panel_genes <- rownames(so_all)
panel_title <- tools::toTitleCase(tolower(panel_genes))
panel_map <- setNames(panel_genes, panel_title)

harm_intersect <- function(sig_genes) {
  sig_title <- tools::toTitleCase(tolower(sig_genes))
  hits <- sig_title[sig_title %in% panel_title]
  matched <- panel_map[hits]
  unique(as.character(matched[!is.na(matched)]))
}

up_genes <- harm_intersect(up_raw)
down_genes <- harm_intersect(down_raw)

if (length(up_genes) == 0 || length(down_genes) == 0) {
  stop("DAM signature genes could not be mapped to expression panel")
}

message(sprintf("Upregulated DAM genes used: %d", length(up_genes)))
message(sprintf("Downregulated DAM genes used: %d", length(down_genes)))

message("Computing AddModuleScore...")
so_all <- AddModuleScore(
  so_all,
  features = list(Upregulated_DAM = up_genes, Downregulated_DAM = down_genes),
  name = "",
  ctrl = 5,
  verbose = FALSE
)

score_df <- tibble(
  cell = colnames(so_all),
  sample = as.character(so_all$sample),
  niche_status = as.character(so_all$niche_status),
  Upregulated_DAM = so_all@meta.data[, "1"],
  Downregulated_DAM = so_all@meta.data[, "2"]
) %>%
  dplyr::filter(sample %in% TIMEPOINTS) %>%
  dplyr::mutate(
    timepoint = factor(TIME_LABELS[sample], levels = c("LCMV 1 wpi", "LCMV 6 wpi"))
  )

n_by_time <- table(score_df$timepoint)
message(sprintf("All microglia counts: 1wpi=%d, 6wpi=%d", n_by_time[["LCMV 1 wpi"]], n_by_time[["LCMV 6 wpi"]]))

p_up <- tryCatch(
  wilcox.test(
    score_df$Upregulated_DAM[score_df$timepoint == "LCMV 6 wpi"],
    score_df$Upregulated_DAM[score_df$timepoint == "LCMV 1 wpi"],
    exact = FALSE
  )$p.value,
  error = function(e) NA_real_
)

p_down <- tryCatch(
  wilcox.test(
    score_df$Downregulated_DAM[score_df$timepoint == "LCMV 6 wpi"],
    score_df$Downregulated_DAM[score_df$timepoint == "LCMV 1 wpi"],
    exact = FALSE
  )$p.value,
  error = function(e) NA_real_
)

stars_up <- pval_stars(p_up)
stars_down <- pval_stars(p_down)

make_panel <- function(df, score_col, pval, stars, title_txt) {
  y_max <- max(df[[score_col]], na.rm = TRUE)
  y_min <- min(df[[score_col]], na.rm = TRUE)
  y_rng <- y_max - y_min
  y_txt <- y_max + 0.08 * y_rng

  ggplot(df, aes(x = timepoint, y = .data[[score_col]], fill = timepoint)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.6, color = NA) +
    geom_boxplot(width = 0.2, outlier.size = 0.3, outlier.alpha = 0.25,
                 color = "grey20", fill = "white") +
    annotate("text", x = 1.5, y = y_txt, label = stars, size = 5, fontface = "bold") +
    scale_fill_manual(values = TIME_COLORS, guide = "none") +
    labs(
      title = title_txt,
      subtitle = sprintf("Wilcoxon p=%s", ifelse(is.na(pval), "NA", formatC(pval, digits = 3, format = "e"))),
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

p1 <- make_panel(score_df, "Upregulated_DAM", p_up, stars_up, "Upregulated DAM")
p2 <- make_panel(score_df, "Downregulated_DAM", p_down, stars_down, "Downregulated DAM")

p_final <- p1 + p2 +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "All microglia (C1qa + P2ry12): DAM score 6wpi vs 1wpi",
    subtitle = sprintf("n(1wpi)=%d, n(6wpi)=%d", n_by_time[["LCMV 1 wpi"]], n_by_time[["LCMV 6 wpi"]]),
    theme = theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 9, color = "grey35", hjust = 0.5)
    )
  )

pdf_out <- file.path(OUT_DIR, "fig_dam_score_all_microglia_6wpi_vs_1wpi.pdf")
jpg_out <- file.path(OUT_DIR, "fig_dam_score_all_microglia_6wpi_vs_1wpi.jpg")

ggsave(pdf_out, p_final, width = 10.5, height = 4.8, device = cairo_pdf)
ggsave(jpg_out, p_final, width = 10.5, height = 4.8, dpi = 300)

stats_out <- tibble(
  comparison = c("Upregulated DAM", "Downregulated DAM"),
  n_1wpi = c(n_by_time[["LCMV 1 wpi"]], n_by_time[["LCMV 1 wpi"]]),
  n_6wpi = c(n_by_time[["LCMV 6 wpi"]], n_by_time[["LCMV 6 wpi"]]),
  median_1wpi = c(
    median(score_df$Upregulated_DAM[score_df$timepoint == "LCMV 1 wpi"], na.rm = TRUE),
    median(score_df$Downregulated_DAM[score_df$timepoint == "LCMV 1 wpi"], na.rm = TRUE)
  ),
  median_6wpi = c(
    median(score_df$Upregulated_DAM[score_df$timepoint == "LCMV 6 wpi"], na.rm = TRUE),
    median(score_df$Downregulated_DAM[score_df$timepoint == "LCMV 6 wpi"], na.rm = TRUE)
  ),
  p_value_wilcoxon = c(p_up, p_down),
  stars = c(stars_up, stars_down)
)

readr::write_csv(stats_out, file.path(OUT_DIR, "dam_score_all_microglia_6wpi_vs_1wpi_stats.csv"))

message("Saved figure:")
message(pdf_out)
message(jpg_out)
message("Saved stats:")
message(file.path(OUT_DIR, "dam_score_all_microglia_6wpi_vs_1wpi_stats.csv"))
