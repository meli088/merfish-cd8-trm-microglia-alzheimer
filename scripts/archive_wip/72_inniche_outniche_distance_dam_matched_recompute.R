#!/usr/bin/env Rscript
# ============================================================================
# Script 72 — Recomputed matched IN/OUT figures
# Distance violin + DAM barplot by distance bins (same style for both)
#
# Outputs:
#   outputs/banksy/microglia_dam_niche/
#     fig_combined_inniche_distance_dam_matched.pdf/.jpg
#     fig_combined_outniche_distance_dam_matched.pdf/.jpg
#     summary_inniche_distance_dam_matched.csv
#     summary_outniche_distance_dam_matched.csv
# ============================================================================

set.seed(1997)

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(Banksy)
  library(FNN)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(Cairo)
})

source("scripts/00_palette.R")

out_dir <- "outputs/banksy/microglia_dam_niche"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

SAMPLE_ORDER <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi", "mock_6wpi")
SAMPLE_LABELS <- c(
  mock_6wpi = "Mock 6 wpi",
  LCMV_1wpi = "LCMV 1 wpi",
  LCMV_3wpi = "LCMV 3 wpi",
  LCMV_6wpi = "LCMV 6 wpi"
)
TIMEPOINTS <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")

# Use the same displayed bins as the existing figure panel
DIST_BINS <- c(0, 50, 100, 200)
DIST_LABELS <- c("<50um", "50-100um", "100-200um")

cond_palette <- c(
  "Mock 6 wpi" = "grey60",
  "LCMV 1 wpi" = "#56B4E9",
  "LCMV 3 wpi" = "#E69F00",
  "LCMV 6 wpi" = "#D55E00"
)

pval_stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "n.s.",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "n.s."
  )
}

save_fig <- function(p, fname, width, height) {
  CairoPDF(file.path(out_dir, paste0(fname, ".pdf")), width = width, height = height)
  print(p)
  dev.off()

  CairoJPEG(file.path(out_dir, paste0(fname, ".jpg")),
            width = width * 150, height = height * 150, res = 150)
  print(p)
  dev.off()
}

find_cl_col <- function(se, lam, res) {
  all_cols <- Banksy::clusterNames(se)
  lam_str <- gsub("\\.", "\\\\.", as.character(lam))
  lam_cols <- all_cols[grep(paste0("lam", lam_str), all_cols)]
  if (length(lam_cols) == 0) return(NULL)
  res_vals <- suppressWarnings(as.numeric(sub(".*_res", "", lam_cols)))
  idx <- which(!is.na(res_vals) & abs(res_vals - res) < 1e-8)
  if (length(idx) == 0) return(NULL)
  lam_cols[idx[1]]
}

message("Loading objects and annotations...")
se_global <- readRDS("objects/04_banksy_joint_lam08_after_bloc3.rds")
se_niche <- readRDS("objects/08_immune_annotated_lam02_res03.rds")

cl_col <- find_cl_col(se_global, lam = 0.2, res = 0.9)
if (is.null(cl_col)) stop("Could not find BANKSY cluster column")

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

domain_labels <- paste0("Domain_", as.character(SummarizedExperiment::colData(se_global)[[cl_col]]))
anno_lookup <- setNames(annotation_map$annotation, annotation_map$banksy_domain)
annotation_global <- ifelse(
  !is.na(anno_lookup[domain_labels]) & anno_lookup[domain_labels] != "",
  anno_lookup[domain_labels],
  "Non annote"
)

# Out-niche mask from global annotation
out_mask <- annotation_global == "Microglia (P2ry12)"

# In-niche mask: prefer explicit global annotation, fallback to IDs from niche object
in_label_candidates <- c("Activated microglia (C1qa)", "Microglia (C1qa)")
niche_ct <- as.character(SummarizedExperiment::colData(se_niche)$cell_type)
in_label <- in_label_candidates[in_label_candidates %in% unique(niche_ct)][1]
if (is.na(in_label)) stop("In-niche C1qa label not found in object 08")
in_ids <- colnames(se_niche)[niche_ct == in_label]

in_mask_global <- annotation_global == "Microglia (C1qa)"
if (!any(in_mask_global)) {
  in_mask_global <- colnames(se_global) %in% in_ids
}

# T-cell / immune mask
tcell_patterns <- c("Immune \\(Acod1\\)", "T cells")
tcell_mask <- grepl(paste(tcell_patterns, collapse = "|"), annotation_global)

samples_global <- as.character(SummarizedExperiment::colData(se_global)$sample)
coords_global <- as.data.frame(SpatialExperiment::spatialCoords(se_global))
colnames(coords_global) <- c("x", "y")

compute_dist <- function(sample_id, mg_mask, tc_mask, niche_status) {
  mg_idx <- which(mg_mask & samples_global == sample_id)
  tc_idx <- which(tc_mask & samples_global == sample_id)

  if (length(mg_idx) == 0) return(NULL)
  if (length(tc_idx) == 0) {
    return(data.frame(
      sample = sample_id,
      cell_id = colnames(se_global)[mg_idx],
      dist_um = NA_real_,
      niche_status = niche_status,
      stringsAsFactors = FALSE
    ))
  }

  mg_coords <- as.matrix(coords_global[mg_idx, c("x", "y")])
  tc_coords <- as.matrix(coords_global[tc_idx, c("x", "y")])
  knn_res <- FNN::get.knnx(data = tc_coords, query = mg_coords, k = 1)

  data.frame(
    sample = sample_id,
    cell_id = colnames(se_global)[mg_idx],
    dist_um = as.numeric(knn_res$nn.dist[, 1]),
    niche_status = niche_status,
    stringsAsFactors = FALSE
  )
}

message("Computing OUT distances...")
dist_out <- lapply(SAMPLE_ORDER, function(s) {
  compute_dist(s, out_mask, tcell_mask, "Out niche")
}) %>% bind_rows()

message("Computing IN distances (exclude IN microglia from target immune pool)...")
tcell_mask_in <- tcell_mask & !in_mask_global

dist_in <- lapply(SAMPLE_ORDER, function(s) {
  compute_dist(s, in_mask_global, tcell_mask_in, "In niche")
}) %>% bind_rows()

# Read existing module scores table from script 39
score_df <- readr::read_csv(file.path(out_dir, "module_scores_per_cell.csv"), show_col_types = FALSE) %>%
  mutate(cell_id_orig = sub("^(out|in)_", "", cell_id))

build_figure <- function(dist_sub, niche_status_value, mg_label, out_stub, title_prefix) {
  dist_sub <- dist_sub %>%
    filter(!is.na(dist_um)) %>%
    mutate(sample_label = factor(SAMPLE_LABELS[sample], levels = SAMPLE_LABELS[SAMPLE_ORDER]))

  # Violin stats vs mock
  wilcox_df <- lapply(TIMEPOINTS, function(s) {
    d1 <- dist_sub$dist_um[dist_sub$sample == s]
    d2 <- dist_sub$dist_um[dist_sub$sample == "mock_6wpi"]
    pv <- tryCatch(wilcox.test(d1, d2, exact = FALSE)$p.value, error = function(e) NA_real_)
    data.frame(
      sample_label = factor(SAMPLE_LABELS[s], levels = SAMPLE_LABELS[SAMPLE_ORDER]),
      stars = pval_stars(pv),
      pval = pv,
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()

  y_max <- max(dist_sub$dist_um, na.rm = TRUE)
  y_rng <- diff(range(dist_sub$dist_um, na.rm = TRUE))
  y_txt <- y_max + 0.05 * y_rng

  p_violin <- ggplot(dist_sub, aes(x = sample_label, y = dist_um, fill = sample_label)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.55, colour = NA) +
    geom_boxplot(width = 0.18, outlier.size = 0.3, outlier.alpha = 0.3,
                 colour = "grey20", fill = "white") +
    geom_text(data = wilcox_df,
              aes(x = sample_label, y = y_txt, label = stars),
              inherit.aes = FALSE, size = 3.6, fontface = "bold", colour = "grey20") +
    scale_fill_manual(values = cond_palette, guide = "none") +
    scale_y_continuous(labels = function(x) paste0(round(x), " um")) +
    labs(
      title = paste0("Distance Microglia (", mg_label, ") -> T cells / Immune niche"),
      subtitle = "Wilcoxon vs Mock",
      x = NULL,
      y = "Distance to nearest T cell (um)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8, colour = "grey40"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )

  # DAM by bins (1wpi vs 6wpi)
  dam_sub <- score_df %>%
    filter(niche_status == niche_status_value, sample %in% c("LCMV_1wpi", "LCMV_6wpi")) %>%
    mutate(sample_label = factor(SAMPLE_LABELS[sample], levels = c("LCMV 1 wpi", "LCMV 6 wpi"))) %>%
    select(cell_id_orig, sample, sample_label, dam_score = Upregulated_DAM)

  dist_1_6 <- dist_sub %>%
    filter(sample %in% c("LCMV_1wpi", "LCMV_6wpi")) %>%
    select(cell_id, sample, dist_um)

  joined <- inner_join(dam_sub, dist_1_6, by = c("cell_id_orig" = "cell_id", "sample"))
  joined$dist_bin <- cut(joined$dist_um,
                         breaks = DIST_BINS,
                         labels = DIST_LABELS,
                         right = FALSE,
                         include.lowest = TRUE)

  bin_stats <- joined %>%
    filter(!is.na(dist_bin)) %>%
    group_by(sample_label, dist_bin) %>%
    summarise(
      mean_dam = mean(dam_score, na.rm = TRUE),
      sem_dam = sd(dam_score, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    )

  dodge <- position_dodge(width = 0.8)

  p_bins <- ggplot(bin_stats,
                   aes(x = dist_bin, y = mean_dam,
                       fill = sample_label, group = sample_label)) +
    geom_bar(stat = "identity", position = dodge, width = 0.7, alpha = 0.85) +
    geom_errorbar(aes(ymin = mean_dam - sem_dam, ymax = mean_dam + sem_dam),
                  position = dodge, width = 0.25, linewidth = 0.5) +
    geom_text(aes(label = paste0("n=", n), y = -Inf),
              position = dodge, vjust = -0.3, size = 2.5, colour = "grey40") +
    scale_fill_manual(values = c("LCMV 1 wpi" = "#56B4E9", "LCMV 6 wpi" = "#D55E00"),
                      name = "Time point") +
    labs(
      title = paste0(title_prefix, " DAM score by distance bin"),
      subtitle = "Upregulated_DAM | 43 panel genes",
      x = "Distance to nearest T cell / Immune niche",
      y = "Mean DAM module score +- SEM"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8, colour = "grey40"),
      axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )

  p_combined <- p_violin + p_bins +
    plot_layout(widths = c(1.1, 1.4)) +
    plot_annotation(
      title = paste0(title_prefix, " microglia distance and DAM summary"),
      theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
    )

  save_fig(p_combined, out_stub, 12, 5)

  dist_summary <- dist_sub %>%
    group_by(sample) %>%
    summarise(n = n(), median_dist = median(dist_um, na.rm = TRUE), .groups = "drop")

  bin_summary <- bin_stats %>%
    mutate(sample = recode(as.character(sample_label),
                           "LCMV 1 wpi" = "LCMV_1wpi",
                           "LCMV 6 wpi" = "LCMV_6wpi")) %>%
    select(sample, dist_bin, n, mean_dam, sem_dam)

  out_summary <- bind_rows(
    dist_summary %>%
      mutate(metric = "distance", dist_bin = NA_character_, value = median_dist, sem = NA_real_) %>%
      select(metric, sample, dist_bin, n, value, sem),
    bin_summary %>%
      mutate(metric = "dam_by_bin", value = mean_dam, sem = sem_dam) %>%
      select(metric, sample, dist_bin, n, value, sem)
  )

  wilcox_summary <- wilcox_df %>%
    mutate(
      metric = "distance_wilcox_vs_mock",
      sample = recode(as.character(sample_label),
                      "LCMV 1 wpi" = "LCMV_1wpi",
                      "LCMV 3 wpi" = "LCMV_3wpi",
                      "LCMV 6 wpi" = "LCMV_6wpi"),
      dist_bin = NA_character_,
      n = NA_integer_,
      value = pval,
      sem = NA_real_
    ) %>%
    select(metric, sample, dist_bin, n, value, sem)

  out_summary <- bind_rows(out_summary, wilcox_summary)

  readr::write_csv(out_summary, file.path(out_dir, paste0("summary_", ifelse(niche_status_value == "In niche", "inniche", "outniche"), "_distance_dam_matched.csv")))

  list(dist_summary = dist_summary, bin_summary = bin_summary, wilcox = wilcox_df)
}

res_in <- build_figure(
  dist_sub = dist_in,
  niche_status_value = "In niche",
  mg_label = "C1qa",
  out_stub = "fig_combined_inniche_distance_dam_matched",
  title_prefix = "IN-niche"
)

res_out <- build_figure(
  dist_sub = dist_out,
  niche_status_value = "Out niche",
  mg_label = "P2ry12",
  out_stub = "fig_combined_outniche_distance_dam_matched",
  title_prefix = "OUT-of-niche"
)

message("Done: generated matched IN and OUT figures with recomputed distances.")
