#!/usr/bin/env Rscript
# ============================================================================
# Script 71 — Matched IN/OUT niche figures: distance violin + DAM bin barplot
#
# Uses existing outputs from script 39:
#   - outputs/banksy/microglia_dam_niche/distance_microglia_to_tcell_per_cell.csv
#   - outputs/banksy/microglia_dam_niche/module_scores_per_cell.csv
#
# Produces two figures with identical layout/style:
#   - fig_combined_inniche_distance_dam_matched.pdf/.jpg
#   - fig_combined_outniche_distance_dam_matched.pdf/.jpg
#
# Plus summary tables:
#   - summary_inniche_distance_dam_matched.csv
#   - summary_outniche_distance_dam_matched.csv
#
# ============================================================================

set.seed(1997)

suppressPackageStartupMessages({
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

# Keep bins exactly like the existing displayed IN figure
DIST_BINS <- c(0, 50, 100, 200)
DIST_LABELS <- c("<50µm", "50-100µm", "100-200µm")

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

message("Loading distance and module score tables...")
dist_df <- readr::read_csv(
  file.path(out_dir, "distance_microglia_to_tcell_per_cell.csv"),
  show_col_types = FALSE
)
score_df <- readr::read_csv(
  file.path(out_dir, "module_scores_per_cell.csv"),
  show_col_types = FALSE
)

score_df <- score_df %>%
  mutate(cell_id_orig = sub("^(out|in)_", "", cell_id))

build_matched_figure <- function(niche_status_value, mg_label, out_stub) {
  message("Building figure for: ", niche_status_value)

  # -------- Violin panel (distance over conditions) --------
  dist_sub <- dist_df %>%
    filter(niche_status == niche_status_value, !is.na(dist_um)) %>%
    mutate(sample_label = factor(SAMPLE_LABELS[sample], levels = SAMPLE_LABELS[SAMPLE_ORDER]))

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
    scale_y_continuous(labels = function(x) paste0(x, " µm")) +
    labs(
      title = paste0("Distance Microglia (", mg_label, ") → T cells / Immune niche"),
      subtitle = "Wilcoxon vs Mock",
      x = NULL,
      y = "Distance to nearest T cell (µm)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8, colour = "grey40"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )

  # -------- Bar panel (DAM by distance bins, 1wpi vs 6wpi) --------
  dam_sub <- score_df %>%
    filter(niche_status == niche_status_value, sample %in% c("LCMV_1wpi", "LCMV_6wpi")) %>%
    mutate(sample_label = factor(SAMPLE_LABELS[sample], levels = c("LCMV 1 wpi", "LCMV 6 wpi"))) %>%
    select(cell_id_orig, sample, sample_label, dam_score = Upregulated_DAM)

  dist_sub_1_6 <- dist_sub %>%
    filter(sample %in% c("LCMV_1wpi", "LCMV_6wpi")) %>%
    select(cell_id, sample, dist_um)

  joined <- inner_join(dam_sub, dist_sub_1_6, by = c("cell_id_orig" = "cell_id", "sample"))
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
      title = paste0(niche_status_value, " DAM score by distance bin"),
      subtitle = "Upregulated_DAM | 43 panel genes",
      x = "Distance to nearest T cell / Immune niche",
      y = "Mean DAM module score ± SEM"
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
      title = paste0(niche_status_value, " microglia distance and DAM summary"),
      theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
    )

  save_fig(p_combined, out_stub, 12, 5)

  # -------- Summary table --------
  dist_summary <- dist_sub %>%
    group_by(sample) %>%
    summarise(
      n = n(),
      median_dist = median(dist_um, na.rm = TRUE),
      .groups = "drop"
    )

  dam_summary <- bin_stats %>%
    mutate(sample = recode(as.character(sample_label),
                           "LCMV 1 wpi" = "LCMV_1wpi",
                           "LCMV 6 wpi" = "LCMV_6wpi")) %>%
    select(sample, dist_bin, n, mean_dam, sem_dam)

  out_summary <- dist_summary %>%
    mutate(metric = "distance") %>%
    rename(value = median_dist) %>%
    mutate(dist_bin = NA_character_, sem = NA_real_) %>%
    select(metric, sample, dist_bin, n, value, sem)

  out_summary <- bind_rows(
    out_summary,
    dam_summary %>%
      mutate(metric = "dam_by_bin", value = mean_dam, sem = sem_dam) %>%
      select(metric, sample, dist_bin, n, value, sem)
  )

  readr::write_csv(out_summary, file.path(out_dir, paste0("summary_", out_stub, ".csv")))

  list(
    wilcox_vs_mock = wilcox_df,
    dist_summary = dist_summary,
    bin_stats = bin_stats
  )
}

res_in <- build_matched_figure(
  niche_status_value = "In niche",
  mg_label = "C1qa",
  out_stub = "fig_combined_inniche_distance_dam_matched"
)

res_out <- build_matched_figure(
  niche_status_value = "Out niche",
  mg_label = "P2ry12",
  out_stub = "fig_combined_outniche_distance_dam_matched"
)

message("Done. Generated matched IN and OUT figures.")
