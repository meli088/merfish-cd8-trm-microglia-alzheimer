#!/usr/bin/env Rscript
# ============================================================================
# Script 78 — OUT-of-niche microglia (P2ry12): simple violin 6wpi vs 1wpi
# Style aligned with LOCKED in-niche figure
#
# Inputs:
#   outputs/banksy/microglia_dam_niche/module_scores_per_cell.csv
#   outputs/banksy/immune_acod1/analysis/figures/dam_module_score_inniche_vs_outniche_stats_6wpi_vs_1wpi.csv
#
# Outputs:
#   outputs/banksy/immune_acod1/analysis/figures/
#     fig_dam_score_outniche_microglia_6wpi_vs_1wpi_simple_violin.pdf
#     fig_dam_score_outniche_microglia_6wpi_vs_1wpi_simple_violin.jpg
#     dam_score_outniche_microglia_6wpi_vs_1wpi_simple_violin_stats.csv
# ============================================================================

set.seed(1997)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

source("scripts/00_palette.R")

scores_file <- "outputs/banksy/microglia_dam_niche/module_scores_per_cell.csv"
stats_file <- "outputs/banksy/immune_acod1/analysis/figures/dam_module_score_inniche_vs_outniche_stats_6wpi_vs_1wpi.csv"
out_dir <- "outputs/banksy/immune_acod1/analysis/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(scores_file)) stop("Missing scores file: ", scores_file)
if (!file.exists(stats_file)) stop("Missing stats file: ", stats_file)

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

# OUT-of-niche per-cell scores
score_df <- readr::read_csv(scores_file, show_col_types = FALSE) %>%
  filter(niche_status == "Out niche", sample %in% c("LCMV_1wpi", "LCMV_6wpi")) %>%
  transmute(
    sample,
    timepoint = factor(TIME_LABELS[sample], levels = c("LCMV 1 wpi", "LCMV 6 wpi")),
    Upregulated_DAM,
    Downregulated_DAM
  )

# Pull exact p-values from authoritative in/out CSV
stats_out <- readr::read_csv(stats_file, show_col_types = FALSE) %>%
  filter(niche_type == "OUT-niche (P2ry12)") %>%
  select(signature, n_1wpi, n_6wpi, median_1wpi, median_6wpi, delta, p_value_wilcoxon, sig)

p_up <- stats_out$p_value_wilcoxon[stats_out$signature == "Upregulated DAM"]
p_down <- stats_out$p_value_wilcoxon[stats_out$signature == "Downregulated DAM"]

make_panel <- function(df, col_name, title_txt, pval) {
  y_max <- max(df[[col_name]], na.rm = TRUE)
  y_min <- min(df[[col_name]], na.rm = TRUE)
  y_txt <- y_max + 0.08 * (y_max - y_min)

  ggplot(df, aes(x = timepoint, y = .data[[col_name]], fill = timepoint)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.6, color = NA) +
    geom_boxplot(width = 0.2, outlier.size = 0.3, outlier.alpha = 0.25,
                 color = "grey20", fill = "white") +
    annotate("text", x = 1.5, y = y_txt, label = pval_stars(pval),
             size = 5, fontface = "bold") +
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
    title = "Out-of-niche microglia (P2ry12): DAM score 6wpi vs 1wpi",
    subtitle = sprintf("n(1wpi)=%d, n(6wpi)=%d", n_by_time[["LCMV 1 wpi"]], n_by_time[["LCMV 6 wpi"]]),
    theme = theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 9, color = "grey35", hjust = 0.5)
    )
  )

pdf_out <- file.path(out_dir, "fig_dam_score_outniche_microglia_6wpi_vs_1wpi_simple_violin.pdf")
jpg_out <- file.path(out_dir, "fig_dam_score_outniche_microglia_6wpi_vs_1wpi_simple_violin.jpg")

# cairo_pdf requested
ggsave(pdf_out, p_final, width = 10.5, height = 4.8, device = cairo_pdf)
ggsave(jpg_out, p_final, width = 10.5, height = 4.8, dpi = 300)

readr::write_csv(stats_out, file.path(out_dir, "dam_score_outniche_microglia_6wpi_vs_1wpi_simple_violin_stats.csv"))

message("Saved:")
message(pdf_out)
message(jpg_out)
