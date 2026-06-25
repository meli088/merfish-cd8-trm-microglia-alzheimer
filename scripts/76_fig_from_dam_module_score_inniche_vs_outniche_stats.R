#!/usr/bin/env Rscript
# ============================================================================
# Script 76 — Summary figure from dam_module_score_inniche_vs_outniche_stats
# ============================================================================

set.seed(1997)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(Cairo)
})

source("scripts/00_palette.R")

in_file <- "outputs/banksy/immune_acod1/analysis/figures/dam_module_score_inniche_vs_outniche_stats_6wpi_vs_1wpi.csv"
out_dir <- "outputs/banksy/immune_acod1/analysis/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_file)) stop("Missing input CSV: ", in_file)

df <- readr::read_csv(in_file, show_col_types = FALSE) %>%
  mutate(
    niche_type = factor(niche_type, levels = c("IN-niche (C1qa)", "OUT-niche (P2ry12)")),
    signature = factor(signature, levels = c("Upregulated DAM", "Downregulated DAM")),
    y_star = pmax(median_1wpi, median_6wpi) + 0.12
  )

plot_df <- bind_rows(
  df %>% transmute(niche_type, signature, timepoint = "LCMV 1 wpi", median = median_1wpi, n = n_1wpi, sig, y_star),
  df %>% transmute(niche_type, signature, timepoint = "LCMV 6 wpi", median = median_6wpi, n = n_6wpi, sig, y_star)
) %>%
  mutate(timepoint = factor(timepoint, levels = c("LCMV 1 wpi", "LCMV 6 wpi")))

p <- ggplot(plot_df, aes(x = niche_type, y = median, color = timepoint, group = timepoint)) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  geom_line(aes(group = timepoint), alpha = 0.5, linewidth = 0.6,
            position = position_dodge(width = 0.4)) +
  geom_text(aes(label = paste0("n=", n)),
            position = position_dodge(width = 0.4), vjust = -1.1, size = 3, color = "grey35") +
  geom_text(
    data = df,
    aes(x = niche_type, y = y_star, label = sig),
    inherit.aes = FALSE,
    size = 4.5,
    fontface = "bold",
    color = "black"
  ) +
  facet_wrap(~signature, nrow = 1, scales = "free_y") +
  scale_color_manual(values = c("LCMV 1 wpi" = "#56B4E9", "LCMV 6 wpi" = "#F28E2B"), name = "Timepoint") +
  labs(
    title = "DAM module score (6wpi vs 1wpi): In-niche vs Out-of-niche",
    subtitle = "Medians from Wilcoxon summary table; stars indicate p-value significance",
    x = NULL,
    y = "Median DAM module score"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    strip.text = element_text(face = "bold")
  )

pdf_out <- file.path(out_dir, "fig_dam_module_score_inniche_vs_outniche_6wpi_vs_1wpi_summary.pdf")
jpg_out <- file.path(out_dir, "fig_dam_module_score_inniche_vs_outniche_6wpi_vs_1wpi_summary.jpg")

CairoPDF(pdf_out, width = 10.5, height = 5)
print(p)
dev.off()

CairoJPEG(jpg_out, width = 10.5 * 170, height = 5 * 170, res = 170)
print(p)
dev.off()

message("Saved summary figure:")
message(pdf_out)
message(jpg_out)
