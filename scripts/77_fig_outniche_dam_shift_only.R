#!/usr/bin/env Rscript
# ============================================================================
# Script 77 — OUT-of-niche DAM shift only (1wpi vs 6wpi)
#
# Goal: single figure illustrating statement:
#   - Upregulated DAM decreases in OUT-niche (P2ry12)
#   - Downregulated DAM increases in OUT-niche (P2ry12)
#   - with exact p-values from authoritative CSV
#
# Input:
#   outputs/banksy/immune_acod1/analysis/figures/
#     dam_module_score_inniche_vs_outniche_stats_6wpi_vs_1wpi.csv
#
# Output:
#   outputs/banksy/immune_acod1/analysis/figures/
#     fig_dam_outniche_shift_6wpi_vs_1wpi.pdf
#     fig_dam_outniche_shift_6wpi_vs_1wpi.jpg
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

fmt_p <- function(p) {
  formatC(p, format = "e", digits = 2)
}

df <- readr::read_csv(in_file, show_col_types = FALSE) %>%
  filter(niche_type == "OUT-niche (P2ry12)") %>%
  mutate(
    signature = factor(signature, levels = c("Upregulated DAM", "Downregulated DAM")),
    direction = ifelse(delta < 0, "decrease", "increase"),
    p_label = paste0("p = ", fmt_p(p_value_wilcoxon)),
    n_label = paste0("n1=", n_1wpi, " | n6=", n_6wpi)
  )

plot_df <- bind_rows(
  df %>% transmute(signature, timepoint = "LCMV 1 wpi", median = median_1wpi, p_label, n_label),
  df %>% transmute(signature, timepoint = "LCMV 6 wpi", median = median_6wpi, p_label, n_label)
) %>%
  mutate(timepoint = factor(timepoint, levels = c("LCMV 1 wpi", "LCMV 6 wpi")))

ann_df <- df %>%
  mutate(
    y_top = pmax(median_1wpi, median_6wpi) + 0.08,
    arrow_label = ifelse(direction == "decrease", "decrease", "increase")
  )

p <- ggplot(plot_df, aes(x = timepoint, y = median, color = timepoint, group = 1)) +
  geom_line(linewidth = 1.0, alpha = 0.7) +
  geom_point(size = 3.8) +
  geom_text(aes(label = sprintf("%.3f", median)), vjust = -0.9, size = 3.2, color = "grey25") +
  geom_text(
    data = ann_df,
    aes(x = 1.5, y = y_top, label = p_label),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4.1,
    color = "black"
  ) +
  geom_text(
    data = ann_df,
    aes(x = 1.5, y = y_top - 0.045, label = arrow_label),
    inherit.aes = FALSE,
    size = 3.2,
    color = "grey30"
  ) +
  geom_text(
    data = ann_df,
    aes(x = 1.5, y = y_top - 0.09, label = n_label),
    inherit.aes = FALSE,
    size = 3.0,
    color = "grey35"
  ) +
  facet_wrap(~signature, nrow = 1, scales = "free_y") +
  scale_color_manual(values = c("LCMV 1 wpi" = "#56B4E9", "LCMV 6 wpi" = "#F28E2B"), guide = "none") +
  labs(
    title = "OUT-of-niche microglia (P2ry12): DAM shift from 1wpi to 6wpi",
    subtitle = "Highly significant shift toward homeostasis",
    x = NULL,
    y = "Median DAM module score"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

pdf_out <- file.path(out_dir, "fig_dam_outniche_shift_6wpi_vs_1wpi.pdf")
jpg_out <- file.path(out_dir, "fig_dam_outniche_shift_6wpi_vs_1wpi.jpg")

CairoPDF(pdf_out, width = 9.2, height = 4.6)
print(p)
dev.off()

CairoJPEG(jpg_out, width = 9.2 * 170, height = 4.6 * 170, res = 170)
print(p)
dev.off()

message("Saved:")
message(pdf_out)
message(jpg_out)
