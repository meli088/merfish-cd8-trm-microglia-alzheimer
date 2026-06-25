#!/usr/bin/env Rscript
# ============================================================================
# Script 68 — DAM temporal changes: Summary barplot 6wpi vs 1wpi
#             Shows delta changes (Δ = median_6wpi - median_1wpi) per niche
#
# Outputs:
#   - fig_dam_temporal_changes_6wpi_vs_1wpi.pdf/.jpg
#   - dam_temporal_summary_6wpi_vs_1wpi.csv
#
# ============================================================================

set.seed(1997)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(Cairo)
})

# Load stats from script 67
stats_path <- "outputs/banksy/immune_acod1/analysis/figures/dam_module_score_inniche_vs_outniche_stats_6wpi_vs_1wpi.csv"
stats <- read_csv(stats_path, show_col_types = FALSE)

# Prepare for visualization
stats_clean <- stats %>%
  mutate(
    niche_short = case_when(
      grepl("IN-niche", niche_type) ~ "IN-niche\n(C1qa)",
      grepl("OUT-niche", niche_type) ~ "OUT-niche\n(P2ry12)",
      TRUE ~ niche_type
    ),
    sig_short = case_when(
      grepl("Upregulated", signature) ~ "Upregulated\nDAM",
      grepl("Downregulated", signature) ~ "Downregulated\nDAM",
      TRUE ~ signature
    ),
    # Color by significance
    sig_color = case_when(
      sig == "***" ~ "#D55E00",  # Orange for *** (OUT-niche)
      sig == "n.s." ~ "grey70",  # Grey for n.s. (IN-niche)
      TRUE ~ "grey50"
    ),
    # Direction
    direction = ifelse(delta > 0, "Increase", "Decrease")
  ) %>%
  select(niche_short, sig_short, delta, p_value_wilcoxon, sig, sig_color, direction, n_1wpi, n_6wpi)

# Save summary table
dir_out <- "outputs/banksy/immune_acod1/analysis/figures"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

summary_out <- stats_clean %>%
  select(niche_short, sig_short, n_1wpi, n_6wpi, delta, p_value_wilcoxon, sig)

write_csv(summary_out, file.path(dir_out, "dam_temporal_summary_6wpi_vs_1wpi.csv"))
message(sprintf("✓ Saved summary: dam_temporal_summary_6wpi_vs_1wpi.csv"))

# Create figure
p <- ggplot(stats_clean, aes(x = niche_short, y = delta, fill = sig_color)) +
  geom_col(width = 0.7, alpha = 0.85, colour = "grey30", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "solid", colour = "black", linewidth = 0.3) +
  facet_wrap(~sig_short, scales = "fixed") +
  
  # Add text labels with values and significance
  geom_text(
    aes(label = sprintf("Δ=%.3f\np%s", delta, sig)),
    vjust = ifelse(stats_clean$delta > 0, -0.5, 1.5),
    size = 3.5,
    fontface = "bold",
    colour = "grey20"
  ) +
  
  scale_fill_identity(guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  
  labs(
    title = "DAM Module Score Changes: 6wpi vs 1wpi",
    subtitle = "Δ = median(6wpi) - median(1wpi) | Wilcoxon test",
    x = NULL,
    y = "Δ DAM Score"
  ) +
  
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 9, colour = "grey50", hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "white", colour = "grey30"),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 10, face = "bold")
  )

# Save figure
pdf_path <- file.path(dir_out, "fig_dam_temporal_changes_6wpi_vs_1wpi.pdf")
jpg_path <- file.path(dir_out, "fig_dam_temporal_changes_6wpi_vs_1wpi.jpg")

CairoPDF(pdf_path, width = 8, height = 5)
print(p)
dev.off()

CairoJPEG(jpg_path, width = 8 * 150, height = 5 * 150, res = 150)
print(p)
dev.off()

message(sprintf("✓ Saved figure: fig_dam_temporal_changes_6wpi_vs_1wpi"))

# Print text summary
cat("\n=== DAM Temporal Changes (6wpi vs 1wpi) ===\n\n")
for (i in 1:nrow(stats_clean)) {
  row <- stats_clean[i, ]
  cat(sprintf(
    "%s | %s\n  Δ = %.4f (p=%s, %s) | n_1wpi=%d, n_6wpi=%d\n\n",
    row$niche_short, row$sig_short,
    row$delta, 
    ifelse(row$p_value_wilcoxon < 0.0001, "<0.0001", sprintf("%.4f", row$p_value_wilcoxon)),
    row$sig, row$n_1wpi, row$n_6wpi
  ))
}
