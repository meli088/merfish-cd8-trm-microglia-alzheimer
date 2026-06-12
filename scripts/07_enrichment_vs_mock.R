#!/usr/bin/env Rscript
# ============================================================================
# STEP: Enrichment vs mock (log2 fold-change of annotation % over time)
# ============================================================================
# Purpose:
#  - Calculate enrichment (log2 fold change) of annotation percentages in
#    infected samples (LCMV_1wpi, LCMV_3wpi, LCMV_6wpi) vs mock_6wpi.
#  - Produce a CSV table and two publication-style plots (all annotations + top N).
#
# Inputs:
#  - outputs/banksy/umap_annotated/composition_by_sample_long.csv
#  - ncells_by_sample_lam02_res09_joint_long.csv (for consistent annotation labels)
#
# Outputs:
#  - outputs/banksy/enrichment/enrichment_vs_mock_table.csv
#  - outputs/banksy/enrichment/enrichment_vs_mock_all.pdf/jpg
#  - outputs/banksy/enrichment/enrichment_vs_mock_top10.pdf/jpg
# ============================================================================

library(tidyverse)
library(ggplot2)

base_path <- normalizePath(".")  # Run this script from the project root directory
setwd(base_path)

# Palette et ordre des annotations (fichier partagé)
source("scripts/00_palette.R")

# Output dir
out_dir <- file.path("outputs","banksy","enrichment_vs_mock")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Read composition table (generated in previous step)
comp_csv <- file.path("outputs","banksy","umap_annotated","composition_by_sample_long.csv")
if (!file.exists(comp_csv)) {
  stop("Fichier composition introuvable : ", comp_csv,
       "\nVeuillez d'abord lancer le script 06_composition_by_sample.R")
}
comp <- read_csv(comp_csv, show_col_types = FALSE)
if (!"sample" %in% colnames(comp)) {
  stop("Colonne 'sample' absente du CSV. Colonnes pr\u00e9sentes : ",
       paste(colnames(comp), collapse = ", "))
}

# Ensure sample ordering and mapping to numeric timepoints
sample_map <- tibble(
  sample = c("mock_6wpi","LCMV_1wpi","LCMV_3wpi","LCMV_6wpi"),
  timepoint = c(0,1,3,6)
)

comp2 <- comp %>%
  left_join(sample_map, by = "sample")

# Use mock as reference: get pct per annotation in mock
eps <- 0.01  # pseudocount in percentage points to avoid division by zero
mock_pct <- comp2 %>%
  filter(sample == "mock_6wpi") %>%
  select(annotation, pct_mock = pct_cells)

# Calculate enrichment for infected samples
infected_samples <- c("LCMV_1wpi","LCMV_3wpi","LCMV_6wpi")
res_tbl <- comp2 %>%
  filter(sample %in% infected_samples) %>%
  left_join(mock_pct, by = "annotation") %>%
  mutate(
    pct_cells = pct_cells,
    pct_mock = ifelse(is.na(pct_mock), 0, pct_mock),
    pct_cells_adj = pct_cells + eps,
    pct_mock_adj = pct_mock + eps,
    log2_fc_vs_mock = log2(pct_cells_adj / pct_mock_adj),
    timepoint = as.integer(timepoint)
  ) %>%
  select(annotation, sample, timepoint, pct_cells, pct_mock, log2_fc_vs_mock) %>%
  arrange(annotation, timepoint)

# Save results
csv_out <- file.path(out_dir, "enrichment_vs_mock_table.csv")
write.csv(res_tbl, csv_out, row.names = FALSE)
cat("Saved enrichment table to:", csv_out, "\n")

# Ordre des annotations pour la légende
anno_ordered <- order_annotations(unique(res_tbl$annotation))
res_tbl$annotation <- factor(res_tbl$annotation, levels = anno_ordered)

# Plot: all annotations
p_all <- ggplot(res_tbl, aes(x = timepoint, y = log2_fc_vs_mock, group = annotation, color = annotation)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_line(size = 0.6, alpha = 0.9) +
  geom_point(size = 1.8) +
  scale_x_continuous(breaks = c(1,3,6), labels = c("1","3","6")) +
  scale_color_manual(values = GLOBAL_PALETTE, na.value = "grey70") +
  labs(x = "w.p.i.", y = "Log2FC of annotation percentage infected vs mock",
       title = "Annotation enrichment over time vs mock (log2FC)") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

pdf_all <- file.path(out_dir, "enrichment_vs_mock_all.pdf")
jpg_all <- file.path(out_dir, "enrichment_vs_mock_all.jpg")
ggsave(pdf_all, p_all, width = 10, height = 6)
ggsave(jpg_all, p_all, width = 10, height = 6, dpi = 300)
cat("Saved plots (all annotations) to:", pdf_all, jpg_all, "\n")

# Select top N annotations by max absolute log2FC across timepoints
top_n <- 10
top_annos <- res_tbl %>%
  group_by(annotation) %>%
  summarise(max_abs = max(abs(log2_fc_vs_mock), na.rm = TRUE)) %>%
  arrange(desc(max_abs)) %>%
  slice_head(n = top_n) %>%
  pull(annotation)

res_top <- res_tbl %>% filter(annotation %in% top_annos)

p_top <- ggplot(res_top, aes(x = timepoint, y = log2_fc_vs_mock, group = annotation, color = annotation)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = c(1,3,6), labels = c("1","3","6")) +
  scale_color_manual(values = GLOBAL_PALETTE, na.value = "grey70") +
  labs(x = "w.p.i.", y = "Log2FC of annotation percentage infected vs mock",
       title = paste0("Top ", top_n, " annotations by absolute Log2FC")) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

pdf_top <- file.path(out_dir, "enrichment_vs_mock_top10.pdf")
jpg_top <- file.path(out_dir, "enrichment_vs_mock_top10.jpg")
ggsave(pdf_top, p_top, width = 8, height = 5)
ggsave(jpg_top, p_top, width = 8, height = 5, dpi = 300)
cat("Saved plots (top annotations) to:", pdf_top, jpg_top, "\n")

cat("\nSTEP COMPLETE: enrichment vs mock.\n")
cat("Outputs saved in:", out_dir, "\n")
