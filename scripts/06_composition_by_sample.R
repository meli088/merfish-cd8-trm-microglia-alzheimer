#!/usr/bin/env Rscript
# ============================================================================
# STEP 1: Composition Analysis by Sample
# ============================================================================
# Purpose:
#   Calculate and visualize annotation composition per sample using the 
#   annotated BANKSY clustering from lambda=0.2, resolution=0.9.
#   Produces a long-format CSV table and publication-style stacked barplot.
#
# Inputs:
#   - objects/04_banksy_joint_lam08_after_bloc3.rds
#   - ncells_by_sample_lam02_res09_joint_long.csv (annotation mapping)
#
# Outputs:
#   - composition_by_sample_long.csv (table with counts and percentages)
#   - composition_by_sample_stackedbar.pdf (publication style)
#   - composition_by_sample_stackedbar.jpg (publication style)
# ============================================================================

library(Banksy)
library(SpatialExperiment)
library(Seurat)
library(tidyverse)
library(ggplot2)

# Set up paths
base_path <- normalizePath(".")  # Run this script from the project root directory
setwd(base_path)

# Palette et ordre des annotations (fichier partagé)
source("scripts/00_palette.R")

# Create output directory if it doesn't exist
output_dir <- "outputs/banksy/umap_annotated"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================================
# STEP 1: Load SpatialExperiment object and extract clustering
# ============================================================================
cat("\n[Step 1] Loading SpatialExperiment object...\n")

# Load the object (with both Harmony embeddings and UMAP)
obj_path <- "objects/04_banksy_joint_lam08_after_bloc3.rds"
se_joint <- readRDS(obj_path)
cat("  - Loaded:", obj_path, "\n")
cat("  - Dimensions:", nrow(se_joint), "genes x", ncol(se_joint), "cells\n")

# Extract clustering column for lambda=0.2, resolution=0.9
# Using the helper function to find the correct column name
find_cl_col <- function(se, lam, res) {
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  cn <- colnames(cd)
  lam_pat <- paste0("lam", gsub("\\.", "\\\\.", as.character(lam)))
  res_pat <- paste0("_res", gsub("\\.", "\\\\.", as.character(res)), "$")
  cols <- cn[grepl(lam_pat, cn) & grepl(res_pat, cn)]
  if (length(cols) == 0) {
    stop("Clustering column not found for lambda=", lam, " resolution=", res)
  }
  return(cols[1])
}

lam <- 0.2
res <- 0.9
cluster_col <- find_cl_col(se_joint, lam, res)
cat("  - Found clustering column:", cluster_col, "\n")

# Extract metadata
cd <- as.data.frame(SummarizedExperiment::colData(se_joint))
banksy_clusters <- as.numeric(cd[[cluster_col]])  # Convert to numeric to handle the mapping
samples <- as.character(cd$sample)
cell_ids <- colnames(se_joint)

cat("  - Unique samples:", paste(unique(samples), collapse = ", "), "\n")
cat("  - Unique clusters:", length(unique(banksy_clusters)), "\n")

# ============================================================================
# STEP 2: Load annotation mapping from CSV
# ============================================================================
cat("\n[Step 2] Loading annotation mapping from CSV...\n")

csv_path <- "ncells_by_sample_lam02_res09_joint_long.csv"
anno_data <- read.delim(csv_path, sep = ";", stringsAsFactors = FALSE)
cat("  - Loaded CSV:", csv_path, "\n")
cat("  - Dimensions:", nrow(anno_data), "rows x", ncol(anno_data), "columns\n")
cat("  - Columns:", paste(colnames(anno_data), collapse = ", "), "\n")

# Create annotation mapping: banksy_domain -> annotation
# Extract unique domain-annotation pairs (exclude empty annotations)
anno_map <- anno_data %>%
  filter(annotation != "" & !is.na(annotation)) %>%
  mutate(annotation = trimws(annotation)) %>%  # Remove leading/trailing whitespace
  select(banksy_domain, annotation) %>%
  distinct() %>%
  arrange(banksy_domain)

cat("  - Unique annotations:", nrow(anno_map), "\n")
cat("  - Annotation labels:\n")
for (i in seq_len(nrow(anno_map))) {
  cat("    -", anno_map$banksy_domain[i], "->", anno_map$annotation[i], "\n")
}

# Create a lookup table using numeric cluster IDs
# Extract numeric part from "Domain_X" format
anno_map$cluster_id <- as.numeric(gsub("Domain_", "", anno_map$banksy_domain))
anno_lookup <- setNames(anno_map$annotation, anno_map$cluster_id)

# ============================================================================
# STEP 3: Create cell metadata with annotations
# ============================================================================
cat("\n[Step 3] Creating cell metadata with annotations...\n")

# Create domain names from numeric cluster IDs
banksy_domain_names <- paste0("Domain_", banksy_clusters)

cell_meta <- tibble(
  cell_id = cell_ids,
  sample = samples,
  banksy_cluster_id = banksy_clusters,
  banksy_domain = banksy_domain_names,
  annotation = trimws(anno_lookup[as.character(banksy_clusters)])  # Trim whitespace from lookup
)

# Check for unmapped clusters (these will have NA annotations)
n_unmapped <- sum(is.na(cell_meta$annotation))
if (n_unmapped > 0) {
  cat("  WARNING: ", n_unmapped, "cells have unmapped clusters\n")
  unmapped_domains <- unique(cell_meta$banksy_domain[is.na(cell_meta$annotation)])
  cat("  Unmapped domains:", paste(unmapped_domains, collapse = ", "), "\n")
}

cat("  - Total cells:", nrow(cell_meta), "\n")
cat("  - Cells with annotation:", sum(!is.na(cell_meta$annotation)), "\n")

# ============================================================================
# STEP 4: Calculate composition by sample and annotation
# ============================================================================
cat("\n[Step 4] Calculating composition by sample...\n")

# Define sample order for consistent visualization
sample_order <- c("mock_6wpi", "LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")

# Count cells per sample per annotation
composition <- cell_meta %>%
  filter(!is.na(annotation)) %>%  # Exclude unmapped clusters
  mutate(annotation = trimws(annotation)) %>%  # Final trim to ensure consistency
  group_by(sample, annotation) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  # Calculate percentage within each sample
  group_by(sample) %>%
  mutate(
    total_cells_in_sample = sum(n_cells),
    pct_cells = 100 * n_cells / total_cells_in_sample
  ) %>%
  ungroup() %>%
  # Ensure consistent sample order
  mutate(sample = factor(sample, levels = sample_order)) %>%
  arrange(sample, desc(n_cells)) %>%
  # Create output format: sample, annotation, n_cells, pct_cells
  select(sample, annotation, n_cells, pct_cells)

## (debug prints removed)

# Display composition summary by sample
cat("\n  Composition summary by sample:\n")
for (s in sample_order) {
  sample_comp <- composition %>% filter(sample == s)
  total_cells <- sum(sample_comp$n_cells)
  cat("\n  Sample:", s, "(", total_cells, "cells total )\n")
  
  # Show top 5 most abundant annotations
  top5 <- sample_comp %>% slice_head(n = 5)
  for (i in seq_len(nrow(top5))) {
    row <- top5[i, ]
    cat("    -", row$annotation, ":", row$n_cells, "cells", 
      "(", round(row$pct_cells, 2), "%)\n")
  }

  # Always show Htr2c if present and not in top 5
  htr2c_row <- sample_comp %>% filter(grepl("htr2c", tolower(annotation)))
  htr2c_in_top5 <- any(grepl("htr2c", tolower(top5$annotation)))
  
  if (nrow(htr2c_row) > 0 && !htr2c_in_top5) {
    cat("    -", htr2c_row$annotation[1], ":", htr2c_row$n_cells[1], "cells",
      "(", round(htr2c_row$pct_cells[1], 2), "%) [not in top5]\n")
  }
  
  if (nrow(htr2c_row) == 0) {
    cat("    - WARNING: Htr2c not found in composition for", s, "\n")
  }

}

# ============================================================================
# STEP 5: Save composition table as long-format CSV
# ============================================================================
cat("\n[Step 5] Saving composition table...\n")

composition_long <- composition %>%
  mutate(
    sample = as.character(sample),
    pct_cells = round(pct_cells, 4)
  )

csv_output <- file.path(output_dir, "composition_by_sample_long.csv")
write.csv(composition_long, csv_output, row.names = FALSE)
cat("  - Saved to:", csv_output, "\n")

# ============================================================================
# STEP 6: Create stacked barplot
# ============================================================================
cat("\n[Step 6] Creating stacked barplot...\n")

# Create color palette for annotations
# Use a consistent set of colors for each unique annotation
annotations_unique <- unique(composition_long$annotation)
n_annos <- length(annotations_unique)

cat("  - Annotations in composition_long for plot:", n_annos, "unique annotations\n")

# Utilise GLOBAL_PALETTE ; fallback grey70 pour annotations non référencées
# Ordonner les annotations selon ANNOTATION_ORDER
anno_ordered <- order_annotations(annotations_unique)

# Préparer les données avec le bon ordre de facteur
composition_plot <- composition_long %>%
  mutate(
    sample     = factor(sample, levels = sample_order),
    annotation = factor(annotation, levels = anno_ordered)
  ) %>%
  arrange(sample, annotation)

## (debug prints removed)

# Créer le stacked barplot
p <- ggplot(composition_plot, aes(x = sample, y = pct_cells, fill = annotation)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(
    values = GLOBAL_PALETTE,
    na.value = "grey70",
    name = "Annotation",
    guide = guide_legend(
      ncol = 1,
      title.position = "top",
      label.hjust = 0
    )
  ) +
  labs(
    x = "Sample",
    y = "Percentage of cells (%)",
    title = "Cell type composition by sample (BANKSY λ=0.2, res=0.9)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.box.margin = margin(l = 10),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
  ) +
  coord_cartesian(ylim = c(0, 100))

# ============================================================================
# STEP 7: Save plot as PDF and JPG
# ============================================================================
cat("\n[Step 7] Saving figures...\n")

# PDF output
pdf_output <- file.path(output_dir, "composition_by_sample_stackedbar.pdf")
ggsave(
  pdf_output,
  plot = p,
  width = 10,
  height = 7,
  dpi = 300,
  device = "pdf"
)
cat("  - Saved PDF to:", pdf_output, "\n")

# JPG output
jpg_output <- file.path(output_dir, "composition_by_sample_stackedbar.jpg")
ggsave(
  jpg_output,
  plot = p,
  width = 10,
  height = 7,
  dpi = 300,
  device = "jpg"
)
cat("  - Saved JPG to:", jpg_output, "\n")

# ============================================================================
# STEP 8: Print final summary
# ============================================================================
cat("\n")
cat(paste(c("", rep("=", 78)), collapse = ""), "\n")
cat("STEP 1 COMPLETE: Composition Analysis by Sample\n")
cat(paste(c("", rep("=", 78)), collapse = ""), "\n\n")

cat("Summary Statistics:\n")
cat("  - Total cells analyzed:", sum(composition_long$n_cells), "\n")
cat("  - Number of samples:", length(unique(composition_long$sample)), "\n")
cat("  - Number of unique annotations:", length(unique(composition_long$annotation)), "\n")
cat("  - Average cells per annotation per sample:", 
    round(mean(composition_long$n_cells), 1), "\n")

cat("\nOutputs generated:\n")
cat("  1.", csv_output, "\n")
cat("  2.", pdf_output, "\n")
cat("  3.", jpg_output, "\n")

cat("\nNext: Run enrichment analysis and/or spatial distribution analysis.\n")
cat("=" %+% paste(rep("=", 78), collapse = ""), "\n\n")

cat("\n✓ Step 1 finished successfully.\n")
