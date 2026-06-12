#!/usr/bin/env Rscript
# =============================================================
# Script: 24_spatial_grid_100um_global.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-05
#
# Correction of script 23:
#   Script 23 used the immune subclustered object
#   (08_immune_annotated_lam02_res03.rds), which no longer contains
#   the broad annotation labels "Immune (Acod1)" and
#   "IFN responsive (Ifit1)". Those labels only exist in the global
#   BANKSY object (04_banksy_joint_lam08_after_bloc3.rds).
#
# Goal:
#   Perform 100 µm x 100 µm grid-based spatial co-occurrence analysis
#   between "Immune (Acod1)" and "IFN responsive (Ifit1)" using the
#   correct globally-annotated BANKSY object.
#
# Inputs:
#   objects/04_banksy_joint_lam08_after_bloc3.rds
#     — global BANKSY SpatialExperiment object
#     — lambda=0.8 joint object; contains lambda=0.2 / res=0.9 cluster col
#     — spatialCoords: sdimx / sdimy (µm)
#   ncells_by_sample_lam02_res09_joint_long.csv
#     — annotation mapping: banksy_domain -> annotation label
#
# Definitions:
#   Immune-positive zone : 100 µm bin with >= 1 "Immune (Acod1)" cell
#   IFN-positive zone    : 100 µm bin with >= 1 "IFN responsive (Ifit1)" cell
#   Immune-enriched zone : Immune-positive bin with n_immune strictly >
#                          mean(n_immune) among Immune-positive bins in
#                          that sample
#   IFN-enriched zone    : IFN-positive bin with n_ifn strictly >
#                          mean(n_ifn) among IFN-positive bins in that sample
#
# Outputs (folder: outputs/banksy/inflammatory_niche_step2_grid_100um_global/):
#   diagnostics.txt                         — object / label / cell count check
#   grid_per_cell.csv                       — per-cell grid assignment
#   grid_per_bin.csv                        — per-bin total/immune/ifn counts
#   grid_per_bin_enrichment_flags.csv       — enrichment flags per bin
#   grid_positive_zones_summary.csv         — positive zone counts per sample
#   grid_enriched_zones_summary.csv         — enriched zone counts + overlap
#   fig_grid_overlay_<sample>.pdf/jpg       — grid overlay (one sample)
#   fig_grid_immune_<sample>.pdf/jpg        — Immune heatmap per sample
#   fig_grid_ifn_<sample>.pdf/jpg           — IFN heatmap per sample
# =============================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(Banksy)
  library(tidyverse)
  library(ggplot2)
})

base_path <- normalizePath(".")  # Run this script from the project root directory
setwd(base_path)

out_dir <- file.path("outputs", "banksy",
                     "inflammatory_niche_step2_grid_100um_global")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE,
                                      showWarnings = FALSE)

# =============================================================
# Parameters
# =============================================================

GRID_SIZE_UM  <- 100
LAM           <- 0.2
RES_TARGET    <- 0.9

# Exact annotation labels to match (script stops if either is absent)
IMMUNE_LABEL  <- "Immune (Acod1)"
IFN_LABEL     <- "IFN responsive (Ifit1)"

SAMPLE_ORDER  <- c("mock_6wpi", "LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")
sample_labels <- c(mock_6wpi = "Mock 6 wpi", LCMV_1wpi = "LCMV 1 wpi",
                   LCMV_3wpi = "LCMV 3 wpi", LCMV_6wpi = "LCMV 6 wpi")

# Colour palette for global annotations
PALETTE <- c(
  "Immune (Acod1)"               = "#E41A1C",
  "Microglia (P2ry12)"           = "#377EB8",
  "Astrocytes (Fgfr3)"           = "#4DAF4A",
  "Astrocytes (Gfap)"            = "#984EA3",
  "IFN responsive (Ifit1)"       = "#A65628",
  "Oligodendrocytes (Plp1)"      = "#F781BF",
  "OPC (Pdgfra)"                 = "#FF7F00",
  "Vascular (Cldn5)"             = "#66C2A5",
  "Fibroblasts/VSMC (Col1a1)"    = "#8DA0CB",
  "Choroid plexus Epi (Enpp2)"   = "#FFD92F",
  "Prolif neural/glial (Ccdc153)"= "#A6D854",
  "Neurons (Adora2a)"            = "#1B9E77",
  "Neurons (Arhgap36)"           = "#D95F02",
  "Neurons (Crhbp)"              = "#7570B3",
  "Neurons (Dkk3)"               = "#E7298A",
  "Neurons (Htr2c)"              = "#66A61E",
  "Neurons (Nefm)"               = "#E6AB02",
  "Neurons (Neurod1)"            = "#A6761D",
  "Neurons (Rxfp1)"              = "#125988",
  "Neurons (Slc17a8)"            = "#FB9A99",
  "Non annote"                   = "#AAAAAA"
)

# =============================================================
# Theme / save helpers
# =============================================================

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1,
                                      hjust = 0, margin = margin(b = 4)),
      plot.subtitle    = element_text(size = base_size - 2, color = "grey40",
                                      hjust = 0, lineheight = 1.3,
                                      margin = margin(b = 8)),
      axis.text        = element_text(size = base_size - 2),
      axis.title       = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size - 1, face = "bold"),
      legend.text      = element_text(size = base_size - 1.5),
      plot.margin      = margin(10, 14, 10, 10)
    )
}

save_fig <- function(plot, name, width, height, dpi_jpg = 220) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")),
         plot = plot, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(out_dir, paste0(name, ".jpg")),
         plot = plot, width = width, height = height, dpi = dpi_jpg)
  message("  Saved: ", name, ".pdf/jpg")
}

# =============================================================
# 1. Load the global BANKSY object
# =============================================================

obj_file <- file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds")
stopifnot("Object file not found" = file.exists(obj_file))

message("Loading: ", obj_file)
se <- readRDS(obj_file)
message("  Loaded: ", ncol(se), " cells")
message("  Class : ", class(se)[1])
message("  reducedDims: ", paste(reducedDimNames(se), collapse = ", "))

cd <- as.data.frame(colData(se))

# =============================================================
# 2. Find BANKSY cluster column (lambda=0.2, res=0.9)
# =============================================================

find_cl_col <- function(se, lam, res) {
  all_cols <- clusterNames(se)
  lam_str  <- gsub("\\.", "\\\\.", as.character(lam))
  lam_cols <- all_cols[grep(paste0("lam", lam_str), all_cols)]
  if (length(lam_cols) == 0) return(NULL)
  res_vals <- suppressWarnings(as.numeric(sub(".*_res", "", lam_cols)))
  idx <- which(!is.na(res_vals) & abs(res_vals - res) < 1e-8)
  if (length(idx) == 0) return(NULL)
  lam_cols[idx[1]]
}

cl_col <- find_cl_col(se, LAM, RES_TARGET)
if (is.null(cl_col)) {
  stop("BANKSY cluster column not found for lambda=", LAM,
       " res=", RES_TARGET,
       "\nAvailable columns: ", paste(clusterNames(se), collapse = ", "))
}
message("  BANKSY cluster column: ", cl_col)

# =============================================================
# 3. Reconstruct annotation mapping from CSV
# =============================================================

csv_path <- "ncells_by_sample_lam02_res09_joint_long.csv"
stopifnot("Annotation CSV not found" = file.exists(csv_path))

message("Loading annotation CSV: ", csv_path)
annot_long <- read_delim(
  csv_path,
  delim      = ";",
  locale     = locale(decimal_mark = "."),
  show_col_types = FALSE,
  trim_ws    = TRUE
) %>%
  select(-matches("^Unnamed")) %>%
  mutate(
    banksy_domain = as.character(banksy_domain),
    annotation    = trimws(as.character(annotation))
  )

message("  CSV rows: ", nrow(annot_long))

# Build unique banksy_domain -> annotation mapping
annotation_map <- annot_long %>%
  filter(!is.na(annotation), annotation != "") %>%
  distinct(banksy_domain, annotation)

# Check for ambiguous mappings (one domain -> multiple annotations)
ambig <- annotation_map %>% count(banksy_domain) %>% filter(n > 1)
if (nrow(ambig) > 0) {
  cat("Ambiguous domain mappings:\n")
  print(annotation_map %>% semi_join(ambig, by = "banksy_domain"))
  stop("Ambiguous annotation mapping — fix the CSV before continuing.")
}

message("  Unique banksy_domain -> annotation mappings: ",
        nrow(annotation_map))
cat("  Annotations in map:\n")
cat(paste0("    - ", sort(unique(annotation_map$annotation)), "\n"),
    sep = "")

# =============================================================
# 4. Assign annotation to each cell
# =============================================================

domain_labels <- paste0("Domain_", as.character(colData(se)[[cl_col]]))

anno_lookup <- setNames(annotation_map$annotation,
                        annotation_map$banksy_domain)

annotation <- ifelse(
  !is.na(anno_lookup[domain_labels]) & anno_lookup[domain_labels] != "",
  anno_lookup[domain_labels],
  "Non annote"
)

cat("\nAnnotation distribution (all cells):\n")
print(sort(table(annotation), decreasing = TRUE))

# =============================================================
# 5. Verify that both target labels are present — stop if not
# =============================================================

if (!IMMUNE_LABEL %in% annotation) {
  stop(
    "CRITICAL: Label '", IMMUNE_LABEL, "' was NOT found in annotations.\n",
    "Labels present: ", paste(sort(unique(annotation)), collapse = ", ")
  )
}
if (!IFN_LABEL %in% annotation) {
  stop(
    "CRITICAL: Label '", IFN_LABEL, "' was NOT found in annotations.\n",
    "Labels present: ", paste(sort(unique(annotation)), collapse = ", ")
  )
}
message("  Confirmed: '", IMMUNE_LABEL, "' found (",
        sum(annotation == IMMUNE_LABEL), " cells)")
message("  Confirmed: '", IFN_LABEL,    "' found (",
        sum(annotation == IFN_LABEL),    " cells)")

# =============================================================
# 6. Extract spatial coordinates
# =============================================================

xy <- as.data.frame(SpatialExperiment::spatialCoords(se))
if (!all(c("sdimx", "sdimy") %in% colnames(xy))) {
  colnames(xy)[1:2] <- c("sdimx", "sdimy")
}

cat("\nSpatial coordinate ranges:\n")
cat("  sdimx:", round(min(xy$sdimx), 1), "to",
    round(max(xy$sdimx), 1), "µm\n")
cat("  sdimy:", round(min(xy$sdimy), 1), "to",
    round(max(xy$sdimy), 1), "µm\n")

# =============================================================
# 7. Build master cell table
# =============================================================

cell_df <- data.frame(
  cell_id    = colnames(se),
  sample     = factor(as.character(cd$sample), levels = SAMPLE_ORDER),
  annotation = annotation,
  sdimx      = xy$sdimx,
  sdimy      = xy$sdimy,
  is_immune  = annotation == IMMUNE_LABEL,
  is_ifn     = annotation == IFN_LABEL,
  stringsAsFactors = FALSE
)

cat("\nCells per sample:\n")
print(table(cell_df$sample))
cat("\nImmune and IFN cells per sample:\n")
print(
  cell_df %>%
    group_by(sample) %>%
    summarise(n_immune = sum(is_immune), n_ifn = sum(is_ifn), .groups = "drop")
)

# =============================================================
# 8. Write diagnostics text file
# =============================================================

diag_lines <- c(
  "=============================================================",
  "DIAGNOSTICS: 24_spatial_grid_100um_global.R",
  paste0("Date: ", Sys.Date()),
  "=============================================================",
  "",
  "--- Object ---",
  paste0("File : ", obj_file),
  paste0("Class: ", class(se)[1]),
  paste0("Total cells: ", ncol(se)),
  "",
  "--- Clustering column ---",
  paste0("lambda    : ", LAM),
  paste0("resolution: ", RES_TARGET),
  paste0("Column    : ", cl_col),
  "",
  "--- Annotation CSV ---",
  paste0("File: ", csv_path),
  paste0("Unique domain->annotation mappings: ", nrow(annotation_map)),
  "",
  "--- Target label confirmation ---",
  paste0("'", IMMUNE_LABEL, "' : FOUND — ",
         sum(annotation == IMMUNE_LABEL), " cells"),
  paste0("'", IFN_LABEL, "' : FOUND — ",
         sum(annotation == IFN_LABEL), " cells"),
  "",
  "--- All annotation labels present ---",
  paste0("  ", sort(unique(annotation))),
  "",
  "--- Grid ---",
  paste0("Bin size: ", GRID_SIZE_UM, " µm x ", GRID_SIZE_UM, " µm"),
  "Origin: per-sample (min sdimx, min sdimy of that sample)",
  "grid_x = floor((sdimx - min(sdimx)) / 100)  [0-based]",
  "grid_y = floor((sdimy - min(sdimy)) / 100)  [0-based]",
  "",
  "--- Enrichment rule ---",
  paste0("Immune-enriched: Immune-positive bin with"),
  paste0("  n_immune > mean(n_immune) across Immune-positive bins in that sample"),
  paste0("IFN-enriched: IFN-positive bin with"),
  paste0("  n_ifn > mean(n_ifn) across IFN-positive bins in that sample")
)
writeLines(diag_lines, file.path(out_dir, "diagnostics.txt"))
message("Saved: diagnostics.txt")

# =============================================================
# 9. Assign cells to 100 µm grid bins (per-sample origin)
# =============================================================

message("\nAssigning cells to ", GRID_SIZE_UM, " µm bins...")

cell_df <- cell_df %>%
  group_by(sample) %>%
  mutate(
    grid_x  = as.integer(floor((sdimx - min(sdimx)) / GRID_SIZE_UM)),
    grid_y  = as.integer(floor((sdimy - min(sdimy)) / GRID_SIZE_UM)),
    grid_id = paste0(grid_x, "_", grid_y)
  ) %>%
  ungroup()

cat("\nOccupied bins per sample:\n")
print(cell_df %>%
  group_by(sample) %>%
  summarise(n_bins = n_distinct(grid_id), .groups = "drop"))

write.csv(
  cell_df %>% select(cell_id, sample, annotation,
                     sdimx, sdimy, grid_x, grid_y, grid_id,
                     is_immune, is_ifn),
  file.path(out_dir, "grid_per_cell.csv"),
  row.names = FALSE
)
message("Saved: grid_per_cell.csv")

# =============================================================
# 10. Per-bin counts
# =============================================================

message("Computing per-bin counts...")

bin_df <- cell_df %>%
  group_by(sample, grid_id, grid_x, grid_y) %>%
  summarise(
    n_total  = n(),
    n_immune = sum(is_immune),
    n_ifn    = sum(is_ifn),
    .groups  = "drop"
  )

cat("\nTotal bins per sample:\n")
print(table(bin_df$sample))

write.csv(bin_df, file.path(out_dir, "grid_per_bin.csv"),
          row.names = FALSE)
message("Saved: grid_per_bin.csv")

# =============================================================
# 11. Positive-zone summary
# =============================================================

message("Computing positive-zone summaries...")

pos_summary <- bin_df %>%
  group_by(sample) %>%
  summarise(
    n_bins_total            = n(),
    n_immune_pos_zones      = sum(n_immune > 0),
    n_ifn_pos_zones         = sum(n_ifn > 0),
    mean_immune_in_imm_pos  = round(mean(n_immune[n_immune > 0]), 3),
    mean_ifn_in_ifn_pos     = round(mean(n_ifn[n_ifn > 0]), 3),
    max_immune              = max(n_immune),
    max_ifn                 = max(n_ifn),
    .groups = "drop"
  ) %>%
  mutate(
    pct_bins_immune_pos = round(100 * n_immune_pos_zones / n_bins_total, 2),
    pct_bins_ifn_pos    = round(100 * n_ifn_pos_zones    / n_bins_total, 2)
  )

cat("\nPositive-zone summary:\n")
print(as.data.frame(pos_summary), row.names = FALSE)
write.csv(pos_summary,
          file.path(out_dir, "grid_positive_zones_summary.csv"),
          row.names = FALSE)
message("Saved: grid_positive_zones_summary.csv")

# =============================================================
# 12. Enriched-zone summary
# =============================================================

message("Computing enriched-zone tables...")

bin_enriched <- bin_df %>%
  group_by(sample) %>%
  mutate(
    mean_immune_pos = mean(n_immune[n_immune > 0]),
    mean_ifn_pos    = mean(n_ifn[n_ifn > 0]),
    immune_pos      = n_immune > 0,
    ifn_pos         = n_ifn    > 0,
    immune_enrich   = immune_pos & (n_immune > mean_immune_pos),
    ifn_enrich      = ifn_pos    & (n_ifn    > mean_ifn_pos)
  ) %>%
  ungroup()

enrich_summary <- bin_enriched %>%
  group_by(sample) %>%
  summarise(
    n_immune_pos_zones         = sum(immune_pos),
    n_ifn_pos_zones            = sum(ifn_pos),
    n_immune_enrich_zones      = sum(immune_enrich),
    n_ifn_enrich_zones         = sum(ifn_enrich),
    n_both_enrich              = sum(immune_enrich & ifn_enrich),
    pct_immune_enrich_also_ifn = round(
      100 * sum(immune_enrich & ifn_enrich) /
        max(sum(immune_enrich), 1), 2),
    pct_ifn_enrich_also_immune = round(
      100 * sum(immune_enrich & ifn_enrich) /
        max(sum(ifn_enrich), 1), 2),
    ratio_immune_enrich_vs_pos = round(
      sum(immune_enrich) / max(sum(immune_pos), 1), 4),
    ratio_ifn_enrich_vs_pos    = round(
      sum(ifn_enrich)    / max(sum(ifn_pos),    1), 4),
    .groups = "drop"
  )

cat("\nEnriched-zone summary:\n")
print(as.data.frame(enrich_summary), row.names = FALSE)
write.csv(enrich_summary,
          file.path(out_dir, "grid_enriched_zones_summary.csv"),
          row.names = FALSE)
message("Saved: grid_enriched_zones_summary.csv")

write.csv(
  bin_enriched %>%
    select(sample, grid_id, grid_x, grid_y,
           n_total, n_immune, n_ifn,
           immune_pos, ifn_pos, immune_enrich, ifn_enrich),
  file.path(out_dir, "grid_per_bin_enrichment_flags.csv"),
  row.names = FALSE
)
message("Saved: grid_per_bin_enrichment_flags.csv")

# =============================================================
# 13. Helper figures
# =============================================================

message("\nGenerating helper figures...")

# Shared fill scales
immune_scale <- scale_fill_gradientn(
  name     = paste0("Immune\n(", IMMUNE_LABEL, ")\ncount"),
  colors   = c("grey96", "#fee8c8", "#fdbb84", "#e34a33", "#8b0000"),
  na.value = "grey96",
  guide    = guide_colorbar(barwidth = 0.7, barheight = 4)
)
ifn_scale <- scale_fill_gradientn(
  name     = paste0("IFN\n(", IFN_LABEL, ")\ncount"),
  colors   = c("grey96", "#edf8b1", "#7fcdbb", "#1d91c0", "#081d58"),
  na.value = "grey96",
  guide    = guide_colorbar(barwidth = 0.7, barheight = 4)
)

for (samp in SAMPLE_ORDER) {

  s_label  <- sample_labels[samp]
  cells_s  <- cell_df      %>% filter(sample == samp)
  bins_s   <- bin_enriched %>% filter(sample == samp)

  if (nrow(cells_s) == 0) {
    message("  Skipping ", samp, " (no cells)")
    next
  }

  min_x <- min(cells_s$sdimx)
  min_y <- min(cells_s$sdimy)

  bins_plot <- bins_s %>%
    mutate(
      tile_x = min_x + (grid_x + 0.5) * GRID_SIZE_UM,
      tile_y = min_y + (grid_y + 0.5) * GRID_SIZE_UM
    )

  # -----------------------------------------------------------
  # 13a. Grid overlay with cell-type colours — all samples
  # -----------------------------------------------------------

    all_anno      <- sort(unique(cell_df$annotation))
    palette_sub   <- PALETTE[all_anno]
    palette_sub[is.na(palette_sub)] <- "grey70"
    names(palette_sub) <- all_anno

    x_breaks <- seq(min_x, max(cells_s$sdimx) + GRID_SIZE_UM, GRID_SIZE_UM)
    y_breaks <- seq(min_y, max(cells_s$sdimy) + GRID_SIZE_UM, GRID_SIZE_UM)

    p_overlay <- ggplot(cells_s,
        aes(x = sdimx, y = sdimy, color = annotation)) +
      geom_point(size = 0.25, alpha = 0.6, shape = 16) +
      geom_vline(xintercept = x_breaks, color = "grey35",
                 linewidth = 0.18, alpha = 0.45) +
      geom_hline(yintercept = y_breaks, color = "grey35",
                 linewidth = 0.18, alpha = 0.45) +
      scale_color_manual(values = palette_sub, name = "Annotation") +
      coord_equal() +
      labs(
        title    = paste0("100 µm grid overlay — ", s_label),
        subtitle = paste0(
          "Grid = ", GRID_SIZE_UM, " µm × ", GRID_SIZE_UM, " µm | ",
          nrow(cells_s), " cells | ",
          n_distinct(cells_s$grid_id), " occupied bins\n",
          "Red = ", IMMUNE_LABEL,
          " | Brown = ", IFN_LABEL
        ),
        x = "X (µm)", y = "Y (µm)"
      ) +
      theme_pub() +
      theme(
        axis.text       = element_text(size = 7),
        legend.text     = element_text(size = 7),
        legend.key.size = unit(0.32, "cm")
      ) +
      guides(color = guide_legend(
        override.aes = list(size = 2, alpha = 1), ncol = 1))

    save_fig(p_overlay,
             paste0("fig_grid_overlay_", samp),
             width = 13, height = 9)

  # -----------------------------------------------------------
  # 13b. Heatmap: Immune count per bin
  # -----------------------------------------------------------
  n_imm_pos  <- sum(bins_s$immune_pos)
  mean_imm   <- if (n_imm_pos > 0)
                  round(mean(bins_s$n_immune[bins_s$immune_pos]), 2)
                else NA_real_

  p_immune <- ggplot(bins_plot,
      aes(x = tile_x, y = tile_y, fill = n_immune)) +
    geom_tile(width  = GRID_SIZE_UM, height = GRID_SIZE_UM,
              color  = "grey80", linewidth = 0.08) +
    immune_scale +
    coord_equal() +
    labs(
      title    = paste0("Immune (Acod1) density — ", s_label),
      subtitle = paste0(
        GRID_SIZE_UM, " µm bins | Immune-positive bins: ", n_imm_pos,
        " | mean in positive bins: ", mean_imm
      ),
      x = "X (µm)", y = "Y (µm)"
    ) +
    theme_pub() +
    theme(axis.text = element_text(size = 7))

  save_fig(p_immune,
           paste0("fig_grid_immune_", samp),
           width = 8, height = 6.5)

  # -----------------------------------------------------------
  # 13c. Heatmap: IFN count per bin
  # -----------------------------------------------------------
  n_ifn_pos  <- sum(bins_s$ifn_pos)
  mean_ifn   <- if (n_ifn_pos > 0)
                  round(mean(bins_s$n_ifn[bins_s$ifn_pos]), 2)
                else NA_real_

  p_ifn <- ggplot(bins_plot,
      aes(x = tile_x, y = tile_y, fill = n_ifn)) +
    geom_tile(width  = GRID_SIZE_UM, height = GRID_SIZE_UM,
              color  = "grey80", linewidth = 0.08) +
    ifn_scale +
    coord_equal() +
    labs(
      title    = paste0("IFN responsive (Ifit1) density — ", s_label),
      subtitle = paste0(
        GRID_SIZE_UM, " µm bins | IFN-positive bins: ", n_ifn_pos,
        " | mean in positive bins: ", mean_ifn
      ),
      x = "X (µm)", y = "Y (µm)"
    ) +
    theme_pub() +
    theme(axis.text = element_text(size = 7))

  save_fig(p_ifn,
           paste0("fig_grid_ifn_", samp),
           width = 8, height = 6.5)
}

# =============================================================
# 14. Faceted heatmaps with SHARED color scale across samples
# =============================================================
# These complement the per-sample heatmaps above. A single fill
# scale is fixed from 0 to the global maximum across all samples,
# so visual intensity is directly comparable between conditions.
# =============================================================

message("\nGenerating faceted shared-scale heatmaps...")

# Build a combined plotting table across all samples.
# For each sample we recompute the tile centroid in µm using the
# per-sample grid origin (same logic as the per-sample loop above).
all_bins_plot <- lapply(SAMPLE_ORDER, function(samp) {
  cells_s <- cell_df %>% filter(sample == samp)
  if (nrow(cells_s) == 0) return(NULL)
  min_x <- min(cells_s$sdimx)
  min_y <- min(cells_s$sdimy)
  bin_enriched %>%
    filter(sample == samp) %>%
    mutate(
      tile_x        = min_x + (grid_x + 0.5) * GRID_SIZE_UM,
      tile_y        = min_y + (grid_y + 0.5) * GRID_SIZE_UM,
      sample_label  = factor(sample_labels[as.character(sample)],
                             levels = sample_labels[SAMPLE_ORDER])
    )
}) %>%
  bind_rows()

# Global maxima (used as shared upper limits)
global_max_immune <- max(all_bins_plot$n_immune, na.rm = TRUE)
global_max_ifn    <- max(all_bins_plot$n_ifn,    na.rm = TRUE)

message("  Global max n_immune across samples: ", global_max_immune)
message("  Global max n_ifn    across samples: ", global_max_ifn)

# -----------------------------------------------------------
# 14a. Faceted Immune heatmap — shared scale
# -----------------------------------------------------------

p_imm_facet <- ggplot(all_bins_plot,
    aes(x = tile_x, y = tile_y, fill = n_immune)) +
  geom_tile(width  = GRID_SIZE_UM, height = GRID_SIZE_UM,
            color  = "grey85", linewidth = 0.07) +
  scale_fill_gradientn(
    name   = paste0("Immune\n(Acod1)\ncount"),
    colors = c("grey96", "#fee8c8", "#fdbb84", "#e34a33", "#8b0000"),
    limits = c(0, global_max_immune),
    na.value = "grey96",
    guide  = guide_colorbar(barwidth = 0.7, barheight = 5)
  ) +
  facet_wrap(~ sample_label, nrow = 2) +
  coord_equal() +
  labs(
    title    = "Immune (Acod1) density — shared scale across samples",
    subtitle = paste0(
      GRID_SIZE_UM, " µm bins | fill scale: 0 – ", global_max_immune,
      " cells/bin (global maximum) | scale is identical across all panels"
    ),
    x = "X (µm)", y = "Y (µm)"
  ) +
  theme_pub() +
  theme(
    axis.text        = element_text(size = 6.5),
    strip.text       = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.spacing    = unit(0.6, "cm")
  )

save_fig(p_imm_facet,
         "fig_grid_immune_faceted_shared_scale",
         width = 11, height = 11)

# -----------------------------------------------------------
# 14b. Faceted IFN heatmap — shared scale
# -----------------------------------------------------------

p_ifn_facet <- ggplot(all_bins_plot,
    aes(x = tile_x, y = tile_y, fill = n_ifn)) +
  geom_tile(width  = GRID_SIZE_UM, height = GRID_SIZE_UM,
            color  = "grey85", linewidth = 0.07) +
  scale_fill_gradientn(
    name   = paste0("IFN\n(Ifit1)\ncount"),
    colors = c("grey96", "#edf8b1", "#7fcdbb", "#1d91c0", "#081d58"),
    limits = c(0, global_max_ifn),
    na.value = "grey96",
    guide  = guide_colorbar(barwidth = 0.7, barheight = 5)
  ) +
  facet_wrap(~ sample_label, nrow = 2) +
  coord_equal() +
  labs(
    title    = "IFN responsive (Ifit1) density — shared scale across samples",
    subtitle = paste0(
      GRID_SIZE_UM, " µm bins | fill scale: 0 – ", global_max_ifn,
      " cells/bin (global maximum) | scale is identical across all panels"
    ),
    x = "X (µm)", y = "Y (µm)"
  ) +
  theme_pub() +
  theme(
    axis.text        = element_text(size = 6.5),
    strip.text       = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.spacing    = unit(0.6, "cm")
  )

save_fig(p_ifn_facet,
         "fig_grid_ifn_faceted_shared_scale",
         width = 11, height = 11)

# =============================================================
# Done
# =============================================================

message("\n=== Done. Outputs in: ", out_dir, " ===\n")
message("Files produced:")
message("  diagnostics.txt")
message("  grid_per_cell.csv")
message("  grid_per_bin.csv")
message("  grid_per_bin_enrichment_flags.csv")
message("  grid_positive_zones_summary.csv")
message("  grid_enriched_zones_summary.csv")
message("  fig_grid_overlay_LCMV_1wpi.pdf/jpg")
message("  fig_grid_immune_<sample>.pdf/jpg  (x", length(SAMPLE_ORDER), ")")
message("  fig_grid_ifn_<sample>.pdf/jpg     (x", length(SAMPLE_ORDER), ")")
message("  fig_grid_immune_faceted_shared_scale.pdf/jpg")
message("  fig_grid_ifn_faceted_shared_scale.pdf/jpg")

# =============================================================
# 15. Global-threshold enrichment (pooled mean across all samples)
# =============================================================
# Complementary to section 12 (per-sample thresholds).
# Here the enrichment threshold is derived from all samples pooled:
#   global_mean_immune_pos = mean(n_immune) across ALL Immune-positive
#                            bins from every sample combined
#   global_mean_ifn_pos    = mean(n_ifn)    across ALL IFN-positive
#                            bins from every sample combined
#
# New flags:
#   immune_enrich_global = n_immune > global_mean_immune_pos
#   ifn_enrich_global    = n_ifn    > global_mean_ifn_pos
#
# Existing per-sample enrichment flags are not modified.
# =============================================================

message("\n--- Section 15: global-threshold enrichment ---")

# ------------------------------------------------------------------
# 15-1. Compute pooled global means
# ------------------------------------------------------------------

global_mean_immune_pos <- mean(bin_df$n_immune[bin_df$n_immune > 0],
                               na.rm = TRUE)
global_mean_ifn_pos    <- mean(bin_df$n_ifn[bin_df$n_ifn > 0],
                               na.rm = TRUE)

message("  Pooled global mean Immune (positive bins): ",
        round(global_mean_immune_pos, 3))
message("  Pooled global mean IFN    (positive bins): ",
        round(global_mean_ifn_pos, 3))

# ------------------------------------------------------------------
# 15-2. Add global enrichment flags to bin table
# ------------------------------------------------------------------

bin_global <- bin_enriched %>%           # reuse existing flags too
  mutate(
    immune_enrich_global = n_immune > global_mean_immune_pos,
    ifn_enrich_global    = n_ifn    > global_mean_ifn_pos
  )

# ------------------------------------------------------------------
# 15-3. Save per-bin table with global flags
# ------------------------------------------------------------------

write.csv(
  bin_global %>%
    select(sample, grid_id, grid_x, grid_y,
           n_total, n_immune, n_ifn,
           immune_pos, ifn_pos,
           immune_enrich_global, ifn_enrich_global),
  file.path(out_dir, "grid_per_bin_global_threshold_flags.csv"),
  row.names = FALSE
)
message("  Saved: grid_per_bin_global_threshold_flags.csv")

# ------------------------------------------------------------------
# 15-4. Per-sample summary with global thresholds
# ------------------------------------------------------------------

global_threshold_summary <- bin_global %>%
  group_by(sample) %>%
  summarise(
    n_bins_total                          = n(),
    n_immune_pos_zones                    = sum(immune_pos),
    n_ifn_pos_zones                       = sum(ifn_pos),
    n_immune_enrich_global                = sum(immune_enrich_global),
    n_ifn_enrich_global                   = sum(ifn_enrich_global),
    pct_immune_enrich_global_among_immune_pos =
      round(100 * sum(immune_enrich_global & immune_pos) /
              max(sum(immune_pos), 1), 2),
    pct_ifn_enrich_global_among_ifn_pos =
      round(100 * sum(ifn_enrich_global & ifn_pos) /
              max(sum(ifn_pos), 1), 2),
    .groups = "drop"
  )

cat("\nGlobal-threshold summary by sample:\n")
print(as.data.frame(global_threshold_summary), row.names = FALSE)
write.csv(global_threshold_summary,
          file.path(out_dir, "grid_global_threshold_summary.csv"),
          row.names = FALSE)
message("  Saved: grid_global_threshold_summary.csv")
