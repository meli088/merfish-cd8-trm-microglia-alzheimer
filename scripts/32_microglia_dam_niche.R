#!/usr/bin/env Rscript
# =============================================================
# Script: 32_microglia_dam_niche.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-05
#
# Goal:
#   Compare microglia in vs out of the immune niche:
#   - DAM module score (AddModuleScore) on combined microglia population
#   - DEG in niche vs out niche per timepoint
#   - DEG across spatial distance bins from immune niche
#
# "In niche"  = Microglia (C1qa)  — from object 08 (immune niche)
# "Out niche" = Microglia (P2ry12) — from object 04 (global BANKSY)
#
# Inputs:
#   objects/04_banksy_joint_lam08_after_bloc3.rds
#   objects/08_immune_annotated_lam02_res03.rds
#   ncells_by_sample_lam02_res09_joint_long.csv
#   outputs/banksy/dam_signature_curation/Upregulated_DAM.csv
#   outputs/banksy/dam_signature_curation/Downregulated_DAM.csv
#   outputs/banksy/dam_signature_curation/Microglia_signature_union.csv
#   outputs/banksy/nearest_immune_distance/nearest_immune_distance_per_cell_lam02_res09.csv
#
# Outputs (outputs/banksy/microglia_dam_niche/):
#   [Sec 2] module_scores_per_cell.csv
#           fig_module_score_[sig]_violin.pdf/jpg
#           fig_module_score_[sig]_spatial.pdf/jpg
#   [Sec 3] DEG_microglia_in_vs_out_[timepoint].csv
#           fig_volcano_in_vs_out_[timepoint].pdf/jpg
#   [Sec 4] DEG_by_distance_[timepoint].csv
#           fig_ndeg_by_distance_[timepoint].pdf/jpg
#           fig_volcano_distance_max_[timepoint].pdf/jpg
# =============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(Banksy)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(Cairo)
})

set.seed(1997)

base_path <- normalizePath(".")   # Run from project root
setwd(base_path)

# Palette et ordre des annotations (fichier partagé)
source("scripts/00_palette.R")

# =============================================================
# Global parameters
# =============================================================

SAMPLE_ORDER           <- c("mock_6wpi", "LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")
MICROGLIA_GLOBAL_LABEL <- "Microglia (P2ry12)"   # out niche — global object
MICROGLIA_NICHE_LABEL  <- "Microglia (C1qa)"      # in niche  — object 08
LAM                    <- 0.2
RES_TARGET             <- 0.9
DISTANCE_THRESHOLDS    <- c(100, 200, 300)          # seuils µm pour Section 4
FDR_CUTOFF             <- 0.05
FC_CUTOFF              <- 0.25
TOP_N_LABEL            <- 15

TIMEPOINTS <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")

out_dir <- "outputs/banksy/microglia_dam_niche"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

SAMPLE_LABELS <- c(
  mock_6wpi = "Mock 6 wpi",
  LCMV_1wpi = "LCMV 1 wpi",
  LCMV_3wpi = "LCMV 3 wpi",
  LCMV_6wpi = "LCMV 6 wpi"
)

COND_PALETTE <- c(
  "Mock 6 wpi" = "grey60",
  "LCMV 1 wpi" = "#56B4E9",
  "LCMV 3 wpi" = "#E69F00",
  "LCMV 6 wpi" = "#D55E00"
)

direction_colors <- c(
  "up"   = "#B2182B",
  "down" = "#2166AC",
  "ns"   = "grey75"
)

niche_palette <- c(
  "Out niche" = "grey60",
  "In niche"  = "#B2182B"
)

# FindMarkers settings
MIN_PCT   <- 0.05
FC_THRESH <- 0.1
TEST_USE  <- "wilcox"

# =============================================================
# Helper: find BANKSY cluster column (from scripts 28/29/31)
# =============================================================

find_cl_col <- function(se, lam, res) {
  all_cols <- clusterNames(se)
  lam_str  <- gsub("\\.", "\\\\.", as.character(lam))
  lam_cols <- all_cols[grep(paste0("lam", lam_str), all_cols)]
  if (length(lam_cols) == 0) return(NULL)
  res_vals <- suppressWarnings(as.numeric(sub(".*_res", "", lam_cols)))
  idx      <- which(!is.na(res_vals) & abs(res_vals - res) < 1e-8)
  if (length(idx) == 0) return(NULL)
  lam_cols[idx[1]]
}

# Helper: save PDF + JPG with cairo_pdf
save_fig <- function(p, fname, width, height) {
  CairoPDF(file.path(out_dir, paste0(fname, ".pdf")), width = width, height = height)
  print(p)
  dev.off()
  CairoJPEG(file.path(out_dir, paste0(fname, ".jpg")),
            width = width * 150, height = height * 150, res = 150)
  print(p)
  dev.off()
  message("  Saved: ", fname, ".pdf / .jpg")
}

# =============================================================
# SECTION 1 — Load objects and reconstruct annotations
# =============================================================

message("\n=== SECTION 1: Loading objects ===\n")

# --- 1a. Global BANKSY object ---
obj_file <- file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds")
stopifnot("Global object not found" = file.exists(obj_file))
message("Loading: ", obj_file)
se_global <- readRDS(obj_file)
message("  ", ncol(se_global), " cells | class: ", class(se_global)[1])

cl_col <- find_cl_col(se_global, LAM, RES_TARGET)
if (is.null(cl_col)) {
  stop("BANKSY cluster column not found for lambda=", LAM, " res=", RES_TARGET,
       "\nAvailable: ", paste(clusterNames(se_global), collapse = ", "))
}
message("  Cluster column: ", cl_col)

# --- 1b. Reconstruct annotation from CSV (exact pattern from scripts 28/29/31) ---
csv_path <- "ncells_by_sample_lam02_res09_joint_long.csv"
stopifnot("Annotation CSV not found" = file.exists(csv_path))

annot_long <- read_delim(
  csv_path,
  delim          = ";",
  locale         = locale(decimal_mark = "."),
  show_col_types = FALSE,
  trim_ws        = TRUE
) %>%
  select(-matches("^Unnamed")) %>%
  mutate(
    banksy_domain = as.character(banksy_domain),
    annotation    = trimws(as.character(annotation))
  )

annotation_map <- annot_long %>%
  filter(!is.na(annotation), annotation != "") %>%
  distinct(banksy_domain, annotation)

message("  Annotation mappings: ", nrow(annotation_map))

domain_labels <- paste0("Domain_", as.character(colData(se_global)[[cl_col]]))
anno_lookup   <- setNames(annotation_map$annotation, annotation_map$banksy_domain)

annotation_global <- ifelse(
  !is.na(anno_lookup[domain_labels]) & anno_lookup[domain_labels] != "",
  anno_lookup[domain_labels],
  "Non annote"
)

cat("\nAnnotation distribution (global, all cells):\n")
print(sort(table(annotation_global), decreasing = TRUE))

if (!MICROGLIA_GLOBAL_LABEL %in% annotation_global) {
  stop("Label '", MICROGLIA_GLOBAL_LABEL, "' not found in global annotations.\n",
       "Present: ", paste(sort(unique(annotation_global)), collapse = ", "))
}

# --- 1c. Extract Out-niche microglia (P2ry12) from global object ---
out_idx <- which(annotation_global == MICROGLIA_GLOBAL_LABEL)
se_out  <- se_global[, out_idx]
colData(se_out)$niche_status <- "Out niche"
message("\nOut-niche [", MICROGLIA_GLOBAL_LABEL, "]: ", ncol(se_out), " cells")
cat("  Per sample:\n")
print(table(as.character(colData(se_out)$sample)))

# --- 1d. Niche object (08) and extract In-niche microglia (C1qa) ---
niche_file <- file.path("objects", "08_immune_annotated_lam02_res03.rds")
stopifnot("Niche object not found" = file.exists(niche_file))
message("\nLoading: ", niche_file)
se_niche <- readRDS(niche_file)
message("  ", ncol(se_niche), " cells | class: ", class(se_niche)[1])

cell_type_niche <- as.character(colData(se_niche)$cell_type)
if (!MICROGLIA_NICHE_LABEL %in% cell_type_niche) {
  stop("Label '", MICROGLIA_NICHE_LABEL, "' not found in se_niche$cell_type.\n",
       "Present: ", paste(sort(unique(cell_type_niche)), collapse = ", "))
}

in_idx <- which(cell_type_niche == MICROGLIA_NICHE_LABEL)
se_in  <- se_niche[, in_idx]
colData(se_in)$niche_status <- "In niche"
message("\nIn-niche [", MICROGLIA_NICHE_LABEL, "]: ", ncol(se_in), " cells")
cat("  Per sample:\n")
print(table(as.character(colData(se_in)$sample)))

# --- 1e. Summary: n cells per niche status × sample ---
cat("\n--- Cell count summary (in / out niche) ---\n")
for (s in SAMPLE_ORDER) {
  n_in  <- sum(as.character(colData(se_in)$sample)  == s, na.rm = TRUE)
  n_out <- sum(as.character(colData(se_out)$sample) == s, na.rm = TRUE)
  cat(sprintf("  %-12s  In niche: %4d   Out niche: %4d\n", s, n_in, n_out))
}

# =============================================================
# SECTION 2 — AddModuleScore: DAM signature per cell
# =============================================================

message("\n=== SECTION 2: AddModuleScore ===\n")

# --- 2a. Read and case-harmonise DAM signatures ---
sig_dir <- file.path("outputs", "banksy", "dam_signature_curation")

read_sig <- function(fname) {
  f <- file.path(sig_dir, fname)
  stopifnot("Signature file not found" = file.exists(f))
  read.csv(f, stringsAsFactors = FALSE)$gene
}

sig_raw <- list(
  "Upregulated_DAM"           = read_sig("Upregulated_DAM.csv"),
  "Downregulated_DAM"         = read_sig("Downregulated_DAM.csv"),
  "Microglia_signature_union" = read_sig("Microglia_signature_union.csv")
)
for (nm in names(sig_raw)) {
  message("  ", nm, ": ", length(sig_raw[[nm]]), " genes (raw)")
}

# Panel uses mouse Titlecase; signature files use human ALLCAPS
panel_genes    <- rownames(se_out)
panel_title    <- tools::toTitleCase(tolower(panel_genes))
panel_name_map <- setNames(panel_genes, panel_title)

harm_intersect <- function(sig_genes) {
  sig_title <- tools::toTitleCase(tolower(sig_genes))
  hits      <- sig_title[sig_title %in% panel_title]
  matched   <- panel_name_map[hits]
  unique(as.character(matched[!is.na(matched)]))
}

sig_filtered <- lapply(sig_raw, harm_intersect)
for (nm in names(sig_filtered)) {
  message("  ", nm, ": ", length(sig_filtered[[nm]]),
          " panel genes (from ", length(sig_raw[[nm]]), " raw)")
}

sig_modules <- sig_filtered[sapply(sig_filtered, length) > 0]
if (length(sig_modules) == 0) stop("No signature genes matched the panel.")

# --- 2b. Extract XY coords before Seurat conversion ---
coords_out <- as.data.frame(spatialCoords(se_out))
colnames(coords_out) <- c("x", "y")
coords_out$orig_id <- colnames(se_out)
coords_out$sample  <- as.character(colData(se_out)$sample)

coords_in <- as.data.frame(spatialCoords(se_in))
colnames(coords_in) <- c("x", "y")
coords_in$orig_id <- colnames(se_in)
coords_in$sample  <- as.character(colData(se_in)$sample)

# --- 2c. Convert each to Seurat and normalise ---
assay_out <- if ("counts" %in% assayNames(se_out)) "counts" else assayNames(se_out)[1]
assay_in  <- if ("counts" %in% assayNames(se_in))  "counts" else assayNames(se_in)[1]

message("Converting Out-niche to Seurat...")
so_out <- suppressWarnings(as.Seurat(se_out, counts = assay_out, data = NULL))
dassay_out <- if ("RNA" %in% SeuratObject::Assays(so_out)) "RNA" else SeuratObject::Assays(so_out)[1]
DefaultAssay(so_out) <- dassay_out
so_out <- NormalizeData(so_out, assay = dassay_out, verbose = FALSE)
so_out$niche_status <- "Out niche"
so_out$sample       <- as.character(colData(se_out)$sample)
message("  ", ncol(so_out), " cells")

message("Converting In-niche to Seurat...")
so_in <- suppressWarnings(as.Seurat(se_in, counts = assay_in, data = NULL))
dassay_in <- if ("RNA" %in% SeuratObject::Assays(so_in)) "RNA" else SeuratObject::Assays(so_in)[1]
DefaultAssay(so_in) <- dassay_in
so_in <- NormalizeData(so_in, assay = dassay_in, verbose = FALSE)
so_in$niche_status <- "In niche"
so_in$sample       <- as.character(colData(se_in)$sample)
message("  ", ncol(so_in), " cells")

# Restrict to common genes
common_genes <- intersect(rownames(so_out), rownames(so_in))
message("  Common genes: ", length(common_genes))
so_out_sub <- subset(so_out, features = common_genes)
so_in_sub  <- subset(so_in,  features = common_genes)

# --- 2d. Merge and re-normalise ---
message("Merging Seurat objects...")
so_combined <- merge(
  so_out_sub, so_in_sub,
  add.cell.ids = c("out", "in"),
  project      = "microglia_niche"
)

dassay_comb <- if ("RNA" %in% SeuratObject::Assays(so_combined)) "RNA" else SeuratObject::Assays(so_combined)[1]
DefaultAssay(so_combined) <- dassay_comb
so_combined <- NormalizeData(so_combined, assay = dassay_comb, verbose = FALSE)
message("  Combined: ", ncol(so_combined), " cells")

# Build coord lookup keyed by merged cell name ("out_<orig>" / "in_<orig>")
coords_combined <- bind_rows(
  data.frame(
    x            = coords_out$x,
    y            = coords_out$y,
    merged_id    = paste0("out_", coords_out$orig_id),
    sample       = coords_out$sample,
    niche_status = "Out niche",
    stringsAsFactors = FALSE
  ),
  data.frame(
    x            = coords_in$x,
    y            = coords_in$y,
    merged_id    = paste0("in_", coords_in$orig_id),
    sample       = coords_in$sample,
    niche_status = "In niche",
    stringsAsFactors = FALSE
  )
)

# --- 2e. AddModuleScore for each signature ---
# ctrl must be < (n_panel_genes - n_sig_genes); cap at 5 for small MERFISH panels
n_panel <- nrow(so_combined)
ctrl_n  <- max(1L, min(5L, floor((n_panel - max(sapply(sig_modules, length))) / 2)))
message("Running AddModuleScore (ctrl=", ctrl_n, ")...")
for (nm in names(sig_modules)) {
  message("  ", nm, " (", length(sig_modules[[nm]]), " genes)")
  set.seed(1997)
  so_combined <- AddModuleScore(
    so_combined,
    features = list(sig_modules[[nm]]),
    name     = nm,
    ctrl     = ctrl_n,
    seed     = 1997
  )
}

# Collect scores into a data frame (Seurat appends "1" to the name)
scores_df <- data.frame(
  cell_id      = colnames(so_combined),
  sample       = so_combined$sample,
  niche_status = so_combined$niche_status,
  stringsAsFactors = FALSE
)
for (nm in names(sig_modules)) {
  score_col <- paste0(nm, "1")
  scores_df[[nm]] <- so_combined[[score_col]][, 1]
}
scores_df$sample_label <- factor(
  SAMPLE_LABELS[scores_df$sample],
  levels = SAMPLE_LABELS[SAMPLE_ORDER]
)
scores_df$niche_status <- factor(scores_df$niche_status,
                                  levels = c("Out niche", "In niche"))

write.csv(scores_df, file.path(out_dir, "module_scores_per_cell.csv"), row.names = FALSE)
message("Saved: module_scores_per_cell.csv")

# --- 2f. Violin + spatial plots per signature ---
for (nm in names(sig_modules)) {
  sig_slug    <- gsub("[^a-z0-9]", "_", tolower(nm))
  sig_display <- gsub("_", " ", nm)
  message("\nPlotting: ", nm)

  # Wilcoxon test per sample (computed manually to avoid dplyr scoping issues)
  wilcox_res <- do.call(rbind, lapply(SAMPLE_ORDER, function(s) {
    d <- scores_df[scores_df$sample == s, ]
    d_in  <- d[[nm]][d$niche_status == "In niche"]
    d_out <- d[[nm]][d$niche_status == "Out niche"]
    pval  <- if (length(d_in) < 3 || length(d_out) < 3) {
      NA_real_
    } else {
      tryCatch(
        wilcox.test(d_in, d_out, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )
    }
    data.frame(sample = s, pval = pval, stringsAsFactors = FALSE)
  }))
  wilcox_res$pval_label <- with(wilcox_res, case_when(
    is.na(pval)  ~ "n.s.",
    pval < 0.001 ~ "***",
    pval < 0.01  ~ "**",
    pval < 0.05  ~ "*",
    TRUE         ~ "n.s."
  ))
  wilcox_res$sample_label <- factor(
    SAMPLE_LABELS[wilcox_res$sample],
    levels = SAMPLE_LABELS[SAMPLE_ORDER]
  )

  # y_pos for Wilcoxon annotation: above the highest violin per facet
  y_pos_df <- do.call(rbind, lapply(SAMPLE_ORDER, function(s) {
    vals <- scores_df[[nm]][scores_df$sample == s]
    y_max <- if (length(vals) == 0) 0 else max(vals, na.rm = TRUE)
    y_rng <- if (length(vals) < 2)  0.1 else diff(range(vals, na.rm = TRUE))
    data.frame(
      sample_label = factor(SAMPLE_LABELS[s], levels = SAMPLE_LABELS[SAMPLE_ORDER]),
      y_pos        = y_max + 0.1 * y_rng,
      stringsAsFactors = FALSE
    )
  }))
  y_pos_df <- merge(y_pos_df,
                    wilcox_res[, c("sample_label", "pval_label")],
                    by = "sample_label", all.x = TRUE)

  # ---- Violin plot ----
  scores_plot  <- scores_df %>% filter(!is.na(sample_label))
  # geom_violin requires >= 2 data points per group; filter separately to avoid warnings
  violin_data  <- scores_plot %>%
    group_by(sample_label, niche_status) %>%
    filter(n() >= 2) %>%
    ungroup()

  p_violin <- ggplot(
    scores_plot,
    aes(x = niche_status, y = .data[[nm]], fill = niche_status)
  ) +
    geom_violin(data = violin_data,
                trim = TRUE, scale = "width", alpha = 0.5, colour = NA) +
    geom_boxplot(width = 0.18, outlier.size = 0.4, outlier.alpha = 0.4,
                 colour = "grey20", fill = "white") +
    geom_text(
      data        = y_pos_df,
      aes(x       = 1.5, y = y_pos, label = pval_label),
      inherit.aes = FALSE,
      size        = 3.5, fontface = "bold", colour = "grey20"
    ) +
    facet_wrap(~ sample_label, nrow = 1) +
    scale_fill_manual(values = niche_palette, guide = "none") +
    labs(
      title    = paste0("Module score: ", sig_display),
      subtitle = paste0(length(sig_modules[[nm]]),
                        " panel genes | Wilcoxon per condition"),
      x        = NULL,
      y        = "Module score (AddModuleScore)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title         = element_text(face = "bold", size = 11, hjust = 0),
      plot.subtitle      = element_text(size = 8, colour = "grey40", hjust = 0),
      strip.background   = element_rect(fill = "grey95"),
      strip.text         = element_text(size = 9),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank()
    )

  save_fig(p_violin, paste0("fig_module_score_", sig_slug, "_violin"),
           width = 9, height = 4.5)

  # ---- Spatial feature plot ----
  # Join scores back to XY coordinates via merged cell ID
  score_lookup <- data.frame(
    merged_id = scores_df$cell_id,
    score     = scores_df[[nm]],
    stringsAsFactors = FALSE
  )
  spatial_df <- coords_combined %>%
    left_join(score_lookup, by = "merged_id") %>%
    filter(!is.na(sample), sample %in% SAMPLE_ORDER) %>%
    mutate(
      sample_label = factor(SAMPLE_LABELS[sample], levels = SAMPLE_LABELS[SAMPLE_ORDER])
    )

  score_range <- quantile(spatial_df$score, c(0.01, 0.99), na.rm = TRUE)

  p_spatial <- ggplot(
    spatial_df,
    aes(x      = x,
        y      = y,
        colour = pmin(pmax(score, score_range[1]), score_range[2]))
  ) +
    geom_point(size = 0.3, alpha = 0.7, stroke = 0) +
    facet_wrap(~ sample_label, nrow = 1, scales = "free") +
    scale_colour_gradientn(
      colours = c("#2166AC", "grey90", "#B2182B"),
      limits  = score_range,
      name    = "Score",
      oob     = scales::squish
    ) +
    labs(
      title    = paste0("Spatial — Module score: ", sig_display),
      subtitle = "Microglia (in niche + out niche combined)",
      x        = "X (\u00b5m)",
      y        = "Y (\u00b5m)"
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title        = element_text(face = "bold", size = 10, hjust = 0),
      plot.subtitle     = element_text(size = 8, colour = "grey40", hjust = 0),
      strip.background  = element_rect(fill = "grey95"),
      strip.text        = element_text(size = 8),
      panel.grid        = element_blank(),
      axis.text         = element_text(size = 6),
      legend.key.height = unit(0.5, "cm")
    )

  save_fig(p_spatial, paste0("fig_module_score_", sig_slug, "_spatial"),
           width = 14, height = 4)
}

# =============================================================
# SECTION 3 — DEG: microglia In niche vs Out niche
# =============================================================

message("\n=== SECTION 3: DEG In niche vs Out niche ===\n")

# Set identity to niche_status for FindMarkers
so3 <- SetIdent(so_combined, value = "niche_status")

make_volcano_niche <- function(mk, tp_slug, tp_label) {
  mk$neg_log10_fdr <- -log10(pmax(mk$p_val_adj, .Machine$double.eps))
  x_cap <- max(3.5, quantile(abs(mk$avg_log2FC), 0.999, na.rm = TRUE) * 1.05)
  y_cap <- max(5,   quantile(mk$neg_log10_fdr,   0.999, na.rm = TRUE) * 1.05)
  mk$x_plot <- pmin(pmax(mk$avg_log2FC,   -x_cap), x_cap)
  mk$y_plot <- pmin(mk$neg_log10_fdr, y_cap)

  sig_mk <- mk %>% filter(direction != "ns")
  if (nrow(sig_mk) > 0) {
    lab_genes <- sig_mk %>%
      arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
      slice_head(n = TOP_N_LABEL) %>%
      pull(gene)
  } else {
    lab_genes <- mk %>%
      arrange(p_val_adj) %>%
      slice_head(n = min(TOP_N_LABEL, 10)) %>%
      pull(gene)
  }
  mk$label <- ifelse(mk$gene %in% lab_genes, mk$gene, NA_character_)

  n_up   <- sum(mk$direction == "up",   na.rm = TRUE)
  n_down <- sum(mk$direction == "down", na.rm = TRUE)

  ggplot(mk, aes(x = x_plot, y = y_plot, colour = direction)) +
    geom_point(size = 1.2, alpha = 0.7, stroke = 0) +
    geom_text_repel(
      data          = mk %>% filter(!is.na(label)),
      aes(label     = label),
      size          = 2.6,
      fontface      = "italic",
      max.overlaps  = 20,
      box.padding   = 0.3,
      point.padding = 0.2,
      segment.size  = 0.3,
      segment.color = "grey50",
      show.legend   = FALSE
    ) +
    geom_hline(yintercept = -log10(FDR_CUTOFF),
               linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF),
               linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    scale_colour_manual(
      values = direction_colors,
      labels = c(
        up   = paste0("Up in In niche (n=", n_up,   ")"),
        down = paste0("Down in In niche (n=", n_down, ")"),
        ns   = "Not significant"
      ),
      name = NULL
    ) +
    scale_x_continuous(limits = c(-x_cap, x_cap),
                       expand = expansion(mult = 0.02)) +
    scale_y_continuous(limits = c(0, y_cap),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(
      title    = paste0("Microglia: In niche vs Out niche \u2014 ", tp_label),
      subtitle = paste0(nrow(mk), " genes tested  |  FDR \u2264 ", FDR_CUTOFF,
                        "  |  |log2FC| > ", FC_CUTOFF),
      x        = "log2FC (In niche / Out niche)",
      y        = expression(-log[10](FDR))
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title       = element_text(face = "bold", size = 11, hjust = 0),
      plot.subtitle    = element_text(size = 8, colour = "grey40", hjust = 0),
      legend.position  = "bottom",
      legend.text      = element_text(size = 8),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey93", linewidth = 0.3)
    )
}

for (tp in TIMEPOINTS) {
  tp_slug  <- gsub("[^a-z0-9]", "_", tolower(tp))
  tp_label <- SAMPLE_LABELS[tp]
  n_in  <- sum(so_combined$sample == tp & so_combined$niche_status == "In niche",
               na.rm = TRUE)
  n_out <- sum(so_combined$sample == tp & so_combined$niche_status == "Out niche",
               na.rm = TRUE)
  message("\n--- ", tp, " : In (n=", n_in, ") vs Out (n=", n_out, ") ---")

  if (n_in < 5 || n_out < 5) {
    message("  Skipping: insufficient cells (< 5 in one group).")
    next
  }

  cells_tp <- colnames(so3)[so3$sample == tp]
  so3_tp   <- subset(so3, cells = cells_tp)

  mk <- tryCatch(
    FindMarkers(
      so3_tp,
      ident.1         = "In niche",
      ident.2         = "Out niche",
      only.pos        = FALSE,
      min.pct         = MIN_PCT,
      logfc.threshold = FC_THRESH,
      return.thresh   = 1,
      test.use        = TEST_USE,
      verbose         = FALSE
    ),
    error = function(e) {
      message("  ERROR: ", conditionMessage(e)); NULL
    }
  )

  if (is.null(mk) || nrow(mk) == 0) {
    message("  No results returned — skipping.")
    next
  }

  mk$gene      <- rownames(mk)
  mk$timepoint <- tp
  mk$direction <- case_when(
    mk$p_val_adj <= FDR_CUTOFF & mk$avg_log2FC >  FC_CUTOFF ~ "up",
    mk$p_val_adj <= FDR_CUTOFF & mk$avg_log2FC < -FC_CUTOFF ~ "down",
    TRUE ~ "ns"
  )
  mk <- mk %>%
    arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
    select(gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj, direction, timepoint)

  n_up   <- sum(mk$direction == "up")
  n_down <- sum(mk$direction == "down")
  message("  ", nrow(mk), " genes | ", n_up, " up | ", n_down, " down",
          " (FDR<=", FDR_CUTOFF, ", |log2FC|>", FC_CUTOFF, ")")

  csv_out <- file.path(out_dir, paste0("DEG_microglia_in_vs_out_", tp_slug, ".csv"))
  write.csv(mk, csv_out, row.names = FALSE)
  message("  Saved: DEG_microglia_in_vs_out_", tp_slug, ".csv")

  p <- make_volcano_niche(mk, tp_slug, tp_label)
  save_fig(p, paste0("fig_volcano_in_vs_out_", tp_slug), width = 6, height = 5.5)
}

# =============================================================
# SECTION 4 — DEG microglie (P2ry12) proche vs loin de la niche
#             Seuils : 100, 200, 300 µm
# =============================================================

message("\n=== SECTION 4: DEG by distance threshold ===\n")

# Supprimer d'éventuels fichiers de l'ancienne version (bin pairs)
old_patterns <- c(
  file.path(out_dir, "DEG_by_distance_lcmv_1wpi.csv"),
  file.path(out_dir, "DEG_by_distance_lcmv_6wpi.csv")
)
for (f in old_patterns) {
  if (file.exists(f)) { file.remove(f); message("  Removed old file: ", f) }
}

dist_file <- file.path(
  "outputs", "banksy", "nearest_immune_distance",
  "nearest_immune_distance_per_cell_lam02_res09.csv"
)
stopifnot("Distance CSV not found" = file.exists(dist_file))

dist_df <- read.csv(dist_file, stringsAsFactors = FALSE) %>%
  mutate(annotation = trimws(annotation))
message("  Distance CSV: ", nrow(dist_df), " cells")

# Keep only Microglia (P2ry12)
dist_mg_all <- dist_df %>%
  filter(annotation == MICROGLIA_GLOBAL_LABEL)
message("  Microglia (P2ry12) rows: ", nrow(dist_mg_all))

if (nrow(dist_mg_all) < 50) {
  stop("Too few Microglia rows in distance CSV (", nrow(dist_mg_all),
       "). Check annotation label / trimws.")
}

# Match cell_ids to se_out colnames
matched_all <- dist_mg_all$cell_id[dist_mg_all$cell_id %in% colnames(se_out)]
message("  Matched to se_out: ", length(matched_all), " / ", nrow(dist_mg_all))
if (length(matched_all) < 50) {
  stop("Too few cells matched (", length(matched_all), "). Check cell_id format.")
}

# Seurat conversion of all matched microglia
se_dist_all <- se_out[, matched_all]
assay_dist  <- if ("counts" %in% assayNames(se_dist_all)) "counts" else assayNames(se_dist_all)[1]
so_dist_all <- suppressWarnings(as.Seurat(se_dist_all, counts = assay_dist, data = NULL))
dassay_dist <- if ("RNA" %in% SeuratObject::Assays(so_dist_all)) "RNA" else SeuratObject::Assays(so_dist_all)[1]
DefaultAssay(so_dist_all) <- dassay_dist
so_dist_all <- NormalizeData(so_dist_all, assay = dassay_dist, verbose = FALSE)
so_dist_all$sample <- as.character(colData(se_dist_all)$sample)
# Attach distance values
dist_lookup_all <- setNames(dist_mg_all$nearest_immune_distance_um, dist_mg_all$cell_id)
so_dist_all$nearest_dist_um <- dist_lookup_all[colnames(so_dist_all)]

for (tp in c("LCMV_1wpi", "LCMV_6wpi")) {
  tp_slug  <- gsub("[^a-z0-9]", "_", tolower(tp))
  tp_label <- SAMPLE_LABELS[tp]
  message("\n--- Distance threshold analysis: ", tp, " ---")

  cells_tp <- colnames(so_dist_all)[so_dist_all$sample == tp]
  if (length(cells_tp) < 20) {
    message("  Skipping: too few cells (", length(cells_tp), ").")
    next
  }
  so_tp <- subset(so_dist_all, cells = cells_tp)

  dist_vals <- so_tp$nearest_dist_um
  cat("  Distance range: ", round(min(dist_vals, na.rm = TRUE), 1),
      " - ", round(max(dist_vals, na.rm = TRUE), 1), " µm\n")
  cat("  Quantiles (25/50/75%):",
      round(quantile(dist_vals, c(0.25, 0.5, 0.75), na.rm = TRUE), 1), "\n")

  # ---- Partie 1 : courbe n DEG vs seuil ----
  threshold_results <- lapply(DISTANCE_THRESHOLDS, function(thr) {
    cells_proche <- colnames(so_tp)[so_tp$nearest_dist_um <  thr]
    cells_loin   <- colnames(so_tp)[so_tp$nearest_dist_um >= thr]
    n_proche <- length(cells_proche)
    n_loin   <- length(cells_loin)
    message("  Seuil ", thr, " µm — proche: ", n_proche, " | loin: ", n_loin)

    if (n_proche < 5 || n_loin < 5) {
      message("  Skipping: insufficient cells.")
      return(list(
        threshold = thr, n_proche = n_proche, n_loin = n_loin,
        n_sig = 0L, n_up = 0L, n_down = 0L, deg = NULL
      ))
    }

    # Assign group identity
    so_thr <- so_tp
    so_thr$dist_group <- ifelse(
      so_thr$nearest_dist_um < thr, "proche", "loin"
    )
    so_thr <- SetIdent(so_thr, value = "dist_group")

    mk_thr <- tryCatch(
      FindMarkers(
        so_thr,
        ident.1         = "proche",
        ident.2         = "loin",
        only.pos        = FALSE,
        min.pct         = MIN_PCT,
        logfc.threshold = FC_THRESH,
        return.thresh   = 1,
        test.use        = TEST_USE,
        verbose         = FALSE
      ),
      error = function(e) {
        message("  ERROR: ", conditionMessage(e)); NULL
      }
    )

    if (is.null(mk_thr) || nrow(mk_thr) == 0) {
      return(list(
        threshold = thr, n_proche = n_proche, n_loin = n_loin,
        n_sig = 0L, n_up = 0L, n_down = 0L, deg = NULL
      ))
    }

    mk_thr$gene      <- rownames(mk_thr)
    mk_thr$threshold <- thr
    mk_thr$direction <- case_when(
      mk_thr$p_val_adj <= FDR_CUTOFF & mk_thr$avg_log2FC >  FC_CUTOFF ~ "up",
      mk_thr$p_val_adj <= FDR_CUTOFF & mk_thr$avg_log2FC < -FC_CUTOFF ~ "down",
      TRUE ~ "ns"
    )
    mk_thr <- mk_thr %>%
      arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
      select(gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj, direction, threshold)

    n_sig  <- sum(mk_thr$direction != "ns")
    n_up   <- sum(mk_thr$direction == "up")
    n_down <- sum(mk_thr$direction == "down")
    message("  ", nrow(mk_thr), " genes tested | sig=", n_sig,
            " (", n_up, " up / ", n_down, " down)")

    list(
      threshold = thr, n_proche = n_proche, n_loin = n_loin,
      n_sig = as.integer(n_sig),
      n_up  = as.integer(n_up),
      n_down = as.integer(n_down),
      deg   = mk_thr
    )
  })

  # Curve data frame
  curve_df <- do.call(rbind, lapply(threshold_results, function(r) {
    data.frame(
      threshold = r$threshold,
      n_proche  = r$n_proche,
      n_loin    = r$n_loin,
      n_sig     = r$n_sig,
      n_up      = r$n_up,
      n_down    = r$n_down
    )
  }))

  cond_color <- unname(COND_PALETTE[tp_label])
  if (is.na(cond_color)) cond_color <- "steelblue"

  p_curve <- ggplot(curve_df, aes(x = threshold, y = n_sig)) +
    geom_line(colour = cond_color, linewidth = 1) +
    geom_point(colour = cond_color, size = 3) +
    geom_text(aes(label = n_sig), vjust = 1.5, size = 3.2, colour = "grey20") +
    scale_x_continuous(
      breaks = DISTANCE_THRESHOLDS,
      labels = paste0(DISTANCE_THRESHOLDS, " \u00b5m")
    ) +
    labs(
      title    = paste0("N DEG significatifs vs seuil de distance \u2014 ", tp_label),
      subtitle = paste0(
        "Microglia (P2ry12) | proche (<seuil) vs loin (\u2265seuil)\n",
        "FDR \u2264 ", FDR_CUTOFF, " | |log2FC| > ", FC_CUTOFF
      ),
      x = "Seuil de distance (\u00b5m)",
      y = "N DEG significatifs"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title       = element_text(face = "bold", size = 11, hjust = 0),
      plot.subtitle    = element_text(size = 8, colour = "grey40", hjust = 0),
      panel.grid.minor = element_blank()
    )

  save_fig(p_curve, paste0("fig_ndeg_by_distance_", tp_slug),
           width = 5, height = 4)

  # ---- Partie 2 : volcano au seuil avec le plus de DEG ----
  has_deg <- sapply(threshold_results, function(r) !is.null(r$deg) && r$n_sig > 0)

  if (any(has_deg)) {
    best_idx <- which.max(sapply(threshold_results, function(r) r$n_sig))
    best     <- threshold_results[[best_idx]]
    thr_best <- best$threshold
    message("  Best threshold: ", thr_best, " µm (n_sig=", best$n_sig, ")")

    mk_max <- best$deg
    write.csv(
      mk_max,
      file.path(out_dir, paste0("DEG_distance_max_", tp_slug, ".csv")),
      row.names = FALSE
    )
    message("  Saved: DEG_distance_max_", tp_slug, ".csv")

    mk_max$neg_log10_fdr <- -log10(pmax(mk_max$p_val_adj, .Machine$double.eps))
    x_cap <- max(3.5, quantile(abs(mk_max$avg_log2FC), 0.999, na.rm = TRUE) * 1.05)
    y_cap <- max(5,   quantile(mk_max$neg_log10_fdr,   0.999, na.rm = TRUE) * 1.05)
    mk_max$x_plot <- pmin(pmax(mk_max$avg_log2FC, -x_cap), x_cap)
    mk_max$y_plot <- pmin(mk_max$neg_log10_fdr, y_cap)

    sig_max <- mk_max %>% filter(direction != "ns")
    if (nrow(sig_max) > 0) {
      lab_genes_max <- sig_max %>%
        arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
        slice_head(n = TOP_N_LABEL) %>%
        pull(gene)
    } else {
      lab_genes_max <- mk_max %>%
        arrange(p_val_adj) %>%
        slice_head(n = min(TOP_N_LABEL, 10)) %>%
        pull(gene)
    }
    mk_max$label <- ifelse(mk_max$gene %in% lab_genes_max, mk_max$gene, NA_character_)

    n_up_m   <- sum(mk_max$direction == "up",   na.rm = TRUE)
    n_down_m <- sum(mk_max$direction == "down", na.rm = TRUE)

    p_volc_max <- ggplot(mk_max,
                         aes(x = x_plot, y = y_plot, colour = direction)) +
      geom_point(size = 1.2, alpha = 0.7, stroke = 0) +
      geom_text_repel(
        data          = mk_max %>% filter(!is.na(label)),
        aes(label     = label),
        size          = 2.6,
        fontface      = "italic",
        max.overlaps  = 20,
        box.padding   = 0.3,
        point.padding = 0.2,
        segment.size  = 0.3,
        segment.color = "grey50",
        show.legend   = FALSE
      ) +
      geom_hline(yintercept = -log10(FDR_CUTOFF),
                 linetype = "dashed", colour = "grey40", linewidth = 0.4) +
      geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF),
                 linetype = "dashed", colour = "grey40", linewidth = 0.4) +
      scale_colour_manual(
        values = direction_colors,
        labels = c(
          up   = paste0("Up in proche (<", thr_best, " \u00b5m) (n=", n_up_m,   ")"),
          down = paste0("Down in proche (<", thr_best, " \u00b5m) (n=", n_down_m, ")"),
          ns   = "Not significant"
        ),
        name = NULL
      ) +
      scale_x_continuous(limits = c(-x_cap, x_cap),
                         expand = expansion(mult = 0.02)) +
      scale_y_continuous(limits = c(0, y_cap),
                         expand = expansion(mult = c(0, 0.05))) +
      labs(
        title    = paste0("Microglia \u2014 proche (<", thr_best,
                          " \u00b5m) vs loin (\u2265", thr_best,
                          " \u00b5m) \u2014 ", tp_label),
        subtitle = paste0(nrow(mk_max), " genes tested  |  FDR \u2264 ", FDR_CUTOFF,
                          "  |  |log2FC| > ", FC_CUTOFF),
        x        = paste0("log2FC (proche / loin, seuil=", thr_best, " \u00b5m)"),
        y        = expression(-log[10](FDR))
      ) +
      theme_bw(base_size = 10) +
      theme(
        plot.title       = element_text(face = "bold", size = 11, hjust = 0),
        plot.subtitle    = element_text(size = 8, colour = "grey40", hjust = 0),
        legend.position  = "bottom",
        legend.text      = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey93", linewidth = 0.3)
      )

    save_fig(p_volc_max, paste0("fig_volcano_distance_max_", tp_slug),
             width = 6, height = 5.5)
  } else {
    message("  No significant DEGs at any threshold — skipping volcano.")
  }
}

# =============================================================
# Done
# =============================================================

message("\n=== Done. All outputs in: ", out_dir, " ===\n")
cat("Expected outputs:\n")
for (nm in names(sig_modules)) {
  sig_slug <- gsub("[^a-z0-9]", "_", tolower(nm))
  cat("  fig_module_score_", sig_slug, "_violin.pdf/jpg\n", sep = "")
  cat("  fig_module_score_", sig_slug, "_spatial.pdf/jpg\n", sep = "")
}
for (tp in TIMEPOINTS) {
  tp_slug <- gsub("[^a-z0-9]", "_", tolower(tp))
  cat("  DEG_microglia_in_vs_out_", tp_slug, ".csv\n", sep = "")
  cat("  fig_volcano_in_vs_out_", tp_slug, ".pdf/jpg\n", sep = "")
}
for (tp in c("LCMV_1wpi", "LCMV_6wpi")) {
  tp_slug <- gsub("[^a-z0-9]", "_", tolower(tp))
  cat("  fig_ndeg_by_distance_", tp_slug, ".pdf/jpg\n", sep = "")
  cat("  DEG_distance_max_", tp_slug, ".csv\n", sep = "")
  cat("  fig_volcano_distance_max_", tp_slug, ".pdf/jpg\n", sep = "")
}
cat("  module_scores_per_cell.csv\n")
