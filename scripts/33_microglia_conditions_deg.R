#!/usr/bin/env Rscript
# =============================================================
# Script: 33_microglia_conditions_deg.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-05
#
# Goal:
#   DEG analysis within the global Microglia (P2ry12) population,
#   comparing selected conditions, then cross-referencing with the
#   in-niche vs out-niche DEG produced by script 32.
#
# Contrasts (Section 2):
#   LCMV_6wpi  vs mock_6wpi
#   LCMV_6wpi  vs LCMV_1wpi
#   LCMV_1wpi  vs mock_6wpi
#
# Overlap (Section 3):
#   For 1wpi: conditions DEG (LCMV_1wpi vs mock) ∩ niche DEG (in vs out, 1wpi)
#   For 6wpi: conditions DEG (LCMV_6wpi vs mock) ∩ niche DEG (in vs out, 6wpi)
#
# Inputs:
#   objects/04_banksy_joint_lam08_after_bloc3.rds
#   ncells_by_sample_lam02_res09_joint_long.csv
#   outputs/banksy/microglia_dam_niche/DEG_microglia_in_vs_out_lcmv_1wpi.csv
#   outputs/banksy/microglia_dam_niche/DEG_microglia_in_vs_out_lcmv_6wpi.csv
#
# Outputs (outputs/banksy/microglia_conditions_deg/):
#   [Sec 2] DEG_microglia_[contrast].csv
#           fig_volcano_microglia_[contrast].pdf/jpg
#   [Sec 3] common_genes_lcmv_1wpi.csv
#           common_genes_lcmv_6wpi.csv
#           fig_overlap_lcmv_1wpi.pdf/jpg
#           fig_overlap_lcmv_6wpi.pdf/jpg
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
MICROGLIA_GLOBAL_LABEL <- "Microglia (P2ry12)"
LAM                    <- 0.2
RES_TARGET             <- 0.9
FDR_CUTOFF             <- 0.05
FC_CUTOFF              <- 0.25
TOP_N_LABEL            <- 15
SEED                   <- 1997

out_dir <- "outputs/banksy/microglia_conditions_deg"
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

# Contrasts: list of (ident.1, ident.2) pairs
CONTRASTS <- list(
  list(id1 = "LCMV_6wpi", id2 = "mock_6wpi"),
  list(id1 = "LCMV_6wpi", id2 = "LCMV_1wpi"),
  list(id1 = "LCMV_1wpi", id2 = "mock_6wpi")
)

# FindMarkers settings
MIN_PCT   <- 0.05
FC_THRESH <- 0.1
TEST_USE  <- "wilcox"

# =============================================================
# Helpers
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

contrast_slug <- function(id1, id2) {
  gsub("[^a-z0-9]", "_", tolower(paste0(id1, "_vs_", id2)))
}

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
# SECTION 1 — Load object and extract Microglia (P2ry12)
# =============================================================

message("\n=== SECTION 1: Loading object and extracting microglia ===\n")

obj_file <- file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds")
stopifnot("Object file not found" = file.exists(obj_file))

message("Loading: ", obj_file)
se_global <- readRDS(obj_file)
message("  ", ncol(se_global), " cells | class: ", class(se_global)[1])

# --- Find BANKSY cluster column ---
cl_col <- find_cl_col(se_global, LAM, RES_TARGET)
if (is.null(cl_col)) {
  stop("BANKSY cluster column not found for lambda=", LAM, " res=", RES_TARGET,
       "\nAvailable: ", paste(clusterNames(se_global), collapse = ", "))
}
message("  Cluster column: ", cl_col)

# --- Reconstruct annotation from CSV (same pattern as scripts 29/31/32) ---
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

cat("\nAnnotation distribution (all cells):\n")
print(sort(table(annotation_global), decreasing = TRUE))

if (!MICROGLIA_GLOBAL_LABEL %in% annotation_global) {
  stop("Label '", MICROGLIA_GLOBAL_LABEL, "' not found in annotations.\n",
       "Present: ", paste(sort(unique(annotation_global)), collapse = ", "))
}

# --- Subset to Microglia (P2ry12) ---
mg_idx <- which(annotation_global == MICROGLIA_GLOBAL_LABEL)
se_mg  <- se_global[, mg_idx]
message("\n  Microglia (P2ry12) cells: ", ncol(se_mg))

sample_vec <- as.character(colData(se_mg)$sample)
cat("\nCells per condition:\n")
print(table(sample_vec))

# --- Convert to Seurat and normalize ---
assay_name    <- if ("counts" %in% assayNames(se_mg)) "counts" else assayNames(se_mg)[1]
so            <- suppressWarnings(as.Seurat(se_mg, counts = assay_name, data = NULL))
default_assay <- if ("RNA" %in% SeuratObject::Assays(so)) "RNA" else SeuratObject::Assays(so)[1]
DefaultAssay(so) <- default_assay

needs_norm <- tryCatch({
  dm <- GetAssayData(so, assay = default_assay, layer = "data")
  nrow(dm) == 0 || ncol(dm) == 0
}, error = function(e) TRUE)
if (needs_norm) {
  so <- NormalizeData(so, assay = default_assay, verbose = FALSE)
  message("  NormalizeData done")
}

so <- SetIdent(so, value = "sample")

# =============================================================
# SECTION 2 — DEG between conditions: volcano plots
# =============================================================

message("\n=== SECTION 2: DEG between conditions ===\n")

make_volcano <- function(mk, id1, id2) {
  slug  <- contrast_slug(id1, id2)
  lab1  <- SAMPLE_LABELS[id1]
  lab2  <- SAMPLE_LABELS[id2]

  mk$neg_log10_fdr <- -log10(pmax(mk$p_val_adj, .Machine$double.eps))
  x_cap <- max(3.5, quantile(abs(mk$avg_log2FC), 0.999, na.rm = TRUE) * 1.05)
  y_cap <- max(5,   quantile(mk$neg_log10_fdr,   0.999, na.rm = TRUE) * 1.05)
  mk$x_plot <- pmin(pmax(mk$avg_log2FC, -x_cap), x_cap)
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

  p <- ggplot(mk, aes(x = x_plot, y = y_plot, colour = direction)) +
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
        up   = paste0("Up in ", lab1, " (n=", n_up,   ")"),
        down = paste0("Down in ", lab1, " (n=", n_down, ")"),
        ns   = "Not significant"
      ),
      name = NULL
    ) +
    scale_x_continuous(limits = c(-x_cap, x_cap),
                       expand = expansion(mult = 0.02)) +
    scale_y_continuous(limits = c(0, y_cap),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(
      title    = paste0("Microglia (P2ry12): ", lab1, " vs ", lab2),
      subtitle = paste0(nrow(mk), " genes tested  |  FDR \u2264 ", FDR_CUTOFF,
                        "  |  |log2FC| > ", FC_CUTOFF),
      x        = paste0("log2FC (", lab1, " / ", lab2, ")"),
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

  save_fig(p, paste0("fig_volcano_microglia_", slug), width = 6, height = 5.5)
}

# Run all contrasts
all_deg <- list()   # keyed by "id1_vs_id2" (raw, not slugified)

for (ct in CONTRASTS) {
  id1  <- ct$id1
  id2  <- ct$id2
  slug <- contrast_slug(id1, id2)
  n1   <- sum(sample_vec == id1, na.rm = TRUE)
  n2   <- sum(sample_vec == id2, na.rm = TRUE)
  message("\n--- ", id1, " vs ", id2, "  (n=", n1, " vs n=", n2, ") ---")

  if (n1 < 5 || n2 < 5) {
    message("  Skipping: insufficient cells (< 5 in one group).")
    next
  }

  mk <- tryCatch(
    FindMarkers(
      so,
      ident.1         = id1,
      ident.2         = id2,
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
  mk$contrast  <- paste0(id1, "_vs_", id2)
  mk$direction <- case_when(
    mk$p_val_adj <= FDR_CUTOFF & mk$avg_log2FC >  FC_CUTOFF ~ "up",
    mk$p_val_adj <= FDR_CUTOFF & mk$avg_log2FC < -FC_CUTOFF ~ "down",
    TRUE ~ "ns"
  )
  mk <- mk %>%
    arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
    select(gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj, direction, contrast)

  n_up   <- sum(mk$direction == "up")
  n_down <- sum(mk$direction == "down")
  message("  ", nrow(mk), " genes tested | ", n_up, " up | ", n_down,
          " down (FDR<=", FDR_CUTOFF, ", |log2FC|>", FC_CUTOFF, ")")

  write.csv(mk,
            file.path(out_dir, paste0("DEG_microglia_", slug, ".csv")),
            row.names = FALSE)
  message("  Saved: DEG_microglia_", slug, ".csv")

  all_deg[[paste0(id1, "_vs_", id2)]] <- mk

  make_volcano(mk, id1, id2)
}

if (length(all_deg) == 0) stop("No DEG results for any contrast.")

# =============================================================
# SECTION 3 — Cross-comparison with in niche vs out niche DEGs
# =============================================================

message("\n=== SECTION 3: Overlap conditions DEG vs niche DEG ===\n")

niche_dir <- file.path("outputs", "banksy", "microglia_dam_niche")

# Mapping: timepoint key → (conditions contrast key, niche DEG file)
overlap_map <- list(
  lcmv_1wpi = list(
    tp_label     = "LCMV 1 wpi",
    cond_key     = "LCMV_1wpi_vs_mock_6wpi",
    niche_file   = file.path(niche_dir, "DEG_microglia_in_vs_out_lcmv_1wpi.csv")
  ),
  lcmv_6wpi = list(
    tp_label     = "LCMV 6 wpi",
    cond_key     = "LCMV_6wpi_vs_mock_6wpi",
    niche_file   = file.path(niche_dir, "DEG_microglia_in_vs_out_lcmv_6wpi.csv")
  )
)

for (tp_slug in names(overlap_map)) {
  om         <- overlap_map[[tp_slug]]
  tp_label   <- om$tp_label
  cond_key   <- om$cond_key
  niche_file <- om$niche_file

  message("\n--- ", tp_label, " ---")

  # --- Conditions DEG (from Section 2) ---
  if (!cond_key %in% names(all_deg)) {
    message("  Conditions contrast '", cond_key, "' not available — skipping.")
    next
  }
  cond_deg  <- all_deg[[cond_key]]
  cond_sig  <- cond_deg %>%
    filter(direction != "ns") %>%
    pull(gene)
  message("  Conditions DEG sig (FDR<=", FDR_CUTOFF, ", |log2FC|>", FC_CUTOFF,
          "): ", length(cond_sig))

  # --- Niche DEG (from script 32) ---
  if (!file.exists(niche_file)) {
    message("  Niche DEG file not found: ", niche_file, " — skipping.")
    next
  }
  niche_deg  <- read.csv(niche_file, stringsAsFactors = FALSE)
  niche_sig  <- niche_deg %>%
    filter(p_val_adj <= FDR_CUTOFF, abs(avg_log2FC) > FC_CUTOFF) %>%
    pull(gene)
  message("  Niche DEG sig: ", length(niche_sig))

  # --- Compute overlap ---
  genes_only_cond  <- setdiff(cond_sig,  niche_sig)
  genes_common     <- intersect(cond_sig, niche_sig)
  genes_only_niche <- setdiff(niche_sig, cond_sig)

  message("  Only conditions: ", length(genes_only_cond))
  message("  Common:          ", length(genes_common))
  message("  Only niche:      ", length(genes_only_niche))

  if (length(genes_common) > 0) {
    # Build common gene table with both log2FCs
    common_cond  <- cond_deg  %>%
      filter(gene %in% genes_common) %>%
      select(gene,
             log2FC_conditions  = avg_log2FC,
             fdr_conditions     = p_val_adj,
             dir_conditions     = direction)
    common_niche <- niche_deg %>%
      filter(gene %in% genes_common) %>%
      select(gene,
             log2FC_niche       = avg_log2FC,
             fdr_niche          = p_val_adj,
             dir_niche          = direction)
    common_df <- common_cond %>%
      left_join(common_niche, by = "gene") %>%
      arrange(fdr_conditions, desc(abs(log2FC_conditions)))

    write.csv(common_df,
              file.path(out_dir, paste0("common_genes_", tp_slug, ".csv")),
              row.names = FALSE)
    message("  Saved: common_genes_", tp_slug, ".csv")

    if (nrow(common_df) > 0) {
      cat("  Common genes:\n")
      print(common_df[, c("gene", "log2FC_conditions", "dir_conditions",
                          "log2FC_niche", "dir_niche")],
            row.names = FALSE)
    }
  } else {
    message("  No common genes — skipping common_genes CSV.")
  }

  # --- Overlap bar plot ---
  # Horizontal bar showing 3 groups: only conditions | common | only niche
  bar_df <- data.frame(
    group = factor(
      c(
        paste0("Only conditions\n(", tp_label, " vs Mock)"),
        "Common\n(both)",
        "Only niche\n(In vs Out)"
      ),
      levels = c(
        paste0("Only conditions\n(", tp_label, " vs Mock)"),
        "Common\n(both)",
        "Only niche\n(In vs Out)"
      )
    ),
    n = c(length(genes_only_cond), length(genes_common), length(genes_only_niche)),
    fill_col = c("#56B4E9", "#7B2D8B", "#B2182B")
  )

  p_bar <- ggplot(bar_df, aes(x = group, y = n, fill = group)) +
    geom_col(width = 0.6, colour = "grey30", linewidth = 0.3) +
    geom_text(aes(label = n), vjust = -0.5, size = 3.5, colour = "grey20") +
    scale_fill_manual(values = setNames(bar_df$fill_col, bar_df$group),
                      guide  = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title    = paste0("DEG overlap — ", tp_label),
      subtitle = paste0(
        "Conditions: ", tp_label, " vs Mock  |  Niche: In vs Out\n",
        "FDR \u2264 ", FDR_CUTOFF, "  |  |log2FC| > ", FC_CUTOFF
      ),
      x = NULL,
      y = "N significant DEGs"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title       = element_text(face = "bold", size = 11, hjust = 0),
      plot.subtitle    = element_text(size = 8, colour = "grey40", hjust = 0),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x      = element_text(size = 9)
    )

  save_fig(p_bar, paste0("fig_overlap_", tp_slug), width = 5, height = 4.5)

  # --- Scatter: log2FC conditions vs log2FC niche for common genes ---
  if (length(genes_common) >= 3) {
    scatter_df <- common_cond %>%
      left_join(common_niche, by = "gene")

    p_scatter <- ggplot(scatter_df,
                        aes(x = log2FC_conditions, y = log2FC_niche,
                            label = gene)) +
      geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
      geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
      geom_point(colour = "#7B2D8B", size = 2.5, alpha = 0.8) +
      geom_text_repel(
        size          = 2.8,
        fontface      = "italic",
        max.overlaps  = 20,
        box.padding   = 0.3,
        point.padding = 0.2,
        segment.size  = 0.3,
        segment.color = "grey50"
      ) +
      labs(
        title    = paste0("Common DEGs — ", tp_label),
        subtitle = paste0("x: log2FC (", tp_label, " / Mock)  |  ",
                          "y: log2FC (In niche / Out niche)"),
        x        = paste0("log2FC conditions (", tp_label, " vs Mock)"),
        y        = "log2FC niche (In vs Out)"
      ) +
      theme_bw(base_size = 10) +
      theme(
        plot.title       = element_text(face = "bold", size = 11, hjust = 0),
        plot.subtitle    = element_text(size = 8, colour = "grey40", hjust = 0),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey93", linewidth = 0.3)
      )

    save_fig(p_scatter, paste0("fig_overlap_scatter_", tp_slug),
             width = 5, height = 5)
  }
}

# =============================================================
# Done
# =============================================================

message("\n=== Done. All outputs in: ", out_dir, " ===\n")
cat("Expected outputs:\n")
for (ct in CONTRASTS) {
  slug <- contrast_slug(ct$id1, ct$id2)
  cat("  DEG_microglia_", slug, ".csv\n", sep = "")
  cat("  fig_volcano_microglia_", slug, ".pdf/jpg\n", sep = "")
}
for (tp_slug in names(overlap_map)) {
  cat("  common_genes_", tp_slug, ".csv\n", sep = "")
  cat("  fig_overlap_", tp_slug, ".pdf/jpg\n", sep = "")
  cat("  fig_overlap_scatter_", tp_slug, ".pdf/jpg  (if >= 3 common genes)\n", sep = "")
}
