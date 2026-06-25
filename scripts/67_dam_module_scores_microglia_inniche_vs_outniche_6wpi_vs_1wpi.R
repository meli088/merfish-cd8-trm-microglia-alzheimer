#!/usr/bin/env Rscript
# ============================================================================
# Script 67 — DAM module scores: Microglia IN-niche (C1qa) vs OUT-niche (P2ry12)
#             Comparison 6wpi vs 1wpi with both Upregulated and Downregulated DAM
#
# GOAL: Violin plots + Wilcoxon stats for:
#   - C1qa IN-niche (from object 08)
#   - P2ry12 OUT-niche (from object 04)
#   Stratified by timepoint (1wpi vs 6wpi) and signature (Up/Down DAM)
#
# Outputs:
#   - fig_dam_up_microglia_inniche_vs_outniche_6wpi_vs_1wpi.pdf/.jpg
#   - fig_dam_down_microglia_inniche_vs_outniche_6wpi_vs_1wpi.pdf/.jpg
#   - dam_module_score_inniche_vs_outniche_stats_6wpi_vs_1wpi.csv
#
# ============================================================================

set.seed(1997)

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(Banksy)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(Cairo)
})

# Utility function for p-value stars
pval_stars <- function(p) {
  case_when(
    is.na(p)  ~ "n.s.",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "n.s."
  )
}

# Save figure function
save_fig <- function(p, name, w = 7, h = 6) {
  dir_out <- "outputs/banksy/immune_acod1/analysis/figures"
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  
  pdf_path <- file.path(dir_out, paste0(name, ".pdf"))
  jpg_path <- file.path(dir_out, paste0(name, ".jpg"))
  
  CairoPDF(pdf_path, width = w, height = h)
  print(p)
  dev.off()
  
  CairoJPEG(jpg_path, width = w * 150, height = h * 150, res = 150)
  print(p)
  dev.off()
  
  message(sprintf("  ✓ Saved: %s", name))
}

# ============================================================================
# SECTION 1: Load objects and DAM signatures
# ============================================================================

message("\n=== SECTION 1: Loading objects and signatures ===\n")

# Load immune niche object (IN-niche C1qa)
se_immune <- readRDS("objects/08_immune_annotated_lam02_res03.rds")
message(sprintf("  Loaded immune object: %d cells", ncol(se_immune)))

# Load global object (OUT-niche P2ry12) — extract annotation via BANKSY clusters
se_global <- readRDS("objects/04_banksy_joint_lam08_after_bloc3.rds")
message(sprintf("  Loaded global object: %d cells", ncol(se_global)))

# Helper: find BANKSY cluster column (lambda=0.2, res=0.9)
find_cl_col <- function(se, lam, res) {
  all_cols <- Banksy::clusterNames(se)
  lam_str  <- gsub("\\.", "\\\\.", as.character(lam))
  lam_cols <- all_cols[grep(paste0("lam", lam_str), all_cols)]
  if (length(lam_cols) == 0) return(NULL)
  res_vals <- suppressWarnings(as.numeric(sub(".*_res", "", lam_cols)))
  idx      <- which(!is.na(res_vals) & abs(res_vals - res) < 1e-8)
  if (length(idx) == 0) return(NULL)
  lam_cols[idx[1]]
}

LAM <- 0.2
RES <- 0.9
cl_col <- find_cl_col(se_global, LAM, RES)
if (is.null(cl_col)) {
  stop("BANKSY cluster column not found for lambda=", LAM, " res=", RES)
}
message(sprintf("  Using cluster column: %s", cl_col))

# Reconstruct BANKSY annotations from CSV
csv_path <- "ncells_by_sample_lam02_res09_joint_long.csv"
annot_long <- read_delim(
  csv_path, delim = ";", locale = locale(decimal_mark = "."),
  show_col_types = FALSE, trim_ws = TRUE
) %>%
  select(-matches("^Unnamed")) %>%
  mutate(
    banksy_domain = as.character(banksy_domain),
    annotation = trimws(as.character(annotation))
  )

annotation_map <- annot_long %>%
  filter(!is.na(annotation), annotation != "") %>%
  distinct(banksy_domain, annotation)

domain_labels <- paste0("Domain_", as.character(SummarizedExperiment::colData(se_global)[[cl_col]]))
anno_lookup <- setNames(annotation_map$annotation, annotation_map$banksy_domain)
annotation_global <- ifelse(
  !is.na(anno_lookup[domain_labels]) & anno_lookup[domain_labels] != "",
  anno_lookup[domain_labels],
  "Non annotated"
)

message(sprintf("  Reconstructed BANKSY annotations (%d mappings)", nrow(annotation_map)))

# Load DAM signatures
upregulated_dam <- read_csv("outputs/banksy/dam_signature_curation/Upregulated_DAM.csv", show_col_types = FALSE)
downregulated_dam <- read_csv("outputs/banksy/dam_signature_curation/Downregulated_DAM.csv", show_col_types = FALSE)

# Panel uses mouse Titlecase; signature files use human ALLCAPS
panel_genes <- rownames(se_global)
panel_title <- tools::toTitleCase(tolower(panel_genes))
panel_name_map <- setNames(panel_genes, panel_title)

# Harmonise and match signatures
harm_intersect <- function(sig_genes) {
  sig_title <- tools::toTitleCase(tolower(sig_genes))
  hits <- sig_title[sig_title %in% panel_title]
  matched <- panel_name_map[hits]
  unique(as.character(matched[!is.na(matched)]))
}

up_genes_orig <- harm_intersect(upregulated_dam$gene)
down_genes_orig <- harm_intersect(downregulated_dam$gene)

message(sprintf("  Upregulated DAM: %d genes matched", length(up_genes_orig)))
message(sprintf("  Downregulated DAM: %d genes matched", length(down_genes_orig)))

# ============================================================================
# SECTION 2: Extract and process populations
# ============================================================================

message("\n=== SECTION 2: Extract populations ===\n")

# --- IN-NICHE: C1qa microglies from immune object (08)
cell_type_niche <- as.character(SummarizedExperiment::colData(se_immune)$cell_type)
inniche_label_candidates <- c("Activated microglia (C1qa)", "Microglia (C1qa)")
inniche_label <- inniche_label_candidates[inniche_label_candidates %in% unique(cell_type_niche)][1]

if (is.na(inniche_label)) {
  stop("IN-niche labels not found. Present: ", paste(sort(unique(cell_type_niche)), collapse = ", "))
}

inniche_idx <- which(cell_type_niche == inniche_label)
se_inniche <- se_immune[, inniche_idx]
message(sprintf("  IN-niche (%s): %d cells", inniche_label, length(inniche_idx)))

# --- OUT-NICHE: P2ry12 microglies from global object (04) via BANKSY annotations
outniche_label <- "Microglia (P2ry12)"
if (!outniche_label %in% annotation_global) {
  stop("OUT-niche label '", outniche_label, "' not found. Present: ",
       paste(sort(unique(annotation_global)), collapse = ", "))
}

outniche_idx <- which(annotation_global == outniche_label)
se_outniche <- se_global[, outniche_idx]
message(sprintf("  OUT-niche (%s): %d cells", outniche_label, length(outniche_idx)))

# ============================================================================
# SECTION 3: Extract 1wpi and 6wpi subsets
# ============================================================================

message("\n=== SECTION 3: Subset by timepoint ===\n")

# IN-niche subsets
sample_inniche <- as.character(SummarizedExperiment::colData(se_inniche)$sample)
inniche_1wpi_idx <- which(sample_inniche == "LCMV_1wpi")
inniche_6wpi_idx <- which(sample_inniche == "LCMV_6wpi")

se_inniche_1wpi <- se_inniche[, inniche_1wpi_idx]
se_inniche_6wpi <- se_inniche[, inniche_6wpi_idx]

# OUT-niche subsets
sample_outniche <- as.character(SummarizedExperiment::colData(se_outniche)$sample)
outniche_1wpi_idx <- which(sample_outniche == "LCMV_1wpi")
outniche_6wpi_idx <- which(sample_outniche == "LCMV_6wpi")

se_outniche_1wpi <- se_outniche[, outniche_1wpi_idx]
se_outniche_6wpi <- se_outniche[, outniche_6wpi_idx]

message(sprintf("  IN-niche: 1wpi=%d cells, 6wpi=%d cells", ncol(se_inniche_1wpi), ncol(se_inniche_6wpi)))
message(sprintf("  OUT-niche: 1wpi=%d cells, 6wpi=%d cells", ncol(se_outniche_1wpi), ncol(se_outniche_6wpi)))

# ============================================================================
# SECTION 4: Compute DAM scores for IN-niche
# ============================================================================

message("\n=== SECTION 4: Compute DAM scores (IN-niche) ===\n")

inniche_1wpi_scores <- tibble()
inniche_6wpi_scores <- tibble()

# Process IN-niche 1wpi
if (ncol(se_inniche_1wpi) >= 10) {
  so_1wpi <- suppressWarnings(as.Seurat(se_inniche_1wpi, counts = assayNames(se_inniche_1wpi)[1], data = NULL))
  dassay <- if ("RNA" %in% SeuratObject::Assays(so_1wpi)) "RNA" else SeuratObject::Assays(so_1wpi)[1]
  DefaultAssay(so_1wpi) <- dassay
  so_1wpi <- NormalizeData(so_1wpi, assay = dassay, verbose = FALSE)
  so_1wpi <- AddModuleScore(
    so_1wpi,
    features = list(Upregulated_DAM = up_genes_orig, Downregulated_DAM = down_genes_orig),
    name = "", ctrl = 5, verbose = FALSE
  )
  inniche_1wpi_scores <- tibble(
    timepoint = "1wpi",
    niche_status = "IN-niche",
    cell_id = colnames(so_1wpi),
    Upregulated_DAM = so_1wpi@meta.data[, "1"],
    Downregulated_DAM = so_1wpi@meta.data[, "2"]
  )
  message(sprintf("  ✓ IN-niche 1wpi: computed scores (%d cells)", nrow(inniche_1wpi_scores)))
} else {
  message("  ✗ IN-niche 1wpi: insufficient cells")
}

# Process IN-niche 6wpi
if (ncol(se_inniche_6wpi) >= 10) {
  so_6wpi <- suppressWarnings(as.Seurat(se_inniche_6wpi, counts = assayNames(se_inniche_6wpi)[1], data = NULL))
  dassay <- if ("RNA" %in% SeuratObject::Assays(so_6wpi)) "RNA" else SeuratObject::Assays(so_6wpi)[1]
  DefaultAssay(so_6wpi) <- dassay
  so_6wpi <- NormalizeData(so_6wpi, assay = dassay, verbose = FALSE)
  so_6wpi <- AddModuleScore(
    so_6wpi,
    features = list(Upregulated_DAM = up_genes_orig, Downregulated_DAM = down_genes_orig),
    name = "", ctrl = 5, verbose = FALSE
  )
  inniche_6wpi_scores <- tibble(
    timepoint = "6wpi",
    niche_status = "IN-niche",
    cell_id = colnames(so_6wpi),
    Upregulated_DAM = so_6wpi@meta.data[, "1"],
    Downregulated_DAM = so_6wpi@meta.data[, "2"]
  )
  message(sprintf("  ✓ IN-niche 6wpi: computed scores (%d cells)", nrow(inniche_6wpi_scores)))
} else {
  message("  ✗ IN-niche 6wpi: insufficient cells")
}

inniche_scores <- bind_rows(inniche_1wpi_scores, inniche_6wpi_scores)

# ============================================================================
# SECTION 5: Compute DAM scores for OUT-niche
# ============================================================================

message("\n=== SECTION 5: Compute DAM scores (OUT-niche) ===\n")

outniche_1wpi_scores <- tibble()
outniche_6wpi_scores <- tibble()

# Process OUT-niche 1wpi
if (ncol(se_outniche_1wpi) >= 10) {
  so_1wpi <- suppressWarnings(as.Seurat(se_outniche_1wpi, counts = assayNames(se_outniche_1wpi)[1], data = NULL))
  dassay <- if ("RNA" %in% SeuratObject::Assays(so_1wpi)) "RNA" else SeuratObject::Assays(so_1wpi)[1]
  DefaultAssay(so_1wpi) <- dassay
  so_1wpi <- NormalizeData(so_1wpi, assay = dassay, verbose = FALSE)
  so_1wpi <- AddModuleScore(
    so_1wpi,
    features = list(Upregulated_DAM = up_genes_orig, Downregulated_DAM = down_genes_orig),
    name = "", ctrl = 5, verbose = FALSE
  )
  outniche_1wpi_scores <- tibble(
    timepoint = "1wpi",
    niche_status = "OUT-niche",
    cell_id = colnames(so_1wpi),
    Upregulated_DAM = so_1wpi@meta.data[, "1"],
    Downregulated_DAM = so_1wpi@meta.data[, "2"]
  )
  message(sprintf("  ✓ OUT-niche 1wpi: computed scores (%d cells)", nrow(outniche_1wpi_scores)))
} else {
  message("  ✗ OUT-niche 1wpi: insufficient cells")
}

# Process OUT-niche 6wpi
if (ncol(se_outniche_6wpi) >= 10) {
  so_6wpi <- suppressWarnings(as.Seurat(se_outniche_6wpi, counts = assayNames(se_outniche_6wpi)[1], data = NULL))
  dassay <- if ("RNA" %in% SeuratObject::Assays(so_6wpi)) "RNA" else SeuratObject::Assays(so_6wpi)[1]
  DefaultAssay(so_6wpi) <- dassay
  so_6wpi <- NormalizeData(so_6wpi, assay = dassay, verbose = FALSE)
  so_6wpi <- AddModuleScore(
    so_6wpi,
    features = list(Upregulated_DAM = up_genes_orig, Downregulated_DAM = down_genes_orig),
    name = "", ctrl = 5, verbose = FALSE
  )
  outniche_6wpi_scores <- tibble(
    timepoint = "6wpi",
    niche_status = "OUT-niche",
    cell_id = colnames(so_6wpi),
    Upregulated_DAM = so_6wpi@meta.data[, "1"],
    Downregulated_DAM = so_6wpi@meta.data[, "2"]
  )
  message(sprintf("  ✓ OUT-niche 6wpi: computed scores (%d cells)", nrow(outniche_6wpi_scores)))
} else {
  message("  ✗ OUT-niche 6wpi: insufficient cells")
}

outniche_scores <- bind_rows(outniche_1wpi_scores, outniche_6wpi_scores)

# Combine all scores
all_scores <- bind_rows(inniche_scores, outniche_scores) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("1wpi", "6wpi")),
    niche_status = factor(niche_status, levels = c("OUT-niche", "IN-niche"))
  )

message(sprintf("\n  Total cells: %d", nrow(all_scores)))

# ============================================================================
# SECTION 6: Wilcoxon tests and statistics
# ============================================================================

message("\n=== SECTION 6: Wilcoxon tests ===\n")

stats_list <- list()

# Upregulated DAM, 6wpi vs 1wpi, IN-niche
if (nrow(inniche_1wpi_scores) > 0 & nrow(inniche_6wpi_scores) > 0) {
  p_val <- wilcox.test(
    inniche_6wpi_scores$Upregulated_DAM,
    inniche_1wpi_scores$Upregulated_DAM,
    exact = FALSE
  )$p.value
  
  stats_list[[1]] <- tibble(
    niche_type = "IN-niche (C1qa)",
    signature = "Upregulated DAM",
    n_1wpi = nrow(inniche_1wpi_scores),
    n_6wpi = nrow(inniche_6wpi_scores),
    median_1wpi = median(inniche_1wpi_scores$Upregulated_DAM, na.rm = TRUE),
    median_6wpi = median(inniche_6wpi_scores$Upregulated_DAM, na.rm = TRUE),
    delta = median(inniche_6wpi_scores$Upregulated_DAM, na.rm = TRUE) - 
            median(inniche_1wpi_scores$Upregulated_DAM, na.rm = TRUE),
    p_value_wilcoxon = p_val,
    sig = pval_stars(p_val)
  )
  message(sprintf("  IN-niche Upregulated DAM: p=%.4f %s", p_val, pval_stars(p_val)))
}

# Downregulated DAM, 6wpi vs 1wpi, IN-niche
if (nrow(inniche_1wpi_scores) > 0 & nrow(inniche_6wpi_scores) > 0) {
  p_val <- wilcox.test(
    inniche_6wpi_scores$Downregulated_DAM,
    inniche_1wpi_scores$Downregulated_DAM,
    exact = FALSE
  )$p.value
  
  stats_list[[2]] <- tibble(
    niche_type = "IN-niche (C1qa)",
    signature = "Downregulated DAM",
    n_1wpi = nrow(inniche_1wpi_scores),
    n_6wpi = nrow(inniche_6wpi_scores),
    median_1wpi = median(inniche_1wpi_scores$Downregulated_DAM, na.rm = TRUE),
    median_6wpi = median(inniche_6wpi_scores$Downregulated_DAM, na.rm = TRUE),
    delta = median(inniche_6wpi_scores$Downregulated_DAM, na.rm = TRUE) - 
            median(inniche_1wpi_scores$Downregulated_DAM, na.rm = TRUE),
    p_value_wilcoxon = p_val,
    sig = pval_stars(p_val)
  )
  message(sprintf("  IN-niche Downregulated DAM: p=%.4f %s", p_val, pval_stars(p_val)))
}

# Upregulated DAM, 6wpi vs 1wpi, OUT-niche
if (nrow(outniche_1wpi_scores) > 0 & nrow(outniche_6wpi_scores) > 0) {
  p_val <- wilcox.test(
    outniche_6wpi_scores$Upregulated_DAM,
    outniche_1wpi_scores$Upregulated_DAM,
    exact = FALSE
  )$p.value
  
  stats_list[[3]] <- tibble(
    niche_type = "OUT-niche (P2ry12)",
    signature = "Upregulated DAM",
    n_1wpi = nrow(outniche_1wpi_scores),
    n_6wpi = nrow(outniche_6wpi_scores),
    median_1wpi = median(outniche_1wpi_scores$Upregulated_DAM, na.rm = TRUE),
    median_6wpi = median(outniche_6wpi_scores$Upregulated_DAM, na.rm = TRUE),
    delta = median(outniche_6wpi_scores$Upregulated_DAM, na.rm = TRUE) - 
            median(outniche_1wpi_scores$Upregulated_DAM, na.rm = TRUE),
    p_value_wilcoxon = p_val,
    sig = pval_stars(p_val)
  )
  message(sprintf("  OUT-niche Upregulated DAM: p=%.4f %s", p_val, pval_stars(p_val)))
}

# Downregulated DAM, 6wpi vs 1wpi, OUT-niche
if (nrow(outniche_1wpi_scores) > 0 & nrow(outniche_6wpi_scores) > 0) {
  p_val <- wilcox.test(
    outniche_6wpi_scores$Downregulated_DAM,
    outniche_1wpi_scores$Downregulated_DAM,
    exact = FALSE
  )$p.value
  
  stats_list[[4]] <- tibble(
    niche_type = "OUT-niche (P2ry12)",
    signature = "Downregulated DAM",
    n_1wpi = nrow(outniche_1wpi_scores),
    n_6wpi = nrow(outniche_6wpi_scores),
    median_1wpi = median(outniche_1wpi_scores$Downregulated_DAM, na.rm = TRUE),
    median_6wpi = median(outniche_6wpi_scores$Downregulated_DAM, na.rm = TRUE),
    delta = median(outniche_6wpi_scores$Downregulated_DAM, na.rm = TRUE) - 
            median(outniche_1wpi_scores$Downregulated_DAM, na.rm = TRUE),
    p_value_wilcoxon = p_val,
    sig = pval_stars(p_val)
  )
  message(sprintf("  OUT-niche Downregulated DAM: p=%.4f %s", p_val, pval_stars(p_val)))
}

stats_df <- bind_rows(stats_list)

# Save stats
dir_out <- "outputs/banksy/immune_acod1/analysis/figures"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
stats_path <- file.path(dir_out, "dam_module_score_inniche_vs_outniche_stats_6wpi_vs_1wpi.csv")
write_csv(stats_df, stats_path)
message(sprintf("\n  ✓ Saved stats to: %s", stats_path))

# ============================================================================
# SECTION 7: Violin plots for Upregulated DAM
# ============================================================================

message("\n=== SECTION 7: Create violin plots ===\n")

# Color palette
pal_timepoint <- c("1wpi" = "#56B4E9", "6wpi" = "#F28E2B")

# Compute y-position for annotations
compute_y_annotation <- function(df, col_name) {
  y_max <- max(df[[col_name]], na.rm = TRUE)
  y_min <- min(df[[col_name]], na.rm = TRUE)
  y_range <- y_max - y_min
  y_max + 0.08 * y_range
}

# Upregulated DAM
y_pos_up <- compute_y_annotation(all_scores, "Upregulated_DAM")

# Add significance annotations
p_inniche_up <- 0.38423488749279683  # IN-niche Upregulated
p_outniche_up <- 4.280054395379414e-130  # OUT-niche Upregulated

stars_inniche_up <- pval_stars(p_inniche_up)
stars_outniche_up <- pval_stars(p_outniche_up)

p_up <- ggplot(all_scores, aes(x = niche_status, y = Upregulated_DAM, fill = timepoint)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.2, outlier.size = 0.5, outlier.alpha = 0.3,
               colour = "grey30", fill = "white", position = position_dodge(0.9)) +
  annotate("text", x = 1, y = y_pos_up, label = stars_inniche_up,
           size = 5, fontface = "bold", colour = "grey20") +
  annotate("text", x = 2, y = y_pos_up, label = stars_outniche_up,
           size = 5, fontface = "bold", colour = "grey20") +
  scale_fill_manual(values = pal_timepoint, name = "Timepoint") +
  labs(
    title = "Upregulated DAM: 6wpi vs 1wpi",
    subtitle = "Microglia IN-niche (C1qa) vs OUT-niche (P2ry12)",
    x = NULL,
    y = "DAM module score (AddModuleScore)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, colour = "grey60"),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

save_fig(p_up, "fig_dam_up_microglia_inniche_vs_outniche_6wpi_vs_1wpi", 6.5, 5)

# ============================================================================
# SECTION 8: Violin plots for Downregulated DAM
# ============================================================================

y_pos_down <- compute_y_annotation(all_scores, "Downregulated_DAM")

p_inniche_down <- 0.5791845521768408  # IN-niche Downregulated
p_outniche_down <- 1.8309345052148906e-48  # OUT-niche Downregulated

stars_inniche_down <- pval_stars(p_inniche_down)
stars_outniche_down <- pval_stars(p_outniche_down)

p_down <- ggplot(all_scores, aes(x = niche_status, y = Downregulated_DAM, fill = timepoint)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.2, outlier.size = 0.5, outlier.alpha = 0.3,
               colour = "grey30", fill = "white", position = position_dodge(0.9)) +
  annotate("text", x = 1, y = y_pos_down, label = stars_inniche_down,
           size = 5, fontface = "bold", colour = "grey20") +
  annotate("text", x = 2, y = y_pos_down, label = stars_outniche_down,
           size = 5, fontface = "bold", colour = "grey20") +
  scale_fill_manual(values = pal_timepoint, name = "Timepoint") +
  labs(
    title = "Downregulated DAM: 6wpi vs 1wpi",
    subtitle = "Microglia IN-niche (C1qa) vs OUT-niche (P2ry12)",
    x = NULL,
    y = "DAM module score (AddModuleScore)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, colour = "grey60"),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

save_fig(p_down, "fig_dam_down_microglia_inniche_vs_outniche_6wpi_vs_1wpi", 6.5, 5)

message("\n=== Script completed successfully ===\n")
