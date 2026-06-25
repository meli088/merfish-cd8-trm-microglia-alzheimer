#!/usr/bin/env Rscript
# ============================================================================
# Script 39 — Microglia niche analysis (revised)
#
# SECTION 1 : Volcano DEG microglie In vs Out niche, tous timepoints LCMV
#             mergés en un seul contraste
#             Output : outputs/banksy/microglia_dam_niche/
#                      fig_volcano_in_vs_out_all_timepoints_merged.pdf/jpg
#
# SECTION 2 : DAM module score In vs Out niche, tous timepoints LCMV mergés
#             Violin plot global + Wilcoxon
#             Output : outputs/banksy/microglia_dam_niche/
#                      fig_module_score_merged_violin.pdf/jpg
#
# SECTION 3 : Distance Microglia (P2ry12) -> T cells / Immune (Acod1)
#             au cours du temps — violin par condition
#             Output : outputs/banksy/microglia_dam_niche/
#                      fig_distance_microglia_to_tcell_overtime.pdf/jpg
#
# SECTION 4 : DAM signature (Upregulated_DAM) en fonction de bins de
#             distance aux T cells — 1wpi et 6wpi séparément
#             Output : outputs/banksy/microglia_dam_niche/
#                      fig_dam_score_by_distance_bins.pdf/jpg
#
# SECTION 5 : Violin + barplot DAM score In niche vs Hors niche
#             (Upregulated_DAM, LCMV merged)
#             Output : outputs/banksy/microglia_dam_niche/
#                      fig_dam_score_violin_in_vs_out_niche.pdf/jpg
#                      fig_dam_score_barplot_in_vs_out_niche.pdf/jpg
#                      fig_violin_barplot_in_vs_out_niche.pdf/jpg
#
# SECTION 6 : IN-niche only (Microglia C1qa)
#             Distance to nearest T cell / Immune (Acod1), puis DAM par bins
#             Figure combinée : violin par condition + barplot bins 1wpi vs 6wpi
#             Output : outputs/banksy/microglia_dam_niche/
#                      fig_distance_microglia_INNICHE_to_tcell.pdf/jpg
#
# Inputs :
#   objects/04_banksy_joint_lam08_after_bloc3.rds
#   objects/08_immune_annotated_lam02_res03.rds
#   ncells_by_sample_lam02_res09_joint_long.csv
#   outputs/banksy/dam_signature_curation/Upregulated_DAM.csv
#   outputs/banksy/dam_signature_curation/Downregulated_DAM.csv
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
  library(ggrepel)
  library(FNN)
  library(Cairo)
  library(patchwork)
})

base_path <- normalizePath(".")
setwd(base_path)
# Restaurer %in% base après chargement Banksy (qui le masque parfois)
`%in%` <- base::`%in%`
# Forcer Assays de Seurat (Banksy peut le masquer)
Assays <- SeuratObject::Assays
source("scripts/00_palette.R")

# ------------------------------------------------------------------
# Paramètres globaux
# ------------------------------------------------------------------
SAMPLE_ORDER  <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi", "mock_6wpi")
SAMPLE_LABELS <- c(
  mock_6wpi = "Mock 6 wpi",
  LCMV_1wpi = "LCMV 1 wpi",
  LCMV_3wpi = "LCMV 3 wpi",
  LCMV_6wpi = "LCMV 6 wpi"
)
TIMEPOINTS <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")

MICROGLIA_GLOBAL_LABEL <- "Microglia (P2ry12)"
MICROGLIA_NICHE_LABEL  <- "Microglia (C1qa)"

# Labels T cells dans l'objet global (pattern OR)
TCELL_PATTERNS <- c("Immune \\(Acod1\\)", "T cells")

FDR_CUTOFF  <- 0.05
FC_CUTOFF   <- 0.25
TOP_N_LABEL <- 15
MIN_PCT     <- 0.05
FC_THRESH   <- 0.10

DIST_BINS <- c(0, 50, 100, 200, 300, Inf)
DIST_LABELS <- c("<50µm", "50-100µm", "100-200µm", "200-300µm", ">300µm")

direction_colors <- c("up" = "#B2182B", "down" = "#2166AC", "ns" = "grey75")
niche_palette    <- c("Out niche" = "grey60", "In niche" = "#B2182B")

cond_palette <- c(
  "Mock 6 wpi"  = "grey60",
  "LCMV 1 wpi"  = "#56B4E9",
  "LCMV 3 wpi"  = "#E69F00",
  "LCMV 6 wpi"  = "#D55E00"
)

out_dir <- file.path("outputs", "banksy", "microglia_dam_niche")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
find_cl_col <- function(se, lam, res) {
  all_cols <- clusterNames(se)
  lam_str  <- gsub("\\.", "\\\\.", as.character(lam))
  lam_cols <- all_cols[grep(paste0("lam", lam_str), all_cols)]
  if (length(lam_cols) == 0) stop("Cluster column not found for lam=", lam)
  res_vals <- suppressWarnings(as.numeric(sub(".*_res", "", lam_cols)))
  idx      <- which(!is.na(res_vals) & abs(res_vals - res) < 1e-8)
  if (length(idx) == 0) stop("Cluster column not found for res=", res)
  lam_cols[idx[1]]
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

pval_stars <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "n.s.",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "n.s."
  )
}

# ==================================================================
# CHARGEMENT DES OBJETS ET RECONSTRUCTION DES ANNOTATIONS
# ==================================================================
message("\n=== Chargement des objets ===\n")

# --- Objet global (04) ---
obj_file <- file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds")
stopifnot(file.exists(obj_file))
se_global <- readRDS(obj_file)
message("Global: ", ncol(se_global), " cells")

cl_col <- find_cl_col(se_global, 0.2, 0.9)
message("  Cluster column: ", cl_col)

# Reconstruction annotations depuis CSV
csv_path <- "ncells_by_sample_lam02_res09_joint_long.csv"
stopifnot(file.exists(csv_path))
annot_long <- read_delim(csv_path, delim = ";", show_col_types = FALSE,
                          trim_ws = TRUE) %>%
  select(-matches("^Unnamed")) %>%
  mutate(banksy_domain = as.character(banksy_domain),
         annotation    = trimws(as.character(annotation)))

annotation_map <- annot_long %>%
  filter(!is.na(annotation), annotation != "") %>%
  distinct(banksy_domain, annotation)

domain_labels      <- paste0("Domain_", as.character(colData(se_global)[[cl_col]]))
anno_lookup        <- setNames(annotation_map$annotation, annotation_map$banksy_domain)
annotation_global  <- ifelse(!is.na(anno_lookup[domain_labels]) &
                               anno_lookup[domain_labels] != "",
                             anno_lookup[domain_labels], "Non annote")

message("Global annotations:\n")
print(sort(table(annotation_global), decreasing = TRUE))

# Cellules out-niche (Microglia P2ry12)
out_idx  <- which(annotation_global == MICROGLIA_GLOBAL_LABEL)
se_out   <- se_global[, out_idx]
colData(se_out)$niche_status <- "Out niche"
message("\nOut-niche [", MICROGLIA_GLOBAL_LABEL, "]: ", ncol(se_out), " cells")
print(table(as.character(colData(se_out)$sample)))

# --- Objet niche (08) ---
niche_file <- file.path("objects", "08_immune_annotated_lam02_res03.rds")
stopifnot(file.exists(niche_file))
se_niche       <- readRDS(niche_file)
cell_type_col  <- if ("cell_type" %in% colnames(colData(se_niche))) "cell_type" else
                  stop("'cell_type' column not found in se_niche colData")
cell_type_niche <- as.character(colData(se_niche)[[cell_type_col]])

if (!MICROGLIA_NICHE_LABEL %in% cell_type_niche) {
  stop("Label '", MICROGLIA_NICHE_LABEL, "' not found. Present: ",
       paste(sort(unique(cell_type_niche)), collapse = ", "))
}
in_idx   <- which(cell_type_niche == MICROGLIA_NICHE_LABEL)
se_in    <- se_niche[, in_idx]
colData(se_in)$niche_status <- "In niche"
message("\nIn-niche [", MICROGLIA_NICHE_LABEL, "]: ", ncol(se_in), " cells")
print(table(as.character(colData(se_in)$sample)))

# Coordonnées spatiales (avant conversion Seurat)
coords_out <- as.data.frame(spatialCoords(se_out))
colnames(coords_out) <- c("x", "y")
coords_out$orig_id     <- colnames(se_out)
coords_out$sample      <- as.character(colData(se_out)$sample)
coords_out$niche_status <- "Out niche"

# ==================================================================
# Construction de l'objet Seurat fusionné (LCMV uniquement)
# ==================================================================
message("\n=== Construction objet Seurat LCMV merged ===\n")

assay_out_names <- tryCatch(as.character(assayNames(se_out)), error = function(e) character(0))
assay_out <- if ("counts" %in% assay_out_names) "counts" else assay_out_names[1]
assay_in_names  <- tryCatch(as.character(assayNames(se_in)),  error = function(e) character(0))
assay_in  <- if ("counts" %in% assay_in_names)  "counts" else assay_in_names[1]
message("  assay_out='", assay_out, "' | assay_in='", assay_in, "'")

message("Converting Out-niche to Seurat...")
# Sanitiser colData : supprimer les colonnes de type liste non supportées par as.Seurat
cd_out <- colData(se_out)
bad_cols_out <- sapply(colnames(cd_out), function(x) is.list(cd_out[[x]]))
if (any(bad_cols_out)) {
  message("  Dropping list columns from se_out colData: ", paste(names(which(bad_cols_out)), collapse = ", "))
  colData(se_out) <- cd_out[, !bad_cols_out, drop = FALSE]
}
# Conversion robuste : créer Seurat directement depuis la matrice de counts
message("  Extracting counts assay...")
cnt_out_raw <- assay(se_out, assay_out)
message("  class cnt_out_raw: ", class(cnt_out_raw))
cnt_out <- methods::as(cnt_out_raw, "dgCMatrix")
message("  dim cnt_out: ", nrow(cnt_out), " x ", ncol(cnt_out))
rownames(cnt_out) <- rownames(se_out)
colnames(cnt_out) <- colnames(se_out)
message("  CreateSeuratObject out...")
so_out <- suppressWarnings(CreateSeuratObject(counts = cnt_out, project = "out_niche"))
message("  so_out created: ", ncol(so_out), " cells")
rm(cnt_out, cnt_out_raw)
da_out <- if ("RNA" %in% Assays(so_out)) "RNA" else Assays(so_out)[1]
message("  da_out: ", da_out)
DefaultAssay(so_out) <- da_out
message("  NormalizeData out...")
so_out <- NormalizeData(so_out, verbose = FALSE)
message("  NormalizeData out done")
so_out$niche_status <- "Out niche"
so_out$sample       <- as.character(colData(se_out)$sample)
message("  so_out metadata set")

message("Converting In-niche to Seurat...")
cd_in <- colData(se_in)
bad_cols_in <- sapply(colnames(cd_in), function(x) is.list(cd_in[[x]]))
if (any(bad_cols_in)) {
  message("  Dropping list columns from se_in colData: ", paste(names(which(bad_cols_in)), collapse = ", "))
  colData(se_in) <- cd_in[, !bad_cols_in, drop = FALSE]
}
cnt_in_raw <- assay(se_in, assay_in)
cnt_in <- methods::as(cnt_in_raw, "dgCMatrix")
rownames(cnt_in) <- rownames(se_in)
colnames(cnt_in) <- colnames(se_in)
so_in <- suppressWarnings(CreateSeuratObject(counts = cnt_in, project = "in_niche"))
rm(cnt_in, cnt_in_raw)
da_in <- if ("RNA" %in% Assays(so_in)) "RNA" else Assays(so_in)[1]
DefaultAssay(so_in) <- da_in
so_in <- NormalizeData(so_in, verbose = FALSE)
so_in$niche_status <- "In niche"
so_in$sample       <- as.character(colData(se_in)$sample)

common_genes <- intersect(rownames(so_out), rownames(so_in))
message("  Common genes: ", length(common_genes))
so_out_sub <- subset(so_out, features = common_genes)
so_in_sub  <- subset(so_in,  features = common_genes)

message("Merging...")
so_combined <- merge(so_out_sub, so_in_sub,
                     add.cell.ids = c("out", "in"),
                     project      = "microglia_niche")
da_comb <- if ("RNA" %in% Assays(so_combined)) "RNA" else Assays(so_combined)[1]
DefaultAssay(so_combined) <- da_comb
so_combined <- NormalizeData(so_combined, verbose = FALSE)
# Seurat v5 : JoinLayers requis après merge avant FindMarkers
so_combined <- tryCatch(JoinLayers(so_combined), error = function(e) so_combined)

# Filtrer sur LCMV uniquement
so_lcmv <- subset(so_combined, subset = sample %in% TIMEPOINTS)
message("  After LCMV filter: ", ncol(so_lcmv), " cells (", ncol(so_combined), " total)")

# ==================================================================
# SECTION 1 — Volcano DEG : tous timepoints LCMV mergés
# ==================================================================
message("\n=== SECTION 1 : Volcano DEG in vs out (tous LCMV mergés) ===\n")

Idents(so_lcmv) <- "niche_status"
deg_merged <- FindMarkers(
  so_lcmv,
  ident.1    = "In niche",
  ident.2    = "Out niche",
  test.use   = "wilcox",
  min.pct    = MIN_PCT,
  logfc.threshold = FC_THRESH,
  verbose    = TRUE
)

deg_merged$gene <- rownames(deg_merged)
deg_merged <- deg_merged %>%
  mutate(
    direction = case_when(
      avg_log2FC >  FC_CUTOFF & p_val_adj < FDR_CUTOFF ~ "up",
      avg_log2FC < -FC_CUTOFF & p_val_adj < FDR_CUTOFF ~ "down",
      TRUE ~ "ns"
    ),
    neg_log10_fdr = -log10(p_val_adj + 1e-300)
  )

write.csv(deg_merged,
          file.path(out_dir, "DEG_microglia_in_vs_out_all_timepoints_merged.csv"),
          row.names = FALSE)
message("  Saved DEG CSV")

# Top labels (symétriques, n par direction)
n_up   <- sum(deg_merged$direction == "up")
n_down <- sum(deg_merged$direction == "down")
lab_genes <- bind_rows(
  deg_merged %>% filter(direction == "up")   %>%
    slice_max(order_by = neg_log10_fdr, n = min(TOP_N_LABEL, n_up)),
  deg_merged %>% filter(direction == "down")  %>%
    slice_max(order_by = neg_log10_fdr, n = min(TOP_N_LABEL, n_down))
)

n_up_sig   <- sum(deg_merged$direction == "up")
n_down_sig <- sum(deg_merged$direction == "down")

p_volcano <- ggplot(deg_merged,
                    aes(x = avg_log2FC, y = neg_log10_fdr, color = direction)) +
  geom_point(size = 1.2, alpha = 0.7, stroke = 0) +
  geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF), linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_text_repel(data = lab_genes,
                  aes(label = gene),
                  size = 2.8, max.overlaps = 25, min.segment.length = 0.2,
                  segment.color = "grey60", segment.size = 0.3) +
  scale_color_manual(
    values = direction_colors,
    labels = c(up   = paste0("Up in niche (", n_up_sig, ")"),
               down = paste0("Down in niche (", n_down_sig, ")"),
               ns   = "n.s."),
    name   = NULL
  ) +
  labs(
    title    = "Microglia: In niche vs Out niche (all LCMV timepoints merged)",
    subtitle = sprintf("Wilcoxon | FDR < %.2f | log2FC > %.2f | n_in=%d n_out=%d",
                       FDR_CUTOFF, FC_CUTOFF,
                       sum(so_lcmv$niche_status == "In niche"),
                       sum(so_lcmv$niche_status == "Out niche")),
    x = "avg log2 FC (In / Out niche)",
    y = "-log10(FDR)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle = element_text(size = 8, color = "grey40", hjust = 0),
    legend.position = "top"
  )

save_fig(p_volcano, "fig_volcano_in_vs_out_all_timepoints_merged", 7, 6)

# ==================================================================
# SECTION 2 — Module score DAM : tous timepoints LCMV mergés
# ==================================================================
message("\n=== SECTION 2 : Module score DAM (LCMV mergés) ===\n")

sig_dir <- file.path("outputs", "banksy", "dam_signature_curation")
read_sig <- function(fname) {
  f <- file.path(sig_dir, fname)
  stopifnot(file.exists(f))
  read.csv(f, stringsAsFactors = FALSE)$gene
}

sig_raw <- list(
  Upregulated_DAM   = read_sig("Upregulated_DAM.csv"),
  Downregulated_DAM = read_sig("Downregulated_DAM.csv")
)

panel_genes    <- rownames(so_lcmv)
panel_title    <- tools::toTitleCase(tolower(panel_genes))
panel_name_map <- setNames(panel_genes, panel_title)

harm_intersect <- function(sig) {
  hits <- tools::toTitleCase(tolower(sig))
  hits <- hits[hits %in% panel_title]
  unique(as.character(panel_name_map[hits]))
}
sig_modules <- lapply(sig_raw, harm_intersect)
sig_modules  <- sig_modules[sapply(sig_modules, length) > 0]
message("Signature genes matched:")
for (nm in names(sig_modules)) message("  ", nm, ": ", length(sig_modules[[nm]]))

ctrl_n <- max(1L, min(5L, floor((nrow(so_lcmv) - max(sapply(sig_modules, length))) / 2)))
for (nm in names(sig_modules)) {
  set.seed(1997)
  so_lcmv <- AddModuleScore(so_lcmv,
                             features = list(sig_modules[[nm]]),
                             name     = nm, ctrl = ctrl_n, seed = 1997)
}

scores_df2 <- data.frame(
  cell_id      = colnames(so_lcmv),
  sample       = so_lcmv$sample,
  niche_status = factor(so_lcmv$niche_status, levels = c("Out niche", "In niche")),
  stringsAsFactors = FALSE
)
for (nm in names(sig_modules)) {
  scores_df2[[nm]] <- so_lcmv[[paste0(nm, "1")]][, 1]
}

# Violin pour chaque signature
for (nm in names(sig_modules)) {
  sig_display <- gsub("_", " ", nm)

  d_in  <- scores_df2[[nm]][scores_df2$niche_status == "In niche"]
  d_out <- scores_df2[[nm]][scores_df2$niche_status == "Out niche"]
  pval  <- tryCatch(wilcox.test(d_in, d_out, exact = FALSE)$p.value,
                    error = function(e) NA_real_)
  stars <- pval_stars(pval)

  y_max  <- max(scores_df2[[nm]], na.rm = TRUE)
  y_rng  <- diff(range(scores_df2[[nm]], na.rm = TRUE))
  y_text <- y_max + 0.08 * y_rng

  p_vio <- ggplot(scores_df2,
                  aes(x = niche_status, y = .data[[nm]], fill = niche_status)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.5, colour = NA) +
    geom_boxplot(width = 0.18, outlier.size = 0.4, outlier.alpha = 0.4,
                 colour = "grey20", fill = "white") +
    annotate("text", x = 1.5, y = y_text, label = stars,
             size = 4, fontface = "bold", colour = "grey20") +
    scale_fill_manual(values = niche_palette, guide = "none") +
    labs(
      title    = paste0("Module score: ", sig_display, " (all LCMV merged)"),
      subtitle = sprintf("Wilcoxon p %s | n_in=%d n_out=%d",
                         ifelse(is.na(pval), "NA", formatC(pval, digits = 3, format = "e")),
                         sum(scores_df2$niche_status == "In niche"),
                         sum(scores_df2$niche_status == "Out niche")),
      x = NULL, y = "Module score (AddModuleScore)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title         = element_text(face = "bold", size = 12),
      plot.subtitle      = element_text(size = 8, colour = "grey40"),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank()
    )

  sig_slug <- gsub("[^a-z0-9]", "_", tolower(nm))
  save_fig(p_vio, paste0("fig_module_score_merged_violin_", sig_slug), 5, 5)
}

# ==================================================================
# SECTION 3 — Distance Microglia → T cells (Hors niche vs In niche)
# ==================================================================
message("\n=== SECTION 3 : Distance Microglia -> T cells (Hors vs In niche) ===\n")

# Préparer le data frame global (tous échantillons)
samples_global <- as.character(colData(se_global)$sample)
spatial_global <- as.data.frame(spatialCoords(se_global))
colnames(spatial_global) <- c("x", "y")

# Identifier les labels T cells / Immune (Acod1) dans les annotations globales
tcell_mask <- grepl(paste(TCELL_PATTERNS, collapse = "|"), annotation_global)

# Hors niche: microglia P2ry12 dans l'objet global
out_mask <- annotation_global == MICROGLIA_GLOBAL_LABEL
# In niche: microglia C1qa de l'objet niche, retrouvées dans l'objet global
in_mask <- colnames(se_global) %in% colnames(se_in)

message("  T-cell / Immune (Acod1) cells across all samples: ", sum(tcell_mask))
message("  Microglia Hors niche (P2ry12): ", sum(out_mask))
message("  Microglia In niche (C1qa): ", sum(in_mask))
cat("  T-cell labels found:\n")
print(table(annotation_global[tcell_mask]))

compute_dist_by_group <- function(sample_id, mg_mask, group_label) {
  mg_idx <- which(mg_mask & samples_global == sample_id)
  tc_idx <- which(tcell_mask & samples_global == sample_id)

  if (length(mg_idx) == 0) {
    message("  ", sample_id, " | ", group_label, ": aucune microglie — ignoré")
    return(NULL)
  }
  if (length(tc_idx) == 0) {
    message("  ", sample_id, " | ", group_label, ": aucune cellule T / Immune (Acod1) — distance NA")
    return(data.frame(
      sample       = sample_id,
      cell_id      = colnames(se_global)[mg_idx],
      dist_um      = NA_real_,
      niche_status = group_label,
      stringsAsFactors = FALSE
    ))
  }

  mg_coords <- as.matrix(spatial_global[mg_idx, c("x", "y")])
  tc_coords <- as.matrix(spatial_global[tc_idx, c("x", "y")])

  knn_res <- FNN::get.knnx(data = tc_coords, query = mg_coords, k = 1)
  d_um <- as.numeric(knn_res$nn.dist[, 1])

  message(sprintf("  %-12s | %-9s  microglia: %d  T cells: %d  median dist: %.1f µm",
                  sample_id, group_label, length(mg_idx), length(tc_idx), median(d_um)))

  data.frame(
    sample       = sample_id,
    cell_id      = colnames(se_global)[mg_idx],
    dist_um      = d_um,
    niche_status = group_label,
    stringsAsFactors = FALSE
  )
}

dist_list <- lapply(SAMPLE_ORDER, function(s) {
  bind_rows(
    compute_dist_by_group(s, out_mask, "Out niche"),
    compute_dist_by_group(s, in_mask, "In niche")
  )
})

dist_df <- bind_rows(dist_list) %>%
  mutate(
    sample_label = factor(SAMPLE_LABELS[sample], levels = SAMPLE_LABELS[SAMPLE_ORDER]),
    niche_group = case_when(
      niche_status == "In niche" ~ "In niche",
      TRUE ~ "Hors niche"
    ),
    niche_group = factor(niche_group, levels = c("Hors niche", "In niche"))
  )

write.csv(dist_df,
          file.path(out_dir, "distance_microglia_to_tcell_per_cell.csv"),
          row.names = FALSE)
message("  Saved distance CSV")

# --- Figure par condition : Microglia (P2ry12) → T cells, toutes conditions ---
set.seed(1997)
dist_p2ry12 <- dist_df %>%
  filter(niche_status == "Out niche", !is.na(dist_um))

# Wilcoxon vs Mock pour chaque timepoint LCMV
wilcox_p2ry12 <- lapply(TIMEPOINTS, function(s) {
  d1 <- dist_p2ry12$dist_um[dist_p2ry12$sample == s]
  d2 <- dist_p2ry12$dist_um[dist_p2ry12$sample == "mock_6wpi"]
  pv <- tryCatch(wilcox.test(d1, d2, exact = FALSE)$p.value,
                 error = function(e) NA_real_)
  data.frame(
    sample_label = factor(SAMPLE_LABELS[s], levels = SAMPLE_LABELS[SAMPLE_ORDER]),
    stars        = pval_stars(pv),
    pval         = pv,
    stringsAsFactors = FALSE
  )
})
wilcox_p2ry12_df <- bind_rows(wilcox_p2ry12)

ymax_dist <- max(dist_p2ry12$dist_um, na.rm = TRUE)
yrng_dist <- diff(range(dist_p2ry12$dist_um, na.rm = TRUE))
ytext_dist <- ymax_dist + 0.05 * yrng_dist

niche_palette_blue_orange <- c(
  "Hors niche" = "#2C7FB8",
  "In niche" = "#F28E2B"
)

p_dist <- ggplot(dist_p2ry12,
                 aes(x = sample_label, y = dist_um, fill = sample_label)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.55, colour = NA) +
  geom_boxplot(width = 0.18, outlier.size = 0.3, outlier.alpha = 0.3,
               colour = "grey20", fill = "white") +
  geom_text(data = wilcox_p2ry12_df,
            aes(x = sample_label, y = ytext_dist, label = stars),
            inherit.aes = FALSE, size = 3.6, fontface = "bold", colour = "grey20") +
  scale_fill_manual(values = cond_palette, guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x, " \u00b5m")) +
  labs(
    title    = "Distance Microglia (P2ry12) \u2192 T cells / Immune niche",
    subtitle = "Wilcoxon vs Mock",
    x        = NULL,
    y        = "Distance to nearest T cell (\u00b5m)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 8, colour = "grey40"),
    axis.text.x      = element_text(angle = 30, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

save_fig(p_dist, "fig_distance_microglia_to_tcell_overtime", 6, 5)

# ==================================================================
# SECTION 4 — DAM score par bins de distance (Hors niche vs In niche)
# ==================================================================
message("\n=== SECTION 4 : DAM score par bins de distance (Hors vs In niche) ===\n")

if (!"Upregulated_DAM" %in% colnames(scores_df2)) {
  message("  Upregulated_DAM score not found in scores_df2; section 4 skipped.")
} else {
  dam_scores_inout <- scores_df2 %>%
    filter(sample %in% TIMEPOINTS) %>%
    mutate(
      cell_id_orig = sub("^(out|in)_", "", cell_id),
      niche_group = case_when(
        as.character(niche_status) == "In niche" ~ "In niche",
        TRUE ~ "Hors niche"
      ),
      niche_group = factor(niche_group, levels = c("Hors niche", "In niche"))
    ) %>%
    dplyr::select(cell_id_orig, sample, niche_group, dam_score = Upregulated_DAM)

  dist_sub <- dist_df %>%
    filter(sample %in% TIMEPOINTS, !is.na(dist_um)) %>%
    dplyr::select(cell_id, sample, niche_group, dist_um)

  joined_df <- inner_join(
    dist_sub,
    dam_scores_inout,
    by = c("cell_id" = "cell_id_orig", "sample", "niche_group")
  )
  message("  Cells joined (distance x DAM score): ", nrow(joined_df))

  if (nrow(joined_df) == 0) {
    message("  WARNING: no cells matched between distance and DAM score tables.")
    message("  Distance cell_ids (first 3): ", paste(head(dist_sub$cell_id, 3), collapse = ", "))
    message("  DAM score cell_ids (first 3): ", paste(head(dam_scores_inout$cell_id_orig, 3), collapse = ", "))
  }

  if (nrow(joined_df) > 0) {
    joined_df$dist_bin <- cut(joined_df$dist_um,
                              breaks = DIST_BINS,
                              labels = DIST_LABELS,
                              right = FALSE, include.lowest = TRUE)

    bin_stats <- joined_df %>%
      filter(!is.na(dist_bin)) %>%
      group_by(niche_group, dist_bin) %>%
      summarise(
        mean_dam = mean(dam_score, na.rm = TRUE),
        sem_dam = sd(dam_score, na.rm = TRUE) / sqrt(n()),
        n = n(),
        .groups = "drop"
      )

    dodge <- position_dodge(width = 0.8)

    p_bins <- ggplot(bin_stats,
                     aes(x = dist_bin, y = mean_dam,
                         fill = niche_group, group = niche_group)) +
      geom_bar(stat = "identity", position = dodge, width = 0.7, alpha = 0.85) +
      geom_errorbar(aes(ymin = mean_dam - sem_dam, ymax = mean_dam + sem_dam),
                    position = dodge, width = 0.25, linewidth = 0.5) +
      geom_text(aes(label = paste0("n=", n), y = -Inf),
                position = dodge, vjust = -0.3, size = 2.5, colour = "grey40") +
      scale_fill_manual(values = niche_palette_blue_orange, name = "Niche status") +
      labs(
        title = "Upregulated DAM score by distance to T cells",
        subtitle = paste0("LCMV merged | Hors niche vs In niche | ",
                          length(sig_modules[["Upregulated_DAM"]]), " panel genes"),
        x = "Distance to nearest T cell / Immune niche",
        y = "Mean DAM module score ± SEM"
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8, colour = "grey40"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
        panel.grid.minor = element_blank(),
        legend.position = "top"
      )

    save_fig(p_bins, "fig_dam_score_by_distance_bins", 7, 5)
  } else {
    message("  Section 4 skipped (no matched cells). Check cell ID format.")
  }
}

# ==================================================================
# SECTION 5 — Violin + barplot DAM score In niche vs Out niche
# ==================================================================
message("\n=== SECTION 5 : Violin + barplot DAM in/out niche ===\n")

if (!"Upregulated_DAM" %in% colnames(scores_df2)) {
  message("  Upregulated_DAM score not found in scores_df2; section 5 skipped.")
} else {
  dam_df <- scores_df2 %>%
    dplyr::select(niche_status, Upregulated_DAM) %>%
    dplyr::mutate(
      niche_status = factor(as.character(niche_status),
                            levels = c("Out niche", "In niche")),
      niche_group = case_when(
        as.character(niche_status) == "In niche" ~ "In niche",
        TRUE ~ "Hors niche"
      ),
      niche_group = factor(niche_group, levels = c("Hors niche", "In niche")),
      dam_positive = Upregulated_DAM > 0
    )

  niche_palette_blue_orange <- c(
    "Hors niche" = "#2C7FB8",
    "In niche"  = "#F28E2B"
  )

  # Violin DAM score by niche status
  pval_dam <- tryCatch(
    wilcox.test(
      dam_df$Upregulated_DAM[dam_df$niche_status == "In niche"],
      dam_df$Upregulated_DAM[dam_df$niche_status == "Out niche"],
      exact = FALSE
    )$p.value,
    error = function(e) NA_real_
  )
  stars_dam <- pval_stars(pval_dam)

  y_max_dam <- max(dam_df$Upregulated_DAM, na.rm = TRUE)
  y_rng_dam <- diff(range(dam_df$Upregulated_DAM, na.rm = TRUE))
  y_txt_dam <- y_max_dam + 0.08 * y_rng_dam

  p_vio_dam <- ggplot(dam_df,
                      aes(x = niche_group, y = Upregulated_DAM, fill = niche_group)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.55, colour = NA) +
    geom_boxplot(width = 0.18, outlier.size = 0.3, outlier.alpha = 0.35,
                 colour = "grey20", fill = "white") +
    annotate("text", x = 1.5, y = y_txt_dam, label = stars_dam,
             size = 4.2, fontface = "bold", colour = "grey20") +
    scale_fill_manual(values = niche_palette_blue_orange, guide = "none") +
    labs(
      title = "Upregulated DAM score: In niche vs Hors niche",
      subtitle = sprintf("Wilcoxon p %s | n_in=%d n_out=%d",
                         ifelse(is.na(pval_dam), "NA", formatC(pval_dam, digits = 3, format = "e")),
                         sum(dam_df$niche_group == "In niche"),
                         sum(dam_df$niche_group == "Hors niche")),
      x = NULL,
      y = "DAM module score (AddModuleScore)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8, colour = "grey40"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )

  # Barplot of % cells with DAM score > 0
  dam_bar <- dam_df %>%
    group_by(niche_group) %>%
    summarise(
      n_cells = n(),
      n_pos = sum(dam_positive, na.rm = TRUE),
      pct_pos = 100 * n_pos / n_cells,
      .groups = "drop"
    )

  p_bar_dam <- ggplot(dam_bar,
                      aes(x = niche_group, y = pct_pos, fill = niche_group)) +
    geom_col(width = 0.65, alpha = 0.9, colour = "grey30", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%\nn=%d", pct_pos, n_cells)),
              vjust = -0.35, size = 3) +
    scale_fill_manual(values = niche_palette_blue_orange, guide = "none") +
    expand_limits(y = max(dam_bar$pct_pos, na.rm = TRUE) * 1.15) +
    labs(
      title = "% cellules avec DAM score > 0",
      x = NULL,
      y = "% cellules"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )

  save_fig(p_vio_dam, "fig_dam_score_violin_in_vs_out_niche", 5.2, 4.8)
  save_fig(p_bar_dam, "fig_dam_score_barplot_in_vs_out_niche", 4.2, 4.8)

  p_sec5 <- p_vio_dam + p_bar_dam +
    plot_layout(widths = c(1.5, 1)) +
    plot_annotation(
      title = "DAM niche status summary",
      theme = theme(
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
      )
    )

  save_fig(p_sec5, "fig_violin_barplot_in_vs_out_niche", 9, 4.8)
}

# ==================================================================
# SECTION 6 — IN-niche microglia: distance + DAM bins (figure combinée)
# ==================================================================
message("\n=== SECTION 6 : IN-niche distance + DAM bins ===\n")

# In-niche microglia dans l'objet global.
# Priorité au label explicite MICROGLIA_NICHE_LABEL si présent.
in_mask_global <- annotation_global == MICROGLIA_NICHE_LABEL
if (!any(in_mask_global)) {
  message("  Label '", MICROGLIA_NICHE_LABEL, "' absent de l'annotation globale; fallback via IDs de se_in.")
  in_mask_global <- colnames(se_global) %in% colnames(se_in)
}

# Exclure les C1qa microglia du pool T cells (sinon distance = 0 car elles sont
# dans le domaine "Immune (Acod1)" de l'objet global)
tcell_mask_sec6 <- grepl(paste(TCELL_PATTERNS, collapse = "|"), annotation_global) &
                   !in_mask_global
samples_global_sec6 <- as.character(colData(se_global)$sample)
coords_global_sec6 <- as.data.frame(spatialCoords(se_global))
colnames(coords_global_sec6) <- c("x", "y")

dist_in_list <- lapply(SAMPLE_ORDER, function(s) {
  mg_idx <- which(in_mask_global & samples_global_sec6 == s)
  tc_idx <- which(tcell_mask_sec6 & samples_global_sec6 == s)

  if (length(mg_idx) == 0) {
    message("  ", s, ": aucune microglie IN-niche — ignoré")
    return(NULL)
  }
  if (length(tc_idx) == 0) {
    message("  ", s, ": aucune cellule T / Immune (Acod1) — distance NA")
    return(data.frame(
      sample = s,
      cell_id = colnames(se_global)[mg_idx],
      dist_um = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  mg_coords <- as.matrix(coords_global_sec6[mg_idx, c("x", "y")])
  tc_coords <- as.matrix(coords_global_sec6[tc_idx, c("x", "y")])
  knn_res <- FNN::get.knnx(data = tc_coords, query = mg_coords, k = 1)
  d_um <- as.numeric(knn_res$nn.dist[, 1])

  message(sprintf("  %-12s  in-niche microglia: %d  T cells: %d  median dist: %.1f µm",
                  s, length(mg_idx), length(tc_idx), median(d_um, na.rm = TRUE)))

  data.frame(
    sample = s,
    cell_id = colnames(se_global)[mg_idx],
    dist_um = d_um,
    stringsAsFactors = FALSE
  )
})

dist_in_df <- bind_rows(dist_in_list)
dist_in_df <- dist_in_df %>%
  mutate(sample_label = factor(SAMPLE_LABELS[sample], levels = SAMPLE_LABELS[SAMPLE_ORDER]))

wilcox_in <- lapply(TIMEPOINTS, function(s) {
  d1 <- dist_in_df$dist_um[dist_in_df$sample == s & !is.na(dist_in_df$dist_um)]
  d2 <- dist_in_df$dist_um[dist_in_df$sample == "mock_6wpi" & !is.na(dist_in_df$dist_um)]
  pv <- tryCatch(wilcox.test(d1, d2, exact = FALSE)$p.value,
                 error = function(e) NA_real_)
  data.frame(
    sample_label = factor(SAMPLE_LABELS[s], levels = SAMPLE_LABELS[SAMPLE_ORDER]),
    stars = pval_stars(pv),
    pval = pv,
    stringsAsFactors = FALSE
  )
})
wilcox_in_df <- bind_rows(wilcox_in)

ymax_in <- max(dist_in_df$dist_um, na.rm = TRUE)
yrng_in <- diff(range(dist_in_df$dist_um, na.rm = TRUE))
ytext_in <- ymax_in + 0.05 * yrng_in

p_in_violin <- ggplot(dist_in_df %>% filter(!is.na(dist_um)),
                      aes(x = sample_label, y = dist_um, fill = sample_label)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.55, colour = NA) +
  geom_boxplot(width = 0.18, outlier.size = 0.3, outlier.alpha = 0.3,
               colour = "grey20", fill = "white") +
  geom_text(data = wilcox_in_df,
            aes(x = sample_label, y = ytext_in, label = stars),
            inherit.aes = FALSE, size = 3.6, fontface = "bold", colour = "grey20") +
  scale_fill_manual(values = cond_palette, guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x, " µm")) +
  labs(
    title    = "Distance Microglia (C1qa) → T cells / Immune niche",
    subtitle = "Wilcoxon vs Mock",
    x        = NULL,
    y        = "Distance to nearest T cell (µm)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, colour = "grey40"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

if (!"Upregulated_DAM" %in% colnames(scores_df2)) {
  stop("Upregulated_DAM score not found in scores_df2 for Section 6.")
}

dam_in_df <- scores_df2 %>%
  filter(niche_status == "In niche", sample %in% c("LCMV_1wpi", "LCMV_6wpi")) %>%
  mutate(
    cell_id_orig = sub("^(out|in)_", "", cell_id),
    sample_label = factor(SAMPLE_LABELS[sample], levels = c("LCMV 1 wpi", "LCMV 6 wpi"))
  ) %>%
  dplyr::select(cell_id_orig, sample, sample_label, dam_score = Upregulated_DAM)

dist_in_1_6 <- dist_in_df %>%
  filter(sample %in% c("LCMV_1wpi", "LCMV_6wpi"), !is.na(dist_um)) %>%
  dplyr::select(cell_id, sample, dist_um)

joined_in <- inner_join(dam_in_df, dist_in_1_6,
                        by = c("cell_id_orig" = "cell_id", "sample"))

joined_in$dist_bin <- cut(joined_in$dist_um,
                          breaks = DIST_BINS,
                          labels = DIST_LABELS,
                          right = FALSE, include.lowest = TRUE)

bin_stats_in <- joined_in %>%
  filter(!is.na(dist_bin)) %>%
  group_by(sample_label, dist_bin) %>%
  summarise(
    mean_dam = mean(dam_score, na.rm = TRUE),
    sem_dam = sd(dam_score, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

dodge_sec6 <- position_dodge(width = 0.8)

p_in_bins <- ggplot(bin_stats_in,
                    aes(x = dist_bin, y = mean_dam,
                        fill = sample_label, group = sample_label)) +
  geom_bar(stat = "identity", position = dodge_sec6, width = 0.7, alpha = 0.85) +
  geom_errorbar(aes(ymin = mean_dam - sem_dam, ymax = mean_dam + sem_dam),
                position = dodge_sec6, width = 0.25, linewidth = 0.5) +
  geom_text(aes(label = paste0("n=", n), y = -Inf),
            position = dodge_sec6, vjust = -0.3, size = 2.5, colour = "grey40") +
  scale_fill_manual(
    values = c("LCMV 1 wpi" = "#56B4E9", "LCMV 6 wpi" = "#D55E00"),
    name = "Time point"
  ) +
  labs(
    title = "IN-niche DAM score by distance bin",
    subtitle = paste0("Upregulated_DAM | ", length(sig_modules[["Upregulated_DAM"]]), " panel genes"),
    x = "Distance to nearest T cell / Immune niche",
    y = "Mean DAM module score ± SEM"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, colour = "grey40"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

save_fig(p_in_violin, "fig_distance_microglia_INNICHE_to_tcell", 6, 5)

p_sec6 <- p_in_violin + p_in_bins +
  plot_layout(widths = c(1.1, 1.4)) +
  plot_annotation(
    title = "IN-niche microglia distance and DAM summary",
    theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
  )

save_fig(p_sec6, "fig_combined_inniche_distance_dam", 12, 5)

message("\n=== Script 39 terminé avec succès ===\n")
