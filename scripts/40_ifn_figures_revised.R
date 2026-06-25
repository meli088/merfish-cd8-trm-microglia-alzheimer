#!/usr/bin/env Rscript
# ============================================================================
# Script 40 — IFN figures revised (objet global)
#
# SECTION 1 : Dotplot Ifng par annotation (tous types cellulaires,
#             tous échantillons mergés, objet global 04)
#             Output : outputs/banksy/ifn_immune_overlay/
#                      fig_ifng_dotplot_all_celltypes.pdf/jpg
#
# SECTION 2 : Volcano Ifng+ vs Ifng- (tous types confondus, mergé)
#             FindMarkers : quels gènes distinguent les cellules Ifng+ ?
#             Output : outputs/banksy/ifn_immune_overlay/
#                      fig_volcano_ifng_pos_vs_neg_merged.pdf/jpg
#
# SECTION 3 : XY spatial map de l'échantillon le plus riche en Ifng+
#             + Lineplot % Ifng+ dans Immune (Acod1) overtime
#             Output : outputs/banksy/ifn_immune_overlay/
#                      fig_ifng_spatial_focus.pdf/jpg
#                      fig_ifng_tcell_overtime.pdf/jpg
#
# Inputs :
#   objects/04_banksy_joint_lam08_after_bloc3.rds
#   ncells_by_sample_lam02_res09_joint_long.csv
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
  library(scales)
  library(Cairo)
})

base_path <- normalizePath(".")
setwd(base_path)
source("scripts/00_palette.R")

# ------------------------------------------------------------------
# Paramètres
# ------------------------------------------------------------------
SAMPLE_ORDER  <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi", "mock_6wpi")
SAMPLE_LABELS <- c(
  mock_6wpi = "Mock 6 wpi",
  LCMV_1wpi = "LCMV 1 wpi",
  LCMV_3wpi = "LCMV 3 wpi",
  LCMV_6wpi = "LCMV 6 wpi"
)
IMMUNE_ACOD1_LABEL <- "Immune (Acod1)"
FDR_CUTOFF         <- 0.05
FC_CUTOFF          <- 0.25
TOP_N_LABEL        <- 15
MIN_PCT            <- 0.05
FC_THRESH          <- 0.10

direction_colors <- c("up" = "#B2182B", "down" = "#2166AC", "ns" = "grey75")

out_dir <- file.path("outputs", "banksy", "ifn_immune_overlay")
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

# ==================================================================
# Chargement de l'objet global + reconstruction annotations
# ==================================================================
message("\n=== Chargement de l'objet global ===\n")

obj_file <- file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds")
stopifnot(file.exists(obj_file))
se <- readRDS(obj_file)
message("Cells: ", ncol(se))

cl_col <- find_cl_col(se, 0.2, 0.9)
message("Cluster column: ", cl_col)

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

domain_labels     <- paste0("Domain_", as.character(colData(se)[[cl_col]]))
anno_lookup       <- setNames(annotation_map$annotation, annotation_map$banksy_domain)
annotation_global <- ifelse(!is.na(anno_lookup[domain_labels]) &
                              anno_lookup[domain_labels] != "",
                            anno_lookup[domain_labels], "Non annote")

samples_vec   <- as.character(colData(se)$sample)
spatial_mat   <- as.data.frame(spatialCoords(se))
colnames(spatial_mat) <- c("x", "y")

# Extraction de l'expression d'Ifng
avail_assays <- assayNames(se)
assay_use <- if ("logcounts" %in% avail_assays) "logcounts" else
             if ("normcounts" %in% avail_assays) "normcounts" else avail_assays[1]
message("Assay: ", assay_use)

expr_mat <- assay(se, assay_use)

if (!"Ifng" %in% rownames(expr_mat)) {
  stop("Ifng not found in rownames. Available (first 20): ",
       paste(head(rownames(expr_mat), 20), collapse = ", "))
}

ifng_expr <- as.numeric(expr_mat["Ifng", ])
ifng_pos  <- ifng_expr > 0
message("Ifng+ cells: ", sum(ifng_pos), " / ", length(ifng_pos),
        sprintf(" (%.2f%%)", 100 * mean(ifng_pos)))

# Data frame de base
base_df <- data.frame(
  x          = spatial_mat$x,
  y          = spatial_mat$y,
  sample     = samples_vec,
  annotation = annotation_global,
  ifng_expr  = ifng_expr,
  ifng_pos   = ifng_pos,
  stringsAsFactors = FALSE
)

base_df$annotation <- dplyr::recode(
  base_df$annotation,
  "Prolif neural/glial (Ccdc153)" = "Ependymal (Ccdc153)"
)

message("Annotations présentes:\n")
print(sort(table(base_df$annotation), decreasing = TRUE))

# ==================================================================
# SECTION 1 — Dotplot Ifng par annotation (tous échantillons mergés)
# ==================================================================
message("\n=== SECTION 1 : Dotplot Ifng par annotation ===\n")

ann_order <- order_annotations(unique(base_df$annotation), extended = TRUE)
ann_order <- ann_order[ann_order %in% unique(base_df$annotation)]

dot_data <- do.call(rbind, lapply(ann_order, function(ann) {
  idx <- base_df$annotation == ann
  data.frame(
    annotation = ann,
    mean_expr  = mean(base_df$ifng_expr[idx], na.rm = TRUE),
    pct_pos    = mean(base_df$ifng_pos[idx],  na.rm = TRUE) * 100,
    n_cells    = sum(idx),
    stringsAsFactors = FALSE
  )
}))

# Plafonnement à q99 pour la couleur
cap_val <- quantile(dot_data$mean_expr, 0.99, na.rm = TRUE)
if (cap_val == 0) cap_val <- max(dot_data$mean_expr, na.rm = TRUE)
dot_data$mean_expr_capped <- pmin(dot_data$mean_expr, cap_val)
dot_data$annotation <- factor(dot_data$annotation, levels = rev(ann_order))

p_dot <- ggplot(dot_data, aes(x = "Ifng", y = annotation)) +
  geom_point(aes(size = pct_pos, color = mean_expr_capped)) +
  scale_color_gradient(
    low  = "grey90",
    high = "#CC0000",
    name = "Mean expr.\n(logcounts, q99 cap)"
  ) +
  scale_size_continuous(range = c(0.5, 8), name = "% cells\nIfng > 0") +
  labs(
    title    = "Ifng expression by cell type (all samples merged)",
    subtitle = sprintf("Global BANKSY annotations | assay: %s | n = %d cells",
                       assay_use, nrow(base_df)),
    x = NULL, y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle    = element_text(size = 8, colour = "grey40"),
    axis.text.y      = element_text(size = 9),
    axis.text.x      = element_text(size = 10),
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
  )

n_ann  <- length(ann_order)
fig_h  <- max(4, n_ann * 0.32 + 2)

save_fig(p_dot, "fig_ifng_dotplot_all_celltypes", 6, fig_h)

# SECTION 1b — Dotplot Ifng par annotation et timepoint (LCMV uniquement)
message("\n=== SECTION 1b : Dotplot Ifng par annotation et timepoint (1/3/6 wpi) ===\n")

tp_order <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")

# Same normalization strategy as bubble plots:
# 1) divide by per-cell detected features
# 2) rescale by median nFeature within each sample
if ("nFeature_RNA" %in% colnames(as.data.frame(colData(se)))) {
  nfeature_vec <- as.numeric(colData(se)$nFeature_RNA)
} else {
  counts_assay <- if ("counts" %in% avail_assays) "counts" else assay_use
  counts_mat <- assay(se, counts_assay)
  nfeature_vec <- Matrix::colSums(counts_mat > 0)
}

depth_df <- data.frame(
  sample = samples_vec,
  nfeature = nfeature_vec,
  stringsAsFactors = FALSE
) %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    med_nfeature = median(nfeature[nfeature > 0], na.rm = TRUE),
    .groups = "drop"
  )

norm_tp_df <- base_df %>%
  dplyr::mutate(nfeature = nfeature_vec) %>%
  dplyr::left_join(depth_df, by = "sample") %>%
  dplyr::mutate(
    ifng_norm = dplyr::if_else(
      is.na(nfeature) | nfeature <= 0 | is.na(med_nfeature) | med_nfeature <= 0,
      NA_real_,
      ifng_expr / nfeature * med_nfeature
    )
  )

dot_tp <- norm_tp_df %>%
  dplyr::filter(sample %in% tp_order) %>%
  dplyr::group_by(annotation, sample) %>%
  dplyr::summarise(
    median_expr_pos_norm = ifelse(
      sum(ifng_pos, na.rm = TRUE) > 0,
      median(ifng_norm[ifng_pos], na.rm = TRUE),
      NA_real_
    ),
    pct_pos   = mean(ifng_pos,  na.rm = TRUE) * 100,
    n_cells   = dplyr::n(),
    .groups = "drop"
  ) %>%
  tidyr::complete(
    annotation = ann_order,
    sample = tp_order,
    fill = list(median_expr_pos_norm = NA_real_, pct_pos = 0, n_cells = 0)
  ) %>%
  dplyr::mutate(
    annotation = factor(annotation, levels = rev(ann_order)),
    sample = factor(sample, levels = tp_order, labels = SAMPLE_LABELS[tp_order])
  )

cap_tp <- quantile(dot_tp$median_expr_pos_norm, 0.90, na.rm = TRUE)
if (!is.finite(cap_tp) || cap_tp == 0) cap_tp <- max(dot_tp$median_expr_pos_norm, na.rm = TRUE)
if (!is.finite(cap_tp)) cap_tp <- 1
dot_tp$median_expr_capped <- pmin(dot_tp$median_expr_pos_norm, cap_tp)

p_dot_tp <- ggplot() +
  geom_point(
    data = dot_tp %>% dplyr::filter(pct_pos == 0),
    aes(x = sample, y = annotation),
    shape = 4, size = 1.8, color = "grey75", stroke = 0.7
  ) +
  geom_point(
    data = dot_tp %>% dplyr::filter(pct_pos > 0),
    aes(x = sample, y = annotation, size = pct_pos, color = median_expr_capped),
    alpha = 0.95
  ) +
  scale_color_gradient(
    low  = "grey90",
    high = "#CC0000",
    name = "Median expr.\n(norm among Ifng+, q90 cap)"
  ) +
  scale_size_continuous(range = c(0.5, 8), name = "% cells\nIfng > 0") +
  labs(
    title = "Ifng expression by cell type (LCMV time points)",
    subtitle = "1wpi, 3wpi, 6wpi | normalized by per-cell nFeature and sample median depth",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle    = element_text(size = 8, colour = "grey40"),
    axis.text.y      = element_text(size = 9),
    axis.text.x      = element_text(size = 10),
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
  )

save_fig(p_dot_tp, "fig_ifng_dotplot_all_celltypes_timepoints", 7.2, fig_h)

# ==================================================================
# SECTION 2 — Volcano Ifng+ vs Ifng- (tous types cellulaires, mergé)
# ==================================================================
message("\n=== SECTION 2 : Volcano Ifng+ vs Ifng- ===\n")

# Conversion en Seurat (counts bruts)
assay_counts <- if ("counts" %in% avail_assays) "counts" else avail_assays[1]
message("  Using counts assay: ", assay_counts, " for Seurat conversion")

so <- suppressWarnings(as.Seurat(se, counts = assay_counts, data = NULL))
da <- if ("RNA" %in% Assays(so)) "RNA" else Assays(so)[1]
DefaultAssay(so) <- da
so <- NormalizeData(so, verbose = FALSE)

# AddMetaData avec data.frame nommé (rownames = colnames(so))
so <- AddMetaData(so, metadata = data.frame(
  ifng_status = ifelse(ifng_pos, "Ifng_pos", "Ifng_neg"),
  annotation  = annotation_global,
  sample      = samples_vec,
  row.names   = colnames(so),
  stringsAsFactors = FALSE
))

n_pos <- sum(so$ifng_status == "Ifng_pos")
n_neg <- sum(so$ifng_status == "Ifng_neg")
message("  Ifng_pos: ", n_pos, "  Ifng_neg: ", n_neg)

if (n_pos < 10) {
  message("  WARNING: fewer than 10 Ifng+ cells — volcano skipped.")
} else {
  Idents(so) <- "ifng_status"
  deg_ifng <- FindMarkers(
    so,
    ident.1         = "Ifng_pos",
    ident.2         = "Ifng_neg",
    test.use        = "wilcox",
    min.pct         = MIN_PCT,
    logfc.threshold = FC_THRESH,
    verbose         = TRUE
  )
  deg_ifng$gene <- rownames(deg_ifng)
  deg_ifng <- deg_ifng %>%
    mutate(
      direction = case_when(
        avg_log2FC >  FC_CUTOFF & p_val_adj < FDR_CUTOFF ~ "up",
        avg_log2FC < -FC_CUTOFF & p_val_adj < FDR_CUTOFF ~ "down",
        TRUE ~ "ns"
      ),
      neg_log10_fdr = -log10(p_val_adj + 1e-300)
    )

  write.csv(deg_ifng,
            file.path(out_dir, "DEG_ifng_pos_vs_neg_merged.csv"),
            row.names = FALSE)

  n_up   <- sum(deg_ifng$direction == "up")
  n_down <- sum(deg_ifng$direction == "down")

  lab_genes <- bind_rows(
    deg_ifng %>% filter(direction == "up")  %>%
      slice_max(order_by = neg_log10_fdr, n = min(TOP_N_LABEL, n_up)),
    deg_ifng %>% filter(direction == "down") %>%
      slice_max(order_by = neg_log10_fdr, n = min(TOP_N_LABEL, n_down))
  )

  p_volc <- ggplot(deg_ifng,
                   aes(x = avg_log2FC, y = neg_log10_fdr, color = direction)) +
    geom_point(size = 1.2, alpha = 0.7, stroke = 0) +
    geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF), linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    geom_text_repel(data = lab_genes, aes(label = gene),
                    size = 2.8, max.overlaps = 25, min.segment.length = 0.2,
                    segment.color = "grey60", segment.size = 0.3) +
    scale_color_manual(
      values = direction_colors,
      labels = c(up   = paste0("Up in Ifng+ (", n_up, ")"),
                 down = paste0("Down in Ifng+ (", n_down, ")"),
                 ns   = "n.s."),
      name   = NULL
    ) +
    labs(
      title    = "Which markers distinguish Ifng+ cells? (all annotations merged)",
      subtitle = sprintf("Wilcoxon | FDR < %.2f | log2FC > %.2f | n_pos=%d n_neg=%d",
                         FDR_CUTOFF, FC_CUTOFF, n_pos, n_neg),
      x = "avg log2 FC (Ifng+ / Ifng-)",
      y = "-log10(FDR)"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 12, hjust = 0),
      plot.subtitle   = element_text(size = 8, color = "grey40", hjust = 0),
      legend.position = "top"
    )

  save_fig(p_volc, "fig_volcano_ifng_pos_vs_neg_merged", 7, 6)
}

# ==================================================================
# SECTION 3A — XY spatial map : échantillon avec le + de cellules Ifng+
# ==================================================================
message("\n=== SECTION 3A : XY spatial focus (échantillon Ifng+) ===\n")

# % Ifng+ par échantillon
pct_by_sample <- base_df %>%
  group_by(sample) %>%
  summarise(pct_pos = mean(ifng_pos) * 100,
            n_pos   = sum(ifng_pos),
            n_total = n(),
            .groups = "drop") %>%
  arrange(desc(pct_pos))

message("% Ifng+ par échantillon:\n")
print(pct_by_sample)

focus_sample <- pct_by_sample$sample[1]
message("  Echantillon focus: ", focus_sample,
        sprintf(" (%.2f%% Ifng+, n=%d)", pct_by_sample$pct_pos[1], pct_by_sample$n_pos[1]))

df_focus <- base_df %>% filter(sample == focus_sample)
cap_focus <- quantile(df_focus$ifng_expr, 0.99, na.rm = TRUE)
if (cap_focus == 0) cap_focus <- max(df_focus$ifng_expr[df_focus$ifng_expr > 0], na.rm = TRUE)
df_focus$ifng_capped <- pmin(df_focus$ifng_expr, cap_focus)

# Ordonner par expression croissante (cellules exprimantes au 1er plan)
df_focus <- df_focus[order(df_focus$ifng_capped), ]

pt_size_f <- if (nrow(df_focus) > 20000) 0.2 else if (nrow(df_focus) > 5000) 0.4 else 0.6

p_spatial <- ggplot(df_focus, aes(x = x, y = y, color = ifng_capped)) +
  geom_point(size = pt_size_f, stroke = 0, alpha = 0.85) +
  scale_color_gradient(
    low   = "grey92",
    high  = "#CC0000",
    name  = sprintf("Ifng\n(q99 clip: %.2f)", cap_focus),
    limits = c(0, cap_focus),
    oob   = scales::squish
  ) +
  coord_fixed() +
  labs(
    title    = paste0("Ifng spatial expression — ", SAMPLE_LABELS[focus_sample]),
    subtitle = sprintf("%d Ifng+ cells (%.2f%% of %d total)",
                       pct_by_sample$n_pos[1], pct_by_sample$pct_pos[1],
                       pct_by_sample$n_total[1]),
    x = "X (µm)", y = "Y (µm)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title   = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 8, colour = "grey40", hjust = 0.5),
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8),
    axis.text    = element_text(size = 7),
    panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.4)
  )

save_fig(p_spatial, "fig_ifng_spatial_focus", 7, 7)

# ==================================================================
# SECTION 3B — Lineplot % Ifng+ dans Immune (Acod1) overtime
# ==================================================================
message("\n=== SECTION 3B : Lineplot % Ifng+ dans Immune (Acod1) overtime ===\n")

# Restriction aux cellules Immune (Acod1)
immune_df <- base_df %>%
  filter(annotation == IMMUNE_ACOD1_LABEL)

if (nrow(immune_df) == 0) {
  message("  WARNING: aucune cellule '", IMMUNE_ACOD1_LABEL,
          "' trouvée — section 3B ignorée.")
  message("  Annotations présentes: ",
          paste(sort(unique(base_df$annotation)), collapse = ", "))
} else {
  pct_immune <- immune_df %>%
    group_by(sample) %>%
    summarise(
      n_total  = n(),
      n_pos    = sum(ifng_pos),
      pct_pos  = mean(ifng_pos) * 100,
      .groups  = "drop"
    )

  # Compléter avec les échantillons manquants (0 cellules)
  pct_immune <- pct_immune %>%
    complete(sample = SAMPLE_ORDER, fill = list(n_total = 0, n_pos = 0, pct_pos = 0))

  pct_immune$sample_label <- factor(SAMPLE_LABELS[pct_immune$sample],
                                     levels = SAMPLE_LABELS[SAMPLE_ORDER])

  message("% Ifng+ dans Immune (Acod1) par échantillon:\n")
  print(pct_immune)

  # Palette conditions
  cond_palette <- c(
    "Mock 6 wpi" = "grey60",
    "LCMV 1 wpi" = "#56B4E9",
    "LCMV 3 wpi" = "#E69F00",
    "LCMV 6 wpi" = "#D55E00"
  )

  # Convertir en position numérique pour le lineplot (style script 16)
  # mock = 9 wpi (position référence à droite), LCMV 1/3/6 sur axe gauche
  xpos_map <- c(
    "LCMV 1 wpi" = 1,
    "LCMV 3 wpi" = 3,
    "LCMV 6 wpi" = 6,
    "Mock 6 wpi" = 9
  )
  pct_immune$x_pos <- xpos_map[as.character(pct_immune$sample_label)]

  p_line <- ggplot(pct_immune,
                   aes(x = x_pos, y = pct_pos,
                       color = sample_label, group = 1)) +
    geom_line(color = "grey50", linewidth = 0.7) +
    geom_point(aes(color = sample_label, size = n_total), stroke = 0) +
    geom_text(aes(label = sprintf("n=%d", n_total)),
              vjust = -1.2, size = 2.8, color = "grey40") +
    scale_color_manual(values = cond_palette, name = "Condition") +
    scale_size_continuous(range = c(2, 6), guide = "none") +
    scale_x_continuous(
      breaks = c(1, 3, 6, 9),
      labels = c("1 wpi", "3 wpi", "6 wpi", "mock\n6 wpi")
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
    annotate("segment", x = 7.5, xend = 7.5, y = -Inf, yend = Inf,
             linetype = "dashed", color = "grey60", linewidth = 0.4) +
    annotate("text", x = 4.0, y = Inf, label = "rLCMV",
             vjust = 1.5, size = 3.5, color = "grey30") +
    annotate("text", x = 9.0, y = Inf, label = "Mock",
             vjust = 1.5, size = 3.5, color = "grey30") +
    labs(
      title    = paste0("% Ifng-expressing cells in ", IMMUNE_ACOD1_LABEL),
      subtitle = "Global BANKSY annotations | Ifng > 0 threshold",
      x        = "Time point",
      y        = "% Ifng+ cells"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 12, hjust = 0),
      plot.subtitle   = element_text(size = 8, colour = "grey40"),
      axis.text       = element_text(size = 10),
      legend.position = "right",
      legend.text     = element_text(size = 9)
    )

  save_fig(p_line, "fig_ifng_tcell_overtime", 6.5, 4.5)
}

# ==================================================================
# SECTION 3C — XY spatial Ifng dans la niche Immune (Acod1) uniquement
# ==================================================================
message("\n=== SECTION 3C : XY spatial Ifng — niche Immune (Acod1) ===\n")

df_niche <- base_df %>% filter(annotation == IMMUNE_ACOD1_LABEL)

if (nrow(df_niche) == 0) {
  message("  WARNING: aucune cellule '", IMMUNE_ACOD1_LABEL, "' — section 3C ignorée.")
} else {
  samples_present_niche <- unique(df_niche$sample)
  sample_levels_niche   <- c(
    SAMPLE_ORDER[SAMPLE_ORDER %in% samples_present_niche],
    setdiff(samples_present_niche, SAMPLE_ORDER)
  )
  df_niche$sample_label <- factor(
    SAMPLE_LABELS[df_niche$sample],
    levels = SAMPLE_LABELS[sample_levels_niche]
  )

  message("  Cellules Immune (Acod1): ", nrow(df_niche))
  print(table(df_niche$sample))

  # Clipping q99 global sur les cellules de la niche
  cap_niche <- quantile(df_niche$ifng_expr, 0.99, na.rm = TRUE)
  if (cap_niche == 0) {
    pos_vals  <- df_niche$ifng_expr[df_niche$ifng_expr > 0]
    cap_niche <- if (length(pos_vals) > 0) max(pos_vals) else 1
  }
  df_niche$ifng_capped <- pmin(df_niche$ifng_expr, cap_niche)

  # Tri par expression croissante -> cellules Ifng+ au premier plan
  df_niche <- df_niche[order(df_niche$ifng_capped), ]

  pt_size_n <- if (nrow(df_niche) > 20000) 0.3 else if (nrow(df_niche) > 5000) 0.6 else 1.0

  p_niche <- ggplot(df_niche, aes(x = x, y = y, color = ifng_capped)) +
    geom_point(size = pt_size_n, stroke = 0, alpha = 0.9) +
    scale_color_gradient(
      low    = "grey92",
      high   = "#CC0000",
      name   = sprintf("Ifng\n(q99: %.2f)", cap_niche),
      limits = c(0, cap_niche),
      oob    = scales::squish
    ) +
    facet_wrap(~sample_label, nrow = 1, scales = "free") +
    labs(
      title    = paste0("Ifng spatial expression — ", IMMUNE_ACOD1_LABEL, " cells only"),
      subtitle = sprintf("n = %d cells | Ifng+ = %d (%.1f%%) | q99 clip",
                         nrow(df_niche), sum(df_niche$ifng_pos),
                         100 * mean(df_niche$ifng_pos)),
      x = "X (µm)", y = "Y (µm)"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle    = element_text(size = 8, colour = "grey40", hjust = 0.5),
      strip.text       = element_text(face = "bold", size = 9),
      strip.background = element_blank(),
      axis.text        = element_text(size = 6),
      legend.title     = element_text(size = 8),
      legend.text      = element_text(size = 7),
      panel.border     = element_rect(color = "grey60", fill = NA, linewidth = 0.4),
      aspect.ratio     = 1
    )

  save_fig(p_niche, "fig_ifng_spatial_immune_niche", 14, 5)
}

message("\n=== Script 40 terminé avec succès ===\n")
