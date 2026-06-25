#!/usr/bin/env Rscript
# ============================================================================
# Script 42 — DAM figures urgent
#
# SECTION 1 : DAM module score Microglia (C1qa), 6wpi vs 1wpi (in-niche)
#             Output : outputs/banksy/microglia_dam_niche/
#                      fig_dam_score_inniche_6wpi_vs_1wpi_violin.pdf/.jpg
#
# SECTION 2 : DAM module score Mac (Ctss), 6wpi vs 1wpi
#             Output : outputs/banksy/microglia_dam_niche/
#                      fig_dam_score_mac_6wpi_vs_1wpi_violin.pdf/.jpg
#
# SECTION 3 : Figure combinée (Microglia + Macrophages)
#             Output : outputs/banksy/microglia_dam_niche/
#                      fig_dam_score_combined_6wpi_vs_1wpi.pdf/.jpg
#
# SECTION 4 : nDEG summary par cluster (6wpi vs 1wpi)
#             Input  : outputs/banksy/immune_niche_volcano_by_celltype/
#                      ndeg_summary_table.csv
#             Output : outputs/banksy/immune_niche_volcano_by_celltype/
#                      fig_ndeg_summary_by_celltype_updated.pdf/.jpg
#
# SECTION 5 : Zoom XY niche highlighted (LCMV_1wpi)
#             Output : outputs/banksy/spatial_annotations/
#                      xy_zoom_niche_highlighted_lcmv_1wpi.pdf/.jpg
#
# Seed = 1997 | source scripts/00_palette.R | cairo_pdf
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

base_path <- normalizePath(".")
setwd(base_path)

# Restaurer operators/fonctions parfois masqués
`%in%` <- base::`%in%`
Assays <- SeuratObject::Assays
source("scripts/00_palette.R")

# ------------------------------------------------------------------
# Global params
# ------------------------------------------------------------------
TIMEPOINTS <- c("LCMV_1wpi", "LCMV_6wpi")
TIMEPOINT_LABELS <- c(
  LCMV_1wpi = "LCMV 1 wpi",
  LCMV_6wpi = "LCMV 6 wpi"
)
TIMEPOINT_COLORS <- c(
  LCMV_1wpi = "#56B4E9",  # bleu
  LCMV_6wpi = "#F28E2B"   # orange
)

MICROGLIA_LABEL <- "Microglia (C1qa)"
MAC_LABEL <- "Mac (Ctss)"

obj_08 <- file.path("objects", "08_immune_annotated_lam02_res03.rds")
obj_04 <- file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds")
sig_file <- file.path("outputs", "banksy", "dam_signature_curation", "Upregulated_DAM.csv")
csv_annot <- "ncells_by_sample_lam02_res09_joint_long.csv"
ndeg_file <- file.path("outputs", "banksy", "immune_niche_volcano_by_celltype", "ndeg_summary_table.csv")

out_dam <- file.path("outputs", "banksy", "microglia_dam_niche")
out_ndeg <- file.path("outputs", "banksy", "immune_niche_volcano_by_celltype")
out_xy <- file.path("outputs", "banksy", "spatial_annotations")

if (!dir.exists(out_dam)) dir.create(out_dam, recursive = TRUE)
if (!dir.exists(out_ndeg)) dir.create(out_ndeg, recursive = TRUE)
if (!dir.exists(out_xy)) dir.create(out_xy, recursive = TRUE)

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
save_fig <- function(p, out_dir, fname, width, height) {
  pdf_file <- file.path(out_dir, paste0(fname, ".pdf"))
  jpg_file <- file.path(out_dir, paste0(fname, ".jpg"))

  cairo_pdf(pdf_file, width = width, height = height)
  print(p)
  dev.off()

  jpeg(jpg_file,
       width = round(width * 150),
       height = round(height * 150),
       res = 150,
       quality = 95)
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

read_upregulated_dam <- function(fpath) {
  stopifnot(file.exists(fpath))
  sig <- read.csv(fpath, stringsAsFactors = FALSE)
  if (!"gene" %in% colnames(sig)) {
    stop("Column 'gene' not found in ", fpath)
  }
  unique(as.character(sig$gene))
}

# Harmonize signature names to object rownames using title-case matching.
match_signature_genes <- function(sig_genes, panel_genes) {
  panel_title <- tools::toTitleCase(tolower(panel_genes))
  map_panel <- setNames(panel_genes, panel_title)
  sig_title <- tools::toTitleCase(tolower(sig_genes))
  sig_title <- sig_title[sig_title %in% panel_title]
  unique(as.character(map_panel[sig_title]))
}

# Build a Seurat object from selected cells and compute DAM score.
compute_dam_score_df <- function(se_obj, cell_type_label, timepoints, sig_genes) {
  stopifnot("cell_type" %in% colnames(colData(se_obj)))
  stopifnot("sample" %in% colnames(colData(se_obj)))

  ct_vec <- as.character(colData(se_obj)$cell_type)
  smp_vec <- as.character(colData(se_obj)$sample)

  idx <- which(ct_vec == cell_type_label & smp_vec %in% timepoints)
  if (length(idx) == 0) {
    stop("No cells found for label '", cell_type_label, "' in selected timepoints")
  }

  se_sub <- se_obj[, idx]

  assay_names <- tryCatch(as.character(assayNames(se_sub)), error = function(e) character(0))
  assay_use <- if ("counts" %in% assay_names) "counts" else assay_names[1]
  if (is.na(assay_use) || assay_use == "") {
    stop("No assay found in subset object for ", cell_type_label)
  }

  cnt_raw <- assay(se_sub, assay_use)
  cnt <- methods::as(cnt_raw, "dgCMatrix")
  rownames(cnt) <- rownames(se_sub)
  colnames(cnt) <- colnames(se_sub)

  so <- suppressWarnings(CreateSeuratObject(counts = cnt, project = "dam_urgent"))
  da <- if ("RNA" %in% Assays(so)) "RNA" else Assays(so)[1]
  DefaultAssay(so) <- da
  so <- NormalizeData(so, verbose = FALSE)
  so$sample <- as.character(colData(se_sub)$sample)
  so$cell_type <- as.character(colData(se_sub)$cell_type)

  sig_hits <- match_signature_genes(sig_genes, rownames(so))
  if (length(sig_hits) == 0) {
    stop("No Upregulated_DAM genes matched rownames for ", cell_type_label)
  }

  ctrl_n <- max(1L, min(5L, floor((nrow(so) - length(sig_hits)) / 2)))
  set.seed(1997)
  so <- AddModuleScore(
    so,
    features = list(sig_hits),
    name = "Upregulated_DAM",
    ctrl = ctrl_n,
    seed = 1997
  )

  data.frame(
    cell_id = colnames(so),
    sample = factor(so$sample, levels = timepoints),
    sample_label = factor(TIMEPOINT_LABELS[so$sample], levels = TIMEPOINT_LABELS[timepoints]),
    dam_score = so$Upregulated_DAM1,
    cell_type = cell_type_label,
    stringsAsFactors = FALSE
  )
}

make_violin_dam <- function(df, title_text) {
  d1 <- df$dam_score[df$sample == "LCMV_1wpi"]
  d2 <- df$dam_score[df$sample == "LCMV_6wpi"]
  pval <- tryCatch(wilcox.test(d2, d1, exact = FALSE)$p.value, error = function(e) NA_real_)
  stars <- pval_stars(pval)

  y_max <- max(df$dam_score, na.rm = TRUE)
  y_rng <- diff(range(df$dam_score, na.rm = TRUE))
  if (!is.finite(y_rng) || y_rng == 0) y_rng <- 1
  y_txt <- y_max + 0.08 * y_rng

  ggplot(df, aes(x = sample_label, y = dam_score, fill = sample)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.55, colour = NA) +
    geom_boxplot(width = 0.18, outlier.size = 0.35, outlier.alpha = 0.35,
                 colour = "grey20", fill = "white") +
    annotate("text", x = 1.5, y = y_txt, label = stars,
             size = 4.2, fontface = "bold", colour = "grey20") +
    scale_fill_manual(values = TIMEPOINT_COLORS, guide = "none") +
    labs(
      title = title_text,
      subtitle = sprintf("Wilcoxon 6wpi vs 1wpi: p %s | n(1wpi)=%d n(6wpi)=%d",
                         ifelse(is.na(pval), "NA", formatC(pval, digits = 3, format = "e")),
                         sum(df$sample == "LCMV_1wpi"),
                         sum(df$sample == "LCMV_6wpi")),
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
}

# ------------------------------------------------------------------
# Load shared inputs
# ------------------------------------------------------------------
message("\n=== Loading shared inputs ===")

stopifnot(file.exists(obj_08))
stopifnot(file.exists(obj_04))
stopifnot(file.exists(sig_file))
stopifnot(file.exists(ndeg_file))
stopifnot(file.exists(csv_annot))

se_08 <- readRDS(obj_08)
se_04 <- readRDS(obj_04)
up_dam_genes <- read_upregulated_dam(sig_file)
message("  Upregulated_DAM genes loaded: ", length(up_dam_genes))

# ==================================================================
# SECTION 1 — DAM in-niche Microglia (C1qa): 6wpi vs 1wpi
# ==================================================================
message("\n=== SECTION 1: DAM Microglia (C1qa) 6wpi vs 1wpi ===")

df_micro <- compute_dam_score_df(
  se_obj = se_08,
  cell_type_label = MICROGLIA_LABEL,
  timepoints = TIMEPOINTS,
  sig_genes = up_dam_genes
)

p_micro <- make_violin_dam(
  df_micro,
  "DAM module score in-niche microglia (C1qa): 6wpi vs 1wpi"
)

save_fig(
  p_micro,
  out_dir = out_dam,
  fname = "fig_dam_score_inniche_6wpi_vs_1wpi_violin",
  width = 5.4,
  height = 5
)

# ==================================================================
# SECTION 2 — DAM Macrophages (Mac Ctss): 6wpi vs 1wpi
# ==================================================================
message("\n=== SECTION 2: DAM Mac (Ctss) 6wpi vs 1wpi ===")

df_mac <- compute_dam_score_df(
  se_obj = se_08,
  cell_type_label = MAC_LABEL,
  timepoints = TIMEPOINTS,
  sig_genes = up_dam_genes
)

p_mac <- make_violin_dam(
  df_mac,
  "DAM module score macrophages (Mac Ctss): 6wpi vs 1wpi"
)

save_fig(
  p_mac,
  out_dir = out_dam,
  fname = "fig_dam_score_mac_6wpi_vs_1wpi_violin",
  width = 5.4,
  height = 5
)

# ==================================================================
# SECTION 3 — Combined figure (microglia + macrophages)
# ==================================================================
message("\n=== SECTION 3: Combined DAM figure ===")

p_combined <- p_micro + p_mac +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "DAM module score at chronic vs acute stage : microglia and macrophages",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5)
    )
  )

save_fig(
  p_combined,
  out_dir = out_dam,
  fname = "fig_dam_score_combined_6wpi_vs_1wpi",
  width = 10.8,
  height = 5
)

# ==================================================================
# SECTION 4 — nDEG summary by cluster (updated)
# ==================================================================
message("\n=== SECTION 4: nDEG summary updated ===")

ndeg <- read.csv(ndeg_file, stringsAsFactors = FALSE)

# Flexible column mapping.
ct_col <- if ("cell_type" %in% names(ndeg)) "cell_type" else names(ndeg)[grep("cell", names(ndeg), ignore.case = TRUE)[1]]
up_col <- if ("n_up" %in% names(ndeg)) "n_up" else names(ndeg)[grep("up", names(ndeg), ignore.case = TRUE)[1]]
down_col <- if ("n_down" %in% names(ndeg)) "n_down" else names(ndeg)[grep("down", names(ndeg), ignore.case = TRUE)[1]]
contrast_col <- if ("contrast_label" %in% names(ndeg)) {
  "contrast_label"
} else if ("contrast" %in% names(ndeg)) {
  "contrast"
} else {
  NA_character_
}

if (any(is.na(c(ct_col, up_col, down_col)))) {
  stop("Could not identify required columns in ndeg_summary_table.csv")
}

ndeg2 <- ndeg %>%
  mutate(
    cell_type_plot = as.character(.data[[ct_col]]),
    n_up_plot = as.numeric(.data[[up_col]]),
    n_down_plot = as.numeric(.data[[down_col]]),
    contrast_plot = if (is.na(contrast_col)) "All" else as.character(.data[[contrast_col]])
  ) %>%
  filter(!is.na(cell_type_plot)) %>%
  mutate(
    n_up_plot = ifelse(is.na(n_up_plot), 0, n_up_plot),
    n_down_plot = ifelse(is.na(n_down_plot), 0, n_down_plot),
    n_total = n_up_plot + n_down_plot
  )

# Order by total DEGs across contrasts.
ct_order <- ndeg2 %>%
  group_by(cell_type_plot) %>%
  summarise(total_all = sum(n_total, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_all)) %>%
  pull(cell_type_plot)

ndeg2$cell_type_plot <- factor(ndeg2$cell_type_plot, levels = rev(ct_order))

ndeg_long <- bind_rows(
  ndeg2 %>% transmute(cell_type_plot, contrast_plot, direction = "Up", n_deg = n_up_plot),
  ndeg2 %>% transmute(cell_type_plot, contrast_plot, direction = "Down", n_deg = n_down_plot)
)

p_ndeg <- ggplot(ndeg_long, aes(x = n_deg, y = cell_type_plot, fill = direction)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65,
           colour = "grey35", linewidth = 0.2) +
  scale_fill_manual(values = c("Up" = "#D62728", "Down" = "#2C7FB8"), name = NULL) +
  labs(
    title = "nDEG summary by cell type (updated)",
    subtitle = "Barres rouges = up | barres bleues = down",
    x = "n DEGs",
    y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, colour = "grey40"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  )

if (length(unique(ndeg_long$contrast_plot)) > 1) {
  p_ndeg <- p_ndeg + facet_wrap(~ contrast_plot)
}

fig_h <- max(5.2, 0.34 * length(unique(ndeg_long$cell_type_plot)) + 1.7)
save_fig(
  p_ndeg,
  out_dir = out_ndeg,
  fname = "fig_ndeg_summary_by_celltype_updated",
  width = 10,
  height = fig_h
)

# ==================================================================
# SECTION 5 — XY zoom niche highlighted (LCMV_1wpi)
# ==================================================================
message("\n=== SECTION 5: XY zoom niche highlighted (LCMV_1wpi) ===")

# Rebuild annotations from cluster mapping CSV.
cl_col <- find_cl_col(se_04, 0.2, 0.9)
cd04 <- as.data.frame(SummarizedExperiment::colData(se_04))
xy04 <- as.data.frame(SpatialExperiment::spatialCoords(se_04))
colnames(xy04) <- c("x", "y")

anno_long <- read.delim(csv_annot, sep = ";", stringsAsFactors = FALSE) %>%
  dplyr::filter(annotation != "" & !is.na(annotation)) %>%
  dplyr::select(banksy_domain, annotation) %>%
  dplyr::distinct()
anno_long$cluster_id <- as.numeric(gsub("Domain_", "", anno_long$banksy_domain))
anno_lookup <- setNames(trimws(anno_long$annotation), anno_long$cluster_id)

clusters04 <- as.numeric(cd04[[cl_col]])
annot04 <- anno_lookup[as.character(clusters04)]
annot04[is.na(annot04)] <- "Non annote"

plot_df <- data.frame(
  x = xy04$x,
  y = xy04$y,
  sample = as.character(cd04$sample),
  annotation = annot04,
  stringsAsFactors = FALSE
)

df_1wpi <- plot_df %>% dplyr::filter(sample == "LCMV_1wpi")
if (nrow(df_1wpi) == 0) stop("No LCMV_1wpi cells found in object 04")

immune_1wpi <- df_1wpi %>% dplyr::filter(annotation == "Immune (Acod1)")
if (nrow(immune_1wpi) == 0) stop("No Immune (Acod1) cells found in LCMV_1wpi")

cx <- mean(immune_1wpi$x, na.rm = TRUE)
cy <- mean(immune_1wpi$y, na.rm = TRUE)
half_window <- 800
xlim_zoom <- c(cx - half_window, cx + half_window)
ylim_zoom <- c(cy - half_window, cy + half_window)

df_zoom <- df_1wpi %>%
  dplyr::filter(
    x >= xlim_zoom[1], x <= xlim_zoom[2],
    y >= ylim_zoom[1], y <= ylim_zoom[2]
  )

# Group into 5 highlighted classes + background.
is_neuron <- grepl("Neurons|Excitatory|Inhibitory|Neuron", df_zoom$annotation)
is_astro <- grepl("Astrocytes", df_zoom$annotation)

df_zoom$group <- dplyr::case_when(
  df_zoom$annotation == "Immune (Acod1)" ~ "Immune (Acod1)",
  df_zoom$annotation == "IFN responsive (Ifit1)" ~ "IFN responsive",
  df_zoom$annotation == "Microglia (P2ry12)" ~ "Microglia (P2ry12)",
  is_neuron ~ "Neurons",
  is_astro ~ "Astrocytes",
  TRUE ~ "Other"
)

group_levels <- c(
  "Immune (Acod1)",
  "IFN responsive",
  "Microglia (P2ry12)",
  "Neurons",
  "Astrocytes",
  "Other"
)

df_zoom$group <- factor(df_zoom$group, levels = group_levels)

group_pal <- c(
  "Immune (Acod1)" = "#D62728",  # red
  "IFN responsive" = "#F28E2B",   # orange
  "Microglia (P2ry12)" = "#2C7FB8", # blue
  "Neurons" = "#D9D9D9",          # light gray
  "Astrocytes" = "#4DAF4A",       # green
  "Other" = "#EEEEEE"             # very pale gray
)

# Plot background first, then key populations.
df_other <- df_zoom %>% dplyr::filter(group == "Other")
df_key <- df_zoom %>% dplyr::filter(group != "Other")

p_xy <- ggplot() +
  geom_point(
    data = df_other,
    aes(x = x, y = y),
    color = group_pal[["Other"]],
    size = 0.3,
    alpha = 0.2,
    stroke = 0
  ) +
  geom_point(
    data = df_key,
    aes(x = x, y = y, color = group),
    size = 1.5,
    alpha = 0.95,
    stroke = 0
  ) +
  scale_color_manual(values = group_pal[group_levels[group_levels != "Other"]], drop = FALSE) +
  coord_fixed(xlim = xlim_zoom, ylim = ylim_zoom, expand = FALSE) +
  labs(
    title = "XY zoom niche highlighted — LCMV 1wpi",
    subtitle = sprintf("Fenêtre ±%dµm autour du centroïde Immune (Acod1)", half_window),
    x = "X (µm)",
    y = "Y (µm)",
    color = "Population"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 8.5, hjust = 0.5, colour = "grey40"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.4)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

save_fig(
  p_xy,
  out_dir = out_xy,
  fname = "xy_zoom_niche_highlighted_lcmv_1wpi",
  width = 7.4,
  height = 7
)

message("\n=== Script 42 terminé avec succès ===\n")
#!/usr/bin/env Rscript
# ============================================================================
# Script 42 — DAM figures urgent
#
# SECTION 1 : DAM module score microglia in-niche (Microglia (C1qa))
#             6wpi vs 1wpi (violin + Wilcoxon)
#             Output : outputs/banksy/microglia_dam_niche/
#                      fig_dam_score_inniche_6wpi_vs_1wpi_violin.pdf/jpg
#
# SECTION 2 : DAM module score macrophages (Mac (Ctss))
#             6wpi vs 1wpi (violin + Wilcoxon)
#             Output : outputs/banksy/microglia_dam_niche/
#                      fig_dam_score_mac_6wpi_vs_1wpi_violin.pdf/jpg
#
# SECTION 3 : Figure combinée microglia + macrophages (patchwork)
#             Output : outputs/banksy/microglia_dam_niche/
#                      fig_dam_score_combined_6wpi_vs_1wpi.pdf/jpg
#
# SECTION 4 : nDEG summary par cluster (6wpi vs 1wpi)
#             Input  : outputs/banksy/immune_niche_volcano_by_celltype/
#                      ndeg_summary_table.csv
#             Output : outputs/banksy/immune_niche_volcano_by_celltype/
#                      fig_ndeg_summary_by_celltype_updated.pdf/jpg
#
# SECTION 5 : Zoom XY niche highlighted (LCMV_1wpi)
#             Output : outputs/banksy/spatial_annotations/
#                      xy_zoom_niche_highlighted_lcmv_1wpi.pdf/jpg
#
# Seed = 1997, cairo_pdf, source scripts/00_palette.R
# ============================================================================

set.seed(1997)

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(tidyverse)
  library(ggplot2)
  library(Cairo)
  library(patchwork)
})

base_path <- normalizePath(".")
setwd(base_path)
source("scripts/00_palette.R")

# ------------------------------------------------------------------
# Paths and constants
# ------------------------------------------------------------------
obj_immune_path <- file.path("objects", "08_immune_annotated_lam02_res03.rds")
obj_global_path <- file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds")
sig_up_path <- file.path("outputs", "banksy", "dam_signature_curation", "Upregulated_DAM.csv")

out_dam_dir <- file.path("outputs", "banksy", "microglia_dam_niche")
out_ndeg_dir <- file.path("outputs", "banksy", "immune_niche_volcano_by_celltype")
out_spatial_dir <- file.path("outputs", "banksy", "spatial_annotations")

if (!dir.exists(out_dam_dir)) dir.create(out_dam_dir, recursive = TRUE)
if (!dir.exists(out_ndeg_dir)) dir.create(out_ndeg_dir, recursive = TRUE)
if (!dir.exists(out_spatial_dir)) dir.create(out_spatial_dir, recursive = TRUE)

TIMEPOINTS_KEEP <- c("LCMV_1wpi", "LCMV_6wpi")
TP_LABELS <- c("LCMV_1wpi" = "1 wpi", "LCMV_6wpi" = "6 wpi")
TP_COLORS <- c("1 wpi" = "#56B4E9", "6 wpi" = "#F28E2B")

MICRO_LABEL <- "Microglia (C1qa)"
MAC_LABEL <- "Mac (Ctss)"

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
save_fig <- function(p, out_dir, fname, width, height) {
  cairo_pdf(file.path(out_dir, paste0(fname, ".pdf")), width = width, height = height)
  print(p)
  dev.off()

  jpeg(file.path(out_dir, paste0(fname, ".jpg")),
       width = width * 150, height = height * 150, res = 150, quality = 95)
  print(p)
  dev.off()

  message("  Saved: ", file.path(out_dir, paste0(fname, ".pdf")), " / .jpg")
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

read_upregulated_dam <- function(path_csv, panel_genes) {
  if (!file.exists(path_csv)) stop("Signature file not found: ", path_csv)
  sig_raw <- read.csv(path_csv, stringsAsFactors = FALSE)$gene

  # Harmonize case as in script 39
  panel_title <- tools::toTitleCase(tolower(panel_genes))
  panel_name_map <- setNames(panel_genes, panel_title)

  sig_title <- tools::toTitleCase(tolower(sig_raw))
  sig_title <- sig_title[sig_title %in% panel_title]
  sig_genes <- unique(as.character(panel_name_map[sig_title]))

  sig_genes
}

find_cl_col <- function(se, lam, res) {
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  cn <- colnames(cd)
  lam_pat <- paste0("lam", gsub("\\.", "\\\\.", as.character(lam)))
  res_pat <- paste0("_res", gsub("\\.", "\\\\.", as.character(res)), "$")
  cols <- cn[grepl(lam_pat, cn) & grepl(res_pat, cn)]
  if (length(cols) == 0) stop("Clustering column not found for lam=", lam, " res=", res)
  cols[1]
}

build_violin <- function(df, y_col, title_txt, subtitle_prefix) {
  d1 <- df[[y_col]][df$time_label == "6 wpi"]
  d2 <- df[[y_col]][df$time_label == "1 wpi"]
  pval <- tryCatch(wilcox.test(d1, d2, exact = FALSE)$p.value,
                   error = function(e) NA_real_)
  stars <- pval_stars(pval)

  y_max <- max(df[[y_col]], na.rm = TRUE)
  y_rng <- diff(range(df[[y_col]], na.rm = TRUE))
  y_txt <- y_max + 0.08 * y_rng

  p <- ggplot(df, aes(x = time_label, y = .data[[y_col]], fill = time_label)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.55, colour = NA) +
    geom_boxplot(width = 0.18, outlier.size = 0.3, outlier.alpha = 0.35,
                 colour = "grey20", fill = "white") +
    annotate("text", x = 1.5, y = y_txt, label = stars,
             size = 4.2, fontface = "bold", colour = "grey20") +
    scale_fill_manual(values = TP_COLORS, guide = "none") +
    labs(
      title = title_txt,
      subtitle = sprintf(
        "%s | Wilcoxon p %s | n(1wpi)=%d n(6wpi)=%d",
        subtitle_prefix,
        ifelse(is.na(pval), "NA", formatC(pval, digits = 3, format = "e")),
        sum(df$time_label == "1 wpi"),
        sum(df$time_label == "6 wpi")
      ),
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

  list(plot = p, pval = pval)
}

# ==================================================================
# SECTION 1 + 2 — DAM module score in immune object (Microglia/Mac)
# ==================================================================
message("\n=== SECTION 1+2: DAM module score microglia/macrophages (1wpi vs 6wpi) ===\n")

stopifnot(file.exists(obj_immune_path))
se_immune <- readRDS(obj_immune_path)

if (!"cell_type" %in% colnames(colData(se_immune))) {
  stop("'cell_type' column not found in immune object colData")
}
if (!"sample" %in% colnames(colData(se_immune))) {
  stop("'sample' column not found in immune object colData")
}

cell_type <- as.character(colData(se_immune)$cell_type)
sample_vec <- as.character(colData(se_immune)$sample)

idx_keep <- which(cell_type %in% c(MICRO_LABEL, MAC_LABEL) & sample_vec %in% TIMEPOINTS_KEEP)
if (length(idx_keep) == 0) {
  stop("No cells found for requested groups: ",
       paste(c(MICRO_LABEL, MAC_LABEL), collapse = " / "),
       " at 1wpi/6wpi")
}

se_sub <- se_immune[, idx_keep]

# Robust conversion to Seurat from counts
assay_names <- tryCatch(as.character(assayNames(se_sub)), error = function(e) character(0))
assay_use <- if ("counts" %in% assay_names) "counts" else assay_names[1]

cd_sub <- colData(se_sub)
bad_cols <- sapply(colnames(cd_sub), function(x) is.list(cd_sub[[x]]))
if (any(bad_cols)) {
  colData(se_sub) <- cd_sub[, !bad_cols, drop = FALSE]
}

cnt_raw <- assay(se_sub, assay_use)
cnt <- methods::as(cnt_raw, "dgCMatrix")
rownames(cnt) <- rownames(se_sub)
colnames(cnt) <- colnames(se_sub)

so_sub <- suppressWarnings(CreateSeuratObject(counts = cnt, project = "dam_urgent"))
rm(cnt, cnt_raw)

DefaultAssay(so_sub) <- if ("RNA" %in% SeuratObject::Assays(so_sub)) "RNA" else SeuratObject::Assays(so_sub)[1]
so_sub <- NormalizeData(so_sub, verbose = FALSE)

so_sub$cell_type <- as.character(colData(se_sub)$cell_type)
so_sub$sample <- as.character(colData(se_sub)$sample)

# DAM signature and module score
sig_up <- read_upregulated_dam(sig_up_path, rownames(so_sub))
if (length(sig_up) == 0) {
  stop("No Upregulated_DAM genes matched the panel genes in this object")
}
message("Upregulated_DAM genes matched: ", length(sig_up))

ctrl_n <- max(1L, min(5L, floor((nrow(so_sub) - length(sig_up)) / 2)))
set.seed(1997)
so_sub <- AddModuleScore(
  so_sub,
  features = list(sig_up),
  name = "Upregulated_DAM",
  ctrl = ctrl_n,
  seed = 1997
)

score_df <- data.frame(
  cell_id = colnames(so_sub),
  cell_type = so_sub$cell_type,
  sample = so_sub$sample,
  dam_score = so_sub$Upregulated_DAM1,
  stringsAsFactors = FALSE
) %>%
  mutate(
    time_label = factor(TP_LABELS[sample], levels = c("1 wpi", "6 wpi"))
  )

# Section 1: Microglia (C1qa)
message("\n--- Section 1: Microglia (C1qa) ---")
df_micro <- score_df %>% filter(cell_type == MICRO_LABEL)
if (nrow(df_micro) == 0) {
  stop("No cells found for ", MICRO_LABEL, " at 1wpi/6wpi")
}

micro_res <- build_violin(
  df = df_micro,
  y_col = "dam_score",
  title_txt = "DAM module score in-niche microglia: 6wpi vs 1wpi",
  subtitle_prefix = "Microglia (C1qa)"
)
p_micro <- micro_res$plot
save_fig(p_micro, out_dam_dir, "fig_dam_score_inniche_6wpi_vs_1wpi_violin", 5.2, 5)

# Section 2: Mac (Ctss)
message("\n--- Section 2: Macrophages (Mac (Ctss)) ---")
df_mac <- score_df %>% filter(cell_type == MAC_LABEL)
if (nrow(df_mac) == 0) {
  stop("No cells found for ", MAC_LABEL, " at 1wpi/6wpi")
}

mac_res <- build_violin(
  df = df_mac,
  y_col = "dam_score",
  title_txt = "DAM module score macrophages: 6wpi vs 1wpi",
  subtitle_prefix = "Mac (Ctss)"
)
p_mac <- mac_res$plot
save_fig(p_mac, out_dam_dir, "fig_dam_score_mac_6wpi_vs_1wpi_violin", 5.2, 5)

# ==================================================================
# SECTION 3 — Combined DAM figure (microglia + macrophages)
# ==================================================================
message("\n=== SECTION 3: Combined DAM score figure ===\n")

p_combined <- p_micro + p_mac +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "DAM module score at chronic vs acute stage : microglia and macrophages",
    theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
  )

save_fig(p_combined, out_dam_dir, "fig_dam_score_combined_6wpi_vs_1wpi", 10.2, 5)

# ==================================================================
# SECTION 4 — Updated nDEG summary by cell type
# ==================================================================
message("\n=== SECTION 4: nDEG summary by cluster (updated) ===\n")

ndeg_path <- file.path(out_ndeg_dir, "ndeg_summary_table.csv")
if (!file.exists(ndeg_path)) {
  stop("nDEG summary file not found: ", ndeg_path)
}

ndeg <- read.csv(ndeg_path, stringsAsFactors = FALSE)
required_cols <- c("cell_type", "n_up", "n_down")
miss <- setdiff(required_cols, colnames(ndeg))
if (length(miss) > 0) {
  stop("Missing required columns in ndeg_summary_table.csv: ", paste(miss, collapse = ", "))
}

if (!"n_total" %in% colnames(ndeg)) {
  ndeg$n_total <- ndeg$n_up + ndeg$n_down
}

if ("skipped" %in% colnames(ndeg)) {
  ndeg <- ndeg %>% filter(!skipped)
}

contrast_col <- if ("contrast_label" %in% colnames(ndeg)) {
  "contrast_label"
} else if ("contrast" %in% colnames(ndeg)) {
  "contrast"
} else {
  NA_character_
}

if (is.na(contrast_col)) {
  ndeg$contrast_plot <- "All contrasts"
} else {
  ndeg$contrast_plot <- as.character(ndeg[[contrast_col]])
}

ct_order <- ndeg %>%
  group_by(cell_type) %>%
  summarise(total = sum(n_total, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(cell_type)

ndeg_long <- ndeg %>%
  mutate(cell_type = factor(cell_type, levels = rev(ct_order))) %>%
  select(cell_type, contrast_plot, n_up, n_down) %>%
  pivot_longer(cols = c("n_up", "n_down"),
               names_to = "direction", values_to = "n_deg") %>%
  mutate(
    direction = recode(direction, n_up = "up", n_down = "down"),
    n_signed = ifelse(direction == "down", -n_deg, n_deg)
  )

p_ndeg <- ggplot(ndeg_long, aes(x = cell_type, y = n_signed, fill = direction)) +
  geom_col(width = 0.72, colour = "grey30", linewidth = 0.2) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
  coord_flip() +
  scale_fill_manual(
    values = c("up" = "#B2182B", "down" = "#2166AC"),
    labels = c("up" = "Up", "down" = "Down"),
    name = NULL
  ) +
  scale_y_continuous(labels = function(x) abs(x), expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title = "nDEG summary by cell type (6wpi vs 1wpi and available contrasts)",
    subtitle = "Bars: up (red) and down (blue) significant DEGs",
    x = NULL,
    y = "n DEGs"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, colour = "grey40"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  )

if (length(unique(ndeg_long$contrast_plot)) > 1) {
  p_ndeg <- p_ndeg + facet_wrap(~contrast_plot, ncol = 1, scales = "free_y")
}

n_ct <- length(unique(ndeg_long$cell_type))
h_ndeg <- max(5.0, 0.32 * n_ct + 2.2)
save_fig(p_ndeg, out_ndeg_dir, "fig_ndeg_summary_by_celltype_updated", 8.5, h_ndeg)

# ==================================================================
# SECTION 5 — XY zoom niche highlighted (LCMV_1wpi)
# ==================================================================
message("\n=== SECTION 5: XY zoom niche highlighted (LCMV_1wpi) ===\n")

stopifnot(file.exists(obj_global_path))
se_global <- readRDS(obj_global_path)

cluster_col <- find_cl_col(se_global, 0.2, 0.9)
cd_global <- as.data.frame(SummarizedExperiment::colData(se_global))
spatial_mat <- as.data.frame(SpatialExperiment::spatialCoords(se_global))
if (!all(c("sdimx", "sdimy") %in% colnames(spatial_mat))) {
  stop("Spatial coordinates must contain columns sdimx and sdimy")
}

# Reconstruct annotations from CSV mapping
csv_annot <- "ncells_by_sample_lam02_res09_joint_long.csv"
if (!file.exists(csv_annot)) stop("Annotation CSV not found: ", csv_annot)
anno_data <- read.delim(csv_annot, sep = ";", stringsAsFactors = FALSE)
anno_map <- anno_data %>%
  filter(annotation != "" & !is.na(annotation)) %>%
  select(banksy_domain, annotation) %>%
  distinct()
anno_map$cluster_id <- as.numeric(gsub("Domain_", "", anno_map$banksy_domain))
anno_lookup <- setNames(trimws(anno_map$annotation), anno_map$cluster_id)

banksy_clusters <- as.numeric(cd_global[[cluster_col]])
cell_annotations <- anno_lookup[as.character(banksy_clusters)]
cell_annotations[is.na(cell_annotations)] <- "Non annote"

plot_df <- data.frame(
  x = spatial_mat$sdimx,
  y = spatial_mat$sdimy,
  sample = as.character(cd_global$sample),
  annotation = cell_annotations,
  stringsAsFactors = FALSE
)

df_1wpi <- plot_df %>% filter(sample == "LCMV_1wpi")
if (nrow(df_1wpi) == 0) {
  stop("No LCMV_1wpi cells found in global object")
}

immune_1wpi <- df_1wpi %>% filter(annotation == "Immune (Acod1)")
if (nrow(immune_1wpi) == 0) {
  stop("No Immune (Acod1) cells found in LCMV_1wpi")
}

cx <- mean(immune_1wpi$x, na.rm = TRUE)
cy <- mean(immune_1wpi$y, na.rm = TRUE)
half_win <- 800
xlim_zoom <- c(cx - half_win, cx + half_win)
ylim_zoom <- c(cy - half_win, cy + half_win)

df_zoom <- df_1wpi %>%
  filter(x >= xlim_zoom[1], x <= xlim_zoom[2],
         y >= ylim_zoom[1], y <= ylim_zoom[2])

# Highlight classes requested by user
is_neuron <- grepl("Neurons", df_zoom$annotation)
is_astro <- grepl("Astrocytes", df_zoom$annotation)

df_zoom$display <- case_when(
  df_zoom$annotation == "Immune (Acod1)" ~ "Immune (Acod1)",
  df_zoom$annotation == "IFN responsive (Ifit1)" ~ "IFN responsive (Ifit1)",
  df_zoom$annotation == "Microglia (P2ry12)" ~ "Microglia (P2ry12)",
  is_neuron ~ "Neurons",
  is_astro ~ "Astrocytes",
  TRUE ~ "Other"
)

display_levels <- c(
  "Other",
  "Neurons",
  "Astrocytes",
  "Microglia (P2ry12)",
  "IFN responsive (Ifit1)",
  "Immune (Acod1)"
)
df_zoom$display <- factor(df_zoom$display, levels = display_levels)

display_colors <- c(
  "Immune (Acod1)" = "#D62728",
  "IFN responsive (Ifit1)" = "#F28E2B",
  "Microglia (P2ry12)" = "#1F77B4",
  "Neurons" = "#D9D9D9",
  "Astrocytes" = "#2CA02C",
  "Other" = "#EFEFEF"
)

# Plot other cells first, then highlighted populations
bg_df <- df_zoom %>% filter(display == "Other")
fg_df <- df_zoom %>% filter(display != "Other")

p_zoom <- ggplot() +
  geom_point(data = bg_df, aes(x = x, y = y, color = display),
             size = 0.3, alpha = 0.2, stroke = 0) +
  geom_point(data = fg_df, aes(x = x, y = y, color = display),
             size = 1.5, alpha = 0.95, stroke = 0) +
  scale_color_manual(values = display_colors, drop = FALSE, name = "Annotation") +
  coord_fixed(xlim = xlim_zoom, ylim = ylim_zoom, expand = FALSE) +
  labs(
    title = "LCMV 1wpi niche zoom with highlighted annotations",
    subtitle = sprintf("Fenetre +/- %d um autour du centroide des cellules Immune (Acod1)", half_win),
    x = "X (um)",
    y = "Y (um)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.4)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

save_fig(p_zoom, out_spatial_dir, "xy_zoom_niche_highlighted_lcmv_1wpi", 7.2, 7)

message("\n=== Script 42 terminé avec succès ===\n")
