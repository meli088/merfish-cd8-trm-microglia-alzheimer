#!/usr/bin/env Rscript

set.seed(1997)

suppressPackageStartupMessages({
  library(Seurat)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(Cairo)
})

base_path <- normalizePath(".")
setwd(base_path)

source("scripts/00_palette.R")

OUT_DIR <- file.path("outputs", "banksy", "rfp_plots_readable")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SAMPLE_6WPI <- "LCMV_6wpi"
SAMPLE_MOCK <- "mock_6wpi"

NEURON_SUBTYPES <- c(
  "Neurons (Adora2a)", "Neurons (Arhgap36)", "Neurons (Crhbp)",
  "Neurons (Dkk3)", "Neurons (Htr2c)", "Neurons (Nefm)",
  "Neurons (Neurod1)", "Neurons (Rxfp1)", "Neurons (Slc17a8)"
)

rename_annotation <- function(x) {
  x <- as.character(x)
  x[x %in% NEURON_SUBTYPES] <- "Neurons"
  x[x == "Prolif neural/glial (Ccdc153)"] <- "Ependymal (Ccdc153)"
  x
}

find_cl_col <- function(se, lam, res) {
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  cn <- colnames(cd)
  lam_pat <- paste0("lam", gsub("\\.", "\\\\.", as.character(lam)))
  res_pat <- paste0("_res", gsub("\\.", "\\\\.", as.character(res)), "$")
  cols <- cn[grepl(lam_pat, cn) & grepl(res_pat, cn)]
  if (length(cols) == 0) {
    stop("Clustering column not found for lambda=", lam, " and resolution=", res)
  }
  cols[1]
}

save_cairo_pair <- function(plot_obj, out_stem, width, height, jpg_res = 170) {
  cairo_pdf(paste0(out_stem, ".pdf"), width = width, height = height)
  print(plot_obj)
  dev.off()

  CairoJPEG(
    paste0(out_stem, ".jpg"),
    width = round(width * jpg_res),
    height = round(height * jpg_res),
    res = jpg_res,
    quality = 95
  )
  print(plot_obj)
  dev.off()

  message("Saved: ", paste0(out_stem, ".pdf"))
  message("Saved: ", paste0(out_stem, ".jpg"))
}

# =============================================================
# SECTION 1 — RFP enrichment (6wpi vs mock) for all global clusters
# =============================================================
message("\n=== Section 1: RFP enrichment 6wpi vs mock ===")

obj_rfp <- readRDS(file.path("objects", "03_rfp_analysis.rds"))
md <- as.data.frame(obj_rfp@meta.data)

required_cols <- c("sample", "Anti.RFP_raw", "Anti.RFP_high_pass")
missing_cols <- setdiff(required_cols, colnames(md))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns in 03_rfp_analysis.rds: ",
       paste(missing_cols, collapse = ", "))
}

if ("CellType" %in% colnames(md)) {
  md$annotation <- as.character(md$CellType)
} else if ("cell_type" %in% colnames(md)) {
  md$annotation <- as.character(md$cell_type)
} else {
  stop("No CellType/cell_type column found in 03_rfp_analysis.rds")
}

md <- md %>%
  mutate(
    sample = as.character(sample),
    annotation = rename_annotation(annotation),
    Anti.RFP_raw = as.numeric(Anti.RFP_raw),
    Anti.RFP_high_pass = as.numeric(Anti.RFP_high_pass)
  ) %>%
  filter(sample %in% c(SAMPLE_6WPI, SAMPLE_MOCK),
         !is.na(annotation), annotation != "", annotation != "Non annote")

md_long <- md %>%
  select(annotation, sample, Anti.RFP_raw, Anti.RFP_high_pass) %>%
  pivot_longer(
    cols = c(Anti.RFP_raw, Anti.RFP_high_pass),
    names_to = "measure",
    values_to = "rfp_value"
  ) %>%
  mutate(
    measure = recode(measure,
                     "Anti.RFP_raw" = "raw",
                     "Anti.RFP_high_pass" = "highpass")
  ) %>%
  filter(!is.na(rfp_value), is.finite(rfp_value))

rfp_summary <- md_long %>%
  group_by(measure, annotation, sample) %>%
  summarise(
    median_rfp = median(rfp_value, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

rfp_wide <- rfp_summary %>%
  select(measure, annotation, sample, median_rfp, n_cells) %>%
  pivot_wider(
    names_from = sample,
    values_from = c(median_rfp, n_cells),
    values_fill = list(median_rfp = NA_real_, n_cells = 0)
  )

# Robust pseudocount to avoid division by zero when medians are 0.
eps <- 1e-9
rfp_wide <- rfp_wide %>%
  mutate(
    log2FC = log2((median_rfp_LCMV_6wpi + eps) / (median_rfp_mock_6wpi + eps))
  ) %>%
  arrange(measure, desc(log2FC))

# Wilcoxon test per annotation and measure: LCMV_6wpi vs mock_6wpi
wilcox_df <- md_long %>%
  group_by(measure, annotation) %>%
  summarise(
    p_wilcox = tryCatch(
      wilcox.test(
        rfp_value[sample == SAMPLE_6WPI],
        rfp_value[sample == SAMPLE_MOCK],
        alternative = "two.sided",
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  group_by(measure) %>%
  mutate(
    p_adj_bh = p.adjust(p_wilcox, method = "BH"),
    sig_label = case_when(
      is.na(p_adj_bh) ~ "ns",
      p_adj_bh < 1e-4 ~ "****",
      p_adj_bh < 1e-3 ~ "***",
      p_adj_bh < 1e-2 ~ "**",
      p_adj_bh < 5e-2 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  ungroup()

rfp_wide <- rfp_wide %>%
  left_join(wilcox_df, by = c("measure", "annotation"))

if (nrow(rfp_wide) == 0) {
  stop("No rows available for 6wpi vs mock comparison in Section 1")
}

palette_local <- c(
  GLOBAL_PALETTE,
  "Ependymal (Ccdc153)" = GLOBAL_PALETTE[["Prolif neural/glial (Ccdc153)"]]
)

for (m in c("raw", "highpass")) {
  rfp_m <- rfp_wide %>%
    filter(measure == m) %>%
    arrange(desc(log2FC))

  if (nrow(rfp_m) == 0) next

  rfp_m$annotation <- factor(rfp_m$annotation, levels = rev(rfp_m$annotation))
  rfp_m$bar_color <- unname(palette_local[as.character(rfp_m$annotation)])
  rfp_m$bar_color[is.na(rfp_m$bar_color)] <- "grey75"

  xmax_abs <- max(abs(rfp_m$log2FC), na.rm = TRUE)
  x_text <- ifelse(rfp_m$log2FC >= 0,
                   rfp_m$log2FC + 0.04 * xmax_abs,
                   rfp_m$log2FC - 0.04 * xmax_abs)
  rfp_m$hjust_sig <- ifelse(rfp_m$log2FC >= 0, 0, 1)
  rfp_m$x_text_sig <- x_text

  measure_label <- ifelse(m == "raw", "Anti.RFP_raw", "Anti.RFP_high_pass")

  p_fc <- ggplot(rfp_m, aes(x = log2FC, y = annotation, fill = bar_color)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.35, color = "grey40") +
    scale_fill_identity() +
    labs(
      title = "RFP signal enrichment at 6wpi vs mock",
      subtitle = paste0("log2FC median RFP intensity (LCMV 6wpi / mock 6wpi) | ", measure_label),
      x = paste0("log2FC median ", measure_label),
      y = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey35"),
      axis.text.y = element_text(size = 9)
    )

  p_fc_wilcox <- p_fc +
    geom_text(
      data = rfp_m,
      aes(x = x_text_sig, y = annotation, label = sig_label, hjust = hjust_sig),
      inherit.aes = FALSE,
      size = 3,
      color = "black"
    ) +
    labs(
      subtitle = paste0(
        "log2FC median RFP intensity (LCMV 6wpi / mock 6wpi) | ",
        measure_label,
        " + Wilcoxon (BH): **** <1e-4, *** <1e-3, ** <1e-2, * <0.05, ns"
      )
    )

  save_cairo_pair(
    p_fc,
    file.path(OUT_DIR, paste0("fig_rfp_fc_6wpi_vs_mock_allclusters_", m, "_no_wilcox")),
    width = 10,
    height = max(6.5, 0.35 * nrow(rfp_m) + 1.8)
  )

  save_cairo_pair(
    p_fc_wilcox,
    file.path(OUT_DIR, paste0("fig_rfp_fc_6wpi_vs_mock_allclusters_", m, "_wilcox")),
    width = 10,
    height = max(6.5, 0.35 * nrow(rfp_m) + 1.8)
  )
}

write.csv(
  rfp_wide %>% select(measure, annotation, log2FC,
                      median_rfp_LCMV_6wpi, median_rfp_mock_6wpi,
                      n_cells_LCMV_6wpi, n_cells_mock_6wpi,
                      p_wilcox, p_adj_bh, sig_label),
  file.path(OUT_DIR, "rfp_fc_6wpi_vs_mock_wilcox_table.csv"),
  row.names = FALSE
)
message("Saved: ", file.path(OUT_DIR, "rfp_fc_6wpi_vs_mock_wilcox_table.csv"))

# =============================================================
# SECTION 2 — XY map of Anti.RFP_raw in 6wpi with niche zoom
# =============================================================
message("\n=== Section 2: XY RFP map 6wpi + niche zoom ===")

se_joint <- readRDS(file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds"))
cd <- as.data.frame(SummarizedExperiment::colData(se_joint))

if (!"sample" %in% colnames(cd)) {
  stop("Column 'sample' not found in colData of 04_banksy_joint_lam08_after_bloc3.rds")
}

# Build global annotation from BANKSY domains + CSV mapping.
cluster_col <- find_cl_col(se_joint, lam = 0.2, res = 0.9)
annot_long <- read.delim("ncells_by_sample_lam02_res09_joint_long.csv", sep = ";", stringsAsFactors = FALSE)
anno_map <- annot_long %>%
  filter(!is.na(annotation), annotation != "") %>%
  transmute(
    banksy_domain = as.character(banksy_domain),
    annotation = rename_annotation(trimws(as.character(annotation)))
  ) %>%
  distinct()
anno_map$cluster_id <- as.numeric(gsub("Domain_", "", anno_map$banksy_domain))
anno_lookup <- setNames(anno_map$annotation, anno_map$cluster_id)

banksy_clusters <- as.numeric(cd[[cluster_col]])
annotation <- anno_lookup[as.character(banksy_clusters)]
annotation[is.na(annotation)] <- "Non annote"
annotation <- rename_annotation(annotation)

xy <- as.data.frame(SpatialExperiment::spatialCoords(se_joint))
if (all(c("sdimx", "sdimy") %in% colnames(xy))) {
  x <- as.numeric(xy$sdimx)
  y <- as.numeric(xy$sdimy)
} else {
  x <- as.numeric(xy[[1]])
  y <- as.numeric(xy[[2]])
}

# Retrieve Anti.RFP_raw from se_joint if available, fallback to 03_rfp_analysis.rds metadata.
if ("Anti.RFP_raw" %in% colnames(cd)) {
  rfp_vals <- as.numeric(cd$Anti.RFP_raw)
} else {
  md_rfp <- as.data.frame(obj_rfp@meta.data)
  if (!"Anti.RFP_raw" %in% colnames(md_rfp)) {
    stop("Anti.RFP_raw not found in both se_joint colData and 03_rfp_analysis.rds meta.data")
  }
  ids_joint <- colnames(se_joint)
  ids_rfp <- rownames(md_rfp)
  rfp_vals <- as.numeric(md_rfp[ids_joint, "Anti.RFP_raw"])
}

plot_df <- data.frame(
  cell_id = colnames(se_joint),
  sample = as.character(cd$sample),
  annotation = annotation,
  sdimx = x,
  sdimy = y,
  Anti.RFP_raw = rfp_vals,
  stringsAsFactors = FALSE
) %>%
  filter(sample == SAMPLE_6WPI, !is.na(Anti.RFP_raw), is.finite(Anti.RFP_raw))

if (nrow(plot_df) == 0) {
  stop("No LCMV_6wpi cells with Anti.RFP_raw available for Section 2")
}

imm_df <- plot_df %>% filter(annotation == "Immune (Acod1)")
if (nrow(imm_df) == 0) {
  message("No Immune (Acod1) cells found in LCMV_6wpi; using full tissue centroid for zoom")
  cx <- median(plot_df$sdimx, na.rm = TRUE)
  cy <- median(plot_df$sdimy, na.rm = TRUE)
} else {
  cx <- median(imm_df$sdimx, na.rm = TRUE)
  cy <- median(imm_df$sdimy, na.rm = TRUE)
}

zoom_half <- 800
xmin <- cx - zoom_half
xmax <- cx + zoom_half
ymin <- cy - zoom_half
ymax <- cy + zoom_half

zoom_df <- plot_df %>%
  filter(sdimx >= xmin, sdimx <= xmax, sdimy >= ymin, sdimy <= ymax)

# Shared color scale (white -> red) across global and zoom.
rfp_max <- quantile(plot_df$Anti.RFP_raw, 0.995, na.rm = TRUE)
rfp_max <- max(rfp_max, max(plot_df$Anti.RFP_raw, na.rm = TRUE) * 0.01)

p_global <- ggplot(plot_df, aes(x = sdimx, y = sdimy, color = Anti.RFP_raw)) +
  geom_point(size = 0.23, alpha = 0.85, stroke = 0) +
  scale_color_gradient(
    low = "#FFFFFF",
    high = "#D7191C",
    limits = c(0, rfp_max),
    oob = scales::squish,
    name = "Anti.RFP_raw"
  ) +
  labs(
    title = "LCMV 6wpi — Anti.RFP_raw (global)",
    x = "X (µm)",
    y = "Y (µm)"
  ) +
  coord_fixed() +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10.5),
    legend.position = "right"
  )

bar_len <- 200
bar_margin <- 60
bar_x0 <- xmax - bar_margin - bar_len
bar_x1 <- xmax - bar_margin
bar_y <- ymin + bar_margin

p_zoom <- ggplot(zoom_df, aes(x = sdimx, y = sdimy, color = Anti.RFP_raw)) +
  geom_point(size = 0.45, alpha = 0.9, stroke = 0) +
  scale_color_gradient(
    low = "#FFFFFF",
    high = "#D7191C",
    limits = c(0, rfp_max),
    oob = scales::squish,
    name = "Anti.RFP_raw"
  ) +
  coord_fixed(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  annotate("segment", x = bar_x0, xend = bar_x1, y = bar_y, yend = bar_y,
           linewidth = 1.2, color = "black") +
  annotate("text", x = (bar_x0 + bar_x1) / 2, y = bar_y + 45,
           label = "200 µm", size = 3, color = "black") +
  labs(
    title = "Niche-centered zoom (±800 µm)",
    x = "X (µm)",
    y = "Y (µm)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10.5),
    legend.position = "none"
  )

p_xy <- p_global + p_zoom +
  plot_layout(widths = c(1.25, 1)) +
  plot_annotation(
    title = "Anti.RFP_raw spatial intensity in LCMV 6wpi",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12))
  )

save_cairo_pair(
  p_xy,
  file.path(OUT_DIR, "fig_rfp_xy_6wpi_zoom_niche"),
  width = 14,
  height = 6.8
)

message("\nDone: script 59 completed successfully.")
