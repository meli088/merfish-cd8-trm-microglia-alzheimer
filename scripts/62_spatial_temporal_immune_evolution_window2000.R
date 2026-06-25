#!/usr/bin/env Rscript

set.seed(1997)

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Cairo)
})

base_path <- normalizePath(".")
setwd(base_path)

source("scripts/00_palette.R")

OBJ_FILE <- file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds")
IMMUNE_SUB_FILE <- file.path("objects", "08_immune_annotated_lam02_res03.rds")
ANNOT_CSV <- file.path("ncells_by_sample_lam02_res09_joint_long.csv")
OUT_DIR <- file.path("outputs", "banksy", "immune_acod1", "analysis", "figures")
OUT_BASE_MAP <- file.path(OUT_DIR, "fig_spatial_temporal_immune_evolution_window2000")
OUT_BASE_BAR <- file.path(OUT_DIR, "fig_spatial_temporal_immune_counts_window2000")
OUT_COUNTS <- file.path(OUT_DIR, "spatial_temporal_immune_window2000_counts.csv")

TIME_ORDER <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")
WINDOW_HALF_UM <- 1500
SCALE_BAR_UM <- 500

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
}

find_cl_col <- function(cd, lam = 0.2, res = 0.9) {
  cn <- colnames(cd)
  lam_pat <- paste0("lam", gsub("\\.", "\\\\.", as.character(lam)))
  res_pat <- paste0("_res", gsub("\\.", "\\\\.", as.character(res)), "$")
  cols <- cn[grepl(lam_pat, cn) & grepl(res_pat, cn)]
  if (length(cols) == 0) {
    stop("Could not find clustering column for lam=", lam, " and res=", res)
  }
  cols[1]
}

if (!file.exists(OBJ_FILE)) stop("Missing object: ", OBJ_FILE)
if (!file.exists(IMMUNE_SUB_FILE)) stop("Missing immune subcluster object: ", IMMUNE_SUB_FILE)
if (!file.exists(ANNOT_CSV)) stop("Missing annotation CSV: ", ANNOT_CSV)

message("Loading: ", OBJ_FILE)
se <- readRDS(OBJ_FILE)
se_immune <- readRDS(IMMUNE_SUB_FILE)
cd <- as.data.frame(SummarizedExperiment::colData(se))
cd_immune <- as.data.frame(SummarizedExperiment::colData(se_immune))
xy <- as.data.frame(SpatialExperiment::spatialCoords(se))

if (!all(c("sample") %in% colnames(cd))) {
  stop("colData must contain 'sample'")
}
if (!all(c("cell_type") %in% colnames(cd_immune))) {
  stop("immune subcluster object colData must contain 'cell_type'")
}
if (!all(c("sdimx", "sdimy") %in% colnames(xy))) {
  stop("spatialCoords must contain sdimx/sdimy")
}

resolve_immune_color <- function(label) {
  if (label %in% names(IMMUNE_PALETTE)) {
    return(unname(IMMUNE_PALETTE[[label]]))
  }
  "#999999"  # Fallback gray
}

cl_col <- find_cl_col(cd, lam = 0.2, res = 0.9)
cluster_ids <- as.numeric(cd[[cl_col]])

annot_long <- read.delim(ANNOT_CSV, sep = ";", stringsAsFactors = FALSE)
annot_lookup <- annot_long %>%
  filter(!is.na(annotation), annotation != "") %>%
  transmute(
    cluster_id = as.numeric(gsub("Domain_", "", as.character(banksy_domain))),
    annotation = trimws(as.character(annotation))
  ) %>%
  distinct(cluster_id, .keep_all = TRUE)

annotation <- annot_lookup$annotation[match(cluster_ids, annot_lookup$cluster_id)]
annotation[is.na(annotation)] <- "Non annote"

immune_map <- setNames(as.character(cd_immune$cell_type), colnames(se_immune))

df <- data.frame(
  cell_id = colnames(se),
  sample = as.character(cd$sample),
  annotation = annotation,
  x = as.numeric(xy$sdimx),
  y = as.numeric(xy$sdimy),
  stringsAsFactors = FALSE
) %>%
  filter(sample %in% TIME_ORDER, !is.na(x), !is.na(y))

if (nrow(df) == 0) stop("No cells for requested timepoints")

# Shared centered coordinates for same absolute spatial window.
df <- df %>%
  group_by(sample) %>%
  mutate(
    x_center = x - median(x, na.rm = TRUE),
    y_center = y - median(y, na.rm = TRUE)
  ) %>%
  ungroup()

df_win <- df %>%
  filter(
    x_center >= -WINDOW_HALF_UM, x_center <= WINDOW_HALF_UM,
    y_center >= -WINDOW_HALF_UM, y_center <= WINDOW_HALF_UM
  ) %>%
  mutate(
    is_immune = annotation == "Immune (Acod1)",
    immune_subcluster = immune_map[cell_id],
    immune_subcluster = ifelse(is_immune, immune_subcluster, NA_character_)
  )

immune_levels <- sort(unique(df_win$immune_subcluster[!is.na(df_win$immune_subcluster)]))
immune_levels <- order_annotations(immune_levels, extended = TRUE)
immune_levels <- immune_levels[immune_levels %in% unique(df_win$immune_subcluster[!is.na(df_win$immune_subcluster)])]
immune_palette <- setNames(vapply(immune_levels, resolve_immune_color, character(1)), immune_levels)

counts <- df_win %>%
  group_by(sample) %>%
  summarise(
    n_total = n(),
    n_immune = sum(is_immune, na.rm = TRUE),
    pct_immune = round(100 * n_immune / pmax(1, n_total), 2),
    .groups = "drop"
)

counts$sample <- factor(counts$sample, levels = TIME_ORDER)
counts <- counts %>% arrange(sample)
write.csv(counts, OUT_COUNTS, row.names = FALSE)
message("Saved: ", OUT_COUNTS)

df_win$sample <- factor(df_win$sample, levels = TIME_ORDER)

sb_x1 <- WINDOW_HALF_UM - 180
sb_x0 <- sb_x1 - SCALE_BAR_UM
sb_y <- -WINDOW_HALF_UM + 180

sb_df <- data.frame(
  sample = factor(TIME_ORDER, levels = TIME_ORDER),
  x0 = sb_x0,
  x1 = sb_x1,
  y = sb_y,
  label_x = (sb_x0 + sb_x1) / 2,
  label_y = sb_y + 120,
  label = paste0(SCALE_BAR_UM, " um"),
  stringsAsFactors = FALSE
)

labs_df <- counts %>%
  transmute(
    sample,
    label = paste0("n=", format(n_total, big.mark = ",", scientific = FALSE)),
    x = -WINDOW_HALF_UM + 120,
    y = WINDOW_HALF_UM - 120
  )

p_map <- ggplot() +
  geom_point(
    data = df_win,
    aes(x = x_center, y = y_center),
    color = "#CECECE",
    size = 0.25,
    alpha = 0.72,
    stroke = 0
  ) +
  geom_point(
    data = df_win %>% filter(is_immune),
    aes(x = x_center, y = y_center, color = immune_subcluster),
    size = 0.48,
    alpha = 0.95,
    stroke = 0
  ) +
  scale_color_manual(values = immune_palette, breaks = immune_levels, drop = FALSE, name = "Immune subclusters") +
  geom_segment(
    data = sb_df,
    aes(x = x0, xend = x1, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.9,
    color = "black"
  ) +
  geom_text(
    data = sb_df,
    aes(x = label_x, y = label_y, label = label),
    inherit.aes = FALSE,
    size = 2.8,
    color = "black"
  ) +
  geom_text(
    data = labs_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3,
    color = "black"
  ) +
  facet_wrap(~ sample, nrow = 1) +
  coord_fixed(
    xlim = c(-WINDOW_HALF_UM, WINDOW_HALF_UM),
    ylim = c(-WINDOW_HALF_UM, WINDOW_HALF_UM),
    expand = FALSE
  ) +
  labs(
    title = "Temporal evolution of immune spatial distribution",
    subtitle = "Shared window: +/-2000 um | all cells in gray, Immune (Acod1) subclusters colored",
    x = "X (um, centered per sample)",
    y = "Y (um, centered per sample)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey35"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.65, "cm"),
    panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.3)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3.2, alpha = 1)))

p_bar <- ggplot(counts, aes(x = sample, y = n_total)) +
  geom_col(fill = "#BDBDBD", width = 0.7) +
  geom_text(aes(label = format(n_total, big.mark = ",", scientific = FALSE)), vjust = -0.35, size = 3.1) +
  labs(
    title = "Cells inside +/-2000 um window",
    x = NULL,
    y = "Number of cells"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
    axis.text.x = element_text(face = "bold")
  )

cairo_pdf(paste0(OUT_BASE_MAP, ".pdf"), width = 14.2, height = 5.8)
print(p_map)
dev.off()

CairoJPEG(
  paste0(OUT_BASE_MAP, ".jpg"),
  width = round(14.2 * 180),
  height = round(5.8 * 180),
  res = 180,
  quality = 95
)
print(p_map)
dev.off()

cairo_pdf(paste0(OUT_BASE_BAR, ".pdf"), width = 7.2, height = 4.8)
print(p_bar)
dev.off()

CairoJPEG(
  paste0(OUT_BASE_BAR, ".jpg"),
  width = round(7.2 * 180),
  height = round(4.8 * 180),
  res = 180,
  quality = 95
)
print(p_bar)
dev.off()

message("Saved: ", paste0(OUT_BASE_MAP, ".pdf"))
message("Saved: ", paste0(OUT_BASE_MAP, ".jpg"))
message("Saved: ", paste0(OUT_BASE_BAR, ".pdf"))
message("Saved: ", paste0(OUT_BASE_BAR, ".jpg"))
