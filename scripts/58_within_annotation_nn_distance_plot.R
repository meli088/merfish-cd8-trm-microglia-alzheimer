#!/usr/bin/env Rscript

set.seed(1997)

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(FNN)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(Cairo)
})

base_path <- normalizePath(".")
setwd(base_path)

source("scripts/00_palette.R")

summary_csv <- file.path("outputs", "banksy", "within_annotation_nn_summary.csv")
out_dir <- file.path("outputs", "banksy", "spatial_annotations")
out_base <- file.path(out_dir, "fig_within_annotation_nn_distance_by_condition")
sample_order <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi", "mock_6wpi")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

find_cl_col <- function(se, lam, res) {
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  cn <- colnames(cd)
  lam_pat <- paste0("lam", gsub("\\.", "\\\\.", as.character(lam)))
  res_pat <- paste0("_res", gsub("\\.", "\\\\.", as.character(res)), "$")
  cols <- cn[grepl(lam_pat, cn) & grepl(res_pat, cn)]
  if (length(cols) == 0) {
    stop("Clustering column not found for lambda=", lam, " resolution=", res)
  }
  cols[1]
}

weighted_median <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  if (length(x) == 0) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= 0.5)[1]]
}

compute_global_summary <- function() {
  message("Building global within-annotation NN summary from object...")

  se <- readRDS(file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds"))
  cluster_col <- find_cl_col(se, lam = 0.2, res = 0.9)

  cd <- as.data.frame(SummarizedExperiment::colData(se))
  samples <- as.character(cd$sample)
  banksy_clusters <- as.numeric(cd[[cluster_col]])

  annot_long <- read.delim("ncells_by_sample_lam02_res09_joint_long.csv", sep = ";", stringsAsFactors = FALSE)
  anno_map <- annot_long[!is.na(annot_long$annotation) & annot_long$annotation != "",
                         c("banksy_domain", "annotation")]
  anno_map$banksy_domain <- as.character(anno_map$banksy_domain)
  anno_map$annotation <- trimws(as.character(anno_map$annotation))
  anno_map$annotation[anno_map$annotation == "Prolif neural/glial (Ccdc153)"] <- "Ependymal (Ccdc153)"
  anno_map <- unique(anno_map)

  anno_map$cluster_id <- as.numeric(gsub("Domain_", "", anno_map$banksy_domain))
  anno_lookup <- setNames(anno_map$annotation, anno_map$cluster_id)
  annotation <- anno_lookup[as.character(banksy_clusters)]
  annotation[is.na(annotation)] <- "Non annote"

  xy <- as.data.frame(SpatialExperiment::spatialCoords(se))
  if (all(c("sdimx", "sdimy") %in% colnames(xy))) {
    x <- xy$sdimx
    y <- xy$sdimy
  } else {
    x <- xy[[1]]
    y <- xy[[2]]
  }

  df <- data.frame(
    sample = samples,
    annotation = annotation,
    x = x,
    y = y,
    stringsAsFactors = FALSE
  ) %>%
    filter(annotation != "Non annote")

  split_key <- paste(df$sample, df$annotation, sep = "||")
  idx_groups <- split(seq_len(nrow(df)), split_key)

  summary_rows <- lapply(idx_groups, function(idx) {
    sub <- df[idx, , drop = FALSE]
    n_cells <- nrow(sub)

    if (n_cells < 2) {
      nn <- rep(NA_real_, n_cells)
    } else {
      coords <- as.matrix(sub[, c("x", "y")])
      knn <- FNN::get.knnx(coords, coords, k = 2)
      nn <- knn$nn.dist[, 2]
    }

    data.frame(
      annotation = sub$annotation[1],
      sample = sub$sample[1],
      n_cells = n_cells,
      n_valid = sum(is.finite(nn)),
      median_nn_um = median(nn, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  summary_df <- bind_rows(summary_rows) %>%
    arrange(annotation, factor(sample, levels = sample_order))

  write.csv(summary_df, summary_csv, row.names = FALSE)
  message("Saved: ", summary_csv)
}

if (!file.exists(summary_csv)) {
  compute_global_summary()
}

message("Reading: ", summary_csv)
nn_summary <- read.csv(summary_csv, stringsAsFactors = FALSE)

if (!"annotation" %in% colnames(nn_summary)) {
  if ("cell_type" %in% colnames(nn_summary)) {
    nn_summary$annotation <- nn_summary$cell_type
  } else {
    stop("No annotation column found in summary CSV.")
  }
}

if (!"median_nn_um" %in% colnames(nn_summary)) {
  stop("Column median_nn_um is required in summary CSV.")
}

if ("sample" %in% colnames(nn_summary)) {
  nn_summary$sample <- factor(nn_summary$sample,
    levels = sample_order
  )
}

nn_summary$annotation <- dplyr::recode(
  as.character(nn_summary$annotation),
  "Prolif neural/glial (Ccdc153)" = "Ependymal (Ccdc153)"
)

if (!"n_valid" %in% colnames(nn_summary)) {
  nn_summary$n_valid <- 1
}

neuron_subtypes <- c(
  "Neurons (Adora2a)", "Neurons (Arhgap36)", "Neurons (Crhbp)",
  "Neurons (Dkk3)", "Neurons (Htr2c)", "Neurons (Nefm)",
  "Neurons (Neurod1)", "Neurons (Rxfp1)", "Neurons (Slc17a8)"
)

plot_df <- nn_summary %>%
  filter(is.finite(median_nn_um)) %>%
  mutate(
    annotation = ifelse(annotation %in% neuron_subtypes, "Neurons", annotation),
    sample = factor(as.character(sample), levels = sample_order)
  ) %>%
  group_by(sample, annotation) %>%
  summarise(
    median_nn_um = weighted_median(median_nn_um, n_valid),
    n_valid = sum(n_valid, na.rm = TRUE),
    .groups = "drop"
  )

if (nrow(plot_df) == 0) {
  stop("No valid data points available to plot.")
}

order_ref <- plot_df %>%
  filter(sample == "LCMV_1wpi") %>%
  arrange(median_nn_um) %>%
  pull(annotation)

if (length(order_ref) == 0) {
  order_ref <- character(0)
}

order_other <- plot_df %>%
  filter(!annotation %in% order_ref) %>%
  group_by(annotation) %>%
  summarise(m = weighted_median(median_nn_um, n_valid), .groups = "drop") %>%
  arrange(m) %>%
  pull(annotation)

annotation_levels <- unique(c(order_ref, order_other))
plot_df$annotation <- factor(plot_df$annotation, levels = annotation_levels)

plot_df <- plot_df %>% arrange(sample, annotation)

palette_local <- c(
  GLOBAL_PALETTE,
  "Ependymal (Ccdc153)" = GLOBAL_PALETTE[["Prolif neural/glial (Ccdc153)"]]
)
plot_df$bar_color <- ifelse(
  as.character(plot_df$annotation) == "Immune (Acod1)",
  palette_local[["Immune (Acod1)"]],
  ifelse(
    as.character(plot_df$annotation) == "IFN responsive (Ifit1)",
    palette_local[["IFN responsive (Ifit1)"]],
    "grey75"
  )
)

p <- ggplot(plot_df, aes(x = median_nn_um, y = annotation, fill = bar_color)) +
  geom_col(width = 0.75) +
  scale_fill_identity() +
  facet_wrap(~sample, ncol = 1, scales = "fixed") +
  labs(
    title = "Within-cluster nearest neighbor distance by cell type",
    x = "median_nn_um",
    y = "annotation"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

pdf_file <- paste0(out_base, ".pdf")
jpg_file <- paste0(out_base, ".jpg")

cairo_pdf(pdf_file, width = 10, height = 14)
print(p)
dev.off()

CairoJPEG(jpg_file, width = 1500, height = 2100, res = 150, quality = 95)
print(p)
dev.off()

message("Saved: ", pdf_file)
message("Saved: ", jpg_file)
