#!/usr/bin/env Rscript

set.seed(1997)

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(dplyr)
  library(ggplot2)
  library(Cairo)
})

base_path <- normalizePath(".")
setwd(base_path)

source("scripts/00_palette.R")

OBJ_FILE <- file.path("objects", "08_immune_annotated_lam02_res03.rds")
UP_SIG_FILE <- file.path("outputs", "banksy", "dam_signature_curation", "Upregulated_DAM.csv")
DOWN_SIG_FILE <- file.path("outputs", "banksy", "dam_signature_curation", "Downregulated_DAM.csv")
OUT_DIR <- file.path("outputs", "banksy", "immune_acod1", "analysis", "figures")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TIME_REF <- "LCMV_1wpi"
TIME_CASE <- "LCMV_6wpi"
TIME_LABELS <- c(LCMV_1wpi = "LCMV 1 wpi", LCMV_6wpi = "LCMV 6 wpi")
TIME_COLORS <- c("LCMV 1 wpi" = "#56B4E9", "LCMV 6 wpi" = "#F28E2B")

POPULATIONS <- c(
  "Activated microglia (C1qa)",
  "Antigen-presenting myeloid cells (Cd74)"
)

read_sig <- function(path) {
  if (!file.exists(path)) stop("Missing signature file: ", path)
  x <- read.csv(path, stringsAsFactors = FALSE)
  if (!"gene" %in% colnames(x)) stop("Column 'gene' not found in: ", path)
  unique(as.character(x$gene))
}

match_signature_genes <- function(sig_genes, panel_genes) {
  panel_title <- tools::toTitleCase(tolower(panel_genes))
  map_panel <- setNames(panel_genes, panel_title)
  sig_title <- tools::toTitleCase(tolower(sig_genes))
  hits <- sig_title[sig_title %in% panel_title]
  out <- unique(as.character(map_panel[hits]))
  out[!is.na(out)]
}

pval_stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "n.s.",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "n.s."
  )
}

save_fig <- function(p, out_base, width = 6.4, height = 4.9) {
  cairo_pdf(paste0(out_base, ".pdf"), width = width, height = height)
  print(p)
  dev.off()

  CairoJPEG(
    paste0(out_base, ".jpg"),
    width = width * 150,
    height = height * 150,
    res = 150,
    quality = 95
  )
  print(p)
  dev.off()
}

compute_plot_df <- function(se, pop_label, sig_genes, sig_label) {
  cd <- as.data.frame(colData(se))
  stopifnot("cell_type missing in object" = "cell_type" %in% colnames(cd))
  stopifnot("sample missing in object" = "sample" %in% colnames(cd))

  idx <- which(
    as.character(cd$cell_type) == pop_label &
      as.character(cd$sample) %in% c(TIME_REF, TIME_CASE)
  )

  if (length(idx) == 0) stop("No cells found for population: ", pop_label)

  se_sub <- se[, idx]
  assay_use <- if ("counts" %in% assayNames(se_sub)) "counts" else assayNames(se_sub)[1]

  so <- suppressWarnings(as.Seurat(se_sub, counts = assay_use, data = NULL))
  da <- if ("RNA" %in% SeuratObject::Assays(so)) "RNA" else SeuratObject::Assays(so)[1]
  DefaultAssay(so) <- da
  so <- NormalizeData(so, assay = da, verbose = FALSE)

  so$sample <- as.character(colData(se_sub)$sample)

  sig_hits <- match_signature_genes(sig_genes, rownames(so))
  if (length(sig_hits) == 0) {
    stop("No matched genes in panel for signature ", sig_label, " and population ", pop_label)
  }

  ctrl_n <- max(1L, min(5L, floor((nrow(so) - length(sig_hits)) / 2)))
  score_name <- paste0(gsub("[^A-Za-z0-9]", "_", sig_label), "_score")

  so <- AddModuleScore(
    so,
    features = list(sig_hits),
    name = score_name,
    ctrl = ctrl_n,
    seed = 1997
  )

  score_col <- paste0(score_name, "1")

  data.frame(
    sample = as.character(so$sample),
    module_score = so[[score_col, drop = TRUE]],
    stringsAsFactors = FALSE
  ) %>%
    filter(sample %in% c(TIME_REF, TIME_CASE)) %>%
    mutate(
      sample_label = factor(TIME_LABELS[sample], levels = c("LCMV 1 wpi", "LCMV 6 wpi"))
    )
}

make_violin <- function(df, pop_label, sig_title_text) {
  d1 <- df$module_score[df$sample_label == "LCMV 1 wpi"]
  d2 <- df$module_score[df$sample_label == "LCMV 6 wpi"]

  pval <- tryCatch(
    wilcox.test(d2, d1, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )

  stars <- pval_stars(pval)
  y_max <- max(df$module_score, na.rm = TRUE)
  y_rng <- diff(range(df$module_score, na.rm = TRUE))
  if (!is.finite(y_rng) || y_rng == 0) y_rng <- 1
  y_txt <- y_max + 0.08 * y_rng

  ggplot(df, aes(x = sample_label, y = module_score, fill = sample_label)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.58, colour = NA) +
    geom_boxplot(
      width = 0.18,
      outlier.size = 0.35,
      outlier.alpha = 0.35,
      colour = "grey20",
      fill = "white"
    ) +
    annotate(
      "text",
      x = 1.5,
      y = y_txt,
      label = stars,
      size = 4.3,
      fontface = "bold",
      colour = "grey20"
    ) +
    scale_fill_manual(values = TIME_COLORS, guide = "none") +
    labs(
      title = paste0(sig_title_text, " score in ", pop_label, ": 6wpi vs 1wpi"),
      subtitle = paste0(
        "Wilcoxon p = ", ifelse(is.na(pval), "NA", formatC(pval, digits = 3, format = "e")),
        " | n(1wpi)=", sum(df$sample_label == "LCMV 1 wpi"),
        " n(6wpi)=", sum(df$sample_label == "LCMV 6 wpi")
      ),
      x = NULL,
      y = "Module score (AddModuleScore)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8, colour = "grey40"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

message("Loading object: ", OBJ_FILE)
if (!file.exists(OBJ_FILE)) stop("Missing object file: ", OBJ_FILE)
se <- readRDS(OBJ_FILE)

up_sig <- read_sig(UP_SIG_FILE)
down_sig <- read_sig(DOWN_SIG_FILE)

run_one <- function(pop_label, sig_genes, sig_title_text, out_stub) {
  df <- compute_plot_df(se, pop_label, sig_genes, sig_title_text)
  p <- make_violin(df, pop_label, sig_title_text)
  save_fig(p, file.path(OUT_DIR, out_stub))
  message("Saved: ", out_stub, ".pdf/.jpg")
}

# 1) Activated microglia (C1qa)
run_one(
  pop_label = "Activated microglia (C1qa)",
  sig_genes = up_sig,
  sig_title_text = "Upregulated DAM",
  out_stub = "fig_dam_up_score_microglia_6wpi_vs_1wpi"
)
run_one(
  pop_label = "Activated microglia (C1qa)",
  sig_genes = down_sig,
  sig_title_text = "Downregulated DAM",
  out_stub = "fig_dam_down_score_microglia_6wpi_vs_1wpi"
)

# 2) Antigen-presenting myeloid cells (Cd74)
run_one(
  pop_label = "Antigen-presenting myeloid cells (Cd74)",
  sig_genes = up_sig,
  sig_title_text = "Upregulated DAM",
  out_stub = "fig_dam_up_score_mac_6wpi_vs_1wpi"
)
run_one(
  pop_label = "Antigen-presenting myeloid cells (Cd74)",
  sig_genes = down_sig,
  sig_title_text = "Downregulated DAM",
  out_stub = "fig_dam_down_score_mac_6wpi_vs_1wpi"
)

message("Done. Outputs in: ", OUT_DIR)
