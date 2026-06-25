#!/usr/bin/env Rscript

set.seed(1997)

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(Cairo)
})

base_path <- normalizePath(".")
setwd(base_path)

source("scripts/00_palette.R")

OBJ_FILE <- file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds")
ANNOT_CSV <- "ncells_by_sample_lam02_res09_joint_long.csv"
DAM_SIG_FILE <- file.path("outputs", "banksy", "dam_signature_curation", "Upregulated_DAM.csv")

OUT_DIR <- file.path("outputs", "banksy", "immune_niche_volcano_by_celltype")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TARGET_CLUSTERS <- c("Mac (Ctss)", "T cells (Gzmb)", "Microglia (C1qa)")
ID1 <- "LCMV_6wpi"
ID2 <- "LCMV_1wpi"

MIN_CELLS <- 10
MIN_PCT <- 0.05
FC_THRESH <- 0.1
FDR_CUTOFF <- 0.05
FC_CUTOFF <- 0.25
TOP_N_LABEL <- 15

find_cl_col <- function(cd, lam = 0.2, res = 0.9) {
  cn <- colnames(cd)
  lam_pat <- paste0("lam", gsub("\\.", "\\\\.", as.character(lam)))
  res_pat <- paste0("_res", gsub("\\.", "\\\\.", as.character(res)), "$")
  cols <- cn[grepl(lam_pat, cn) & grepl(res_pat, cn)]
  if (length(cols) == 0) stop("Could not find clustering column for lam=", lam, " res=", res)
  cols[1]
}

slug <- function(x) gsub("[^a-z0-9]+", "_", tolower(trimws(x)))

save_fig <- function(p, out_base, width = 7, height = 5.5) {
  CairoPDF(paste0(out_base, ".pdf"), width = width, height = height)
  print(p)
  dev.off()
  CairoJPEG(paste0(out_base, ".jpg"), width = width * 150, height = height * 150, res = 150)
  print(p)
  dev.off()
}

build_volcano <- function(mk, ct_label) {
  mk$neg_log10_fdr <- -log10(pmax(mk$p_val_adj, .Machine$double.eps))
  x_cap <- max(3.5, quantile(abs(mk$avg_log2FC), 0.999, na.rm = TRUE) * 1.05)
  y_cap <- max(5, quantile(mk$neg_log10_fdr, 0.999, na.rm = TRUE) * 1.05)
  mk$x_plot <- pmin(pmax(mk$avg_log2FC, -x_cap), x_cap)
  mk$y_plot <- pmin(mk$neg_log10_fdr, y_cap)

  sig_mk <- mk %>% filter(direction != "ns")
  lab_genes <- if (nrow(sig_mk) > 0) {
    sig_mk %>% arrange(p_val_adj, desc(abs(avg_log2FC))) %>% slice_head(n = TOP_N_LABEL) %>% pull(gene)
  } else {
    mk %>% arrange(p_val_adj) %>% slice_head(n = min(TOP_N_LABEL, 10)) %>% pull(gene)
  }
  mk$label <- ifelse(mk$gene %in% lab_genes, mk$gene, NA_character_)

  n_up <- sum(mk$direction == "up", na.rm = TRUE)
  n_down <- sum(mk$direction == "down", na.rm = TRUE)

  ggplot(mk, aes(x = x_plot, y = y_plot, colour = direction)) +
    geom_point(size = 1.15, alpha = 0.72, stroke = 0) +
    geom_text_repel(
      data = mk %>% filter(!is.na(label)),
      aes(label = label),
      size = 2.5,
      fontface = "italic",
      max.overlaps = 20,
      box.padding = 0.28,
      point.padding = 0.15,
      segment.size = 0.28,
      segment.color = "grey50",
      show.legend = FALSE
    ) +
    geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF), linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    scale_colour_manual(
      values = c(up = "#B2182B", down = "#2166AC", ns = "grey75"),
      labels = c(
        up = paste0("Up in 6wpi (n=", n_up, ")"),
        down = paste0("Down in 6wpi (n=", n_down, ")"),
        ns = "Not significant"
      ),
      name = NULL
    ) +
    scale_x_continuous(limits = c(-x_cap, x_cap), expand = expansion(mult = 0.02)) +
    scale_y_continuous(limits = c(0, y_cap), expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = paste0(ct_label, ": LCMV 6 wpi vs LCMV 1 wpi"),
      subtitle = paste0(nrow(mk), " genes tested | Wilcoxon | FDR<=", FDR_CUTOFF, " | |log2FC|>", FC_CUTOFF),
      x = "log2FC (6wpi / 1wpi)",
      y = expression(-log[10](FDR))
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8, colour = "grey40"),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey93", linewidth = 0.3)
    )
}

message("Loading object: ", OBJ_FILE)
se <- readRDS(OBJ_FILE)
cd <- as.data.frame(colData(se))

if (!all(c("sample") %in% colnames(cd))) stop("Missing sample column in colData")

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
sample_vec <- as.character(cd$sample)

deg_summary <- list()

for (ct in TARGET_CLUSTERS) {
  idx <- which(annotation == ct & sample_vec %in% c(ID1, ID2))
  n1 <- sum(sample_vec[idx] == ID1)
  n2 <- sum(sample_vec[idx] == ID2)

  message("\n--- ", ct, " | ", ID1, " vs ", ID2, " | n=", n1, " vs ", n2, " ---")

  if (n1 < MIN_CELLS || n2 < MIN_CELLS) {
    message("Skipping ", ct, " (insufficient cells)")
    deg_summary[[ct]] <- data.frame(cell_type = ct, n_up = NA_integer_, n_down = NA_integer_, n_total_plot = NA_integer_)
    next
  }

  se_sub <- se[, idx]
  assay_use <- if ("counts" %in% assayNames(se_sub)) "counts" else assayNames(se_sub)[1]

  so <- suppressWarnings(as.Seurat(se_sub, counts = assay_use, data = NULL))
  da <- if ("RNA" %in% SeuratObject::Assays(so)) "RNA" else SeuratObject::Assays(so)[1]
  DefaultAssay(so) <- da

  needs_norm <- tryCatch({
    dm <- GetAssayData(so, assay = da, layer = "data")
    nrow(dm) == 0 || ncol(dm) == 0
  }, error = function(e) TRUE)
  if (needs_norm) so <- NormalizeData(so, assay = da, verbose = FALSE)

  so$sample <- as.character(colData(se_sub)$sample)
  so <- SetIdent(so, value = "sample")

  mk <- FindMarkers(
    so,
    ident.1 = ID1,
    ident.2 = ID2,
    only.pos = FALSE,
    min.pct = MIN_PCT,
    logfc.threshold = FC_THRESH,
    return.thresh = 1,
    test.use = "wilcox",
    verbose = FALSE
  )

  if (is.null(mk) || nrow(mk) == 0) {
    message("No DEG result for ", ct)
    deg_summary[[ct]] <- data.frame(cell_type = ct, n_up = NA_integer_, n_down = NA_integer_, n_total_plot = NA_integer_)
    next
  }

  mk$gene <- rownames(mk)
  mk$direction <- case_when(
    mk$p_val_adj <= FDR_CUTOFF & mk$avg_log2FC > FC_CUTOFF ~ "up",
    mk$p_val_adj <= FDR_CUTOFF & mk$avg_log2FC < -FC_CUTOFF ~ "down",
    TRUE ~ "ns"
  )
  mk <- mk %>% arrange(p_val_adj, desc(abs(avg_log2FC)))

  n_up <- sum(mk$direction == "up", na.rm = TRUE)
  n_down <- sum(mk$direction == "down", na.rm = TRUE)

  deg_summary[[ct]] <- data.frame(cell_type = ct, n_up = n_up, n_down = n_down, n_total_plot = n_up + n_down)

  ct_slug <- slug(ct)
  contrast_slug <- paste0(slug(ID1), "_vs_", slug(ID2))

  out_csv <- file.path(OUT_DIR, paste0("DEG_", ct_slug, "__", contrast_slug, ".csv"))
  write.csv(mk, out_csv, row.names = FALSE)

  p_vol <- build_volcano(mk, ct)
  out_base <- file.path(OUT_DIR, paste0("fig_volcano_", ct_slug, "__", contrast_slug))
  save_fig(p_vol, out_base, width = 7, height = 5.5)

  message("Saved DEG + volcano for ", ct)
}

summary_df <- bind_rows(deg_summary)
write.csv(summary_df, file.path(OUT_DIR, "ndeg_summary_6wpi_vs_1wpi_selected_clusters.csv"), row.names = FALSE)

# nDEG summary figure
sum_plot_df <- summary_df %>%
  tidyr::pivot_longer(cols = c("n_up", "n_down"), names_to = "direction", values_to = "n_deg") %>%
  mutate(
    direction = factor(direction, levels = c("n_up", "n_down"), labels = c("Up in 6wpi", "Down in 6wpi")),
    cell_type = factor(cell_type, levels = rev(TARGET_CLUSTERS))
  )

p_ndeg <- ggplot(sum_plot_df, aes(x = n_deg, y = cell_type, fill = direction)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.72, colour = "white") +
  scale_fill_manual(values = c("Up in 6wpi" = "#B2182B", "Down in 6wpi" = "#2166AC"), name = NULL) +
  labs(
    title = "Number of DEGs per cluster (LCMV 6 wpi vs LCMV 1 wpi)",
    subtitle = "Wilcoxon test | FDR<=0.05 and |log2FC|>0.25",
    x = "Number of DEGs",
    y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, colour = "grey40"),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

save_fig(p_ndeg, file.path(OUT_DIR, "fig_ndeg_summary_6wpi_vs_1wpi_selected_clusters"), width = 7.4, height = 4.4)

# DAM signature in Microglia (C1qa): 6wpi vs 1wpi
if (file.exists(DAM_SIG_FILE)) {
  dam_sig <- read.csv(DAM_SIG_FILE, stringsAsFactors = FALSE)
  dam_genes <- unique(as.character(dam_sig$gene))

  mg_label <- "Microglia (C1qa)"
  mg_idx <- which(annotation == mg_label & sample_vec %in% c(ID1, ID2))

  if (length(mg_idx) > 0) {
    se_mg <- se[, mg_idx]
    assay_use <- if ("counts" %in% assayNames(se_mg)) "counts" else assayNames(se_mg)[1]

    so_mg <- suppressWarnings(as.Seurat(se_mg, counts = assay_use, data = NULL))
    da_mg <- if ("RNA" %in% SeuratObject::Assays(so_mg)) "RNA" else SeuratObject::Assays(so_mg)[1]
    DefaultAssay(so_mg) <- da_mg
    so_mg <- NormalizeData(so_mg, assay = da_mg, verbose = FALSE)

    # Match signature genes to panel in title-case to handle species casing differences
    panel_genes <- rownames(so_mg)
    panel_title <- tools::toTitleCase(tolower(panel_genes))
    panel_map <- setNames(panel_genes, panel_title)
    sig_title <- tools::toTitleCase(tolower(dam_genes))
    sig_hits <- unique(as.character(panel_map[sig_title[sig_title %in% panel_title]]))
    sig_hits <- sig_hits[!is.na(sig_hits)]

    if (length(sig_hits) > 0) {
      ctrl_n <- max(1L, min(5L, floor((nrow(so_mg) - length(sig_hits)) / 2)))
      so_mg <- AddModuleScore(
        so_mg,
        features = list(sig_hits),
        name = "Upregulated_DAM",
        ctrl = ctrl_n,
        seed = 1997
      )

      dam_df <- data.frame(
        sample = as.character(colData(se_mg)$sample),
        dam_score = so_mg$Upregulated_DAM1,
        stringsAsFactors = FALSE
      ) %>%
        filter(sample %in% c(ID1, ID2)) %>%
        mutate(sample = factor(sample, levels = c(ID2, ID1), labels = c("LCMV 1 wpi", "LCMV 6 wpi")))

      pval <- tryCatch(
        wilcox.test(
          dam_df$dam_score[dam_df$sample == "LCMV 6 wpi"],
          dam_df$dam_score[dam_df$sample == "LCMV 1 wpi"],
          exact = FALSE
        )$p.value,
        error = function(e) NA_real_
      )

      p_dam <- ggplot(dam_df, aes(x = sample, y = dam_score, fill = sample)) +
        geom_violin(trim = TRUE, scale = "width", alpha = 0.58, colour = NA) +
        geom_boxplot(width = 0.18, outlier.size = 0.35, outlier.alpha = 0.35, colour = "grey20", fill = "white") +
        scale_fill_manual(values = c("LCMV 1 wpi" = "#56B4E9", "LCMV 6 wpi" = "#F28E2B"), guide = "none") +
        labs(
          title = "DAM signature score in Microglia (C1qa): LCMV 6 wpi vs 1 wpi",
          subtitle = paste0("Wilcoxon p = ", ifelse(is.na(pval), "NA", formatC(pval, digits = 3, format = "e"))),
          x = NULL,
          y = "DAM module score (AddModuleScore)"
        ) +
        theme_bw(base_size = 10) +
        theme(
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 8, colour = "grey40"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank()
        )

      save_fig(
        p_dam,
        file.path("outputs", "banksy", "microglia_dam_niche", "fig_dam_signature_microglia_c1qa_6wpi_vs_1wpi_violin"),
        width = 6.2,
        height = 4.8
      )

      write.csv(
        dam_df,
        file.path("outputs", "banksy", "microglia_dam_niche", "dam_signature_microglia_c1qa_6wpi_vs_1wpi_scores.csv"),
        row.names = FALSE
      )
    }
  }
}

message("\nDone. Outputs in:")
message(" - ", OUT_DIR)
message(" - outputs/banksy/microglia_dam_niche")
