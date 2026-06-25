#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(tidyverse)
})

base_path <- normalizePath(".")  # Run this script from the project root directory
setwd(base_path)

# Palette et ordre des annotations (fichier partagé)
source("scripts/00_palette.R")

# -------------------------------------------------------
# Default parameters
# -------------------------------------------------------
IMMUNE_LABEL <- "Immune (Acod1)"
FOLDER_NAME  <- NULL
RES          <- 0.3
TOP_N_UNION  <- 10
FC_THRESH    <- 0.05
MIN_PCT      <- 0.10
ANNOT_FILE   <- NULL   # optional: CSV with columns 'cluster' and 'annotation'

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  for (i in seq_along(args)) {
    a <- args[i]
    if (a %in% c("--label", "-l") && (i + 1) <= length(args)) IMMUNE_LABEL <- args[i + 1]
    if (a == "--folder"           && (i + 1) <= length(args)) FOLDER_NAME  <- args[i + 1]
    if (a == "--res"              && (i + 1) <= length(args)) RES          <- as.numeric(args[i + 1])
    if (a == "--top"              && (i + 1) <= length(args)) TOP_N_UNION  <- as.integer(args[i + 1])
    if (a == "--annot"            && (i + 1) <= length(args)) ANNOT_FILE   <- args[i + 1]
  }
}

slugify   <- function(x) gsub("^_+|_+$", "", gsub("[^a-z0-9]+", "_", tolower(trimws(x))))
safe_label <- if (!is.null(FOLDER_NAME)) tolower(trimws(FOLDER_NAME)) else slugify(IMMUNE_LABEL)
res_id    <- gsub("\\.", "", sprintf("%.1f", RES))

# Auto-detect annotation file if not provided
# Priorité : CSV de mapping propre exporté par 17_annotate_immune_object.R
if (is.null(ANNOT_FILE)) {
  clean_map <- file.path("outputs", "banksy", safe_label, "annotation_map_immune_res03.csv")
  if (file.exists(clean_map)) {
    ANNOT_FILE <- clean_map
    message("  Auto-detected annotation map: ", ANNOT_FILE)
  } else {
    # Fallback : DEGs CSV annoté manuellement (colonnes cluster;annotation)
    auto_annot <- file.path("outputs", "banksy", safe_label, "analysis", "DEGs",
                             paste0("DEGs_all_res", res_id, ".csv"))
    if (file.exists(auto_annot)) {
      ANNOT_FILE <- auto_annot
      message("  Auto-detected annotation file: ", ANNOT_FILE)
    }
  }
}

out_dir   <- file.path("outputs", "banksy", safe_label, "domain_DEG_evolution", paste0("res", res_id))
csv_dir   <- file.path(out_dir, "DEG_csv")
fig_dir   <- file.path(out_dir, "heatmaps")
for (d in c(csv_dir, fig_dir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

message("=== 16_domain_DEG_evolution.R | ", safe_label, " | res=", RES, " ===")

# -------------------------------------------------------
# Load annotation map (optional --annot)
# -------------------------------------------------------
annot_map <- NULL
if (!is.null(ANNOT_FILE) && file.exists(ANNOT_FILE)) {
  annot_df  <- read.csv(ANNOT_FILE, sep = ";", stringsAsFactors = FALSE)
  annot_map <- setNames(trimws(annot_df$annotation), trimws(annot_df$cluster))
  annot_map <- annot_map[!duplicated(names(annot_map))]
  annot_map <- annot_map[nchar(annot_map) > 0]  # drop empty annotations
  message("  Annotation map loaded: ", length(annot_map), " domains")
} else if (!is.null(ANNOT_FILE)) {
  warning("Annotation file not found: ", ANNOT_FILE)
}

# -------------------------------------------------------
# Load BANKSY object
# -------------------------------------------------------
cands <- c(
  file.path("objects", paste0("07_immune_banksy_lam02_", safe_label, ".rds")),
  file.path("objects", "07_immune_banksy_lam02_after_bloc1.rds")
)
banksy_obj <- cands[file.exists(cands)][1]
if (is.na(banksy_obj)) stop("BANKSY object not found for: ", safe_label)
se <- readRDS(banksy_obj)
message("  ", ncol(se), " cells loaded from: ", banksy_obj)

assay_name    <- if ("counts" %in% assayNames(se)) "counts" else assayNames(se)[1]
so            <- as.Seurat(se, counts = assay_name, data = NULL)
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

# -------------------------------------------------------
# Find cluster column for requested resolution
# -------------------------------------------------------
cd      <- SummarizedExperiment::colData(se)
res_txt <- sprintf("%.1f", RES)
res_pat <- paste0("_res", gsub("\\.", "\\\\.", res_txt), "$")
lam_pat <- "lam0\\.2|lam02"
cl_cols <- colnames(cd)[grepl(lam_pat, colnames(cd)) & grepl(res_pat, colnames(cd))]
if (length(cl_cols) == 0) {
  stop("No cluster column for res=", RES, " | available: ",
       paste(colnames(cd)[grepl("lam", colnames(cd))], collapse = ", "))
}
cl_col <- cl_cols[1]
message("  Cluster column: ", cl_col)

so@meta.data[["banksy_domain"]] <- factor(
  paste0("Domain_", as.character(cd[[cl_col]])),
  levels = paste0("Domain_", sort(unique(as.integer(cd[[cl_col]]))))
)

# -------------------------------------------------------
# Sample info
# -------------------------------------------------------
sample_vec    <- as.character(so@meta.data$sample)
avail_samples <- unique(sample_vec)
if (!"LCMV_6wpi" %in% avail_samples) stop("LCMV_6wpi not found in object")

sample_order <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi", "mock_6wpi")
present      <- sample_order[sample_order %in% avail_samples]
comparisons  <- setdiff(present, c("LCMV_6wpi", "mock_6wpi"))
comp_tags    <- gsub("LCMV_", "", comparisons)
comp_tags    <- gsub("_", " ", comp_tags)
comp_slugs   <- gsub("[^a-z0-9]", "_", tolower(comp_tags))

domain_levels <- levels(so@meta.data[["banksy_domain"]])
message("  Domains: ", paste(domain_levels, collapse = ", "))
message("  Comparisons: 6wpi vs ", paste(comparisons, collapse = ", "))

# -------------------------------------------------------
# Create combined ident: domain__sample
# -------------------------------------------------------
so@meta.data[["sample_domain"]] <- paste0(so@meta.data[["banksy_domain"]], "__", sample_vec)
so <- SetIdent(so, value = "sample_domain")

# -------------------------------------------------------
# Accumulateur pour la section résumé nDEG
summary_records <- list()

# -------------------------------------------------------
# Loop over domains
# -------------------------------------------------------
for (dom in domain_levels) {
  message("\n=== Domain: ", dom, " ===")

  cells_6wpi <- which(so@meta.data[["banksy_domain"]] == dom & sample_vec == "LCMV_6wpi")
  if (length(cells_6wpi) < 1) {
    message("  Skipping: only ", length(cells_6wpi), " cells in LCMV_6wpi")
    next
  }

  dom_slug      <- gsub("[^a-z0-9]", "_", tolower(dom))
  dom_label     <- if (!is.null(annot_map) && dom %in% names(annot_map)) annot_map[[dom]] else dom
  dom_label_slug <- gsub("[^a-z0-9]", "_", tolower(dom_label))
  all_comp_deg  <- list()

  # --- PART 1: FindMarkers per comparison ---
  for (i in seq_along(comparisons)) {
    comp      <- comparisons[i]
    comp_slug <- comp_slugs[i]
    comp_tag  <- comp_tags[i]

    cells_comp <- which(so@meta.data[["banksy_domain"]] == dom & sample_vec == comp)
    if (length(cells_comp) < 10) {
      message("  Skip vs ", comp, ": only ", length(cells_comp), " cells")
      next
    }

    id1 <- paste0(dom, "__LCMV_6wpi")
    id2 <- paste0(dom, "__", comp)

    mk <- tryCatch(
      FindMarkers(so, ident.1 = id1, ident.2 = id2,
                  only.pos = FALSE, min.pct = MIN_PCT,
                  logfc.threshold = FC_THRESH, return.thresh = 1,
                  test.use = "wilcox", verbose = FALSE),
      error = function(e) { message("  ERROR vs ", comp, ": ", conditionMessage(e)); NULL }
    )
    if (is.null(mk) || nrow(mk) == 0) {
      message("  No markers vs ", comp)
      next
    }

    mk$gene       <- rownames(mk)
    mk$direction  <- ifelse(mk$avg_log2FC > 0, "up_6wpi", "down_6wpi")
    mk$comparison <- paste0("6wpi_vs_", comp_slug)
    mk$domain     <- dom
    mk <- mk %>%
      filter(p_val_adj < 0.05) %>%
      arrange(desc(avg_log2FC)) %>%
      select(domain, gene, avg_log2FC, direction, pct.1, pct.2, p_val, p_val_adj, comparison)

    all_comp_deg[[comp]] <- mk
    n_up   <- sum(mk$avg_log2FC > 0)
    n_down <- sum(mk$avg_log2FC < 0)
    message("  vs ", comp_tag, ": ", nrow(mk), " DEGs (", n_up, " up / ", n_down, " down)")

    fname <- file.path(csv_dir, paste0("DEG_", dom_slug, "_", dom_label_slug, "_6wpi_vs_", comp_slug, ".csv"))
    tryCatch(
      write.csv(mk, fname, row.names = FALSE),
      error = function(e) message("  WARNING write: ", conditionMessage(e))
    )
  }

  # Accumulate nDEG counts for summary figure
  for (i in seq_along(comparisons)) {
    comp_i <- comparisons[i]
    if (comp_i %in% names(all_comp_deg)) {
      mk_i <- all_comp_deg[[comp_i]]
      summary_records[[length(summary_records) + 1]] <- data.frame(
        annotation = dom_label,
        comparison = paste0("6wpi vs ", comp_tags[i]),
        n_up       = sum(mk_i$avg_log2FC > 0),
        n_down     = sum(mk_i$avg_log2FC < 0),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(all_comp_deg) == 0) {
    message("  No significant DEGs for any comparison, skip heatmap")
    next
  }

  # --- PART 2: Heatmap evolution ---
  all_deg_df <- do.call(rbind, all_comp_deg)

  top_genes_union <- all_deg_df %>%
    arrange(desc(abs(avg_log2FC))) %>%
    distinct(gene, .keep_all = TRUE) %>%
    slice_head(n = TOP_N_UNION * 2) %>%
    pull(gene)

  message("  Union top genes for heatmap: ", length(top_genes_union))
  if (length(top_genes_union) < 2) {
    message("  Too few genes, skip heatmap")
    next
  }

  # Build log2FC matrix from DEG results (6wpi vs each comparison)
  fc_col_names <- paste0("6wpi vs ", comp_tags)
  fc_mat <- matrix(0, nrow = length(top_genes_union), ncol = length(comparisons),
                   dimnames = list(top_genes_union, fc_col_names))
  for (j in seq_along(comparisons)) {
    cj <- comparisons[j]
    if (cj %in% names(all_comp_deg)) {
      df_j  <- all_comp_deg[[cj]]
      idx   <- match(top_genes_union, df_j$gene)
      valid <- !is.na(idx)
      fc_mat[valid, j] <- df_j$avg_log2FC[idx[valid]]
    }
  }

  # Row direction (kept for summary CSV only)
  row_dir <- setNames(sapply(rownames(fc_mat), function(g) {
    n_up   <- sum(sapply(all_comp_deg, function(df) g %in% df$gene[df$direction == "up_6wpi"]))
    n_down <- sum(sapply(all_comp_deg, function(df) g %in% df$gene[df$direction == "down_6wpi"]))
    if (n_up > n_down) "up_6wpi"
    else if (n_down > n_up) "down_6wpi"
    else {
      fc_vals <- sapply(all_comp_deg, function(df) {
        idx <- which(df$gene == g)
        if (length(idx) == 0) NA_real_ else df$avg_log2FC[idx[1]]
      })
      fc_vals <- fc_vals[!is.na(fc_vals)]
      if (length(fc_vals) == 0 || fc_vals[which.max(abs(fc_vals))] < 0) "down_6wpi" else "up_6wpi"
    }
  }), rownames(fc_mat))

  # Order rows: up_6wpi block first, then down_6wpi, then NS
  gene_order <- c(
    names(row_dir)[row_dir == "up_6wpi"],
    names(row_dir)[row_dir == "down_6wpi"],
    names(row_dir)[row_dir == "NS"]
  )
  fc_mat_ord <- fc_mat[gene_order, , drop = FALSE]

  # --- Summary table: log2FC per comparison + significance ---
  summary_wide <- tibble(gene = gene_order, direction_majority = row_dir[gene_order])
  for (j in seq_along(comparisons)) {
    cj <- comparisons[j]; sj <- comp_slugs[j]
    if (cj %in% names(all_comp_deg)) {
      df_j <- all_comp_deg[[cj]] %>% select(gene, avg_log2FC, p_val_adj)
      colnames(df_j)[2:3] <- c(paste0("log2FC_vs_", sj), paste0("padj_vs_", sj))
      summary_wide <- left_join(summary_wide, df_j, by = "gene")
    } else {
      summary_wide[[paste0("log2FC_vs_", sj)]] <- NA_real_
      summary_wide[[paste0("padj_vs_", sj)]]   <- NA_real_
    }
  }
  padj_cols <- grep("^padj_", colnames(summary_wide), value = TRUE)
  summary_wide$sig_in_comparisons <- apply(
    summary_wide[, padj_cols, drop = FALSE], 1,
    function(x) {
      idx <- which(!is.na(x) & x < 0.05)
      if (length(idx) == 0) "none" else paste(comp_slugs[idx], collapse = "|")
    }
  )
  summary_wide$domain      <- dom
  summary_wide$annotation  <- dom_label
  write.csv(summary_wide,
            file.path(csv_dir, paste0("summary_genes_", dom_slug, "_", dom_label_slug, ".csv")),
            row.names = FALSE)
  message("  Summary table: ", nrow(summary_wide), " genes saved")

  # Une heatmap par comparaison
  for (j in seq_along(comparisons)) {
    cj      <- comparisons[j]
    cj_tag  <- comp_tags[j]
    cj_slug <- comp_slugs[j]

    if (!cj %in% names(all_comp_deg)) next

    deg_j <- all_comp_deg[[cj]]

    # Top gènes pour cette comparaison uniquement, triés par |log2FC|
    top_genes_j <- deg_j %>%
      arrange(desc(abs(avg_log2FC))) %>%
      distinct(gene, .keep_all = TRUE) %>%
      slice_head(n = TOP_N_UNION) %>%
      pull(gene)

    if (length(top_genes_j) < 2) next

    # Barplot log2FC (ggplot2)
    plot_df <- deg_j %>%
      filter(gene %in% top_genes_j) %>%
      mutate(
        gene    = factor(gene, levels = gene[order(avg_log2FC)]),
        sig     = case_when(
          p_val_adj < 0.001 ~ "***",
          p_val_adj < 0.01  ~ "**",
          p_val_adj < 0.05  ~ "*",
          TRUE              ~ ""
        ),
        label_x = ifelse(avg_log2FC >= 0, avg_log2FC + 0.05, avg_log2FC - 0.05),
        hjust   = ifelse(avg_log2FC >= 0, 0, 1)
      )

    p_bar <- ggplot(plot_df, aes(x = avg_log2FC, y = gene, fill = direction)) +
      geom_col(width = 0.7) +
      geom_vline(xintercept = 0, linewidth = 0.4, color = "grey30") +
      geom_text(aes(x = label_x, label = sig, hjust = hjust),
                size = 3, vjust = 0.75) +
      scale_fill_manual(
        values = c(up_6wpi = "#B2182B", down_6wpi = "#2166AC"),
        labels = c(up_6wpi = "Up at 6wpi", down_6wpi = "Down at 6wpi"),
        name   = NULL
      ) +
      labs(
        x        = bquote(log[2]*"FC  (6wpi vs "*.(cj_tag)*")"),
        y        = NULL,
        title    = dom_label,
        subtitle = paste0("6wpi vs ", cj_tag, "  |  top ", length(top_genes_j), " DEGs by |log2FC|")
      ) +
      theme_classic(base_size = 11) +
      theme(
        plot.title      = element_text(face = "bold", size = 12),
        plot.subtitle   = element_text(color = "grey50", size = 9),
        axis.text.y     = element_text(size = 9, face = "italic"),
        axis.text.x     = element_text(size = 9),
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"),
        legend.text     = element_text(size = 9)
      )

    fig_h      <- max(4, length(top_genes_j) * 0.4 + 2)
    fname_base <- file.path(fig_dir, paste0("barplot_", dom_slug, "_", dom_label_slug, "_vs_", cj_slug))
    ggsave(paste0(fname_base, ".pdf"), p_bar, width = 5, height = fig_h)
    ggsave(paste0(fname_base, ".jpg"), p_bar, width = 5, height = fig_h, dpi = 300)
    message("  Barplot saved: barplot_", dom_slug, "_", dom_label_slug, "_vs_", cj_slug)
  }
}

# -------------------------------------------------------
# Section finale : figure récapitulative nDEG par domaine
# -------------------------------------------------------
if (length(summary_records) > 0) {
  message("\n=== Section finale: fig_ndeg_summary_by_domain ===")

  summary_df <- bind_rows(summary_records)

  # Ordre par nDEG total croissant → après coord_flip, le plus haut est en haut
  annot_order <- summary_df %>%
    group_by(annotation) %>%
    summarise(total = sum(n_up + n_down), .groups = "drop") %>%
    arrange(total) %>%
    pull(annotation)

  # Couleurs GLOBAL_PALETTE pour les labels d'axe
  label_colors <- GLOBAL_PALETTE[annot_order]
  label_colors[is.na(label_colors)] <- "grey40"

  # Long format : n_up positif, n_down négatif
  plot_ndeg <- summary_df %>%
    mutate(annotation = factor(annotation, levels = annot_order)) %>%
    pivot_longer(cols = c(n_up, n_down),
                 names_to = "direction", values_to = "n") %>%
    mutate(
      count     = ifelse(direction == "n_down", -n, n),
      direction = factor(direction, levels = c("n_up", "n_down"),
                         labels = c("Up at 6wpi", "Down at 6wpi"))
    )

  p_ndeg <- ggplot(plot_ndeg, aes(x = annotation, y = count, fill = direction)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey30") +
    facet_wrap(~ comparison, ncol = length(unique(plot_ndeg$comparison))) +
    scale_fill_manual(
      values = c("Up at 6wpi" = "#B2182B", "Down at 6wpi" = "#2166AC"),
      name   = NULL
    ) +
    scale_y_continuous(
      labels = function(x) abs(x),
      name   = "Number of DEGs (FDR \u2264 0.05)"
    ) +
    labs(x = NULL,
         title = paste0("DEG evolution by domain  |  ", safe_label)) +
    coord_flip() +
    theme_bw(base_size = 10) +
    theme(
      plot.title         = element_text(face = "bold", size = 11),
      axis.text.y        = element_text(size = 9),
      strip.text         = element_text(face = "bold", size = 10),
      legend.position    = "bottom",
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank()
    )

  # Labels d'axe colorés par GLOBAL_PALETTE via ggtext (si disponible)
  if (requireNamespace("ggtext", quietly = TRUE)) {
    md_labels <- setNames(
      vapply(as.character(annot_order), function(a) {
        col <- if (!is.na(label_colors[a])) label_colors[a] else "grey40"
        sprintf("<span style='color:%s'>%s</span>", col, a)
      }, character(1)),
      as.character(annot_order)
    )
    p_ndeg <- p_ndeg +
      scale_x_discrete(labels = md_labels) +
      theme(axis.text.y = ggtext::element_markdown(size = 9))
  }

  n_comps  <- length(unique(plot_ndeg$comparison))
  ndeg_w   <- max(5, n_comps * 3.5)
  ndeg_h   <- max(3, length(annot_order) * 0.5 + 1.5)

  fname_ndeg <- file.path(fig_dir, "fig_ndeg_summary_by_domain")
  ggsave(paste0(fname_ndeg, ".pdf"), p_ndeg, width = ndeg_w, height = ndeg_h)
  ggsave(paste0(fname_ndeg, ".jpg"), p_ndeg, width = ndeg_w, height = ndeg_h, dpi = 300)
  message("  Saved: fig_ndeg_summary_by_domain.pdf / .jpg")

  write.csv(summary_df, file.path(fig_dir, "ndeg_summary_by_domain.csv"), row.names = FALSE)
  message("  Saved: ndeg_summary_by_domain.csv")
} else {
  message("\n  No DEGs found for any domain — summary plot skipped.")
}

message("\nDone. Outputs in: ", out_dir)
