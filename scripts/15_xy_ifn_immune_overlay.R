#!/usr/bin/env Rscript
# =============================================================
# Script: 15_xy_ifn_immune_overlay.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Description:
#   Figure 1 — XY spatial maps (one 2×2 grid, one panel per sample)
#     with three layers:
#       - Grey background : all other cells  (size 0.1, alpha 0.3)
#       - Red             : Immune (Acod1)   = Domain_13  (size 0.3, alpha 0.9)
#       - Blue            : IFN responsive (Ifit1) = Domain_16  (size 0.3, alpha 0.9)
#     Source: 04_banksy_joint_lam08_after_bloc3.rds
#             ncells_by_sample_lam02_res09_joint_long.csv
#
#   Figure 2 — % Ifng+ cells by cell type and condition
#     Check whether Ifng is in the MERFISH panel.
#     If yes → grouped barplot (x = cell_type, y = % Ifng+, color = sample).
#     If no  → report + barplot of IFN-γ proximal genes (Ifit1, Cxcl10,
#               Irf7, Stat1, Gbp2) expressed in any cell (pct_expressing).
#     Source: 03_all_clustered.rds
#
# Output: outputs/banksy/ifn_immune_overlay/
#           fig1_xy_ifn_immune_overlay.pdf / .jpg
#           fig2_ifng_pct_by_celltype.pdf  / .jpg   (or fig2_ifng_proximal.pdf/.jpg)
# =============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(tidyverse)
  library(patchwork)
  library(ggplot2)
})

base_path <- normalizePath(".")  # Run this script from the project root directory
setwd(base_path)

out_dir <- file.path("outputs", "banksy", "ifn_immune_overlay")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

SAMPLE_ORDER <- c("mock_6wpi", "LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")

# Domains of interest (BANKSY lam=0.2, res=0.9)
DOMAIN_IMMUNE <- "Domain_13"   # Immune (Acod1)
DOMAIN_IFN    <- "Domain_16"   # IFN responsive (Ifit1)

COLOR_IMMUNE <- "#D7191C"   # red
COLOR_IFN    <- "#2C7BB6"   # blue
COLOR_OTHER  <- "#CCCCCC"   # light grey

# Gènes IFN à tester — tous confirmés présents dans le panel MERFISH
# Ordonnés par importance biologique (facteur de transcription > source > ISGs)
# Stat1  : facteur de transcription central de la voie IFN
# Ifng   : source de l'IFN-γ (producteurs)
# Ifit1  : ISG canonique, nomme Domain_16
# Irf7   : facteur de transcription ISG
# Cxcl10 : chimiokine induite par IFN-γ
# Gbp2   : GTPase induite par IFN-γ
IFN_PROXIMAL <- c("Stat1", "Ifng", "Ifit1", "Irf7", "Cxcl10", "Gbp2")

# -------------------------------------------------------
# Helpers
# -------------------------------------------------------
find_cl_col <- function(se, lam, res) {
  cn      <- colnames(SummarizedExperiment::colData(se))
  lam_pat <- paste0("lam", gsub("\\.", "\\\\.", as.character(lam)))
  res_pat <- paste0("_res", gsub("\\.", "\\\\.", sprintf("%.1f", res)), "$")
  cols    <- cn[grepl(lam_pat, cn) & grepl(res_pat, cn)]
  if (length(cols) == 0) stop("Clustering column not found for lam=", lam, " res=", res)
  cols[1]
}

save_plot <- function(p, base, w, h) {
  pdf_file <- paste0(base, ".pdf")
  jpg_file <- paste0(base, ".jpg")
  ggsave(pdf_file, plot = p, width = w, height = h, device = cairo_pdf)
  ggsave(jpg_file, plot = p, width = w, height = h, dpi = 450)
  message("  Saved: ", pdf_file)
  message("  Saved: ", jpg_file)
}

# ================================================================
# FIGURE 1 — XY spatial overlay: IFN responsive + Immune (Acod1)
# ================================================================
message("\n=== Figure 1: XY spatial overlay ===")
message("Loading: objects/04_banksy_joint_lam08_after_bloc3.rds")
se <- readRDS(file.path("objects", "04_banksy_joint_lam08_after_bloc3.rds"))
message("  ", ncol(se), " cells")

cluster_col <- find_cl_col(se, 0.2, 0.9)
message("  Cluster column: ", cluster_col)

cd      <- as.data.frame(SummarizedExperiment::colData(se))
xy      <- as.data.frame(SpatialExperiment::spatialCoords(se))

if (!all(c("sdimx", "sdimy") %in% colnames(xy))) {
  stop("Spatial coordinates must contain columns 'sdimx' and 'sdimy'")
}

# Build per-cell domain label
domain_labels <- paste0("Domain_", as.character(cd[[cluster_col]]))

# Assign layer: immune / ifn / other
layer <- dplyr::case_when(
  domain_labels == DOMAIN_IMMUNE ~ "Immune (Acod1)",
  domain_labels == DOMAIN_IFN    ~ "IFN responsive (Ifit1)",
  TRUE                           ~ "Other"
)

plot_df <- tibble(
  x      = xy$sdimx,
  y      = xy$sdimy,
  sample = as.character(cd$sample),
  layer  = factor(layer, levels = c("Other", "Immune (Acod1)", "IFN responsive (Ifit1)"))
)

# Keep only known samples and order them
samples_present <- unique(plot_df$sample)
sample_levels   <- c(
  SAMPLE_ORDER[SAMPLE_ORDER %in% samples_present],
  setdiff(samples_present, SAMPLE_ORDER)
)
plot_df$sample <- factor(plot_df$sample, levels = sample_levels)

# Counts summary
cat("\nCell counts per layer:\n")
print(
  plot_df %>%
    filter(layer != "Other") %>%
    count(sample, layer) %>%
    spread(layer, n, fill = 0)
)

# ------------------------------------------------------------------
# Convex hull per sample — tissue silhouette (like script 10)
# ------------------------------------------------------------------
hull_df <- do.call(rbind, lapply(levels(plot_df$sample), function(s) {
  d <- plot_df %>% filter(sample == s)
  if (nrow(d) < 3) return(NULL)
  h <- d[chull(d$x, d$y), c("x", "y")]
  h$sample <- s
  h
}))
if (!is.null(hull_df) && nrow(hull_df) > 0) {
  hull_df$sample <- factor(hull_df$sample, levels = levels(plot_df$sample))
} else {
  hull_df <- tibble(x = numeric(), y = numeric(),
                    sample = factor(character(), levels = levels(plot_df$sample)))
}

# ------------------------------------------------------------------
# One ggplot panel per sample
# ------------------------------------------------------------------
make_sample_panel <- function(sname, df) {
  d_samp  <- df %>% filter(sample == sname)
  d_other <- d_samp %>% filter(layer == "Other")
  d_imm   <- d_samp %>% filter(layer == "Immune (Acod1)")
  d_ifn   <- d_samp %>% filter(layer == "IFN responsive (Ifit1)")

  xr   <- range(d_samp$x, na.rm = TRUE)
  yr   <- range(d_samp$y, na.rm = TRUE)
  xpad <- max(diff(xr) * 0.08, 50)
  ypad <- max(diff(yr) * 0.08, 50)

  nice_name <- dplyr::recode(
    sname,
    mock_6wpi  = "Mock 6 wpi",
    LCMV_1wpi  = "LCMV 1 wpi",
    LCMV_3wpi  = "LCMV 3 wpi",
    LCMV_6wpi  = "LCMV 6 wpi"
  )

  hull_sample <- hull_df %>% filter(sample == sname)

  ggplot() +
    # tissue silhouette
    (if (nrow(hull_sample) > 0)
      geom_polygon(
        data    = hull_sample,
        mapping = aes(x = x, y = y),
        fill    = "#f8f8f8",
        color   = NA
      )
    else NULL) +
    # grey background cells
    geom_point(
      data    = d_other,
      mapping = aes(x = x, y = y),
      color   = COLOR_OTHER,
      size    = 0.1,
      alpha   = 0.3,
      stroke  = 0,
      show.legend = FALSE
    ) +
    # Immune (Acod1) — red
    geom_point(
      data    = d_imm,
      mapping = aes(x = x, y = y, color = layer),
      size    = 0.3,
      alpha   = 0.9,
      stroke  = 0
    ) +
    # IFN responsive — blue
    geom_point(
      data    = d_ifn,
      mapping = aes(x = x, y = y, color = layer),
      size    = 0.3,
      alpha   = 0.9,
      stroke  = 0
    ) +
    scale_color_manual(
      name   = NULL,
      values = c(
        "Immune (Acod1)"         = COLOR_IMMUNE,
        "IFN responsive (Ifit1)" = COLOR_IFN
      ),
      drop = FALSE
    ) +
    coord_fixed(
      xlim   = c(xr[1] - xpad, xr[2] + xpad),
      ylim   = c(yr[1] - ypad, yr[2] + ypad),
      expand = FALSE
    ) +
    labs(
      title = nice_name,
      x     = "X coordinate (um)",
      y     = "Y coordinate (um)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title  = element_text(face = "bold", size = 10, hjust = 0.5),
      axis.title  = element_text(size = 10),
      axis.text   = element_text(size = 8),
      axis.line   = element_line(linewidth = 0.5, color = "black"),
      axis.ticks  = element_line(linewidth = 0.4, color = "black"),
      plot.margin = margin(6, 8, 6, 8)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
}

panels <- lapply(sample_levels, make_sample_panel, df = plot_df)
panels <- panels[!vapply(panels, is.null, logical(1))]

fig1 <- wrap_plots(panels, ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Spatial distribution of IFN-responsive and immune cells during LCMV infection",
    subtitle = "BANKSY spatial domains (\u03bb = 0.2, res = 0.9) | Immune (Acod1): Domain 13, IFN-responsive (Ifit1): Domain 16",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40")
    )
  ) &
  
theme(legend.position = "bottom")

save_plot(fig1, file.path(out_dir, "fig1_xy_ifn_immune_overlay"), w = 9.2, h = 8.0)

# ================================================================
# FIGURE 2 — % Ifng+ cells by BANKSY domain and condition
# Source: 04_banksy_joint_lam08_after_bloc3.rds (déjà chargé)
# ================================================================
message("\n=== Figure 2: % Ifng+ / proxy genes by BANKSY domain ===")

# se est déjà chargé depuis Figure 1 — on le réutilise
genes_panel <- rownames(se)

# Récupérer les domaines BANKSY (lam=0.2, res=0.9) comme "types cellulaires"
cd_fig2       <- as.data.frame(SummarizedExperiment::colData(se))
domain_labels <- paste0("Domain_", as.character(cd_fig2[[cluster_col]]))
samples_fig2  <- as.character(cd_fig2$sample)

# Restreindre à la niche (Domain_13 = Immune Acod1, Domain_16 = IFN responsive)
niche_mask    <- domain_labels %in% c(DOMAIN_IMMUNE, DOMAIN_IFN)
domain_labels <- domain_labels[niche_mask]
samples_fig2  <- samples_fig2[niche_mask]
message("  Cellules dans la niche : ", sum(niche_mask), " / ", length(niche_mask))

# Mapper les domaines vers leurs annotations
anno_data <- read.delim("ncells_by_sample_lam02_res09_joint_long.csv",
                        sep = ";", stringsAsFactors = FALSE)
anno_map <- anno_data %>%
  filter(annotation != "" & !is.na(annotation)) %>%
  mutate(annotation = trimws(annotation)) %>%
  select(banksy_domain, annotation) %>%
  distinct()
anno_lookup <- setNames(anno_map$annotation, anno_map$banksy_domain)
# Labels propres : annotation seule (sans préfixe Domain_XX)
domain_annot <- ifelse(
  !is.na(anno_lookup[domain_labels]) & anno_lookup[domain_labels] != "",
  anno_lookup[domain_labels],
  domain_labels
)

# Gènes à tester
genes_to_test <- IFN_PROXIMAL[IFN_PROXIMAL %in% genes_panel]
genes_absent  <- IFN_PROXIMAL[!IFN_PROXIMAL %in% genes_panel]
message("  Gènes trouvés dans le panel : ", paste(genes_to_test, collapse = ", "))
if (length(genes_absent) > 0)
  message("  Absents du panel           : ", paste(genes_absent, collapse = ", "))
if (length(genes_to_test) == 0) stop("Aucun gène proxy dans le panel.")

# Calculer % expressing par domaine × échantillon
counts_mat <- SummarizedExperiment::assay(se, "counts")

pct_list <- lapply(genes_to_test, function(g) {
  expr_vec <- as.numeric(counts_mat[g, niche_mask])
  tibble(
    domain    = domain_labels,
    annot     = domain_annot,
    sample    = samples_fig2,
    gene      = g,
    expressed = expr_vec > 0
  ) %>%
    group_by(gene, domain, annot, sample) %>%
    summarise(
      n_cells  = n(),
      n_pos    = sum(expressed),
      pct_expr = 100 * n_pos / n_cells,
      .groups  = "drop"
    )
})
pct_df <- bind_rows(pct_list)

# Forcer l'ordre des gènes (même ordre que IFN_PROXIMAL)
pct_df$gene <- factor(pct_df$gene, levels = genes_to_test)

pct_df$sample <- factor(
  pct_df$sample,
  levels = c(SAMPLE_ORDER[SAMPLE_ORDER %in% pct_df$sample],
             setdiff(unique(pct_df$sample), SAMPLE_ORDER))
)

sample_colors <- c(
  mock_6wpi  = "#999999",
  LCMV_1wpi  = "#56B4E9",
  LCMV_3wpi  = "#E69F00",
  LCMV_6wpi  = "#D55E00"
)

# Ordre des conditions dans la légende
sample_labels <- c(
  mock_6wpi  = "Mock 6 wpi",
  LCMV_1wpi  = "LCMV 1 wpi",
  LCMV_3wpi  = "LCMV 3 wpi",
  LCMV_6wpi  = "LCMV 6 wpi"
)

fig2 <- ggplot(pct_df,
       aes(x = annot, y = pct_expr, fill = sample)) +
  geom_bar(
    stat     = "identity",
    position = position_dodge(width = 0.75),
    width    = 0.65,
    color    = "white",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = sample_colors,
    labels = sample_labels,
    name   = NULL
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08)),
    labels = function(x) paste0(x, "%")
  ) +
  facet_wrap(~ gene, ncol = 2, scales = "fixed") +
  coord_flip() +
  labs(
    title    = "Expression of IFN-\u03b3 pathway genes within the immune niche",
    subtitle = "Proportion of cells with counts > 0 | Immune niche: Acod1+ immune cells (Domain 13) and IFN-responsive cells (Domain 16)",
    x = NULL,
    y = "% expressing cells"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle    = element_text(size = 9, hjust = 0, color = "grey50"),
    axis.text.y      = element_text(size = 10),
    axis.text.x      = element_text(size = 9),
    axis.line.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold.italic", size = 11),
    panel.spacing.x  = unit(1.8, "cm"),
    legend.position  = "top",
    legend.key.size  = unit(0.45, "cm"),
    legend.text      = element_text(size = 9),
    plot.margin      = margin(8, 16, 8, 8)
  )

save_plot(fig2, file.path(out_dir, "fig2_ifng_proxy_by_banksy_domain"),
          w = 10, h = 9)
write.csv(pct_df, file.path(out_dir, "fig2_ifng_proxy_by_banksy_domain.csv"),
          row.names = FALSE)

message("\n=== Done. Outputs in: ", out_dir, " ===\n")
