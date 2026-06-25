#!/usr/bin/env Rscript
# =============================================================
# Script: 45_rfp_fc_all_clusters.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-06
#
# Goal:
#   Pour CHAQUE annotation de l'objet global (tous les clusters lam02 res09),
#   calculer le log2FC d'intensité RFP (raw + high_pass) :
#     cluster focal  vs  toutes autres cellules annotées (intra-échantillon)
#
#   Deux mesures retenues (vol_norm supprimé) :
#     raw      : Anti.RFP_raw
#     highpass : Anti.RFP_high_pass
#
# Input :
#   objects/03_rfp_analysis.rds
#     — Seurat avec colonnes protéiques corrigées (EntityID join) et
#       CellType issu de 05_joint_annotated_lam02_res09.rds
#       → construit par 34_rebuild_protein_metadata.R
#
# Outputs :  outputs/banksy/rfp_plots_readable/
#   rfp_fc_all_clusters.csv
#   fig_rfp_fc_all_clusters_barplot.pdf/jpg  — facet par annotation
#   fig_rfp_fc_all_clusters_heatmap.pdf/jpg  — heatmap annotations × conditions
# =============================================================

set.seed(1997)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(Cairo)
})

setwd(normalizePath("."))
source("scripts/00_palette.R")

# =============================================================
# 1. Paramètres
# =============================================================

OBJ_IN       <- file.path("objects", "03_rfp_analysis.rds")
OUT_DIR      <- file.path("outputs", "banksy", "rfp_plots_readable")

SAMPLE_ORDER  <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi", "mock_6wpi")
SAMPLE_LABELS <- c(
  mock_6wpi = "Mock 6 wpi",
  LCMV_1wpi = "LCMV 1 wpi",
  LCMV_3wpi = "LCMV 3 wpi",
  LCMV_6wpi = "LCMV 6 wpi"
)

UNANNOTATED <- "Non annote"
MIN_CELLS   <- 5

MEASURES <- c(raw = "Anti.RFP_raw", highpass = "Anti.RFP_high_pass")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================
# 2. Helper : sauvegarde cairo_pdf + JPEG
# =============================================================

save_fig <- function(p, stem, width, height) {
  cairo_pdf(paste0(stem, ".pdf"), width = width, height = height)
  print(p)
  dev.off()
  CairoJPEG(paste0(stem, ".jpg"),
            width = width * 150, height = height * 150, res = 150)
  print(p)
  dev.off()
  message("  Saved: ", basename(stem), ".pdf / .jpg")
}

# =============================================================
# 3. Chargement de l'objet
# =============================================================

message("\n=== Chargement de l'objet ===")
obj <- readRDS(OBJ_IN)
message("  Cellules  : ", ncol(obj))

required_cols <- c("Anti.RFP_raw", "Anti.RFP_high_pass", "CellType", "sample")
missing_cols  <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_cols) > 0) {
  stop("Colonnes manquantes : ", paste(missing_cols, collapse = ", "),
       "\n  -> Lancer 34_rebuild_protein_metadata.R d'abord.")
}

obj$sample <- factor(obj$sample,
                     levels = intersect(SAMPLE_ORDER, unique(as.character(obj$sample))))

ana_df <- as.data.frame(obj@meta.data) %>%
  filter(CellType != UNANNOTATED, !is.na(CellType)) %>%
  mutate(sample      = factor(sample, levels = levels(obj$sample)),
         condition   = factor(SAMPLE_LABELS[as.character(sample)],
                              levels = SAMPLE_LABELS))

message("  Cellules annotées : ", nrow(ana_df))

all_annotations <- sort(unique(as.character(ana_df$CellType)))
message("  Annotations : ", length(all_annotations))
for (a in all_annotations) message("    ", a)

# =============================================================
# 4. Calcul du log2FC (médiane) par annotation × sample × mesure
# =============================================================

message("\n=== Calcul log2FC intensité RFP ===")

fc_rows <- list()

for (samp in levels(ana_df$sample)) {
  df_s <- ana_df %>% filter(sample == samp)

  for (anno in all_annotations) {
    cells_in  <- df_s %>% filter(CellType == anno)
    cells_out <- df_s %>% filter(CellType != anno)

    n_in  <- nrow(cells_in)
    n_out <- nrow(cells_out)

    if (n_in < MIN_CELLS) {
      message("  SKIP (n_in=", n_in, ") : ", anno, " / ", samp)
      next
    }

    for (mname in names(MEASURES)) {
      col   <- MEASURES[[mname]]
      v_in  <- cells_in[[col]][!is.na(cells_in[[col]])]
      v_out <- cells_out[[col]][!is.na(cells_out[[col]])]

      med_in  <- median(v_in)
      med_out <- median(v_out)

      log2fc <- if (!is.na(med_out) && med_out > 0) {
        log2(med_in / med_out)
      } else {
        NA_real_
      }

      fc_rows[[length(fc_rows) + 1]] <- data.frame(
        sample      = samp,
        condition   = SAMPLE_LABELS[[samp]],
        annotation  = anno,
        measure     = mname,
        median_in   = med_in,
        median_out  = med_out,
        log2FC      = log2fc,
        n_in        = n_in,
        n_out       = n_out,
        stringsAsFactors = FALSE
      )
    }
  }
}

fc_df <- bind_rows(fc_rows) %>%
  mutate(
    sample     = factor(sample,     levels = levels(ana_df$sample)),
    condition  = factor(condition,  levels = SAMPLE_LABELS),
    annotation = factor(annotation, levels = intersect(ANNOTATION_ORDER,
                                                        all_annotations)),
    measure    = factor(measure,    levels = names(MEASURES))
  )

message("  Lignes : ", nrow(fc_df))
print(fc_df %>% select(sample, annotation, measure, log2FC, n_in) %>%
      arrange(annotation, sample, measure))

csv_out <- file.path(OUT_DIR, "rfp_fc_all_clusters.csv")
write.csv(fc_df, csv_out, row.names = FALSE)
message("  Saved: rfp_fc_all_clusters.csv")

# =============================================================
# 5. Figure A — Barplot facetté par annotation
#    x = condition, y = log2FC, fill = mesure
# =============================================================

message("\n=== Figure A : barplot facetté par annotation ===")

MEASURE_COLORS <- c(raw = "#4393C3", highpass = "#D6604D")
MEASURE_LABELS <- c(raw = "Anti.RFP_raw", highpass = "Anti.RFP_high_pass")

n_anno <- length(levels(fc_df$annotation))
ncol_facet <- min(4L, n_anno)

p_bar <- ggplot(fc_df,
                aes(x = condition, y = log2FC,
                    fill = measure, group = measure)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.65, color = "white", linewidth = 0.2,
           na.rm = TRUE) +
  geom_hline(yintercept = 0, linewidth = 0.45, color = "grey30") +
  facet_wrap(~ annotation, ncol = ncol_facet,
             labeller = label_wrap_gen(width = 26)) +
  scale_fill_manual(values = MEASURE_COLORS,
                    labels = MEASURE_LABELS,
                    name = "Mesure RFP") +
  scale_x_discrete(labels = function(x) sub(" wpi", "\nwpi", x)) +
  labs(
    title    = "log2FC intensité RFP — cluster focal vs toutes autres cellules annotées",
    subtitle = paste0("log2(médiane_in / médiane_out) | intra-échantillon | ",
                      "hors cluster = toutes cellules annotées sauf cluster focal"),
    x = NULL,
    y = "log2FC (médiane in / out)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    plot.title       = element_text(face = "bold", size = 11),
    plot.subtitle    = element_text(size = 8, color = "grey40"),
    axis.text.x      = element_text(angle = 0, hjust = 0.5, size = 8),
    axis.text.y      = element_text(size = 8),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold", size = 8),
    panel.spacing    = unit(0.5, "cm"),
    legend.position  = "bottom"
  )

nrow_facet <- ceiling(n_anno / ncol_facet)
save_fig(p_bar,
         file.path(OUT_DIR, "fig_rfp_fc_all_clusters_barplot"),
         width = 14, height = 3.5 * nrow_facet + 1.5)

# =============================================================
# 6. Figure B — Heatmap annotations × conditions (facet par mesure)
# =============================================================

message("\n=== Figure B : heatmap annotations × conditions ===")

# Symétrie de l'échelle autour de 0
max_abs <- max(abs(fc_df$log2FC), na.rm = TRUE)
max_abs <- ceiling(max_abs * 10) / 10   # arrondi au 0.1 supérieur

p_heat <- ggplot(fc_df,
                 aes(x = condition, y = annotation, fill = log2FC)) +
  geom_tile(color = "white", linewidth = 0.35) +
  facet_wrap(~ measure,
             labeller = as_labeller(MEASURE_LABELS)) +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    limits   = c(-max_abs, max_abs),
    na.value = "grey88",
    name     = "log2FC\n(médiane)"
  ) +
  scale_x_discrete(labels = function(x) sub(" wpi", "\nwpi", x)) +
  scale_y_discrete(limits = rev) +
  labs(
    title    = "log2FC intensité RFP par cluster et condition",
    subtitle = "log2(médiane_in / médiane_out) | intra-échantillon | NA = effectif < 5",
    x = NULL, y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    plot.title       = element_text(face = "bold", size = 11),
    plot.subtitle    = element_text(size = 8, color = "grey40"),
    axis.text.x      = element_text(size = 9),
    axis.text.y      = element_text(size = 9),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold", size = 10),
    panel.spacing    = unit(1, "cm"),
    legend.position  = "right"
  )

save_fig(p_heat,
         file.path(OUT_DIR, "fig_rfp_fc_all_clusters_heatmap"),
         width = 12, height = max(6, 0.38 * n_anno + 2))

# =============================================================
# 7. Résumé console
# =============================================================

# =============================================================
# 7b. Figure C — Immune + IFN seuls (2 panneaux)
# =============================================================

message("\n=== Figure C : Immune + IFN seuls ===")

IMMUNE_IFN <- c("Immune (Acod1)", "IFN responsive (Ifit1)")

fc_immune_ifn <- fc_df %>%
  filter(annotation %in% IMMUNE_IFN) %>%
  mutate(annotation = factor(annotation, levels = IMMUNE_IFN))

make_barplot <- function(df, title_str, ncol_f) {
  ggplot(df, aes(x = condition, y = log2FC, fill = measure, group = measure)) +
    geom_col(position = position_dodge(width = 0.7),
             width = 0.65, color = "white", linewidth = 0.2,
             na.rm = TRUE) +
    geom_hline(yintercept = 0, linewidth = 0.45, color = "grey30") +
    facet_wrap(~ annotation, ncol = ncol_f,
               labeller = label_wrap_gen(width = 28)) +
    scale_fill_manual(values = MEASURE_COLORS, labels = MEASURE_LABELS,
                      name = "Mesure RFP") +
    scale_x_discrete(labels = function(x) sub(" wpi", "\nwpi", x)) +
    labs(
      title    = title_str,
      subtitle = "log2(médiane_in / médiane_out) | intra-échantillon | out = toutes cellules annotées hors groupe",
      x = NULL, y = "log2FC (médiane in / out)"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(size = 8.5, color = "grey40"),
      axis.text.x      = element_text(angle = 0, hjust = 0.5, size = 9),
      axis.text.y      = element_text(size = 9),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(face = "bold", size = 10),
      panel.spacing    = unit(0.8, "cm"),
      legend.position  = "bottom"
    )
}

p_immune_ifn <- make_barplot(
  fc_immune_ifn,
  "log2FC intensité RFP — Immune (Acod1) et IFN responsive (Ifit1)",
  ncol_f = 2
)

save_fig(p_immune_ifn,
         file.path(OUT_DIR, "fig_rfp_fc_immune_ifn"),
         width = 8, height = 5.5)

# =============================================================
# 7c. Figure D — groupes biologiques poolés
#     Immune | IFN | Neurons (poolés) | Microglia (poolés)
#     Astrocytes (poolés) | Oligodendrocytes
# =============================================================

message("\n=== Figure D : groupes biologiques poolés ===")

# Définition des groupes d'appartenance
NEURON_TYPES <- grep("^Neurons?\\s*\\(", all_annotations, value = TRUE)
MICRO_TYPES  <- grep("^Microglia\\s*\\(", all_annotations, value = TRUE)
ASTRO_TYPES  <- grep("^Astrocytes?\\s*\\(", all_annotations, value = TRUE)
OLIGO_TYPES  <- grep("^Oligodendrocytes?\\s*\\(", all_annotations, value = TRUE)

message("  Neurons  : ", paste(NEURON_TYPES, collapse = ", "))
message("  Microglia: ", paste(MICRO_TYPES,  collapse = ", "))
message("  Astro    : ", paste(ASTRO_TYPES,  collapse = ", "))
message("  Oligo    : ", paste(OLIGO_TYPES,  collapse = ", "))

# Tous les sous-clusters individuels (aucun pooling)
indiv_neuron <- setNames(as.list(NEURON_TYPES), NEURON_TYPES)
indiv_micro  <- setNames(as.list(MICRO_TYPES),  MICRO_TYPES)
indiv_astro  <- setNames(as.list(ASTRO_TYPES),  ASTRO_TYPES)
indiv_oligo  <- setNames(as.list(OLIGO_TYPES),  OLIGO_TYPES)

GROUP_MAP <- c(
  list(
    "Immune (Acod1)"         = c("Immune (Acod1)"),
    "IFN responsive (Ifit1)" = c("IFN responsive (Ifit1)")
  ),
  indiv_neuron,
  indiv_micro,
  indiv_astro,
  indiv_oligo
)
GROUP_MAP <- Filter(function(x) length(x) > 0, GROUP_MAP)

GROUP_ORDER <- names(GROUP_MAP)

fc_group_rows <- list()

for (samp in levels(ana_df$sample)) {
  df_s <- ana_df %>% filter(sample == samp)

  for (grp_name in GROUP_ORDER) {
    in_types  <- GROUP_MAP[[grp_name]]
    cells_in  <- df_s %>% filter(CellType %in% in_types)
    cells_out <- df_s %>% filter(!CellType %in% in_types)

    n_in  <- nrow(cells_in)
    n_out <- nrow(cells_out)

    if (n_in < MIN_CELLS) {
      message("  SKIP (n_in=", n_in, ") : ", grp_name, " / ", samp)
      next
    }

    for (mname in names(MEASURES)) {
      col   <- MEASURES[[mname]]
      v_in  <- cells_in[[col]][!is.na(cells_in[[col]])]
      v_out <- cells_out[[col]][!is.na(cells_out[[col]])]

      med_in  <- median(v_in)
      med_out <- median(v_out)

      log2fc <- if (!is.na(med_out) && med_out > 0) log2(med_in / med_out) else NA_real_

      fc_group_rows[[length(fc_group_rows) + 1]] <- data.frame(
        sample     = samp,
        condition  = SAMPLE_LABELS[[samp]],
        annotation = grp_name,
        measure    = mname,
        median_in  = med_in,
        median_out = med_out,
        log2FC     = log2fc,
        n_in       = n_in,
        n_out      = n_out,
        stringsAsFactors = FALSE
      )
    }
  }
}

fc_group_df <- bind_rows(fc_group_rows) %>%
  mutate(
    sample     = factor(sample,     levels = levels(ana_df$sample)),
    condition  = factor(condition,  levels = SAMPLE_LABELS),
    annotation = factor(annotation, levels = GROUP_ORDER),
    measure    = factor(measure,    levels = names(MEASURES))
  )

csv_group <- file.path(OUT_DIR, "rfp_fc_grouped_clusters.csv")
write.csv(fc_group_df, csv_group, row.names = FALSE)
message("  Saved: rfp_fc_grouped_clusters.csv")

p_grouped <- make_barplot(
  fc_group_df,
  "log2FC intensité RFP — Immune, IFN, Neurons, Microglia, Astrocytes, Oligodendrocytes",
  ncol_f = 3
)

n_grp <- length(GROUP_ORDER)
save_fig(p_grouped,
         file.path(OUT_DIR, "fig_rfp_fc_grouped_clusters"),
         width = 13, height = 3.5 * ceiling(n_grp / 3) + 1.5)

# =============================================================
# 7d. Figure E — 5 groupes biologiques poolés (résumé compact)
#     Immune | IFN | Neurons | Microglia | Oligodendrocytes
# =============================================================

message("\n=== Figure E : 5 groupes poolés (résumé compact) ===")

SUMMARY_MAP <- list(
  "Immune (Acod1)"         = c("Immune (Acod1)"),
  "IFN responsive (Ifit1)" = c("IFN responsive (Ifit1)"),
  "Neurons"                = NEURON_TYPES,
  "Microglia"              = MICRO_TYPES,
  "Oligodendrocytes"       = OLIGO_TYPES
)
SUMMARY_MAP   <- Filter(function(x) length(x) > 0, SUMMARY_MAP)
SUMMARY_ORDER <- names(SUMMARY_MAP)

fc_summary_rows <- list()

for (samp in levels(ana_df$sample)) {
  df_s <- ana_df %>% filter(sample == samp)

  for (grp_name in SUMMARY_ORDER) {
    in_types  <- SUMMARY_MAP[[grp_name]]
    cells_in  <- df_s %>% filter(CellType %in% in_types)
    cells_out <- df_s %>% filter(!CellType %in% in_types)

    n_in  <- nrow(cells_in)
    n_out <- nrow(cells_out)

    if (n_in < MIN_CELLS) next

    for (mname in names(MEASURES)) {
      col    <- MEASURES[[mname]]
      v_in   <- cells_in[[col]][!is.na(cells_in[[col]])]
      v_out  <- cells_out[[col]][!is.na(cells_out[[col]])]
      med_in  <- median(v_in)
      med_out <- median(v_out)
      log2fc  <- if (!is.na(med_out) && med_out > 0) log2(med_in / med_out) else NA_real_

      fc_summary_rows[[length(fc_summary_rows) + 1]] <- data.frame(
        sample     = samp,
        condition  = SAMPLE_LABELS[[samp]],
        annotation = grp_name,
        measure    = mname,
        median_in  = med_in,
        median_out = med_out,
        log2FC     = log2fc,
        n_in       = n_in,
        n_out      = n_out,
        stringsAsFactors = FALSE
      )
    }
  }
}

fc_summary_df <- bind_rows(fc_summary_rows) %>%
  mutate(
    sample     = factor(sample,     levels = levels(ana_df$sample)),
    condition  = factor(condition,  levels = SAMPLE_LABELS),
    annotation = factor(annotation, levels = SUMMARY_ORDER),
    measure    = factor(measure,    levels = names(MEASURES))
  )

p_summary <- make_barplot(
  fc_summary_df,
  "log2FC intensité RFP — Immune, IFN, Neurons, Microglia, Oligodendrocytes",
  ncol_f = 5
)

save_fig(p_summary,
         file.path(OUT_DIR, "fig_rfp_fc_summary_5groups"),
         width = 14, height = 5)

# =============================================================
# 8. Résumé console
# =============================================================

message("\n", strrep("=", 60))
message("SUMMARY — 45_rfp_fc_all_clusters")
message(strrep("=", 60))
message("  Cellules annotées    : ", nrow(ana_df))
message("  Annotations          : ", n_anno)
message("  Lignes résultat      : ", nrow(fc_df))
message("  Output folder        : ", normalizePath(OUT_DIR))
message("  rfp_fc_all_clusters.csv")
message("  fig_rfp_fc_all_clusters_barplot.pdf/jpg")
message("  fig_rfp_fc_all_clusters_heatmap.pdf/jpg")
message("  fig_rfp_fc_immune_ifn.pdf/jpg")
message("  rfp_fc_grouped_clusters.csv")
  message("  fig_rfp_fc_grouped_clusters.pdf/jpg")
  message("  fig_rfp_fc_summary_5groups.pdf/jpg")
message(strrep("=", 60))
