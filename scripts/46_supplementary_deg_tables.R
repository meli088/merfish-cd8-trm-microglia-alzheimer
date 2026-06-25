#!/usr/bin/env Rscript
# =============================================================
# Script: 46_supplementary_deg_tables.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-06
#
# Goal:
#   Compiler toutes les tables DEG du pipeline en un seul
#   fichier Excel multi-onglets (Supplementary_DEG_tables.xlsx).
#
#   Pour chaque onglet :
#     - Colonne "significant" ajoutée (p_val_adj < 0.05 & |log2FC| > 0.25)
#     - Trié par avg_log2FC décroissant
#     - En-têtes en gras, lignes signif. surlignées en jaune
#     - Largeurs de colonnes auto-ajustées
#     - Si CSV absent : onglet vide avec message "File not found"
#
# Output : outputs/banksy/supplementary/Supplementary_DEG_tables.xlsx
# =============================================================

set.seed(1997)

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
})

setwd(normalizePath("."))

OUT_DIR  <- file.path("outputs", "banksy", "supplementary")
OUT_FILE <- file.path(OUT_DIR, "Supplementary_DEG_tables.xlsx")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================
# 1. Définition des onglets (ordre, nom, chemin)
# =============================================================

TABS <- list(
  list(
    sheet = "DEG_global_clusters",
    path  = file.path("outputs", "banksy", "umap_annotated_OK",
                      "DEG_annotated_lam02_res09.csv")
  ),
  list(
    sheet = "DEG_immune_subclusters",
    path  = file.path("outputs", "banksy", "immune_acod1",
                      "analysis", "DEGs", "DEGs_all_res03.csv")
  ),
  list(
    sheet = "DEG_microglia_1wpi_vs_mock",
    path  = file.path("outputs", "banksy", "microglia_conditions_deg",
                      "DEG_microglia_lcmv_1wpi_vs_mock_6wpi.csv")
  ),
  list(
    sheet = "DEG_microglia_6wpi_vs_mock",
    path  = file.path("outputs", "banksy", "microglia_conditions_deg",
                      "DEG_microglia_lcmv_6wpi_vs_mock_6wpi.csv")
  ),
  list(
    sheet = "DEG_microglia_6wpi_vs_1wpi",
    path  = file.path("outputs", "banksy", "microglia_conditions_deg",
                      "DEG_microglia_lcmv_6wpi_vs_lcmv_1wpi.csv")
  ),
  list(
    sheet = "DEG_microglia_in_vs_out_1wpi",
    path  = file.path("outputs", "banksy", "microglia_dam_niche",
                      "DEG_microglia_in_vs_out_lcmv_1wpi.csv")
  ),
  list(
    sheet = "DEG_microglia_in_vs_out_3wpi",
    path  = file.path("outputs", "banksy", "microglia_dam_niche",
                      "DEG_microglia_in_vs_out_lcmv_3wpi.csv")
  ),
  list(
    sheet = "DEG_microglia_in_vs_out_6wpi",
    path  = file.path("outputs", "banksy", "microglia_dam_niche",
                      "DEG_microglia_in_vs_out_lcmv_6wpi.csv")
  ),
  list(
    sheet = "DEG_microglia_in_vs_out_merged",
    path  = file.path("outputs", "banksy", "microglia_dam_niche",
                      "DEG_microglia_in_vs_out_all_timepoints_merged.csv")
  ),
  list(
    sheet = "DEG_tcells_6wpi_vs_1wpi",
    path  = file.path("outputs", "banksy", "immune_niche_volcano_by_celltype",
                      "DEG_t_cells_gzmb__lcmv_6wpi_vs_lcmv_1wpi.csv")
  ),
  list(
    sheet = "DEG_mac_6wpi_vs_1wpi",
    path  = file.path("outputs", "banksy", "immune_niche_volcano_by_celltype",
                      "DEG_mac_ctss__lcmv_6wpi_vs_lcmv_1wpi.csv")
  ),
  list(
    sheet = "DEG_microglia_c1qa_6wpi_vs_1wpi",
    path  = file.path("outputs", "banksy", "immune_niche_volcano_by_celltype",
                      "DEG_microglia_c1qa__lcmv_6wpi_vs_lcmv_1wpi.csv")
  ),
  list(
    sheet = "DEG_ifng_pos_vs_neg",
    path  = file.path("outputs", "banksy", "ifn_immune_overlay",
                      "DEG_ifng_pos_vs_neg_merged.csv")
  )
)

# =============================================================
# 2. Seuils significativité
# =============================================================

FDR_CUTOFF <- 0.05
FC_CUTOFF  <- 0.25

# =============================================================
# 3. Styles openxlsx
# =============================================================

style_header <- createStyle(
  fontName  = "Arial",
  fontSize  = 10,
  textDecoration = "bold",
  fgFill    = "#D9D9D9",
  border    = "Bottom",
  borderColour = "#666666",
  wrapText  = FALSE
)

style_signif <- createStyle(
  fgFill = "#FFFACD"   # jaune citron clair
)

style_body <- createStyle(
  fontName = "Arial",
  fontSize = 10
)

# =============================================================
# 4. Création du workbook
# =============================================================

wb <- createWorkbook()

for (tab in TABS) {
  sname <- tab$sheet
  fpath <- tab$path

  message("\n--- Onglet : ", sname)
  addWorksheet(wb, sheetName = sname)

  # ── Fichier absent ──────────────────────────────────────────
  if (!file.exists(fpath)) {
    warning("  Fichier introuvable : ", fpath)
    msg_df <- data.frame(message = paste("File not found:", fpath),
                         stringsAsFactors = FALSE)
    writeData(wb, sheet = sname, x = msg_df, startRow = 1,
              headerStyle = style_header)
    addStyle(wb, sheet = sname, style = style_body,
             rows = 2, cols = 1, gridExpand = TRUE)
    setColWidths(wb, sheet = sname, cols = 1, widths = "auto")
    next
  }

  # ── Lecture ─────────────────────────────────────────────────
  df <- tryCatch(
    read.csv(fpath, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      warning("  Erreur lecture : ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(df) || nrow(df) == 0) {
    msg_df <- data.frame(message = paste("Empty or unreadable file:", fpath),
                         stringsAsFactors = FALSE)
    writeData(wb, sheet = sname, x = msg_df, startRow = 1,
              headerStyle = style_header)
    next
  }

  message("  Lignes lues : ", nrow(df), " | Colonnes : ",
          paste(names(df), collapse = ", "))

  # ── Colonne "significant" ────────────────────────────────────
  has_fdr <- "p_val_adj"     %in% names(df)
  has_fc  <- "avg_log2FC"    %in% names(df)

  if (has_fdr && has_fc) {
    df <- df %>%
      mutate(
        significant = !is.na(p_val_adj) & !is.na(avg_log2FC) &
                      p_val_adj < FDR_CUTOFF &
                      abs(avg_log2FC) > FC_CUTOFF
      ) %>%
      arrange(desc(avg_log2FC))
  } else {
    df$significant <- NA
    message("  WARN : colonnes p_val_adj / avg_log2FC absentes — significant = NA")
  }

  # ── Écriture ─────────────────────────────────────────────────
  writeData(wb, sheet = sname, x = df, startRow = 1,
            headerStyle = style_header, borders = "none")

  n_rows <- nrow(df)
  n_cols <- ncol(df)

  # Style corps
  if (n_rows > 0) {
    addStyle(wb, sheet = sname, style = style_body,
             rows = 2:(n_rows + 1), cols = 1:n_cols,
             gridExpand = TRUE)
  }

  # Surlignage lignes significatives
  if (has_fdr && has_fc && any(df$significant, na.rm = TRUE)) {
    sig_rows <- which(df$significant) + 1L   # +1 pour header
    addStyle(wb, sheet = sname, style = style_signif,
             rows = sig_rows, cols = 1:n_cols,
             gridExpand = TRUE, stack = TRUE)
    message("  Lignes surlignées : ", length(sig_rows))
  }

  # Largeurs auto
  setColWidths(wb, sheet = sname, cols = 1:n_cols, widths = "auto")

  message("  OK — ", n_rows, " lignes écrites")
}

# =============================================================
# 5. Sauvegarde
# =============================================================

saveWorkbook(wb, file = OUT_FILE, overwrite = TRUE)
message("\n", strrep("=", 60))
message("Fichier sauvegardé : ", normalizePath(OUT_FILE))
message("Onglets : ", length(TABS))
message(strrep("=", 60))
