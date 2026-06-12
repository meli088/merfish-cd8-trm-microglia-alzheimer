#!/usr/bin/env Rscript
# =============================================================
# Script: 25_dam_signature_curation.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-05
#
# Goal:
#   Read the microglia/DAM signature Excel file compiled by Zoe
#   (Martin-Ferrera 2025), retain only the biologically meaningful
#   curated signature groups, and export clean gene lists for
#   downstream dotplots.
#
# Input:
#   outputs/Microglia.sig. Martin-ferrera2025.xlsx
#     — 3 columns: Gene | Population | Study
#     — 16 344 rows across all groups
#
# Retained signature groups:
#   - Upregulated DAM      (Keren-Shaul 2017)
#   - Downregulated DAM    (Keren-Shaul 2017)
#   - Microglia signature  (Gosselin 2017 + Galatro → per study + union)
#
# Excluded:
#   - Anonymous cluster labels: Mic0-3, Micro0-1, MG0-MG12, a-f, 0-12
#     (study-specific clusters without cross-study interpretable labels)
#   - Human DIMs / DAMs / YAMs (Silvin) — analysis is in mouse
#   - Aging (Galatro) — not the focus of this analysis
#
# Output folder: outputs/banksy/dam_signature_curation/
#
# Outputs:
#   raw_population_study_summary.csv
#   curated_summary.csv
#   Upregulated_DAM.csv / .txt
#   Downregulated_DAM.csv / .txt
#   Microglia_signature_Gosselin2017.csv / .txt
#   Microglia_signature_Galatro.csv / .txt
#   Microglia_signature_union.csv / .txt
#   wide_format_readable.csv
#   curation_report.txt
# =============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
})

base_path <- normalizePath(".")  # Run this script from the project root directory
setwd(base_path)

out_dir <- file.path("outputs", "banksy", "dam_signature_curation")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE,
                                      showWarnings = FALSE)

# =============================================================
# 1. Read the Excel file
# =============================================================

excel_path <- file.path("outputs", "Microglia.sig. Martin-ferrera2025.xlsx")
stopifnot("Excel file not found" = file.exists(excel_path))

message("Reading: ", excel_path)
raw <- read_excel(excel_path, sheet = 1) %>%
  rename_with(trimws) %>%
  mutate(
    Gene       = trimws(as.character(Gene)),
    Population = trimws(as.character(Population)),
    Study      = trimws(as.character(Study))
  ) %>%
  filter(!is.na(Gene), Gene != "")

message("  Total rows: ", nrow(raw))
message("  Unique Population x Study combinations: ",
        nrow(distinct(raw, Population, Study)))

# =============================================================
# 2. Raw summary table (all groups before curation)
# =============================================================

raw_summary <- raw %>%
  group_by(Population, Study) %>%
  summarise(n_genes_raw = n(), .groups = "drop") %>%
  arrange(Study, Population)

cat("\nAll Population × Study groups in the file:\n")
print(as.data.frame(raw_summary), row.names = FALSE)

write.csv(raw_summary,
          file.path(out_dir, "raw_population_study_summary.csv"),
          row.names = FALSE)
message("\nSaved: raw_population_study_summary.csv")

# =============================================================
# 3. Curation — define retained groups
# =============================================================

# Retain only these three Population labels:
RETAIN_POPULATIONS <- c(
  "Upregulated DAM",
  "Downregulated DAM",
  "Microglia signature"
)

curated <- raw %>%
  filter(Population %in% RETAIN_POPULATIONS)

message("\nRows retained after curation: ", nrow(curated))
cat("Retained groups:\n")
print(as.data.frame(count(curated, Population, Study)), row.names = FALSE)

# =============================================================
# 4. Gene name normalisation
# =============================================================
# Gene names in the file are a mix of human (ALL CAPS) and mouse
# (Title case, e.g. Trem2). We keep them as-is but remove exact
# duplicates within each signature group (case-sensitive).
# The downstream filtering to MERFISH panel will handle
# species-level harmonisation.

# =============================================================
# 5. Upregulated DAM
# =============================================================

up_dam <- curated %>%
  filter(Population == "Upregulated DAM") %>%
  pull(Gene) %>%
  unique() %>%
  sort()

message("\nUpregulated DAM: ", length(up_dam), " unique genes")

write.csv(data.frame(gene = up_dam),
          file.path(out_dir, "Upregulated_DAM.csv"),
          row.names = FALSE)
writeLines(up_dam, file.path(out_dir, "Upregulated_DAM.txt"))
message("  Saved: Upregulated_DAM.csv / .txt")

# =============================================================
# 6. Downregulated DAM
# =============================================================

down_dam <- curated %>%
  filter(Population == "Downregulated DAM") %>%
  pull(Gene) %>%
  unique() %>%
  sort()

message("\nDownregulated DAM: ", length(down_dam), " unique genes")

write.csv(data.frame(gene = down_dam),
          file.path(out_dir, "Downregulated_DAM.csv"),
          row.names = FALSE)
writeLines(down_dam, file.path(out_dir, "Downregulated_DAM.txt"))
message("  Saved: Downregulated_DAM.csv / .txt")

# =============================================================
# 7. Microglia signature — per study + union
# =============================================================

mg_sig <- curated %>% filter(Population == "Microglia signature")

mg_studies <- sort(unique(mg_sig$Study))
message("\nMicroglia signature studies found: ",
        paste(mg_studies, collapse = ", "))

mg_per_study <- list()

for (st in mg_studies) {
  genes_st <- mg_sig %>%
    filter(Study == st) %>%
    pull(Gene) %>%
    unique() %>%
    sort()

  # Build a safe filename (remove spaces and slashes)
  st_safe <- gsub("[^A-Za-z0-9_]", "", gsub(" ", "_", st))

  fname <- paste0("Microglia_signature_", st_safe)
  write.csv(data.frame(gene = genes_st),
            file.path(out_dir, paste0(fname, ".csv")),
            row.names = FALSE)
  writeLines(genes_st, file.path(out_dir, paste0(fname, ".txt")))
  message("  ", st, ": ", length(genes_st), " unique genes → ",
          fname, ".csv/.txt")

  mg_per_study[[st]] <- genes_st
}

# Union across all Microglia signature studies
mg_union <- sort(unique(unlist(mg_per_study)))
message("  Union (all studies): ", length(mg_union), " unique genes")

write.csv(data.frame(gene = mg_union),
          file.path(out_dir, "Microglia_signature_union.csv"),
          row.names = FALSE)
writeLines(mg_union, file.path(out_dir, "Microglia_signature_union.txt"))
message("  Saved: Microglia_signature_union.csv / .txt")

# =============================================================
# 8. Curated summary table
# =============================================================

# Per-study rows for Microglia signature
mg_study_rows <- lapply(names(mg_per_study), function(st) {
  data.frame(
    signature_name = "Microglia signature",
    study          = st,
    n_raw          = nrow(mg_sig %>% filter(Study == st)),
    n_unique       = length(mg_per_study[[st]]),
    notes          = "per-study list",
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

# Union row for Microglia signature
mg_union_row <- data.frame(
  signature_name = "Microglia signature",
  study          = paste(mg_studies, collapse = " + "),
  n_raw          = nrow(mg_sig),
  n_unique       = length(mg_union),
  notes          = "union across studies",
  stringsAsFactors = FALSE
)

# DAM rows
dam_rows <- data.frame(
  signature_name = c("Upregulated DAM", "Downregulated DAM"),
  study = c(
    unique(curated$Study[curated$Population == "Upregulated DAM"]),
    unique(curated$Study[curated$Population == "Downregulated DAM"])
  ),
  n_raw = c(
    nrow(curated %>% filter(Population == "Upregulated DAM")),
    nrow(curated %>% filter(Population == "Downregulated DAM"))
  ),
  n_unique = c(length(up_dam), length(down_dam)),
  notes    = "single study",
  stringsAsFactors = FALSE
)

curated_summary <- bind_rows(dam_rows, mg_study_rows, mg_union_row)

cat("\nCurated summary:\n")
print(as.data.frame(curated_summary), row.names = FALSE)

write.csv(curated_summary,
          file.path(out_dir, "curated_summary.csv"),
          row.names = FALSE)
message("\nSaved: curated_summary.csv")

# =============================================================
# 9. Wide-format readable CSV
# =============================================================

wide_lists <- list(
  Upregulated_DAM              = up_dam,
  Downregulated_DAM            = down_dam,
  Microglia_signature_union    = mg_union
)

# Add per-study Microglia signature columns
for (st in names(mg_per_study)) {
  st_safe <- gsub("[^A-Za-z0-9_]", "", gsub(" ", "_", st))
  wide_lists[[paste0("Microglia_signature_", st_safe)]] <- mg_per_study[[st]]
}

max_len <- max(sapply(wide_lists, length))
wide_df <- as.data.frame(lapply(wide_lists, function(g) {
  c(g, rep(NA_character_, max_len - length(g)))
}), check.names = FALSE)

write.csv(wide_df,
          file.path(out_dir, "wide_format_readable.csv"),
          row.names = FALSE, na = "")
message("Saved: wide_format_readable.csv")

# =============================================================
# 10. Curation report
# =============================================================

excluded_pop <- sort(setdiff(unique(raw$Population), RETAIN_POPULATIONS))

report_lines <- c(
  "=============================================================",
  "CURATION REPORT: DAM / Microglia signature gene lists",
  paste0("Script : 25_dam_signature_curation.R"),
  paste0("Date   : ", Sys.Date()),
  paste0("Input  : outputs/Microglia.sig. Martin-ferrera2025.xlsx"),
  paste0("Output : outputs/banksy/dam_signature_curation/"),
  "=============================================================",
  "",
  "--- Source file overview ---",
  paste0("Total rows : ", nrow(raw)),
  paste0("Unique Population x Study combinations: ",
         nrow(distinct(raw, Population, Study))),
  "",
  "--- Retained signature groups ---",
  paste0("  - Upregulated DAM      (Keren-Shaul 2017) : ",
         length(up_dam), " unique genes"),
  paste0("  - Downregulated DAM    (Keren-Shaul 2017) : ",
         length(down_dam), " unique genes"),
  paste0("  - Microglia signature  (per study):"),
  paste0(sapply(names(mg_per_study), function(st)
    paste0("      ", st, " : ", length(mg_per_study[[st]]), " unique genes"))),
  paste0("  - Microglia signature  (union)     : ",
         length(mg_union), " unique genes"),
  "",
  "--- Excluded population labels ---",
  paste0("  ", excluded_pop),
  "",
  "--- Rationale for exclusions ---",
  "  * Anonymous cluster labels (Mic0, MG0, Micro0, 0, 1, a, b, …):",
  "    Study-specific cluster IDs without cross-study interpretable names.",
  "    Including these would inflate gene lists with cluster-level markers",
  "    that do not correspond to a defined, named biological population.",
  "",
  "  * Human DIMs / DAMs / YAMs (Silvin):",
  "    This analysis is performed on mouse MERFISH data.",
  "    Human signatures are intentionally excluded to avoid cross-species",
  "    gene name mismatches (human ALLCAPS vs mouse Title case).",
  "",
  "  * Aging (Galatro):",
  "    Not the focus of the current analysis.",
  "    Can be re-added if aging comparisons are required later.",
  "",
  "--- Downstream usage notes ---",
  "  Full CSVs / TXTs contain all unique curated genes.",
  "  Downstream dotplot figures will additionally:",
  "    - filter to the MERFISH gene panel",
  "    - select top 10–15 genes per signature",
  "  These filtering steps are NOT performed here.",
  "",
  "--- Output files ---",
  "  raw_population_study_summary.csv  — all groups before curation",
  "  curated_summary.csv               — retained groups with n_unique",
  "  Upregulated_DAM.csv / .txt",
  "  Downregulated_DAM.csv / .txt",
  paste0(sapply(names(mg_per_study), function(st) {
    st_safe <- gsub("[^A-Za-z0-9_]", "", gsub(" ", "_", st))
    paste0("  Microglia_signature_", st_safe, ".csv / .txt")
  })),
  "  Microglia_signature_union.csv / .txt",
  "  wide_format_readable.csv",
  "  curation_report.txt",
  "============================================================="
)

writeLines(report_lines, file.path(out_dir, "curation_report.txt"))
message("Saved: curation_report.txt")

# =============================================================
# Done
# =============================================================

message("\n=== Done. Outputs in: ", out_dir, " ===\n")
