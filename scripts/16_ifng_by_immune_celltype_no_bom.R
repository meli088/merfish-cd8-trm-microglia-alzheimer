#!/usr/bin/env Rscript
# =============================================================
# Script: 16_ifng_by_immune_celltype.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
# Author: Mélina Farshchi
# Date: 2026-05
#
# Biological question:
#   Do T cells, especially T CD8, account for most of the Ifng signal
#   in the immune niche during LCMV infection?
#
# Input:
#   objects/08_immune_annotated_lam02_res03.rds
#   BANKSY sub-clustering on Domain_13 (Immune Acod1), lambda=0.2, res=0.3
#   Manual annotations in colData(se)$cell_type (17 sub-domains)
#
# Strategy:
#   Ifng expression is extremely sparse (few positive cells, few distinct
#   expression levels). Line/polygon plots are not appropriate.
#   We use:
#     Fig 1 -- Bubble plot (global overview, all cell types)
#              x = condition, y = cell type
#              size = % Ifng+ cells, color = median Ifng expr among positives
#     Fig 2 -- Strip plot (intensity in Ifng+ cells, focus populations only)
#              restricted to cell types with >= MIN_POS_FOR_STRIP Ifng+ cells
#
# Outputs:
#   CSV:
#     ifng_per_cell_table.csv             -- one row per cell
#     ifng_summary_by_celltype.csv        -- per cell type, all conditions
#     ifng_summary_by_celltype_sample.csv -- per cell type x condition
#   Figures (PDF + JPG):
#     fig1_ifng_bubble                    -- bubble plot, all cell types
#     fig2_ifng_strip_positives           -- strip plot, focus populations
# =============================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(tidyverse)
  library(ggplot2)
})

base_path <- normalizePath(".")  # Run this script from the project root directory
setwd(base_path)

out_dir <- file.path("outputs", "banksy", "ifn_immune_overlay")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================
# Parameters
# =============================================================

SAMPLE_ORDER <- c("LCMV_1wpi", "LCMV_3wpi", "LCMV_6wpi")

sample_labels <- c(
  LCMV_1wpi = "LCMV\n1 wpi",
  LCMV_3wpi = "LCMV\n3 wpi",
  LCMV_6wpi = "LCMV\n6 wpi"
)

sample_colors <- c(
  LCMV_1wpi = "#56B4E9",
  LCMV_3wpi = "#E69F00",
  LCMV_6wpi = "#D55E00"
)

# Minimum number of Ifng+ cells (across all conditions) required to include
# a cell type in the strip intensity plot (Fig 2).
MIN_POS_FOR_STRIP <- 3L

# Patterns for focus populations (flexible matching against actual labels)
FOCUS_PATTERNS <- c("T CD8", "T CD4", "Mono", "Mac", "Microglia",
                    "IFN", "Ifit", "Interferon")

# Preferred biological order for y-axis labels
PREFERRED_ORDER <- c("T CD8", "T CD4", "Mono", "Mac", "Microglia")

# =============================================================
# 1. Load annotated object and inspect Ifng assay
# =============================================================

message("Loading: objects/08_immune_annotated_lam02_res03.rds")
se_imm <- readRDS(file.path("objects", "08_immune_annotated_lam02_res03.rds"))
# Exclude mock_6wpi — analysis focuses on LCMV time course
se_imm <- se_imm[, as.character(SummarizedExperiment::colData(se_imm)$sample) != "mock_6wpi"]
message("  ", ncol(se_imm), " cells loaded (mock_6wpi excluded)")

cd_imm     <- as.data.frame(colData(se_imm))
counts_mat <- SummarizedExperiment::assay(se_imm, "counts")

# Verify expression scale.
# If values are continuous non-integers (0.693, 1.099, 1.386 ...),
# the assay stores log1p-like normalized expression, not raw integer counts.
# The stored values are used as-is; Ifng_positive = stored value > 0.
ifng_raw <- as.numeric(counts_mat["Ifng", ])

cat("\n--- Ifng assay inspection ---\n")
cat("  n cells         :", length(ifng_raw), "\n")
cat("  Min             :", min(ifng_raw), "\n")
cat("  Max             :", round(max(ifng_raw), 4), "\n")
cat("  n > 0           :", sum(ifng_raw > 0), "\n")
pos_vals <- sort(unique(round(ifng_raw[ifng_raw > 0], 4)))
cat("  Unique non-zero :", paste(head(pos_vals, 12), collapse = ", "), "\n")
cat("  Scale           :",
    if (all(pos_vals == round(pos_vals))) "integer counts (raw)"
    else "log1p-like (continuous, not raw counts)", "\n\n")

ifng_expr <- round(ifng_raw, 4)  # 4 dp to avoid floating-point key mismatches

cat("Cells per cell type (all conditions combined):\n")
print(sort(table(cd_imm$cell_type), decreasing = TRUE))

# =============================================================
# 2. Per-cell table
# =============================================================

message("\nBuilding per-cell table...")

cell_df <- data.frame(
  cell_id       = colnames(se_imm),
  sample        = factor(as.character(cd_imm$sample), levels = SAMPLE_ORDER),
  cell_type     = as.character(cd_imm$cell_type),
  Ifng_expr     = ifng_expr,     # stored log1p-like expression value
  Ifng_positive = ifng_expr > 0, # TRUE if any Ifng expression detected
  stringsAsFactors = FALSE
)

cat("\nIfng+ cells per condition:\n")
print(
  cell_df %>%
    group_by(sample) %>%
    summarise(n_cells = n(), n_pos = sum(Ifng_positive),
              pct_pos = round(100 * n_pos / n_cells, 2), .groups = "drop")
)

write.csv(cell_df,
          file.path(out_dir, "ifng_per_cell_table.csv"),
          row.names = FALSE)
message("Saved: ifng_per_cell_table.csv")

# =============================================================
# 3. Summary tables
# =============================================================

message("\nComputing summary tables...")

## 3a. Per cell type x condition
summary_ct_samp <- cell_df %>%
  group_by(cell_type, sample) %>%
  summarise(
    n_cells         = n(),
    n_pos           = sum(Ifng_positive),
    pct_pos         = round(100 * n_pos / n_cells, 3),
    mean_expr_all   = round(mean(Ifng_expr), 5),
    median_expr_all = round(median(Ifng_expr), 5),
    max_expr        = round(max(Ifng_expr), 5),
    mean_expr_pos   = ifelse(n_pos > 0,
                             round(mean(Ifng_expr[Ifng_positive]), 5), NA_real_),
    median_expr_pos = ifelse(n_pos > 0,
                             round(median(Ifng_expr[Ifng_positive]), 5), NA_real_),
    .groups = "drop"
  )

write.csv(summary_ct_samp,
          file.path(out_dir, "ifng_summary_by_celltype_sample.csv"),
          row.names = FALSE)
message("Saved: ifng_summary_by_celltype_sample.csv")

## 3b. Per cell type (all conditions collapsed)
summary_ct <- cell_df %>%
  group_by(cell_type) %>%
  summarise(
    n_cells          = n(),
    n_pos            = sum(Ifng_positive),
    pct_pos          = round(100 * n_pos / n_cells, 3),
    mean_expr_all    = round(mean(Ifng_expr), 5),
    median_expr_pos  = ifelse(sum(Ifng_positive) > 0,
                              round(median(Ifng_expr[Ifng_positive]), 5), NA_real_),
    max_expr         = round(max(Ifng_expr), 5),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_pos))

write.csv(summary_ct,
          file.path(out_dir, "ifng_summary_by_celltype.csv"),
          row.names = FALSE)
cat("\nSummary by cell type (all conditions):\n")
print(summary_ct)
message("Saved: ifng_summary_by_celltype.csv")

# =============================================================
# 4. Cell type ordering for figures
# =============================================================

all_ct <- unique(cell_df$cell_type)

# Biologically meaningful order:
#   PREFERRED_ORDER patterns first (in listed order), then remaining alphabetically
preferred_matched <- unique(unlist(lapply(PREFERRED_ORDER, function(pat) {
  grep(pat, all_ct, value = TRUE, ignore.case = TRUE)
})))
remaining    <- sort(setdiff(all_ct, preferred_matched))
ct_order_bio <- c(preferred_matched, remaining)

cat("\nCell type ordering for figures (top to bottom):\n")
print(ct_order_bio)

# =============================================================
# 5. Publication theme
# =============================================================

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1,
                                      hjust = 0, margin = margin(b = 4)),
      plot.subtitle    = element_text(size = base_size - 2, color = "grey40",
                                      hjust = 0, lineheight = 1.3,
                                      margin = margin(b = 10)),
      axis.text        = element_text(size = base_size - 1.5),
      axis.title       = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size - 1, face = "bold"),
      legend.text      = element_text(size = base_size - 1.5),
      legend.key.size  = unit(0.45, "cm"),
      panel.grid.major = element_line(color = "grey93", linewidth = 0.3),
      plot.margin      = margin(12, 18, 12, 12)
    )
}

# =============================================================
# 6. Figure 1 -- Bubble plot (global overview, all cell types)
#
#   x    = condition
#   y    = cell type (reversed biological order so T CD8 is at top)
#   size = % Ifng-positive cells
#   color = median Ifng expression among positive cells only (log1p-like scale)
#
#   Grey cross = 0 Ifng+ cells detected in that condition x cell type.
#   Bubble size encodes frequency; color encodes intensity among expressors.
# =============================================================

message("\nFigure 1: bubble plot (all cell types, % Ifng+, median expr among positives)...")

bubble_df <- summary_ct_samp %>%
  mutate(
    cell_type = factor(cell_type, levels = rev(ct_order_bio)),
    sample    = factor(sample, levels = SAMPLE_ORDER)
  )

bubble_df_zero <- bubble_df %>% filter(pct_pos == 0)
bubble_df_pos  <- bubble_df %>% filter(pct_pos >  0)

expr_range <- range(bubble_df_pos$median_expr_pos, na.rm = TRUE)

fig1 <- ggplot() +
  # Grey cross for zero-positive conditions (shows the condition was measured)
  geom_point(
    data  = bubble_df_zero,
    aes(x = sample, y = cell_type),
    shape = 4, size = 1.8, color = "grey75", stroke = 0.7
  ) +
  # Filled bubble for conditions with >= 1 Ifng+ cell
  geom_point(
    data  = bubble_df_pos,
    aes(x = sample, y = cell_type,
        size  = pct_pos,
        color = median_expr_pos),
    alpha = 0.88
  ) +
  scale_x_discrete(labels = sample_labels) +
  scale_size_continuous(
    name   = "% Ifng\u207a cells",
    range  = c(2, 14),
    breaks = c(0.5, 1, 2, 5, 10),
    labels = function(x) paste0(x, "%")
  ) +
  scale_color_gradientn(
    name   = "Median Ifng expression\namong Ifng\u207a cells\n(log1p-like scale)",
    colors = c("#fecc5c", "#fd8d3c", "#e31a1c", "#800026"),
    na.value = "grey80",
    limits = if (nrow(bubble_df_pos) == 0 ||
                 all(is.na(bubble_df_pos$median_expr_pos))) NULL else expr_range
  ) +
  labs(
    title    = "Ifng expression in immune niche cell populations during LCMV infection",
    subtitle = paste0(
      "Bubble size = % of Ifng-positive cells | ",
      "Bubble color = median Ifng expression among positive cells (log1p-like scale)\n",
      "Grey cross (\u00d7) = 0 Ifng\u207a cells in that condition | ",
      "BANKSY sub-clustering on Immune (Acod1) domain (\u03bb\u00a0=\u00a00.2, res\u00a0=\u00a00.3)"
    ),
    x = "Condition", y = NULL
  ) +
  theme_pub() +
  theme(
    legend.position = "right",
    axis.text.x     = element_text(hjust = 0.5, size = 9.5),
    axis.text.y     = element_text(size = 9.5)
  )

n_ct <- length(ct_order_bio)
h1   <- max(5, n_ct * 0.42 + 3)
w1   <- 9.5

ggsave(file.path(out_dir, "fig1_ifng_bubble.pdf"),
       plot = fig1, width = w1, height = h1, device = cairo_pdf)
ggsave(file.path(out_dir, "fig1_ifng_bubble.jpg"),
       plot = fig1, width = w1, height = h1, dpi = 300)
message("  Saved: fig1_ifng_bubble  (", w1, " x ", round(h1, 1), " in)")

# =============================================================
# 7. Figure 2A and 2B -- Strip plots (Ifng+ cells only)
#
#   2A — All annotated cell types with >= 1 Ifng-positive cell.
#   2B — Restricted to focus biological populations:
#        T CD8, T CD4, Mono, Mac, Microglia, Glials, Oligo, Vascular
#        (flexible pattern matching against actual labels in object)
#
#   Each point = one Ifng-positive cell.
#   x = condition, y = Ifng expression value (log1p-like scale).
#   Horizontal bar = median per condition.
#   Right y-axis: approximate transcript count = round(exp(y) - 1)
# =============================================================

# Focus populations for Fig 2B (flexible pattern matching)
FOCUS_B_PATTERNS <- c("T CD8")

# ---------------------------------------------------------------------------
# Helper: build a strip plot of Ifng+ cells for a given set of cell types.
# Arguments:
#   pos_df       — data frame of Ifng-positive cells (cell_type, sample, Ifng_expr)
#   all_cells_df — full cell_df (to compute n_total for facet labels)
#   plot_order   — character vector: desired panel order (left-to-right, top-to-bottom)
#   title_str    — figure title
#   subtitle_str — figure subtitle
# Returns: a ggplot object
# ---------------------------------------------------------------------------
build_strip_plot <- function(pos_df, all_cells_df, plot_order,
                              title_str, subtitle_str) {

  # y-axis ticks from unique expression values present in this subset
  yb    <- sort(unique(pos_df$Ifng_expr))
  ylab1 <- sprintf("%.2f", yb)
  ylab2 <- as.character(round(expm1(yb)))

  # Facet label: "cell_type\n(n_pos Ifng+ / n_total cells)"
  lbl_df <- all_cells_df %>%
    filter(cell_type %in% plot_order) %>%
    group_by(cell_type) %>%
    summarise(n_total = n(), n_pos = sum(Ifng_positive), .groups = "drop") %>%
    mutate(
      strip_label = paste0(cell_type, "\n(",
                           n_pos, " Ifng\u207a / ", n_total, " cells)"),
      cell_type   = factor(cell_type, levels = plot_order)
    )

  lbl_levels <- lbl_df$strip_label[match(plot_order, lbl_df$cell_type)]

  plot_df <- pos_df %>%
    left_join(lbl_df %>% select(cell_type, strip_label), by = "cell_type") %>%
    mutate(
      cell_type   = factor(cell_type,   levels = plot_order),
      sample      = factor(sample,      levels = SAMPLE_ORDER),
      strip_label = factor(strip_label, levels = lbl_levels)
    )

  med_df <- plot_df %>%
    group_by(cell_type, sample, strip_label) %>%
    summarise(med_expr = median(Ifng_expr), .groups = "drop") %>%
    mutate(
      cell_type   = factor(cell_type,   levels = plot_order),
      sample      = factor(sample,      levels = SAMPLE_ORDER),
      strip_label = factor(strip_label, levels = lbl_levels)
    )

  ggplot(plot_df, aes(x = sample, y = Ifng_expr, color = sample)) +
    geom_jitter(width = 0.22, height = 0,
                size = 2.5, alpha = 0.72, shape = 16) +
    geom_crossbar(
      data        = med_df,
      aes(x = sample, y = med_expr, ymin = med_expr, ymax = med_expr),
      width       = 0.55, linewidth = 1.0, color = "grey15",
      middle.linewidth = 0, inherit.aes = FALSE
    ) +
    scale_color_manual(values = sample_colors, labels = sample_labels,
                       name = "Condition") +
    scale_x_discrete(labels = sample_labels) +
    scale_y_continuous(
      name     = "Ifng expression (log1p-like scale)",
      breaks   = yb,
      labels   = ylab1,
      sec.axis = sec_axis(
        transform = ~ expm1(.),
        name      = "Approx. transcript count  [round(exp(x)\u22121)]",
        breaks    = expm1(yb),
        labels    = ylab2
      )
    ) +
    facet_wrap(~ strip_label, ncol = min(3L, length(plot_order)),
               scales = "fixed") +
    labs(title = title_str, subtitle = subtitle_str, x = NULL) +
    theme_pub() +
    theme(
      legend.position  = "top",
      axis.text.x      = element_text(size = 8.5, hjust = 0.5),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text       = element_text(face = "bold", size = 9, lineheight = 1.2)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1)))
}

# ---------------------------------------------------------------------------
# Figure 2A — All annotated cell types with >= 1 Ifng+ cell
# ---------------------------------------------------------------------------
message("\nFigure 2A: strip plot (all cell types with >= 1 Ifng+ cell)...")

ct_2a <- ct_order_bio[ct_order_bio %in%
                        (cell_df %>% filter(Ifng_positive) %>%
                           pull(cell_type) %>% unique())]

cat("\nCell types included in Figure 2A:\n")
print(
  cell_df %>% filter(cell_type %in% ct_2a) %>%
    group_by(cell_type) %>%
    summarise(n_total = n(), n_pos = sum(Ifng_positive), .groups = "drop") %>%
    arrange(match(cell_type, ct_2a))
)

pos_df_2a <- cell_df %>%
  filter(cell_type %in% ct_2a, Ifng_positive) %>%
  mutate(cell_type = factor(cell_type, levels = ct_2a),
         sample    = factor(sample,    levels = SAMPLE_ORDER))

n_2a    <- length(ct_2a)
n_col_a <- min(3L, n_2a)
h_2a    <- max(5, ceiling(n_2a / n_col_a) * 3.8 + 2.5)
w_2a    <- n_col_a * 3.8 + 2.5

fig2a <- build_strip_plot(
  pos_df       = pos_df_2a,
  all_cells_df = cell_df,
  plot_order   = ct_2a,
  title_str    = "Ifng expression intensity in Ifng-positive cells (all annotated cell types)",
  subtitle_str = paste0(
    "Each point = one Ifng-positive cell | Horizontal bar = median per condition\n",
    "Expression: stored log1p-like scale | Right y-axis: approx. transcript count\n",
    "All cell types with \u2265 1 Ifng\u207a cell | BANKSY (\u03bb\u00a0=\u00a00.2, res\u00a0=\u00a00.3)"
  )
)

ggsave(file.path(out_dir, "fig2A_ifng_strip_positives_all_celltypes.pdf"),
       plot = fig2a, width = w_2a, height = h_2a, device = cairo_pdf)
ggsave(file.path(out_dir, "fig2A_ifng_strip_positives_all_celltypes.jpg"),
       plot = fig2a, width = w_2a, height = h_2a, dpi = 300)
message("  Saved: fig2A_ifng_strip_positives_all_celltypes  (",
        w_2a, " x ", round(h_2a, 1), " in)")

# ---------------------------------------------------------------------------
# Figure 2B — Focus biological populations only
# ---------------------------------------------------------------------------
message("\nFigure 2B: strip plot (focus biological populations)...")

focus_b_matched <- unique(unlist(lapply(FOCUS_B_PATTERNS, function(pat) {
  grep(pat, all_ct, value = TRUE, ignore.case = TRUE)
})))

# Keep biological order from ct_order_bio; restrict to those with >= 1 Ifng+ cell
ct_2b <- ct_order_bio[ct_order_bio %in% focus_b_matched &
                        ct_order_bio %in%
                        (cell_df %>% filter(Ifng_positive) %>%
                           pull(cell_type) %>% unique())]

cat("\nFocus populations evaluated for Figure 2B:\n")
print(
  cell_df %>% filter(cell_type %in% focus_b_matched) %>%
    group_by(cell_type) %>%
    summarise(n_total = n(), n_pos = sum(Ifng_positive), .groups = "drop") %>%
    mutate(included = ifelse(cell_type %in% ct_2b,
                             "yes (>= 1 Ifng+)", "no (0 Ifng+)")) %>%
    arrange(match(cell_type, ct_order_bio))
)

if (length(ct_2b) == 0) {

  message("  No focus population has any Ifng+ cells. Skipping Figure 2B.")

} else {

  pos_df_2b <- cell_df %>%
    filter(cell_type %in% ct_2b, Ifng_positive) %>%
    mutate(cell_type = factor(cell_type, levels = ct_2b),
           sample    = factor(sample,    levels = SAMPLE_ORDER))

  n_2b    <- length(ct_2b)
  n_col_b <- min(3L, n_2b)
  h_2b    <- max(5, ceiling(n_2b / n_col_b) * 3.8 + 2.5)
  w_2b    <- n_col_b * 3.8 + 2.5

  fig2b <- build_strip_plot(
    pos_df       = pos_df_2b,
    all_cells_df = cell_df,
    plot_order   = ct_2b,
    title_str    = "Ifng expression intensity in Ifng-positive T CD8 cells",
    subtitle_str = paste0(
      "Each point = one Ifng-positive T CD8 (Gzmb) cell | Horizontal bar = median per condition\n",
      "Expression: stored log1p-like scale | Right y-axis: approx. transcript count\n",
      "BANKSY sub-clustering on Immune (Acod1) domain (\u03bb\u00a0=\u00a00.2, res\u00a0=\u00a00.3)"
    )
  )

  ggsave(file.path(out_dir, "fig2B_ifng_strip_positives_selected_celltypes.pdf"),
         plot = fig2b, width = w_2b, height = h_2b, device = cairo_pdf)
  ggsave(file.path(out_dir, "fig2B_ifng_strip_positives_selected_celltypes.jpg"),
         plot = fig2b, width = w_2b, height = h_2b, dpi = 300)
  message("  Saved: fig2B_ifng_strip_positives_selected_celltypes  (",
          w_2b, " x ", round(h_2b, 1), " in)")
}

# =============================================================
# 8. Figure 3 -- Line plot closely matching the reference style
#
#   x = pseudo-time layout:
#       1, 3, 6 for rLCMV
#       8 for mock 6 wpi (visually separated)
#   y = % Ifng-positive cells
#   one color per cell type
#   infected samples connected
#   mock shown as isolated points only
#
#   3A — all cell types
#   3B — focus populations
# =============================================================

message("\nFigure 3: reference-style line plot (% Ifng-expressing cells)...")

# ---------------------------------------------------------------------------
# X positions to mimic the reference layout
# rLCMV = 1, 3, 6
# mock  = separate point at x = 8
# ---------------------------------------------------------------------------
x_pos_map <- c(
  LCMV_1wpi = 1,
  LCMV_3wpi = 3,
  LCMV_6wpi = 6,
  mock_6wpi = 8
)

line_df <- summary_ct_samp %>%
  mutate(
    x_pos = x_pos_map[as.character(sample)],
    group_type = ifelse(sample == "mock_6wpi", "mock", "rLCMV")
  )

# ---------------------------------------------------------------------------
# Color palette
# ---------------------------------------------------------------------------
n_ct_all <- length(ct_order_bio)
ct_pal <- setNames(
  colorRampPalette(c(
    "#4E79A7", "#A0CBE8", "#B07AA1", "#F1B555", "#EDC948",
    "#E15759", "#76B7B2", "#59A14F", "#FF9DA7", "#17BECF",
    "#D4A6C8", "#C49A6C", "#86BCB6", "#8CD17D", "#FABFD2",
    "#9C755F", "#BAB0AC"
  ))(n_ct_all),
  ct_order_bio
)

# Emphasize T cells if present
if ("T CD8 (Gzmb)" %in% names(ct_pal)) ct_pal["T CD8 (Gzmb)"] <- "#E15759"
if ("T CD4 (Foxp3)" %in% names(ct_pal)) ct_pal["T CD4 (Foxp3)"] <- "#F28E8B"

# Focus populations for Fig 3B
FOCUS_LINE_PATTERNS <- c("T CD8", "T CD4", "Mono", "Mac", "Microglia",
                         "IFN", "Ifit", "Interferon")
focus_line_ct <- ct_order_bio[ct_order_bio %in%
  unique(unlist(lapply(FOCUS_LINE_PATTERNS, function(pat)
    grep(pat, all_ct, value = TRUE, ignore.case = TRUE))))]

# ---------------------------------------------------------------------------
# Helper to build a plot visually close to the example figure
# ---------------------------------------------------------------------------
build_reference_line_plot <- function(df, ct_subset, title_str, y_max = NULL) {

  df_sub <- df %>%
    filter(cell_type %in% ct_subset) %>%
    mutate(cell_type = factor(cell_type, levels = ct_subset))

  df_inf <- df_sub %>% filter(group_type == "rLCMV") %>% arrange(cell_type, x_pos)
  df_mock <- df_sub %>% filter(group_type == "mock")

  pal_sub <- ct_pal[ct_subset]

  p <- ggplot() +
    geom_line(
      data = df_inf,
      aes(x = x_pos, y = pct_pos, color = cell_type, group = cell_type),
      linewidth = 0.55,
      alpha = 0.95
    ) +
    geom_point(
      data = df_inf,
      aes(x = x_pos, y = pct_pos, color = cell_type),
      size = 1.8,
      alpha = 0.95
    ) +
    geom_point(
      data = df_mock,
      aes(x = x_pos, y = pct_pos, color = cell_type),
      size = 1.8,
      alpha = 0.95
    ) +
    scale_color_manual(values = pal_sub, name = NULL) +
    scale_x_continuous(
      breaks = c(1, 3, 6, 8),
      labels = c("1", "3", "6", "6"),
      limits = c(0.6, 8.4),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title = title_str,
      x = "w.p.i.",
      y = "% IFNg expressing cells"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 12, hjust = 0),
      axis.title.x    = element_text(size = 11),
      axis.title.y    = element_text(size = 11),
      axis.text       = element_text(size = 10),
      legend.position = "right",
      legend.text     = element_text(size = 9),
      legend.key.size = unit(0.45, "cm"),
      plot.margin     = margin(10, 18, 18, 10)
    ) +
    guides(color = guide_legend(override.aes = list(linewidth = 1.2, size = 3)))

  if (!is.null(y_max)) {
    p <- p + coord_cartesian(ylim = c(0, y_max))
  }

  # add bottom group labels exactly like the reference logic
  p <- p +
    annotate("text", x = 3.3, y = -0.9, label = "rLCMV", size = 4) +
    annotate("text", x = 8.0, y = -0.9, label = "mock", size = 4)

  p
}

# ---------------------------------------------------------------------------
# 3A — all cell types
# ---------------------------------------------------------------------------
fig3a <- build_reference_line_plot(
  df = line_df,
  ct_subset = ct_order_bio,
  title_str = "Percentage of Ifng-expressing cells by cell type"
)

w3a <- 7.2
h3a <- 4.8

ggsave(file.path(out_dir, "fig3_ifng_pct_lineplot_all.pdf"),
       plot = fig3a, width = w3a, height = h3a, device = cairo_pdf)
ggsave(file.path(out_dir, "fig3_ifng_pct_lineplot_all.jpg"),
       plot = fig3a, width = w3a, height = h3a, dpi = 300)
message("  Saved: fig3_ifng_pct_lineplot_all")

# ---------------------------------------------------------------------------
# 3B — focus populations
# ---------------------------------------------------------------------------
if (length(focus_line_ct) > 0) {
  fig3b <- build_reference_line_plot(
    df = line_df,
    ct_subset = focus_line_ct,
    title_str = "Percentage of Ifng-expressing cells in key immune populations"
  )

  w3b <- 6.5
  h3b <- 4.6

  ggsave(file.path(out_dir, "fig3_ifng_pct_lineplot_focus.pdf"),
         plot = fig3b, width = w3b, height = h3b, device = cairo_pdf)
  ggsave(file.path(out_dir, "fig3_ifng_pct_lineplot_focus.jpg"),
         plot = fig3b, width = w3b, height = h3b, dpi = 300)
  message("  Saved: fig3_ifng_pct_lineplot_focus")
} else {
  message("  No focus cell types found — Fig 3B skipped.")
}

message("\n=== Done. Outputs in: ", out_dir, " ===\n")