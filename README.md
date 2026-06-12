# MERFISH — CD8 TRM / Microglia Niche Analysis in Alzheimer's Disease

> **Status : unpublished — work in progress**

## Overview

This repository contains the analysis pipeline for a spatial transcriptomics study
using **MERFISH** (Multiplexed Error-Robust Fluorescence In Situ Hybridization) data.

The project investigates the spatial organization of **CD8 tissue-resident memory T cells (TRM)**
and **microglia** in a mouse model of Alzheimer's disease (LCMV system), with a focus on
inflammatory niche characterization and DAM (Disease-Associated Microglia) signatures.

Conditions:
- `Mock_6wpi` — uninfected control
- `LCMV_1wpi`, `LCMV_3wpi`, `LCMV_6wpi` — LCMV-infected at 1, 3, and 6 weeks post-infection

## Input data

Raw MERFISH data are **not included** in this repository (large files, not yet public).

Expected folder structure:

```
data/
  slide4/
    region_R2/    ← Mock_6wpi   (cell_by_gene.csv, cell_metadata.csv, detected_transcripts.csv, cell_boundaries.parquet)
    region_R3/    ← LCMV_1wpi
  slide2/
    region_R1/    ← LCMV_3wpi
    region_R2/    ← LCMV_6wpi
```

Data were acquired with the **Vizgen MERSCOPE** platform and loaded via `Seurat::LoadVizgen`.
A custom fork of Seurat (`alikhuseynov/seurat`, branch `vizgen_seurat5`) is required for initial loading
(see `scripts/00_load_data.R`, flag `INSTALL_SEURAT`).

## Analysis pipeline

Scripts are numbered in execution order and run from the **project root**.

| Phase | Scripts | Description |
|---|---|---|
| 1 — QC | `00` → `02` | Load raw data, quality filtering, log1p normalization |
| 2 — Global clustering | `03`, `04a` | Seurat clustering + BANKSY spatial clustering (λ=0.2, res=0.9) |
| 3 — Global annotation | `05c`, `05d`, `06`, `07` | UMAP annotation from `ncells_by_sample_lam02_res09_joint_long.csv`, composition, enrichment vs mock |
| 4 — Spatial QC | `09`, `10` | Nearest immune-cell distances, XY spatial annotation maps |
| 5 — Immune sub-clustering | `11` → `14` | Extract Acod1 immune subset, re-normalize, BANKSY re-cluster |
| 6 — Immune annotation | `17` | Annotate immune sub-clusters from manual DEG review |
| 7 — Niche analyses | `15`, `16_*`, `18` → `20`, `22`, `24` | IFN overlay, DEG evolution, inflammatory grid, RFP signal |
| 8 — Microglia | `25`, `26`, `28`, `29`, `30`, `31`, `32`, `33` | DAM signature, dotplots, DEGs by condition |
| 9 — RFP | `34_rebuild_protein_metadata`, `35`, `36` | Anti-RFP signal, fold-change by niche |

## Key outputs

| Output | Description |
|---|---|
| `outputs/banksy/umap_annotated/` | Global BANKSY UMAP figures and annotation tables |
| `outputs/banksy/enrichment_vs_mock/` | Log2FC enrichment tables and plots per domain |
| `outputs/banksy/nearest_immune_distance/` | Per-cell distances to nearest immune cell |
| `outputs/banksy/inflammatory_niche_step2_grid_100um_global/` | 100 µm spatial grid enrichment (script 24) |
| `outputs/banksy/immune_acod1/` | Immune sub-cluster figures, DEGs, domain evolution |
| `outputs/rfp_by_celltype/` | Anti-RFP signal by cell type |
| `outputs/rfp_foldchange_niche/` | RFP fold-change by spatial niche |
| `objects/` | Serialized R objects at each pipeline step (not versioned) |

## Dependencies

- R ≥ 4.3
- Seurat ≥ 5.4.0 / SeuratObject ≥ 5.4.0
- Banksy (Bioconductor)
- SpatialExperiment
- ggplot2, dplyr, data.table, arrow, hdf5r

## Contact

Mélina Farshchi
