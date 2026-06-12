# =============================================================
# 00_palette.R
# Project: LCMV MERFISH — TRM-Microglia niche analysis
#
# Palette globale partagée par tous les scripts de figures.
# Usage : source("scripts/00_palette.R")
#
# GLOBAL_PALETTE  — named vector de couleurs (annotation -> hex)
# ANNOTATION_ORDER — vecteur ordonné pour les légendes
# =============================================================

GLOBAL_PALETTE <- c(
  # Populations d'intérêt — couleurs vives et contrastées
  "Immune (Acod1)"                      = "#D62728",
  "IFN responsive (Ifit1)"              = "#FF7F0E",
  "Microglia (P2ry12)"                  = "#1F77B4",

  # Sous-clusters de la niche immune
  "T cells (Gzmb)"                       = "#E41A1C",
  "Mac (Ctss)"                          = "#984EA3",
  "Mono (Lyz2)"                         = "#4DAF4A",
  "Microglia (C1qa)"                    = "#377EB8",
  "B cells (Cd19)"                      = "#A65628",
  "T cell / Neuron doublet / Cycling 1" = "#F781BF",

  # Populations gliales/neuronales — couleurs atténuées
  "Astrocytes (Fgfr3)"                  = "#98DF8A",
  "Astrocytes (Gfap)"                   = "#2CA02C",
  "Oligodendrocytes (Plp1)"             = "#17BECF",
  "OPC (Pdgfra)"                        = "#005F73",
  "Vascular (Cldn5)"                    = "#CEA2F7",
  "Fibroblasts/VSMC (Col1a1)"           = "#9467BD",
  "Choroid plexus Epi (Enpp2)"          = "#FFD700",
  "Prolif neural/glial (Ccdc153)"       = "#BCBD22",
  "Neurons (Adora2a)"                   = "#AEC7E8",
  "Neurons (Arhgap36)"                  = "#C7C7C7",
  "Neurons (Crhbp)"                     = "#DBDB8D",
  "Neurons (Dkk3)"                      = "#9EDAE5",
  "Neurons (Htr2c)"                     = "#F7B6D2",
  "Neurons (Nefm)"                      = "#C49C94",
  "Neurons (Neurod1)"                   = "#FFBB78",
  "Neurons (Rxfp1)"                     = "#D9D9D9",
  "Neurons (Slc17a8)"                   = "#E7CB94",
  "Non annote"                          = "#AAAAAA",

  # Sous-clusters IFN et microglie
  "Excitatory neurons (Satb2)"          = "#AEC7E8",
  "Inhibitory neurons (Htr2c)"          = "#F7B6D2",
  "Neurons (Fam107a)"                   = "#C7C7C7",
  "Neurons (Rbfox3)"                    = "#DBDB8D",
  "Glials (Gja1)"                       = "#98DF8A",
  "Oligo"                               = "#17BECF",
  "Oligo ? (Nkx2-9)"                   = "#9EDAE5",
  "OPC (Sox10) / Cycling 2"             = "#005F73",
  "Vascular"                            = "#CEA2F7",
  "Vascular (Igfbp2)"                   = "#9467BD",
  "Unknown"                             = "#D9D9D9"
)

ANNOTATION_ORDER <- c(
  "Immune (Acod1)",
  "IFN responsive (Ifit1)",
  "Microglia (P2ry12)",
  "T cells (Gzmb)", "Mac (Ctss)",
  "Mono (Lyz2)", "Microglia (C1qa)", "B cells (Cd19)",
  "T cell / Neuron doublet / Cycling 1",
  "Astrocytes (Fgfr3)", "Astrocytes (Gfap)",
  "Oligodendrocytes (Plp1)", "OPC (Pdgfra)",
  "Vascular (Cldn5)", "Fibroblasts/VSMC (Col1a1)",
  "Choroid plexus Epi (Enpp2)",
  "Prolif neural/glial (Ccdc153)",
  "Neurons (Adora2a)", "Neurons (Arhgap36)",
  "Neurons (Crhbp)", "Neurons (Dkk3)",
  "Neurons (Htr2c)", "Neurons (Nefm)",
  "Neurons (Neurod1)", "Neurons (Rxfp1)",
  "Neurons (Slc17a8)",
  "Non annote"
)

# Ordre étendu (inclut sous-clusters IFN / immune)
ANNOTATION_ORDER_EXTENDED <- c(
  ANNOTATION_ORDER,
  "Excitatory neurons (Satb2)",
  "Inhibitory neurons (Htr2c)",
  "Neurons (Fam107a)",
  "Neurons (Rbfox3)",
  "Glials (Gja1)",
  "Oligo",
  "Oligo ? (Nkx2-9)",
  "OPC (Sox10) / Cycling 2",
  "Vascular",
  "Vascular (Igfbp2)",
  "Unknown"
)

# Helper : ordonne un vecteur d'annotations selon ANNOTATION_ORDER(_EXTENDED)
# puis place en fin les annotations inconnues, avec na.value = "grey70" pour
# les couleurs manquantes dans GLOBAL_PALETTE.
order_annotations <- function(annots, extended = FALSE) {
  ref <- if (extended) ANNOTATION_ORDER_EXTENDED else ANNOTATION_ORDER
  known   <- intersect(ref, annots)
  unknown <- setdiff(annots, ref)
  c(known, sort(unknown))
}
