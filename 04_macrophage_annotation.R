# =====================================
# Annotate macrophage subclusters.
#
# Workflow:
#   1. Assess low-resolution (0.1) clustering and annotate
#   2. Assess high-resolution (1.0) clustering and annotate
#   3. Create downsampled object
#   4. Save outputs
#
# Inputs:
#   results/macrophages/macs_subclustered.rds
#
# Outputs:
#   results/macrophages/macs_annotated.rds
#   results/macrophages/macs_markers.rds
#
# =====================================

# load packages
library(Seurat)
library(scCustomize)
library(dplyr)
library(ggplot2)
library(patchwork)

# Helper Functions
genecheck <- function(
    object,
    gene,
    reduction = "umap.cca"
){

    p1 <- FeaturePlot_scCustom(
        object,
        features = gene,
        label = TRUE,
        reduction = reduction
    )

    p2 <- VlnPlot(
        object,
        features = gene
    )

    p1 / p2
}

# load macrophage object and markers
macs <- readRDS(
    "results/macrophages/macs_subclustered.rds"
)
 macs.markers <- readRDS(
	"results/macrophage/macs_markers.rds"
)

# =====================================
# Low-resolution (0.1) annotation
# =====================================

Idents(macs) <- "RNA_snn_res.0.1"

DimPlot_scCustom(
    macs,
    reduction = "umap.cca",
    label = TRUE
)

FeaturePlot_scCustom(
    macs,
    reduction = "umap.cca",
    label = TRUE,
    features = c(
        "FABP4",
        "MARCO",
        "LGMN",
        "FCN1",
        "CD14",
        "MKI67"
    )
)

macs_res_0.1_ids <- c(
    "alveolar macrophage",
    "interstitial macrophage",
    "monocyte/recently derived macrophage",
    "cycling macrophages"
)

names(macs_res_0.1_ids) <- levels(macs)

macs <- RenameIdents(
    macs,
    macs_res_0.1_ids
)

macs$MacsLabels_0.1 <- Idents(macs)

# Cell cycle scoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

macs <- CellCycleScoring(
    macs,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = TRUE
)

# =====================================
# High-resolution (1.0) annotation
# =====================================

# Return to final clustering
Idents(macs) <- "cca_clusters"

# Heatmap of cluster markers
macs.TOP5 <- macs.markers %>%
    group_by(cluster) %>%
    slice_max(
        avg_log2FC,
        n = 5
    )

DoHeatmap(
    macs,
    features = macs.TOP5$gene
) +
    NoLegend() +
    scale_fill_gradientn(
        colours = c(
            "blue",
            "white",
            "red"
        )
    )

#use genecheck helper function to assess distribution of expression
genecheck(macs, "IFI27")
genecheck(macs, "SPP1")
genecheck(macs, "CXCL8")
genecheck(macs, "LILRA5")

# Marker genes used for macrophage subtype annotation

# IFI27+ AMs
# IFI27

# FN1+ AMs
# FN1
# LYAR

# FOLR2+ IMs
# FOLR2
# CCL8

# F13A1+ IMs
# F13A1
# FCER2
# CCL13

# Monocytes
# S100A8
# S100A9
# FCN1
# CD14

# CXCL8+ Monocytes
# CXCL8
# FCGR2B
# FCN1
# S100A8
# S100A9

# CLEC10A+ macrophages
# CLEC10A

# Interferon-stimulated macrophages (ISGs)
# CXCL10
# CCL2
# CCL8
# ISG15
# IFIT1
# IL1RN

# Galectin-2 macrophages
# LGALS2

# LILR macrophages
# LILRA2
# LILRA5
# LILRB2

# SPP1 macrophages
# SPP1
# MMP9
# CHI3L1

# Proliferating macrophages
# MKI67
# TOP2A
# UBE2C

# Main annotation
macs_res_1.0_ids <- c(
  "IFI27+ AMs",        # 0
  "IFI27+ AMs",        # 1
  "FN1+ AMs",          # 2
  "IFI27+ AMs",        # 3
  "FOLR2+ IMs",        # 4
  "Monocytes",         # 5
  "CLEC10A+",          # 6
  "ISGs",              # 7
  "FN1+ AMs",          # 8
  "Galectin-2",        # 9
  "LILR",              # 10
  "CXCL8+ Monocytes",  # 11
  "FN1+ AMs",          # 12
  "Prolif",            # 13
  "F13A1+ IMs",        # 14
  "SPP1",              # 15
  "SPP1",              # 16
  "CXCL8+ Monocytes"   # 17 - singleton cluster; merged with cluster 11
)

names(macs_res_1.0_ids) <- levels(macs)

macs <- RenameIdents(
    macs,
    macs_res_1.0_ids
)

macs$MacsLabels_1.0 <- Idents(macs)

# Visualize annotations
DimPlot_scCustom(
    macs,
    reduction = "umap.cca",
    label = TRUE
)

# Markers by macrophage subtype
Idents(macs) <- "MacsLabels_1.0"

macs.subtype.markers <- FindAllMarkers(
    macs,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
)

macs.subtype.TOP5 <- macs.subtype.markers %>%
    group_by(cluster) %>%
    slice_max(
        avg_log2FC,
        n = 5
    )

# =====================================
# Create downsampled object for balanced downstream visualization
# =====================================

Idents(macs) <- "Condition"

macs_downsample <- subset(
    macs,
    downsample = 4000
)

Idents(macs) <- "MacsLabels_1.0"
Idents(macs_downsample) <- "MacsLabels_1.0"

# =====================================
# Save outputs
# =====================================

dir.create(
    "results/macrophages",
    recursive = TRUE,
    showWarnings = FALSE
)

saveRDS(
    macs,
    "results/macrophages/macs_annotated.rds"
)

saveRDS(
    macs_downsample,
    "results/macrophages/macs_annotated_downsample.rds"
)


saveRDS(
    macs.subtype.markers,
    "results/macrophages/macs_subtype_markers.rds"
)

saveRDS(
    macs.TOP5,
    "results/macrophages/macs_top5_markers.rds"
)

writeLines(
    capture.output(sessionInfo()),
    "results/macrophages/sessionInfo_04.txt"
)
