# =====================================
# Integration of CLAD and CIT lung scRNA-seq datasets
#
# Input:
#   - Clustered Seurat objects generated from
#     preprocessing pipeline
#
# Output:
#   - Integrated Seurat object
#   - Marker tables
#   - Annotated cell types
#
# 1. Load clustered Seurat objects and merge
# 2. Run pre-integration workup
# 3. Perform CCA integration
# 4. Manual cluster annotation
# 5. Save objects
# =====================================

# Load packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

options(future.globals.maxSize = 4000 * 1024^3)

# =====================================
# Load clustered Seurat objects and merge
# =====================================

sample_ids <- c(
  "CLAD040lung",
  "CLAD046lung",
  "CLAD053lung",
  "CLAD062lung",
  "CLAD063lung",
  "CLAD064lung",
  "CLAD068lung",
  "CLAD070lung",
  "CLAD073lung",
  "CLAD080lung",
  "CLAD081lung",
  "CLAD082lung",
  "CLAD087lung",
  "CLAD093lung",
  "CIT2",
  "CIT3",
  "CIT4",
  "CIT5",
  "CIT6",
  "CIT7"
)

objects <- setNames(
  lapply(
    sample_ids,
    function(id) {
      readRDS(
        file.path(
          "results",
          "clustered",
          paste0(id, "_clustered.rds")
        )
      )
    }
  ),
  sample_ids
)

# Merge all samples
combined <- merge(
  x = objects[[1]],
  y = objects[-1]
)

# Prepare assay layers for Seurat v5 integration
combined <- JoinLayers(combined)

combined[["RNA"]] <- split(
  combined[["RNA"]],
  f = combined$orig.ident
)


# =====================================
# Pre-integration workup analysis
# =====================================
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- FindNeighbors(
  combined,
  dims = 1:30
)
combined <- FindClusters(
  combined,
  resolution = 0.5,
  cluster.name = "unintegrated_clusters"
)
combined <- RunUMAP(
  combined,
  dims = 1:30,
  reduction.name = "umap.unintegrated"
)

# Unintegrated QC plots
DimPlot(
  combined,
  reduction = "umap.unintegrated",
  group.by = "orig.ident"
)

DimPlot(
  combined,
  reduction = "umap.unintegrated",
  group.by = "Condition"
)

# =====================================
# CCA integration
# =====================================
combined <- IntegrateLayers(
  object = combined,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca"
)

combined <- FindNeighbors(
  combined,
  reduction = "integrated.cca",
  dims = 1:30
)

combined <- FindClusters(
  combined,
  resolution = 1.0,
  cluster.name = "cca_clusters"
)

combined <- RunUMAP(
  combined,
  reduction = "integrated.cca",
  dims = 1:30,
  reduction.name = "umap.cca"
)

combined <- JoinLayers(combined)

# Integrated QC plots
DimPlot(
  combined,
  reduction = "umap.cca",
  group.by = "orig.ident"
)

DimPlot(
  combined,
  reduction = "umap.cca",
  group.by = "Condition"
)

# Marker discovery
markers <- FindAllMarkers(
  combined,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top10_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(
    avg_log2FC,
    n = 10
  )
# =====================================
# Manual clusters annotation
# =====================================

# Gene inspection helper
genecheck <- function(object, gene) {

  p1 <- FeaturePlot(
    object,
    features = gene,
    reduction = "umap.cca",
    label = TRUE
  )

  p2 <- VlnPlot(
    object,
    features = gene
  ) +
    NoLegend() +
    ggtitle(NULL)

  p1 / p2
}

# ============================================================
# Marker genes used for annotation

# Pan-Immune
# PTPRC

# T cell
# CD3E
# CD4
# CD8A

# NK Cell
# NCR1
# KIR2DL4

# B cell and plasma cell
# CD79A
# JCHAIN

# Myeloid cells
# FCN1
# MARCO
# CD163
# CD14

# Dendritic cells
# CD74
# FCER1A
# CLEC10A

# Mast cell
# TPSAB1

# Pan-Epithelial
# EPCAM

# Type 1 alveolar cell
# AGER
# CAV1

# Type 2 alveolar cell
# SFTPC
# SFTPA1
# LAMP3
# NAPSA

# Club cells
# SCGB1A1
# SCGB3A2
# CYP2F1

# Ciliated epithelial cell
# FOXJ1
# TUBB4B
# RSPH1
# DNAH5

# Fibroblasts/Stromal
# COL1A1
# COL1A2
# LUM
# DCN

# Endothelial
# PECAM1
# CLDN5
# VWF

# Manual cell type annotation

new.cluster.ids <- c(
  "monocyte/macrophage",  # 0
  "monocyte/macrophage",  # 1
  "T cell",               # 2
  "AT2",                  # 3
  "AT2",                  # 4
  "monocyte/macrophage",  # 5
  "NK cell",              # 6
  "monocyte/macrophage",  # 7
  "T cell",               # 8
  "AT1",                  # 9
  "AT2",                  # 10
  "ciliated cell",        # 11
  "ciliated cell",        # 12
  "endothelial cell",     # 13
  "endothelial cell",     # 14
  "cycling",              # 15
  "B cell",               # 16
  "monocyte/macrophage",  # 17
  "club cell",            # 18
  "fibroblast",           # 19
  "endothelial cell",     # 20
  "endothelial cell",     # 21
  "dendritic cell",       # 22
  "mast cell",            # 23
  "B cell",               # 24
  "monocyte/macrophage"   # 25
)

names(new.cluster.ids) <- levels(combined)

combined <- RenameIdents(
  combined,
  new.cluster.ids
)

combined$CellType <- Idents(combined)

# ============================================================
# Save outputs
# ============================================================

dir.create(
  "results/integrated",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  combined,
  file = "results/integrated/CLAD_CIT_integrated_annotated.rds"
)

saveRDS(
  markers,
  file = "results/integrated/CLAD_CIT_markers.rds"
)

saveRDS(
  top10_markers,
  file = "results/integrated/CLAD_CIT_top10_markers.rds"
)

writeLines(
  capture.output(sessionInfo()),
  "results/integrated/sessionInfo.txt"
)