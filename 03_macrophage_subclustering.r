# =====================================
# Subcluster macrophages from the integrated CLAD/CIT object.
#
# Workflow:
#   1. Subset macrophages and perform workup
#   2. Remove contaminating cells
#   3. Recluster macrophages
#   4. Resolution optimization using clustree
#   5. Identify cluster markers
#
# Input:
#   results/integrated/CLAD_CIT_integrated_annotated.rds
#
# Output:
#   results/macrophages/macs_subclustered.rds
#   results/macrophages/macs_markers.rds
#   results/macrophages/macs_top10_markers.rds
#
# =====================================

# load packages
library(Seurat)
library(scCustomize)
library(dplyr)
library(ggplot2)
library(clustree)
library(patchwork)
options(future.globals.maxSize = 4000 * 1024^3)

# Helper functions
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

# Load integrated object
lung <- readRDS(
    "results/integrated/CLAD_CIT_integrated_annotated.rds"
)

# =====================================
# Subset macrophage from lung dataset and perform reclustering
# =====================================

macs <- subset(
    lung,
    ident = c(
        "monocyte/macrophage",
        "cycling"
    )
)

rm(lung)

# SCTransform normalization
macs <- SCTransform(
    macs,
    vars.to.regress = c(
        "nFeature_RNA",
        "nCount_RNA",
        "percent.mt"
    )
)

DefaultAssay(macs) <- "RNA"

# Initial reclustering
macs <- FindNeighbors(
    macs,
    reduction = "integrated.cca",
    dims = 1:30
)

macs <- FindClusters(
    macs,
    resolution = 0.5,
    cluster.name = "cca_clusters"
)

macs <- RunUMAP(
    macs,
    reduction = "integrated.cca",
    dims = 1:30,
    reduction.name = "umap.cca"
)

# Remove DoubletFinder metadata
remove_cols <- grep(
    "DF.classification|pANN",
    colnames(macs@meta.data),
    value = TRUE
)
if(length(remove_cols) > 0){

    macs@meta.data <- macs@meta.data[
        ,
        !colnames(macs@meta.data) %in% remove_cols
    ]

}
# Reorder disease groups
macs$Condition <- factor(
    macs$Condition,
    levels = c(
        "BOS",
        "RAS",
        "CIT"
    )
)

# =====================================
# Remove contaminating cell types
# =====================================

#check for contaminating clusters
genecheck(macs, "SFTPC")
genecheck(macs, "EPCAM")
genecheck(macs, "PTPRC")

DimPlot_scCustom(
    macs,
    group.by = "CellAnnotation",
    reduction = "umap.cca"
)


# retain macrophage/cycling labels
# remove EPCAM+, CD3E+, CD79B+ cells

Idents(macs) <- "CellAnnotation"
macs <- subset(
    macs,
    ident = c(
        "macrophage",
        "cycling"
    )
)
contaminating_genes <- c(
    "EPCAM",
    "CD3E",
    "CD79B"
)
for(gene in contaminating_genes){

    macs <- subset(
        macs,
        subset = get(gene) > 0,
        invert = TRUE
    )

}

# Recluster after contaminant removal
macs <- FindNeighbors(
    macs,
    reduction = "integrated.cca",
    dims = 1:30
)

macs <- FindClusters(
    macs,
    resolution = 0.8
)

macs <- RunUMAP(
    macs,
    reduction = "integrated.cca",
    dims = 1:30,
    reduction.name = "umap.cca"
)

# Remove cluster 19 due to off-target contamination
macs <- subset(
    macs,
    ident = "19",
    invert = TRUE
)

# =====================================
# Recluster after contamination removal and clean up meta.data
# =====================================

macs <- FindNeighbors(
    macs,
    reduction = "integrated.cca",
    dims = 1:30
)

macs <- RunUMAP(
    macs,
    reduction = "integrated.cca",
    dims = 1:30,
    reduction.name = "umap.cca"
)

# Include a CLADstatus metadata (CLAD versus CIT)
macs$CLADstatus <- ifelse(
    macs$Condition %in% c(
        "BOS",
        "RAS"
    ),
    "CLAD",
    "CIT"
)

# Remove old clustering resolution columns
resolution_cols <- grep(
    "^RNA_snn_res",
    colnames(macs@meta.data),
    value = TRUE
)

if(length(resolution_cols) > 0){

    macs@meta.data <- macs@meta.data[
        ,
        !colnames(macs@meta.data) %in%
            resolution_cols
    ]

}

# =====================================
# Resolution optimization using Clustree
# =====================================

# Resolution sweep
resolution_range <- seq(
    from = 0.1,
    to = 1.4,
    by = 0.1
)

macs <- FindClusters(
    macs,
    resolution = resolution_range
)

# Clustree evaluation
clustree(
    macs@meta.data,
    prefix = "RNA_snn_res."
)

# Visual inspection
p1 <- DimPlot(
    macs,
    reduction = "umap.cca",
    group.by = "RNA_snn_res.1",
    label = TRUE
) +
    NoLegend()

p2 <- DimPlot(
    macs,
    reduction = "umap.cca",
    split.by = "RNA_snn_res.1",
    group.by = "RNA_snn_res.1",
    label = TRUE,
    ncol = 5
) +
    NoLegend()

p1 | p2

DimPlot_scCustom(
    macs,
    reduction = "umap.cca",
    group.by = "RNA_snn_res.1",
    split.by = "Condition"
) +
    NoLegend()

# Generate UMAPs for all tested resolutions
for(res in resolution_range){

    res_col <- paste0(
        "RNA_snn_res.",
        res
    )

    print(
        DimPlot(
            macs,
            reduction = "umap.cca",
            group.by = res_col,
            label = TRUE
        ) +
        ggtitle(
            paste(
                "Resolution",
                res
            )
        ) +
        NoLegend()
    )
}

# Based on clustree + biological inspection,
# resolution = 1.0 was selected.
macs <- FindClusters(
    macs,
    resolution = 1.0,
    cluster.name = "cca_clusters"
)

# =====================================
# Cluster identification through markers
# =====================================

Idents(macs) <- "cca_clusters"

macs.markers <- FindAllMarkers(
    macs,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
)

macs.TOP10 <- macs.markers %>%
    group_by(cluster) %>%
    slice_max(
        avg_log2FC,
        n = 10
    )

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
    "results/macrophages/macs_subclustered.rds"
)

saveRDS(
    macs.markers,
    "results/macrophages/macs_markers.rds"
)

saveRDS(
    macs.TOP10,
    "results/macrophages/macs_top10_markers.rds"
)

writeLines(
    capture.output(sessionInfo()),
    "results/macrophages/sessionInfo.txt"
)