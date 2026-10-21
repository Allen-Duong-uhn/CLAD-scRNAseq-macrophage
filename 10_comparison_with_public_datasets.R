# =====================================
#   Differential expression analyses of macrophages.
#
# Workflow:
#   1. Preparation of Yan et al. and Khatri et al. CLAD lung datasets
#   2. Subsetting of public CLAD macrophages and analysis
#   3. Integration of public_CLAD macrophages with Toronto_CLAD macrophages
#   4. Preparation of Reyfman et al. and Morse et al. IPF lung datasets
#   5. Subsetting of public IPF macrophages and analysis
#   6. Integration of IPF macrophages with Toronto_CLAD macrophages
#   7. Save outputs
#
# Inputs:
#   results/macrophages/macs_annotated.rds
#
# Outputs:
#   results/
#
# =====================================

# load packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(scCustomize)

# Helper Functions
genecheck <- function(object, gene) {
  
  reduction <- if ("umap.cca" %in% Reductions(object)) {
    "umap.cca"
  } else if ("umap.harmony" %in% Reductions(object)) {
    "umap.harmony"
  } else {
    stop("Neither 'umap.cca' nor 'umap.harmony' is available in the object.")
  }
  
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

qc_filter <- function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 &
                 nCount_RNA > 1000 & nCount_RNA < 30000 &
                 percent.mt < 15)
}

process_seurat <- function(obj) {
  obj %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30) %>%
    FindNeighbors(dims = 1:30) %>%
    FindClusters(resolution = 0.5)
}


# =====================================
# Prepare Yan et al. and Khatri et al. CLAD lung datasets
# =====================================

# File paths (EDIT THESE)
data_dir <- "path/to/data"

# Load and prepare Yan et al. datasets
yan_paths <- c("CLAD1", "CLAD2", "CLAD3", "CLAD5")
Yan_list <- lapply(yan_paths, function(sample_name) {
  counts <- Read10X(
    file.path(data_dir, sample_name)
  )
  obj <- CreateSeuratObject(
    counts,
    project = sample_name
  )
  obj$SampleID <- obj$orig.ident
  obj$SampleType <- "Lung"
  obj <- qc_filter(obj)
  obj <- process_seurat(obj)
  return(obj)
})
names(Yan_list) <- yan_paths

# Load and prepare Khatri et al. dataset
norm_mat <- read.csv(
  file.path(data_dir, "GSE224210_Normalized_Count_Matrix.csv.gz"),
  row.names = 1
)

Khatri_lung <- CreateSeuratObject(counts = as.matrix(norm_mat))

Khatri_lung <- SetAssayData(
  object = Khatri_lung,
  assay = "RNA",
  layer = "data",
  new.data = as.matrix(norm_mat)
)

# preprocessing
Khatri_lung <- Khatri_lung %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.5, cluster.name = "unintegrated_clusters") %>%
  RunUMAP(dims = 1:30,reduction.name = "umap.unintegrated")

# Add metadata
IDs <- read.csv(
  file.path(data_dir, "GSE224210_Cell_Identities.csv.gz"),
  row.names = 1,
  header = TRUE,
  check.names = FALSE
)
rownames(IDs) <- gsub("-", ".", rownames(IDs))
common <- intersect(
  colnames(Khatri_lung),
  rownames(IDs)
)
Khatri_lung <- AddMetaData(Khatri_lung, IDs[common, ])

#copy 'sample ID' to "SampleID" metadata
Khatri_lung$SampleID <- Khatri_lung$`sample ID`
Khatri_lung$`sample ID` <- NULL

# remove non-CLAD samples
Khatri_lung <- SetIdent(Khatri_lung, value = "SampleID")
Khatri_lung <- subset(Khatri_lung, ident = c("20-469_DT", "21-287_DT", "21-361_DT"), invert = TRUE)
table(Khatri_lung$SampleID)

# Merge and integrate Yan and Khatri data sets
public_CLAD <- merge(Khatri_lung, y = Yan_list)
public_CLAD <- JoinLayers(public_CLAD)
public_CLAD[["RNA"]] <- split(
  public_CLAD[["RNA"]],
  f = public_CLAD$SampleID
)

# preprocessing before integration
public_CLAD <- public_CLAD %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.5, cluster.name = "unintegrated_clusters") %>%
  RunUMAP(dims = 1:30, reduction.name = "umap.unintegrated")

# Harmony integration
public_CLAD <- IntegrateLayers(public_CLAD, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "integrated.harmony")

public_CLAD <- public_CLAD %>% 
  FindNeighbors(reduction = "integrated.harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.5, cluster.name = "harmony_clusters") %>% 
  RunUMAP(reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony") %>% 
  JoinLayers()

# Add "StudyID" meta.data
public_CLAD$StudyID <- ifelse(
  public_CLAD$SampleID %in% c("18-577", "20-469", "21-013", "21-361"),
  "Khatri et al.",
  "Yan et al."
)

#public_CLAD_merge annotation
public_CLAD.markers <- FindAllMarkers(
  public_CLAD,
	only.pos = TRUE,
	min.pct = 0.25,
	logfc.threshold = 0.25
)

public_CLAD.markers %>% 
	group_by(cluster) %>% 
	top_n(n = 10, wt = avg_log2FC
)
public_CLAD.TOP10 <- public_CLAD.markers %>% 
	group_by(cluster) %>% 
	top_n(n = 10, wt = avg_log2FC
)

public_CLAD.TOP10 <- public_CLAD.markers %>%
	group_by(cluster) %>%
	slice_max(order_by = avg_log2FC, n = 10)

#Manual cluster annotations
#Cluster annotations were assigned based on canonical marker expression.
#Marker genes used for selected annotations are indicated in comments.

public_CLAD.cluster.ids <- c(
"T cell", # 0: CD3E
"macrophage/monocyte", # 1: MARCO
"T cell", # 2
"AT2", # 3: SFTPC
"neutrophil", # 4
"endothelial cell", # 5: low VWF, high PECAM1
"endothelial cell", # 6: VWF
"AT1", # 7: AGER
"fibroblast", # 8: COL1A1
"macrophage/monocyte", # 9: monocyte-like, FCN1
"NK cell", # 10: NCR1
"macrophage/monocyte", # 11
"red blood cell", # 12: HBB
"fibroblast", # 13
"AT2", # 14
"AT1", # 15
"dendritic cell", # 16: HLA-DRA; lack of MARCO/FCN1
"smooth muscle cell", # 17: ACTA2
"endothelial cell", # 18
"endothelial cell", # 19
"B cell", # 20: plasma-cell-like, JCHAIN
"ciliated cell", # 21: DNAH5
"B cell", # 22
"endothelial cell", # 23: lymphatic, PROX1
"mast cell", # 24
"club cell", # 25
"cycling", # 26
"macrophage/monocyte", # 27
"macrophage/monocyte", # 28
"macrophage/monocyte", # 29
"AT2", # 30
"T cell" # 31
)

names(public_CLAD.cluster.ids) <- levels(public_CLAD_merge)

public_CLAD_merge <- RenameIdents(
	public_CLAD_merge,
	public_CLAD.cluster.ids
)

public_CLAD_merge$CellType <- Idents(public_CLAD_merge)

# =====================================
# Identifying and subsetting public CLAD macrophages 
# =====================================

#identify and subset myeloid cells
VlnPlot(public_CLAD, features = c("MARCO", "CD14", "CD68", "FCN1", "LGMN"))
FeaturePlot(public_CLAD, reduction = "umap.harmony", features = c("MARCO", "CD14", "CD68", "FCN1", "LGMN"), label = TRUE)
macs_publicCLAD <- subset(public_CLAD, ident = c("1", "9", "11", "27", "28", "29"))

macs_publicCLAD <- FindNeighbors(macs_publicCLAD, reduction = "integrated.harmony", dims = 1:30)
macs_publicCLAD <- FindClusters(macs_publicCLAD, resolution = 1.2, cluster.name = "harmony_clusters")
macs_publicCLAD <- RunUMAP(macs_publicCLAD, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony")

#identifying and removing contaminating clusters
#contamination was assessed using canonical markers
genecheck(macs_publicCLAD, "CD3E") #cluster 12 is high
genecheck(macs_publicCLAD, "EPCAM") #cluster 17 and 22 is high
genecheck(macs_publicCLAD, "CD79A") #none
genecheck(macs_publicCLAD, "PTPRC") #cluster 22 does not express CD45 at all. 
                                    #Cluster 3, 8, and 19 had some expression but majority of cells have none.
macs_publicCLAD <- subset(macs_publicCLAD, ident = c("12", "17", "22"), invert = TRUE)

# recalculate macrophage clustering after contamination removal
macs_publicCLAD <- FindNeighbors(macs_publicCLAD, reduction = "integrated.harmony", dims = 1:30)
macs_publicCLAD <- FindClusters(macs_publicCLAD, resolution = 1.0, cluster.name = "harmony_clusters")
macs_publicCLAD <- RunUMAP(macs_publicCLAD, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony")
DimPlot(macs_publicCLAD, reduction = "umap.harmony", label = TRUE)

# =====================================
# Integration of public_CLAD macrophages with Toronto_CLAD macrophages
# =====================================

# load macs
macs <- readRDS(
  "results/macrophages/macs_subclustered.rds"
)

# subset only CLAD samples
macs_TorontoCLAD <- subset(
  macs,
  subset = Condition %in% c(
    "BOS",
    "RAS"
  )
)
rm(macs)
#add StudyID
macs_TorontoCLAD$StudyID <- "Duong et al."

# merge datasets
macs_mergeCLAD<- merge(macs_TorontoCLAD, y = c(macs_publicCLAD))

# prepare RNA layers for integration
macs_mergeCLAD <- JoinLayers(
  macs_mergeCLAD
)

macs_mergeCLAD[["RNA"]] <- split(
  macs_mergeCLAD[["RNA"]],
  f = macs_mergeCLAD$SampleID
)

# Unintegrated analysis
macs_mergeCLAD <- macs_mergeCLAD %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(
    dims = 1:30
  ) %>%
  FindClusters(
    resolution = 0.5,
    cluster.name = "unintegrated_clusters"
  ) %>%
  RunUMAP(
    dims = 1:30,
    reduction.name = "umap.unintegrated"
  )

# Harmony integration
macs_mergeCLAD <- IntegrateLayers(
  object = macs_mergeCLAD,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.harmony"
)

macs_mergeCLAD <- macs_mergeCLAD %>%
  FindNeighbors(
    reduction = "integrated.harmony",
    dims = 1:30
  ) %>%
  FindClusters(
    resolution = 0.5,
    cluster.name = "harmony_clusters"
  ) %>%
  RunUMAP(
    reduction = "integrated.harmony",
    dims = 1:30,
    reduction.name = "umap.harmony"
  ) %>%
  JoinLayers()

# remove legacy reduction if present
if ("umap.cca" %in% Reductions(
  macs_mergeCLAD
)) {
  macs_mergeCLAD[["umap.cca"]] <- NULL
}


# Visualize integration by sample
DimPlot(
  macs_mergeCLAD,
  reduction = "umap.harmony",
  group.by = c("StudyID"))
DimPlot(macs_mergeCLAD, reduction = 'umap.harmony', label = TRUE)

# =====================================
# Prepare Refyman et al. and Morse et al. IPF lung datasets
# =====================================

## File paths (EDIT THESE)
data_dir <- "path/to/GSE122960"

# Load and prepare Reyfman et al. datasets
reyfman_paths <- c(
  IPF1 = file.path(data_dir,"GSM3489183_IPF_01_filtered_gene_bc_matrices_h5.h5"),
  IPF2 = file.path(data_dir,"GSM3489184_IPF_02_filtered_gene_bc_matrices_h5.h5"),
  IPF3 = file.path(data_dir,"GSM3489188_IPF_03_filtered_gene_bc_matrices_h5.h5"),
  IPF4 = file.path(data_dir,"GSM3489190_IPF_04_filtered_gene_bc_matrices_h5.h5")
                   )

Reyfman_list <- lapply(names(reyfman_paths), function(sample_name) {
	counts <- Read10X_h5(
		reyfman_paths[[sample_name]],
		use.names = TRUE,
		unique.features = TRUE
	)

	obj <- CreateSeuratObject(
		counts = counts,
		project = sample_name
	)

	obj$SampleID <- obj$orig.ident
	obj$SampleType <- "Lung"
	obj$Condition <- "IPF"
	obj <- qc_filter(obj)
	obj <- process_seurat(obj)

	return(obj)
})
names(Reyfman_list) <- names(reyfman_paths)

# Load and prepare Morse et al. datasets
data_dir <- "path/to/GSE128033"

morse_paths <- c(
"SC87",
"SC88",
"SC89",
"SC93",
"SC94",
"SC95",
"SC153",
"SC154"
)

Morse_list <- lapply(morse_paths, function(sample_name) {

	counts <- Read10X(
	file.path(data_dir, paste0(sample_name, "_IPF"))
	)

	obj <- CreateSeuratObject(
	counts = counts,
	project = sample_name
	)

	obj$SampleID <- obj$orig.ident
	obj$SampleType <- "Lung"
	obj$Condition <- "IPF"
	obj <- qc_filter(obj)
	obj <- process_seurat(obj)

	return(obj)
})

names(Morse_list) <- morse_paths


## Integration of IPF
IPF <- merge(Reyfman_list[[1]], y = c(Reyfman_list[-1], Morse_list))
IPF <- JoinLayers(IPF)
IPF[["RNA"]] <- split(
IPF[["RNA"]],
f = IPF$SampleID
)

# preprocessing before integration
IPF <- IPF %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.5, cluster.name = "unintegrated_clusters") %>%
  RunUMAP(dims = 1:30, reduction.name = "umap.unintegrated")

# Harmony integration
IPF <- IntegrateLayers(IPF, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "integrated.harmony")

IPF <- IPF %>% 
  FindNeighbors(reduction = "integrated.harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.5, cluster.name = "harmony_clusters") %>% 
  RunUMAP(reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony") %>% 
  JoinLayers()

# Visualize integrated dataset
DimPlot(IPF, reduction = "umap.harmony", group.by = "SampleID")
DimPlot(IPF, reduction = "umap.harmony", label = TRUE)

# Identify top markers to aid in IPF annotation
IPF.markers <- FindAllMarkers(IPF, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
IPF.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
IPF.TOP10 <- IPF.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#Manual cluster annotations
#Cluster annotations were assigned based on canonical marker expression.
#Marker genes used for selected annotations are indicated in comments.
ipf.cluster.ids <- c("macrophage/monocyte" #0
                    ,"endothelial cell" #1 
                    , "macrophage/monocyte" #2
                    , "Alveolar epi cells" #3
                    , "T cell" #4
                    , "fibroblast" #5
                    , "T cell" #6
                    , "macrophage/monocyte" #7
                    , "dendritic cell" #8
                    , "mast cell" #9
                    , "ciliated cell" #10
                    , "basal cell" #11
                    , "club cell" #12
                    , "cycling" #13
                    , "B cell" #14
                    , "B cell" #15
                    , "smooth muscle cell" #16 vascular
                    , "endothelial cell" #17 lymphatic 
                    , "B cell" #18
)
names(ipf.cluster.ids) <- levels(IPF)
IPF <- RenameIdents(IPF, ipf.cluster.ids)
IPF$CellType <- Idents(IPF)

DimPlot(IPF, label = TRUE, reduction = "umap.harmony")

# =====================================
# Identifying and subsetting IPF macrophages 
# =====================================

#subset myeloid cells
macs_IPF <- subset(IPF, ident = c("macrophage/monocyte", "cycling"))

macs_IPF <- FindNeighbors(macs_IPF, reduction = "integrated.harmony", dims = 1:30)
macs_IPF <- FindClusters(macs_IPF, resolution = 1.0, cluster.name = "harmony_clusters")
macs_IPF <- RunUMAP(macs_IPF, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony")

#identifying and removing contaminating clusters
#contamination was assessed using canonical markers
genecheck(macs_IPF, "CD3E") #
genecheck(macs_IPF, "EPCAM") #
genecheck(macs_IPF, "CD79A") #
genecheck(macs_IPF, "PTPRC") # 
genecheck(macs_IPF, "MARCO") #
genecheck(macs_IPF, "CD14") #
genecheck(macs_IPF, "MKI67") # 
macs_IPF <- subset(macs_IPF, ident = c("12", "13", "14"), invert = TRUE)

# recalculate macrophage clustering after contamination removal
macs_IPF <- FindNeighbors(macs_IPF, reduction = "integrated.harmony", dims = 1:30)
macs_IPF <- FindClusters(macs_IPF, resolution = 1.0, cluster.name = "harmony_clusters")
macs_IPF <- RunUMAP(macs_IPF, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony")
DimPlot(macs_IPF, reduction = "umap.harmony", label = TRUE)

# =====================================
# Integration of IPF macrophages with Toronto_CLAD macrophages
# =====================================

# merge datasets
macs_merge_CLAD_IPF<- merge(macs_TorontoCLAD, y = c(macs_IPF))

# prepare layers for integration
class(macs_merge_CLAD_IPF[["RNA"]])
macs_merge_CLAD_IPF <- JoinLayers(macs_merge_CLAD_IPF)
macs_merge_CLAD_IPF[["RNA"]] <- split(macs_merge_CLAD_IPF[["RNA"]], f = macs_merge_CLAD_IPF$SampleID)
class(macs_merge_CLAD_IPF[["RNA"]])
Layers(macs_merge_CLAD_IPF[["RNA"]])

# unintegrated analysis
macs_merge_CLAD_IPF <- macs_merge_CLAD_IPF %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.5, cluster.name = "unintegrated_clusters") %>%
  RunUMAP(dims = 1:30, reduction.name = "umap.unintegrated")
DimPlot(macs_merge_CLAD_IPF, group.by = "SampleID")

# Harmony integration
macs_merge_CLAD_IPF <- IntegrateLayers(macs_merge_CLAD_IPF, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "integrated.harmony")

macs_merge_CLAD_IPF <- macs_merge_CLAD_IPF %>% 
  FindNeighbors(reduction = "integrated.harmony", dims = 1:30) %>% 
  FindClusters(resolution = 1.0, cluster.name = "harmony_clusters") %>% 
  RunUMAP(reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony") %>% 
  JoinLayers()

# -------------------------------------
# Save outputs
# -------------------------------------

dir.create(
  "results/public_datasets",
  recursive = TRUE,
  showWarnings = FALSE
)

# Public CLAD datasets
saveRDS(
  public_CLAD,
  "results/public_datasets/public_CLAD.rds"
)

saveRDS(
  macs_publicCLAD,
  "results/public_datasets/macs_publicCLAD.rds"
)

saveRDS(
  macs_mergeCLAD,
  "results/public_datasets/macs_mergeCLAD.rds"
)

# IPF datasets
saveRDS(
  IPF,
  "results/public_datasets/IPF.rds"
)

saveRDS(
  macs_IPF,
  "results/public_datasets/macs_IPF.rds"
)

saveRDS(
  macs_merge_CLAD_IPF,
  "results/public_datasets/macs_merge_CLAD_IPF.rds"
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  "results/public_datasets/sessionInfo.txt"
)