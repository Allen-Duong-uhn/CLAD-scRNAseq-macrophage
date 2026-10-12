# =====================================
# Lung scRNA-seq preprocessing pipeline
#
# Workflow:
# 1. Read Cell Ranger outputs
# 2. SoupX ambient RNA correction
# 3. Quality control filtering
# 4. DoubletFinder
# 5. Clustering
# 6. Marker identification

# Input:
#   - filtered matrices
#   - raw matrices
#
# Output:
#   - Clustered Seurat objects
#
# =====================================

set.seed(12345)


# Libraries

library(Seurat)
library(dplyr)
library(purrr)
library(SoupX)
library(SingleR)
library(celldex)
library(DoubletFinder)
library(future)

options(future.globals.maxSize = 4000 * 1024^3)

hpca.se <- HumanPrimaryCellAtlasData()

# Sample metadata

sample_info <- tibble::tribble(
  ~sample_id, ~condition, ~time, ~sample_type, ~filtered_path, ~raw_path,

  "CLAD01", "BOS", "EndStage", "Lung",
  "data/CLAD01/filtered_feature_bc_matrix",
  "data/CLAD01/raw_feature_bc_matrix",

  "CLAD02", "BOS", "EndStage", "Lung",
  "data/CLAD02/filtered_feature_bc_matrix",
  "data/CLAD02/raw_feature_bc_matrix",

  "CLAD03", "BOS", "EndStage", "Lung",
  "data/CLAD03/filtered_feature_bc_matrix",
  "data/CLAD03/raw_feature_bc_matrix",

  "CLAD04", "BOS", "EndStage", "Lung",
  "data/CLAD04/filtered_feature_bc_matrix",
  "data/CLAD04/raw_feature_bc_matrix",

  "CLAD05", "RAS", "EndStage", "Lung",
  "data/CLAD05/filtered_feature_bc_matrix",
  "data/CLAD05/raw_feature_bc_matrix",

  "CLAD06", "RAS", "EndStage", "Lung",
  "data/CLAD06/filtered_feature_bc_matrix",
  "data/CLAD06/raw_feature_bc_matrix",

  "CLAD07", "BOS", "EndStage", "Lung",
  "data/CLAD07/filtered_feature_bc_matrix",
  "data/CLAD07/raw_feature_bc_matrix",

  "CLAD08", "BOS", "EndStage", "Lung",
  "data/CLAD08/filtered_feature_bc_matrix",
  "data/CLAD08/raw_feature_bc_matrix",

  "CLAD09", "RAS", "EndStage", "Lung",
  "data/CLAD09/filtered_feature_bc_matrix",
  "data/CLAD09/raw_feature_bc_matrix",

  "CLAD10", "RAS", "EndStage", "Lung",
  "data/CLAD10/filtered_feature_bc_matrix",
  "data/CLAD10/raw_feature_bc_matrix",

  "CLAD11", "RAS", "EndStage", "Lung",
  "data/CLAD11/filtered_feature_bc_matrix",
  "data/CLAD11/raw_feature_bc_matrix",

  "CLAD12", "RAS", "EndStage", "Lung",
  "data/CLAD12/filtered_feature_bc_matrix",
  "data/CLAD12/raw_feature_bc_matrix",

  "CLAD13", "RAS", "EndStage", "Lung",
  "data/CLAD13/filtered_feature_bc_matrix",
  "data/CLAD13/raw_feature_bc_matrix",

  "CLAD14", "BOS", "EndStage", "Lung",
  "data/CLAD14/filtered_feature_bc_matrix",
  "data/CLAD14/raw_feature_bc_matrix"
)

# Helper functions

remove_ribosomal <- function(count_matrix) {

  rpl <- grep("^RPL", rownames(count_matrix))
  rps <- grep("^RPS", rownames(count_matrix))

  count_matrix[-c(rpl, rps), ]
}

soupx_precluster <- function(obj) {

  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)

  obj <- RunPCA(
    obj,
    features = VariableFeatures(obj)
  )

  obj <- FindNeighbors(obj)
  obj <- FindClusters(obj)

  obj <- RunUMAP(
    obj,
    dims = 1:10
  )

  obj
}

run_soupx <- function(obj,
                      filtered_path,
                      raw_path) {

  filtered_counts <- Read10X(filtered_path)
  raw_counts <- Read10X(raw_path)

  filtered_counts <- remove_ribosomal(filtered_counts)
  raw_counts <- remove_ribosomal(raw_counts)

  sc <- SoupChannel(
    tod = raw_counts,
    toc = filtered_counts
  )

  sc <- setClusters(
    sc,
    obj$seurat_clusters
  )

  sc <- setDR(
    sc,
    obj@reductions$umap@cell.embeddings
  )

  sc <- autoEstCont(sc)

  corrected_counts <- adjustCounts(sc)

  colnames(corrected_counts) <- colnames(obj)

  CreateSeuratObject(
    counts = corrected_counts,
    project = obj@project.name
  )
}

run_qc <- function(obj) {

  obj[["percent.mt"]] <- PercentageFeatureSet(
    obj,
    pattern = "^MT-"
  )

  subset(
    obj,
    subset =
      nFeature_RNA > 200 &
      nFeature_RNA < 5000 &
      nCount_RNA > 1000 &
      nCount_RNA < 30000 &
      percent.mt < 15
  )
}

run_singler <- function(obj) {

  annotations <- SingleR(
    test = GetAssayData(obj),
    ref = hpca.se,
    labels = hpca.se$label.main,
    fine.tune = TRUE,
    genes = "de",
    tune.thresh = 0.05
  )

  obj$SingleR.Mainlabels <- annotations$labels

  obj
}

prepare_doubletfinder <- function(obj) {

  obj <- NormalizeData(obj)

  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 2000
  )

  obj <- ScaleData(obj)

  obj <- RunPCA(obj)

  obj <- RunUMAP(
    obj,
    dims = 1:10
  )

  obj <- FindNeighbors(
    obj,
    dims = 1:10
  )

  obj <- FindClusters(
    obj,
    resolution = 0.5
  )

  obj
}

paramSweep_v3_aw<-function(seu, PCs=1:10, sct = FALSE, num.cores=1) {
  require(Seurat); require(fields);
  ## Set pN-pK param sweep ranges
  pK <- c(0.0005, 0.001, 0.005, seq(0.01,0.3,by=0.01))
  pN <- seq(0.15,0.35,by=0.05)
  
  ## Remove pK values with too few cells
  min.cells <- round(nrow(seu@meta.data)/(1-0.05) - nrow(seu@meta.data))
  pK.test <- round(pK*min.cells)
  pK <- pK[which(pK.test >= 1)]
  
  
  real.cells <- rownames(seu@meta.data)
  data <- seu@assays$RNA3@counts
  n.real.cells <- ncol(data)
  
  ## Extract pre-processing parameters from original data analysis workflow
  orig.commands <- seu@commands
  
  ## Iterate through pN, computing pANN vectors at varying pK
  #no_cores <- detectCores()-1
  if(num.cores>1){
    require(parallel)
    cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)),
                        FUN = paramSweep,
                        n.real.cells,
                        real.cells,
                        pK,
                        pN,
                        data,
                        orig.commands,
                        PCs,
                        sct,mc.cores=num.cores)
    stopCluster(cl)
  }else{
    output2 <- lapply(as.list(1:length(pN)),
                      FUN = paramSweep,
                      n.real.cells,
                      real.cells,
                      pK,
                      pN,
                      data,
                      orig.commands,
                      PCs,
                      sct)
  }
  
  ## Write parallelized output into list
  sweep.res.list <- list()
  list.ind <- 0
  for(i in 1:length(output2)){
    for(j in 1:length(output2[[i]])){
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }
  
  ## Assign names to list of results
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep = "_" ))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
  
}


run_doubletfinder_pipeline <- function(
    obj,
    pcs = 10,
    expected_rate = 0.075,
    sct = FALSE
) {
  # Parameter sweep using the installed DoubletFinder implementation
  sweep.res <- DoubletFinder::paramSweep(
    obj,
    PCs = 1:pcs,
    sct = sct,
    num.cores = 1
  )
  
  sweep.stats <- DoubletFinder::summarizeSweep(
    sweep.res,
    GT = FALSE
  )
  
  bcmvn <- DoubletFinder::find.pK(sweep.stats)
  
  best_pk <- as.numeric(
    as.character(
      bcmvn$pK[which.max(bcmvn$BCmetric)]
    )
  )
  
  # Estimate homotypic doublets
  homotypic.prop <- DoubletFinder::modelHomotypic(
    obj$SingleR.Mainlabels
  )
  
  nExp <- round(expected_rate * ncol(obj))
  nExp.adj <- round(nExp * (1 - homotypic.prop))
  
  # Run DoubletFinder
  obj <- DoubletFinder::doubletFinder(
    obj,
    PCs = 1:pcs,
    pN = 0.25,
    pK = best_pk,
    nExp = nExp.adj,
    reuse.pANN = NULL,
    sct = sct
  )
  
  # Identify the newest classification column
  df_col <- grep(
    "^DF.classifications",
    colnames(obj@meta.data),
    value = TRUE
  )
  
  if (length(df_col) == 0L) {
    stop("No DoubletFinder classification column was created.")
  }
  
  df_col <- tail(df_col, 1)
  
  # Retain singlets
  singlet_cells <- rownames(obj@meta.data)[
    obj@meta.data[[df_col]] == "Singlet"
  ]
  
  obj <- subset(obj, cells = singlet_cells)
  
  return(obj)
}

run_clustering <- function(obj) {

  obj <- RunUMAP(
    obj,
    dims = 1:30
  )

  obj <- FindNeighbors(
    obj,
    dims = 1:30
  )

  obj <- FindClusters(
    obj,
    resolution = 0.5
  )

  obj
}

add_metadata <- function(obj,
                         metadata_row) {

  obj$Condition <- metadata_row$condition
  obj$SampleID <- metadata_row$sample_id
  obj$SampleType <- metadata_row$sample_type
  obj$Time <- metadata_row$time

  obj
}

find_markers <- function(obj) {

  FindAllMarkers(
    obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
}


# =====================================
# Load samples
# =====================================

samples <- map2(
  sample_info$filtered_path,
  sample_info$sample_id,
  ~CreateSeuratObject(
    counts = Read10X(.x),
    project = .y
  )
)

names(samples) <- sample_info$sample_id

# =====================================
# SoupX
# =====================================
samples <- imap(
  samples,
  function(obj, sample_name) {

    row <- sample_info[
      sample_info$sample_id == sample_name,
    ]

    obj <- soupx_precluster(obj)

    run_soupx(
      obj,
      filtered_path = row$filtered_path,
      raw_path = row$raw_path
    )
  }
)
# =====================================
# QC Filtering
# =====================================
samples <- map(samples, run_qc)


# =====================================
# DoubletFinder
# =====================================

# SingleR
samples <- map(
  samples,
  prepare_doubletfinder
)

samples <- map(
  samples,
  run_singler
)

# Doublet Finder
samples <- map(
  samples,
  run_doubletfinder_pipeline
)

# =====================================
# Clustering
# =====================================
samples <- map(
  samples,
  run_clustering
)

# Metadata
samples <- imap(
  samples,
  function(obj, sample_name) {

    row <- sample_info[
      sample_info$sample_id == sample_name,
    ]

    add_metadata(obj, row)
  }
)

# =====================================
# Marker genes
# =====================================
markers <- map(
  samples,
  find_markers
)

top10_markers <- map(
  markers,
  ~ .x %>%
    group_by(cluster) %>%
    slice_max(
      avg_log2FC,
      n = 10
    )
)

# =====================================
# Save outputs
# =====================================
dir.create(
  "results/clustered",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/markers",
  recursive = TRUE,
  showWarnings = FALSE
)

iwalk(
  samples,
  ~saveRDS(
    .x,
    file.path(
      "results/clustered",
      paste0(.y, "_clustered.rds")
    )
  )
)

iwalk(
  markers,
  ~saveRDS(
    .x,
    file.path(
      "results/markers",
      paste0(.y, "_markers.rds")
    )
  )
)

saveRDS(
  sample_info,
  "results/sample_metadata.rds"
)

writeLines(
  capture.output(sessionInfo()),
  "results/sessionInfo_01.txt"
)

