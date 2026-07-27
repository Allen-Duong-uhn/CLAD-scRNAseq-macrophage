# =====================================
# Lung scRNA-seq preprocessing pipeline
#
# Steps:
# 1. Read Cell Ranger outputs
# 2. SoupX ambient RNA correction
# 3. Quality control filtering
# 4. DoubletFinder
# 5. Clustering
# 6. Marker identification
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

# Sample metadata

sample_info <- tibble::tribble(
  ~sample_id, ~condition, ~time, ~sample_type, ~filtered_path, ~raw_path,

  "CLAD040", "BOS", "EndStage", "Lung",
  "data/CLAD040/filtered_feature_bc_matrix",
  "data/CLAD040/raw_feature_bc_matrix",

  "CLAD046", "BOS", "EndStage", "Lung",
  "data/CLAD046/filtered_feature_bc_matrix",
  "data/CLAD046/raw_feature_bc_matrix",

  "CLAD053", "BOS", "EndStage", "Lung",
  "data/CLAD053/filtered_feature_bc_matrix",
  "data/CLAD053/raw_feature_bc_matrix",

  "CLAD062", "BOS", "EndStage", "Lung",
  "data/CLAD062/filtered_feature_bc_matrix",
  "data/CLAD062/raw_feature_bc_matrix",

  "CLAD063", "RAS", "EndStage", "Lung",
  "data/CLAD063/filtered_feature_bc_matrix",
  "data/CLAD063/raw_feature_bc_matrix",

  "CLAD064", "RAS", "EndStage", "Lung",
  "data/CLAD064/filtered_feature_bc_matrix",
  "data/CLAD064/raw_feature_bc_matrix",

  "CLAD068", "BOS", "EndStage", "Lung",
  "data/CLAD068/filtered_feature_bc_matrix",
  "data/CLAD068/raw_feature_bc_matrix",

  "CLAD070", "BOS", "EndStage", "Lung",
  "data/CLAD070/filtered_feature_bc_matrix",
  "data/CLAD070/raw_feature_bc_matrix",

  "CLAD073", "RAS", "EndStage", "Lung",
  "data/CLAD073/filtered_feature_bc_matrix",
  "data/CLAD073/raw_feature_bc_matrix",

  "CLAD080", "RAS", "EndStage", "Lung",
  "data/CLAD080/filtered_feature_bc_matrix",
  "data/CLAD080/raw_feature_bc_matrix",

  "CLAD081", "RAS", "EndStage", "Lung",
  "data/CLAD081/filtered_feature_bc_matrix",
  "data/CLAD081/raw_feature_bc_matrix",

  "CLAD082", "RAS", "EndStage", "Lung",
  "data/CLAD082/filtered_feature_bc_matrix",
  "data/CLAD082/raw_feature_bc_matrix",

  "CLAD087", "RAS", "EndStage", "Lung",
  "data/CLAD087/filtered_feature_bc_matrix",
  "data/CLAD087/raw_feature_bc_matrix",

  "CLAD093", "BOS", "EndStage", "Lung",
  "data/CLAD093/filtered_feature_bc_matrix",
  "data/CLAD093/raw_feature_bc_matrix"
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

run_doubletfinder_pipeline <- function(obj,
                                       pcs = 10,
                                       expected_rate = 0.075) {

  obj[["RNA3"]] <- as(
    obj[["RNA"]],
    Class = "Assay"
  )
  sweep.res <- paramSweep_v3_aw(
    obj,
    PCs = 1:pcs,
    sct = TRUE
  )
  sweep.stats <- summarizeSweep(
    sweep.res,
    GT = FALSE
  )
  bcmvn <- find.pK(sweep.stats)
  best_pk <- as.numeric(
    as.character(
      bcmvn$pK[
        which.max(
          bcmvn$BCmetric
        )
      ]
    )
  )
  homotypic.prop <- modelHomotypic(
    obj$SingleR.Mainlabels
  )
  nExp <- round(
    expected_rate *
      ncol(obj)
  )
  nExp.adj <- round(
    nExp * (1 - homotypic.prop)
  )
  obj <- doubletFinder(
    obj,
    PCs = 1:pcs,
    pN = 0.25,
    pK = best_pk,
    nExp = nExp.adj,
    reuse.pANN = FALSE,
    sct = TRUE
  )
  df_col <- grep(
    "^DF.classifications",
    colnames(obj@meta.data),
    value = TRUE
  )
  obj <- subset(
    obj,
    subset = get(df_col) == "Singlet"
  )
  obj[["RNA3"]] <- NULL
  obj
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
  "results/sessionInfo.txt"
)