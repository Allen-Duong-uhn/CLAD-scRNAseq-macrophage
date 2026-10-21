# =====================================
# Script: 09_scenic.R
#
# Workflow:
#   1. Load packages and initialize SCENIC
#   2. Load CLAD macrophages
#   3. Gene filtering
#   4. GENIE3
#   5. SCENIC
#   6. Regulon activity matrix at cluster level
#
# Inputs:
#   results/macrophages/macs_annotated.rds
#
# Outputs:
#   results/scenic/
#
# =====================================

# =====================================
# Libraries
# =====================================

library(Seurat)
library(SCENIC)
library(AUCell)
library(dplyr)
library(ComplexHeatmap)
library(Matrix)

# =====================================
# Set SCENIC Parameters and initialize SCENIC
# =====================================

org <- "hgnc"

dataset_title <- "CLAD_macrophages"

db_dir <- "resources/SCENIC_databases"

output_dir <- "results/scenic"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

data(defaultDbNames)

dbs <- defaultDbNames[[org]]

#ensure hg19-500bp-upstream-7species.mc9nr.feather, hg19-tss-centered-10kb-7species.mc9nr.feather are located in dbDir and dbs
scenicOptions <- initializeScenic(
  org = org,
  dbDir = db_dir,
  dbs = dbs,
  datasetTitle = dataset_title,
  nCores = parallel::detectCores()
)
motifAnnotations_hgnc <- motifAnnotations
scenicOptions <- initializeScenic(
  org = org,
  dbDir = db_dir,
  dbs = dbs,
  datasetTitle = dataset_title,
  nCores = parallel::detectCores()
)

rm(motifAnnotations)


# =====================================
# Load CLAD macrophages and downsample
# =====================================

macs <- readRDS(
  "results/macrophages/macs_annotated.rds"
)

Idents(macs) <- "Condition"

clad_macs <- subset(
  macs,
  ident = c(
    "BOS",
    "RAS"
  )
)

clad_macs_ds <- subset(
  clad_macs,
  downsample = 5000
)

# Add cell metadata
cellInfo <- data.frame(
  cluster =
    clad_macs_ds$MacsLabels_1.0
)

rownames(cellInfo) <- colnames(
  clad_macs_ds
)

# acquire expression matrices to be used for SCENIC
expr_counts <- tryCatch(
  GetAssayData(
    clad_macs_ds,
    assay = "RNA",
    layer = "counts"
  ),
  error = function(e) {
    GetAssayData(
      clad_macs_ds,
      assay = "RNA",
      layer = "data"
    )
  }
)

expr_norm <- GetAssayData(
  clad_macs_ds,
  assay = "RNA",
  layer = "data"
)

# =====================================
# Gene filtering
# =====================================

minCountsPerGene <- 10 * 0.05 * ncol(expr_counts)
minSamples <- ncol(expr_counts) * 0.05

gene_counts <- Matrix::rowSums(expr_counts)
gene_samples <- Matrix::rowSums(expr_counts > 0)

genes_kept <- names(gene_counts)[
  gene_counts >= minCountsPerGene &
    gene_samples >= minSamples
]

exprMat_filtered <- expr_counts[
  genes_kept,
  ,
  drop = FALSE
]
exprMat_filtered <- as.matrix(exprMat_filtered)

# =====================================
# Run GENIE3
# =====================================

runCorrelation(
  exprMat_filtered,
  scenicOptions
)

runGenie3(
  exprMat_filtered,
  scenicOptions
)

exprMat_log <- log2(
  exprMat_filtered + 1
)

# =====================================
# Run SCENIC regulon inference
# =====================================

scenicOptions <-
  runSCENIC_1_coexNetwork2modules(
    scenicOptions
  )

scenicOptions <-
  runSCENIC_2_createRegulons(
    scenicOptions
  )

scenicOptions <-
  runSCENIC_3_scoreCells(
    scenicOptions,
    exprMat_log
  )

regulonAUC <- loadInt(
  scenicOptions,
  "aucell_regulonAUC"
)

regulonAUC <-
  regulonAUC[
    onlyNonDuplicatedExtended(
      rownames(regulonAUC)
    ),
  ]

regulonAUC_mat <- getAUC(
  regulonAUC
)

# =====================================
# Cluster-level regulon activity
# =====================================

regulon_activity_by_cluster <- sapply(

  split(
    rownames(cellInfo),
    cellInfo$cluster
  ),

  function(cells){

    valid_cells <- intersect(
      cells,
      colnames(
        regulonAUC_mat
      )
    )

    if(length(valid_cells) == 0){

      return(
        rep(
          NA,
          nrow(
            regulonAUC_mat
          )
        )
      )

    }

    rowMeans(
      regulonAUC_mat[
        ,
        valid_cells,
        drop = FALSE
      ]
    )

  }

)

regulon_activity_scaled <-
  t(
    scale(
      t(
        regulon_activity_by_cluster
      ),
      center = TRUE,
      scale = TRUE
    )
  )

# =====================================
# Save outputs
# =====================================

saveRDS(
  clad_macs_ds,
  file.path(
    output_dir,
    "clad_macs_downsampled.rds"
  )
)

saveRDS(
  genes_kept,
  file.path(
    output_dir,
    "genes_kept.rds"
  )
)

saveRDS(
  scenicOptions,
  file.path(
    output_dir,
    "scenicOptions.rds"
  )
)

saveRDS(
  regulonAUC,
  file.path(
    output_dir,
    "regulonAUC.rds"
  )
)

saveRDS(
  regulonAUC_mat,
  file.path(
    output_dir,
    "regulonAUC_matrix.rds"
  )
)

saveRDS(
  regulon_activity_by_cluster,
  file.path(
    output_dir,
    "regulon_activity_by_cluster.rds"
  )
)

saveRDS(
  regulon_activity_scaled,
  file.path(
    output_dir,
    "regulon_activity_scaled.rds"
  )
)

saveRDS(
  cellInfo,
  file.path(
    output_dir,
    "cellInfo.rds"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    output_dir,
    "sessionInfo_09.txt"
  )
)
