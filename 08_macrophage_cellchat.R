# =====================================
# Script: 08_cellchat.R
#
# Purpose:
#   CellChat analysis of CLAD lungs.
#
# Inputs:
#   results/integrated/CLAD_CIT_integrated_annotated.rds
#   results/macrophages/macs_annotated.rds
#
# Outputs:
#   results/cellchat/
#
# =====================================

# =====================================
# Libraries
# =====================================

library(Seurat)
library(CellChat)

# =====================================
# Helper function
# =====================================

run_cellchat <- function(
  object,
  group.by
){

  cc <- createCellChat(
    object = object,
    group.by = group.by,
    assay = "RNA"
  )

  cc@DB <- CellChatDB.human

  cc <- subsetData(cc)

  cc <- identifyOverExpressedGenes(cc)

  cc <- identifyOverExpressedInteractions(cc)

  cc <- computeCommunProb(
    cc,
    type = "truncatedMean",
    trim = 0.1
  )

  cc <- filterCommunication(
    cc,
    min.cells = 10
  )

  cc <- computeCommunProbPathway(cc)

  cc <- aggregateNet(cc)

  cc <- netAnalysis_computeCentrality(
    cc,
    slot.name = "netP"
  )

  return(cc)

}

# =====================================
# Load objects
# =====================================

lung <- readRDS(
  "results/integrated/CLAD_CIT_integrated_annotated.rds"
)

macs <- readRDS(
  "results/macrophages/macs_annotated.rds"
)

# =====================================
# Add SPP1 and LILR labels to lung object
# =====================================

spp1_barcodes <- colnames(
  macs
)[
  macs$MacsLabels_1.0 == "SPP1"
]

lilr_barcodes <- colnames(
  macs
)[
  macs$MacsLabels_1.0 == "LILR"
]

lung$CellType_CellChat <- as.character(
  lung$CellType
)

lung$CellType_CellChat[
  spp1_barcodes
] <- "SPP1_macrophage"

lung$CellType_CellChat[
  lilr_barcodes
] <- "LILR_macrophage"

lung$CellType_CellChat <- factor(
  lung$CellType_CellChat
)

# =====================================
# Subset CLAD lungs
#
# CLAD = BOS + RAS
# =====================================

clad_lung <- subset(
  lung,
  subset = Condition %in% c(
    "BOS",
    "RAS"
  )
)

# =====================================
# Run CellChat
# =====================================

clad_cellchat <- run_cellchat(
  clad_lung,
  group.by = "CellType_CellChat"
)

# =====================================
# Communication table
# =====================================

clad_communications <- subsetCommunication(
  clad_cellchat
)

# =====================================
# Macrophage-restricted CellChat
# =====================================

selected_groups <- c(
  "monocyte/macrophage",
  "LILR_macrophage",
  "SPP1_macrophage"
)

clad_cellchat_macs <- subsetCellChat(
  clad_cellchat,
  idents = selected_groups
)

clad_cellchat_macs <- netAnalysis_computeCentrality(
  clad_cellchat_macs,
  slot.name = "netP"
)

# =====================================
# Save outputs
# =====================================
dir.create(
  "results/cellchat",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  lung,
  "results/cellchat/CLAD_CIT_lung_cellchat_metadata.rds"
)

saveRDS(
  clad_cellchat,
  "results/cellchat/CLAD_cellchat.rds"
)

saveRDS(
  clad_cellchat_macs,
  "results/cellchat/CLAD_cellchat_macrophages.rds"
)

saveRDS(
  clad_communications,
  "results/cellchat/CLAD_communications.rds"
)

write.csv(
  clad_communications,
  "results/cellchat/CLAD_communications.csv",
  row.names = FALSE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  "results/cellchat/sessionInfo_08.txt"
)
