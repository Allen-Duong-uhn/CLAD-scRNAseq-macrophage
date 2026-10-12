# =====================================
#   Differential expression analyses of macrophages.
#
# Workflow:
#   1. Compare CLAD versus CIT with Wilcoxon rank-sum and MAST analysis
#   2. Compare BOS vs RAS with Wilcoxon rank-sum and MAST analysis
#   3. Save outputs
#
# Inputs:
#   results/macrophages/macs_annotated.rds
#
# Outputs:
#   results/deg/mast_cit_vs_clad.rds
#   results/mast_cit_vs_clad.rds
#
# =====================================

# load libraries
library(Seurat)
library(dplyr)

# Load macrophage object

macs <- readRDS(
  "results/macrophages/macs_annotated.rds"
)

# =====================================
# Comparing CLAD versus CIT
# =====================================

Idents(macs) <- "CLADstatus"

markers_clad_vs_ctrl <- FindAllMarkers(
  macs,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
)

mast_clad_vs_ctrl <- FindMarkers(
  macs,
  ident.1 = "Control",
  ident.2 = "CLAD",
  test.use = "MAST"
)

mast_clad_vs_ctrl <- mast_clad_vs_ctrl %>%
  tibble::rownames_to_column(
    var = "gene"
  ) %>%
  arrange(
    p_val_adj
  )

# =====================================
# Comparing BOS versus RAS
# =====================================

Idents(macs) <- "Condition"

macs_bos_ras <- subset(
  macs,
  ident = c(
    "BOS",
    "RAS"
  )
)

markers_bos_vs_ras <- FindAllMarkers(
  macs_bos_ras,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
)

mast_bos_vs_ras <- FindMarkers(
  macs_bos_ras,
  ident.1 = "BOS",
  ident.2 = "RAS",
  test.use = "MAST"
)

mast_bos_vs_ras <- mast_bos_vs_ras %>%
  tibble::rownames_to_column(
    var = "gene"
  ) %>%
  arrange(
    p_val_adj
  )

# =====================================
# Save outputs
# =====================================
dir.create(
  "results/deg",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  markers_clad_vs_ctrl,
  "results/deg/markers_CLAD_vs_Ctrl.rds"
)

saveRDS(
  mast_clad_vs_ctrl,
  "results/deg/mast_CLAD_vs_Ctrl.rds"
)

saveRDS(
  markers_bos_vs_ras,
  "results/deg/markers_BOS_vs_RAS.rds"
)

saveRDS(
  mast_bos_vs_ras,
  "results/deg/mast_BOS_vs_RAS.rds"
)

write.csv(
  mast_ctrl_vs_clad,
  "results/deg/mast_CLAD_vs_Ctrl.csv",
  row.names = FALSE
)

write.csv(
  mast_bos_vs_ras,
  "results/deg/mast_BOS_vs_RAS.csv",
  row.names = FALSE
)

# =====================================
# Session information
# =====================================
writeLines(
  capture.output(
    sessionInfo()
  ),
  "results/deg/sessionInfo05.txt"
)
