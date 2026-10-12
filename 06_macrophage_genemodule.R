# =====================================
#   Calculate macrophage functional module scores.
#
# Workflow:
#   1. Input gene modules as lists
#   2. Calculate module score for all gene sets on macrophage object
#   3. Average the gene module expression at the sample level for pro-inflammatory and fibrosis modules
#   4. Save outputs
#
# Inputs:
#   results/macrophages/macs_annotated.rds
#
# Outputs:
#   results/module_scores/proinflam_avg_scores_wide.rds
#   results/module_scores/fibrosis_avg_scores_wide.rds
#   results/module_scores/proinflam_avg_scores_wide.csv
#   results/module_scores/fibrosis_avg_scores_wide.csv
#   results/macrophages/macs_annotated_module.rds
#
# =====================================

# load packages
library(Seurat)
library(dplyr)
library(tidyr)

# load macrophages
macs <- readRDS(
  "results/macrophages/macs_annotated.rds"
)

# =====================================
# Input gene modules
# =====================================

proinflam_genes <- c(
  "ACOD1","PLA2G7","FMNL2","FYB1","AOAH",
  "F13A1","ARG1","HP","LCN2","ARID5B",
  "HLA-DRA","SASH1","PRDX6","GBP2",
  "CTSS","CIITA","FAM126B","AIF1",
  "CXCL2","CCL4","CCL3","IL1B",
  "LILRA5","CD274","IGSF6",
  "IL1RN","SLFN5","RNF213","CCRL2",
  "IL36G","GADD45A","BASP1","CLEC4E",
  "HCAR2","IRGM"
)

fibrosis_genes <- c(
  "CTSD","FABP5","PSAP","CD63",
  "GPNMB","CTSS","CSTB","LGALS3",
  "SGK1","SPP1","CTSL","CD9",
  "FTH1","MSRB1","ATP6V0D2",
  "TMSB10","TXN1","CHI3L1"
)

engulfment_genes <- c(
  "SYK","FYN","PIK3CG","PIK3CD",
  "INPP5D","RAC1","CDC42","RHOA",
  "ABI3","DOCK5","ARAP2","PSTPIP2"
)

M1_genes <- c(
  "TYMP","XAF1","PSMA2","OAS2","NAMPT",
  "GADD45G","SLC31A2","PFKP","PSMB9","INHBA",
  "APOL6","SLC2A6","SPHK1","PFKFB3","PSME2",
  "EDN1","APOL2","TNFSF10","IL6","CCL20",
  "SLCO4A1","AK3","APOL1","PDGFA","IRF7",
  "BIRC3","IRF1","ATF3","IL15","HSD11B1",
  "IL12B","SLC7A5","IGFBP4","CHI3L2","APOL3",
  "HESX1","IL7R","FAS","OASL","CCL15",
  "IL15RA","VCAN","CCL5","TNF","IL2RA",
  "PTX3","PLA1A","BCL2A1","CXCL9","CXCL10",
  "CCR7","IDO1","CCL19","CXCL11"
)

M2_genes <- c(
  "MAF","MS4A4A","MRC1","FGL2",
  "CCL13","CD36","MS4A6A","CLEC4A",
  "CD209","P2RY14","FN1","SLC38A6",
  "CLEC7A","CCL18","CXCR4","MSR1",
  "CCL23","HS3ST1","TLR5","CHN2",
  "HRH1","ALOX15","EGR2","TGFBR2",
  "LPAR6","SLC4A7","SEPP1","CA2",
  "TGFBI","SLCO2B1","LTA4H","IGF1",
  "TPST2","CERK","HS3ST2","GPR86",
  "ADK","HNMT","GAS7","CTSC",
  "HEXB","LIPA"
)

gene_modules <- list(
  Proinflam = proinflam_genes,
  Fibrosis = fibrosis_genes,
  Engulfment = engulfment_genes,
  M1 = M1_genes,
  M2 = M2_genes
)

# =====================================
# Calculate module scores
# =====================================

for(module_name in names(gene_modules)){

  macs <- AddModuleScore(
    object = macs,
    features = list(
      gene_modules[[module_name]]
    ),
    name = paste0(
      module_name,
      "_score"
    ),
    ctrl = 5
  )

}

# =====================================
# Averaging module scores at the sample-level
# =====================================

meta <- macs@meta.data

proinflam_avg_scores <- meta %>%
  group_by(
    orig.ident,
    MacsLabels_1.0
  ) %>%
  summarise(
    mean_module_score = mean(
      Proinflam_score1
    ),
    n_cells = n(),
    .groups = "drop"
  )

fibrosis_avg_scores <- meta %>%
  group_by(
    orig.ident,
    MacsLabels_1.0
  ) %>%
  summarise(
    mean_module_score = mean(
      Fibrosis_score1
    ),
    n_cells = n(),
    .groups = "drop"
  )

sample_condition <- meta %>%
  select(
    orig.ident,
    Condition
  ) %>%
  distinct()

proinflam_avg_scores <- proinflam_avg_scores %>%
  left_join(
    sample_condition,
    by = "orig.ident"
  )

fibrosis_avg_scores <- fibrosis_avg_scores %>%
  left_join(
    sample_condition,
    by = "orig.ident"
  )

# convert to wide-format tables for easier transfer of data to GraphPad Prism
proinflam_avg_scores_wide <- proinflam_avg_scores %>%
  select(
    orig.ident,
    MacsLabels_1.0,
    mean_module_score
  ) %>%
  pivot_wider(
    names_from = MacsLabels_1.0,
    values_from = mean_module_score
  )

fibrosis_avg_scores_wide <- fibrosis_avg_scores %>%
  select(
    orig.ident,
    MacsLabels_1.0,
    mean_module_score
  ) %>%
  pivot_wider(
    names_from = MacsLabels_1.0,
    values_from = mean_module_score
  )
  
# =====================================
# Save outputs
# =====================================
dir.create(
  "results/module_scores",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  proinflam_avg_scores_wide,
  "results/module_scores/proinflam_avg_scores_wide.rds"
)

saveRDS(
  fibrosis_avg_scores_wide,
  "results/module_scores/fibrosis_avg_scores_wide.rds"
)

write.csv(
  proinflam_avg_scores_wide,
  "results/module_scores/proinflam_avg_scores_wide.csv",
  row.names = FALSE
)

write.csv(
  fibrosis_avg_scores_wide,
  "results/module_scores/fibrosis_avg_scores_wide.csv",
  row.names = FALSE
)

saveRDS(
   macs,
   "results/macrophages/macs_annotated_module.rds"
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  "results/module_scores/sessionInfo_06.txt"
)
