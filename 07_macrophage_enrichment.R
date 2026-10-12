# =====================================
# Script: 07_macrophage_enrichment.R
#
# Workflow:
#   1. Load annotated macrophages
#   2. Generate subtype-specific marker lists
#   3. Build ranked gene lists
#   4. Perform GO Biological Process GSEA
#   5. Save outputs
#
# Inputs:
#   results/macrophages/macs_annotated.rds
#
# Outputs:
#   results/enrichment/marker_lists.rds
#   results/enrichment/ranked_gene_lists.rds
#   results/enrichment/gsea_results.rds
#
# Notes:
#   Marker discovery is performed here directly using
#   FindMarkers() for each macrophage subtype.
#
# =====================================

# =====================================
# Libraries
# =====================================

library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# =====================================
# Load annotated macrophages
# =====================================

macs <- readRDS(
  "results/macrophages/macs_annotated.rds"
)

Idents(macs) <- "MacsLabels_1.0"

# =====================================
# Generate subtype-specific markers
# =====================================
subtypes <- levels(macs)

marker_lists <- lapply(
  subtypes,
  function(x){

    FindMarkers(
      object = macs,
      ident.1 = x,
      only.pos = FALSE
    )

  }
)

names(marker_lists) <- subtypes

# =====================================
# Build ranked gene lists
#
# Uses avg_log2FC as ranking metric.
# =====================================

ranked_gene_lists <- lapply(
  marker_lists,
  function(df){

    gene_list <- df$avg_log2FC

    names(gene_list) <- rownames(df)

    gene_list <- gene_list[
      !is.na(gene_list)
    ]

    gene_list <- sort(
      gene_list,
      decreasing = TRUE
    )

    gene_list

  }
)

# =====================================
# GO Biological Process GSEA
# =====================================

gsea_results <- lapply(
  ranked_gene_lists,
  function(gene_list){

    gseGO(
      geneList = gene_list,
      ont = "BP",
      keyType = "SYMBOL",
      OrgDb = org.Hs.eg.db,
      minGSSize = 10,
      maxGSSize = 800,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = FALSE
    )

  }
)

names(gsea_results) <- names(
  ranked_gene_lists
)

# =====================================
# Save outputs
# =====================================
dir.create(
  "results/enrichment",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  marker_lists,
  "results/enrichment/marker_lists.rds"
)

saveRDS(
  ranked_gene_lists,
  "results/enrichment/ranked_gene_lists.rds"
)

saveRDS(
  gsea_results,
  "results/enrichment/gsea_results.rds"
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  "results/enrichment/sessionInfo_07.txt"
)
