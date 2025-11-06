###############################################################
# 📘 Single-cell RNA-seq QC and Normalization Workflow
# Author: Deeksha Sharma
# Purpose: Perform QC filtering, visualize metrics, and run SCTransform
###############################################################

##--------------------------------------------------------------
## 1. Load Required Libraries
##--------------------------------------------------------------

# Core data manipulation & visualization
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(ggpubr)

# Seurat ecosystem
library(Seurat)
library(SeuratObject)
library(glmGamPoi)

# Functional and enrichment analysis
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(enrichR)

# Utilities and helper packages
library(tibble)
library(stringr)
library(openxlsx)
library(writexl)
library(parallel)
library(future)
library(CellChat)
library(monocle3)
library(Nebulosa)
library(pheatmap)

# Optional: for troubleshooting or extensions
library(devtools)

##--------------------------------------------------------------
## 2. Load Pre-Merged Seurat Object
##--------------------------------------------------------------

# Load your previously merged object
obj <- readRDS(file = "D2A1plus_PYMT_1-merged-object.RDS")

cat("✅ Merged Seurat object successfully loaded.\n")

# Inspect object metadata
print(obj)
head(obj@meta.data)

##--------------------------------------------------------------
## 3. Quality Control (QC)
##--------------------------------------------------------------

cat("🔍 Running QC metrics...\n")

# 3.1 Record number of cells before filtering
ncells.before.qc.filter <- ncol(obj)
cat("Cells before QC filtering:", ncells.before.qc.filter, "\n")

# 3.2 Compute mitochondrial gene percentage
# Note: Mouse mitochondrial genes often start with "mt-"
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

# 3.3 Visualize QC metrics (Feature count, UMI count, MT%)
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 3.4 Scatter plots to explore relationships between QC features
plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

##--------------------------------------------------------------
## 4. Apply QC Filters
##--------------------------------------------------------------

cat("🧹 Filtering low-quality cells...\n")

# Define filtering thresholds
# Adjust thresholds depending on your tissue or sequencing depth
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 100000 & percent.mt < 10)

# 4.1 Record cell counts after filtering
ncells.after.qc.filter <- ncol(obj)

# 4.2 Calculate how many and what percentage were removed
ncells.filtered <- ncells.before.qc.filter - ncells.after.qc.filter
percent.filtered <- round((ncells.filtered / ncells.before.qc.filter) * 100, 2)

cat("Cells after filtering:", ncells.after.qc.filter, "\n")
cat("Filtered out:", ncells.filtered, "cells (", percent.filtered, "%)\n")

# 4.3 Store QC summary as a table
qc_results <- data.frame(
  Metric = c("Number of cells before QC",
             "Number of cells after QC",
             "Number of cells filtered",
             "Percent filtered"),
  Value = c(ncells.before.qc.filter,
            ncells.after.qc.filter,
            ncells.filtered,
            percent.filtered)
)

print(qc_results)

# 4.4 Save QC summary table
write.csv(
  qc_results,
  file = file.path(getwd(), "results", "tables", "1-qc_filtering_results.csv"),
  row.names = FALSE
)

# 4.5 Save filtered Seurat object
saveRDS(obj, file = "2-QCfiltered-object.RDS")
cat("✅ QC-filtered Seurat object saved.\n")

##--------------------------------------------------------------
## 5. Normalization via SCTransform
##--------------------------------------------------------------

cat("⚙️ Running SCTransform normalization (this may take a while)...\n")

# 5.1 Load filtered object (optional if not already in memory)
obj <- readRDS(file = "2-QCfiltered-object.RDS")

# 5.2 Increase memory limit for large datasets (prevent session abort)
options(future.globals.maxSize = 8000 * 1024^2)  # 8 GB memory cap (adjust as needed)

# 5.3 Run SCTransform
# This step normalizes data, removes technical noise, and regresses mitochondrial content
obj <- SCTransform(
  obj,
  vars.to.regress = "percent.mt",
  verbose = TRUE
)

cat("✅ SCTransform normalization complete.\n")

# 5.4 Save normalized object immediately — this step can take hours
saveRDS(obj, file = "3-SCtransform-object.RDS")

# 5.5 Confirm save and reload if needed
obj <- readRDS("3-SCtransform-object.RDS")

cat("📦 SCTransformed Seurat object loaded and ready for downstream analysis.\n")

##--------------------------------------------------------------
## ✅ Next Steps (Recommended)
##--------------------------------------------------------------
# - Run PCA: obj <- RunPCA(obj)
# - Find clusters: obj <- FindClusters(obj, resolution = 0.5)
# - Run UMAP: obj <- RunUMAP(obj, dims = 1:20)
# - Visualize clusters: DimPlot(obj, reduction = "umap", group.by = "treatment")
#
# Save these steps in a separate script to maintain modular workflows.
###############################################################
