# Simple-R-coding-for-single-cell-RNA-sequencing
###############################################################
# 📘 Single-cell RNA-seq Data Integration and Metadata Setup
# Author: [Deeksha Sharma]
# Purpose: Load, merge, and annotate 10X Genomics scRNA-seq data
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

# Seurat ecosystem for scRNA-seq analysis
library(Seurat)
library(SeuratObject)
library(glmGamPoi)

# Functional analysis
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(enrichR)

# Misc utilities
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

# Optional: For debugging or development
library(devtools)

##--------------------------------------------------------------
## 2. Set Up Directory Structure
##--------------------------------------------------------------

# Define and create output directories for results
tables_dir <- file.path(getwd(), "results", "tables")
figures_dir <- file.path(getwd(), "results", "figures")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

##Check directory
cat("Working directory:", getwd(), "\n")
cat("Results directories initialized.\n")

##--------------------------------------------------------------
## 3. Load and Create Seurat Objects
##--------------------------------------------------------------

# 🔹 Define base data paths for clarity
base_path <- "/home/deekshas/R/7.R DATA and Scripts/Burness Data"

# PYMT samples (Vehicle or control)
pymt_paths <- file.path(base_path, "PYMT DATA", paste0("11529-CH-", 1:6, "_sample_filtered_feature_bc_matrix"))

# D2A1 samples (Treatment and control)
d2a1_paths <- file.path(base_path, "D2A1_CT1", paste0("11529-CH-", 7:12, "_sample_filtered_feature_bc_matrix"))

# Function to read and create Seurat objects
create_obj <- function(path) {
  CreateSeuratObject(counts = Read10X(data.dir = path))
}

# Read all 12 datasets
all_paths <- c(pymt_paths, d2a1_paths)
seurat_list <- lapply(all_paths, create_obj)

# Assign sample names for tracking
sample_names <- paste0("obj", seq_along(seurat_list))
names(seurat_list) <- sample_names

##--------------------------------------------------------------
## 4. Merge Seurat Objects
##--------------------------------------------------------------

merged_samples <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = sample_names
)

cat("✅ Merged", length(seurat_list), "samples into one Seurat object.\n")
print(merged_samples)

# Save merged object
saveRDS(merged_samples, file = "all_1-merged-object.RDS")

# (Optional) Reload to confirm
merged_samples <- readRDS("all_1-merged-object.RDS")

##--------------------------------------------------------------
## 5. Join Layers and Cleanup
##--------------------------------------------------------------

obj <- JoinLayers(merged_samples)

# Remove intermediate objects to free memory
rm(list = sample_names)
gc()  # Run garbage collection

# Save joint object
saveRDS(obj, file = "all_1-joint-object.RDS")

# Reload to confirm
obj <- readRDS("all_1-joint-object.RDS")

##--------------------------------------------------------------
## 6. Annotate Metadata (Sample, Cell Line, Treatment)
##--------------------------------------------------------------

# Extract barcodes and sample IDs from cell names
rownames(obj[[]])  # View structure
sample.by.barcode <- word(rownames(obj[[]]), start = 1, end = 1, sep = fixed("_"))

# Add sample information to metadata
obj$sample <- sample.by.barcode

# Define cell line assignments
cell.line <- rep(NA, length(sample.by.barcode))
cell.line[sample.by.barcode %in% paste0("11529-CH-", 1:6)] <- "PYMT"
cell.line[sample.by.barcode %in% paste0("11529-CH-", 7:12)] <- "D2A1"

obj$cell_line <- cell.line

# Define treatment conditions
treatment <- rep(NA, length(sample.by.barcode))
treatment[sample.by.barcode %in% c("11529-CH-1", "11529-CH-2", "11529-CH-3", 
                                   "11529-CH-7", "11529-CH-8", "11529-CH-9")] <- "Vehicle"
treatment[sample.by.barcode %in% c("11529-CH-4", "11529-CH-5", "11529-CH-6", 
                                   "11529-CH-10", "11529-CH-11", "11529-CH-12")] <- "ZBC260"

obj$treatment <- treatment

##--------------------------------------------------------------
## 7. Save Final Annotated Object
##--------------------------------------------------------------

saveRDS(obj, file = "D2A1plus_PYMT_1-merged-object.RDS")

cat("✅ Final Seurat object saved with metadata annotations.\n")

# Inspect metadata
print(obj)
head(obj@meta.data)
                                   
                                   
                                   
                                   
                                   
                                   
                                   
                                   
