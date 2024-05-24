library(SingleCellExperiment)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)

# QUALITY CONTROL SETUP
# Load the PBMC dataset and create a Seurat object for each sample
ctrl_raw_feature_bc_matrix.data <- Read10X(data.dir = "C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\GSE262446_RAW\\ctrl_raw_features_bc_matrix\\")
ctrl_raw_feature_bc_matrix <- CreateSeuratObject(counts = ctrl_raw_feature_bc_matrix.data, project = "ctrl_raw_features_bc_matrix", min.features = 200)

stim_raw_feature_bc_matrix.data <- Read10X(data.dir = "C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\GSE262446_RAW\\stim_raw_features_bc_matrix\\")
stim_raw_feature_bc_matrix <- CreateSeuratObject(counts = stim_raw_feature_bc_matrix.data, project = "stim_raw_features_bc_matrix", min.features = 200)

# Check the metadata in the new Seurat objects
head(ctrl_raw_feature_bc_matrix@meta.data)
head(stim_raw_feature_bc_matrix@meta.data)

# Create a merged Seurat object
merged_seurat <- merge(x = ctrl_raw_feature_bc_matrix, y = stim_raw_feature_bc_matrix, add.cell.id = c("ctrl", "stim"))
# Check that the merged object has appropriate sample-specific prefixes
head(merged_seurat@meta.data)