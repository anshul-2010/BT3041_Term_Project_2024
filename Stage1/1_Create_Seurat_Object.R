# Importing all the necessary libraries

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

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\GSE262446_\\GSM8169898\\")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

## An object of class Seurat
## 20658 features across 9124 samples within 1 assay
## Active assay: RNA (20658 features, 0 variable features)
## 1 layer present: counts

# Let's examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "KRT8", "CAV1"), 1:30]

## 4 x 30 sparse Matrix of class "dgCMatrix"
##   [[ suppressing 30 column names ‘AAACCCAAGGTAAGGA-1’, ‘AAACCCAAGGTCCTGC-1’, ‘AAACCCACATGGGTCC-1’ ... ]]
                                                                               
## CD3D  . . .  . . .  .  . .  .  .  .  .  . . . .  . .  . .  . . .  . .  . . .  .
## TCL1A . . .  . . .  .  . .  .  .  .  .  . . . .  . .  . .  . . .  . .  . . .  .
## KRT8  4 5 3 31 8 . 44 23 . 18 16 10 31 10 4 8 . 20 4  3 8 18 9 7 16 4 17 . 3 37
## CAV1  . . 6  6 7 1  1  . .  2  .  .  .  . . . .  3 1 11 1  . . 1  1 .  . . .  1

dense.size <- object.size(as.matrix(pbmc.data))
dense.size
## 2631808424 bytes

sparse.size <- object.size(pbmc.data)
sparse.size
## 304564664 bytes

dense.size / sparse.size
## 8.6 bytes