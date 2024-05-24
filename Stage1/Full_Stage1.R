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

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Show the QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

##                      orig.ident nCount_RNA nFeature_RNA percent.mt
## AAACCCAAGGTAAGGA-1     pbmc3k      13920         3252  26.918103
## AAACCCAAGGTCCTGC-1     pbmc3k      12524         2984  11.609709
## AAACCCACATGGGTCC-1     pbmc3k       9196         2393   5.491518
## AAACCCAGTCACCGCA-1     pbmc3k      31777         5208   7.077446
## AAACCCAGTCTGTAGT-1     pbmc3k       9006         2284   4.696869

# Visualize QC metrics as a violin plot
violin_plot <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\violinplot.jpeg", width=600, height=350)
violin_plot
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e., columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
feature_scatter <- plot1 + plot2
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\featurescatter.jpeg", width=1000, height=350)
feature_scatter
dev.off()

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
## Normalizing layer: counts
## Performing log-normalization
## 0%   10   20   30   40   50   60   70   80   90   100%
## [----|----|----|----|----|----|----|----|----|----|
## **************************************************|

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
## Finding variable features for layer counts
## Calculating gene variances
## 0%   10   20   30   40   50   60   70   80   90   100%
## [----|----|----|----|----|----|----|----|----|----|
## **************************************************|
## Calculating feature variances of standardized and clipped values
## 0%   10   20   30   40   50   60   70   80   90   100%
## [----|----|----|----|----|----|----|----|----|----|
## **************************************************|

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
feature_plot <- plot1 + plot2
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\variablefeature.jpeg",  width=1000, height=350)
feature_plot
dev.off()

all.genes <- rownames(pbmc)
# pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
## PC_ 1 
## Positive:  HSPB1, KRT19, S100A2, DEFB1, TACSTD2 
## Negative:  CD52, PTPRC, HCST, CXCR4, CD48 
## PC_ 2 
## Positive:  SPRR1B, KRT6A, S100A8, KRT13, KRT6C 
## Negative:  MGP, SERPING1, CCDC80, CXCL12, FN1 
## PC_ 3 
## Positive:  KRT6C, SPRR1B, KRT6B, SPRR2A, S100A8 
## Negative:  FTH1, KRT8, KRT19, DEFB1, KRT18 
## PC_ 4 
## Positive:  CLDN7, CST3, SMIM22, WFDC2, TSPAN1 
## Negative:  KRT15, IGFBP6, DST, COL17A1, TNC 
## PC_ 5 
## Positive:  TSPAN1, CYP4B1, WFDC2, CLDN7, S100A4 
## Negative:  IL3RA, VWF, AQP1, EGFL7, TSPAN7

vizdim <- VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\pcdimloading.jpeg",  width=600, height=350)
vizdim
dev.off()

dim_plot <- DimPlot(pbmc, reduction = "pca") + NoLegend()
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\dimensionplot.jpeg",  width=600, height=350)
dim_plot
dev.off()

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

elbow_plot <- ElbowPlot(pbmc)
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\elbowplot.jpeg",  width=600, height=350)
elbow_plot
dev.off()

pbmc <- FindNeighbors(pbmc, dims = 1:7)
## Computing nearest neighbor graph
## Computing SNN
pbmc <- FindClusters(pbmc, resolution = 0.5)
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

## Number of nodes: 948
## Number of edges: 25403

## Running Louvain algorithm...
## 0%   10   20   30   40   50   60   70   80   90   100%
## [----|----|----|----|----|----|----|----|----|----|
## **************************************************|
## Maximum modularity in 10 random starts: 0.8815
## Number of communities: 8
## Elapsed time: 0 seconds

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
## AAACCCAGTCTGTAGT-1 AAACCCATCAGGACGA-1 AAAGGATAGGAGAGTA-1 AAAGGATTCTCTAAGG-1 
##                  2                  2                  1                  3 
## AAAGGTAAGATCCCGC-1 
##                  4 
## Levels: 0 1 2 3 4 5 6 7

pbmc <- RunUMAP(pbmc, dims = 1:10)
# Note that label=TRUE or use LabelClusters function to help label individual clusters
## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
## This message will be shown once per session
## 00:25:56 UMAP embedding parameters a = 0.9922 b = 1.112
## 00:25:56 Read 948 rows and found 10 numeric columns
## 00:25:56 Using Annoy for neighbor search, n_neighbors = 30
## 00:25:56 Building Annoy index with metric = cosine, n_trees = 50
## 0%   10   20   30   40   50   60   70   80   90   100%
## [----|----|----|----|----|----|----|----|----|----|
## **************************************************|
## 00:25:56 Writing NN index file to temp file C:\Users\Dell\AppData\Local\Temp\RtmpKGBii6\file40904d6d5075
## 00:25:56 Searching Annoy index using 1 thread, search_k = 3000
## 00:25:56 Annoy recall = 100%
## 00:25:57 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
## 00:25:58 Initializing from normalized Laplacian + noise (using RSpectra)
## 00:25:58 Commencing optimization for 500 epochs, with 34494 positive edges
## Using method 'umap'
## 0%   10   20   30   40   50   60   70   80   90   100%
## [----|----|----|----|----|----|----|----|----|----|
## **************************************************|
## 00:26:01 Optimization finished

umap <- DimPlot(pbmc, reduction = "umap")
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\umap.jpeg", width=600, height=350)
umap
dev.off()

saveRDS(pbmc, file = "C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\pbmc_tutorial.rds")

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)
##                p_val avg_log2FC pct.1 pct.2    p_val_adj
## DST     1.425910e-50   1.852940 0.910 0.335 2.945644e-46
## TNC     2.189956e-48   2.100784 0.677 0.135 4.524012e-44
## CAV1    8.679581e-48   1.804594 0.845 0.251 1.793028e-43
## KRT15   1.686067e-45   1.853981 0.935 0.473 3.483076e-41
## COL17A1 3.494946e-41   1.742382 0.729 0.201 7.219858e-37

pbmc.markers <- FindAllMarkers(pbmc, only.pos=TRUE)
pbmc.markers %>%
group_by(cluster) %>%
dplyr::filter(avg_log2FC > 1)

for (i in 0:7) {
  # Identify markers for each cluster compared to all other cells
  markers <- FindMarkers(pbmc,
                         ident.1 = i,
                         logfc.threshold = 0.5,      # Higher log-fold change threshold
                         test.use = "roc",        # Using Wilcoxon Rank Sum test
       
                  only.pos = TRUE,            # Only consider positive markers
                         min.pct = 0.5,             # Minimum fraction of cells in the cluster expressing the gene
                         min.diff.pct = 0.3)         # Minimum difference in fraction of expressing cells between groups
  
  # Store the result in the list, using the cluster number as the name
  all_markers[[paste0("cluster", i, ".markers")]] <- markers
}

cluster0_markers <- all_markers[["cluster0.markers"]]
cluster1_markers <- all_markers[["cluster1.markers"]]
cluster2_markers <- all_markers[["cluster2.markers"]]
cluster3_markers <- all_markers[["cluster3.markers"]]
cluster4_markers <- all_markers[["cluster4.markers"]]
cluster5_markers <- all_markers[["cluster5.markers"]]
cluster6_markers <- all_markers[["cluster6.markers"]]
cluster7_markers <- all_markers[["cluster7.markers"]]

# Take the first 5 gene expressions from the above obtained clusters of markers.
VlnPlot(pbmc, features = c("KRT8", "TNFSF10", "KRT19", "NUPR1", "AGR2"))
VlnPlot(pbmc, features = c("SPRR1B", "KRT6B", "S100A8", "IL1RN", "SBSN"))
VlnPlot(pbmc, features = c("DST", "TNC", "GPC3", "ADIRF", "CAV1"))
VlnPlot(pbmc, features = c("ADIRF", "KRT15", "CLU", "SERPINF1", "DLK2"))
VlnPlot(pbmc, features = c("S100A2", "KRT5", "WFDC2", "CXCL14", "IGFBP6"))
VlnPlot(pbmc, features = c("BTG1", "CD52", "IL32", "PTPRC", "CD3D"))
VlnPlot(pbmc, features = c("KRT7", "CD24", "CLDN7", "TSPAN1", "KRT8"))
VlnPlot(pbmc, features = c("HLA-DRA", "CD74", "HLA-DPA1", "HLA-DPB1", "HLA-DRB1"))

FeaturePlot(pbmc, features = c("AGR2", "SPRR1B", "TNC", "DLK2", "CXCL14", "IL32", "KRT7", "HLA-DRA"), label=TRUE, label.size=5)

pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

new.cluster.ids <- c("Epithelial_1.1", "Epithelial_1.2", "EndothelialA", "EndothelialB", "Astrocyte", "T_cell", "Keratinocyte", "mDC_cell")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction="umap", label=TRUE, pt.size=1.0) + NoLegend()
