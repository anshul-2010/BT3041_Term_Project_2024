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