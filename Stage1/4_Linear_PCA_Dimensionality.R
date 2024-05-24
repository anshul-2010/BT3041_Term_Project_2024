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