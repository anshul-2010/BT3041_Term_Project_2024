# CLUSTERING CELLS
DimHeatmap(seurat_integrated, dims = 1:9, cells = 500, balanced = TRUE)
print(x = seurat_integrated[["pca"]], dims = 1:10, nfeatures = 5)

jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\umap_integrated_elbow.jpeg", width=600)
ElbowPlot(object = seurat_integrated, ndims = 40)
dev.off()

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, dims = 1:20)
seurat_integrated <- FindClusters(object = seurat_integrated, resolution = c(0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4))

seurat_integration@meta.data %>%
View()

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.2"
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\umap_integrated_2.jpeg", width=600)
# Plot the UMAP
DimPlot(seurat_integrated, reduction = "umap", label = TRUE, label.size = 6)
dev.off()

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\umap_integrated_4.jpeg", width=600)
# Plot the UMAP
DimPlot(seurat_integrated, reduction = "umap", label = TRUE, label.size = 6)
dev.off()

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.5"
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\umap_integrated_5.jpeg", width=600)
# Plot the UMAP
DimPlot(seurat_integrated, reduction = "umap", label = TRUE, label.size = 6)
dev.off()

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.6"
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\umap_integrated_6.jpeg", width=600)
# Plot the UMAP
DimPlot(seurat_integrated, reduction = "umap", label = TRUE, label.size = 6)
dev.off()

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\umap_integrated_8.jpeg", width=600)
# Plot the UMAP
DimPlot(seurat_integrated, reduction = "umap", label = TRUE, label.size = 6)
dev.off()

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.1.0"
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\umap_integrated_10.jpeg", width=600)
# Plot the UMAP
DimPlot(seurat_integrated, reduction = "umap", label = TRUE, label.size = 6)
dev.off()

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.1.2"
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\umap_integrated_12.jpeg", width=600)
# Plot the UMAP
DimPlot(seurat_integrated, reduction = "umap", label = TRUE, label.size = 6)
dev.off()

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.1.4"
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\umap_integrated_14.jpeg", width=600)
# Plot the UMAP
DimPlot(seurat_integrated, reduction = "umap", label = TRUE, label.size = 6)
dev.off()