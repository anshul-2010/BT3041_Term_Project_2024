# INTEGRATION TESTING USING CCA
seura_phase <- RunUMAP(seurat_phase, dims = 1:40,reduction = "pca")
seura_plot <- DimPlot(seura_phase)
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\seura_plot.jpeg", width=600)
seura_plot
dev.off()

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)
# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, anchor.features = integ_features)
# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, normalization.method = "SCT", anchor.features = integ_features)
# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
saveRDS(seurat_integrated, "C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\integrated_seurat.rds")

# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)
# Plot PCA
integrated_pca <- PCAPlot(seurat_integrated, split.by = "sample")
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\pca_plot_integrated.jpeg", width=600)
integrated_pca
dev.off()

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40, reduction = "pca")
# Plot UMAP
integrated_umap <- DimPlot(seurat_integrated)
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\umap_integrated.jpeg", width=600)
integrated_umap
dev.off()