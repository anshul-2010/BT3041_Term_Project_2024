# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat_integrated, vars = c("ident", "sample")) %>%
dplyr::count(ident, sample)

# Barplot of number of cells per cluster by sample
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\clustering_quality_control\\cell_count_distribution_2.jpeg", width=600)
ggplot(n_cells, aes(x=ident, y=n, fill=sample)) +
    geom_bar(position=position_dodge(), stat="identity") +
    geom_text(aes(label=n), vjust = -.2, position=position_dodge(1))
dev.off()

# UMAP of cells in each cluster by sample
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\clustering_quality_control\\cell_count_distribution_umap.jpeg", width=1000)
DimPlot(seurat_integrated, label=TRUE, split.by="sample", label.size=6) +NoLegend()
dev.off()

# Barplot of proportion of cells in each cluster by sample
ggplot(seurat_integrated@meta.data) +
    geom_bar(aes(x=integrated_snn_res.0.2, fill=sample), position=position_fill())

# Explore whether clusters segregate by cell cycle phase
DimPlot(seurat_integrated, label = TRUE, split.by = "Phase")  + NoLegend()

# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")
FeaturePlot(seurat_integrated, reduction = "umap", features = metrics, pt.size = 0.2, order = TRUE, min.cutoff = 'q10', label = TRUE, label.size=6)

# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:16), "ident", "umap_1", "umap_2")
# Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated, vars = columns)

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_integrated, vars = c("ident", "umap_1", "umap_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(umap_1), y=mean(umap_2))
  
# Plotting a UMAP plot for each of the PCs
map(paste0("PC_", 1:16), function(pc){
        ggplot(pc_data, 
               aes(umap_1, umap_2)) +
                geom_point(aes_string(color=pc), 
                           alpha = 0.7) +
                scale_color_gradient(guide = FALSE, 
                                     low = "grey90", 
                                     high = "blue")  +
                geom_text(data=umap_label, 
                          aes(label=ident, x, y)) +
                ggtitle(pc)
}) %>% 
        plot_grid(plotlist = .)

# Examine PCA results 
print(seurat_integrated[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(object = seurat_integrated, reduction = "umap", label = TRUE) + NoLegend()

# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("KRT8", "SBSN", "DST", "MT1G", "GPX2"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

seurat_integrated.markers <- FindAllMarkers(seurat_integrated, only.pos=TRUE)
seurat_integrated %>%
group_by(cluster) %>%
dplyr::filter(avg_log2FC > 1)

markers <- list()

for (i in 0:8) {
  # Identify markers for each cluster compared to all other cells
  markers <- FindMarkers(seurat_integrated,
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
cluster8_markers <- all_markers[["cluster7.markers"]]

# Take the first 5 gene expressions from the above obtained clusters of markers.
VlnPlot(seurat_integratedatures = c("KRT8", "KRT19"))
VlnPlot(seurat_integratedatures = c("SERPINB3", "MAL", "S100A8", "KRT1", "SBSN"))
VlnPlot(seurat_integratedatures = c("DST", "COL17A1", "CAV1"))
VlnPlot(seurat_integratedatures = c("ADIRF", "CLU", "DLK2"))
VlnPlot(seurat_integratedatures = c("S100A2", "MT1G", "STMN1", "KIAA0101"))
VlnPlot(seurat_integratedatures = c("GPX2", "AGR2", "FAM3B"))
VlnPlot(seurat_integratedatures = c("B2M", "CD24"))
VlnPlot(seurat_integratedatures = c("HLA-DRA", "CD74", "HLA-DPA1", "HLA-DPB1", "HLA-DRB1"))
VlnPlot(seurat_integratedatures = c("CAV1", "TIMP2"))
# 8th was very sparse.

FeaturePlot(seurat_integrated, features = c("KRT8", "SBSN", "DST", "KIAA0101", "GPX2", "B2M", "HLA-DPB1", "TIMP2"), label=TRUE, label.size=5)

seurat_integrated.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
DoHeatmap(seurat_integrated, features = top10$gene) + NoLegend()

new.cluster.ids <- c("Epithelial_1.1", "Epithelial_1.2", "EndothelialA", "EndothelialB", "Astrocyte", "T_cell", "Keratinocyte", "mDC_cell")

names(new.cluster.ids) <- levels(seurat_integrated)
seurat_integrated <- RenameIdents(seurat_integrated, new.cluster.ids)
DimPlot(seurat_integrated, reduction="umap", label=TRUE, pt.size=1)

# Create dotplot based on RNA expression
DotPlot(seurat_integrated, markers, assay="RNA")