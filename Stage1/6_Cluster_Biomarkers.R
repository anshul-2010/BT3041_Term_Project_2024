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