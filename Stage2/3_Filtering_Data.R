# FILTERING

# Cell level

# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, subset = (nUMI >= 500) & (nGene >= 250) & (log10GenesPerUMI > 0.80) & (mitoRatio < 0.20))

# Gene level

filtered_seurat_join <- JoinLayers(filtered_seurat)
# Extract counts
counts <- GetAssayData(object = filterd_seurat_join, layer = "counts")
# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)