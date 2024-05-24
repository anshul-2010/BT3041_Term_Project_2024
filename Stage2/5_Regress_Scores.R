# EXPLORE SOURCES OF UNWANTED VARIATION
# Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)

# Load cell cycle markers
load("C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\GSE262446_RAW\\cycle.rda")
# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, g2m.features = g2m_genes, s.features = s_genes)
# View cell cycle scores and phases assigned to cells                                 
View(seurat_phase@meta.data)

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Identify the 15 most highly variable genes
ranked_variable_genes <- VariableFeatures(seurat_phase)
top_genes <- ranked_variable_genes[1:15]
# Plot the average expression and variance of these genes
# With labels to indicate which genes are in the top 15
p <- VariableFeaturePlot(seurat_phase)
variable_features = LabelPoints(plot = p, points = top_genes, repel = TRUE)
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\variablefeature.jpeg", width=600)
variable_features
dev.off()

# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
dimension_plot_separate_phase <- DimPlot(seurat_phase, reduction = "pca", group.by = "Phase", split.by = "Phase")
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\dimensionplot_separatephase.jpeg", width=600)
dimension_plot_separate_phase
dev.off()

dimension_plot_phase <- DimPlot(seurat_phase, reduction = "pca", group.by = "Phase")
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\dimensionplot_phase.jpeg", width=600)
dimension_plot_phase
dev.off()

# Check quartile values
summary(seurat_phase@meta.data$mitoRatio)
# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, breaks=c(-Inf, 0.04398, 0.06875, 0.10182, Inf), labels=c("Low", "Medium", "Medium high", "High"))

dimension_plot_separate_mitoFr <- DimPlot(seurat_phase, reduction = "pca", group.by = "mitoFr", split.by = "mitoFr")
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\dimensionplot_separatemitoFr.jpeg", width=600)
dimension_plot_separate_mitoFr
dev.off()

dimension_plot_mitoFr <- DimPlot(seurat_phase, reduction = "pca", group.by = "mitoFr")
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\dimensionplot_mitoFr.jpeg", width=600)
dimension_plot_mitoFr
dev.off()

# Regress out cell cycle scores
seurat_phase <- ScaleData(seurat_phase, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat_phase))
# Not working in local system, too much space needed

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(seurat_phase, split.by = "sample")
options(future.globals.maxSize = 4000 * 1024^2)
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"), vst.flavor = "v2")
  }
# Save the split seurat object
saveRDS(split_seurat, "C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\split_seurat.rds")